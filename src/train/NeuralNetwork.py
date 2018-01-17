import numpy as np
from keras.callbacks import EarlyStopping, ModelCheckpoint
from keras.models import load_model
from keras.utils import plot_model
from sklearn.metrics import confusion_matrix, precision_score, recall_score, f1_score
from train.CustomLoss import weighted_binary_loss

import train.Divider as div
import train.Util as ut
from train import DataAugmenter as da

from train.JsonParser import JsonParser

class NeuralNetwork:
    def __init__(self, conditions=[], set_division=(0.7, 0.1, 0.2), parameters=JsonParser('src/train/parameters.json'), verbose=True):
        """
        Create deep narrow network
        :param conditions: default: no conditions for loading data
        :param set_division: default: train 70%, test 10%, validate 20%
        :param parameters: load default parameters file
        """
        # Assume no 'image only'
        self.image_only = False
        self.parameters = parameters
        # Determine number of entries used
        self.used_size = self.__determine_used_size()
        self.conditions = conditions
        self.entries = []
        self.names, self.peptides_data, self.mhcs_data, self.results = [], [], [], []
        self.total_size = 0
        self.set_division = set_division
        self.train_indices, self.test_indices, self.val_indices = [], [], []
        self.batch_size = 64
        self.store_weights_file = 'temp.hdf5'
        self.initial_weights = None
        self.peptides, self.mhcs = [], []
        self.augmenter = da.DataAugmenter(peptide_length=9)
        self.verbose = verbose
        if self.verbose:
            print("Created NN of type", self.type())

    def set_verbose(self, v):
        """
        Set verbosity of NN
        """
        self.verbose = v

    def save(self, filename=""):
        """
        Save the NN
        """
        assert(self.model is not None)
        if not filename:
            filename = self.p_model_filename()

        self.model.save(filename)

    def load(self, filename=""):
        """
        Load a saved NN
        """
        if not filename:
            filename = self.p_model_filename()

        self.model = load_model(
            filename, 
            custom_objects={'weighted_binary_loss': weighted_binary_loss}
        )

    def load_data(self, datafile=''):
        """
        Load data from files and return list of (unique) names
        """
        self.__load_raw_data(datafile)
        return self.names

    def set_test_type(self, type='random', filter=None):
        """
        Set the type of testing. Options are:
        - 'random': filter is ignored (default)
        - 'list': filter is list of protein names
        - 'family': filter should be list of families
        - 'type': filter should be list of type (NOTE: for performance, limit to 1 type in filter)

        e.g.
            set_test_type(type='family' filter=['HLA-A*01', 'HLA-A*02'])
            will put all sequences of HLA-A*01 and HLA-A*02 families in test set.
        """
        assert(self.divider is not None)
        self.divider.set_test_type(type, filter)

    def train(self, indices, oa='sigmoid', input_weights_file='', output_weights_file=''):
        (self.train_indices, self.test_indices, self.val_indices) = indices
        assert(self.entries)

        (self.augmented_mhc, self.augmented_peptides, self.augmented_binding_results, self.augmented_indices) = self.__augment_sets()

        # Collect actual data used in training
        self.__collect_peptide_data()
        self.__collect_mhc_data()
        self.__collect_label_data()

        # Load model
        self.model = self.p_load_model(oa)
        self.initial_weights = self.model.get_weights()

        # Optionally load pretrained weights
        if input_weights_file:
            self.model.load_weights(input_weights_file)

        # Train the model if there is a training set
        if self.train_indices:
            mhc_train_array = self.mhcs.take(self.train_indices, 0)
            peptides_train_array = self.peptides.take(self.train_indices, 0)
            training_data = [mhc_train_array, peptides_train_array]

            output_weights = self.store_weights_file
            if output_weights_file:
                output_weights = output_weights_file
            self.model.fit(training_data,
                           self.binding_result.take(self.train_indices, 0),
                           epochs=self.parameters.epochs(),
                           batch_size=self.batch_size,
                           callbacks=self.__callbacks(output_weights),
                           validation_data=([self.mhcs.take(self.val_indices, 0), self.peptides.take(self.val_indices, 0)],
                                       self.binding_result.take(self.val_indices, 0))
                           )

        # Evaluate the model
        return self.__evaluate()

    def predict(self, entry):
        """
        Predict @p entry's label
        """
        assert(self.model is not None)
        (name, mhc, peptide) = entry

        # Convert MHC sequence using BLOSUM matrix
        ndmhc = np.zeros((1, len(mhc), 20))
        blossom_scores = [ut.to_blossom(mhc)]
        ndmhc = self.__fill_ndarray(ndmhc, blossom_scores)

        # Convert peptide sequence using BLOSUM matrix
        mhpeptides = np.zeros((1, len(peptide), 20))
        blossom_scores = [ut.to_blossom(peptide)]
        mhpeptides = self.__fill_ndarray(mhpeptides, blossom_scores)

        # Perform prediction
        return self.model.predict([ndmhc, mhpeptides]).flatten()

    def type(self):
        """
        Return NN type
        """
        return None

    def p_load_model(self, oa):
        """
        Load the model
        """
        return None

    def p_model_filename(self):
        """
        Return the models filename
        :return:
        """
        return None

    def reinitialize_weights(self):
        """
        Reset the models original weights (pre-training)

        """
        if self.initial_weights:
            self.model.set_weights(self.initial_weights)

    def __evaluate(self):
        """
        Evaluate the trained NN
        """
        # Take the test data
        mhc_test_array = self.mhcs.take(self.test_indices, 0)
        peptides_test_array = self.peptides.take(self.test_indices, 0)
        test_data = [mhc_test_array, peptides_test_array]

        if test_data[0].size > 0 and test_data[1].size > 0:
            # Load physchem data if available
            if self.parameters.physchem() > 0:
                test_data = [mhc_test_array[:, :, :self.parameters.aa_size()],
                             mhc_test_array[:, :, self.parameters.aa_size():],
                             peptides_test_array[:, :, :self.parameters.aa_size()],
                             peptides_test_array[:, :, self.parameters.aa_size():]
                             ]

            # Perform keras evaluation (= validation metrics)
            scores = self.model.evaluate(test_data, self.binding_result.take(self.test_indices, 0))
            # Perform prediction on test data
            preds = self.model.predict(test_data).flatten()
            labels = self.binding_result.take(self.test_indices, 0).flatten().astype(np.int)
            prec = precision_score(labels, [custom_round(s) for s in preds])

            # Print performance
            if self.verbose:
                for i in range(len(scores)):
                    print("%s: %.2f%%" % (self.model.metrics_names[i], scores[i] * 100))
                print(np.bincount(labels))
                print("skPrecision: %.2f%%" % (prec * 100))
                print("skRecall: %.2f%%" % (recall_score(labels, [custom_round(s) for s in preds]) * 100))
                print("f1Score: %.2f%%" % (f1_score(labels, [custom_round(s) for s in preds]) * 100))
                print("confMatrix: \n" + str(confusion_matrix(labels, [custom_round(s) for s in preds])))
            return preds, labels

    def __callbacks(self, filepath_name):
        """
        Define NN callbacks during training
        """
        return [
            EarlyStopping(monitor='val_mean_absolute_error', mode='auto', patience=self.parameters.esp(), verbose=1),
            ModelCheckpoint(filepath=filepath_name, monitor='val_mean_absolute_error', verbose=1, save_best_only=True,
                            save_weights_only=False, mode='auto', period=1)
        ]

    def __load_raw_data(self, datafile):
        """
        Load the data from @p datafile
        """
        self.entries = self.p_get_dataset(datafile)
        self.entries = self.__filter_dataset(self.conditions)
        if not self.entries:
            raise Exception("No data to train")
        self.names = [entry[0] for entry in self.entries]
        self.peptides_data = [entry[1] for entry in self.entries]
        self.mhcs_data = [entry[2] for entry in self.entries]
        self.total_size = len(self.names) + self.parameters.num_augmentation()
        self.results = [entry[3] for entry in self.entries]
        if self.verbose:
            print("neg vs pos: " + str(np.bincount([int(x) for x in self.results])))

    def __collect_peptide_data(self):
        """
        Collect the peptide data
        """
        max_len = self.__find_max_peptlength()
        # Create PEPTIDES data
        self.peptides = np.zeros((len(self.peptides_data) + self.parameters.num_augmentation(), max_len, self.parameters.feature_size()))
        # Fill data + augment if requested
        self.peptides = self.__fill_ndarray(self.peptides, [ut.to_blossom(seq, normalize=True) for seq in self.peptides_data] +
                                            [ut.to_blossom(seq, normalize=True) for seq in self.augmented_peptides])
        if self.parameters.physchem():
            self.peptides = self.append_physchem(self.peptides, [ut.to_physchem(seq) for seq in self.peptides_data] +
                                       [ut.to_physchem(seq) for seq in self.augmented_peptides])
        # print("peptides shape: " + str(self.peptides.shape))

    def __collect_mhc_data(self):
        """
        Collect the MHC data
        """
        # Create MHCs data
        self.mhcs = np.zeros((len(self.mhcs_data) + self.parameters.num_augmentation(), len(self.mhcs_data[0]), self.parameters.feature_size()))
        # Fill data + augment if requested
        self.mhcs = self.__fill_ndarray(self.mhcs, [ut.to_blossom(seq, normalize=True) for seq in self.mhcs_data] +
                                        [ut.to_blossom(seq, normalize=True) for seq in self.augmented_mhc])
        if self.parameters.physchem():
            self.mhcs = self.append_physchem(self.mhcs, [ut.to_physchem(seq) for seq in self.mhcs_data] +
                                   [ut.to_physchem(seq, normalize=True) for seq in self.augmented_mhc])
        # print("MHCs shape: " + str(self.mhcs.shape))

    def __collect_label_data(self):
        """
        Collect the labels
        """
        # Get label data
        self.binding_result = np.zeros((len(self.mhcs_data) + self.parameters.num_augmentation(), 1))
        for i in range(len(self.results)):
            self.binding_result[i] = self.results[i]

        # Augment data if needed
        self.train_indices += self.augmented_indices

    def __determine_used_size(self):
        """
        Determine the number of elements used
        """
        preset_size = self.parameters.size()
        # If only an image is to be generated, load only 100 elements
        if self.image_only:
            return 100
        elif preset_size == 'max':
            return self.linecount() - 1
        else:
            return preset_size

    def linecount(self):
        """
        Count the number of element in parameters.datafile
        """
        lc = 0
        with open(self.parameters.datafile()) as f:
            for i, _ in enumerate(f):
                lc += 1
        return lc

    def __find_max_peptlength(self):
        """
        Find max peptidelength used
        """
        # Find max peptide length
        max_len = 0
        for peptide in self.peptides_data:
            length = len(peptide)
            if length > max_len:
                max_len = length
                if self.verbose:
                    print("max pept length = " + str(length))
                    print("peptide = " + str(peptide))
        return max_len

    def __filter_dataset(self, conditions):
        """
        Apply @p conditions to the dataset
        """
        return [el for el in self.entries if self.__validate_conditions(el, conditions)]

    def __validate_conditions(self, entry, conditions):
        """
        Validate @p conditions for @p entry
        """
        for condition in conditions:
            if not condition(entry):
                return False
        return True

    def __fill_ndarray(self, ndarray, data_array, size_limit=0):
        """
        Fill @p ndarray with element from @p data_array
        """
        for i in range(len(data_array)):
            if i > size_limit and size_limit > 0:
                return ndarray
            for ii in range(len(data_array[i])):
                for iii in range(len(data_array[i][ii])):
                    ndarray[i][ii][iii] = data_array[i][ii][iii]

        return ndarray

    def __augment_sets(self):
        """
        Augment the datasets if parameters.num_augmentation() > 0
        """
        # print("len(train_indices): " + str(len(self.train_indices)))
        augmented_start_index = 0

        if len(self.train_indices) != 0:
            augmented_start_index = self.train_indices[-1] + 1

        # print("augmented start index: " + str(augmented_start_index))
        (augmented_mhc, augmented_peptides, augmented_binding_results, augmented_indices) = ([], [], [], [])
        for index in range(self.parameters.num_augmentation()):
            new_entry = self.augmenter.generate()

            # Add mhc sequence to mhc_set
            augmented_mhc.append(new_entry['mhc'])

            # Add peptide sequence to peptide_set
            augmented_peptides.append(new_entry['peptide'])

            # Add binding_result to binding_result_set
            augmented_binding_results.append(new_entry['result'])

            # Add new index
            augmented_indices.append(augmented_start_index + index)

        return augmented_mhc, augmented_peptides, augmented_binding_results, augmented_indices

    def p_get_dataset(self, datafile, col=None):
        """
        Retrieve the dataset
        """
        if datafile == '':
            datafile = self.parameters.datafile()
        data = np.loadtxt(datafile,
                          dtype='str',
                          skiprows=1,
                          delimiter=';',
                          usecols=col
                          )[:self.used_size]
        return data

    def generate_schema(self, filename=''):
        """
        Generate image from NN
        """
        try:
            if not filename:
                filename = ut.current_time() + '_' + self.type()
            plot_model(self.model, show_shapes=True, to_file='model_images/' + filename + '.png')
        except AttributeError:
            print("Can't generate schema, model not initialized.")


def custom_round(val, threshold=0.5):
    """
    Custom round function for values between 0.0 and 1.0
    """
    return int(val >= threshold)
