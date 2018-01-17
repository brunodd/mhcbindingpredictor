import random
random.seed(1)


class Divider:
    def __init__(self, data, custom_division, verbose=False):
        self.__data = data
        self.__filter_type = 'random'
        self.__filter = None
        self.__train_percentage = custom_division[0]
        self.__test_percentage = custom_division[1]
        self.__validation_percentage = custom_division[2]
        self.verbose = verbose

        assert(self.__train_percentage + self.__test_percentage + self.__validation_percentage == 1.0)
        assert(self.__test_percentage >= 0.0
               and self.__validation_percentage >= 0.0
               and self.__train_percentage >= 0.0)


    def __family(self, full_name):
        """
        Get the 'family' name of a sequence

        e.g. HLA-A*01:02 yields HLA-A*01
        """
        return full_name.split(':')[0]

    def families(self):
        """
        Get list of all (unique) family names in @p full_names
        """
        p_families = self.__extract_families(self.__data)
        return p_families

    def __extract_families(self, full_names):
        """
        Get list of all (unique) family names in @p full_names
        """
        family_names = []
        for name in full_names:
            family_names.append(self.__family(name))

        # Add sorting for reproducibility
        return sorted(list(set(family_names)))

    def __type(self, full_name):
        """
        Get the style of a sequence

        e.g. HLA-A*01:02 yields HLA-A
        """
        return full_name.split('*')[0]

    def proteins(self):
        """
        Get list of all proteins
        """
        # Add sorting for reproducibility
        return sorted(list(set(self.__data)))

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
        self.__filter_type = type
        self.__filter = filter

    def size(self):
        """
        Get size of dataset
        """
        return len(self.__data)

    def min_size_training_data(self):
        """
        Get minimum size of the training set based on the train_percentage
        """
        return round(self.size() * self.__train_percentage)

    def validate(self):
        """
        True if validation set is requested
        """
        return self.__validation_percentage > 0

    def __create_all_sets(self, train, test):
        """
        Create train, test and validation set
        """
        if self.validate():
            train_set, val_set = self.__split_train_in_train_val(train)
            return train_set, test, val_set
        return train, test, []

    def __random_families(self):
        """
        Randomly divide the families in train, test, validation sets
        """
        p_families = self.families()
        random.shuffle(p_families)
        train, test = self.extract_test(p_families)
        return self.__create_all_sets(train, test)

    def __split_train_in_train_val(self, train_data):
        """
        Split the train set into train and validation set
        """
        val_size = int(len(train_data) * self.__validation_percentage)
        train, val = train_data[:-val_size], train_data[-val_size:]
        return train, val

    def random_families(self):
        """
        Get train, test and validation sets with random families
        """
        while True:
            test_indices, train_indices, val_indices = [], [], []
            # Generate random family division
            train, test, validation = self.__random_families()

            # Divide indices in relevant datasets
            for index, name in enumerate(self.__data):
                if self.__family(name) in test:
                    test_indices.append(index)
                elif self.__family(name) in validation:
                    val_indices.append(index)
                else:
                    train_indices.append(index)

            # Verify if training data is large enough
            if self.min_size_training_data() <= len(train_indices):
                if self.verbose:
                    print("Test set: " + str(test))
                    if val_indices:
                        print("Validation set; " + str(validation))
                return train_indices, test_indices, val_indices
            elif self.verbose:
                print("min_training_size: " + str(self.min_size_training_data()))
                print("deprecated_train_01 indices: " + str(len(train_indices)))
                print("test indices; " + str(len(test_indices)))
                print("val indices: " + str(len(val_indices)))


    def __random_proteins(self):
        """
        Divide proteins randomly in train, test and validation set
        """
        p_proteins = self.proteins()
        random.shuffle(p_proteins)
        train, test = self.extract_test(p_proteins)
        return self.__create_all_sets(train, test)

    def random_proteins(self):
        """
        Get train, test and validation set with randomly divided proteins
        """
        count = 0
        while True:
            test_indices, train_indices, val_indices = [], [], []
            count += 1
            train, test, validation = self.__random_proteins()
            for index, name in enumerate(self.__data):
                if name in test:
                    test_indices.append(index)
                elif name in validation:
                    val_indices.append(index)
                else:
                    train_indices.append(index)

            if self.min_size_training_data() <= len(train_indices) or count >= 100:
                if self.verbose:
                    print("Test set: " + str(test))
                    if val_indices:
                        print("Validation set; " + str(validation))
                return train_indices, test_indices, val_indices
            elif self.verbose:
                print("min_training_size: " + str(self.min_size_training_data()))
                print("deprecated_train_01n indices: " + str(len(train_indices)))
                print("test indices; " + str(len(test_indices)))
                print("val indices: " + str(len(val_indices)))

    def extract_test(self, mylist):
        """
        Divide data in @p mylist into train and test data according to given test_percentage.
        """
        if self.__test_percentage == 0.0:
            return mylist, []
        elif self.__test_percentage == 1.0:
            return [], mylist
        else:
            test_size = round(len(mylist) * self.__test_percentage)
            return mylist[:-test_size], mylist[-test_size:]

    def get_indices(self):
        """
        Return tuple of (train_indices, test_indices, val_indices)
        """
        train_indices, test_indices, val_indices = [], [], []
        if self.__filter_type == 'random':
            if self.verbose:
                print('type = random')
            train_indices, test_indices, val_indices = self.random_proteins()

        if self.__filter_type == 'list':
            if self.verbose:
                print('type = list')
                print('filter =', self.__filter)
            for index, name in enumerate(self.__data):
                if name in self.__filter:
                    test_indices.append(index)
                else:
                    train_indices.append(index)
        elif self.__filter_type == 'family':
            if self.verbose:
                print('type = family')
                print('filter =', self.__filter)
            if self.__filter:
                for index, name in enumerate(self.__data):
                    if self.__family(name) in self.__filter:
                        test_indices.append(index)
                    else:
                        train_indices.append(index)
            else:
                return self.random_families()

        elif self.__filter_type == 'type':
            if self.verbose:
                print("type = type")
                print('filter = ', self.__filter)
            for index, name in enumerate(self.__data):
                if self.__type(name) in self.__filter:
                    test_indices.append(index)
                else:
                    train_indices.append(index)
        train_indices, val_indices = self.__split_train_in_train_val(train_indices)
        if self.verbose:
            print('train', len(train_indices))
            print('val', len(val_indices))
            print('test', len(test_indices))
        return train_indices, test_indices, val_indices


    def test(self):
        assert(self.__type('HLA-A*02:01') == 'HLA-A')
        assert(self.__family('HLA-A*02:01') == 'HLA-A*02')
        assert('HLA-A*02' in self.__extract_families(self.__data))
        self.set_test_type(type='list', filter=['HLA-A*02:02'])
        train, test = self.get_indices()
        assert(len(test) == 1473)
        self.set_test_type(type='family', filter=['HLA-A*80'])
        train, test = self.get_indices()
        assert(len(test) == 48)
        self.set_test_type(type='type', filter=['HLA-B'])
        train, test = self.get_indices()
        assert(len(test) == 16890)
        self.set_test_type(type='type', filter=['HLA-C'])
        train, test = self.get_indices()
        assert(len(test) == 0)
        return True

def main():
    import numpy as np
    np.random.seed(1)

    # DATA_FILE_NAME = 'TrainingData50k.csv'
    DATA_FILE_NAME = 'TrainingDataAll.csv'
    names = np.loadtxt(DATA_FILE_NAME,
                      dtype='str',
                      skiprows=1,
                      delimiter=';',
                      usecols=0)
    # # print(names[:10])
    # div = Divider(names, custom_division=(0.8, 0.1, 0.1), seed=None)
    # # print(str(div.random_proteins()[1][:10]))
    # # deprecated_train_01, test, val = div.random_families()
    # train, test, val = div.random_families()
    # print("Train: " + str(len(train)))
    # print("Test: " + str(len(test)))
    # print("Validation: " + str(len(val)))
    # div = Divider(names, custom_division=(0.8, 0.1, 0.1), seed=None)
    # # print(str(div.random_proteins()[1][:10]))
    # # deprecated_train_01, test, val = div.random_families()
    # train, test, val = div.random_families()
    # print("Train: " + str(len(train)))
    # print("Test: " + str(len(test)))
    # print("Validation: " + str(len(val)))
    #
    # # print("test: %.2f%%" % + (len(div.random_proteins()[1]) / div.size() * 100) + "%")
    # # print(div.families()[:10])
    # # # assert(div.test())
    # #
    # # print(div.families())
    # # print(str(len(div.families())))
    # # print("total size: " + str(div.size()))
    # # deprecated_train_01, test = div.random_families()
    # # # print("deprecated_train_01: " + str(deprecated_train_01))
    # # # print("test: " + str(test))
    div = Divider(names, custom_division=(0.7, 0.1, 0.2), seed=1)
    # div.set_test_type(type='family', filter=['HLA-A*01', 'HLA-B*51'])
    # div.get_indices()
    # print(" ")
    # div.set_test_type(type='random')
    # div.get_indices()
    # print(" ")
    print('num names:', len(names))
    # div.set_test_type(type='type', filter=['HLA-C'])
    train, test, val = div.get_indices()
    print("train:", names[train])
    print("\tsize:", len(names[train]))
    print("test:", names[test])
    print("\tsize:", len(names[test]))
    print("val:", names[val])
    print("\tsize:", len(names[val]))
    # print(" ")
    # div.set_test_type(type='list', filter=['HLA-A*01:01', 'HLA-B*51:01'])
    # div.get_indices()

    div.set_test_type(type='family')
    train, test, val = div.get_indices()
    print('num indices: train: {:d}\ttest: {:d}\tval: {:d}'.format(len(train), len(test), len(val)))

