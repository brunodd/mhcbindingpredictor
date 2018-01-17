from sklearn.metrics import confusion_matrix, precision_score, recall_score, f1_score

from train.WideShallow import WideShallow as WS
from train.Divider import Divider
from preprocess.Completer import Completer

import sendmail as sm


def show_results(preds, test_labels):
    """
    Show performance of @p preds w.r.t. @p test_labels
    """
    prec = precision_score([custom_round(float(s)) for s in test_labels],
                           [custom_round(s, 0.5) for s in preds])
    print("skPrecision: %.2f%%" % (prec * 100))
    print("skRecall: %.2f%%" % (recall_score(test_labels, 
                                             [custom_round(s, 0.5) for s in preds]) * 100))
    print("f1Score: %.2f%%" % (f1_score(test_labels, 
                                        [custom_round(s, 0.5) for s in preds]) * 100))
    print("confMatrix: \n" + str(confusion_matrix(test_labels, 
                                                  [custom_round(s, 0.5) for s in preds])))


def perform_prediction(model, elements):
    """
    Predict the classification of @p elements using (trained) model @p model
    """
    epsilon = 10e-8
    preds = []
    for e in elements:
        # Predict
        pred = model.predict(e)
        pred = min(pred[0], 1.0-epsilon)    # Don't allow 100% probability

        # Convert to 'real-life' prediction
        pred = convert_probability(pred)

        # Append to predictions
        preds.append(pred)
    return preds


def load_patient_data(df):
    """
    Load the patient data from @p df
    """
    # Initialize variables
    c = Completer()
    hlaa, hlab, hlac = [], [], []
    hlaa_labels, hlab_labels, hlac_labels = [], [], []

    with open(df) as data_file:
        for patient_entry in data_file.readlines():
            # Split into list
            entry = patient_entry.split(',')

            # Store label
            label = custom_round(float(entry[2]), threshold=0.5)
            # Remove affinity score
            entry = entry[:2]

            # Complete data
            entry = c.complete(entry)

            if 'HLA-A' in entry[0]:
                hlaa.append(entry)
                hlaa_labels.append(label)
            elif 'HLA-B' in entry[0]:
                hlab.append(entry)
                hlab_labels.append(label)
            elif 'HLA-C' in entry[0]:
                hlac.append(entry)
                hlac_labels.append(label)
            else:
                raise RuntimeWarning("Unrecognized HLA type found:", entry)

    print("Loaded all data. Read {} entries".format(len(hlaa+hlab+hlac)))
    return [(hlaa, hlaa_labels), (hlab, hlab_labels), (hlac, hlac_labels)]


def avg(l):
    """
    Computer average of Boolean list @p l
    """
    return sum(l) / len(l)


def custom_round(val, threshold=0.5):
    """
    Custom round function for values between 0.0 and 1.0 with variable threshold
    """
    return int(val >= threshold)


def compute_average_performance(model, k=5):
    """
    Compute the average
    """
    model.load_data()
    all_preds, all_labels = [], []
    # Train the network @p k times, store prediction results
    for _ in range(k):
        print("Resetting weights")
        model.reinitialize_weights()
        p, l = model.train()
        assert(len(p) == len(l))
        all_preds.append(p)
        all_labels.append(l)

    assert(len(all_preds) == len(all_labels))

    accuracy = []
    precision = []
    recall = []
    f1 = []

    # Compute performance metrics per run
    for i in range(len(all_preds)):
        preds = all_preds[i]
        labels = all_labels[i]

        rounded_preds = [custom_round(s) for s in preds]
        prec = precision_score([custom_round(float(s)) for s in labels], rounded_preds)
        matches = [True if preds[i] == labels[i] else False for i in range(len(preds))]
        accuracy.append(avg(matches) * 100)
        precision.append(prec * 100)
        recall.append(recall_score(labels, rounded_preds) * 100)
        f1.append(f1_score(labels, rounded_preds) * 100)

    # Show average performance
    print('accuracy: %.2f%%' % (avg(accuracy)))
    print("skPrecision: %.2f%%" % (avg(precision)))
    print("skRecall: %.2f%%" % (avg(recall)))
    print("f1Score: %.2f%%" % (avg(f1)))


def convert_probability(p):
    """
    Conversion function for conversion
    from distribution (31% negative, 69% positive)
    to (98.7% negative vs 1.3% positive)
    """
    alpha = 0.013/0.69
    beta = 0.987/0.31
    return (alpha * p) / ((alpha * p) + (beta * (1-p)))


if __name__ == '__main__':
    VERBOSITY = True
    try:
        print("=== Test WS architecture on patient data ===")
        # Create NN
        nn = WS(
            conditions=[
                lambda entry: len(entry[1]) == 9
            ]
        )
        nn.set_verbose(VERBOSITY)

        # TRAIN MODEL
        names = nn.load_data()
        # Divide the data in 70% training data, 10% test data and 20% validation data
        div = Divider(names, (0.7, 0.1, 0.2), verbose=VERBOSITY)
        nn.train(div.get_indices())
        nn.save()
        # nn.load()

        # LOAD PATIENT DATA
        [(hlaa, hlaa_labels), (hlab, hlab_labels), (hlac, hlac_labels)] = load_patient_data('data.txt')
        test_set = hlaa + hlab + hlac
        test_labels = hlaa_labels + hlab_labels + hlac_labels

        # PERFORM PREDICTIONS
        preds = perform_prediction(nn, test_set)

        # SHOW RESULTS
        show_results(preds, test_labels)

        # Notify by mail if finished successfully
        sm.alert('Finished training successfully')

    except Exception as e:
        # Notify by mail if failed
        sm.alert(str(e.args))
        raise e
