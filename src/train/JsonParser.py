import json


class JsonParser:
    def __init__(self, filename):
        self.params = []
        with open(filename) as js_file:
            self.params = json.load(js_file)

    def size(self):
        return self.params['size']

    def epochs(self):
        return self.params['epochs']

    def bs(self):
        return self.params['batchsize']

    def datafile(self):
        return self.params['datafile']

    def output_act(self):
        return self.params['activations']['output']

    def esp(self):
        return self.params['early_stopping_patience']

    def feature_size(self):
        return self.aa_size() + self.params['physchem_size']

    def aa_size(self):
        return self.params['aa_size']

    def num_augmentation(self):
        return self.params['augment_size']

    def physchem(self):
        return self.params['physchem_size'] > 0
