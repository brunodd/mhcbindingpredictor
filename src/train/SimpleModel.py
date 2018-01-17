import train.NeuralNetwork as nn

from train.Models import simple_model as nn_model


class SimpleModel(nn.NeuralNetwork):
    def p_load_model(self, oa):
        return nn_model(self.peptides[0].shape,
                        self.mhcs[0].shape,
                        output_activation=oa)

    def p_model_filename(self):
        return "simple.model.hdf5"

    def type(self):
        return 'SimpleModel'
