import train.NeuralNetwork as nn

from train.Models import model_wide_shallow as nn_model


class WideShallow(nn.NeuralNetwork):
    def p_load_model(self, oa):
        return nn_model(self.peptides[0].shape,
                        self.mhcs[0].shape,
                        output_activation=oa)

    def p_model_filename(self):
        return "wideshallow.model.hdf5"

    def type(self):
        return 'WideShallow'

