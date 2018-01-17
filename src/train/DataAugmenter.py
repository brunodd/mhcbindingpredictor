import Bio.SeqUtils as bs
import random

random.seed(1)

class DataAugmenter:
    def __init__(self, peptide_length=None):
        # Store all possible aminoacid symbols
        self.aa = bs.IUPACData.protein_letters
        self.peptide_length = peptide_length


    def __generate_aminoacid(self):
        """
        Generate random aminoacid
        """
        return random.choice(self.aa)

    def __generate_mhc(self):
        """
        Generate MHC sequence of length 116
        """
        return self.__generate_sequence(116)

    
    def __generate_peptide(self):
        """
        Generate peptide sequence of length 9 or 10
        """
        if not self.peptide_length:
            return self.__generate_sequence(random.choice([9,10]))
        return self.__generate_sequence(self.peptide_length)

    def __generate_sequence(self, seq_len):
        """
        Generate random sequence of length @p seq_len
        """
        return ''.join([self.__generate_aminoacid() for _ in range(seq_len)])

    def generate(self):
        """
        Generate random entry

        return {'mhc' : MHC_SEQUENCE,
                'peptide' : PEPTIDE_SEQUENCE,
                'result' : BOOLEAN_RESULT
                }
        """
        return {
                'mhc': self.__generate_mhc(),
                'peptide': self.__generate_peptide(),
                'result': 0
               }

if __name__ == '__main__':
    da = DataAugmenter()
    for _ in range(10):
        print(da.generate())
