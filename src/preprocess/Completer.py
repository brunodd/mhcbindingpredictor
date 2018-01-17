import preprocess.Util as ut
from preprocess.MHCCleaner import MhcParser as mp

class Completer:
    """
    Complete given data entry for model
    """
    def __init__(self):
        self.parser = mp()
        self.dictionary = self.parser.get_dictionary()

    def complete(self, entry):
        """
        Complete entry of type: (mhc name, peptide sequence)
        :return: (mhc name, mhc sequence, peptide sequence)

        :precondition: mhc name has to have structure: HLA-A*01:01
        """
        try:
            assert(len(entry) == 2)
        except AssertionError as err:
            print("entry failed assertion: ", err)
            print("entry: ", entry)
        try:
            name = ut.reformat(entry[0])
            return name, self.dictionary[name], entry[1]
        except KeyError as err:
            print("Check spelling or update the MHC dictionary.")
            print("Failed to retrieve sequence of ", err)


if __name__ == '__main__':
    c = Completer()

    print(c.complete(['HLA-A*01:01', 'MASTTPITM']))
    print(c.complete(['HLA-A*02:01', 'MASTTPITM']))
    print(c.complete(['HLA-A02:01', 'MASTTPITM']))
