import json
import os.path

from preprocess.FastaElement import FastaElement

import preprocess.Util as ut

MHCJSON_PATH = 'src/preprocess/mhc.json'

class MhcParser:
    def __init__(self, fc='SHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDA', lc='RAYLEGTCVEWLRRYLENG'):
        self.first_const = fc
        self.last_const = lc
        self.data = []
        self.shortest = 1000
        self.largest = 0
        self.reset()

    def reset(self):
        self.data = []
        self.shortest = 1000
        self.largest = 0


    def __finish_element(self, element):
        """
        Finalize the processing of @p element.
        """
        if element.valid():
            if element.size() < self.shortest:
                self.shortest = element.size()
            if element.size() > self.largest:
                self.largest = element.size()

            element.set_sequence(self.clean_sequence(element.sequence()))
            self.data.append(element.to_file())
            return element

        return None

    def read_hla_file(self, filename):
        """
        Read the FASTA file containing HLA proteins
        """
        with open(filename) as f:
            fe = FastaElement()
            for line in f.readlines():
                # Start reading new element
                if line[0] == '>':
                    if self.size() > 0 and self.size() % 50 == 0:
                        # return  # Quit after 50 elements
                        print("Read " + str(self.size()) + " sequences")

                    # Clean and store element
                    self.__finish_element(fe)

                    fe.reset()
                    fe.set_header(line)
                else:
                    # Append sequence, skip final newline
                    fe.append_sequence(line[:-1])

            # Clean and store element
            self.__finish_element(fe)

    def clean_sequence(self, sequence):
        """
        Clean @p sequence
        """
        return self.__remove_after_last_const(self.__remove_upto_first_const(sequence))

    def __remove_upto_first_const(self, sequence):
        """
        Remove aminoacids from sequence upto the first constant substring (first delimiter)
        """
        (target, index) = ut.find_closest_match(sequence, self.first_const)
        remainder = sequence[index + len(target):]
        return remainder

    def __remove_after_last_const(self, sequence):
        """
        Remove aminoacids from sequence after second constant substring (second delimiter)
        :param sequence:
        :return:
        """
        (target, index) = ut.find_closest_match(sequence, self.last_const)
        remainder = sequence[:index]
        return remainder

    def write_file(self, filename):
        """
        Write formatted data to @p filename
        """
        with open(filename, "w") as f:
            line = ""
            for key in ["HLA name", "Allele Name", "Sequence"]:
                line += key + ","
            f.write(line[:-1] + "\n")

            for entry in self.data:
                f.write(entry + "\n")

    def remove_duplicates(self):
        """
        Remove duplicate data
        """
        self.data = list(set(self.data))

    def __load_dictionary(self, filepath):
        """
        If the dictionary exists, load it.
        """
        with open(filepath) as mhcjson:
            dictionary = json.loads(mhcjson.readlines()[0])
            # print("Loaded MHC dictionary from '" + filepath + "'")
            return dictionary

    def __store_dictionary(self, filepath, dictionary):
        """
        Store the dictionary
        """
        with open(filepath, 'w') as mhcjson:
            mhcjson.write(json.dumps(dictionary))
            mhcjson.close()
        # print("Newly created MHC dictionary at '" + filepath + "'")

    def get_dictionary(self, overwrite=False):
        """
        Get dictionary: {HLA name: HLA sequence}
        Only load if no dictionary exists.
        """
        dictionary = {}
        if (not overwrite) and os.path.isfile(MHCJSON_PATH):
            dictionary = self.__load_dictionary(MHCJSON_PATH)
        else:
            # If it doesn't exist, create, store and return it.
            for entry in self.data:
                (_, name, sequence) = tuple(entry.split(","))
                if name in dictionary:
                    if dictionary[name] != sequence:
                        print(name + " already exists with other sequence: " + str(dictionary[name]))
                        print("current sequence: " + str(sequence))
                dictionary[name] = sequence
            self.__store_dictionary(MHCJSON_PATH, dictionary)
        return dictionary

    def size(self):
        """
        Get the size of the data
        """
        return len(self.data)


def main():
    FIRST_CONST = 'SHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDA'
    LAST_CONST = 'RAYLEGTCVEWLRRYLENG'

    path_prefix = "data/original/"

    hla_files = [
        "A_prot.fasta",
        "B_prot.fasta",
        "C_prot.fasta"
    ]

    parser = MhcParser(FIRST_CONST, LAST_CONST)

    # Initialize data list
    parser.reset()
    for file in hla_files:
        # Read and parse data
        parser.read_hla_file(path_prefix + file)

    parser.remove_duplicates()

    # Write data to new file
    parser.write_file("MhcData.csv")

    print("Generated " + str(parser.size()) + " new MHC entries.")
