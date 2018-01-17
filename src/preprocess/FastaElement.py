class FastaElement:

    def __init__(self):
        self.reset()

    def reset(self):
        """
        Reset all fields
        """
        self._sequence = ""
        self._name = ""
        self._code = ""
        self._length = 0
        self._valid = False

    def __eq__(self, other):
        if isinstance(other, FastaElement) \
                and self.sequence() == other.sequence() \
                and self.full_name() == other.full_name() \
                and self.code() == other.code():
            return True
        return False

    def __hash__(self):
        return hash(self._sequence)

    def valid(self, verbose=False):
        log = verbose
        if self._name == "":
            if log:
                print("Bad name")
            return False
        if self._code == "":
            if log:
                print("Bad code")
            return False
        if self._length == 0:
            if log:
                print("Bad length")
            return False
        if self._sequence == "":
            if log:
                print("Bad sequence")
            return False
        if log:
            print("self._valid = " + str(self._valid))
            print("Self.name = " + str(self._name))
        return self._valid

    def full_name(self):
        """
        return HLA:HLA00001
        """
        return self._name

    def compact_name(self):
        """
        return HLA-A*01:01
        """
        # Extract mhc type and append '-'
        name = self.full_name().split(":")[0] + "-"
        name += ":".join(self.code().split(":")[:2])
        return name

    def sequence(self):
        return self._sequence

    def size(self):
        return self._length

    def code(self):
        return self._code

    def to_file(self):
        return self.full_name() + "," + self.compact_name() + "," + self.sequence()

    def set_header(self, header_string):
        assert(header_string[0] == '>')
        self.parse_header(header_string)

    def set_sequence(self, sequence_string):
        assert('>' not in sequence_string)
        self._sequence = sequence_string

    def append_sequence(self, subsequence):
        assert('>' not in subsequence)
        self._sequence += subsequence

    def parse_header(self, header_sequence):
        parts = header_sequence.split(" ")
        assert(parts[3].strip() == 'bp')
        self._name = parts[0][1:].strip()
        self._code = parts[1].strip()
        self._length = int(parts[2].strip())
        # Any sequence of the form A*00:00:00N is invalid
        # Sequences of the form A*00:00:00 are valid
        self._valid = (self._code[-1] != 'N')

    def __str__(self):
        return (self.compact_name() + "\n" + self.sequence())
if __name__ == '__main__':
    h1 = ">HLA:HLA09458 A*24:232N 53 bp"
    s1 = "-------------------------SHSMRYFSTSVSRPGPAAGSPASSPWATW - -----" \
        "--T - TRSSCGSTATPRARGW - --------SRGRRGX - -----------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "----------------------------------------------------"

    h2 = ">HLA:HLA06816 B*51:118N 112 bp"
    s2 = "-------------------------SHSMRYFYTAMSRPGRGEPRFIAVGYVDD - ----T" \
        "QFVRFDSDAASPRTEPRAPWIEQEGPEYWDRNTDLQ - ------DQHTDLPREPADRA - --" \
        "--PLLQPERGR - -VSHLADDVWLRRGAGRAPPPRAX------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "------------------------------------------------------------" \
        "----------------------------------------------------"

    header = ">HLA:HLA14404 A*26:118:423 337 bp"
    sequence = "MAVMAPRTLVLLLSGALALTQTWAGSHSMRYFYTSVSRPGRGEPRFIAVGYVDHTQFVRF" \
        "DSDAASQRMEPRAPWIEQEGPEYWDRNTRNVKAHSQTDRANLGTLRGYYNQSEDGSHTIQ" \
        "RMYGCDVGPDGRFLRGYQQNAYDGKDYIALNEDLRSWTAADMAAQITQRKWETAHEAEQW" \
        "RAYLEGRCVEWLRRYLENGKETLQRTDAPKTHMTHHAVSDHEATLRCWALSFYPAEITLT" \
        "WQRDGEDQTQDTELVETRPAGDGTFQKWASVVVPSGQEQRYTCHVQHEGLPKPLTLRWEP" \
        "SSQPTIPIVGIIAGLVLFGAVIAGAVVAAVMWRRKSS"
    header2 = ">HLA:HLA14404 A*26:118N 337 bp"
    sequence2 = "MAVMAPRTLVLLLSGALALTQTWAGSHSMRYFYTSVSRPGRGEPRFIAVGYVDHTQFVRF" \
        "DSDAASQRMEPRAPWIEQEGPEYWDRNTRNVKAHSQTDRANLGTLRGYYNQSEDGSHTIQ" \
        "RMYGCDVGPDGRFLRGYQQNAYDGKDYIALNEDLRSWTAADMAAQITQRKWETAHEAEQW" \
        "RAYLEGRCVEWLRRYLENGKETLQRTDAPKTHMTHHAVSDHEATLRCWALSFYPAEITLT" \
        "WQRDGEDQTQDTELVETRPAGDGTFQKWASVVVPSGQEQRYTCHVQHEGLPKPLTLRWEP" \
        "SSQPTIPIVGIIAGLVLFGAVIAGAVVADWERKRKSS"

    fl = FastaElement()
    fl.set_header(header)
    fl.set_sequence(sequence)

    assert(fl.full_name() == 'HLA:HLA14404')
    assert(fl.size() == 337)
    assert(fl.code() == 'A*26:118:423')
    assert(fl.compact_name() == 'HLA-A*26:118')
    assert(fl.valid())

    fl2 = FastaElement()
    fl2.set_header(header2)
    fl2.set_sequence(sequence2)
    # assert(not fl2.valid())

    assert(len(list(set([fl, fl2]))) == 2)

    f1 = FastaElement()
    f1.set_header(h1)
    f1.set_sequence(s1)

    f2 = FastaElement()
    f2.set_header(h2)
    f2.set_sequence(s2)
    assert(len(list(set([f1, f2]))) == 2)
