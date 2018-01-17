import Bio.SeqUtils as bs
import Bio.SubsMat.MatrixInfo as bm
from pyteomics import mass, electrochem
import pyteomics
from datetime import datetime as dt
import copy



def get_min_max(blossum_type=90):
    """
    Get min and max values in BLOSUM matrix for normalization\
    """
    if blossum_type == 90:
        return min(bm.blosum90.values()), max(bm.blosum90.values())


def normalize_blossum_scores(blossum_vector, min_value, max_value):
    """
    Normalize values in @p blossum_vector using @p min_value and @p max_value
    """
    normalized = []
    for el in blossum_vector:
        assert(max_value >= el >= min_value)
        normalized.append(normalize(el, min_value, max_value))
    return normalized


def normalize(x, mi, ma, a=0, b=1):
    """
    Normalize x into range [a,b] with min and max data values @p min and @p max resp.
    return ((b-a)(x-min) / (max-min)) + a
    """
    return ((b-a)*(x-mi) / (ma-mi)) + a


def to_numeric(sequence):
    """
    Convert @p sequence into list of integers

    Each integer corresponds to a canonical AA.
    """
    num_result = []
    for c in sequence:
        num_result.append(aa2int(c))

    return num_result


def to_blossom(sequence, normalize=False):
    """
    Convert @p sequence to a matrix where each AA is represented by a vector of BLOSUM scores.
    """
    min_val, max_val = get_min_max()
    blossom_result = []
    for c in sequence:
        score_vector = aa2blossom(c)
        if normalize:
            score_vector = normalize_blossum_scores(score_vector, min_val, max_val)
        blossom_result.append(score_vector)

    return blossom_result


def aa2blossom(AA):
    """
    Convert @p AA to a vector with BLOSUM scores
    """
    vect = []
    for i in range(20):
        try:
            vect.append(bm.blosum90[AA, bs.IUPACData.protein_letters[i]])
        except KeyError:
            vect.append(bm.blosum90[bs.IUPACData.protein_letters[i], AA])
    return vect


def aa2int(AA):
    """
    Get the index of @p AA
    """
    try:
        if AA == 'X':
            return 21
        return bs.IUPACData.protein_letters.index(AA) + 1

    except:
        print("Failed getting value of " + AA)
        raise ValueError


def int2aa(integer):
    """
    AA corresponding to the index @p integer
    """
    return bs.IUPACData.protein_letters[integer - 1]


def str_rep(integer):
    """
    Get formatted string representation of int @p integer
    """
    return "{:02d}".format(integer)


def current_time():
    """
    Get the current time in string format
    """
    cur_time = dt.now()
    return str_rep(cur_time.month) + str_rep(cur_time.day) \
           + "_" + str_rep(cur_time.hour) + str_rep(cur_time.minute) + str_rep(cur_time.second)


def __read_physchem_line(line):
    """
    Parse physchem properties on @p line.
    """
    return tuple(line[:-1].split(' '))


def read_physchem(filename):
    """
    Read physchem properties from @p filename
    """
    # @p filename contains data for 23 AAs.
    num_aa = 23
    properties = {}
    min_max_values = []

    with open(filename) as f:
        for index, line in enumerate(f.readlines()):
            # Skip physchem property name
            if index % (num_aa + 1) == 0:
                # Append (min,max) tuple to min_max_values array
                min_max_values.append((1000,-1000))
                continue
            # Read and convert data
            (aa, value) = __read_physchem_line(line)
            (aa, value) = (aa, float(value))

            # Update min and max values encountered so far
            if value < min_max_values[-1][0]:
                min_max_values[-1] = (value, min_max_values[-1][1])
            elif value > min_max_values[-1][1]:
                min_max_values[-1] = (min_max_values[-1][0], value)

            if not aa in properties:
                properties[aa] = []
            properties[aa].append(float(value))

    # Normalize values of physchem
    for i, aa in enumerate(properties.keys()):
        props = properties[aa]
        for ii, val in enumerate(props):
            mi = min_max_values[ii][0]
            ma = min_max_values[ii][1]
            properties[aa][ii] = normalize(val, mi, ma)

    return properties
PHYSCHEM_PROPERTIES = read_physchem('physchem.txt')


def aa2physchem(aa):
    """
    Get physchem properties for AA @p aa.
    """
    # copy data or it gets overwritten
    result = copy.deepcopy(PHYSCHEM_PROPERTIES[aa])

    # Collect the mass of aa
    try:
        _mass = mass.calculate_mass(aa)
        result.append(_mass)
    except pyteomics.auxiliary.PyteomicsError:
        result.append(0.0)

    # Collect the pH of aa
    try:
        # Assume neutral pH
        result.append(electrochem.charge(aa, 7))
    except pyteomics.auxiliary.PyteomicsError:
        result.append(0.0)

    return result



def to_physchem(sequence):
    """
    Convert the sequence @p sequence of AAs into a sequence of vectors containing physchem properties of the corresponding AAs.
    """
    physchem_result = []
    for c in sequence:
        physchem_result.append(aa2physchem(c))
    return physchem_result


def main():
    (mi, ma) = get_min_max()
    print((mi, ma))
    assert(aa2int("Y") == 20)
    assert(int2aa(19) == 'W')
    assert(int2aa(aa2int('A')) == "A")
    assert(to_numeric('ACDEFGHIKLMNPQRSTVWY') == list(range(1,21)))

    for i in range(20):
        aa = bs.IUPACData.protein_letters[i]
        print(aa2blossom(aa))
    vec = to_blossom("ACDEFGHIKLMNPQRSTVWY", True)
    print(vec)
    print(current_time())
    print(normalize_blossum_scores([10], mi, ma))
    print("Success")

    print(read_physchem('physchem.txt'))
    print("Y: " + str(aa2physchem('Y')))
    print("YAF" + str(to_physchem("YAF")))
