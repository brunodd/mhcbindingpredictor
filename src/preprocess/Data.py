import argparse

import numpy as np

DATA_KEYS = ["Allele", "Peptide", "MHC", "Binding result"]

def write_all_2_file(filename, data, size_limit=10):
    """
    Write all formatted data to @p filename
    """
    with open(filename, "w") as f:
        line = ""
        for key in DATA_KEYS:
            line += key + ";"
        f.write(line[:-1] + "\n")

        linecount = 0
        for entry in data:
            if linecount >= size_limit and size_limit > 0:
                return
            linecount += 1
            f.write(str(entry[0]) + ";" + str(entry[1]) + ";" + str(entry[2]) \
                    + ";" + str(entry[3]) + "\n")


def write_peplen9_2_file(filename, data, size_limit=10):
    """
    Write formatted data with peptide length 9 to @p filename
    """
    with open(filename, "w") as f:
        line = ""
        for key in DATA_KEYS:
            line += key + ";"
        f.write(line[:-1] + "\n")

        linecount = 0
        for entry in data:
            if linecount >= size_limit and size_limit > 0:
                return
            if len(entry[1]) == 9:
                linecount += 1
                f.write(str(entry[0]) + ";" + str(entry[1]) + ";" + str(entry[2]) \
                        + ";" + str(entry[3]) + "\n")


def write_peplen10_2_file(filename, data, size_limit=10):
    """
    Write formatted data with peptide length 9 to @p filename
    """
    with open(filename, "w") as f:
        line = ""
        for key in DATA_KEYS:
            line += key + ";"
        f.write(line[:-1] + "\n")

        linecount = 0
        for entry in data:
            if linecount >= size_limit and size_limit > 0:
                return
            if len(entry[1]) == 10:
                linecount += 1
                f.write(str(entry[0]) + ";" + str(entry[1]) + ";" + str(entry[2]) \
                        + ";" + str(entry[3]) + "\n")
    print("Wrote " + str(linecount) + " entries to " + filename)


def write_A_2_file(filename, data, size_limit=10):
    """
    Write formatted data with peptide length 9 to @p filename
    """
    with open(filename, "w") as f:
        line = ""
        for key in DATA_KEYS:
            line += key + ";"
        f.write(line[:-1] + "\n")

        linecount = 0
        for entry in data:
            if linecount >= size_limit and size_limit > 0:
                return
            if 'HLA-A' in entry[0]:
                linecount += 1
                f.write(str(entry[0]) + ";" + str(entry[1]) + ";" + str(entry[2]) \
                        + ";" + str(entry[3]) + "\n")
    print("Wrote " + str(linecount) + " entries to " + filename)


def write_B_2_file(filename, data, size_limit=10):
    """
    Write formatted data with peptide length 9 to @p filename
    """
    with open(filename, "w") as f:
        line = ""
        for key in DATA_KEYS:
            line += key + ";"
        f.write(line[:-1] + "\n")

        linecount = 0
        for entry in data:
            if linecount >= size_limit and size_limit > 0:
                return
            if 'HLA-B' in entry[0]:
                linecount += 1
                f.write(str(entry[0]) + ";" + str(entry[1]) + ";" + str(entry[2]) \
                        + ";" + str(entry[3]) + "\n")
    print("Wrote " + str(linecount) + " entries to " + filename)


def write_C_2_file(filename, data, size_limit=10):
    """
    Write formatted data with peptide length 9 to @p filename
    """
    with open(filename, "w") as f:
        line = ""
        for key in DATA_KEYS:
            line += key + ";"
        f.write(line[:-1] + "\n")

        linecount = 0
        for entry in data:
            if linecount >= size_limit and size_limit > 0:
                return
            if 'HLA-C' in entry[0]:
                linecount += 1
                f.write(str(entry[0]) + ";" + str(entry[1]) + ";" + str(entry[2]) \
                        + ";" + str(entry[3]) + "\n")
    print("Wrote " + str(linecount) + " entries to " + filename)


def read_data(file):
    """
    Read data from file
    """
    data = np.genfromtxt(file,
                         dtype='str',
                         skip_header=True,
                         delimiter=',')

    return data


def mhc_data(raw_data):
    """
    Parse raw MHC data.

    returns tuple of lists: (mhc_name_list, mhc_sequence_list)
    """
    print("Checking MHC data")
    mhcs = []
    data = []
    temp_data = []
    for d in raw_data:
        temp_data.append(list(d))

    raw_data = temp_data

    for d in raw_data:
        # entry = list2string(extract_numeric_sequence(d, index=2))
        entry = d[2]
        # Allele name = index 1
        mhcs.append(d[1])
        data.append(entry)

    print("Size of data = " + str(len(raw_data)))
    return (mhcs, data)


def peptide_data(raw_data):
    """
    Parse raw peptide data.

    returns tuple of lists: (mhc_name_list, peptide_sequence_list, boolean_bindings_list)
    """
    print("Checking Peptide data")
    mhc_names = []
    data = []
    bindings = []
    temp_data = []
    for d in raw_data:
        temp_data.append(list(d))

    raw_data = temp_data
    for d in raw_data:
        entry = d[0]
        # Only keep peptides of length 9 or 10
        if filter_size(entry, 9) or filter_size(entry, 10):
            entry = d[0]
            data.append(entry)
            # Allele name = index 1
            mhc_names.append(d[2])
            # Append if binding is True or False
            bindings.append(int(d[1] == 'True'))

    return (mhc_names, data, bindings)


def filter_size(sequence, size):
    """
    Only keep sequence of length @p size.
    """
    if len(sequence) == size:
        return sequence
    return None


def get_valid_mhc_names(peptides_tuple, mhcs_tuple):
    """
    Get list of MHC names that are also used in experiments with peptides.
    """
    # Extract MHC names from peptides and mhcs data.
    mhcs_peptides = peptides_tuple[0]
    mhcs_mhcs = mhcs_tuple[0]

    inter_set = set(mhcs_peptides).intersection(set(mhcs_mhcs))

    return list(inter_set)


def merge_valid_data(peptides_tuple, mhcs_tuple, valid_mhcs):
    """
    Merge the (valid) mhcs with the peptide experiments.
    """
    # Get the indices of the peptide experiments that are done with valid MHCs.
    peptides_indices = [index for index, mhc in enumerate(peptides_tuple[0]) if mhc in valid_mhcs]

    result = []
    for index in peptides_indices:
        # Create data entry and append to the results
        mhc_name = peptides_tuple[0][index]   # MHC name
        full_entry = [
            mhc_name,   # MHC name
            peptides_tuple[1][index],   # Peptide sequence
            mhcs_tuple[1][mhcs_tuple[0].index(mhc_name)],   # MHC sequence
            peptides_tuple[2][index]    # Binding result
        ]

        result.append(full_entry)

    return result


def main(mhc_file, peptide_file):
    """
    Main function
    """
    # Collect and clean MHC data
    raw_mhc_data = read_data(mhc_file)
    mhcs_tuple = mhc_data(raw_mhc_data)

    # Collect and clean experimental peptide data
    raw_binding_data = read_data(peptide_file)
    peptides_tuple = peptide_data(raw_binding_data)

    # Drop all useless data
    useful_mhcs = get_valid_mhc_names(peptides_tuple, mhcs_tuple)

    # Merge the remaining (valid) data
    training_data = merge_valid_data(peptides_tuple, mhcs_tuple, useful_mhcs)

    # Write the data to a file
    # Size limit set to 500k = larger than total data = process all
    write_all_2_file("TrainingDataAll.csv", training_data, 500000)
    # write_peplen10_2_file("Peplen10TrialFile.csv", training_data, 500000)
    # write_peplen9_2_file("Peplen9TrialFile.csv", training_data, 500000)
    # write_A_2_file("ATrialFile.csv", training_data, 500000)
    # write_B_2_file("BTrialFile.csv", training_data, 500000)
    # write_C_2_file("CTrialFile.csv", training_data, 500000)

    print(str(len(training_data)) + " usefull entries found available")
