import numpy as np

DATA_KEYS = ["Allele", "Peptide", "MHC", "Binding result"]


def extract():
    """
    Extract patient data
    """
    filename = 'data/prior_df.txt'
    with open(filename) as file:
        data = []
        for index, line in enumerate(file.readlines()):
            # Skip first line
            if index == 0:
                continue
            line_list = line.split('\t')

            # Collect data fields
            hla = line_list[3]
            mhc_norm = float(line_list[5])
            peptide = line_list[6]
            data.append([hla, peptide, mhc_norm])

        return data


def write_to_file(filename, data):
    """
    Write (patient)data to @p filename
    """
    with open(filename, 'w') as f:
        for el in data:
            f.write(str(el[0]) + ',' + str(el[1]) + ',' + str(el[2]) + '\n')


if __name__ == '__main__':
    data = extract()
    write_to_file('data.txt', data)
