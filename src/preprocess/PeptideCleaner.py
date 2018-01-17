import argparse

KEYS = ['Description', 'Binds', 'Allele Name']
logfile = "PeptideCleaner.log"


def useful_indices():
    """
    Get list of useful indices from csv file
    """
    return [
        11,     # Allele description
        82,     # Qualitative Measure
        93,     # Allele Name
        96      # MHC allele class
    ]


def binding2binary(data):
    """
    Convert all forms of 'positive' qualitative measure into True/False binding.
    """
    for entry in data:
        entry['Binds'] = entry.pop('Qualitative Measure')

        # if binding is some sort of positive, but not positive-low, set
        # binding to True
        entry['Binds'] = 'positive' in entry['Binds'].lower()

    return data


def remove_non_stars(data):
    """
    Remove all entries without a '*' in the Allele Name.
    """
    return [entry for entry in data if '*' in entry['Allele Name']]


def p_split(text):
    """
    Split line
    """
    return [el.replace('"', '') for el in text.split('","')]


def get_header(filename, header_index=0):
    """
    Get list of columns from @p filename.
    """
    with open(filename) as f:
        for i, line in enumerate(f):
            if i == header_index:
                return line[:-1].split(',')


def p_read_raw_file(filename):
    """
    Read and clean file
    """
    header_index = 1
    columns = get_header(filename, header_index)

    data = []

    with open(filename) as f:
        for i, line in enumerate(f):
            # Skip header row(s)
            if i <= header_index:
                continue

            entry = p_split(line)

            # Skip inconsistent entries
            if len(entry) != 98:
                continue
            data_entry = {}
            try:
                for index in useful_indices():
                    data_entry[columns[index]] = entry[index]
                    # Handle empty value for current field
                    if data_entry[columns[index]] == '':
                        p_error(columns[index], entry)
            except IndexError:
                # Skip erroneous entries
                continue
            data.append(data_entry)

    return data


def p_error(field, entry):
    """
    Log and throw error.
    """
    with open(logfile, 'a+') as log:
        log.write("Failed to add " + str(field) + " to " + str(entry))
        log.close()
    raise IndexError


def p_filter_data(data, field, value):
    """
    Only keep entries with field == value from data
    """
    result = []
    for e in data:
        if e[field] == value:
            result.append(e)
    return result


def p_drop_field(data, field):

    """
    Drop the field @p field from each entry in @p data.

    :param data Dictionary containing a column with key @p field.
    """
    for entry in data:
        del entry[field]
    return data


def p_write_file(filename, data):
    """
    Write formatted @p data to @p filename
    """
    with open(filename, "w") as f:
        line = ""
        for key in KEYS:
            line += key + ","
        f.write(line[:-1] + "\n")

        for entry in data:
            entry = fix_sequence(entry)
            if not entry:
                continue
            line = ""
            # for key in entry:
            for key in KEYS:
                line += str(entry[key]) + ","
            f.write(line[:-1] + "\n")


def fix_sequence(entry):
    """
    Check validity of @p entry. Remove the bad entry.
    """
    # Remove all mutant alleles.
    if 'mutant' in entry['Allele Name']:
        # print("Found mutant")
        # print(entry)
        return None
    # Remove all extra data, store only first part
    if ' ' in entry['Description']:
        # print("Found white space")
        # print(entry)
        return None
    # Remove all entries without a sequence
    if not entry['Description'].isupper():
        # print("Found non-sequence")
        # print(entry)
        return None

    return entry


def main(file):
    """
    Main function
    """
    # Clear logfile
    with open(logfile, 'w') as f:
        f.close()

    # Read all data
    data = p_read_raw_file(file)

    # Remove all non-class I entries
    data = p_filter_data(data, 'MHC allele class', 'I')

    # Remove all entries without a *
    data = remove_non_stars(data)

    # Change qualitative measure into yes/no binding data
    data = binding2binary(data)

    # Remove useless fields
    data = p_drop_field(data, 'MHC allele class')
    # for el in data:
    #     print(el['Qualitative Measure'])

    # Generate new data file
    p_write_file("AlleleData.csv", data)

    print("Generated " + str(len(data)) + " peptide entries.")
