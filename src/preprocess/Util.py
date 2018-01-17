import re


def reformat(string):
    """
    Turn MHC name into structure HLA-A*01:01
    """
    # Split all components
    components = re.findall(r"[\w']+", string)

    # Create string without special symbols
    components = ''.join(components)

    # Insert symbols where needed
    result = components[:3] + '-' + components[3] + '*' + components[4:6] + ':' + components[6:]

    assert (len(result) == 11)

    return result


def levenshtein(s1, s2):
    """
    Copied from https://en.wikibooks.org/wiki/Algorithm_Implementation/Strings/Levenshtein_distance#Python
    """
    if len(s1) < len(s2):
        return levenshtein(s2, s1)

    # len(s1) >= len(s2)
    if len(s2) == 0:
        return len(s1)

    previous_row = range(len(s2) + 1)
    for i, c1 in enumerate(s1):
        current_row = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = previous_row[
                             j + 1] + 1  # j+1 instead of j since previous_row and current_row are one character longer
            deletions = current_row[j] + 1  # than s2
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row

    return previous_row[-1]


def find_closest_match(string, target):
    """
    Find closest match to @p target in @p string
    :return: (word, index)

    with word: matched substring of @p string
    with index: starting index of word in @p string
    """
    if len(string) <= len(target):
        return string

    min_dist = 1000
    (best_match, index) = ('', 0)
    for i in range(len(string) - len(target)):
        word = string[i:i+len(target)]
        dist = levenshtein(word, target)
        if dist < min_dist:
            min_dist = dist
            (best_match, index) = (word, i)

    return (best_match, index)
