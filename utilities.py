
def ToSingleLineOfString(array):
    return ' '.join([str(value) for value in array])


def ConvertBaseToIndex(base):
    if base == 'A':
        index = 0
    elif base == 'C':
        index = 1
    elif base == 'G':
        index  = 2
    elif base == 'T':
        index = 3
    else:
        raise Exception('Invalid base type: %s' % base)
    return index


def ConvertIndexToBase(index):
    if index not in range(4):
        raise Exception("Cannot convert index '%d' to nucleic acid." % index)
    bases = ['A', 'C', 'G', 'T']
    return bases[index]


def ConvertTextToMatrix(rows, converter = None):
    """
    :param rows: list or matrix rows
    :return: a matrix parsed from text (i.e. a 2D array)
    """
    matrix = []
    for row in rows:
        row_with_text_values = row.strip().split(" ")
        if converter:
            matrix.append([converter(str) for str in row_with_text_values])
        else:
            matrix.append(row_with_text_values)
    return matrix


def AreStringListsEqual(list1, list2):
    if len(list1) != len(list2):
        return False
    for val in list1:
        if val not in list2:
            return False
    for val in list2:
        if val not in list1:
            return False
    return True


def reduce(function, iterable, initializer=None):
    it = iter(iterable)
    if initializer is None:
        value = next(it)
    else:
        value = initializer
    for element in it:
        value = function(value, element)
    return value