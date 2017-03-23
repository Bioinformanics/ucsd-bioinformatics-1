
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


def ConvertTextToMatrix(text, num_columns):
    text_array = text.split(" ")
    float_array = map(lambda t: float(t), text_array)
    return [float_array[i:i+num_columns] for i in range(0, len(float_array), num_columns)]

def reduce(function, iterable, initializer=None):
    it = iter(iterable)
    if initializer is None:
        value = next(it)
    else:
        value = initializer
    for element in it:
        value = function(value, element)
    return value