import numbers

def ifnumeric(a):
    if not isinstance(a, numbers.Real):
        raise ValueError("invalid numeric field: " + str(a))

def ifpointlist(pp):
    if not isinstance(pp, list):
        raise ValueError("invalid point")
    for p in pp:
        if not isinstance(p, list) or len(p) < 2 or\
                not isinstance(p[0], numbers.Real) or\
                not isinstance(p[1], numbers.Real):
            raise ValueError("invalid point")


def ifincreasing(vec):
    ifnumericlist(vec)
    for i in range(len(vec) - 1):
        if (vec[i + 1] <= vec[i]):
            raise ValueError("vector should be increasing")


def ifnumericlist(vec):
    for x in vec:
        if not isinstance(x, numbers.Real):
            raise ValueError("invalid numeric field: " + str(vec))


def ifintlist(vec):
    for x in vec:
        if not isinstance(x, numbers.Integral):
            raise ValueError("invalid integer values: " + str(vec))

def ifstrlist(vec):
    for x in vec:
        if not isinstance(x, str):
            raise ValueError("invalid string values: " + str(vec))
