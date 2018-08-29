"""
Contains the different possible cartesian and spherical ordering codes.
"""


def row_cartesian_order(L):
    idx = -1
    for i in range(L + 1):
        l = L - i
        for j in range(i + 1):
            m = i - j
            n = j
            idx += 1
            yield (idx, l, m, n)


def molden_cartesian_order(L):
    # http://www.cmbi.ru.nl/molden/molden_format.html
    if L == 0:
        data = [(0, 0, 0, 0)]
    elif L == 1:
        data = [(0, 1, 0, 0), (1, 0, 1, 0), (2, 0, 0, 1)]
    elif L == 2:
        data = [(0, 2, 0, 0), (1, 0, 2, 0), (2, 0, 0, 2), (3, 1, 1, 0), (4, 1, 0, 1), (5, 0, 1, 1)]
    elif L == 3:
        data = [(0, 3, 0, 0), (1, 0, 3, 0), (2, 0, 0, 3), (3, 2, 1, 0), (4, 2, 0, 1), (5, 1, 0, 2), (6, 1, 2, 0),
                (7, 0, 1, 2), (8, 0, 1, 2), (9, 1, 1, 1)]
    elif L == 4:
        data = [(0, 4, 0, 0), (1, 0, 4, 0), (2, 0, 0, 4), (3, 2, 2, 0), (4, 3, 0, 1), (5, 1, 3, 0), (6, 0, 3,
                                                                                                     1), (7, 1, 0, 3),
                (8, 0, 1, 3), (9, 2, 2, 0), (10, 2, 0, 2), (11, 0, 2, 2), (12, 2, 1, 1), (13, 1, 2, 1), (14, 1, 1, 2)]
    else:
        raise KeyError("Molden ordering only goes to 4 (G)")

    for x in data:
        yield x


def cartesian_order_factory(L, order="row"):
    if order.lower() in ["row", "cca"]:
        return row_cartesian_order(L)
    if order.lower() in ["molden", "libint"]:
        return molden_cartesian_order(L)
    else:
        raise KeyError("Cartesian order '%s' not understood" % order)
