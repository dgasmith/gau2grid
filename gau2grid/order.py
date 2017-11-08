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


def cartesian_order_factory(L, order="row"):
    if order == "row":
        return row_cartesian_order(L)
    else:
        raise KeyError("Cartesian order '%s' not understood" % order)
