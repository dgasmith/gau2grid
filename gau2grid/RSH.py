"""
Cartesian to regular solid harmonics conversion code.
"""

import numpy as np

# Arbitrary precision math with 100 decimal places
import mpmath
mpmath.mp.dps = 100

from . import order


class Memoize(object):
    """
    Simple memoize class for RSH_coefs which is quite expensive
    """

    def __init__(self, func):
        self.func = func
        self.mem = {}

    def __call__(self, *args):
        if args not in self.mem:
            self.mem[args] = self.func(*args)
        return self.mem[args]


def quanta_to_string(lx, ly, lz):
    """Pretty print monomials with quanta lx, ly, lz."""
    string = ""
    string += "X" * lx
    string += "Y" * ly
    string += "Z" * lz
    # if lx:
    #     string += 'x'
    # if lx > 1:
    #     string += '^{}'.format(lx)
    # if ly:
    #     string += 'y'
    # if ly > 1:
    #     string += '^{}'.format(ly)
    # if lz:
    #     string += 'z'
    # if lz > 1:
    #     string += '^{}'.format(lz)
    return string


@Memoize
def cart_to_RSH_coeffs(l):
    """
    Generates a coefficients [ coef, x power, y power, z power ] for each component of
    a regular solid harmonic (in terms of raw Cartesians) with angular momentum l.

    See eq. 23 of ACS, F. C. Pickard, H. F. Schaefer and B. R. Brooks, JCP, 140, 184101 (2014)

    Returns coeffs with order 0, +1, -1, +2, -2, ...
    """

    terms = []
    for m in range(l + 1):
        thisterm = {}
        # p1 = mpmath.sqrt(mpmath.fac(l - m / mpmath.fac(l + m))) * (mpmath.fac(m) / (2**l))
        p1 = mpmath.sqrt((mpmath.fac(l - m)) / (mpmath.fac(l + m))) * ((mpmath.fac(m)) / (2**l))
        if m:
            p1 *= mpmath.sqrt(2.0)
        # Loop over cartesian components
        for lz in range(l + 1):
            for ly in range(l - lz + 1):
                lx = l - ly - lz
                xyz = lx, ly, lz
                j = int((lx + ly - m) / 2)
                if (lx + ly - m) % 2 == 1 or j < 0:
                    continue
                p2 = mpmath.mpf(0)
                for i in range(int((l - m) / 2) + 1):
                    if i >= j:
                        p2 += (-1)**i * mpmath.fac(2 * l - 2 * i) / (
                            mpmath.fac(l - i) * mpmath.fac(i - j) * mpmath.fac(l - m - 2 * i))
                p3 = mpmath.mpf(0)
                for k in range(j + 1):
                    if (j >= k) and (lx >= 2 * k) and (m + 2 * k >= lx):
                        p3 += (-1)**k / (
                            mpmath.fac(j - k) * mpmath.fac(k) * mpmath.fac(lx - 2 * k) * mpmath.fac(m - lx + 2 * k))
                p = p1 * p2 * p3
                # print(p)
                if xyz not in thisterm:
                    thisterm[xyz] = [mpmath.mpf(0.0), mpmath.mpf(0.0)]

                if (m - lx) % 2:
                    # imaginary
                    sign = mpmath.mpf(-1.0)**mpmath.mpf((m - lx - 1) / 2.0)
                    thisterm[xyz][1] += sign * p
                else:
                    # real
                    sign = mpmath.mpf(-1.0)**mpmath.mpf((m - lx) / 2.0)
                    thisterm[xyz][0] += sign * p

        tmp_R = []
        tmp_I = []
        for k, v in thisterm.items():
            if abs(v[0]) > 0:
                tmp_R.append((k, v[0]))
            if abs(v[1]) > 0:
                tmp_I.append((k, v[1]))

        if m == 0:
            # name_R = "R_%d%d" % (l, m)
            terms.append(tmp_R)
        else:
            # name_R = "R_%d%dc" % (l, m)
            # name_I = "R_%d%ds" % (l, m)
            terms.append(tmp_R)
            terms.append(tmp_I)
            # terms[name_R] = tmp_R
            # terms[name_I] = tmp_I

        # for k, v in terms.items():
        #     print(k, v)

    return terms


def cart_to_spherical_transform(data, L, cart_order):
    """
    Transforms a cartesian x points matrix into a spherical x points matrix.
    """

    cart_order = {x[1:]: x[0] for x in order.cartesian_order_factory(L, cart_order)}
    RSH_coefs = cart_to_RSH_coeffs(L)

    nspherical = len(RSH_coefs)
    ret = np.zeros((nspherical, data.shape[1]))

    idx = 0
    for spherical in RSH_coefs:
        for cart_index, scale in spherical:
            ret[idx] += float(scale) * data[cart_order[cart_index]]
        idx += 1

    return ret


def transformation_np_generator(cg, L, cart_order, function_name="generate_transformer"):
    """
    Builds a conversion from cartesian to spherical coordinates
    """

    cart_order = {x[1:]: x[0] for x in order.cartesian_order_factory(L, cart_order)}
    RSH_coefs = cart_to_RSH_coeffs(L)

    nspherical = len(RSH_coefs)

    cg.write("def " + function_name + "_%d(data, out=None):" % L)
    cg.indent()
    cg.write("if out is None:")
    cg.write("    out = np.zeros((%d, data.shape[1]))" % nspherical)


    cg.blankline()
    cg.write("# Contraction loops")

    idx = 0
    for spherical in RSH_coefs:
        op = " ="
        for cart_index, scale in spherical:
            if scale != 1.0:
                cg.write("out[%d][:] %s % .16f * data[%d]" % (idx, op, scale, cart_order[cart_index]))
            else:
                cg.write("out[%d][:] %s data[%d]" % (idx, op, cart_order[cart_index]))
            # cg.write("print(np.linalg.norm(out[%d]))" % idx)
            op = "+="
        cg.write("")
        idx += 1

    cg.write("return out")
    cg.dedent()


def transformation_c_generator(cg, L, cart_order, function_name=""):
    """
    Builds a conversion from cartesian to spherical coordinates in C
    """

    if function_name == "":
        function_name = "cart_to_spherical_L%d" % L

    cart_order = {x[1:]: x[0] for x in order.cartesian_order_factory(L, cart_order)}
    RSH_coefs = cart_to_RSH_coeffs(L)

    nspherical = len(RSH_coefs)

    signature = "void %s(const size_t size, const double* __restrict__ cart, const size_t ncart, double* __restrict__ spherical, const size_t nspherical)" % function_name

    # Start function
    cg.start_c_block(signature)

    cg.write("// R_%d0 Transform" % L)
    _c_spherical_trans(cg, 0, RSH_coefs, cart_order)
    cg.blankline()

    for l in range(L):
        cg.write("// R_%d%dc Transform" % (L, l + 1))
        sidx = 2 * l + 1
        _c_spherical_trans(cg, sidx, RSH_coefs, cart_order)

        sidx = 2 * l + 2
        cg.write("// R_%d%ds Transform" % (L, l + 1))
        _c_spherical_trans(cg, sidx, RSH_coefs, cart_order)
        cg.blankline()

    # End function
    cg.close_c_block()
    return signature


def _c_spherical_trans(cg, sidx, RSH_coefs, cart_order):
    # cg.write("#pragma clang loop vectorize(assume_safety)")
    cg.start_c_block("for (size_t i = 0; i < size; i++)")

    # Figure out where we are summing to
    if sidx == 0:
        lhs = "spherical[i]"
    elif sidx == 1:
        lhs = "spherical[nspherical + i]"
    else:
        lhs = "spherical[%d * nspherical + i]" % sidx

    op = " ="
    for cart_index, scale in RSH_coefs[sidx]:

        # Figure out car idx
        idx = cart_order[cart_index]
        if idx == 0:
            rhs = "cart[i]"
        elif idx == 1:
            rhs = "cart[ncart + i]"
        else:
            rhs = "cart[%d * ncart + i]" % idx

        # Scales
        if scale != 1.0:
            cg.write("%s %s % .16f * %s" % (lhs, op, scale, rhs))
        else:
            cg.write("%s %s %s" % (lhs, op, rhs))
        op = "+="
    cg.blankline()

    cg.close_c_block()
