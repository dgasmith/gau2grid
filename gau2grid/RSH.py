import numpy as np
from scipy.misc import factorial

def quanta_to_string(lx, ly, lz):
    """Pretty print monomials with quanta lx, ly, lz."""
    string = ""
    if lx:
        string += 'x'
    if lx > 1:
        string += '^{}'.format(lx)
    if ly:
        string += 'y'
    if ly > 1:
        string += '^{}'.format(ly)
    if lz:
        string += 'z'
    if lz > 1:
        string += '^{}'.format(lz)
    return string


def cart_to_RSH(l):
    """Generates a coefficients [ coef, x power, y power, z power ] for each component of
       a regular solid harmonic (in terms of raw Cartesians) with angular momentum l.

       See eq. 23 of ACS, F. C. Pickard, H. F. Schaefer and B. R. Brooks, JCP, 140, 184101 (2014)"""

    terms = []
    for m in range(l+1):
        thisterm = {}
        p1 = np.sqrt(float(factorial(l-m))/float(factorial(l+m))) * (float(factorial(m)) / float(2**l))
        if m:
            p1 *= np.sqrt(2.0)
        # Loop over cartesian components
        for lz in range(l+1):
            for ly in range(l-lz+1):
                lx  = l - ly - lz
                xyz = lx,ly,lz
                j = int((lx + ly - m) / 2)
                if ((lx + ly - m)%2 == 1 or j < 0): continue
                p2 = 0
                for i in range(int((l-m)/2)+1):
                    if i >= j:
                        p2 += (-1)**i * factorial(2*l-2*i) / float(factorial(l-i)*factorial(i-j)*factorial(l-m-2*i))
                p3 = 0
                for k in range(j+1):
                    if j >= k and lx >= 2*k and m+2*k >= lx:
                        p3 += (-1)**k / float(factorial(j-k)*factorial(k)*factorial(lx-2*k)*factorial(m-lx+2*k)) 
                p = p1 * p2 * p3
                if xyz not in thisterm:
                    thisterm[xyz] = [0.0, 0.0]
                if (m-lx)%2:
                    # imaginary
                    sign = (-1.0)**((m-lx-1)/2)
                    thisterm[xyz][1] += sign * p
                else:
                    # real
                    sign = (-1.0)**((m-lx)/2)
                    thisterm[xyz][0] += sign * p
        # Gather terms
        real_terms = []
        imag_terms = []
        for x,y,z in thisterm.keys():
            real = thisterm[x,y,z][0]
            imag = thisterm[x,y,z][1]
            if real != 0.0:
                sign = ' - ' if real<0 else ' + '
                real_terms.extend([sign, str(np.abs(real)) + '*' + quanta_to_string(x, y, z)])
            if imag != 0.0:
                sign = ' - ' if imag<0 else ' + '
                imag_terms.extend([sign, str(np.abs(imag)) + '*' + quanta_to_string(x, y, z)])
        real_terms[0] = '-' if real_terms[0] == ' - ' else ''
        real_terms = "".join(real_terms)
        if len(imag_terms):
            imag_terms[0] = '-' if imag_terms[0] == ' - ' else ''
            imag_terms = "".join(imag_terms)
        if m == 0:
            # Only real terms here
            terms.append("R_{}{}0 = ".format(l,m) + real_terms)
        else:
            terms.append("R_{}{}c = ".format(l,m) + real_terms)
            terms.append("R_{}{}s = ".format(l,m) + imag_terms)
    return "\n".join(terms)

if __name__ == "__main__":
    for l in range(1,6):
        print(cart_to_RSH(l) + "\n\n")
