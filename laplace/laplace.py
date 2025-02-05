import numpy as np
import scipy.integrate as integrate
from scipy.special import hyp2f1

def lc(a, s, j, method='hyper'):
    """
        Laplace coefficent b_s^j(a).

        Arguments:
            a = alpha, semi-major axis ratio, between 0 and 1.
            s = half-integer power to which the denominator is raised.
            j = appears in argument \cos(j\theta).
            method: - "hyper" uses the hyperbolique function 2F1, see exercice 6.2 of Murray & Dermott.
                    - "brute" uses a slower brute-force integration, see eq. (6.67) of Murray & Dermott.
    """
    if method == 'hyper':
        j = np.abs(j)
        p = np.prod(np.arange(1,j+1)) #--Factorial j
        q = np.prod(np.arange(s,s+j))
        return 2 * (q/p) * np.power(a,j) * hyp2f1(s, s+j, j+1, a*a)

    if method == 'brute':
        db = lambda t: np.cos(j*t)*np.power(1 - 2*a*np.cos(t) + a*a, -s)
        l,err = integrate.quad(db, 0., 2.*np.pi)
        return l/np.pi

def dnlc(a, s, j, n=0, method='hyper'):
    """
        N-th derivative of Laplace coefficent D^n b_s^j(a).

        Arguments:
            a = alpha, semi-major axis ratio, between 0 and 1.
            s = half-integer power to which the denominator is raised.
            j = appears in argument \cos(j\theta).
            n = integer, order of the derivative. n=0 returns the Laplace coeffiicient. n<0 returns 0.
            method: - "hyper" uses the hyperbolique function 2F1, see exercice 6.2 of Murray & Dermott.
                    - "brute" uses a slower brute-force integration, see eq. (6.67) of Murray & Dermott.
    """
    if n<=0:
        return max(0, n+1) * lc(a, s, j, method)
    else:
        return s*( dnlc(a, s+1 ,j-1, n-1, method) - 2*a*dnlc(a, s+1, j, n-1, method) + dnlc(a, s+1, j+1, n-1, method) - 2*(n-1)*dnlc(a, s+1 ,j, n-2, method))


if __name__ == '__main__':
    import prettytable as pt
    # Use Table 6.1 (page 263) of Murray & Dermott to check coefficents:
    a = 0.480597 # alpha near the 3:1 MMR
    A0 = 0.5 * lc(a,0.5,0)
    A1 = 0.125 * ( 2*a*dnlc(a,0.5,0,1) +  a*a*dnlc(a,0.5,0,2) )
    A2 = -0.5 * a * lc(a,1.5,1)
    A3 = 0.25 * ( 2*lc(a,0.5,1) - 2*a*dnlc(a,0.5,1,1) - a*a*dnlc(a,0.5,1,2) )
    A4 = a * lc(a,1.5,1)
    A5 = 0.125 * ( 21*lc(a,0.5,3) + 10*a*dnlc(a,0.5,3,1) + a*a*dnlc(a,0.5,3,2) )
    print()
    print("Code verification using Table 6.1 of Murray & Dermott")
    t = pt.PrettyTable(['i', 'A_i (code)',  'A_i (M&D)'], padding_width=3, junction_char='-')
    t.float_format = f'.6'
    t.add_row(['0', A0, 1.06671] )
    t.add_row(['1', A1, 0.142097] )
    t.add_row(['2', A2, -0.568387] )
    t.add_row(['3', A3, -0.165406] )
    t.add_row(['4', A4, 1.13677] )
    t.add_row(['5', A5, 0.598100] )
    print(t)
    print()
