import numpy as np
import scipy.integrate as integrate
from scipy.special import hyp2f1


def lc(a, s, j, method="hyper"):
    """
    Laplace coefficent b_s^j(a).

    Args:
        a = alpha, semi-major axis ratio, between 0 and 1.
        s = half-integer power to which the denominator is raised.
        j = appears in argument \cos(j\theta).
        method: - "hyper" uses the hyperbolique function 2F1, see exercice 6.2 of Murray & Dermott.
                - "brute" uses a slower brute-force integration, see eq. (6.67) of Murray & Dermott.

    Returns:
        The Laplace coefficent b_s^j(a).

    Note:
        The hyper method is more accurate and faster than the brute method

    Reference:
        Murray & Dermott, Solar System Dynamics, 1999.
    """

    if method not in ["hyper", "brute"]:
        raise ValueError("method must be 'hyper' or 'brute'.")
    if a > 1:
        raise ValueError("semi-major axis ratio a must be between 0 and 1.")
    if s < 0:
        raise ValueError("s must be a positive half-integer.")

    if method == "hyper":
        j = np.abs(j)
        p = np.prod(np.arange(1, j + 1))  # Factorial j
        q = np.prod(np.arange(s, s + j))
        return 2 * (q / p) * np.power(a, j) * hyp2f1(s, s + j, j + 1, a * a)

    if method == "brute":
        db = lambda t: np.cos(j * t) * np.power(1 - 2 * a * np.cos(t) + a * a, -s)
        l, err = integrate.quad(db, 0.0, 2.0 * np.pi)
        return l / np.pi


def dnlc(a, s, j, n=0, method="hyper"):
    """
    N-th derivative of Laplace coefficent D^n b_s^j(a).

    Args:
        a = alpha, semi-major axis ratio, between 0 and 1.
        s = half-integer power to which the denominator is raised.
        j = appears in argument \cos(j\theta).
        n = integer, order of the derivative. n=0 returns the Laplace coeffiicient. n<0 returns 0.
        method: - "hyper" uses the hyperbolique function 2F1, see exercice 6.2 of Murray & Dermott.
                - "brute" uses a slower brute-force integration, see eq. (6.67) of Murray & Dermott.

    Returns:
        The n-th derivative of the Laplace coefficent b_s^j(a).
    """
    if n <= 0:
        return max(0, n + 1) * lc(a, s, j, method)
    else:
        return s * (
            dnlc(a, s + 1, j - 1, n - 1, method)
            - 2 * a * dnlc(a, s + 1, j, n - 1, method)
            + dnlc(a, s + 1, j + 1, n - 1, method)
            - 2 * (n - 1) * dnlc(a, s + 1, j, n - 2, method)
        )
