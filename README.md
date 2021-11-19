# Laplace coefficients
Python script which computes the Laplace coefficients of celestial mechanics and their n-th order derivatives.

Two methods are provided:
- Either a slow brute-force computation of the integral form of the Laplace coefficients (see equation 6.67 of Murray & Dermott 1999).
- Or a faster computation based on the hypergeometric function (see exercice 6.2 of Murray & Dermott 1999, or books by Brouwer & Clemence 1961 or Hagihara 1972).

Laplace coefficients are computed via the `lc` function. Its arguments are:
- a (float): semi-major axis ratio, smaller than unity.
- s (float): half-integer power to which the denominator is raised.
- j (float): Appears in argument \cos(j\theta)
- method (string): 'hyper' uses the hypergeometric function to compute the Laplace coefficients, while 'brute' does a (slower) brute-force integration. If no method is specified, the default method uses the hypergeometric function. 

N-th derivatives are computed via the `dnlc` function. Its arguments are the same as those of the `lc` function, with the addition of:
- n (integer, n>=0): Indicates that the n-th derivative of the Laplace coefficient should be computed. If n=0, it simply returns the Laplace coefficient.
