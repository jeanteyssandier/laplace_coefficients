# Laplace coefficients
Python script which computes the Laplace coefficients of celestial mechanics and their n-th order derivatives.

The Laplace coefficients are defined by the integral:

b_s^j(α) = (2/π) ∫[from 0 to π] cos(jψ) / (1 - 2α cos(ψ) + α²)^s dψ

These are special functions that appear when analyzing gravitational interactions between orbiting bodies, particularly in perturbation theory. They are especially useful when dealing with planetary perturbations and resonances.

## Overview

Two methods are provided:
- Either a slow brute-force computation of the integral form of the Laplace coefficients (see equation 6.67 of Murray & Dermott 1999).
- Or a faster computation based on the hypergeometric function (see exercice 6.2 of Murray & Dermott 1999, or books by Brouwer & Clemence 1961 or Hagihara 1972).

Laplace coefficients are computed via the `lc` function. Its arguments are:
- a (float): semi-major axis ratio, satisfying 0<a<1.
- s (float): half-integer power to which the denominator is raised.
- j (float): Appears in argument \cos(j\theta)
- method (string): 'hyper' uses the hypergeometric function to compute the Laplace coefficients, while 'brute' does a (slower) brute-force integration. If no method is specified, the default method uses the hypergeometric function. 

N-th derivatives are computed via the `dnlc` function. Its arguments are the same as those of the `lc` function, with the addition of:
- n (integer, n>=0): Indicates that the n-th derivative of the Laplace coefficient should be computed. If n=0, it simply returns the Laplace coefficient.

## Usage

```bash
from laplace.coefficients import lc, dnlc

# Calculate basic Laplace coefficient b_{1/2}^1(0.5)
coefficient = lc(a=0.5, s=0.5, j=1)

# Calculate first derivative
derivative = dnlc(a=0.5, s=0.5, j=1, n=1)
```

## Installation

You can use `pyenv` and the `pyproject.toml` file to set up an environment and install dependencies:

```bash
# Install Python 3.8
pyenv install 3.8

# Create a virtual environment
pyenv virtualenv 3.8 laplace-env

# Activate the environment
pyenv activate laplace-env
```
And then install the package:

```bash
pip install -e .
```