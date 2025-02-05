import numpy as np
from scipy.integrate import quad
from scipy.special import hyp2f1


def _validate_inputs(alpha: float, s: float, method: str) -> None:
    """
    Validate input parameters for Laplace coefficient calculations.

    Args:
        alpha: Semi-major axis ratio (must be between 0 and 1)
        s: Power to which the denominator is raised (must be positive)
        method: Calculation method ('hyper' or 'brute')

    Raises:
        ValueError: If any input parameters are invalid
    """
    if method not in ["hyper", "brute"]:
        raise ValueError("method must be 'hyper' or 'brute'")
    if alpha > 1:
        raise ValueError("semi-major axis ratio alpha must be between 0 and 1")
    if s < 0:
        raise ValueError("s must be a positive half-integer")


def lc(alpha: float, s: float, j: int, method: str = "hyper") -> float:
    """
    Calculate the Laplace coefficient b_s^j(alpha).

    Args:
        alpha: Semi-major axis ratio (must be between 0 and 1)
        s: Half-integer power to which the denominator is raised
        j: Integer appearing in the argument cos(j*theta)
        method: Calculation method, either:
            - "hyper": Uses hypergeometric function 2F1 (faster, more accurate)
            - "brute": Uses direct numerical integration

    Returns:
        float: The Laplace coefficient b_s^j(alpha)

    Examples:
        >>> lc(0.5, 0.5, 1)  # Calculate b_{1/2}^1(0.5)
        0.555866197926681
    """
    _validate_inputs(alpha, s, method)

    if method == "hyper":
        j = abs(j)  # Ensure j is non-negative for hypergeometric calculation
        p = np.prod(np.arange(1, j + 1))  # Factorial of j
        q = np.prod(np.arange(s, s + j))
        return 2 * (q / p) * np.power(alpha, j) * hyp2f1(s, s + j, j + 1, alpha * alpha)

    # Brute force integration method
    integrand = lambda t: np.cos(j * t) * np.power(
        1 - 2 * alpha * np.cos(t) + alpha * alpha, -s
    )
    result, _ = quad(integrand, 0.0, 2.0 * np.pi)
    return result / np.pi


def dnlc(alpha: float, s: float, j: int, n: int = 0, method: str = "hyper") -> float:
    """
    Calculate the nth derivative of the Laplace coefficient D^n b_s^j(alpha).

    Args:
        alpha: Semi-major axis ratio (must be between 0 and 1)
        s: Half-integer power to which the denominator is raised
        j: Integer appearing in the argument cos(j*theta)
        n: Order of the derivative (n=0 returns the Laplace coefficient)
        method: Calculation method ('hyper' or 'brute')

    Returns:
        float: The nth derivative of the Laplace coefficient

    Notes:
        - For n <= 0, returns max(0, n+1) * lc(alpha, s, j, method)
        - Uses a recursive formula for n > 0
    """
    _validate_inputs(alpha, s, method)

    if n <= 0:
        return max(0, n + 1) * lc(alpha, s, j, method)

    # Recursive calculation for higher derivatives
    return s * (
        dnlc(alpha, s + 1, j - 1, n - 1, method)
        - 2 * alpha * dnlc(alpha, s + 1, j, n - 1, method)
        + dnlc(alpha, s + 1, j + 1, n - 1, method)
        - 2 * (n - 1) * dnlc(alpha, s + 1, j, n - 2, method)
    )
