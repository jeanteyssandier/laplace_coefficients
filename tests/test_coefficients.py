import pytest
from laplace.coefficients import lc, dnlc

# Use Table 6.1 (page 263) of Murray & Dermott to check coefficents:

# Test fixtures
@pytest.fixture
def alpha():
    """Alpha value near the 3:1 MMR"""
    return 0.480597


@pytest.fixture
def tolerance():
    """Tolerance for floating point comparisons"""
    return 1e-5


@pytest.mark.parametrize("method", ["hyper", "brute"])
def test_coefficient_A0(alpha, tolerance, method):
    """Test coefficient A0 calculation"""
    A0 = 0.5 * lc(alpha, 0.5, 0, method=method)
    assert abs(A0 - 1.06671) < tolerance


@pytest.mark.parametrize("method", ["hyper", "brute"])
def test_coefficient_A1(alpha, tolerance, method):
    """Test coefficient A1 calculation"""
    A1 = 0.125 * (
        2 * alpha * dnlc(alpha, 0.5, 0, 1, method=method)
        + alpha * alpha * dnlc(alpha, 0.5, 0, 2, method=method)
    )
    assert abs(A1 - 0.142097) < tolerance


@pytest.mark.parametrize("method", ["hyper", "brute"])
def test_coefficient_A2(alpha, tolerance, method):
    """Test coefficient A2 calculation"""
    A2 = -0.5 * alpha * lc(alpha, 1.5, 1, method=method)
    assert abs(A2 - (-0.568387)) < tolerance


@pytest.mark.parametrize("method", ["hyper", "brute"])
def test_coefficient_A3(alpha, tolerance, method):
    """Test coefficient A3 calculation"""
    A3 = 0.25 * (
        2 * lc(alpha, 0.5, 1, method=method)
        - 2 * alpha * dnlc(alpha, 0.5, 1, 1, method=method)
        - alpha * alpha * dnlc(alpha, 0.5, 1, 2, method=method)
    )
    assert abs(A3 - (-0.165406)) < tolerance


@pytest.mark.parametrize("method", ["hyper", "brute"])
def test_coefficient_A4(alpha, tolerance, method):
    """Test coefficient A4 calculation"""
    A4 = alpha * lc(alpha, 1.5, 1, method=method)
    assert abs(A4 - 1.13677) < tolerance


@pytest.mark.parametrize("method", ["hyper", "brute"])
def test_coefficient_A5(alpha, tolerance, method):
    """Test coefficient A5 calculation"""
    A5 = 0.125 * (
        21 * lc(alpha, 0.5, 3, method=method)
        + 10 * alpha * dnlc(alpha, 0.5, 3, 1, method=method)
        + alpha * alpha * dnlc(alpha, 0.5, 3, 2, method=method)
    )
    assert abs(A5 - 0.598100) < tolerance


def test_invalid_method():
    """Test that invalid method raises ValueError"""
    with pytest.raises(ValueError, match="method must be 'hyper' or 'brute'"):
        lc(0.5, 0.5, 0, method="invalid")


def test_invalid_alpha():
    """Test that alpha > 1 raises ValueError"""
    with pytest.raises(
        ValueError, match="semi-major axis ratio alpha must be between 0 and 1"
    ):
        lc(1.5, 0.5, 0)


def test_invalid_s():
    """Test that negative s raises ValueError"""
    with pytest.raises(ValueError, match="s must be a positive half-integer"):
        lc(0.5, -0.5, 0)
