import pytest
from laplace.laplace import lc, dnlc

# Test fixtures
@pytest.fixture
def alpha():
    """Alpha value near the 3:1 MMR"""
    return 0.480597


@pytest.fixture
def tolerance():
    """Tolerance for floating point comparisons"""
    return 1e-5


def test_coefficient_A0(alpha, tolerance):
    """Test coefficient A0 calculation"""
    A0 = 0.5 * lc(alpha, 0.5, 0)
    assert abs(A0 - 1.06671) < tolerance


def test_coefficient_A1(alpha, tolerance):
    """Test coefficient A1 calculation"""
    A1 = 0.125 * (
        2 * alpha * dnlc(alpha, 0.5, 0, 1) + alpha * alpha * dnlc(alpha, 0.5, 0, 2)
    )
    assert abs(A1 - 0.142097) < tolerance


def test_coefficient_A2(alpha, tolerance):
    """Test coefficient A2 calculation"""
    A2 = -0.5 * alpha * lc(alpha, 1.5, 1)
    assert abs(A2 - (-0.568387)) < tolerance


def test_coefficient_A3(alpha, tolerance):
    """Test coefficient A3 calculation"""
    A3 = 0.25 * (
        2 * lc(alpha, 0.5, 1)
        - 2 * alpha * dnlc(alpha, 0.5, 1, 1)
        - alpha * alpha * dnlc(alpha, 0.5, 1, 2)
    )
    assert abs(A3 - (-0.165406)) < tolerance


def test_coefficient_A4(alpha, tolerance):
    """Test coefficient A4 calculation"""
    A4 = alpha * lc(alpha, 1.5, 1)
    assert abs(A4 - 1.13677) < tolerance


def test_coefficient_A5(alpha, tolerance):
    """Test coefficient A5 calculation"""
    A5 = 0.125 * (
        21 * lc(alpha, 0.5, 3)
        + 10 * alpha * dnlc(alpha, 0.5, 3, 1)
        + alpha * alpha * dnlc(alpha, 0.5, 3, 2)
    )
    assert abs(A5 - 0.598100) < tolerance


def test_invalid_method():
    """Test that invalid method raises ValueError"""
    with pytest.raises(ValueError, match="method must be 'hyper' or 'brute'"):
        lc(0.5, 0.5, 0, method="invalid")


def test_invalid_alpha():
    """Test that alpha > 1 raises ValueError"""
    with pytest.raises(
        ValueError, match="semi-major axis ratio a must be between 0 and 1"
    ):
        lc(1.5, 0.5, 0)


def test_invalid_s():
    """Test that negative s raises ValueError"""
    with pytest.raises(ValueError, match="s must be a positive half-integer"):
        lc(0.5, -0.5, 0)
