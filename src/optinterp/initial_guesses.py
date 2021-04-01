import numpy as np

__author__ = 'Rusty Gentile'


def chebyshev_exponent(n, alpha):
    """A suggestion for initial guesses."""

    # The extended Chebyshev nodes
    ig_cheb2 = np.polynomial.chebyshev.chebpts1(n)
    ig_cheb2 = ig_cheb2 / ig_cheb2[-1]

    # This should be a slightly worse guess for optimal interpolation
    ig_cheb1 = np.abs(ig_cheb2) ** alpha * np.sign(ig_cheb2)
    return ig_cheb1, ig_cheb2
