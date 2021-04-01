import numpy as np

__author__ = 'Rusty Gentile'


def chebyshev_exp(n: int, alpha: float):
    """
    Used for initial guesses.

    Parameters
    ----------
    n : int
        Number of nodes
    alpha : float
        Exponent

    Returns
    -------
    The extended Chebyshev nodes raised to a small exponent.
    """
    cheb_exp = np.polynomial.chebyshev.chebpts1(n)
    cheb_exp = cheb_exp / cheb_exp[-1]
    return np.abs(cheb_exp) ** alpha * np.sign(cheb_exp)
