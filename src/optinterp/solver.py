import logging
from . import NodeFinder
from . import chebyshev_exponent

__author__ = 'Rusty Gentile'

logger = logging.getLogger(__name__)


def nodes(n, ltol=1e-8, max_iter=100, alpha=1.000001, ig=None, symmetry=1):
    """
    Parameters
    ----------
    n - Number of nodes.
    ltol - Lebesgue function tolerance. All local maxima should be within one another by this value.
    max_iter - Maximum iterations.
    alpha - Exponent for calculating initial guesses.
    ig - Initial guesses.
    symmetry - Force symmetry options: 1 = right, 2 = left, 3 = average, 4 = none

    Returns
    -------
    optimal interpolation nodes based on the infinity norm of the Lebesgue function.
    """
    nf = NodeFinder(n, symmetry=symmetry)
    if ig is None:
        ig1, ig2 = chebyshev_exponent(n, alpha)
    else:
        ig1 = ig[0]
        ig2 = ig[1]

    nf.update_nodes(ig1)
    nf.update_nodes(ig2)

    for i in range(max_iter):
        try:
            dl = nf.next_nodes()
        except ValueError:
            break

        if abs(dl) <= ltol:
            return nf.nodes

    logger.warning(f'Unable to find optimal nodes after {i} iterations. Returning nodes with ltol={nf.best_dl}')
    return nf.best_nodes