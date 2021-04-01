import logging
import numpy as np
from scipy import optimize
from enum import Enum

__author__ = 'Rusty Gentile'

logger = logging.getLogger(__name__)


class SymmetryOptions(Enum):
    LEFT = 1
    RIGHT = 2
    AVERAGE = 3
    NONE = 4


class NodeFinder:
    """
    A system for finding optimal interpolation nodes.
    """
    
    def __init__(self, n, symmetry=1):
         
        self.n = n
        self.symmetry = SymmetryOptions(symmetry)

        self.nodes = np.arange(-1., 1., n)
        self.best_nodes = np.arange(-1., 1., n)
        self.n_guess = 0
        self.best_dl = np.inf

        # Useful for evaluating the Lebesgue function
        self.iden = np.identity(n)
        self.inv_eye = np.ones([n, n]) - self.iden
        self.ones = np.ones([n, n])
        self.lagrange_denominators = np.ones(n)

        # Local maxima of the Lebesgue function
        self.lx1 = np.zeros(n - 1)
        self.lx2 = np.zeros(n - 1)
        self.ly1 = np.zeros(n - 1)
        self.ly2 = np.zeros(n - 1)

        self.dx1 = np.zeros(n - 1)
        self.dx2 = np.zeros(n - 1)
        self.dy1 = np.zeros(n - 1)
        self.dy2 = np.zeros(n - 1)

        # Assume that dx_min will not be less than dx_min for the extended Chebyshev nodes
        cheb = np.polynomial.chebyshev.chebpts1(n)
        cheb = cheb / cheb[-1]
        self.min_dx = np.ones(n - 1) * (cheb[1] - cheb[0]) / 2.

    def force_symmetry(self, arr):
        """
        Prameter
        --------
        arr : arraylike
            A 1D array

        Returns
        -------
        Symmetric 1D array
        """
        
        midpoint = self.n // 2

        if self.symmetry == SymmetryOptions.LEFT:
            arr[-midpoint:] = np.flip(arr[:midpoint])

        elif self.symmetry == SymmetryOptions.RIGHT:
            arr[:midpoint] = np.flip(arr[-midpoint:])

        elif self.symmetry == SymmetryOptions.AVERAGE:
            arr = (arr + np.flip(arr)) / 2

        elif self.symmetry == SymmetryOptions.NONE:
            pass

        else:
            raise NotImplementedError(f'Unsupported symmetry option: '
                                      '{self.symmetry}')

        return arr
    
    def update_lagrange_denominators(self):
        """
        The Lagrange basis polynomials are given by:
        l_i(x) = Pi_{j}^N (x - x_j) / (x_i - x_j), i =/= j

        Here, we calculate the denominators once for a given set of nodes:
        Pi_{j}^N 1 / (x_i - x_j)
        """
        d_m = self.nodes * self.inv_eye
        self.lagrange_denominators = 1. / np.prod(
            ((self.nodes * self.ones)
             .transpose() - d_m) * self.inv_eye
            + self.iden, axis=1)

    def update_nodes(self, nodes):
        """
        Updates the current guess for optimal interpolation

        Parameters
        ----------
        nodes : arraylike
            New nodes.

        Returns
        -------
        None
        """
        self.n_guess += 1
        self.nodes = nodes
        self.update_lagrange_denominators()
        self.lx2 = np.copy(self.lx1)
        self.dx2 = np.copy(self.dx1)
        self.ly2 = np.copy(self.ly1)
        self.dy2 = np.copy(self.dy1)
        self.local_maxima()

    def lebesgue_function(self, x):
        """
        Evaluates the Lebesgue function for the current guess without 
        looping.
        """
        diffs = self.ones * x - self.nodes
        l_bases = self.inv_eye * diffs + np.diag(self.lagrange_denominators)
        return np.sum(np.abs(np.prod(l_bases, axis=1)))

    def local_maxima(self, xtol=1e-12):
        """
        Updates the local maxima of the Lebesgue function for the current guess.
        """
        for i in range(self.n - 1):
            lb = self.nodes[i]
            ub = self.nodes[i + 1]
            res = optimize.fminbound(lambda x: -self.lebesgue_function(x),
                                     lb, ub, xtol=xtol)
            self.lx1[i] = res
            self.ly1[i] = self.lebesgue_function(res)
            self.dx1[i] = ub - lb

        self.dy1 = self.ly1 - np.mean(self.ly1)

    def next_nodes(self):
        """
        Calculates and updates the next guess for optimal interpolation nodes.

        Returns
        -------
        dl_inf - the maximum difference between local optima of the Lebesgue function
        """

        if self.n_guess < 2:
            raise RuntimeError('At least two initial guesses required to start'
                               ' searching for optimal interpolation nodes.')

        # Secant method
        dy = self.dy1 - self.dy2
        dx = self.dx1 - self.dx2
        dxdy = np.divide(dx, dy, out=np.zeros_like(dx), where=abs(dy) > 1e-12)
        roots = self.dx1 - self.dy1 * dxdy

        # Calculate nodes from dx's
        new_dx = self.force_symmetry(roots)
        new_nodes = np.zeros(self.n)
        new_nodes[1:] = np.tril(self.ones[1:, 1:])@new_dx
        new_nodes = 2 * new_nodes / new_nodes[-1] - 1
        self.update_nodes(new_nodes)

        # Check if this is the best guess so far
        dl_inf = np.max(self.ly1) - np.min(self.ly1)
        if dl_inf < self.best_dl:
            self.best_nodes = np.copy(self.nodes)
            self.best_dl = dl_inf

        return dl_inf
