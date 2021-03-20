import logging
import numpy as np
from scipy import optimize
from enum import Enum

logger = logging.getLogger(__name__)


class SymmetryOptions(Enum):
    LEFT = 1
    RIGHT = 2
    AVERAGE = 3
    NONE = 4


class NodeFinder:
    """
    A system for finding optimal interpolation nodes.
    
    Some key fields:
    ----------------
    nodes - The current guess for optimal interpolation.
    best_nodes - The best guess so far for optimal interpolation.
    N - Number of nodes.
    dxs - The distances between each node for each guess.
    dLs - The differences between the Lebesgue function's local and global 
          maxima at each guess.

    Some key methods:
    -----------------
    next_nodes() - Computes the next guess for optimal interpolation.
    ...

    """
    
    def __init__(self, N, symmetry=1):
         
        self.N = N
        self.dxs = []
        self.dLs = []
        self.symmetry = SymmetryOptions(symmetry)

        # Useful for evaluating the Legesgue function
        self.iden = np.identity(N)
        self.inv_eye = np.ones([N, N]) - self.iden
        self.ones = np.ones([N, N])

    def force_symmetry(self, arr):
        """Ensures an array is a palindrome"""
        
        if self.symmetry == SymmetryOptions.NONE:
            return arr
        
        N = self.N - 1
        midpoint = (N - 1.) / 2
        res = np.copy(arr)
        idx = range(np.int(np.ceil(midpoint)))

        if self.symmetry == SymmetryOptions.LEFT:
            for i in idx:
                res[-(i + 1)] = res[i]
        elif self.symmetry == SymmetryOptions.RIGHT:
            for i in idx:
                res[i] = res[-(i + 1)]
        elif self.symmetry == SymmetryOptions.AVERAGE:
            for i in idx:  
                avg = 0.5 * (arr[i] + arr[-(i+1)])
                res[i] = avg
                res[-(i + 1)] = avg

        else:
            raise NotImplementedError(f'Unsupported symmetry opition: '
                                      '{self.symmetry}')

        return res
    
    def update_lagrange_denominators(self):
        """
        The Lagrange basis polynomials are given by:
        l_i(x) = Pi_{j}^N (x - x_j) / (x_i - x_j), i =/= j

        Here, we calculate the denominators once for a given set of nodes:
        Pi_{j}^N 1 / (x_i - x_j)
        """
        d_m = self.nodes * self.inv_eye
        self.lagrange_denominators = 1. / np.prod(((self.nodes * self.ones)
                                                  .transpose() - d_m) * self.inv_eye
                                                  + self.iden, axis=1)

    def update_nodes(self, nodes):
        """Updates the current guess for optimal interpolation"""
        self.nodes = nodes
        self.update_lagrange_denominators()
        _, dxs, _, dLs = self.local_maxima()
        self.dxs.append(dxs)
        self.dLs.append(dLs)

    def lebesgue_function(self, x):
        """
        Evaluates the Lebesgue function for the current guess without 
        looping.
        """
        diffs = xx = self.ones * x - self.nodes
        l = self.inv_eye * diffs + np.diag(self.lagrange_denominators)
        return np.sum(np.abs(np.prod(l, axis=1)))

    def local_maxima(self, xtol=1e-12):
        """
        Returns:
        --------
        xs - x coordinates of the local maxima of the Lebesgue function.
        dxs - Distances bewteen the current set of nodes.
        ls - Local maxima of the Lebesgue function.
        dls - Differences between the local maxima of the Lebesgue function
              and their median.
        """
        xs = []
        dxs = []
        vals = []
        for i in range(1, self.N):
            lb = self.nodes[i - 1]
            ub = self.nodes[i]
            res = optimize.fminbound(lambda x: -self.lebesgue_function(x),
                                     lb, ub, xtol=xtol)
            xs.append(res)
            dxs.append(ub - lb)
            vals.append(self.lebesgue_function(res))
        
        ls = np.array(vals)
        dls = ls - np.median(ls)
        
        return np.array(xs), np.array(dxs), ls, dls

    def next_nodes(self):
        """
        Calculates and updates the next guess for optimal interpolation nodes.
        """
        if len(self.dxs) < 2:
            raise RuntimeError('At least two initial guesses required to start'
                               ' searching for optimal interpolation nodes.')
            
        y1 = self.dLs[-2]
        y2 = self.dLs[-1]
        x1 = self.dxs[-2]
        x2 = self.dxs[-1]
        
        dy = (y2 - y1) / (x2 - x1)
        roots = np.abs(x2 - y2 / dy)
        
        # TODO: get rid of the for loop
        for i in range(self.N - 1):
            if abs(y2[i]) < 1e-3:
                roots[i] = x2[i]
        
        new_dx = self.force_symmetry(roots)
        
        new_nodes = np.zeros(self.N)
        for i in range(1, self.N):
            new_nodes[i] = np.sum(new_dx[:i])
    
        new_nodes = 2 * new_nodes / new_nodes[-1] - 1
        self.update_nodes(new_nodes)
        
        _, _, Ls, _ = self.local_maxima()
        return Ls
    
    def initial_guess(self, N, alpha):
        """A suggestion for initial guesses."""
    
        # The extended Chebyshev nodes
        ig_cheb2 = np.polynomial.chebyshev.chebpts1(N)
        ig_cheb2 = ig_cheb2 / ig_cheb2[-1]
    
        # Perturb the nodes using a small exponent. This should be a slightly 
        # worse guess for optimal interpolation 
        ig_cheb1 = np.abs(ig_cheb2) ** alpha * np.sign(ig_cheb2)
        return ig_cheb1, ig_cheb2
