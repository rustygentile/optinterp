# Optimal Interpolation Nodes

Computes a set of points with an optimal [Lebesgue constant](https://en.wikipedia.org/wiki/Lebesgue_constant). Having an exact analytical solution for these points is an unresolved problem in mathematics. The code here approximates the values numerically. Numerical approximations have been available for many years, but to my knowledge, there are no other open source Python libraries with this functionality.

## Usage

### Installation

```
pip install optinterp
```

### Example

```python
import optinterp

nds = optinterp.nodes(10)
```
This functions similarly to numpy's `chebpts1` but produces points with a slightly improved Lebesgue constant:
```python
import numpy as np

nds = np.polynomial.chebyshev.chebpts1(10)
nds = nds / nds[-1]
```

## Algorithm Description

This solution expoits the following properties:  
* Optimal interpolation points can take values -1 and 1 for their minimum and maximum.
* To mimimize the global maximum of the Lebesgue function, all local maxima should be equal.
* Moving two adjacent nodes closer together reduces the local maximum of the Lebesgue function at the expense of increasing the other local maxima.

Start with two initial guesses for optimal points. The extended Chebyshev nodes and a set of slightly perturbed Chebyshev nodes. Then for each set of nodes define:

```
dx_i = x_{i+1} - x_i
dL_i = L_i - L_{avg}
```
Where `L_i` is the local maximum of the Lebesgue function between `x_{i+1}` and `x_i`. Now assuming each `dL_i` is a function of `dx_i`, use the Secant method to find roots:
```
dx_{i, n+1} = dx_{i, n} - dL_{i, n} * (dx_{i, n} - dx_{i, n-1}) / (dL_{i, n} - dL_{i, n - 1})
```
For the next iteration, calculate each node `x_i` from these roots of `dx_i` and scale the values to be from -1 to 1.
