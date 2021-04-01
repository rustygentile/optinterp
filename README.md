# Optimal Interpolation Nodes

## Background
Say you have a continuous function and you'd like to do a polynomial fit on it. The algorithm for polynomial interpolation is straightforward, but selecting the points to interpolate is not. One natural choice might be evenly spaced nodes. For example, on the interval [a, b]:

<!---TODO: Rewrite this in LaTex if Github ever supports it :(--->
```
x0 = a  
x1 = a + (b - a) / n  
x2 = a + 2 * (b - a) / n  
...
x{n-1} = a + n * (b - a) / n = b  
```

As it turns out, this actually is a poor choice. For some continuous [functions](https://en.wikipedia.org/wiki/Runge%27s_phenomenon), the fit will even get worse as you increase polynomial order and add interpolation points.

One solution to this problem is to use [Chebyshev](https://en.wikipedia.org/wiki/Chebyshev_nodes) nodes, or the so called "extended" Chebyshev nodes which are scaled to be from -1 to 1. The extended Chebyshev nodes are close to optimal, but having an exact solution for optimal interpolation points remains an unresolved problem in mathematics.

This code approximates the optimal interpolation points numerically. This has been done before, but to my knowledge, there are no other open-source python libraries that do it.

## Usage

### Installation

```
pip install optinterp
```

### Example
<!---TODO: Construct a more illustrative example--->
```python
import optinterp

nds = optinterp.nodes(10)
```

This will give a marginal but non-zero improvement over using the extended Chebyshev nodes:
```python
import numpy as np

nds = np.polynomial.chebyshev.chebpts1(10)
nds = nds / nds[-1]
```

## Algorithm Description

First, to define the problem more formally, we want `n` points between -1 and 1 for which the [Lebesgue](https://en.wikipedia.org/wiki/Lebesgue_constant) constant is mimimized. The technique used here expoits three properties:  
1. Optimal interpolation points can take values -1 and 1 for their minimum and maximum.
2. To mimimize the global maximum of the Lebesgue function, all local maxima should be equal.
3. Moving two adjacent nodes closer together reduces the local maximum of the Lebesgue function at the expense of increasing the other local maxima.

Start with two initial guesses for optimal points. The extended Chebyshev nodes and a set of slightly perturbed Chebyshev nodes. Then for each set of nodes define:

```
dx_i = x_{i+1} - x_i
dL_i = L_i - L_{avg}
```
Where `L_i` is the local maximum of the Lebesgue function between `x_{i+1}` and `x_i`. Then we just use the Secant method to find roots for each `dx_i` that make `dL_i = 0`. For the next iteration, calculate your nodes `x_i` from these roots of `dx_i` and scale the values to be from -1 to 1.
