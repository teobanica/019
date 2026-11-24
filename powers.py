"""Finite character and power-trace checks for Section 4."""

from itertools import combinations
from math import prod

import numpy as np


def check_cue(n):
    # Weyl integration on the eigenangles. For k<=2n each exponent in
    # |Delta|^2 |p_k|^2 has absolute value <= 3n-1. A grid of 3n roots
    # integrates these Laurent polynomials exactly in exact arithmetic.
    grid_size = 3*n
    roots = np.exp(2j * np.pi * np.arange(grid_size) / grid_size)
    frequencies = np.arange(1, 2*n + 1)
    mass = 0.
    moments = np.zeros(len(frequencies))
    for points in combinations(roots, n):
        # Unordered subsets absorb the 1/n! in the Weyl density.
        weight = prod(abs(a-b)**2 for a, b in combinations(points, 2)) / grid_size**n
        power_sums = (np.array(points)[:, None]**frequencies).sum(axis=0)
        mass += weight
        moments += weight * abs(power_sums)**2
    np.testing.assert_allclose(mass, 1, rtol=0, atol=1e-12)
    np.testing.assert_allclose(moments, np.minimum(frequencies, n), rtol=0, atol=1e-12)
    print(f"CUE N={n}, moments to k={2*n}: {moments.round(8)}")


def check_flip():
    # On each span of e_i x e_j, e_j x e_i, the flip interchanges the basis.
    flip = np.array([[0., 1.], [1., 0.]])
    for q in (0.1, 0.35, 0.7, 0.95):
        antisymmetric = np.array([1., -q]) / np.sqrt(1 + q*q)
        compression = antisymmetric @ flip @ antisymmetric
        np.testing.assert_allclose(compression, -2*q/(1+q*q))
        assert 1 + compression**2 < 2
        print(f"q={q:.2f}, compressed flip={compression:.8f}, "
              f"ordered second-moment upper bound={1+compression**2:.8f}")


if __name__ == "__main__":
    for size in (2, 3, 4):
        check_cue(size)
    check_flip()
