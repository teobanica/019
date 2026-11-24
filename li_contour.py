"""Independent contour calculation adapted from Claude's referee script."""

from math import factorial, log, pi

import mpmath as mp
import numpy as np
from scipy.special import digamma

from verify import first_three_li


def li_by_contour(count, radius, samples=512):
    # Borwein's accelerated Dirichlet eta series and its derivative.
    terms = 60
    running_sum = 0
    partial_weights = []
    for j in range(terms + 1):
        numerator = terms * factorial(terms + j - 1) * 4**j
        denominator = factorial(terms - j) * factorial(2 * j)
        running_sum += numerator // denominator
        partial_weights.append(running_sum)
    indices = np.arange(terms)
    weights = (-1.)**indices * (1 - np.array(partial_weights[:-1], dtype=float) /
                                float(partial_weights[-1]))
    log_integers = np.log(indices + 1)
    circle = radius * np.exp(2j * pi * np.arange(samples) / samples)
    s = 1 / (1 - circle)
    powers = np.exp(-s[:, None] * log_integers)
    eta = powers @ weights
    eta_prime = powers @ (-weights * log_integers)
    two = np.exp((1 - s) * log(2))
    zeta_log_derivative = eta_prime / eta - log(2) * two / (1 - two)
    xi_log_derivative = (1 / s + 1 / (s - 1) - log(pi) / 2
                         + digamma(s / 2) / 2 + zeta_log_derivative)
    # d/dz log xi(1/(1-z)) = sum_{n>=1} lambda_n z^(n-1).
    values = xi_log_derivative / (1 - circle)**2
    return np.fft.fft(values)[:count] / (samples * radius**np.arange(count))


def main():
    coefficients = li_by_contour(7, 0.5)
    second_radius = li_by_contour(7, 0.6)
    np.testing.assert_allclose(coefficients, second_radius, rtol=0, atol=2e-12)
    with mp.workdps(45):
        taylor = np.array([float(value) for value in first_three_li()])
    np.testing.assert_allclose(coefficients[:3], taylor, rtol=0, atol=5e-14)
    for n, value in enumerate(coefficients.real, start=1):
        print(f"lambda_{n} = {value:.13f}")
    l1, l2, l3 = coefficients.real[:3]
    print(f"S4+ obstruction from the contour = {l3 - 4*l2 + 7*l1:.13f}")


if __name__ == "__main__":
    main()
