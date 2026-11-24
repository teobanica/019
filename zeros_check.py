"""Resolvents and finite pair counts from the first 80 positive zeta ordinates."""

from collections import Counter
from itertools import combinations

import mpmath as mp


def xi(s):
    return s * (s - 1) * mp.pi**(-s / 2) * mp.gamma(s / 2) * mp.zeta(s) / 2


def log_derivative_b(z):
    # z=s(s-1). At the branch point use the even expansion of xi about 1/2.
    half = mp.mpf(1) / 2
    if z == -half**2:
        return mp.diff(xi, half, 2) / (2 * xi(half))
    s = (1 + mp.sqrt(1 + 4*z)) / 2
    return mp.diff(xi, s) / ((2*s - 1) * xi(s))


def check_resolvent(ordinates):
    quarter = mp.mpf(1) / 4
    atoms = [1 / (quarter + y*y) for y in ordinates]
    l1 = 1 + mp.euler / 2 - mp.log(4*mp.pi) / 2
    for z in (-quarter, mp.mpf(1) / 3, mp.mpf(7)):
        # Subtract the value at zero to make each summand decay as gamma^(-4).
        xi_value = (l1 - log_derivative_b(z)) / z
        partial = mp.fsum(s*s / (1 + z*s) for s in atoms)
        residual = xi_value - partial
        print(f"z={mp.nstr(z, 6):>8}  xi={mp.nstr(xi_value, 12)}  "
              f"zero sum={mp.nstr(partial, 12)}  residual={mp.nstr(residual, 5)}")
        assert 0 < residual < mp.mpf('1e-7')


def check_pair_count(ordinates):
    # T is the last included ordinate. Use the paper's L_T and f(x)=exp(-x^2).
    scale = mp.log(max(ordinates)) / (2 * mp.pi)

    def kernel(delta):
        return mp.exp(-(scale * delta)**2) * 4 / (4 + delta**2)

    count = len(ordinates)
    direct = 2 * mp.fsum(kernel(a-b) for a, b in combinations(ordinates, 2)) / count
    multiplicities = Counter(ordinates)
    product_trace = mp.fsum(ma * mb * kernel(a-b)
                            for a, ma in multiplicities.items()
                            for b, mb in multiplicities.items()) / count
    # Each occurrence contributes exactly one self-pair, even at a repeated atom.
    assert mp.almosteq(direct, product_trace - 1, abs_eps=mp.mpf('1e-23'))
    return direct


def main():
    with mp.workdps(30):
        ordinates = [mp.im(mp.zetazero(k)) for k in range(1, 81)]
        print("80 ordinates through", mp.nstr(ordinates[-1], 14))
        check_resolvent(ordinates)
        print("Gaussian pair statistic =", mp.nstr(check_pair_count(ordinates), 15))
        check_pair_count([ordinates[0], ordinates[0], ordinates[1]])


if __name__ == "__main__":
    main()
