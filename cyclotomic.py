"""Gap norms for the ADE indices of Coxeter number at most 120 (Section 3)."""

from math import gcd, lcm, prod

import mpmath as mp


def gap_norms(limit):
    # Evaluating (z^n-1)/(z-1) at 1 gives n = product_{1<d|n} Phi_d(1).
    norms = {2: 2}
    for n in range(3, limit + 1):
        lower = prod(value for d, value in norms.items() if n % d == 0)
        norms[n], remainder = divmod(n, lower)
        assert remainder == 0
    return norms


def main():
    limit = 120
    norms = gap_norms(limit)
    total = 1
    multiple = 2
    with mp.workdps(50):
        for n in range(3, limit + 1):
            # One representative per conjugate pair of primitive roots.
            conjugates = (4 * mp.sin(mp.pi * a / n)**2
                          for a in range(1, (n-1)//2 + 1) if gcd(a, n) == 1)
            assert abs(mp.fprod(conjugates) - norms[n]) < mp.mpf('1e-40')
            total *= norms[n]
            multiple = lcm(multiple, n)
            assert 2 * total == multiple
    # D_{2m} has h=4m-2; E_6 and E_8 have h=12 and 30.
    exceptional = [*range(6, limit + 1, 4), 12, 30]
    assert all(norms[h] == 1 for h in exceptional)
    print("n           ", 3, 4, 5, 6, 7, 8, 9, 10)
    print("gap norm    ", *(norms[n] for n in range(3, 11)))
    print(f"Each cutoff through {limit} agrees with lcm(1,...,X)/2.")
    print(f"Product at X={limit}: {total}")


if __name__ == "__main__":
    main()
