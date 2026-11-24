"""Li coefficients from the Taylor expansion of log xi at 1 (Sections 1 and 2)."""

import mpmath as mp
import numpy as np


def first_three_li():
    """Use the user's current mpmath precision; return lambda_1 through lambda_3."""
    half = mp.mpf(1) / 2
    euler = mp.euler
    # c_k is the coefficient of t^k in log xi(1+t).
    c1 = 1 + euler / 2 - mp.log(4 * mp.pi) / 2
    c2 = -half + mp.polygamma(1, half) / 8 - mp.stieltjes(1) - euler**2 / 2
    c3 = (mp.mpf(1) / 3 + mp.polygamma(2, half) / 48
          + mp.stieltjes(2) / 2 + euler * mp.stieltjes(1) + euler**3 / 3)
    return c1, 2 * (c1 + c2), 3 * (c1 + 2 * c2 + c3)


def check_magic_example():
    # Each flip block has the form [[p,1-p],[1-p,p]].
    # Two such blocks give a four-point magic unitary with X=2p+2r.
    identity = np.eye(2)
    p = np.diag([1., 0.])
    vector = np.array([4., 3.]) / 5
    r = np.outer(vector, vector)
    assert not np.allclose(p @ r, r @ p)
    q = 4 * identity - 2 * (p + r)
    mean = np.trace(q) / 2
    excess = np.trace(q @ q) / 2 - 2 * mean
    np.testing.assert_allclose(excess, 4 * vector[0]**2)
    print(f"Block magic example: tr(Q^2) - 2 tr(Q) = {excess:.6f}")


def main():
    with mp.workdps(45):
        l1, l2, l3 = first_three_li()
        for n, value in enumerate((l1, l2, l3), start=1):
            print(f"lambda_{n} = {mp.nstr(value, 26)}")
        print("lambda_2/lambda_1 =", mp.nstr(l2 / l1, 18))
        print("lambda_3/lambda_1 =", mp.nstr(l3 / l1, 18))
        obstruction = l3 - 4 * l2 + 7 * l1
        print("lambda_3 - 4 lambda_2 + 7 lambda_1 =", mp.nstr(obstruction, 25))
        assert 3 < l2 / l1 < 4 and l3 / l1 < 9
        printed_value = mp.mpf("-0.0000740575950")
        assert abs(obstruction - printed_value) < mp.mpf("0.5e-13")
    check_magic_example()


if __name__ == "__main__":
    main()
