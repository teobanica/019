# Quantum algebra thoughts on the Riemann hypothesis

**Chat and Teo** (November 2026)

This repository contains the paper (22 pages, PDF and LaTeX source) and the 5 Python programs used for its computations.

## Summary

Li's criterion expresses the Riemann hypothesis through the signs of a sequence of real numbers attached to the zeta function. We ask what happens when we read these numbers as averages of squared characters, first on the rotation groups $SU_2$ and $SO_3$, then on their quantum counterparts.

- **Quantum rotations.** The Riemann hypothesis is equivalent to the existence of a central tracial state on $O_2^+$ whose squared character moments are the normalized Li coefficients. This state is unique. We also examine the corresponding moment problem on $S_4^+$.
- **Jones indices.** For the discrete indices $I_n=4\cos^2(\pi/n)$ below $4$, the field norm of $4-I_n$ equals $p$ when $n$ is a power of the prime $p$, and equals $1$ otherwise. Taking products gives a reformulation of von Koch's criterion.
- **Random matrices.** The character polynomial of the quantum unitary group $U_N^q$ has the same law as the characteristic polynomial of a Haar unitary matrix. This gives a quantum group formulation of Montgomery's pair correlation conjecture for the zeros of the zeta function.

## Contents

[019.pdf](019.pdf) is the paper and [019.tex](019.tex) is its LaTeX source. The programs are in [code019/](code019/).

| Programs | Computation |
| --- | --- |
| `verify.py`, `li_contour.py` | Compute Li coefficients by Taylor expansion and contour integration, and check the trace inequality for quantum permutations. |
| `zeros_check.py` | Compare sums over the first 80 positive zero ordinates with a logarithmic derivative of the completed zeta function, and compute a Gaussian pair statistic. |
| `cyclotomic.py` | Compute the field norms of the index gaps and check their product against the least common multiple formula through 120. |
| `powers.py` | Check power-trace moments for Haar unitary matrices in dimensions 2, 3 and 4, and a flip calculation for the deformed quantum unitary group. |

## Running the code

From the repository root, run the following with Python 3.9 or later.

```sh
cd code019
python3 -m pip install mpmath numpy scipy
python3 verify.py
python3 li_contour.py
python3 zeros_check.py
python3 cyclotomic.py
python3 powers.py
```
