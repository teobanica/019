# Code for *Quantum algebra thoughts on the Riemann hypothesis*

These five Python programs reproduce finite calculations in the paper. Run them from the extracted `code019` folder with Python 3.9 or later. All inputs are specified in the scripts, including the number of zeros. The programs print their results and check them as they run.

Install the dependencies once if needed.

```sh
python3 -m pip install mpmath numpy scipy
```

| Program | Where it is used | Calculation |
| --- | --- | --- |
| `verify.py` | Sections 1 and 2, especially Proposition 2.7 and Theorem 2.11 | Computes the first three Li coefficients from the Taylor expansion at 1, checks their ratios and the signed trace obstruction, and evaluates the trace inequality on a four-point magic unitary. |
| `li_contour.py` | Sections 1 and 2 | Uses Claude's independent contour method from the referee report. It computes seven Li coefficients at two radii and compares the first three with `verify.py`. |
| `zeros_check.py` | Theorem 2.8 and the pair statistic in Section 4 | Compares the differentiated xi function with three truncated resolvent sums using the first 80 positive ordinates. Computes the finite Gaussian pair statistic by two methods and checks the diagonal correction with a repeated input. |
| `cyclotomic.py` | Theorem 3.5 | Checks the gap norms using conjugate products and exact integer cyclotomic evaluations. Compares their cumulative product with half the least common multiple at every integer cutoff from 3 to 120 and checks the additional D and E factors. |
| `powers.py` | Section 4, Propositions 4.9 and 4.13 | Integrates the classical trace moments on finite root-of-unity grids for matrix sizes 2, 3 and 4 and checks the compressed flip on the quantum antisymmetric square. |

```sh
python3 verify.py
python3 li_contour.py
python3 zeros_check.py
python3 cyclotomic.py
python3 powers.py
```

The two Li calculations give approximately −0.00007405759501464 for the combination in Theorem 2.11. At the first 80 ordinates, the three regularized resolvent residuals are about 2.47 × 10⁻⁸ and the Gaussian pair statistic is 0.157510531840443. This statistic uses the paper's global scale and Montgomery weight. It is a finite-height illustration.
