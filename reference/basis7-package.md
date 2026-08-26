# basis7: An S7 Framework for Basis Expansions

Basis expansions as objects, built on the S7 object-oriented system. A
basis carries an interval and a number of functions, and answers with
its design matrix at any points, its derivative of any order, its
integral anchored at the lower endpoint, and the Gram matrix of any
derivative against a chosen measure, which is the matrix of a roughness
penalty. B-spline, Fourier and orthogonal polynomial families ship with
exact formulas; one wrapper orthonormalizes a basis, imposes linear
constraints or rotates it to the Demmler-Reinsch form, and a tensor
product builds a basis of several variables. Numerical fallbacks make a
user-defined basis work from its evaluation alone, and a diagnostic
verifies every component against an independent route.

## See also

Useful links:

- <https://statmodels7.github.io/basis7/>

- <https://github.com/statmodels7/basis7>

- Report bugs at <https://github.com/statmodels7/basis7/issues>

## Author

**Maintainer**: Giovanni Tinervia <giovannitinervia9@gmail.com>
