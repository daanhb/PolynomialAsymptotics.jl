# PolynomialAsymptotics
A package that contains expressions for the asymptotic expansions of orthogonal polynomials of large degree

To run the Fortran code for Gauss-Laguerre, one first obtains the code that accompanies [J. Bremer. ``On the numerical calculation of the roots of special functions satisfying second order ordinary differential equations''. SIAM J. Sci. Comput., 39(1):A55–A82, 2016] which has been moved to https://github.com/JamesCBremerJr/Phase-functions. Next, follow the instructions in fortran/explExpa.f90.

The path to the chebfun library (chebfun.org) needs to be included for the Gauss quadrature tests of the matlab/ folder. These also refer to the JACOBI_v3 and Laguerre3 folders, which construct asymptotic expansions of Jacobi-type and Laguerre-type orthogonal polynomials.

The sage/ folder contains Sage worksheets (see sagemath.org) which construct asymptotic expansions of nodes and weights of Gaussian quadrature rules. The procedure is explained in detail in sage/procedureForAsymptotics.pdf. 
