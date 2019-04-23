# Fortran-Library
Nonlinear optimization, linear algebra, and more 

Selected utilities:
1. Nonlinear optimization: (all parameters conveniently tunable)
* Unconstrained optimizer:
Newton-Raphson, BFGS, limited memory BFGS, conjugate gradient (Dai-Yun and Polak-Ribiere+), trust region
* Constrained optimizer: augmented Lagrangian
2. Linear algebra:
* All kinds of vector & matrix & tensor operation
* LAPACK wrapper for linear solver, eigensystem, matrix norm
3. Geometry transformation:
* Standard geometry (also called standard orientaion)
* Cartesian <-> internal coordinate
* Normal mode and vibrational frequency
4. Nonadiabatic:
* Gradient of eigenstates
* Phase fixing
* Conical intersection adapted coordinate

To see what this library is capable of in detail, you may open certain source file and simply fold all: routines are categorized and folded into different sections (VS code is recommended: press ctrl+k+0)

Dependency:
* This library depends on BLAS & LAPACK & MKL
* MKL reverse communication interface is adopted, so need to compile together with MKL_RCI.f90, which can be found in your MKL installation path

test.f90 and makefile are to build a demo program, testing the functionality

Reference:
> 1. J. Nocedal, S. J. Wright, *Numerical Optimization 2nd edition* (Springer, 2006)
> 2. E. B. Wilson, J. C. Decius, P. C. Cross, *Molecular viobrations: the theory of infrared and Raman vibrational spectra* (Dover, 1980)