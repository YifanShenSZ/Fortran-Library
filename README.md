# Fortran-Library
Nonlinear optimization, linear algebra, and more 

Selected utilities:
* Nonlinear optimization: (all parameters conveniently tunable)
1. Unconstrained optimizer:
Newton-Raphson, BFGS, limited memory BFGS, conjugate gradient (Dai-Yun and Polak-Ribiere+), trust region
2. Constrained optimizer: augmented Lagrangian
* Linear algebra:
1. All kinds of vector & matrix & tensor operation
2. LAPACK wrapper for linear solver, eigensystem, matrix norm
* Geometry transformation:
1. Standard geometry (also called standard orientaion)
2. Cartesian <-> internal coordinate
3. Normal mode and vibrational frequency

To see what this library is capable of in detail, you may open certain source file and simply fold all: routines are categorized and folded into different sections (VS code is recommended: press ctrl+k+0)

Dependency:
* This library depends on BLAS & LAPACK & MKL
* Need to compile together with MKL_RCI.f90, which can be found in your MKL installation path

test.f90 and makefile are to build a demo program, testing the functionality

Reference:
> 1. J. Nocedal, S. J. Wright, *Numerical Optimization 2nd edition* (Springer, 2006)
> 2. E. B. Wilson, J. C. Decius, P. C. Cross, *Molecular viobrations: the theory of infrared and Raman vibrational spectra* (Dover, 1980)