# Fortran-Library
Nonlinear optimization, clustering, and more mathematics & chemistry

Featured utilities:
1. Nonlinear optimization: (all parameters conveniently tunable)
* Unconstrained optimizer:
Newton-Raphson, BFGS, limited memory BFGS, conjugate gradient (Dai-Yun and Polak-Ribiere+), trust region
* Constrained optimizer: augmented Lagrangian
2. Clustering:
* K-means
* Gaussian mixture model
3. Mathematics:
* Combinatorics and selected special functions
* Ordinary differential equation
* 1 dimensional integration
4. Integral transform:
* Fourier transform
* Fast Fourier transform
5. Linear algebra:
* All kinds of vector & matrix & tensor operation
* LAPACK wrapper for linear solver, eigensystem, matrix norm
6. Chemistry:
* Gradient of eigenstates
* Phase fixing
* Conical intersection adapted coordinate
7. Geometry transformation:
* Standard geometry (also called standard orientaion)
* Cartesian <-> internal coordinate
* Normal mode and vibrational frequency
8. General:
* Random number
* Sorting
* Some other basic routines

To see what this library is capable of in detail, you may open certain source file and simply fold all: routines are categorized and folded into different sections (VS code is recommended: press ctrl+k+0)

Dependency:
* This library depends on BLAS & LAPACK & MKL
* MKL reverse communication interface is adopted in nonlinear optimization, MKL discrete Fourier transform interface is adopted in integral transform, so this library needs to be compiled together with mkl_rci and mkl_dfti. They can be found in your MKL installation path

Source code level from bottom to top:
1. General, Mathematics, LinearAlgebra
2. (mkl_rci, NonlinearOptimization), (mkl_dfti, IntegralTransform), Clustering, Statistics, Chemistry
3. GeometryTransformation

test.f90 and makefile are to build a demo program, testing the functionality

Reference:
> 1. J. Nocedal, S. J. Wright, *Numerical Optimization 2nd edition* (Springer, 2006)
> 2. E. B. Wilson, J. C. Decius, P. C. Cross, *Molecular viobrations: the theory of infrared and Raman vibrational spectra* (Dover, 1980)