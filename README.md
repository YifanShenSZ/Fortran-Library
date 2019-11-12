# Fortran-Library
Nonlinear optimization, clustering, and more mathematics & chemistry

## Featured utilities
Nonlinear optimization: (all parameters conveniently tunable)
* Unconstrained optimizer:
Newton-Raphson, BFGS, limited memory BFGS, conjugate gradient (Dai-Yun and Polak-Ribiere+), trust region
* Constrained optimizer: augmented Lagrangian

Clustering:
* K-means
* Gaussian mixture model

Mathematics:
* Combinatorics and selected special functions
* Ordinary differential equation
* 1 dimensional integration

Integral transform:
* Fourier transform
* Fast Fourier transform

Linear algebra:
* All kinds of vector & matrix & tensor operation
* LAPACK wrapper for linear solver, eigensystem, matrix norm

Chemistry:
* Gradient of eigenstates
* Phase fixing
* Conical intersection adapted coordinate

Geometry transformation:
* Standard geometry (also called standard orientaion)
* Cartesian <-> internal coordinate
* Normal mode and vibrational frequency

General:
* Random number
* Sorting
* Some other basic routines

## Installation
1. make, make install
2. (optional) make test, then check the log

## Source
To see what this library is capable of in detail, you may open certain source file then fold all: routines are categorized and folded into different sections (VS code is recommended: press ctrl+k+0)

Source code level from bottom to top:
1. General, Mathematics, LinearAlgebra
2. (mkl_rci, NonlinearOptimization), (mkl_dfti, IntegralTransform), Clustering, Statistics, Chemistry
3. GeometryTransformation

test is to build a demo program, testing the functionality

## Dependency
* This library depends on BLAS & LAPACK & MKL
* MKL reverse communication interface is adopted in nonlinear optimization, MKL discrete Fourier transform interface is adopted in integral transform, so this library needs to be compiled together with mkl_rci and mkl_dfti. They can be found in your MKL installation path

## Reference
> 1. J. Nocedal, S. J. Wright, *Numerical Optimization 2nd edition* (Springer, 2006)
> 2. E. B. Wilson, J. C. Decius, P. C. Cross, *Molecular viobrations: the theory of infrared and Raman vibrational spectra* (Dover, 1980)
> 3. D. R. Yarkony, J. Chem. Phys. 112, 2111 (2000)