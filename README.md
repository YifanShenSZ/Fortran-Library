# Fortran-Library
Mathematical & chemical routines like nonlinear optimization etc. with c++ & python interface

## Featured utilities
Geometry transformation:
* Standard geometry (also called standard orientaion)
* Assimilate a geometry to a reference
* Cartesian <-> internal coordinate
* Normal mode and vibrational frequency

Nonlinear optimization: (all parameters conveniently tunable)
* Unconstrained optimizer:
Newton-Raphson, BFGS, limited memory BFGS, conjugate gradient (Dai-Yun and Polak-Ribiere+), trust region
* Constrained optimizer: augmented Lagrangian

Integral transform:
* Fourier transform
* Fast Fourier transform

Clustering:
* K-means
* Gaussian mixture model

Statistics:
* Variance, R^2
* Distribution

Chemistry:
* Gradient of eigenstates
* Phase fixing
* Conical intersection adapted coordinate

General:
* Random number
* Sorting
* Some other basic routines

Mathematics:
* Selected special functions and combinatorics
* Quaternion operation
* Ordinary differential equation
* 1 dimensional integration

Linear algebra:
* All kinds of vector & matrix & tensor operation
* LAPACK wrapper for linear solver, eigensystem, matrix norm

## Installation
1. Copy mkl_rci.f90 & mkl_dfti.f90 from MKL installation path to source
2. `make`, `make install`
3. `export LIBRARY_PATH=Fortran-Library/lib:$LIBRARY_PATH`
4. `export LD_LIBRARY_PATH=Fortran-Library/lib:$LD_LIBRARY_PATH`
5. `export CPATH=Fortran-Library/include:$CPATH`
6. `export PYTHONPATH=Fortran-Library:$PYTHONPATH`
7. (optional) `make test`, then check the logs

## Usage
* fortran `use FortranLibrary`
* c++ `#include "FortranLibrary.h"`
* python `import FortranLibrary`

## Source
To see what this library is capable of in detail, you may open certain source file then fold all: routines are categorized and folded into different sections (VS code is recommended: press ctrl+k+0)

Source code level from bottom to top:
1. General, Mathematics, LinearAlgebra
2. (mkl_rci, NonlinearOptimization), (mkl_dfti, IntegralTransform), Clustering, Statistics, Chemistry
3. GeometryTransformation

## Dependency
* MKL
* mkl_rci.f90 & mkl_dfti.f90 (can be found in MKL installation path)

## Reference
> 1. J. Nocedal, S. J. Wright, *Numerical Optimization 2nd edition* (Springer, 2006)
> 2. E. B. Wilson, J. C. Decius, P. C. Cross, *Molecular viobrations: the theory of infrared and Raman vibrational spectra* (Dover, 1980)
> 3. D. R. Yarkony, J. Chem. Phys. 112, 2111 (2000)