# Fortran-Library
Mathematical & chemical routines e.g. nonlinear optimization etc. with c++ & python interface

## Featured utilities
Geometry transformation:
* Standard geometry (also called standard orientaion)
* Assimilate a geometry to a reference
* Cartesian <-> internal coordinate
* Normal mode and vibrational frequency

Nonlinear optimization: (all parameters conveniently tunable)
* Unconstrained optimizer:
Newton-Raphson, BFGS, limited memory BFGS, conjugate gradient (Dai-Yuan and Polak-Ribiere+), trust region
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

StringUtility:
* Split string
* Converts multiple spaces and tabs to single spaces
* Removes spaces, tabs, and control characters in string str
* Forked from https://gbenthien.net/strings/index.html

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
libFL.a, libFL.so, *.mod, FortranLibrary.hpp, FortranLibrary (python package directory) will be installed to your prefix (default = Fortran-Library)
1. `export LIBRARY_PATH=prefix/lib:$LIBRARY_PATH`
2. `export LD_LIBRARY_PATH=prefix/lib:$LD_LIBRARY_PATH`
3. `export CPATH=prefix/include:$CPATH`
4. `export PYTHONPATH=prefix:$PYTHONPATH`
5. Copy mkl_rci.f90 & mkl_dfti.f90 from MKL installation path to source (for gnu compiler, additionally modify `DJACOBI` in mkl_rci.f90 to specify the dimension of `x`, `fjac`, `f` explicitly: `x(n)`, `fjac(m,n)`, `f(m)`)
6. `make`, `make install`
7. (optional) `make test`, then check the logs

Optional flags to tune make:
* prefix=
* compiler= (can be intel or gnu)
* intelflag=
* gnuflag=

## Usage
* fortran `use FortranLibrary`
* c++ `#include <FortranLibrary.hpp>`
* python `import FortranLibrary as FL`

## Source
To see what this library is capable of in detail, you may open certain source file then fold all: routines are categorized and folded into different sections (VS code is recommended: press ctrl+k+0)

Source code level from bottom to top:
1. StringUtility.f90, General.f90, Mathematics.f90, LinearAlgebra.f90
2. (mkl_rci.f90, NonlinearOptimization.f90), (mkl_dfti.f90, IntegralTransform.f90), Clustering.f90, Statistics.f90, Chemistry.f90
3. GeometryTransformation.f90
4. FortranLibrary.f90

## Dependency
* MKL
* mkl_rci.f90 & mkl_dfti.f90 (can be found in MKL installation path)

## Publications based on this work
> 1. Y. Shen and D. R. Yarkony, J. Phys. Chem. A 2020, 124, 22, 4539–4548 https://doi.org/10.1021/acs.jpca.0c02763
> 2. Y. Shen and D. R. Yarkony, J. Phys. Chem. Lett. 2020, 11, 17, 7245–7252 https://doi.org/10.1021/acs.jpclett.0c02199

## Reference
> 1. J. Nocedal, S. J. Wright, *Numerical Optimization 2nd edition* (Springer, 2006)
> 2. E. B. Wilson, J. C. Decius, P. C. Cross, *Molecular viobrations: the theory of infrared and Raman vibrational spectra* (Dover, 1980)
> 3. D. R. Yarkony, J. Chem. Phys. 112, 2111 (2000)