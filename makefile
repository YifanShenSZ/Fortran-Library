################################################################
#                                                              #
#        Makefile for a test program on Fortran-Library        #
#                                                              #
################################################################

compiler = ifort
MyLibDir = .
MyLib = $(MyLibDir)/General.f90 $(MyLibDir)/Mathematics.f90 $(MyLibDir)/LinearAlgebra.f90 \
$(MyLibDir)/mkl_rci.f90 $(MyLibDir)/NonlinearOptimization.f90 $(MyLibDir)/mkl_dfti.f90 $(MyLibDir)/IntegralTransform.f90 \
$(MyLibDir)/Clustering.f90 $(MyLibDir)/Statistics.f90 $(MyLibDir)/Chemistry.f90 \
$(MyLibDir)/GeometryTransformation.f90
src = test.f90
exe = test.exe
flag = -m64 -xCORE-AVX2 -mtune=core-avx2 -mkl -ipo -O3 -no-prec-div -fp-model fast=2

$(exe): $(MyLib) $(src)
	$(compiler) $(flag) $^ -o $(exe)

clean:
	rm *.mod