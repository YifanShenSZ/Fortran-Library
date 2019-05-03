################################################################
#                                                              #
#        Makefile for a test program on Fortran-Library        #
#                                                              #
################################################################

compiler = ifort
MyLibDir = .
MyLib = $(MyLibDir)/General.f90 $(MyLibDir)/Mathematics.f90 $(MyLibDir)/LinearAlgebra.f90 $(MyLibDir)/MKL_RCI.f90 $(MyLibDir)/NonlinearOptimization.f90 $(MyLibDir)/GeometryTransformation.f90 $(MyLibDir)/Chemistry.f90
src = test.f90
exe = test.exe
flag = -u -mkl -fast -march=core-avx2

$(exe): $(MyLib) $(src)
	$(compiler) $(flag) $^ -o $(exe)

clean:
	rm $(exe)
	rm *.mod
