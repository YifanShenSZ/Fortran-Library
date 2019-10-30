##############################################
#                                            #
#        Makefile for Fortran-Library        #
#                                            #
##############################################

compiler = ifort
src = General.f90 Mathematics.f90 LinearAlgebra.f90 \
mkl_rci.f90 NonlinearOptimization.f90 mkl_dfti.f90 IntegralTransform.f90 \
Clustering.f90 Statistics.f90 Chemistry.f90 \
GeometryTransformation.f90
flag = -m64 -xCORE-AVX2 -mtune=core-avx2 -mkl -O3 -no-prec-div -fp-model fast=2

make: $(src)
ifeq ($(compiler),ifort)
	ifort $(flag) -c $^
	xiar rc libFL.a *.o
	ifort $(flag) -shared -fpic $^ -o libFL.so
else
	gfortran $(flag) -c $^
	ar cr libFL.a *.o
	gfortran $(flag) -shared -pic $^ -o libFL.so
endif
	rm *.o

install:
	mv *.mod /usr/include
	mv *.so /usr/lib
	mv *.a /usr/lib

test:
	$(compiler) $(flag) -ipo test.f90 libFL.a -o test_static.exe
	./test_static.exe > log_static
	$(compiler) $(flag) -ipo test.f90 libFL.so -o test_dynamic.exe
	./test_dynamic.exe > log_dynamic