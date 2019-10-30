##############################################
#                                            #
#        Makefile for Fortran-Library        #
#                                            #
##############################################

compiler = ifort
flag = -m64 -xCORE-AVX2 -mtune=core-avx2 -mkl -O3 -no-prec-div -fp-model fast=2
prefix = /usr/local
src = General.f90 Mathematics.f90 LinearAlgebra.f90 \
mkl_rci.f90 NonlinearOptimization.f90 mkl_dfti.f90 IntegralTransform.f90 \
Clustering.f90 Statistics.f90 Chemistry.f90 \
GeometryTransformation.f90

libFL.a libFL.so: $(src)
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

.PHONY: install
install:
	mv *.mod $(prefix)/include
	mv *.a   $(prefix)/lib
	mv *.so  $(prefix)/lib

.PHONY: test
test:
	$(compiler) $(flag) -I$(prefix)/include -ipo test.f90 $(prefix)/lib/libFL.a -o test_static.exe
	./test_static.exe > log_static
	$(compiler) $(flag) -I$(prefix)/include -L$(prefix)/lib -lFL -ipo test.f90 -o test_dynamic.exe
ifeq (,$(findstring $(prefix)/lib,$(LD_LIBRARY_PATH)))
	# Please set LD_LIBRARY_PATH properly, then run test_dynamic.exe
else
	./test_dynamic.exe > log_dynamic
endif

.PHONY: clean
clean:
	rm *.mod
	rm *.a
	rm *.so
	rm test_*
	rm log_*