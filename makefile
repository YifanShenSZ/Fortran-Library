##############################################
#                                            #
#        Makefile for Fortran-Library        #
#                                            #
##############################################

prefix = . # Default to install to Fortran-Library
f90 = ifort
cpp = icpc
flag = -m64 -xCORE-AVX2 -mtune=core-avx2 -no-prec-div -fp-model fast=2 -mkl -parallel -O3
src = $(addprefix source/, General.f90 Mathematics.f90 LinearAlgebra.f90 \
mkl_rci.f90 NonlinearOptimization.f90 mkl_dfti.f90 IntegralTransform.f90 \
Clustering.f90 Statistics.f90 Chemistry.f90 \
GeometryTransformation.f90 \
FortranLibrary.f90 )
RealPrefix = $(realpath $(prefix))
incdir = $(RealPrefix)/include
libdir = $(RealPrefix)/lib

libFL.a libFL.so: $(src)
ifeq ($(f90),ifort)
	ifort $(flag) -c $^
	xiar cr libFL.a *.o
	ifort $(flag) -shared -fpic $^ -o libFL.so
else
	gfortran $(flag) -c $^
	ar cr libFL.a *.o
	gfortran $(flag) -shared -pic $^ -o libFL.so
endif
	rm *.o

.PHONY: install
install: | $(incdir) $(libdir)
	mv *.mod $(incdir)
ifneq ($(realpath include),$(incdir))
	cp include/*.h  $(incdir)
	cp include/*.py $(incdir)
endif
	mv *.a   $(libdir)
	mv *.so  $(libdir)

$(incdir):
	mkdir $(incdir)

$(libdir):
	mkdir $(libdir)

.PHONY: test
test:
	$(f90) $(flag) -ipo -I$(incdir) test/test.f90 $(libdir)/libFL.a -o test/test_static.exe
	test/test_static.exe > test/log_static
ifneq (,$(findstring $(libdir),$(LIBRARY_PATH)))
ifneq (,$(findstring $(libdir),$(LD_LIBRARY_PATH)))
	$(f90) $(flag) -ipo -lFL test/test.f90 -o test/test_dynamic.exe
	test/test_dynamic.exe > test/log_dynamic
ifneq (,$(findstring $(incdir),$(CPATH)))
	$(cpp) $(flag) -ipo -lFL test/test.cpp -o test/test_cpp.exe
	test/test_cpp.exe > test/log_cpp
else
	# Please set CPATH properly before running test
endif
ifneq (,$(findstring $(incdir),$(PYTHONPATH)))
	python test/test.py > test/log_py
else
	# Please set PYTHONPATH properly before running test
endif
else
	# Please set LD_LIBRARY_PATH properly before running test
endif
else
	# Please set LIBRARY_PATH properly before running test
endif
