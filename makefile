##############################################
#                                            #
#        Makefile for Fortran-Library        #
#                                            #
##############################################

# Default to install to Fortran-Library
prefix = .
# Can be intel or gnu
compiler = gnu
intelflag = -ipo -m64 -xCORE-AVX2 -mtune=core-avx2 -no-prec-div -fp-model fast=2 -parallel -O3 -mkl
gnuflag = -ffree-line-length-0 -fno-range-check -m64 -march=core-avx2 -mtune=core-avx2 -O3 -I${MKLROOT}/include
gnumkl = -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_gf_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl
src = $(addprefix source/, General.f90 Mathematics.f90 LinearAlgebra.f90 \
mkl_rci.f90 NonlinearOptimization.f90 mkl_dfti.f90 IntegralTransform.f90 \
Clustering.f90 Statistics.f90 Chemistry.f90 \
GeometryTransformation.f90 \
FortranLibrary.f90 )
RealPrefix = $(realpath $(prefix))
incdir = $(RealPrefix)/include
libdir = $(RealPrefix)/lib

libFL.a libFL.so: $(src)
ifeq ($(compiler),intel)
	ifort $(intelflag) -c $^
	xiar rcs libFL.a *.o
	ifort $(intelflag) -shared -fpic $^ -o libFL.so
else
	gfortran $(gnuflag) -c $^
	ar rcs libFL.a *.o
	gfortran $(gnuflag) -shared -fpic $^ -o libFL.so
endif
	rm *.o

.PHONY: install
install: | $(incdir) $(libdir)
	mv *.mod $(incdir)
ifneq ($(realpath include),$(incdir))
	cp include/*.h $(incdir)
endif
	mv *.a  $(libdir)
	mv *.so $(libdir)
ifneq ($(realpath .),$(RealPrefix))
	cp -r FortranLibrary $(RealPrefix)
endif

$(incdir):
	mkdir $(incdir)

$(libdir):
	mkdir $(libdir)

.PHONY: test
test:
ifeq ($(compiler),intel)
	ifort $(intelflag) -I$(incdir) test/test.f90 $(libdir)/libFL.a -o test/test_static.exe
else
	gfortran $(gnuflag) -I$(incdir) test/test.f90 -l:libFL.a $(gnumkl) -o test/test_static.exe
endif
	test/test_static.exe > test/log_static

ifeq (,$(findstring $(libdir),$(LIBRARY_PATH)))
$(error Please add prefix/lib to LIBRARY_PATH)
endif
ifeq (,$(findstring $(libdir),$(LD_LIBRARY_PATH)))
$(error Please add prefix/lib to LD_LIBRARY_PATH)
endif

ifeq ($(compiler),intel)
	ifort $(intelflag) test/test.f90 -lFL -o test/test_dynamic.exe
else
	gfortran $(gnuflag) -I$(incdir) test/test.f90 -lFL $(gnumkl) -o test/test_dynamic.exe
endif
	test/test_dynamic.exe > test/log_dynamic

ifeq (,$(findstring $(incdir),$(CPATH)))
$(error Please add prefix/include to CPATH)
endif
ifeq ($(compiler),intel)
	icpc $(intelflag) test/test.cpp -lFL -o test/test_cpp.exe
else
	g++ $(gnuflag) test/test.cpp -lFL -o test/test_cpp.exe
endif
	test/test_cpp.exe > test/log_cpp

ifeq (,$(findstring $(RealPrefix),$(PYTHONPATH)))
$(error Please add prefix to PYTHONPATH)
endif
	python test/test.py > test/log_py
