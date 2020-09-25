// c++ interface to Fortran-Library
// Include this header so that C++ can link to libFL.so

#ifndef FortranLibrary_hpp
#define FortranLibrary_hpp

#include "NonlinearOptimization.hpp"
#include "GeometryTransformation.hpp"

namespace FL {

// No need to fetch constant, directly provide them in c++
// Unit conversion (Multiplying xIny converts x to y)
#define AMUInAU  = 1822.888486192
#define AInAU    = 1.8897261339212517
#define cm_1InAu = 0.000004556335830019422

// C recognizes fortran functions by their compiled names
// fortran function Func contained in module Mod will be renamed
//     ifort: mod_mp_func_
//     gfortran: __mod_MOD_func
// ( You may use command nm to view the compiled names: nm libFL.so )
// For convenience, rename fortran functions back to their origin
#ifdef __INTEL_COMPILER
extern "C" {
    // General
    void general_mp_dscientificnotation_(double & x, int & i);
    // LinearAlgebra
    void linearalgebra_mp_my_dgemm_t_(double * A, double * B, double * C, const int & M, const int & K, const int & N);
    void linearalgebra_mp_my_dsyev_(const char & jobtype, double * A, double * eigval, const int & N);
    // Chemistry
    void chemistry_mp_initializephasefixing_(const int & NStates);
}
namespace General {
    inline void dScientificNotation(double & x, int & i) {general_mp_dscientificnotation_(x, i);}
}
namespace LA { // LinearAlgebra
    inline void My_dgemm_T(double * A, double * B, double * C, const int & M, const int & K, const int & N) {
        linearalgebra_mp_my_dgemm_t_(A, B, C, M, K, N);
    }
    inline void My_dsyev(const char & jobtype, double * A, double * eigval, const int & N) {
        linearalgebra_mp_my_dsyev_(jobtype, A, eigval, N);
    }
}
namespace Chemistry {
    inline void InitializePhaseFixing(const int & NStates) {chemistry_mp_initializephasefixing_(NStates);}
}
#elif __GNUC__
extern "C" {
    // General
    void __general_MOD_dscientificnotation(double & x, int & i);
    // LinearAlgebra
    void __linearalgebra_MOD_my_dgemm_t(double * A, double * B, double * C, const int & M, const int & K, const int & N);
    void __linearalgebra_MOD_my_dsyev(const char & jobtype, double * A, double * eigval, const int & N);
    // Chemistry
    void __chemistry_MOD_initializephasefixing(const int & NStates);
}
namespace General {
    inline void dScientificNotation(double & x, int & i) {__general_MOD_dscientificnotation(x, i);}
}
namespace LA { // LinearAlgebra
    inline void My_dgemm_T(double * A, double * B, double * C, const int & M, const int & K, const int & N) {
        __linearalgebra_MOD_my_dgemm_t(A, B, C, M, K, N);
    }
    inline void My_dsyev(const char & jobtype, double * A, double * eigval, const int & N) {
        __linearalgebra_MOD_my_dsyev(jobtype, A, eigval, N);
    }
}
namespace Chemistry {
    inline void InitializePhaseFixing(const int & NStates) {__chemistry_MOD_initializephasefixing(NStates);}
}
#endif

} // namespace FL

#endif