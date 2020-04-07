// c++ interface to Fortran-Library
// Include this header so that C++ can link to libFL.so
#ifndef FortranLibrary_hpp
#define FortranLibrary_hpp

namespace FL {
// No need to fetch constant, directly provide them in c++
const double AMUInAU  = 1822.888486192;
const double AInAU    = 1.8897261339212517;
const double cm_1InAu = 0.000004556335830019422;

// C recognizes fortran functions by their compiled names
// fortran function Func contained in module Mod will be renamed
//     ifort: mod_mp_func_
//     gfortran: __mod_MOD_func
// ( You may use command nm to view the compiled names: nm libFL.so )
// For convenience, rename fortran functions back to their origin
#ifdef __INTEL_COMPILER
    extern "C" {
        // General
        void general_mp_showtime_();
        void general_mp_dscientificnotation_(double & x, int32_t & i);
        // Chemistry
        void chemistry_mp_initializephasefixing_(const int32_t & NStates);
        // GeometryTransformation
        int32_t geometrytransformation_mp_defineinternalcoordinate_(const char * format, const char * file, int32_t len_format, int32_t len_file);
        void geometrytransformation_mp_sinternalcoordinate_(const float * r, float * q, const int32_t & cartdim, const int32_t & intdim);
        void geometrytransformation_mp_dinternalcoordinate_(const double * r, double * q, const int32_t & cartdim, const int32_t & intdim);
    }
    namespace General {
        inline void ShowTime() {general_mp_showtime_();}
        inline void dScientificNotation(double & x, int32_t & i) {general_mp_dscientificnotation_(x, i);}
    }
    namespace Chemistry {
        inline void InitializePhaseFixing(const int32_t & NStates) {chemistry_mp_initializephasefixing_(NStates);}
    }
    namespace GeometryTransformation {
        inline int32_t DefineInternalCoordinate(std::string format, std::string file) {
            return geometrytransformation_mp_defineinternalcoordinate_(
                format.c_str(), file.c_str(), format.size(), file.size()
            );
        }
        inline void InternalCoordinate(const float * r, float * q, const int32_t & cartdim, const int32_t & intdim) {
            geometrytransformation_mp_sinternalcoordinate_(r, q, cartdim, intdim);
        }
        inline void InternalCoordinate(const double * r, double * q, const int32_t & cartdim, const int32_t & intdim) {
            geometrytransformation_mp_dinternalcoordinate_(r, q, cartdim, intdim);
        }
    }
#elif __GNUC__
    extern "C" {
        // General
        void __general_MOD_showtime();
        void __general_MOD_dscientificnotation(double & x, int32_t & i);
        // Chemistry
        void __chemistry_MOD_initializephasefixing(const int32_t & NStates);
        // GeometryTransformation
        int32_t __geometrytransformation_MOD_defineinternalcoordinate(const char * format, const char * file, int32_t len_format, int32_t len_file);
        void __geometrytransformation_MOD_sinternalcoordinate(const float * r, float * q, const int32_t & cartdim, const int32_t & intdim);
        void __geometrytransformation_MOD_dinternalcoordinate(const double * r, double * q, const int32_t & cartdim, const int32_t & intdim);
    }
    namespace General {
        inline void ShowTime() {__general_MOD_showtime();}
        inline void dScientificNotation(double & x, int32_t & i) {__general_MOD_dscientificnotation(x, i);}
    }
    namespace Chemistry {
        inline void InitializePhaseFixing(const int32_t & NStates) {__chemistry_MOD_initializephasefixing(NStates);}
    }
    namespace GeometryTransformation {
        inline int32_t DefineInternalCoordinate(std::string format, std::string file) {
            return __geometrytransformation_MOD_defineinternalcoordinate(
                format.c_str(), file.c_str(), format.size(), file.size()
            );
        }
        inline void InternalCoordinate(const float * r, float * q, const int32_t & cartdim, const int32_t & intdim) {
            __geometrytransformation_MOD_sinternalcoordinate(r, q, cartdim, intdim);
        }
        inline void InternalCoordinate(const double * r, double * q, const int32_t & cartdim, const int32_t & intdim) {
            __geometrytransformation_MOD_dinternalcoordinate(r, q, cartdim, intdim);
        }
    }
#endif
}

#endif