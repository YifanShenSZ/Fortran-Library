// c++ interface to Fortran-Library
// Include this header so that C++ can link to libFL.so

#ifndef FortranLibrary_hpp
#define FortranLibrary_hpp

namespace FL {

// No need to fetch constant, directly provide them in c++
#define AMUInAU  1822.888486192;
#define AInAU    1.8897261339212517;
#define cm_1InAu 0.000004556335830019422;

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
    int32_t geometrytransformation_mp_defineinternalcoordinate_(
        const char * format, const char * file, int32_t len_format, int32_t len_file
    );
    void geometrytransformation_mp_sinternalcoordinate_(
        const float * r, float * q, const int32_t & cartdim, const int32_t & intdim
    );
    void geometrytransformation_mp_dinternalcoordinate_(
        const double * r, double * q, const int32_t & cartdim, const int32_t & intdim
    );
    void geometrytransformation_mp_scartesian2internal_(
        const float * r, const float * cartgradT, float * q, float * intgradT, const int32_t & cartdim, const int32_t & intdim
    );
    void geometrytransformation_mp_dcartesian2internal_(
        const double * r, const double * cartgradT, double * q, double * intgradT, const int32_t & cartdim, const int32_t & intdim
    );
    void geometrytransformation_mp_swilsonbmatrixandinternalcoordinate_(
        const float * r, float * BT, float * q, const int32_t & cartdim, const int32_t & intdim
    );
    void geometrytransformation_mp_dwilsonbmatrixandinternalcoordinate_(
        const double * r, double * BT, double * q, const int32_t & cartdim, const int32_t & intdim
    );
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
    inline void Cartesian2Internal(const float * r, const float * cartgradT, float * q, float * intgradT, const int32_t & cartdim, const int32_t & intdim) {
        geometrytransformation_mp_scartesian2internal_(r, cartgradT, q, intgradT, cartdim, intdim);
    }
    inline void Cartesian2Internal(const double * r, const double * cartgradT, double * q, double * intgradT, const int32_t & cartdim, const int32_t & intdim) {
        geometrytransformation_mp_dcartesian2internal_(r, cartgradT, q, intgradT, cartdim, intdim);
    }
    inline void WilsonBMatrixAndInternalCoordinate(const float * r, float * BT, float * q, const int32_t & cartdim, const int32_t & intdim) {
        geometrytransformation_mp_swilsonbmatrixandinternalcoordinate_(r, BT, q, cartdim, intdim);
    }
    inline void WilsonBMatrixAndInternalCoordinate(const double * r, double * BT, double * q, const int32_t & cartdim, const int32_t & intdim) {
        geometrytransformation_mp_dwilsonbmatrixandinternalcoordinate_(r, BT, q, cartdim, intdim);
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
    void __geometrytransformation_MOD_sinternalcoordinate(
        const float * r, float * q, const int32_t & cartdim, const int32_t & intdim
    );
    void __geometrytransformation_MOD_dinternalcoordinate(
        const double * r, double * q, const int32_t & cartdim, const int32_t & intdim
    );
    void __geometrytransformation_MOD_scartesian2internal(
        const float * r, const float * cartgradT, float * q, float * intgradT, const int32_t & cartdim, const int32_t & intdim
    );
    void __geometrytransformation_MOD_dcartesian2internal(
        const double * r, const double * cartgradT, double * q, double * intgradT, const int32_t & cartdim, const int32_t & intdim
    );
    void __geometrytransformation_MOD_swilsonbmatrixandinternalcoordinate(
        const float * r, float * BT, float * q, const int32_t & cartdim, const int32_t & intdim
    );
    void __geometrytransformation_MOD_dwilsonbmatrixandinternalcoordinate(
        const double * r, double * BT, double * q, const int32_t & cartdim, const int32_t & intdim
    );
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
    inline void Cartesian2Internal(const float * r, const float * cartgradT, float * q, float * intgradT, const int32_t & cartdim, const int32_t & intdim) {
        __geometrytransformation_MOD_scartesian2internal(r, cartgradT, q, intgradT, cartdim, intdim);
    }
    inline void Cartesian2Internal(const double * r, const double * cartgradT, double * q, double * intgradT, const int32_t & cartdim, const int32_t & intdim) {
        __geometrytransformation_MOD_dcartesian2internal(r, cartgradT, q, intgradT, cartdim, intdim);
    }
    inline void WilsonBMatrixAndInternalCoordinate(const float * r, float * BT, float * q, const int32_t & cartdim, const int32_t & intdim) {
        __geometrytransformation_MOD_swilsonbmatrixandinternalcoordinate(r, BT, q, cartdim, intdim);
    }
    inline void WilsonBMatrixAndInternalCoordinate(const double * r, double * BT, double * q, const int32_t & cartdim, const int32_t & intdim) {
        __geometrytransformation_MOD_dwilsonbmatrixandinternalcoordinate(r, BT, q, cartdim, intdim);
    }
}
#endif

} // namespace FL

#endif