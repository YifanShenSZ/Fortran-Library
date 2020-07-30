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
    void general_mp_dscientificnotation_(double & x, int & i);
    // NonlinearOptimization
    void nonlinearoptimization_mp_conjugategradient_basic_(
        void (*f)(double &, const double *, const int &),
        void (*fd)(double *, const double *, const int &),
        double * x, const int & dim,
        const char * Method,
        const bool & Strong, const bool & Warning,
        const int & MaxIteration,
        const double & Precision, const double & MinStepLength,
        const double & WolfeConst1, const double & WolfeConst2, const double & Increment,
        int len_Method);
    void nonlinearoptimization_mp_conjugategradient_(
        void (*f)(double &, const double *, const int &),
        void (*fd)(double *, const double *, const int &),
        double * x, const int & dim,
        const char * Method,
        int (*f_fd)(double &, double *, const double *, const int &),
        const bool & Strong, const bool & Warning,
        const int & MaxIteration,
        const double & Precision, const double & MinStepLength,
        const double & WolfeConst1, const double & WolfeConst2, const double & Increment,
        int len_Method);
    void nonlinearoptimization_mp_trustregion_basic_(
        void (*fd)(double *, const double *, const int &, const int &),
        void (*Jacobian)(double *, const double *, const int &, const int &),
        double * x, const int & M, const int & N,
        const bool & Warning,
        const int & MaxIteration, const int & MaxStepIteration,
        const double & Precision, const double & MinStepLength);
    // Chemistry
    void chemistry_mp_initializephasefixing_(const int & NStates);
    // GeometryTransformation
    int geometrytransformation_mp_defineinternalcoordinate_(
        const char * format, const char * file, int len_format, int len_file);
    void geometrytransformation_mp_sinternalcoordinate_(
        const float * r, float * q, const int & cartdim, const int & intdim);
    void geometrytransformation_mp_dinternalcoordinate_(
        const double * r, double * q, const int & cartdim, const int & intdim);
    void geometrytransformation_mp_scartesian2internal_(
        const float * r, const float * cartgradT, float * q, float * intgradT, const int & cartdim, const int & intdim);
    void geometrytransformation_mp_dcartesian2internal_(
        const double * r, const double * cartgradT, double * q, double * intgradT, const int & cartdim, const int & intdim);
    void geometrytransformation_mp_swilsonbmatrixandinternalcoordinate_(
        const float * r, float * BT, float * q, const int & cartdim, const int & intdim);
    void geometrytransformation_mp_dwilsonbmatrixandinternalcoordinate_(
        const double * r, double * BT, double * q, const int & cartdim, const int & intdim);
    // FortranLibrary
    void fortranlibrary_mp_testfortranlibrary_();
}
namespace General {
    inline void ShowTime() {general_mp_showtime_();}
    inline void dScientificNotation(double & x, int & i) {general_mp_dscientificnotation_(x, i);}
}
namespace NonlinearOptimization {
    inline void ConjugateGradient(
        void (*f)(double &, const double *, const int &),
        void (*fd)(double *, const double *, const int &),
        double * x, const int & dim,
        std::string Method="DY",
        const bool & Strong=true, const bool & Warning=true,
        const int & MaxIteration=1000, const double & Precision=1e-15, const double & MinStepLength=1e-15,
        const double & WolfeConst1=1e-4, const double & WolfeConst2=0.45, const double & Increment=1.05
    ) {
        nonlinearoptimization_mp_conjugategradient_basic_(
            f, fd, x, dim,
            Method.c_str(),
            Strong, Warning, MaxIteration, Precision, MinStepLength, WolfeConst1, WolfeConst2, Increment,
            Method.size());
    }
    inline void ConjugateGradient(
        void (*f)(double &, const double *, const int &),
        void (*fd)(double *, const double *, const int &),
        int (*f_fd)(double &, double *, const double *, const int &),
        double * x, const int & dim,
        std::string Method="DY",
        const bool & Strong=true, const bool & Warning=true,
        const int & MaxIteration=1000, const double & Precision=1e-15, const double & MinStepLength=1e-15,
        const double & WolfeConst1=1e-4, const double & WolfeConst2=0.45, const double & Increment=1.05
    ) {
        nonlinearoptimization_mp_conjugategradient_(
            f, fd, x, dim,
            Method.c_str(),
            f_fd,
            Strong, Warning, MaxIteration, Precision, MinStepLength, WolfeConst1, WolfeConst2, Increment,
            Method.size());
    }
    inline void TrustRegion(
        void (*fd)(double *, const double *, const int &, const int &),
        void (*Jacobian)(double *, const double *, const int &, const int &),
        double * x, const int & M, const int & N,
        const bool & Warning=true,
        const int & MaxIteration=1000, const int & MaxStepIteration=100,
        const double & Precision=1e-15, const double & MinStepLength=1e-15
    ) {
        nonlinearoptimization_mp_trustregion_basic_(
            fd, Jacobian, x, M, N,
            Warning,
            MaxIteration, MaxStepIteration,
            Precision, MinStepLength);
    }
}
namespace Chemistry {
    inline void InitializePhaseFixing(const int & NStates) {chemistry_mp_initializephasefixing_(NStates);}
}
namespace GeometryTransformation {
    inline int DefineInternalCoordinate(std::string format, std::string file) {
        return geometrytransformation_mp_defineinternalcoordinate_(
            format.c_str(), file.c_str(), format.size(), file.size()
        );
    }
    inline void InternalCoordinate(const float * r, float * q, const int & cartdim, const int & intdim) {
        geometrytransformation_mp_sinternalcoordinate_(r, q, cartdim, intdim);
    }
    inline void InternalCoordinate(const double * r, double * q, const int & cartdim, const int & intdim) {
        geometrytransformation_mp_dinternalcoordinate_(r, q, cartdim, intdim);
    }
    inline void Cartesian2Internal(const float * r, const float * cartgradT, float * q, float * intgradT, const int & cartdim, const int & intdim) {
        geometrytransformation_mp_scartesian2internal_(r, cartgradT, q, intgradT, cartdim, intdim);
    }
    inline void Cartesian2Internal(const double * r, const double * cartgradT, double * q, double * intgradT, const int & cartdim, const int & intdim) {
        geometrytransformation_mp_dcartesian2internal_(r, cartgradT, q, intgradT, cartdim, intdim);
    }
    inline void WilsonBMatrixAndInternalCoordinate(const float * r, float * BT, float * q, const int & cartdim, const int & intdim) {
        geometrytransformation_mp_swilsonbmatrixandinternalcoordinate_(r, BT, q, cartdim, intdim);
    }
    inline void WilsonBMatrixAndInternalCoordinate(const double * r, double * BT, double * q, const int & cartdim, const int & intdim) {
        geometrytransformation_mp_dwilsonbmatrixandinternalcoordinate_(r, BT, q, cartdim, intdim);
    }
}
namespace FortranLibrary {
    inline void TestFortranLibrary() {fortranlibrary_mp_testfortranlibrary_();};
}
#elif __GNUC__
extern "C" {
    // General
    void __general_MOD_showtime();
    void __general_MOD_dscientificnotation(double & x, int & i);
    // Chemistry
    void __chemistry_MOD_initializephasefixing(const int & NStates);
    // GeometryTransformation
    int __geometrytransformation_MOD_defineinternalcoordinate(const char * format, const char * file, int len_format, int len_file);
    void __geometrytransformation_MOD_sinternalcoordinate(
        const float * r, float * q, const int & cartdim, const int & intdim
    );
    void __geometrytransformation_MOD_dinternalcoordinate(
        const double * r, double * q, const int & cartdim, const int & intdim
    );
    void __geometrytransformation_MOD_scartesian2internal(
        const float * r, const float * cartgradT, float * q, float * intgradT, const int & cartdim, const int & intdim
    );
    void __geometrytransformation_MOD_dcartesian2internal(
        const double * r, const double * cartgradT, double * q, double * intgradT, const int & cartdim, const int & intdim
    );
    void __geometrytransformation_MOD_swilsonbmatrixandinternalcoordinate(
        const float * r, float * BT, float * q, const int & cartdim, const int & intdim
    );
    void __geometrytransformation_MOD_dwilsonbmatrixandinternalcoordinate(
        const double * r, double * BT, double * q, const int & cartdim, const int & intdim
    );
    // FortranLibrary
    void __fortranlibrary_MOD_testfortranlibrary();
}
namespace General {
    inline void ShowTime() {__general_MOD_showtime();}
    inline void dScientificNotation(double & x, int & i) {__general_MOD_dscientificnotation(x, i);}
}
namespace Chemistry {
    inline void InitializePhaseFixing(const int & NStates) {__chemistry_MOD_initializephasefixing(NStates);}
}
namespace GeometryTransformation {
    inline int DefineInternalCoordinate(std::string format, std::string file) {
        return __geometrytransformation_MOD_defineinternalcoordinate(
            format.c_str(), file.c_str(), format.size(), file.size()
        );
    }
    inline void InternalCoordinate(const float * r, float * q, const int & cartdim, const int & intdim) {
        __geometrytransformation_MOD_sinternalcoordinate(r, q, cartdim, intdim);
    }
    inline void InternalCoordinate(const double * r, double * q, const int & cartdim, const int & intdim) {
        __geometrytransformation_MOD_dinternalcoordinate(r, q, cartdim, intdim);
    }
    inline void Cartesian2Internal(const float * r, const float * cartgradT, float * q, float * intgradT, const int & cartdim, const int & intdim) {
        __geometrytransformation_MOD_scartesian2internal(r, cartgradT, q, intgradT, cartdim, intdim);
    }
    inline void Cartesian2Internal(const double * r, const double * cartgradT, double * q, double * intgradT, const int & cartdim, const int & intdim) {
        __geometrytransformation_MOD_dcartesian2internal(r, cartgradT, q, intgradT, cartdim, intdim);
    }
    inline void WilsonBMatrixAndInternalCoordinate(const float * r, float * BT, float * q, const int & cartdim, const int & intdim) {
        __geometrytransformation_MOD_swilsonbmatrixandinternalcoordinate(r, BT, q, cartdim, intdim);
    }
    inline void WilsonBMatrixAndInternalCoordinate(const double * r, double * BT, double * q, const int & cartdim, const int & intdim) {
        __geometrytransformation_MOD_dwilsonbmatrixandinternalcoordinate(r, BT, q, cartdim, intdim);
    }
}
namespace FortranLibrary {
    inline void TestFortranLibrary() {
        __fortranlibrary_MOD_testfortranlibrary();
    };
}
#endif

} // namespace FL

#endif