// c++ interface to Fortran-Library
// Include this header so that C++ can link to libFL.so

#ifndef FortranLibrary_hpp
#define FortranLibrary_hpp

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <forward_list>
#include <tuple>
#include <iterator>
#include <cmath>

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
    void general_mp_showtime_();
    void general_mp_dscientificnotation_(double & x, int & i);
    // LinearAlgebra
    void linearalgebra_mp_my_dgemm_t_(double * A, double * B, double * C, const int & M, const int & K, const int & N);
    void linearalgebra_mp_my_dsyev_(const char & jobtype, double * A, double * eigval, const int & N);
    // NonlinearOptimization
    void nonlinearoptimization_mp_steepestdescent_(
        void (*f)(double &, const double *, const int &),
        void (*fd)(double *, const double *, const int &),
        double * x, const int & dim,
        int (*f_fd)(double &, double *, const double *, const int &),
        const bool & Strong, const bool & Warning,
        const int & MaxIteration,
        const double & Precision, const double & MinStepLength,
        const double & WolfeConst1, const double & WolfeConst2, const double & Increment);
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
    void geometrytransformation_mp_defineinternalcoordinate_(
        int & intdim, int & ID, const char * format, const char * file, int len_format, int len_file);
    void geometrytransformation_mp_sinternalcoordinate_(
        const float * r, float * q, const int & cartdim, const int & intdim, const int & ID);
    void geometrytransformation_mp_dinternalcoordinate_(
        const double * r, double * q, const int & cartdim, const int & intdim, const int & ID);
    void geometrytransformation_mp_scartesian2internal_(
        const float * r, const float * cartgradT, float * q, float * intgradT,
        const int & cartdim, const int & intdim, const int & NStates, const int & ID);
    void geometrytransformation_mp_dcartesian2internal_(
        const double * r, const double * cartgradT, double * q, double * intgradT,
        const int & cartdim, const int & intdim, const int & NStates, const int & ID);
    void geometrytransformation_mp_swilsonbmatrixandinternalcoordinate_(
        const float * r, float * BT, float * q, const int & cartdim, const int & intdim, const int & ID);
    void geometrytransformation_mp_dwilsonbmatrixandinternalcoordinate_(
        const double * r, double * BT, double * q, const int & cartdim, const int & intdim, const int & ID);
    void geometrytransformation_mp_cartesiancoordinate_(
        const double * q, double * r, const int & intdim, const int & cartdim, const double * r0, const int & ID);
    // FortranLibrary
    void fortranlibrary_mp_testfortranlibrary_();
}
namespace General {
    inline void ShowTime() {general_mp_showtime_();}
    inline void dScientificNotation(double & x, int & i) {general_mp_dscientificnotation_(x, i);}
}
namespace LA {
    inline void My_dgemm_T(double * A, double * B, double * C, const int & M, const int & K, const int & N) {
        linearalgebra_mp_my_dgemm_t_(A, B, C, M, K, N);
    }
    inline void My_dsyev(const char & jobtype, double * A, double * eigval, const int & N) {
        linearalgebra_mp_my_dsyev_(jobtype, A, eigval, N);
    }
}
namespace NO { // NonlinearOptimization
    inline void SteepestDescent(
        void (*f)(double &, const double *, const int &),
        void (*fd)(double *, const double *, const int &),
        int (*f_fd)(double &, double *, const double *, const int &),
        double * x, const int & dim,
        const bool & Strong=true, const bool & Warning=true,
        const int & MaxIteration=1000, const double & Precision=1e-15, const double & MinStepLength=1e-15,
        const double & WolfeConst1=1e-4, const double & WolfeConst2=0.45, const double & Increment=1.05
    ) {
        nonlinearoptimization_mp_steepestdescent_(
            f, fd, x, dim,
            f_fd,
            Strong, Warning, MaxIteration, Precision, MinStepLength, WolfeConst1, WolfeConst2, Increment);
    }
    inline void ConjugateGradient(
        void (*f)(double &, const double *, const int &),
        void (*fd)(double *, const double *, const int &),
        double * x, const int & dim,
        const std::string & Method="DY",
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
        const std::string & Method="DY",
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
            Warning, MaxIteration, MaxStepIteration, Precision, MinStepLength);
    }
}
namespace Chemistry {
    inline void InitializePhaseFixing(const int & NStates) {chemistry_mp_initializephasefixing_(NStates);}
}
namespace GT { // GeometryTransformation
    inline std::tuple<int, int> DefineInternalCoordinate(const std::string & format, const std::string & file) {
        int intdim, ID;
        geometrytransformation_mp_defineinternalcoordinate_(intdim, ID, format.c_str(), file.c_str(), format.size(), file.size());
        return std::make_tuple(intdim, ID);
    }
    inline void InternalCoordinate(const float * r, float * q, const int & cartdim, const int & intdim, const int & ID=1) {
        geometrytransformation_mp_sinternalcoordinate_(r, q, cartdim, intdim, ID);
    }
    inline void InternalCoordinate(const double * r, double * q, const int & cartdim, const int & intdim, const int & ID=1) {
        geometrytransformation_mp_dinternalcoordinate_(r, q, cartdim, intdim, ID);
    }
    inline void Cartesian2Internal(const float * r, const float * cartgradT, float * q, float * intgradT,
    const int & cartdim, const int & intdim, const int & NStates, const int & ID=1) {
        geometrytransformation_mp_scartesian2internal_(r, cartgradT, q, intgradT, cartdim, intdim, NStates, ID);
    }
    inline void Cartesian2Internal(const double * r, const double * cartgradT, double * q, double * intgradT,
    const int & cartdim, const int & intdim, const int & NStates, const int & ID=1) {
        geometrytransformation_mp_dcartesian2internal_(r, cartgradT, q, intgradT, cartdim, intdim, NStates, ID);
    }
    inline void WilsonBMatrixAndInternalCoordinate(const float * r, float * BT, float * q, const int & cartdim, const int & intdim, const int & ID=1) {
        geometrytransformation_mp_swilsonbmatrixandinternalcoordinate_(r, BT, q, cartdim, intdim, ID);
    }
    inline void WilsonBMatrixAndInternalCoordinate(const double * r, double * BT, double * q, const int & cartdim, const int & intdim, const int & ID=1) {
        geometrytransformation_mp_dwilsonbmatrixandinternalcoordinate_(r, BT, q, cartdim, intdim, ID);
    }
    inline void CartesianCoordinate(const double * q, double * r, const int & intdim, const int & cartdim, const double * r0, const int & ID=1) {
        geometrytransformation_mp_cartesiancoordinate_(q, r, intdim, cartdim, r0, ID);
    }
}
namespace FL { // FortranLibrary
    inline void TestFortranLibrary() {fortranlibrary_mp_testfortranlibrary_();};
}
#elif __GNUC__
extern "C" {
    // General
    void __general_MOD_showtime();
    void __general_MOD_dscientificnotation(double & x, int & i);
    // LinearAlgebra
    void __linearalgebra_MOD_my_dgemm_t(double * A, double * B, double * C, const int & M, const int & K, const int & N);
    void __linearalgebra_MOD_my_dsyev(const char & jobtype, double * A, double * eigval, const int & N);
    // NonlinearOptimization
    void __nonlinearoptimization_MOD_steepestdescent(
        void (*f)(double &, const double *, const int &),
        void (*fd)(double *, const double *, const int &),
        double * x, const int & dim,
        int (*f_fd)(double &, double *, const double *, const int &),
        const bool & Strong, const bool & Warning,
        const int & MaxIteration,
        const double & Precision, const double & MinStepLength,
        const double & WolfeConst1, const double & WolfeConst2, const double & Increment);
    void __nonlinearoptimization_MOD_conjugategradient_basic(
        void (*f)(double &, const double *, const int &),
        void (*fd)(double *, const double *, const int &),
        double * x, const int & dim,
        const char * Method,
        const bool & Strong, const bool & Warning,
        const int & MaxIteration,
        const double & Precision, const double & MinStepLength,
        const double & WolfeConst1, const double & WolfeConst2, const double & Increment,
        int len_Method);
    void __nonlinearoptimization_MOD_conjugategradient(
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
    void __nonlinearoptimization_MOD_trustregion_basic(
        void (*fd)(double *, const double *, const int &, const int &),
        void (*Jacobian)(double *, const double *, const int &, const int &),
        double * x, const int & M, const int & N,
        const bool & Warning,
        const int & MaxIteration, const int & MaxStepIteration,
        const double & Precision, const double & MinStepLength);
    // Chemistry
    void __chemistry_MOD_initializephasefixing(const int & NStates);
    // GeometryTransformation
    int __geometrytransformation_MOD_defineinternalcoordinate(
        int & intdim, int & ID, const char * format, const char * file, int len_format, int len_file);
    void __geometrytransformation_MOD_sinternalcoordinate(
        const float * r, float * q, const int & cartdim, const int & intdim, const int & ID);
    void __geometrytransformation_MOD_dinternalcoordinate(
        const double * r, double * q, const int & cartdim, const int & intdim, const int & ID);
    void __geometrytransformation_MOD_scartesian2internal(
        const float * r, const float * cartgradT, float * q, float * intgradT,
        const int & cartdim, const int & intdim, const int & NStates, const int & ID);
    void __geometrytransformation_MOD_dcartesian2internal(
        const double * r, const double * cartgradT, double * q, double * intgradT,
        const int & cartdim, const int & intdim, const int & NStates, const int & ID);
    void __geometrytransformation_MOD_swilsonbmatrixandinternalcoordinate(
        const float * r, float * BT, float * q, const int & cartdim, const int & intdim, const int & ID);
    void __geometrytransformation_MOD_dwilsonbmatrixandinternalcoordinate(
        const double * r, double * BT, double * q, const int & cartdim, const int & intdim, const int & ID);
    void __geometrytransformation_MOD_cartesiancoordinate(
        const double * q, double * r, const int & intdim, const int & cartdim, const double * r0, const int & ID);
    // FortranLibrary
    void __fortranlibrary_MOD_testfortranlibrary();
}
namespace General {
    inline void ShowTime() {__general_MOD_showtime();}
    inline void dScientificNotation(double & x, int & i) {__general_MOD_dscientificnotation(x, i);}
}
namespace LA {
    inline void My_dgemm_T(double * A, double * B, double * C, const int & M, const int & K, const int & N) {
        __linearalgebra_MOD_my_dgemm_t(A, B, C, M, K, N);
    }
    inline void My_dsyev(const char & jobtype, double * A, double * eigval, const int & N) {
        __linearalgebra_MOD_my_dsyev(jobtype, A, eigval, N);
    }
}
namespace NO { // NonlinearOptimization
    inline void SteepestDescent(
        void (*f)(double &, const double *, const int &),
        void (*fd)(double *, const double *, const int &),
        int (*f_fd)(double &, double *, const double *, const int &),
        double * x, const int & dim,
        const bool & Strong=true, const bool & Warning=true,
        const int & MaxIteration=1000, const double & Precision=1e-15, const double & MinStepLength=1e-15,
        const double & WolfeConst1=1e-4, const double & WolfeConst2=0.45, const double & Increment=1.05
    ) {
        __nonlinearoptimization_MOD_steepestdescent(
            f, fd, x, dim,
            f_fd,
            Strong, Warning, MaxIteration, Precision, MinStepLength, WolfeConst1, WolfeConst2, Increment);
    }
    inline void ConjugateGradient(
        void (*f)(double &, const double *, const int &),
        void (*fd)(double *, const double *, const int &),
        double * x, const int & dim,
        const std::string & Method="DY",
        const bool & Strong=true, const bool & Warning=true,
        const int & MaxIteration=1000, const double & Precision=1e-15, const double & MinStepLength=1e-15,
        const double & WolfeConst1=1e-4, const double & WolfeConst2=0.45, const double & Increment=1.05
    ) {
        __nonlinearoptimization_MOD_conjugategradient_basic(
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
        const std::string & Method="DY",
        const bool & Strong=true, const bool & Warning=true,
        const int & MaxIteration=1000, const double & Precision=1e-15, const double & MinStepLength=1e-15,
        const double & WolfeConst1=1e-4, const double & WolfeConst2=0.45, const double & Increment=1.05
    ) {
        __nonlinearoptimization_MOD_conjugategradient(
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
        __nonlinearoptimization_MOD_trustregion_basic(
            fd, Jacobian, x, M, N,
            Warning,
            MaxIteration, MaxStepIteration,
            Precision, MinStepLength);
    }
}
namespace Chemistry {
    inline void InitializePhaseFixing(const int & NStates) {__chemistry_MOD_initializephasefixing(NStates);}
}
namespace GT { // GeometryTransformation
    inline std::tuple<int, int> DefineInternalCoordinate(const std::string & format, const std::string & file) {
        int intdim, ID;
        __geometrytransformation_MOD_defineinternalcoordinate(intdim, ID, format.c_str(), file.c_str(), format.size(), file.size());
        return std::make_tuple(intdim, ID);
    }
    inline void InternalCoordinate(const float * r, float * q, const int & cartdim, const int & intdim, const int & ID=1) {
        __geometrytransformation_MOD_sinternalcoordinate(r, q, cartdim, intdim, ID);
    }
    inline void InternalCoordinate(const double * r, double * q, const int & cartdim, const int & intdim, const int & ID=1) {
        __geometrytransformation_MOD_dinternalcoordinate(r, q, cartdim, intdim, ID);
    }
    inline void Cartesian2Internal(const float * r, const float * cartgradT, float * q, float * intgradT,
    const int & cartdim, const int & intdim, const int & NStates, const int & ID=1) {
        __geometrytransformation_MOD_scartesian2internal(r, cartgradT, q, intgradT, cartdim, intdim, NStates, ID);
    }
    inline void Cartesian2Internal(const double * r, const double * cartgradT, double * q, double * intgradT,
    const int & cartdim, const int & intdim, const int & NStates, const int & ID=1) {
        __geometrytransformation_MOD_dcartesian2internal(r, cartgradT, q, intgradT, cartdim, intdim, NStates, ID);
    }
    inline void WilsonBMatrixAndInternalCoordinate(const float * r, float * BT, float * q, const int & cartdim, const int & intdim, const int & ID=1) {
        __geometrytransformation_MOD_swilsonbmatrixandinternalcoordinate(r, BT, q, cartdim, intdim, ID);
    }
    inline void WilsonBMatrixAndInternalCoordinate(const double * r, double * BT, double * q, const int & cartdim, const int & intdim, const int & ID=1) {
        __geometrytransformation_MOD_dwilsonbmatrixandinternalcoordinate(r, BT, q, cartdim, intdim, ID);
    }
    inline void CartesianCoordinate(const double * q, double * r, const int & intdim, const int & cartdim, const double * r0, const int & ID=1) {
        __geometrytransformation_MOD_cartesiancoordinate(q, r, intdim, cartdim, r0, ID);
    }
}
namespace FL { // FortranLibrary
    inline void TestFortranLibrary() {
        __fortranlibrary_MOD_testfortranlibrary();
    };
}
#endif

namespace GT { // GeometryTransformation
    struct InvolvedMotion {
        std::string type;
        std::vector<int> atom;
        double coeff;
    
        inline InvolvedMotion(const std::string & type, const std::vector<int> & atom, const double & coeff) {
            this->type  = type;
            this->atom  = atom;
            this->coeff = coeff;
        }
        ~InvolvedMotion() {}
    };
    struct IntCoordDef {
        std::vector<InvolvedMotion> motion;
    
        IntCoordDef() {}
        ~IntCoordDef() {}
    };
    // We provide a function to load internal coordinate format in c++
    // instead of fetching GeometryTransformation_IntCoordDef
    inline std::tuple<int, std::vector<IntCoordDef>> FetchInternalCoordinateDefinition(const std::string & format, const std::string & file) {
        int intdim = 0;
        std::vector<IntCoordDef> def;
        if (format == "Columbus7") {
            // First line is always "TEXAS"
            // New internal coordinate line starts with 'K'
            std::ifstream ifs; ifs.open(file);
            if (! ifs.good()) {ifs.close(); ifs.open("intcfl");}
            std::string line; std::getline(ifs, line);
            while (true) {
                std::string line; std::getline(ifs, line);
                if (! ifs.good()) break;
                if (line[0] == 'K') {
                    intdim++;
                    def.push_back(IntCoordDef());
                }
                std::string type;
                std::vector<int> atom;
                if (line.substr(20,4) == "STRE") {
                    type = "stretching";
                    atom = {std::stoi(line.substr(28, 5)), std::stoi(line.substr(34, 9))};
                } else if (line.substr(20,4) == "BEND") {
                    type = "bending";
                    atom = {std::stoi(line.substr(28, 6)), std::stoi(line.substr(45, 9)), std::stoi(line.substr(35, 9))};
                } else if (line.substr(20,4) == "TORS") {
                    type = "torsion";
                    atom = {std::stoi(line.substr(28, 6)), std::stoi(line.substr(35, 9)), std::stoi(line.substr(45, 9)), std::stoi(line.substr(55, 9))};
                } else if (line.substr(20,3) == "OUT") {
                    type = "OutOfPlane";
                    atom = {std::stoi(line.substr(28, 6)), std::stoi(line.substr(45, 9)), std::stoi(line.substr(55, 9)), std::stoi(line.substr(35, 9))};
                }
                double coeff = 1.0;
                if (line.substr(10, 10) != "          ") coeff = std::stod(line.substr(10, 10));
                def[intdim-1].motion.push_back(InvolvedMotion(type, atom, coeff));
            }
            ifs.close();
        }
        else {
            // First 6 spaces of a line are reserved to indicate the start of new internal coordinate
            // Example:
            //  coor |   coeff   |    type     |      atom
            // --------------------------------------------------
            //      1    1.000000    stretching     1     2          # Comment
            //           1.000000    stretching     1     3
            //      2    1.000000    stretching     1     2
            //          -1.000000    stretching     1     3
            //      3    1.000000       bending     2     1     3
            std::ifstream ifs; ifs.open(file);
            if (! ifs.good()) {ifs.close(); ifs.open("IntCoordDef");}
            while (true) {
                std::string line; std::getline(ifs, line);
                if (! ifs.good()) break;
                std::istringstream iss(line);
                std::forward_list<std::string> strs(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());
                if (line.substr(0, 6) != "      ") {
                    intdim++;
                    def.push_back(IntCoordDef());
                    strs.pop_front();
                }
                double coeff = std::stod(strs.front()); strs.pop_front();
                std::string type = strs.front(); strs.pop_front();
                std::vector<int> atom;
                if (type == "stretching") {
                    atom.resize(2);
                    atom[0] = std::stoi(strs.front()); strs.pop_front();
                    atom[1] = std::stoi(strs.front());
                }
                else if (type == "bending") {
                    atom.resize(3);
                    atom[0] = std::stoi(strs.front()); strs.pop_front();
                    atom[1] = std::stoi(strs.front()); strs.pop_front();
                    atom[2] = std::stoi(strs.front());
                }
                else if (type == "torsion") {
                    atom.resize(4);
                    atom[0] = std::stoi(strs.front()); strs.pop_front();
                    atom[1] = std::stoi(strs.front()); strs.pop_front();
                    atom[2] = std::stoi(strs.front()); strs.pop_front();
                    atom[3] = std::stoi(strs.front());
                }
                else if (type == "OutOfPlane") {
                    atom.resize(4);
                    atom[0] = std::stoi(strs.front()); strs.pop_front();
                    atom[1] = std::stoi(strs.front()); strs.pop_front();
                    atom[2] = std::stoi(strs.front()); strs.pop_front();
                    atom[3] = std::stoi(strs.front());
                }
                def[intdim-1].motion.push_back(InvolvedMotion(type, atom, coeff));
            }
            ifs.close();
        }
        // Normalized linear combination coefficient
        for (size_t i = 0; i < intdim; i++) {
            double norm2 = 0.0;
            for (size_t j = 0; j < def[i].motion.size(); j++) {
                norm2 += def[i].motion[j].coeff * def[i].motion[j].coeff;
            }
            norm2 = sqrt(norm2);
            for (size_t j = 0; j < def[i].motion.size(); j++) {
                def[i].motion[j].coeff /= norm2;
            }
        }
        return std::make_tuple(intdim, def);
    }
}

} // namespace FL

#endif