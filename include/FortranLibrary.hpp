// c++ interface to Fortran-Library
// Include this header so that C++ can link to libFL.so
#ifndef FortranLibrary_hpp
#define FortranLibrary_hpp

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
        void general_mp_dscientificnotation_(double& x, int& i);
        // FortranLibrary
        void fortranlibrary_mp_testfortranlibrary_();
    }
    // General
        void ShowTime() {general_mp_showtime_();}
        void dScientificNotation(double& x, int& i) {general_mp_dscientificnotation_(x, i);}
    // FortranLibrary
        void TestFortranLibrary() {fortranlibrary_mp_testfortranlibrary_();}
#elif __GNUC__
    extern "C" {
        // General
        void __general_MOD_showtime();
        void __general_MOD_dscientificnotation(double& x, int& i);
        // FortranLibrary
        void __fortranlibrary_MOD_testfortranlibrary();
    }
    // General
        void ShowTime() {__general_MOD_showtime();}
        void dScientificNotation(double& x, int& i) {__general_MOD_dscientificnotation(x, i);}
    // FortranLibrary
        void TestFortranLibrary() {__fortranlibrary_MOD_testfortranlibrary();}
#endif

#endif