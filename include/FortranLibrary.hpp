// c++ interface to Fortran-Library
// Include this header so that C++ can link to libFL.so
#ifndef FortranLibrary_HPP__
#define FortranLibrary_HPP__

// No need to fetch constant, directly provide them in c++
const double AMUInAU  = 1822.888486192;
const double AInAU    = 1.8897261339212517;
const double cm_1InAu = 0.000004556335830019422;

// Raw interface
// C recognizes fortran functions by their compiled names:
// fortran function Func contained in module Mod will be renamed as mod_mp_func_
// ( You may use command nm to view the compiled names: nm libFL.so )
extern "C" {
    // General
    void general_mp_showtime_();
    void general_mp_dscientificnotation_(double& x, int& i);
    // FortranLibrary
    void fortranlibrary_mp_testfortranlibrary_();
}

// Rename fortran functions back to their origin
// General
    void ShowTime() {general_mp_showtime_();}
    void dScientificNotation(double& x, int& i) {general_mp_dscientificnotation_(x,i);}
// FortranLibrary
    void TestFortranLibrary() {fortranlibrary_mp_testfortranlibrary_();}

#endif