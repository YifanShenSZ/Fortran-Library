/* c++ interface to Fortran-Library
   Include this header so that your c++ can link to libFL */

/* Raw interface
   C recognizes fortran functions by their object names:
   fortran function Func contained in module Mod will be renamed as mod_mp_func_
   You may use command nm to view the object names, e.g.:
   nm libFL.a, then search for symbol 'T' */
extern "C" {
    // General
    void general_mp_showtime_();
    void general_mp_dscientificnotation_(const double& x, const int& i);
    // FortranLibrary
    void fortranlibrary_mp_testfortranlibrary_();
}

// Rename fortran functions back to their origin
// General
    void ShowTime() {general_mp_showtime_();}
    void dScientificNotation(const double& x, const int& i) {general_mp_dscientificnotation_(x,i);}
// FortranLibrary
    void FortranLibrary() {fortranlibrary_mp_testfortranlibrary_();}