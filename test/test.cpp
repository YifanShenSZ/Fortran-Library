// A test program on the c++ interface of Fortran-Library
#include <iostream>
#include <FortranLibrary.hpp>

extern "C" {
    #ifdef __INTEL_COMPILER
    void fortranlibrary_mp_testfortranlibrary_();
    #elif __GNUC__
    void __fortranlibrary_MOD_testfortranlibrary();
    #endif
}

int main() {
    int i;
    double x;

    std::cout << ">>> Testing calling from c++... >>>\n";
        std::cout << "\nTime display\n";
            FL::General::ShowTime();
        std::cout << "\nScientific notation\n";
            x=3564.1212587;
            FL::General::dScientificNotation(x,i);
            std::cout << x - 3.5641212587 << i - 3 << "\n";
    std::cout << "\n<<< Calling from c++ test passed <<<\n";

    std::cout << "\n>>> Testing calling from fortran... >>>\n";
    #ifdef __INTEL_COMPILER
    fortranlibrary_mp_testfortranlibrary_();
    #elif __GNUC__
    __fortranlibrary_MOD_testfortranlibrary();
    #endif
    std::cout << "<<< Calling from fortran test passed <<<\n";
}