// A test program on the c++ interface of Fortran-Library
#include <iostream>
#include <FortranLibrary.hpp>

int main() {
    int i;
    double x;

    std::cout << ">>> Testing calling from c++... >>>\n";
        std::cout << "\nTime display\n";
            ShowTime();
        std::cout << "\nScientific notation\n";
            x=3564.1212587;
            dScientificNotation(x,i);
            std::cout << x - 3.5641212587 << i - 3 << "\n";
    std::cout << "\n<<< Calling from c++ test passed <<<\n";

    std::cout << "\n>>> Testing calling from fortran... >>>\n";
        TestFortranLibrary();
    std::cout << "<<< Calling from fortran test passed <<<\n";

}