// A test program on the c++ interface of Fortran-Library
#include <iostream>
#include "FortranLibraryCppInterface.h"
using namespace std;
int main() {
    int i;
    double x;

    cout << ">>> Testing calling from Fortran... >>>\n";
        FortranLibrary();
    cout << "<<< Calling from Fortran test passed <<<\n";

    cout << ">>> Testing calling from c++... >>>\n";
        cout << "\nTime display\n";
            ShowTime();
        cout << "\nScientific notation\n";
            x=3564.1212587;
            dScientificNotation(x,i);
            cout << i - 3 << "\n";
    cout << "<<< Calling from c++ test passed <<<\n";
}