// A test program on the c++ interface of Fortran-Library
#include <iostream>
#include <FortranLibrary.hpp>

using namespace std;

int main() {
    int i;
    double x;

    cout << ">>> Testing calling from c++... >>>\n";
        cout << "\nTime display\n";
            ShowTime();
        cout << "\nScientific notation\n";
            x=3564.1212587;
            dScientificNotation(x,i);
            cout << x - 3.5641212587 << i - 3 << "\n";
    cout << "\n<<< Calling from c++ test passed <<<\n";

    cout << "\n>>> Testing calling from fortran... >>>\n";
        TestFortranLibrary();
    cout << "<<< Calling from fortran test passed <<<\n";

}