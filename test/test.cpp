// A test program on the c++ interface of Fortran-Library
#include <iostream>
#include <cmath>
#include <FortranLibrary.hpp>

// Routines for testing nonlinear-optimization solvers
void f(double & fx, const double * x, const int & dim);
void fd(double * fdx, const double * x, const int & dim);
int f_fd(double & fx, double * fdx, const double * x, const int & dim);
void fd_tr(double * fdx, const double * x, const int & M, const int & N);
void fdd_tr(double * fddx, const double * x, const int & M, const int & N);

double norm(const double * x, const int & dim) {
    double norm = 0;
    for (int i = 0; i < dim; i++) norm += x[i]*x[i];
    return sqrt(norm);
}

int main() {
    int i, dim;
    double x;
    double * vec;

    std::cout << ">>> Testing calling from c++... >>>\n";
        std::cout << "\nTime display" << std::endl;
            FL::General::ShowTime();
        std::cout << "\nScientific notation\n";
            x = 3564.1212587;
            FL::General::dScientificNotation(x,i);
            std::cout << x - 3.5641212587 << ' ' << i - 3 << '\n';
        std::cout << "\n!!!!!!!!!! Testing all nonlinear-optimization solvers... !!!!!!!!!!\n";
            dim = 10;
            vec = new double[dim];
            std::cout << "DY_S\n";
            for (i=0; i<dim; i++) vec[i] = (double)rand() / (double)RAND_MAX;
            FL::NO::ConjugateGradient(f, fd, vec, dim);
            std::cout << norm(vec, dim) << "\n\n";
            std::cout << "DY_S__fdwithf\n";
            for (i=0; i<dim; i++) vec[i] = (double)rand() / (double)RAND_MAX;
            FL::NO::ConjugateGradient(f, fd, f_fd, vec, dim);
            std::cout << norm(vec, dim) << "\n\n";
            std::cout << "dtrnlsp\n";
            for (i=0; i<dim; i++) vec[i] = (double)rand() / (double)RAND_MAX;
            FL::NO::TrustRegion(fd_tr, fdd_tr, vec, dim, dim);
            std::cout << norm(vec, dim) << "\n\n";
    std::cout << "\n<<< Calling from c++ test passed <<<\n";

    std::cout << "\n>>> Testing calling from fortran... >>>" << std::endl;
        FL::FL::TestFortranLibrary();
    std::cout << "<<< Calling from fortran test passed <<<\n";
}

// Routines for testing nonlinear-optimization solvers
void f(double & fx, const double * x, const int & dim) {
    fx = 0.0;
    for (int i = 0; i < dim; i++) fx += pow(x[i], 4);
}
void fd(double * fdx, const double * x, const int & dim) {
    for (int i = 0; i < dim; i++) fdx[i] = 4.0 * pow(x[i], 3);
}
int f_fd(double & fx, double * fdx, const double * x, const int & dim) {
    fx = 0.0;
    for (int i = 0; i < dim; i++) {
        fx += pow(x[i], 4);
        fdx[i] = 4.0 * pow(x[i], 3);
    }
    return 0;
}
void fd_tr(double * fdx, const double * x, const int & M, const int & N) {
    for (int i = 0; i < M; i++) fdx[i] = 4.0 * pow(x[i], 3);
}
void fdd_tr(double * fddx, const double * x, const int & M, const int & N) {
    for (int i = 0; i < M; i++)
    for (int j = 0; j < N; j++) {
        if (i == j) {
            fddx[i * N + j] = 12.0 * x[i]*x[i];
        } else {
            fddx[i * N + j] = 0.0;
        }
    }
}