// A test program on the c++ interface of Fortran-Library
#include <iostream>
#include <cmath>
#include <FortranLibrary.hpp>

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

int fdd(double * fddx, const double * x, const int & dim) {
    for (int i = 0; i < dim; i++) 
    for (int j = 0; j < dim; j++) {
        if (i==j) fddx[i * dim + j] = 12.0 * x[i] * x[i];
        else      fddx[i * dim + j] = 0.0;
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
            fddx[i * N + j] = 12.0 * x[i] * x[i];
        } else {
            fddx[i * N + j] = 0.0;
        }
    }
}

void constraint(double * cx, const double * x, const int & M, const int & N) {
    cx[0] = -1.0;
    for (int i = 0; i < N; i++) cx[0] += x[i] * x[i];
}

void constraintd(double * cdx, const double * x, const int & M, const int & N) {
    for (int i = 0; i < N; i++) cdx[i] += 2.0 * x[i];
}

int constraintdd(double * cddx, const double * x, const int & M, const int & N) {
    for (int i = 0; i < N; i++) 
    for (int j = 0; j < N; j++) {
        if (i==j) cddx[i * N + j] = 2.0;
        else      cddx[i * N + j] = 0.0;
    }
    return 0;
}

double norm(const double * x, const int & dim) {
    double norm = 0;
    for (int i = 0; i < dim; i++) norm += x[i]*x[i];
    return sqrt(norm);
}

int main() {
    std::cout << "This is a test program on Fortran-Library c++ interface\n"
              << "Correct routines should print close to 0\n";

    std::cout << "\nScientific notation\n"; {
        double x = 3564.1212587;
        int i;
        FL::General::dScientificNotation(x,i);
        std::cout << x - 3.5641212587 << ' ' << i - 3 << '\n';
    }

    std::cout << "\n!!!!!!!!!! Testing all nonlinear-optimization solvers... !!!!!!!!!!\n"; {
        int dim = 10;
        double * x = new double[dim];

        std::cout << "Steepest descent\n";
        for (int i = 0; i < dim; i++) x[i] = (double)rand() / (double)RAND_MAX;
        FL::NO::SteepestDescent(f, fd, f_fd, x, dim);
        std::cout << norm(x, dim) << "\n\n";

        std::cout << "Dai-Yuan conjugate gradient: basic version, target and gradient are computed separately\n";
        for (int i = 0; i < dim; i++) x[i] = (double)rand() / (double)RAND_MAX;
        FL::NO::ConjugateGradient(f, fd, x, dim);
        std::cout << norm(x, dim) << "\n\n";

        std::cout << "Dai-Yuan conjugate gradient\n";
        for (int i = 0; i < dim; i++) x[i] = (double)rand() / (double)RAND_MAX;
        FL::NO::ConjugateGradient(f, fd, f_fd, x, dim);
        std::cout << norm(x, dim) << "\n\n";
 
        std::cout << "BFGS\n";
        for (int i = 0; i < dim; i++) x[i] = (double)rand() / (double)RAND_MAX;
        FL::NO::BFGS(f, fd, f_fd, fdd, x, dim);
        std::cout << norm(x, dim) << "\n\n";

        std::cout << "Newton\n";
        for (int i = 0; i < dim; i++) x[i] = (double)rand() / (double)RAND_MAX;
        FL::NO::NewtonRaphson(f, fd, f_fd, fdd, x, dim);
        std::cout << norm(x, dim) << "\n\n";

        std::cout << "dtrnlsp\n";
        for (int i = 0; i < dim; i++) x[i] = (double)rand() / (double)RAND_MAX;
        FL::NO::TrustRegion(fd_tr, fdd_tr, x, dim, dim);
        std::cout << norm(x, dim) << "\n\n";

        std::cout << "Augmented Lagrangian\n";
        for (int i = 0; i < dim; i++) x[i] = (double)rand() / (double)RAND_MAX;
        int M = 1;
        FL::NO::AugmentedLagrangian(f, fd, f_fd, fdd, constraint, constraintd, constraintdd, x, dim, M,
                                    "BFGS", {}, 1.0, 20, 10, "DY", true, false);
        std::cout << norm(x, dim) - 1.0 << "\n\n";
        delete [] x;
    }
}