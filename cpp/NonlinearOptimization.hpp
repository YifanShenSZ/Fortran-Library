#ifndef NonlinearOptimization_hpp
#define NonlinearOptimization_hpp

#include <cstdint>
#include <vector>

namespace FL { namespace NO {

#ifdef __INTEL_COMPILER
    extern "C" {
        void nonlinearoptimization_mp_steepestdescent_(
            // Required argument
            void (*f)(double &, const double *, const int &),
            void (*fd)(double *, const double *, const int &),
            double * x, const int & dim,
            int (*f_fd)(double &, double *, const double *, const int &),
            // Common optional argument
            const int32_t & Strong, const int32_t & Warning,
            const int & MaxIteration,
            const double & Precision, const double & MinStepLength,
            const double & WolfeConst1, const double & WolfeConst2, const double & Increment
        );

        void nonlinearoptimization_mp_conjugategradient_basic_(
            // Required argument
            void (*f)(double &, const double *, const int &),
            void (*fd)(double *, const double *, const int &),
            double * x, const int & dim,
            // Conjugate gradient method
            const char * Method,
            // Common optional argument
            const int32_t & Strong, const int32_t & Warning,
            const int & MaxIteration,
            const double & Precision, const double & MinStepLength,
            const double & WolfeConst1, const double & WolfeConst2, const double & Increment,
            // Fortran requires string length
            int len_Method
        );

        void nonlinearoptimization_mp_conjugategradient_(
            // Required argument
            void (*f)(double &, const double *, const int &),
            void (*fd)(double *, const double *, const int &),
            double * x, const int & dim,
            // Conjugate gradient method
            const char * Method,
            // Required argument
            int (*f_fd)(double &, double *, const double *, const int &),
            // Common optional argument
            const int32_t & Strong, const int32_t & Warning,
            const int & MaxIteration,
            const double & Precision, const double & MinStepLength,
            const double & WolfeConst1, const double & WolfeConst2, const double & Increment,
            // Fortran requires string length
            int len_Method
        );

        void nonlinearoptimization_mp_bfgs_(
            // Required argument
            void (*f)(double &, const double *, const int &),
            void (*fd)(double *, const double *, const int &),
            double * x, const int & dim,
            int (*fdd)(double *, const double *, const int &),
            // Every how many steps compute exact hessian
            const int & ExactStep,
            // Required argument
            int (*f_fd)(double &, double *, const double *, const int &),
            // Common optional argument
            const int32_t & Strong, const int32_t & Warning,
            const int & MaxIteration,
            const double & Precision, const double & MinStepLength,
            const double & WolfeConst1, const double & WolfeConst2, const double & Increment
        );

        void nonlinearoptimization_mp_newtonraphson_(
            // Required argument
            void (*f)(double &, const double *, const int &),
            void (*fd)(double *, const double *, const int &),
            double * x, const int & dim,
            int (*fdd)(double *, const double *, const int &),
            int (*f_fd)(double &, double *, const double *, const int &),
            // Common optional argument
            const int32_t & Strong, const int32_t & Warning,
            const int & MaxIteration,
            const double & Precision, const double & MinStepLength,
            const double & WolfeConst1, const double & WolfeConst2, const double & Increment
        );

        void nonlinearoptimization_mp_trustregion_basic_(
            // Required argument
            void (*residue)(double *, const double *, const int &, const int &),
            void (*Jacobian)(double *, const double *, const int &, const int &),
            double * x, const int & M, const int & N,
            // Optional argument
            const int32_t & Warning,
            const int & MaxIteration, const int & MaxStepIteration,
            const double & Precision, const double & MinStepLength
        );

        void nonlinearoptimization_mp_augmentedlagrangian_(
            // Required argument
            void (*f)(double &, const double *, const int &),
            void (*fd)(double *, const double *, const int &),
            void (*c)(double *, const double *, const int &, const int &),
            void (*cd)(double *, const double *, const int &, const int &),
            double * x, const int & N, const int & M,
            // Augmented Lagrangian parameters
            const char * UnconstrainedSolver, const double * lambda0, const double & miu0,
            // Required argument
            int (*fdd)(double *, const double *, const int &),
            int (*cdd)(double *, const double *, const int &, const int &),
            // Unconstrained solver parameters
            const int & ExactStep, const int & Memory, const char * Method,
            // Required argument
            int (*f_fd)(double &, double *, const double *, const int &),
            // Common optional argument
            const int32_t & Strong, const int32_t & Warning,
            const int & MaxIteration,
            const double & Precision, const double & MinStepLength,
            const double & WolfeConst1, const double & WolfeConst2, const double & Increment,
            // Fortran requires string length
            int len_UnconstrainedSolver, int len_Method
        );
    }

    inline void SteepestDescent(
        void (*f)(double &, const double *, const int &),
        void (*fd)(double *, const double *, const int &),
        int (*f_fd)(double &, double *, const double *, const int &),
        double * x, const int & dim,
        const bool & Strong=true, const bool & Warning=true,
        const int & MaxIteration=1000, const double & Precision=1e-15, const double & MinStepLength=1e-15,
        const double & WolfeConst1=1e-4, const double & WolfeConst2=0.9, const double & Increment=1.05
    ) {
        int32_t s, w;
        if (Strong ) s = -1; else s = 0;
        if (Warning) w = -1; else w = 0;
        nonlinearoptimization_mp_steepestdescent_(
            f, fd, x, dim,
            f_fd,
            s, w, MaxIteration, Precision, MinStepLength, WolfeConst1, WolfeConst2, Increment
        );
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
        int32_t s, w;
        if (Strong ) s = -1; else s = 0;
        if (Warning) w = -1; else w = 0;
        nonlinearoptimization_mp_conjugategradient_basic_(
            f, fd, x, dim,
            Method.c_str(),
            s, w, MaxIteration, Precision, MinStepLength, WolfeConst1, WolfeConst2, Increment,
            Method.size()
        );
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
        int32_t s, w;
        if (Strong ) s = -1; else s = 0;
        if (Warning) w = -1; else w = 0;
        nonlinearoptimization_mp_conjugategradient_(
            f, fd, x, dim,
            Method.c_str(),
            f_fd,
            s, w, MaxIteration, Precision, MinStepLength, WolfeConst1, WolfeConst2, Increment,
            Method.size()
        );
    }

    inline void BFGS(
        void (*f)(double &, const double *, const int &),
        void (*fd)(double *, const double *, const int &),
        int (*f_fd)(double &, double *, const double *, const int &),
        int (*fdd)(double *, const double *, const int &),
        double * x, const int & dim,
        const int & ExactStep=20,
        const bool & Strong=true, const bool & Warning=true,
        const int & MaxIteration=1000, const double & Precision=1e-15, const double & MinStepLength=1e-15,
        const double & WolfeConst1=1e-4, const double & WolfeConst2=0.9, const double & Increment=1.05
    ) {
        int32_t s, w;
        if (Strong ) s = -1; else s = 0;
        if (Warning) w = -1; else w = 0;
        nonlinearoptimization_mp_bfgs_(
            f, fd, x, dim,
            fdd, ExactStep,
            f_fd,
            s, w, MaxIteration, Precision, MinStepLength, WolfeConst1, WolfeConst2, Increment
        );
    }

    inline void NewtonRaphson(
        void (*f)(double &, const double *, const int &),
        void (*fd)(double *, const double *, const int &),
        int (*f_fd)(double &, double *, const double *, const int &),
        int (*fdd)(double *, const double *, const int &),
        double * x, const int & dim,
        const bool & Strong=true, const bool & Warning=true,
        const int & MaxIteration=1000, const double & Precision=1e-15, const double & MinStepLength=1e-15,
        const double & WolfeConst1=1e-4, const double & WolfeConst2=0.9, const double & Increment=1.05
    ) {
        int32_t s, w;
        if (Strong ) s = -1; else s = 0;
        if (Warning) w = -1; else w = 0;
        nonlinearoptimization_mp_newtonraphson_(
            f, fd, x, dim,
            fdd, f_fd,
            s, w, MaxIteration, Precision, MinStepLength, WolfeConst1, WolfeConst2, Increment
        );
    }

    inline void TrustRegion(
        void (*residue)(double *, const double *, const int &, const int &),
        void (*Jacobian)(double *, const double *, const int &, const int &),
        double * x, const int & M, const int & N,
        const bool & Warning=true,
        const int & MaxIteration=1000, const int & MaxStepIteration=100,
        const double & Precision=1e-15, const double & MinStepLength=1e-15
    ) {
        int32_t w;
        if (Warning) w = -1; else w = 0;
        nonlinearoptimization_mp_trustregion_basic_(
            residue, Jacobian, x, M, N,
            w, MaxIteration, MaxStepIteration, Precision, MinStepLength
        );
    }

    inline void AugmentedLagrangian(
        void (*f)(double &, const double *, const int &),
        void (*fd)(double *, const double *, const int &),
        int (*f_fd)(double &, double *, const double *, const int &),
        int (*fdd)(double *, const double *, const int &),
        void (*c)(double *, const double *, const int &, const int &),
        void (*cd)(double *, const double *, const int &, const int &),
        int (*cdd)(double *, const double *, const int &, const int &),
        double * x, const int & N, const int & M,
        const std::string & UnconstrainedSolver="BFGS", std::vector<double> lambda0={}, const double & miu0=1.0,
        const int & ExactStep=20, const int & Memory=10, const std::string & Method="DY",
        const bool & Strong=true, const bool & Warning=true,
        const int & MaxIteration=1000, const double & Precision=1e-15, const double & MinStepLength=1e-15,
        const double & WolfeConst1=1e-4, double WolfeConst2=0.9, const double & Increment=1.05
    ) {
        if (lambda0.empty()) {
            lambda0.resize(M);
            std::fill(lambda0.begin(), lambda0.end(), 0.0);
        }
        int32_t s, w;
        if (Strong ) s = -1; else s = 0;
        if (Warning) w = -1; else w = 0;
        if (UnconstrainedSolver == "ConjugateGradient" && WolfeConst2 == 0.9) WolfeConst2 = 0.45;
        nonlinearoptimization_mp_augmentedlagrangian_(
            f, fd, c, cd, x, N, M,
            UnconstrainedSolver.c_str(), lambda0.data(), miu0,
            fdd, cdd, ExactStep, Memory, Method.c_str(),
            f_fd,
            s, w, MaxIteration, Precision, MinStepLength, WolfeConst1, WolfeConst2, Increment,
            UnconstrainedSolver.size(), Method.size()
        );
    }
#elif __GNUC__
    extern "C" {
        void __nonlinearoptimization_MOD_steepestdescent(
            // Required argument
            void (*f)(double &, const double *, const int &),
            void (*fd)(double *, const double *, const int &),
            double * x, const int & dim,
            int (*f_fd)(double &, double *, const double *, const int &),
            // Common optional argument
            const int32_t & Strong, const int32_t & Warning,
            const int & MaxIteration,
            const double & Precision, const double & MinStepLength,
            const double & WolfeConst1, const double & WolfeConst2, const double & Increment
        );

        void __nonlinearoptimization_MOD_conjugategradient_basic(
            // Required argument
            void (*f)(double &, const double *, const int &),
            void (*fd)(double *, const double *, const int &),
            double * x, const int & dim,
            // Conjugate gradient method
            const char * Method,
            // Common optional argument
            const int32_t & Strong, const int32_t & Warning,
            const int & MaxIteration,
            const double & Precision, const double & MinStepLength,
            const double & WolfeConst1, const double & WolfeConst2, const double & Increment,
            // Fortran requires string length
            int len_Method
        );

        void __nonlinearoptimization_MOD_conjugategradient(
            // Required argument
            void (*f)(double &, const double *, const int &),
            void (*fd)(double *, const double *, const int &),
            double * x, const int & dim,
            // Conjugate gradient method
            const char * Method,
            // Required argument
            int (*f_fd)(double &, double *, const double *, const int &),
            // Common optional argument
            const int32_t & Strong, const int32_t & Warning,
            const int & MaxIteration,
            const double & Precision, const double & MinStepLength,
            const double & WolfeConst1, const double & WolfeConst2, const double & Increment,
            // Fortran requires string length
            int len_Method
        );

        void __nonlinearoptimization_MOD_bfgs(
            // Required argument
            void (*f)(double &, const double *, const int &),
            void (*fd)(double *, const double *, const int &),
            double * x, const int & dim,
            int (*fdd)(double *, const double *, const int &),
            // Every how many steps compute exact hessian
            const int & ExactStep,
            // Required argument
            int (*f_fd)(double &, double *, const double *, const int &),
            // Common optional argument
            const int32_t & Strong, const int32_t & Warning,
            const int & MaxIteration,
            const double & Precision, const double & MinStepLength,
            const double & WolfeConst1, const double & WolfeConst2, const double & Increment
        );

        void __nonlinearoptimization_MOD_newtonraphson(
            // Required argument
            void (*f)(double &, const double *, const int &),
            void (*fd)(double *, const double *, const int &),
            double * x, const int & dim,
            int (*fdd)(double *, const double *, const int &),
            // Required argument
            int (*f_fd)(double &, double *, const double *, const int &),
            // Common optional argument
            const int32_t & Strong, const int32_t & Warning,
            const int & MaxIteration,
            const double & Precision, const double & MinStepLength,
            const double & WolfeConst1, const double & WolfeConst2, const double & Increment
        );

        void __nonlinearoptimization_MOD_trustregion_basic(
            // Required argument
            void (*residue)(double *, const double *, const int &, const int &),
            void (*Jacobian)(double *, const double *, const int &, const int &),
            double * x, const int & M, const int & N,
            // Optional argument
            const int32_t & Warning,
            const int & MaxIteration, const int & MaxStepIteration,
            const double & Precision, const double & MinStepLength
        );

        void __nonlinearoptimization_MOD_augmentedlagrangian(
            // Required argument
            void (*f)(double &, const double *, const int &),
            void (*fd)(double *, const double *, const int &),
            void (*c)(double *, const double *, const int &, const int &),
            void (*cd)(double *, const double *, const int &, const int &),
            double * x, const int & N, const int & M,
            // Augmented Lagrangian parameters
            const char * UnconstrainedSolver, const double * lambda0, const double & miu0,
            // Required argument
            int (*fdd)(double *, const double *, const int &),
            int (*cdd)(double *, const double *, const int &, const int &),
            // Unconstrained solver parameters
            const int & ExactStep, const int & Memory, const char * Method,
            // Required argument
            int (*f_fd)(double &, double *, const double *, const int &),
            // Common optional argument
            const int32_t & Strong, const int32_t & Warning,
            const int & MaxIteration,
            const double & Precision, const double & MinStepLength,
            const double & WolfeConst1, const double & WolfeConst2, const double & Increment,
            // Fortran requires string length
            int len_UnconstrainedSolver, int len_Method
        );
    }

    inline void SteepestDescent(
        void (*f)(double &, const double *, const int &),
        void (*fd)(double *, const double *, const int &),
        int (*f_fd)(double &, double *, const double *, const int &),
        double * x, const int & dim,
        const bool & Strong=true, const bool & Warning=true,
        const int & MaxIteration=1000, const double & Precision=1e-15, const double & MinStepLength=1e-15,
        const double & WolfeConst1=1e-4, const double & WolfeConst2=0.9, const double & Increment=1.05
    ) {
        int32_t s, w;
        if (Strong ) s = -1; else s = 0;
        if (Warning) w = -1; else w = 0;
        __nonlinearoptimization_MOD_steepestdescent(
            f, fd, x, dim,
            f_fd,
            s, w, MaxIteration, Precision, MinStepLength, WolfeConst1, WolfeConst2, Increment
        );
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
        int32_t s, w;
        if (Strong ) s = -1; else s = 0;
        if (Warning) w = -1; else w = 0;
        __nonlinearoptimization_MOD_conjugategradient_basic(
            f, fd, x, dim,
            Method.c_str(),
            s, w, MaxIteration, Precision, MinStepLength, WolfeConst1, WolfeConst2, Increment,
            Method.size()
        );
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
        int32_t s, w;
        if (Strong ) s = -1; else s = 0;
        if (Warning) w = -1; else w = 0;
        __nonlinearoptimization_MOD_conjugategradient(
            f, fd, x, dim,
            Method.c_str(),
            f_fd,
            s, w, MaxIteration, Precision, MinStepLength, WolfeConst1, WolfeConst2, Increment,
            Method.size()
        );
    }

    inline void BFGS(
        void (*f)(double &, const double *, const int &),
        void (*fd)(double *, const double *, const int &),
        int (*f_fd)(double &, double *, const double *, const int &),
        int (*fdd)(double *, const double *, const int &),
        double * x, const int & dim,
        const int & ExactStep=20,
        const bool & Strong=true, const bool & Warning=true,
        const int & MaxIteration=1000, const double & Precision=1e-15, const double & MinStepLength=1e-15,
        const double & WolfeConst1=1e-4, const double & WolfeConst2=0.9, const double & Increment=1.05
    ) {
        int32_t s, w;
        if (Strong ) s = -1; else s = 0;
        if (Warning) w = -1; else w = 0;
        __nonlinearoptimization_MOD_bfgs(
            f, fd, x, dim,
            fdd, ExactStep,
            f_fd,
            s, w, MaxIteration, Precision, MinStepLength, WolfeConst1, WolfeConst2, Increment
        );
    }

    inline void NewtonRaphson(
        void (*f)(double &, const double *, const int &),
        void (*fd)(double *, const double *, const int &),
        int (*f_fd)(double &, double *, const double *, const int &),
        int (*fdd)(double *, const double *, const int &),
        double * x, const int & dim,
        const bool & Strong=true, const bool & Warning=true,
        const int & MaxIteration=1000, const double & Precision=1e-15, const double & MinStepLength=1e-15,
        const double & WolfeConst1=1e-4, const double & WolfeConst2=0.9, const double & Increment=1.05
    ) {
        int32_t s, w;
        if (Strong ) s = -1; else s = 0;
        if (Warning) w = -1; else w = 0;
        __nonlinearoptimization_MOD_newtonraphson(
            f, fd, x, dim,
            fdd, f_fd,
            s, w, MaxIteration, Precision, MinStepLength, WolfeConst1, WolfeConst2, Increment
        );
    }

    inline void TrustRegion(
        void (*residue)(double *, const double *, const int &, const int &),
        void (*Jacobian)(double *, const double *, const int &, const int &),
        double * x, const int & M, const int & N,
        const bool & Warning=true,
        const int & MaxIteration=1000, const int & MaxStepIteration=100,
        const double & Precision=1e-15, const double & MinStepLength=1e-15
    ) {
        int32_t w;
        if (Warning) w = -1; else w = 0;
        __nonlinearoptimization_MOD_trustregion_basic(
            residue, Jacobian, x, M, N,
            w, MaxIteration, MaxStepIteration, Precision, MinStepLength
        );
    }

    inline void AugmentedLagrangian(
        void (*f)(double &, const double *, const int &),
        void (*fd)(double *, const double *, const int &),
        int (*f_fd)(double &, double *, const double *, const int &),
        int (*fdd)(double *, const double *, const int &),
        void (*c)(double *, const double *, const int &, const int &),
        void (*cd)(double *, const double *, const int &, const int &),
        int (*cdd)(double *, const double *, const int &, const int &),
        double * x, const int & N, const int & M,
        const std::string & UnconstrainedSolver="BFGS", std::vector<double> lambda0={}, const double & miu0=1.0,
        const int & ExactStep=20, const int & Memory=10, const std::string & Method="DY",
        const bool & Strong=true, const bool & Warning=true,
        const int & MaxIteration=1000, const double & Precision=1e-15, const double & MinStepLength=1e-15,
        const double & WolfeConst1=1e-4, double WolfeConst2=0.9, const double & Increment=1.05
    ) {
        if (lambda0.empty()) {
            lambda0.resize(M);
            std::fill(lambda0.begin(), lambda0.end(), 0.0);
        }
        int32_t s, w;
        if (Strong ) s = -1; else s = 0;
        if (Warning) w = -1; else w = 0;
        if (UnconstrainedSolver == "ConjugateGradient" && WolfeConst2 == 0.9) WolfeConst2 = 0.45;
        __nonlinearoptimization_MOD_augmentedlagrangian(
            f, fd, c, cd, x, N, M,
            UnconstrainedSolver.c_str(), lambda0.data(), miu0,
            fdd, cdd, ExactStep, Memory, Method.c_str(),
            f_fd,
            s, w, MaxIteration, Precision, MinStepLength, WolfeConst1, WolfeConst2, Increment,
            UnconstrainedSolver.size(), Method.size()
        );
    }
#endif

} // namespace NO
} // namespace FL

#endif