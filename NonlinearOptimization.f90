!Nonlinear optimization routines
!
!Instruction:
!If you want to solve nonlinear equation(s), you may:
!    adopt MKL trust region solver
!    or build your own merit function then use other routines to search for its minimum
!This module mainly provides unconstrained local minimizer
!    Only solvers in Augmented Lagrangian section are constrained
!    Only solvers in Heuristic section are global
!Global optimization remains an unsolved problem, since it is NP-complete
!    2-step method is my favourite approach for global optimization:
!        Step 1: A fast but inexact hopper to explore the phase space
!        Step 2: A rigorous local optimizer to polish the best estimation explored
!    How to build a hopper is highly empirical, examples are SurfGen and cluster structure optimization (basin-hopping),
!    so you have to construct your own hopper based on your specific problem
!    Heuristic algorithm is general but naive: pay exponentially for NP-completeness
module NonlinearOptimization
    use LinearAlgebra!NewtonRaphson and BFGS need it
    use mkl_rci_type !BFGS_NumericalHessian and Trust region section need it
    use mkl_rci      !BFGS_NumericalHessian and Trust region section need it
    implicit none

!Parameter
    real*8::jacobiPrecision=1d-8!DO NOT be too small: this might be the finite difference step length
    !Line search
        real*8::LineSearchStepLowerBound=1d-15
        !Newton-Raphson
            logical::NewtonRaphsonWarning=.true.
            integer::MaxNewtonRaphsonIteration=1000
            real*8::NewtonRaphsonTol=1d-30!Tolerence for 2-norm square of gradient
        !Quasi-Newton
            logical::QuasiNewtonWarning=.true.
            integer::MaxQuasiNewtonIteration=1000
            real*8::QuasiNewtonTol=1d-30!Tolerence for 2-norm square of gradient
        !Conjugate gradient
            logical::ConjugateGradientWarning=.true.
            integer::MaxConjugateGradientIteration=1000
            real*8::ConjugateGradientTol=1d-30!Tolerence for 2-norm square of gradient
    !Trust region
    !Heuristic
    !MKL solver
        logical::trnlspWarning=.true.
        integer::MaxMKLTrustRegionIteration=1000,MaxMKLTrialStepIteration=100
        !MKLTrustRegionTol(1:5) contains the stopping criteria for solving f'(x) = 0:
        !    1, trust region radius < MKLTrustRegionTol(1)
        !    2, || f'(x) ||_2 < MKLTrustRegionTol(2)
        !    3, for all j in [1,N], || Jacobian(:,j) ||_2 < MKLTrustRegionTol(3) 
        !    4, || s ||_2 < MKLTrustRegionTol(4), where s is the trial step
        !    5, || f'(x) ||_2 - || f'(x) - Jacobian . s ||_2 < MKLTrustRegionTol(5)
        !MKLTrustRegionTol(6) is the precision of s calculation
        real*8,dimension(6)::MKLTrustRegionTol=[1d-15,1d-15,1d-15,1d-15,1d-15,1d-15]

contains
!-------------- Line search ---------------
    !Nomenclature:
    !    f = the target function to be minimized
    !    a = the line search step length
    !    p = the line search direction
    !    phi(a) = f( x + a * p ), so phi'(a) = f'( x + a * p ) . p
    !Input format: x and f'(x) are dim dimensional vectors, f''(x) is dim order matrix
    !    subroutine f(f(x),x,dim)
    !    subroutine fd(f'(x),x,dim)
    !    subroutine f_fd(f(x),f'(x),x,dim)
    !    subroutine fdd(f''(x),x,dim)

    !=========== Line searcher ============
        !Some direction finders can only converge under strong Wolfe condition
        !Goldstein is only suitable for Newton, but not even for quasi-Newton
        !For Newton and quasi-Newton, the initial guess can always be a = 1, because their direction vector is well scaled
        !However, for methods whose direction is not determined by inverted (approximate) Hessian multiplying -gradient,
        !e.g., conjugate gradient and stochastic gradient descent, user has to come up with a good initial guess

        !Line search for a step length satisfying Wolfe condition
        !This subroutine is designed to minimize gradient computation
        !Input: c1 & c2 are Wolfe constants (0<c1<c2<1), x is current x
        !       a is initial guess of a, p is current p, fx = f(x), phid0 = phi'(0)
        !Output: a harvests the step length satisfying strong Wolfe condition, x = x + a * p, fx = f(x), fdx = f'(x)
        subroutine Wolfe(c1,c2,f,fd,x,a,p,fx,phid0,fdx,dim)
            real*8,intent(in)::c1,c2
            external::f,fd
            integer,intent(in)::dim
            real*8,dimension(dim),intent(inout)::x
            real*8,intent(inout)::a
            real*8,dimension(dim),intent(in)::p
            real*8,intent(inout)::fx
            real*8,intent(in)::phid0
            real*8,dimension(dim),intent(out)::fdx
            real*8,parameter::increment=1.05d0! > 1
            real*8::c2_m_phid0,fx0,atemp,aold,fold,phidx
            real*8,dimension(dim)::x0
            !Initialize
                x0=x
                fx0=fx
                c2_m_phid0=c2*phid0
            !Check whether initial guess satisfies sufficient decrease condition
            x=x0+a*p
            call f(fx,x,dim)
            if(fx<=fx0+c1*a*phid0) then!Satisfied, search for larger a
                do
                    aold=a
                    fold=fx
                    a=aold*increment
                    x=x0+a*p
                    call f(fx,x,dim)
                    if(fx>fx0+c1*a*phid0) then
                        x=x0+aold*p
                        call fd(fdx,x,dim)
                        phidx=dot_product(fdx,p)
                        if(phidx>c2_m_phid0) then
                            a=aold
                            fx=fold
                        else
                            x=x0
                            atemp=a
                            call WolfeZoom(c1,c2,f,fd,x,a,p,fx0,phid0,aold,atemp,fold,fx,phidx,fdx,dim)
                            fx=fx0
                        end if
                        return
                    end if
                end do
            else!Violated, search for smaller a
                do
                    aold=a
                    fold=fx
                    a=aold/increment
                    x=x0+a*p
                    call f(fx,x,dim)
                    if(fx<=fx0+c1*a*phid0) then
                        call fd(fdx,x,dim)
                        phidx=dot_product(fdx,p)
                        if(phidx<c2_m_phid0) then
                            x=x0
                            atemp=a
                            call WolfeZoom(c1,c2,f,fd,x,a,p,fx0,phid0,atemp,aold,fx,fold,phidx,fdx,dim)
                            fx=fx0
                        end if
                        return
                    end if
                    if(a<LineSearchStepLowerBound) then
                        call fd(fdx,x,dim)
                        return
                    end if
                end do
            end if
        end subroutine Wolfe
        !Support Wolfe. low and up must satisfy:
        !    low < up
        !    low satisfies sufficient decrease condition, but up violates
        !    phi'(low) < 0
        !Input: c1 & c2 are Wolfe constants (0<c1<c2<1), x is current x, p is current p, fx = f(x), phid0 = phi'(0)
        !       low & up are explained above, flow/up & phidlow/up corresponds to low/up
        !Output: a harvests the step length satisfying strong Wolfe condition, x = x + a * p, fx = f(x), fdx = f'(x)
        subroutine WolfeZoom(c1,c2,f,fd,x,a,p,fx,phid0,low,up,flow,fup,phidlow,fdx,dim)
            real*8,intent(in)::c1,c2! 0 < c1 < c2 < 1
            external::f,fd
            integer,intent(in)::dim
            real*8,dimension(dim),intent(inout)::x
            real*8,intent(out)::a
            real*8,dimension(dim),intent(in)::p
            real*8,intent(inout)::fx
            real*8,intent(in)::phid0
            real*8,intent(inout)::low,up,flow,fup,phidlow
            real*8,dimension(dim),intent(out)::fdx
            real*8::c2_m_phid0,fx0,phidnew,phidlow_m_a
            real*8,dimension(dim)::x0
            !Initialize
                x0=x
                fx0=fx
                c2_m_phid0=c2*phid0
                phidlow_m_a=phidlow*a
            do
                !Updata a by quadratic interpolation
                    a=phidlow_m_a*a/2d0/(flow+phidlow_m_a-fup)
                    if(.not.(a>up.and.a<low)) a=(low+up)/2d0
                x=x0+a*p
                call f(fx,x,dim)
                if(fx>fx0+c1*a*phid0) then
                    up=a
                    fup=fx
                else
                    call fd(fdx,x,dim)
                    phidnew=dot_product(fdx,p)
                    if(phidnew>c2_m_phid0) return
                    low=a
                    flow=fx
                    phidlow=phidnew
                    phidlow_m_a=phidlow*a
                end if
                if(Abs(up-low)<LineSearchStepLowerBound) then
                    call fd(fdx,x,dim)
                    return
                end if
            end do
        end subroutine WolfeZoom

        !Line search for a step length satisfying strong Wolfe condition
        !Input: c1 & c2 are Wolfe constants (0<c1<c2<1), x is current x
        !       a is initial guess of a, p is current p, fx = f(x), phid0 = phi'(0)
        !Output: a harvests the step length satisfying strong Wolfe condition, x = x + a * p, fx = f(x), fdx = f'(x)
        subroutine StrongWolfe(c1,c2,f,fd,x,a,p,fx,phid0,fdx,dim)
            real*8,intent(in)::c1,c2
            external::f,fd
            integer,intent(in)::dim
            real*8,dimension(dim),intent(inout)::x
            real*8,intent(inout)::a
            real*8,dimension(dim),intent(in)::p
            real*8,intent(inout)::fx
            real*8,intent(in)::phid0
            real*8,dimension(dim),intent(out)::fdx
            real*8,parameter::increment=1.05d0! > 1
            real*8::c2_m_abs_phid0,fx0,atemp,aold,fold,phidnew,phidold
            real*8,dimension(dim)::x0
            !Initialize
                x0=x
                fx0=fx
                c2_m_abs_phid0=c2*Abs(phid0)
            !Check whether initial guess satisfies sufficient decrease condition
            x=x0+a*p
            call f(fx,x,dim)
            if(fx<=fx0+c1*a*phid0) then!Satisfied, try to search for larger a
                call fd(fdx,x,dim)
                phidnew=dot_product(fdx,p)
                if(phidnew>0d0) then!Curve is heading up
                    if(Abs(phidnew)<=c2_m_abs_phid0) return
                    do!Else have to search for smaller a that phi(a) < phi(aold) & phid(a) > 0 is false
                        aold=a
                        fold=fx
                        phidold=phidnew
                        a=aold/increment
                        x=x0+a*p
                        call f(fx,x,dim)
                        call fd(fdx,x,dim)
                        phidnew=dot_product(fdx,p)
                        if(fx>=fold.or.phidnew<=0d0) then
                            x=x0
                            atemp=a
                            call StrongWolfeZoom(c1,c2,f,fd,x,a,p,fx0,phid0,aold,atemp,fold,fx,phidold,phidnew,fdx,dim)
                            fx=fx0
                            return
                        end if
                        if(a<LineSearchStepLowerBound) return
                    end do
                else!Search for larger a
                    do
                        aold=a
                        fold=fx
                        phidold=phidnew
                        a=aold*increment
                        x=x0+a*p
                        call f(fx,x,dim)
                        call fd(fdx,x,dim)
                        phidnew=dot_product(fdx,p)
                        if(fx>fx0+c1*a*phid0.or.fx>=fold) then
                            x=x0
                            atemp=a
                            call StrongWolfeZoom(c1,c2,f,fd,x,a,p,fx0,phid0,aold,atemp,fold,fx,phidold,phidnew,fdx,dim)
                            fx=fx0
                            return
                        end if
                        if(phidnew>0d0) then
                            if(Abs(phidnew)<=c2_m_abs_phid0) then
                                return
                            else
                                x=x0
                                atemp=a
                                call StrongWolfeZoom(c1,c2,f,fd,x,a,p,fx0,phid0,atemp,aold,fx,fold,phidnew,phidold,fdx,dim)
                                fx=fx0
                                return
                            end if
                        end if
                    end do
                end if
            else!Violated, 1st search for smaller a satisfying sufficient decrease condition
                do
                    aold=a
                    fold=fx
                    a=aold/increment
                    x=x0+a*p
                    call f(fx,x,dim)
                    if(fx<=fx0+c1*a*phid0) then!Found, then look at slope
                        call fd(fdx,x,dim)
                        phidnew=dot_product(fdx,p)
                        if(Abs(phidnew)<=c2_m_abs_phid0) return
                        if(phidnew<0d0) then!Within [a, aold]
                            x=x0+aold*p
                            call fd(fdx,x,dim)
                            phidold=dot_product(fdx,p)
                            x=x0
                            atemp=a
                            call StrongWolfeZoom(c1,c2,f,fd,x,a,p,fx0,phid0,atemp,aold,fx,fold,phidnew,phidold,fdx,dim)
                            fx=fx0
                            return
                        else!Search for such an a that phi(a) < phi(aold) & phid(a) > 0 is false
                            do
                                aold=a
                                fold=fx
                                phidold=phidnew
                                a=aold/increment
                                x=x0+a*p
                                call f(fx,x,dim)
                                call fd(fdx,x,dim)
                                phidnew=dot_product(fdx,p)
                                if(fx>=fold.or.phidnew<=0d0) then
                                    x=x0
                                    atemp=a
                                    call StrongWolfeZoom(c1,c2,f,fd,x,a,p,fx0,phid0,aold,atemp,fold,fx,phidold,phidnew,fdx,dim)
                                    fx=fx0
                                    return
                                end if
                                if(a<LineSearchStepLowerBound) return
                            end do
                        end if
                    end if
                    if(a<LineSearchStepLowerBound) then
                        call fd(fdx,x,dim)
                        return
                    end if
                end do
            end if
        end subroutine StrongWolfe
        !Support StrongWolfe. low and up must satisfy:
        !    low satisfies sufficient decrease condition
        !    ( up - low ) * phi'(low) < 0
        !    [ low, up ] (or [ up, low ]) contains a step length satisfying strong Wolfe condition
        !        This means at least 1 of 3 following statements is true:
        !            up violates the sufficient decrease condition
        !            phi(up) >= phi(low)
        !            up < low & phi'(up) <= 0
        !Input: c1 & c2 are Wolfe constants (0<c1<c2<1), x is current x, p is current p, fx = f(x), phid0 = phi'(0)
        !       low & up are explained above, flow/up & phidlow/up corresponds to low/up
        !Output: a harvests the step length satisfying strong Wolfe condition, x = x + a * p, fx = f(x), fdx = f'(x)
        subroutine StrongWolfeZoom(c1,c2,f,fd,x,a,p,fx,phid0,low,up,flow,fup,phidlow,phidup,fdx,dim)
            real*8,intent(in)::c1,c2
            external::f,fd
            integer,intent(in)::dim
            real*8,dimension(dim),intent(inout)::x
            real*8,intent(out)::a
            real*8,dimension(dim),intent(in)::p
            real*8,intent(inout)::fx
            real*8,intent(in)::phid0
            real*8,intent(inout)::low,up,flow,fup,phidlow,phidup
            real*8,dimension(dim),intent(out)::fdx
            real*8::c2_m_abs_phid0,fx0,phidnew,d1,d2
            real*8,dimension(dim)::x0
            !Initialize
                x0=x
                fx0=fx
                c2_m_abs_phid0=c2*Abs(phid0)
            do
                !Updata a by cubic interpolation
                    d1=phidlow+phidup-3d0*(flow-fup)/(low-up)
                    d2=up-low
                    if(d2>0d0) then
                        d2=Sqrt(d1*d1-phidlow*phidup)
                    else
                        d2=-Sqrt(d1*d1-phidlow*phidup)
                    end if
                    a=up-(up-low)*(phidup+d2-d1)/(phidup-phidlow+2d0*d2)
                    if(.not.(a>min(low,up).and.a<max(low,up))) a=(low+up)/2d0
                x=x0+a*p
                call f(fx,x,dim)
                call fd(fdx,x,dim)
                phidnew=dot_product(fdx,p)
                if(fx>fx0+c1*a*phid0.or.fx>=flow) then
                    up=a
                    fup=fx
                    phidup=phidnew
                else
                    if(Abs(phidnew)<=c2_m_abs_phid0) return
                    if(phidnew*(up-low)>=0d0) then
                        up=low
                        fup=flow
                        phidup=phidlow
                    end if
                    low=a
                    flow=fx
                    phidlow=phidnew
                end if
                if(Abs(up-low)<LineSearchStepLowerBound) return
            end do
        end subroutine StrongWolfeZoom
        
        !Almost same to StrongWolfe, but modified for cheap to evaluate f' along with f case
        !f_fd has the form of: subroutine f_fd(f(x),f'(x),x,dim)
        subroutine StrongWolfe_fdwithf(c1,c2,f,fd,f_fd,x,a,p,fx,phid0,fdx,dim)
            real*8,intent(in)::c1,c2
            external::f,fd,f_fd
            integer,intent(in)::dim
            real*8,dimension(dim),intent(inout)::x
            real*8,intent(inout)::a
            real*8,dimension(dim),intent(in)::p
            real*8,intent(inout)::fx
            real*8,intent(in)::phid0
            real*8,dimension(dim),intent(out)::fdx
            real*8,parameter::increment=1.05d0! > 1
            real*8::c2_m_abs_phid0,fx0,atemp,aold,fold,phidnew,phidold
            real*8,dimension(dim)::x0
            !Initialize
                x0=x
                fx0=fx
                c2_m_abs_phid0=c2*Abs(phid0)
            !Check whether initial guess satisfies sufficient decrease condition
            x=x0+a*p
            call f_fd(fx,fdx,x,dim)
            if(fx<=fx0+c1*a*phid0) then!Satisfied, try to search for larger a
                phidnew=dot_product(fdx,p)
                if(phidnew>0d0) then!Curve is heading up
                    if(Abs(phidnew)<=c2_m_abs_phid0) return
                    do!Else have to search for smaller a that phi(a) < phi(aold) & phid(a) > 0 is false
                        aold=a
                        fold=fx
                        phidold=phidnew
                        a=aold/increment
                        x=x0+a*p
                        call f_fd(fx,fdx,x,dim)
                        phidnew=dot_product(fdx,p)
                        if(fx>=fold.or.phidnew<=0d0) then
                            x=x0
                            atemp=a
                            call StrongWolfeZoom_fdwithf(c1,c2,f_fd,x,a,p,fx0,phid0,aold,atemp,fold,fx,phidold,phidnew,fdx,dim)
                            fx=fx0
                            return
                        end if
                        if(a<LineSearchStepLowerBound) return
                    end do
                else!Search for larger a
                    do
                        aold=a
                        fold=fx
                        phidold=phidnew
                        a=aold*increment
                        x=x0+a*p
                        call f_fd(fx,fdx,x,dim)
                        phidnew=dot_product(fdx,p)
                        if(fx>fx0+c1*a*phid0.or.fx>=fold) then
                            x=x0
                            atemp=a
                            call StrongWolfeZoom_fdwithf(c1,c2,f_fd,x,a,p,fx0,phid0,aold,atemp,fold,fx,phidold,phidnew,fdx,dim)
                            fx=fx0
                            return
                        end if
                        if(phidnew>0d0) then
                            if(Abs(phidnew)<=c2_m_abs_phid0) return
                            x=x0
                            atemp=a
                            call StrongWolfeZoom_fdwithf(c1,c2,f_fd,x,a,p,fx0,phid0,atemp,aold,fx,fold,phidnew,phidold,fdx,dim)
                            fx=fx0
                            return
                        end if
                    end do
                end if
            else!Violated, 1st search for smaller a satisfying sufficient decrease condition
                do
                    aold=a
                    fold=fx
                    a=aold/increment
                    x=x0+a*p
                    call f(fx,x,dim)
                    if(fx<=fx0+c1*a*phid0) then!Found, then look at slope
                        call fd(fdx,x,dim)
                        phidnew=dot_product(fdx,p)
                        if(Abs(phidnew)<=c2_m_abs_phid0) return
                        if(phidnew<0d0) then!Within [a, aold]
                            x=x0+aold*p
                            call fd(fdx,x,dim)
                            phidold=dot_product(fdx,p)
                            x=x0
                            atemp=a
                            call StrongWolfeZoom_fdwithf(c1,c2,f_fd,x,a,p,fx0,phid0,atemp,aold,fx,fold,phidnew,phidold,fdx,dim)
                            fx=fx0
                            return
                        else!Search for such an a that phi(a) < phi(aold) & phid(a) > 0 is false
                            do
                                aold=a
                                fold=fx
                                phidold=phidnew
                                a=aold/increment
                                x=x0+a*p
                                call f_fd(fx,fdx,x,dim)
                                phidnew=dot_product(fdx,p)
                                if(fx>=fold.or.phidnew<=0d0) then
                                    x=x0
                                    atemp=a
                                    call StrongWolfeZoom_fdwithf(c1,c2,f_fd,x,a,p,fx0,phid0,aold,atemp,fold,fx,phidold,phidnew,fdx,dim)
                                    fx=fx0
                                    return
                                end if
                                if(a<LineSearchStepLowerBound) return
                            end do
                        end if
                    end if
                    if(a<LineSearchStepLowerBound) then
                        call fd(fdx,x,dim)
                        return
                    end if
                end do
            end if
        end subroutine StrongWolfe_fdwithf
        !Support StrongWolfe_fdwithf
        subroutine StrongWolfeZoom_fdwithf(c1,c2,f_fd,x,a,p,fx,phid0,low,up,flow,fup,phidlow,phidup,fdx,dim)
            real*8,intent(in)::c1,c2
            external::f_fd
            integer,intent(in)::dim
            real*8,dimension(dim),intent(inout)::x
            real*8,intent(out)::a
            real*8,dimension(dim),intent(in)::p
            real*8,intent(inout)::fx
            real*8,intent(in)::phid0
            real*8,intent(inout)::low,up,flow,fup,phidlow,phidup
            real*8,dimension(dim),intent(out)::fdx
            real*8::c2_m_abs_phid0,fx0,phidnew,d1,d2
            real*8,dimension(dim)::x0
            !Initialize
                x0=x
                fx0=fx
                c2_m_abs_phid0=c2*Abs(phid0)
            do
                !Updata a by cubic interpolation
                    d1=phidlow+phidup-3d0*(flow-fup)/(low-up)
                    d2=up-low
                    if(d2>0d0) then
                        d2=Sqrt(d1*d1-phidlow*phidup)
                    else
                        d2=-Sqrt(d1*d1-phidlow*phidup)
                    end if
                    a=up-(up-low)*(phidup+d2-d1)/(phidup-phidlow+2d0*d2)
                    if(.not.(a>min(low,up).and.a<max(low,up))) a=(low+up)/2d0
                x=x0+a*p
                call f_fd(fx,fdx,x,dim)
                phidnew=dot_product(fdx,p)
                if(fx>fx0+c1*a*phid0.or.fx>=flow) then
                    up=a
                    fup=fx
                    phidup=phidnew
                else
                    if(Abs(phidnew)<=c2_m_abs_phid0) return
                    if(phidnew*(up-low)>=0d0) then
                        up=low
                        fup=flow
                        phidup=phidlow
                    end if
                    low=a
                    flow=fx
                    phidlow=phidnew
                end if
                if(Abs(up-low)<LineSearchStepLowerBound) return
            end do
        end subroutine StrongWolfeZoom_fdwithf        

        !I will write someday (flag)
        subroutine Goldstein(c,f,x,p,phid0,dim)
            real*8,intent(in)::c! 0 < c < 0.5
            external::f
            integer,intent(in)::dim
            real*8,dimension(dim),intent(inout)::x
            real*8,dimension(dim),intent(in)::p
            real*8,intent(in)::phid0
        end subroutine Goldstein
    !================ End =================

    !========== Direction finder ==========
        !If dimensionality is low, adopt quasi-Newton (or Newton if Hessian is cheap and initial guess is close)
        !If dimensionality is so high that O(dim^2) memory is unaffordable, adopt conjugate gradient or L-BFGS

        !On input x is an initial guess of the minimum point of f(x), on exit x is the minimum point
        !Newton-Raphson method, requiring Wolfe condition
        subroutine NewtonRaphson(f,fd,fdd,x,dim)
            external::f,fd,fdd
            integer,intent(in)::dim
            real*8,dimension(dim),intent(inout)::x
            real*8,parameter::c1=1d-4,c2=0.9d0! 0 < c1 < c2 < 1: Wolfe constant
            integer::iIteration,info
            real*8::a,fnew,phidnew,phidold
            real*8,dimension(dim)::p,fdnew
            real*8,dimension(dim,dim)::Hessian
            !Initialize
                call f(fnew,x,dim)
                call fd(fdnew,x,dim)
                call fdd(Hessian,x,dim)
                p=-fdnew
                call My_dposv_poQuery(Hessian,p,dim,info)
                if(info==0) then
                    phidnew=dot_product(fdnew,p)
                    a=1d0
                else!Hessian is not positive definite, use steepest descent direction
                    p=-fdnew
                    phidnew=-dot_product(fdnew,fdnew)
                    if(-phidnew<QuasiNewtonTol) return
                    if(fnew==0d0) then
                        a=1d0
                    else
                        a=-fnew/phidnew
                    end if
                end if
            do iIteration=1,MaxNewtonRaphsonIteration
                !Prepare
                phidold=phidnew
                !Line search
                call Wolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim)
                phidnew=dot_product(fdnew,fdnew)
                if(phidnew<NewtonRaphsonTol) return
                if(dot_product(p,p)*a*a<NewtonRaphsonTol) then
                    if(NewtonRaphsonWarning) then
                        write(*,'(1x,A94)')'Newton-Raphson warning: step length has converged, but gradient norm has not met accuracy goal'
                        write(*,'(1x,A56)')'A best estimation rather than exact solution is returned'
                        write(*,*)'Euclidean norm of gradient =',Sqrt(phidnew)
                        NewtonRaphsonWarning=.false.
                    end if
                    return
                end if
                !Determine new direction
                p=-fdnew
                call My_dposv_poQuery(Hessian,p,dim,info)
                if(info==0) then
                    phidnew=dot_product(fdnew,p)
                    a=1d0
                else!Hessian is not positive definite, use steepest descent direction
                    p=-fdnew
                    phidnew=-phidnew
                    a=a*phidold/phidnew
                end if
            end do
            if(iIteration>MaxNewtonRaphsonIteration.and.NewtonRaphsonWarning) then
                write(*,*)'Failed Newton-Raphson: max iteration exceeded!'
                write(*,*)'Euclidean norm of gradient =',Norm2(fdnew)
                NewtonRaphsonWarning=.false.
            end if
        end subroutine NewtonRaphson
        !Strong Wolfe condition usually performs better
        subroutine NewtonRaphson_Strong(f,fd,fdd,x,dim)
            external::f,fd,fdd
            integer,intent(in)::dim
            real*8,dimension(dim),intent(inout)::x
            real*8,parameter::c1=1d-4,c2=0.9d0! 0 < c1 < c2 < 1: Wolfe constant
            integer::iIteration,info
            real*8::a,fnew,phidnew,phidold
            real*8,dimension(dim)::p,fdnew
            real*8,dimension(dim,dim)::Hessian
            !Initialize
                call f(fnew,x,dim)
                call fd(fdnew,x,dim)
                call fdd(Hessian,x,dim)
                p=-fdnew
                call My_dposv_poQuery(Hessian,p,dim,info)
                if(info==0) then
                    phidnew=dot_product(fdnew,p)
                    a=1d0
                else!Hessian is not positive definite, use steepest descent direction
                    p=-fdnew
                    phidnew=-dot_product(fdnew,fdnew)
                    if(-phidnew<QuasiNewtonTol) return
                    if(fnew==0d0) then
                        a=1d0
                    else
                        a=-fnew/phidnew
                    end if
                end if
            do iIteration=1,MaxNewtonRaphsonIteration
                !Prepare
                phidold=phidnew
                !Line search
                call StrongWolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim)
                phidnew=dot_product(fdnew,fdnew)
                if(phidnew<NewtonRaphsonTol) return
                if(dot_product(p,p)*a*a<NewtonRaphsonTol) then
                    if(NewtonRaphsonWarning) then
                        write(*,'(1x,A94)')'Newton-Raphson warning: step length has converged, but gradient norm has not met accuracy goal'
                        write(*,'(1x,A56)')'A best estimation rather than exact solution is returned'
                        write(*,*)'Euclidean norm of gradient =',Sqrt(phidnew)
                        NewtonRaphsonWarning=.false.
                    end if
                    return 
                end if
                !Determine new direction
                p=-fdnew
                call My_dposv_poQuery(Hessian,p,dim,info)
                if(info==0) then
                    phidnew=dot_product(fdnew,p)
                    a=1d0
                else!Hessian is not positive definite, use steepest descent direction
                    p=-fdnew
                    phidnew=-phidnew
                    a=a*phidold/phidnew
                end if
            end do
            if(iIteration>MaxNewtonRaphsonIteration.and.NewtonRaphsonWarning) then
                write(*,*)'Failed Newton-Raphson: max iteration exceeded!'
                write(*,*)'Euclidean norm of gradient =',Norm2(fdnew)
                NewtonRaphsonWarning=.false.
            end if
        end subroutine NewtonRaphson_Strong
        !When it is cheap to evaluate f' along with f
        subroutine NewtonRaphson_Strong_fdwithf(f,fd,f_fd,fdd,x,dim)
            external::f,fd,f_fd,fdd
            integer,intent(in)::dim
            real*8,dimension(dim),intent(inout)::x
            real*8,parameter::c1=1d-4,c2=0.9d0! 0 < c1 < c2 < 1: Wolfe constant
            integer::iIteration,info
            real*8::a,fnew,phidnew,phidold
            real*8,dimension(dim)::p,fdnew
            real*8,dimension(dim,dim)::Hessian
            !Initialize
                call f_fd(fnew,fdnew,x,dim)
                call fdd(Hessian,x,dim)
                p=-fdnew
                call My_dposv_poQuery(Hessian,p,dim,info)
                if(info==0) then
                    phidnew=dot_product(fdnew,p)
                    a=1d0
                else!Hessian is not positive definite, use steepest descent direction
                    p=-fdnew
                    phidnew=-dot_product(fdnew,fdnew)
                    if(-phidnew<QuasiNewtonTol) return
                    if(fnew==0d0) then
                        a=1d0
                    else
                        a=-fnew/phidnew
                    end if
                end if
            do iIteration=1,MaxNewtonRaphsonIteration
                !Prepare
                phidold=phidnew
                !Line search
                call StrongWolfe_fdwithf(c1,c2,f,fd,f_fd,x,a,p,fnew,phidnew,fdnew,dim)
                phidnew=dot_product(fdnew,fdnew)
                if(phidnew<NewtonRaphsonTol) return
                if(dot_product(p,p)*a*a<NewtonRaphsonTol) then
                    if(NewtonRaphsonWarning) then
                        write(*,'(1x,A94)')'Newton-Raphson warning: step length has converged, but gradient norm has not met accuracy goal'
                        write(*,'(1x,A56)')'A best estimation rather than exact solution is returned'
                        write(*,*)'Euclidean norm of gradient =',Sqrt(phidnew)
                        NewtonRaphsonWarning=.false.
                    end if
                    return 
                end if
                !Determine new direction
                p=-fdnew
                call My_dposv_poQuery(Hessian,p,dim,info)
                if(info==0) then
                    phidnew=dot_product(fdnew,p)
                    a=1d0
                else!Hessian is not positive definite, use steepest descent direction
                    p=-fdnew
                    phidnew=-phidnew
                    a=a*phidold/phidnew
                end if
            end do
            if(iIteration>MaxNewtonRaphsonIteration.and.NewtonRaphsonWarning) then
                write(*,*)'Failed Newton-Raphson: max iteration exceeded!'
                write(*,*)'Euclidean norm of gradient =',Norm2(fdnew)
                NewtonRaphsonWarning=.false.
            end if
        end subroutine NewtonRaphson_Strong_fdwithf

        !On entry and every freq steps compute exact Hessian. 9 < freq < 51 is recommended 
        !On input x is an initial guess of the minimum point of f(x), on exit x is the minimum point
        !Broyden–Fletcher–Goldfarb–Shanno (BFGS) quasi-Newton method, requiring Wolfe condition
        subroutine BFGS(f,fd,fdd,x,dim,freq)
            external::f,fd,fdd
            integer,intent(in)::dim,freq
            real*8,dimension(dim),intent(inout)::x
            real*8,parameter::c1=1d-4,c2=0.9d0! 0 < c1 < c2 < 1: Wolfe constant
            integer::iIteration,i
            real*8::a,fnew,phidnew,rho
            real*8,dimension(dim)::p,fdnew,s,y
            real*8,dimension(dim,dim)::U,H!Approximate inverse Hessian
            !Initialize
                call f(fnew,x,dim)
                call fd(fdnew,x,dim)
                call fdd(H,x,dim)
                p=-fdnew
                call My_dpotri_poQuery(H,dim,i)
                if(i==0) then
                    call syL2U(H,dim)
                    p=-matmul(H,fdnew)
                    phidnew=dot_product(fdnew,p)
                    a=1d0
                else!Hessian is not positive definite, use a as initial approximate inverse Hessian
                    p=-fdnew
                    phidnew=-dot_product(fdnew,fdnew)
                    if(-phidnew<QuasiNewtonTol) return
                    if(fnew==0d0) then
                        a=1d0
                    else
                        a=-fnew/phidnew
                    end if
                    s=x
                    y=fdnew
                    call Wolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim)
                    phidnew=dot_product(fdnew,fdnew)
                    if(phidnew<QuasiNewtonTol) return
                    if(dot_product(p,p)*a*a<QuasiNewtonTol) then
                        if(QuasiNewtonWarning) then
                            write(*,'(1x,A84)')'BFGS warning: step length has converged, but gradient norm has not met accuracy goal'
                            write(*,'(1x,A56)')'A best estimation rather than exact solution is returned'
                            write(*,*)'Euclidean norm of gradient =',Sqrt(phidnew)
                            QuasiNewtonWarning=.false.
                        end if
                        return 
                    end if
                    s=x-s
                    y=fdnew-y
                    rho=1d0/dot_product(y,s)
                    U=-rho*vector_direct_product(y,s,dim,dim)
                    forall(i=1:dim)
                        U(i,i)=U(i,i)+1d0
                    end forall
                    H=a*U
                    H=matmul(transpose(U),H)+rho*vector_direct_product(s,s,dim,dim)
                    p=-matmul(H,fdnew)
                    phidnew=dot_product(fdnew,p)
                    a=1d0
                end if
            do iIteration=1,MaxQuasiNewtonIteration
                !Prepare
                s=x
                y=fdnew
                !Line search
                call Wolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim)
                phidnew=dot_product(fdnew,fdnew)
                if(phidnew<QuasiNewtonTol) return
                if(dot_product(p,p)*a*a<QuasiNewtonTol) then
                    if(QuasiNewtonWarning) then
                        write(*,'(1x,A84)')'BFGS warning: step length has converged, but gradient norm has not met accuracy goal'
                        write(*,'(1x,A56)')'A best estimation rather than exact solution is returned'
                        write(*,*)'Euclidean norm of gradient =',Sqrt(phidnew)
                        QuasiNewtonWarning=.false.
                    end if
                    return 
                end if
                !Determine new direction
                if(mod(iIteration,freq)==0) then!Every freq steps compute exact Hessian
                    call fdd(U,x,dim)
                    call My_dpotri_poQuery(U,dim,i)
                    if(i==0) then!Use exact Hessian if positive definite
                        call sycp(H,U,dim)
                        call syL2U(H,dim)
                        p=-matmul(H,fdnew)
                        phidnew=dot_product(fdnew,p)
                        a=1d0
                        cycle
                    end if
                end if
                !Otherwise use approximate Hessian
                s=x-s
                y=fdnew-y
                rho=1d0/dot_product(y,s)
                U=-rho*vector_direct_product(y,s,dim,dim)
                forall(i=1:dim)
                    U(i,i)=U(i,i)+1d0
                end forall
                H=matmul(transpose(U),matmul(H,U))+rho*vector_direct_product(s,s,dim,dim)
                p=-matmul(H,fdnew)
                phidnew=dot_product(fdnew,p)
                a=1d0
            end do
            if(iIteration>MaxQuasiNewtonIteration.and.QuasiNewtonWarning) then
                write(*,*)'Failed BFGS: max iteration exceeded!'
                write(*,*)'Euclidean norm of gradient =',Norm2(fdnew)
                QuasiNewtonWarning=.false.
            end if
        end subroutine BFGS
        !Do not compute exact Hessian at all
        subroutine BFGS_cheap(f,fd,x,dim)
            external::f,fd
            integer,intent(in)::dim
            real*8,dimension(dim),intent(inout)::x
            real*8,parameter::c1=1d-4,c2=0.9d0! 0 < c1 < c2 < 1: Wolfe constant
            integer::iIteration,i
            real*8::a,fnew,phidnew,rho
            real*8,dimension(dim)::p,fdnew,s,y
            real*8,dimension(dim,dim)::U,H!Approximate inverse Hessian
            !Initialize
                call f(fnew,x,dim)
                call fd(fdnew,x,dim)
                p=-fdnew
                phidnew=-dot_product(fdnew,fdnew)
                if(-phidnew<QuasiNewtonTol) return
                if(fnew==0d0) then
                    a=1d0
                else
                    a=-fnew/phidnew
                end if
                s=x
                y=fdnew
                !Initial approximate inverse Hessian = a
                call Wolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim)
                phidnew=dot_product(fdnew,fdnew)
                if(phidnew<QuasiNewtonTol) return
                if(dot_product(p,p)*a*a<QuasiNewtonTol) then
                    if(QuasiNewtonWarning) then
                        write(*,'(1x,A84)')'BFGS warning: step length has converged, but gradient norm has not met accuracy goal'
                        write(*,'(1x,A56)')'A best estimation rather than exact solution is returned'
                        write(*,*)'Euclidean norm of gradient =',Sqrt(phidnew)
                        QuasiNewtonWarning=.false.
                    end if
                    return 
                end if
                s=x-s
                y=fdnew-y
                rho=1d0/dot_product(y,s)
                U=-rho*vector_direct_product(y,s,dim,dim)
                forall(i=1:dim)
                    U(i,i)=U(i,i)+1d0
                end forall
                H=a*matmul(transpose(U),U)+rho*vector_direct_product(s,s,dim,dim)
                p=-matmul(H,fdnew)
                phidnew=dot_product(fdnew,p)
                a=1d0
            do iIteration=1,MaxQuasiNewtonIteration
                !Prepare
                s=x
                y=fdnew
                !Line search
                call Wolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim)
                phidnew=dot_product(fdnew,fdnew)
                if(phidnew<QuasiNewtonTol) return
                if(dot_product(p,p)*a*a<QuasiNewtonTol) then
                    if(QuasiNewtonWarning) then
                        write(*,'(1x,A84)')'BFGS warning: step length has converged, but gradient norm has not met accuracy goal'
                        write(*,'(1x,A56)')'A best estimation rather than exact solution is returned'
                        write(*,*)'Euclidean norm of gradient =',Sqrt(phidnew)
                        QuasiNewtonWarning=.false.
                    end if
                    return 
                end if
                !Determine new direction
                s=x-s
                y=fdnew-y
                rho=1d0/dot_product(y,s)
                U=-rho*vector_direct_product(y,s,dim,dim)
                forall(i=1:dim)
                    U(i,i)=U(i,i)+1d0
                end forall
                H=matmul(transpose(U),matmul(H,U))+rho*vector_direct_product(s,s,dim,dim)
                p=-matmul(H,fdnew)
                phidnew=dot_product(fdnew,p)
                a=1d0
            end do
            if(iIteration>MaxQuasiNewtonIteration.and.QuasiNewtonWarning) then
                write(*,*)'Failed BFGS: max iteration exceeded!'
                write(*,*)'Euclidean norm of gradient =',Norm2(fdnew)
                QuasiNewtonWarning=.false.
            end if
        end subroutine BFGS_cheap
        !Compute Hessian numerically through central difference by calling djacobi
        !djacobi requires fd_j has the form of: subroutine fd_j(dim,dim,x,f'(x))
        subroutine BFGS_NumericalHessian(f,fd,fd_j,x,dim,freq)
            external::f,fd,fd_j
            integer,intent(in)::dim,freq
            real*8,dimension(dim),intent(inout)::x
            real*8,parameter::c1=1d-4,c2=0.9d0! 0 < c1 < c2 < 1: Wolfe constant
            integer::iIteration,i
            real*8::a,fnew,phidnew,rho
            real*8,dimension(dim)::p,fdnew,s,y
            real*8,dimension(dim,dim)::U,H!Approximate inverse Hessian
            !Initialize
                call f(fnew,x,dim)
                call fd(fdnew,x,dim)
                i=djacobi(fd_j,dim,dim,H,x,jacobiPrecision)
                p=-fdnew
                call My_dpotri_poQuery(H,dim,i)
                if(i==0) then
                    call syL2U(H,dim)
                    p=-matmul(H,fdnew)
                    phidnew=dot_product(fdnew,p)
                    a=1d0
                else!Hessian is not positive definite, use a as initial approximate inverse Hessian
                    p=-fdnew
                    phidnew=-dot_product(fdnew,fdnew)
                    if(-phidnew<QuasiNewtonTol) return
                    if(fnew==0d0) then
                        a=1d0
                    else
                        a=-fnew/phidnew
                    end if
                    s=x
                    y=fdnew
                    call Wolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim)
                    phidnew=dot_product(fdnew,fdnew)
                    if(phidnew<QuasiNewtonTol) return
                    if(dot_product(p,p)*a*a<QuasiNewtonTol) then
                        if(QuasiNewtonWarning) then
                            write(*,'(1x,A84)')'BFGS warning: step length has converged, but gradient norm has not met accuracy goal'
                            write(*,'(1x,A56)')'A best estimation rather than exact solution is returned'
                            write(*,*)'Euclidean norm of gradient =',Sqrt(phidnew)
                            QuasiNewtonWarning=.false.
                        end if
                        return 
                    end if
                    s=x-s
                    y=fdnew-y
                    rho=1d0/dot_product(y,s)
                    U=-rho*vector_direct_product(y,s,dim,dim)
                    forall(i=1:dim)
                        U(i,i)=U(i,i)+1d0
                    end forall
                    H=a*U
                    H=matmul(transpose(U),H)+rho*vector_direct_product(s,s,dim,dim)
                    p=-matmul(H,fdnew)
                    phidnew=dot_product(fdnew,p)
                    a=1d0
                end if
            do iIteration=1,MaxQuasiNewtonIteration
                !Prepare
                s=x
                y=fdnew
                !Line search
                call Wolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim)
                phidnew=dot_product(fdnew,fdnew)
                if(phidnew<QuasiNewtonTol) return
                if(dot_product(p,p)*a*a<QuasiNewtonTol) then
                    if(QuasiNewtonWarning) then
                        write(*,'(1x,A84)')'BFGS warning: step length has converged, but gradient norm has not met accuracy goal'
                        write(*,'(1x,A56)')'A best estimation rather than exact solution is returned'
                        write(*,*)'Euclidean norm of gradient =',Sqrt(phidnew)
                        QuasiNewtonWarning=.false.
                    end if
                    return 
                end if
                !Determine new direction
                if(mod(iIteration,freq)==0) then!Every freq steps compute exact Hessian
                    i=djacobi(fd_j,dim,dim,U,x,jacobiPrecision)
                    call My_dpotri_poQuery(U,dim,i)
                    if(i==0) then!Use exact Hessian if positive definite
                        call sycp(H,U,dim)
                        call syL2U(H,dim)
                        p=-matmul(H,fdnew)
                        phidnew=dot_product(fdnew,p)
                        a=1d0
                        cycle
                    end if
                end if
                !Otherwise use approximate Hessian
                s=x-s
                y=fdnew-y
                rho=1d0/dot_product(y,s)
                U=-rho*vector_direct_product(y,s,dim,dim)
                forall(i=1:dim)
                    U(i,i)=U(i,i)+1d0
                end forall
                H=matmul(transpose(U),matmul(H,U))+rho*vector_direct_product(s,s,dim,dim)
                p=-matmul(H,fdnew)
                phidnew=dot_product(fdnew,p)
                a=1d0
            end do
            if(iIteration>MaxQuasiNewtonIteration.and.QuasiNewtonWarning) then
                write(*,*)'Failed BFGS: max iteration exceeded!'
                write(*,*)'Euclidean norm of gradient =',Norm2(fdnew)
                QuasiNewtonWarning=.false.
            end if
        end subroutine BFGS_NumericalHessian
        !Strong Wolfe condition usually performs better
        subroutine BFGS_Strong(f,fd,fdd,x,dim,freq)
            external::f,fd,fdd
            integer,intent(in)::dim,freq
            real*8,dimension(dim),intent(inout)::x
            real*8,parameter::c1=1d-4,c2=0.9d0! 0 < c1 < c2 < 1: Wolfe constant
            integer::iIteration,i
            real*8::a,fnew,phidnew,rho
            real*8,dimension(dim)::p,fdnew,s,y
            real*8,dimension(dim,dim)::U,H!Approximate inverse Hessian
            !Initialize
                call f(fnew,x,dim)
                call fd(fdnew,x,dim)
                call fdd(H,x,dim)
                p=-fdnew
                call My_dpotri_poQuery(H,dim,i)
                if(i==0) then
                    call syL2U(H,dim)
                    p=-matmul(H,fdnew)
                    phidnew=dot_product(fdnew,p)
                    a=1d0
                else!Hessian is not positive definite, use a as initial approximate inverse Hessian
                    p=-fdnew
                    phidnew=-dot_product(fdnew,fdnew)
                    if(-phidnew<QuasiNewtonTol) return
                    if(fnew==0d0) then
                        a=1d0
                    else
                        a=-fnew/phidnew
                    end if
                    s=x
                    y=fdnew
                    call StrongWolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim)
                    phidnew=dot_product(fdnew,fdnew)
                    if(phidnew<QuasiNewtonTol) return
                    if(dot_product(p,p)*a*a<QuasiNewtonTol) then
                        if(QuasiNewtonWarning) then
                            write(*,'(1x,A84)')'BFGS warning: step length has converged, but gradient norm has not met accuracy goal'
                            write(*,'(1x,A56)')'A best estimation rather than exact solution is returned'
                            write(*,*)'Euclidean norm of gradient =',Sqrt(phidnew)
                            QuasiNewtonWarning=.false.
                        end if
                        return 
                    end if
                    s=x-s
                    y=fdnew-y
                    rho=1d0/dot_product(y,s)
                    U=-rho*vector_direct_product(y,s,dim,dim)
                    forall(i=1:dim)
                        U(i,i)=U(i,i)+1d0
                    end forall
                    H=a*U
                    H=matmul(transpose(U),H)+rho*vector_direct_product(s,s,dim,dim)
                    p=-matmul(H,fdnew)
                    phidnew=dot_product(fdnew,p)
                    a=1d0
                end if
            do iIteration=1,MaxQuasiNewtonIteration
                !Prepare
                s=x
                y=fdnew
                !Line search
                call StrongWolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim)
                phidnew=dot_product(fdnew,fdnew)
                if(phidnew<QuasiNewtonTol) return
                if(dot_product(p,p)*a*a<QuasiNewtonTol) then
                    if(QuasiNewtonWarning) then
                        write(*,'(1x,A84)')'BFGS warning: step length has converged, but gradient norm has not met accuracy goal'
                        write(*,'(1x,A56)')'A best estimation rather than exact solution is returned'
                        write(*,*)'Euclidean norm of gradient =',Sqrt(phidnew)
                        QuasiNewtonWarning=.false.
                    end if
                    return 
                end if
                !Determine new direction
                if(mod(iIteration,freq)==0) then!Every freq steps compute exact Hessian
                    call fdd(U,x,dim)
                    call My_dpotri_poQuery(U,dim,i)
                    if(i==0) then!Use exact Hessian if positive definite
                        call sycp(H,U,dim)
                        call syL2U(H,dim)
                        p=-matmul(H,fdnew)
                        phidnew=dot_product(fdnew,p)
                        a=1d0
                        cycle
                    end if
                end if
                !Otherwise use approximate Hessian
                s=x-s
                y=fdnew-y
                rho=1d0/dot_product(y,s)
                U=-rho*vector_direct_product(y,s,dim,dim)
                forall(i=1:dim)
                    U(i,i)=U(i,i)+1d0
                end forall
                H=matmul(transpose(U),matmul(H,U))+rho*vector_direct_product(s,s,dim,dim)
                p=-matmul(H,fdnew)
                phidnew=dot_product(fdnew,p)
                a=1d0
            end do
            if(iIteration>MaxQuasiNewtonIteration.and.QuasiNewtonWarning) then
                write(*,*)'Failed BFGS: max iteration exceeded!'
                write(*,*)'Euclidean norm of gradient =',Norm2(fdnew)
                QuasiNewtonWarning=.false.
            end if
        end subroutine BFGS_Strong
        subroutine BFGS_Strong_cheap(f,fd,x,dim)
            external::f,fd
            integer,intent(in)::dim
            real*8,dimension(dim),intent(inout)::x
            real*8,parameter::c1=1d-4,c2=0.9d0! 0 < c1 < c2 < 1: Wolfe constant
            integer::iIteration,i
            real*8::a,fnew,phidnew,rho
            real*8,dimension(dim)::p,fdnew,s,y
            real*8,dimension(dim,dim)::U,H!Approximate inverse Hessian
            !Initialize
                call f(fnew,x,dim)
                call fd(fdnew,x,dim)
                p=-fdnew
                phidnew=-dot_product(fdnew,fdnew)
                if(-phidnew<QuasiNewtonTol) return
                if(fnew==0d0) then
                    a=1d0
                else
                    a=-fnew/phidnew
                end if
                s=x
                y=fdnew
                !Initial approximate inverse Hessian = a
                call StrongWolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim)
                phidnew=dot_product(fdnew,fdnew)
                if(phidnew<QuasiNewtonTol) return
                if(dot_product(p,p)*a*a<QuasiNewtonTol) then
                    if(QuasiNewtonWarning) then
                        write(*,'(1x,A84)')'BFGS warning: step length has converged, but gradient norm has not met accuracy goal'
                        write(*,'(1x,A56)')'A best estimation rather than exact solution is returned'
                        write(*,*)'Euclidean norm of gradient =',Sqrt(phidnew)
                        QuasiNewtonWarning=.false.
                    end if
                    return 
                end if
                s=x-s
                y=fdnew-y
                rho=1d0/dot_product(y,s)
                U=-rho*vector_direct_product(y,s,dim,dim)
                forall(i=1:dim)
                    U(i,i)=U(i,i)+1d0
                end forall
                H=a*matmul(transpose(U),U)+rho*vector_direct_product(s,s,dim,dim)
                p=-matmul(H,fdnew)
                phidnew=dot_product(fdnew,p)
                a=1d0
            do iIteration=1,MaxQuasiNewtonIteration
                !Prepare
                s=x
                y=fdnew
                !Line search
                call StrongWolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim)
                phidnew=dot_product(fdnew,fdnew)
                if(phidnew<QuasiNewtonTol) return
                if(dot_product(p,p)*a*a<QuasiNewtonTol) then
                    if(QuasiNewtonWarning) then
                        write(*,'(1x,A84)')'BFGS warning: step length has converged, but gradient norm has not met accuracy goal'
                        write(*,'(1x,A56)')'A best estimation rather than exact solution is returned'
                        write(*,*)'Euclidean norm of gradient =',Sqrt(phidnew)
                        QuasiNewtonWarning=.false.
                    end if
                    return 
                end if
                !Determine new direction
                s=x-s
                y=fdnew-y
                rho=1d0/dot_product(y,s)
                U=-rho*vector_direct_product(y,s,dim,dim)
                forall(i=1:dim)
                    U(i,i)=U(i,i)+1d0
                end forall
                H=matmul(transpose(U),matmul(H,U))+rho*vector_direct_product(s,s,dim,dim)
                p=-matmul(H,fdnew)
                phidnew=dot_product(fdnew,p)
                a=1d0
            end do
            if(iIteration>MaxQuasiNewtonIteration.and.QuasiNewtonWarning) then
                write(*,*)'Failed BFGS: max iteration exceeded!'
                write(*,*)'Euclidean norm of gradient =',Norm2(fdnew)
                QuasiNewtonWarning=.false.
            end if
        end subroutine BFGS_Strong_cheap
        subroutine BFGS_Strong_NumericalHessian(f,fd,fd_j,x,dim,freq)
            external::f,fd,fd_j
            integer,intent(in)::dim,freq
            real*8,dimension(dim),intent(inout)::x
            real*8,parameter::c1=1d-4,c2=0.9d0! 0 < c1 < c2 < 1: Wolfe constant
            integer::iIteration,i
            real*8::a,fnew,phidnew,rho
            real*8,dimension(dim)::p,fdnew,s,y
            real*8,dimension(dim,dim)::U,H!Approximate inverse Hessian
            !Initialize
                call f(fnew,x,dim)
                call fd(fdnew,x,dim)
                i=djacobi(fd_j,dim,dim,H,x,jacobiPrecision)
                p=-fdnew
                call My_dpotri_poQuery(H,dim,i)
                if(i==0) then
                    call syL2U(H,dim)
                    p=-matmul(H,fdnew)
                    phidnew=dot_product(fdnew,p)
                    a=1d0
                else!Hessian is not positive definite, use a as initial approximate inverse Hessian
                    p=-fdnew
                    phidnew=-dot_product(fdnew,fdnew)
                    if(-phidnew<QuasiNewtonTol) return
                    if(fnew==0d0) then
                        a=1d0
                    else
                        a=-fnew/phidnew
                    end if
                    s=x
                    y=fdnew
                    call StrongWolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim)
                    phidnew=dot_product(fdnew,fdnew)
                    if(phidnew<QuasiNewtonTol) return
                    if(dot_product(p,p)*a*a<QuasiNewtonTol) then
                        if(QuasiNewtonWarning) then
                            write(*,'(1x,A84)')'BFGS warning: step length has converged, but gradient norm has not met accuracy goal'
                            write(*,'(1x,A56)')'A best estimation rather than exact solution is returned'
                            write(*,*)'Euclidean norm of gradient =',Sqrt(phidnew)
                            QuasiNewtonWarning=.false.
                        end if
                        return 
                    end if
                    s=x-s
                    y=fdnew-y
                    rho=1d0/dot_product(y,s)
                    U=-rho*vector_direct_product(y,s,dim,dim)
                    forall(i=1:dim)
                        U(i,i)=U(i,i)+1d0
                    end forall
                    H=a*U
                    H=matmul(transpose(U),H)+rho*vector_direct_product(s,s,dim,dim)
                    p=-matmul(H,fdnew)
                    phidnew=dot_product(fdnew,p)
                    a=1d0
                end if
            do iIteration=1,MaxQuasiNewtonIteration
                !Prepare
                s=x
                y=fdnew
                !Line search
                call StrongWolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim)
                phidnew=dot_product(fdnew,fdnew)
                if(phidnew<QuasiNewtonTol) return
                if(dot_product(p,p)*a*a<QuasiNewtonTol) then
                    if(QuasiNewtonWarning) then
                        write(*,'(1x,A84)')'BFGS warning: step length has converged, but gradient norm has not met accuracy goal'
                        write(*,'(1x,A56)')'A best estimation rather than exact solution is returned'
                        write(*,*)'Euclidean norm of gradient =',Sqrt(phidnew)
                        QuasiNewtonWarning=.false.
                    end if
                    return 
                end if
                !Determine new direction
                if(mod(iIteration,freq)==0) then!Every freq steps compute exact Hessian
                    i=djacobi(fd_j,dim,dim,U,x,jacobiPrecision)
                    call My_dpotri_poQuery(U,dim,i)
                    if(i==0) then!Use exact Hessian if positive definite
                        call sycp(H,U,dim)
                        call syL2U(H,dim)
                        p=-matmul(H,fdnew)
                        phidnew=dot_product(fdnew,p)
                        a=1d0
                        cycle
                    end if
                end if
                !Otherwise use approximate Hessian
                s=x-s
                y=fdnew-y
                rho=1d0/dot_product(y,s)
                U=-rho*vector_direct_product(y,s,dim,dim)
                forall(i=1:dim)
                    U(i,i)=U(i,i)+1d0
                end forall
                H=matmul(transpose(U),matmul(H,U))+rho*vector_direct_product(s,s,dim,dim)
                p=-matmul(H,fdnew)
                phidnew=dot_product(fdnew,p)
                a=1d0
            end do
            if(iIteration>MaxQuasiNewtonIteration.and.QuasiNewtonWarning) then
                write(*,*)'Failed BFGS: max iteration exceeded!'
                write(*,*)'Euclidean norm of gradient =',Norm2(fdnew)
                QuasiNewtonWarning=.false.
            end if
        end subroutine BFGS_Strong_NumericalHessian
        !When it is cheap to evaluate f' along with f
        subroutine BFGS_Strong_fdwithf(f,fd,f_fd,fdd,x,dim,freq)
            external::f,fd,f_fd,fdd
            integer,intent(in)::dim,freq
            real*8,dimension(dim),intent(inout)::x
            real*8,parameter::c1=1d-4,c2=0.9d0! 0 < c1 < c2 < 1: Wolfe constant
            integer::iIteration,i
            real*8::a,fnew,phidnew,rho
            real*8,dimension(dim)::p,fdnew,s,y
            real*8,dimension(dim,dim)::U,H!Approximate inverse Hessian
            !Initialize
                call f_fd(fnew,fdnew,x,dim)
                call fdd(H,x,dim)
                p=-fdnew
                call My_dpotri_poQuery(H,dim,i)
                if(i==0) then
                    call syL2U(H,dim)
                    p=-matmul(H,fdnew)
                    phidnew=dot_product(fdnew,p)
                    a=1d0
                else!Hessian is not positive definite, use a as initial approximate inverse Hessian
                    p=-fdnew
                    phidnew=-dot_product(fdnew,fdnew)
                    if(-phidnew<QuasiNewtonTol) return
                    if(fnew==0d0) then
                        a=1d0
                    else
                        a=-fnew/phidnew
                    end if
                    s=x
                    y=fdnew
                    call StrongWolfe_fdwithf(c1,c2,f,fd,f_fd,x,a,p,fnew,phidnew,fdnew,dim)
                    phidnew=dot_product(fdnew,fdnew)
                    if(phidnew<QuasiNewtonTol) return
                    if(dot_product(p,p)*a*a<QuasiNewtonTol) then
                        if(QuasiNewtonWarning) then
                            write(*,'(1x,A84)')'BFGS warning: step length has converged, but gradient norm has not met accuracy goal'
                            write(*,'(1x,A56)')'A best estimation rather than exact solution is returned'
                            write(*,*)'Euclidean norm of gradient =',Sqrt(phidnew)
                            QuasiNewtonWarning=.false.
                        end if
                        return 
                    end if
                    s=x-s
                    y=fdnew-y
                    rho=1d0/dot_product(y,s)
                    U=-rho*vector_direct_product(y,s,dim,dim)
                    forall(i=1:dim)
                        U(i,i)=U(i,i)+1d0
                    end forall
                    H=a*U
                    H=matmul(transpose(U),H)+rho*vector_direct_product(s,s,dim,dim)
                    p=-matmul(H,fdnew)
                    phidnew=dot_product(fdnew,p)
                    a=1d0
                end if
            do iIteration=1,MaxQuasiNewtonIteration
                !Prepare
                s=x
                y=fdnew
                !Line search
                call StrongWolfe_fdwithf(c1,c2,f,fd,f_fd,x,a,p,fnew,phidnew,fdnew,dim)
                phidnew=dot_product(fdnew,fdnew)
                if(phidnew<QuasiNewtonTol) return
                if(dot_product(p,p)*a*a<QuasiNewtonTol) then
                    if(QuasiNewtonWarning) then
                        write(*,'(1x,A84)')'BFGS warning: step length has converged, but gradient norm has not met accuracy goal'
                        write(*,'(1x,A56)')'A best estimation rather than exact solution is returned'
                        write(*,*)'Euclidean norm of gradient =',Sqrt(phidnew)
                        QuasiNewtonWarning=.false.
                    end if
                    return 
                end if
                !Determine new direction
                if(mod(iIteration,freq)==0) then!Every freq steps compute exact Hessian
                    call fdd(U,x,dim)
                    call My_dpotri_poQuery(U,dim,i)
                    if(i==0) then!Use exact Hessian if positive definite
                        call sycp(H,U,dim)
                        call syL2U(H,dim)
                        p=-matmul(H,fdnew)
                        phidnew=dot_product(fdnew,p)
                        a=1d0
                        cycle
                    end if
                end if
                !Otherwise use approximate Hessian
                s=x-s
                y=fdnew-y
                rho=1d0/dot_product(y,s)
                U=-rho*vector_direct_product(y,s,dim,dim)
                forall(i=1:dim)
                    U(i,i)=U(i,i)+1d0
                end forall
                H=matmul(transpose(U),matmul(H,U))+rho*vector_direct_product(s,s,dim,dim)
                p=-matmul(H,fdnew)
                phidnew=dot_product(fdnew,p)
                a=1d0
            end do
            if(iIteration>MaxQuasiNewtonIteration.and.QuasiNewtonWarning) then
                write(*,*)'Failed BFGS: max iteration exceeded!'
                write(*,*)'Euclidean norm of gradient =',Norm2(fdnew)
                QuasiNewtonWarning=.false.
            end if
        end subroutine BFGS_Strong_fdwithf
        subroutine BFGS_Strong_cheap_fdwithf(f,fd,f_fd,x,dim)
            external::f,fd,f_fd
            integer,intent(in)::dim
            real*8,dimension(dim),intent(inout)::x
            real*8,parameter::c1=1d-4,c2=0.9d0! 0 < c1 < c2 < 1: Wolfe constant
            integer::iIteration,i
            real*8::a,fnew,phidnew,rho
            real*8,dimension(dim)::p,fdnew,s,y
            real*8,dimension(dim,dim)::U,H!Approximate inverse Hessian
            !Initialize
                call f_fd(fnew,fdnew,x,dim)
                p=-fdnew
                phidnew=-dot_product(fdnew,fdnew)
                if(-phidnew<QuasiNewtonTol) return
                if(fnew==0d0) then
                    a=1d0
                else
                    a=-fnew/phidnew
                end if
                s=x
                y=fdnew
                !Initial approximate inverse Hessian = a
                call StrongWolfe_fdwithf(c1,c2,f,fd,f_fd,x,a,p,fnew,phidnew,fdnew,dim)
                phidnew=dot_product(fdnew,fdnew)
                if(phidnew<QuasiNewtonTol) return
                if(dot_product(p,p)*a*a<QuasiNewtonTol) then
                    if(QuasiNewtonWarning) then
                        write(*,'(1x,A84)')'BFGS warning: step length has converged, but gradient norm has not met accuracy goal'
                        write(*,'(1x,A56)')'A best estimation rather than exact solution is returned'
                        write(*,*)'Euclidean norm of gradient =',Sqrt(phidnew)
                        QuasiNewtonWarning=.false.
                    end if
                    return 
                end if
                s=x-s
                y=fdnew-y
                rho=1d0/dot_product(y,s)
                U=-rho*vector_direct_product(y,s,dim,dim)
                forall(i=1:dim)
                    U(i,i)=U(i,i)+1d0
                end forall
                H=a*matmul(transpose(U),U)+rho*vector_direct_product(s,s,dim,dim)
                p=-matmul(H,fdnew)
                phidnew=dot_product(fdnew,p)
                a=1d0
            do iIteration=1,MaxQuasiNewtonIteration
                !Prepare
                s=x
                y=fdnew
                !Line search
                call StrongWolfe_fdwithf(c1,c2,f,fd,f_fd,x,a,p,fnew,phidnew,fdnew,dim)
                phidnew=dot_product(fdnew,fdnew)
                if(phidnew<QuasiNewtonTol) return
                if(dot_product(p,p)*a*a<QuasiNewtonTol) then
                    if(QuasiNewtonWarning) then
                        write(*,'(1x,A84)')'BFGS warning: step length has converged, but gradient norm has not met accuracy goal'
                        write(*,'(1x,A56)')'A best estimation rather than exact solution is returned'
                        write(*,*)'Euclidean norm of gradient =',Sqrt(phidnew)
                        QuasiNewtonWarning=.false.
                    end if
                    return 
                end if
                !Determine new direction
                s=x-s
                y=fdnew-y
                rho=1d0/dot_product(y,s)
                U=-rho*vector_direct_product(y,s,dim,dim)
                forall(i=1:dim)
                    U(i,i)=U(i,i)+1d0
                end forall
                H=matmul(transpose(U),matmul(H,U))+rho*vector_direct_product(s,s,dim,dim)
                p=-matmul(H,fdnew)
                phidnew=dot_product(fdnew,p)
                a=1d0
            end do
            if(iIteration>MaxQuasiNewtonIteration.and.QuasiNewtonWarning) then
                write(*,*)'Failed BFGS: max iteration exceeded!'
                write(*,*)'Euclidean norm of gradient =',Norm2(fdnew)
                QuasiNewtonWarning=.false.
            end if
        end subroutine BFGS_Strong_cheap_fdwithf
        subroutine BFGS_Strong_NumericalHessian_fdwithf(f,fd,f_fd,fd_j,x,dim,freq)
            external::f,fd,f_fd,fd_j
            integer,intent(in)::dim,freq
            real*8,dimension(dim),intent(inout)::x
            real*8,parameter::c1=1d-4,c2=0.9d0! 0 < c1 < c2 < 1: Wolfe constant
            integer::iIteration,i
            real*8::a,fnew,phidnew,rho
            real*8,dimension(dim)::p,fdnew,s,y
            real*8,dimension(dim,dim)::U,H!Approximate inverse Hessian
            !Initialize
                call f_fd(fnew,fdnew,x,dim)
                i=djacobi(fd_j,dim,dim,H,x,jacobiPrecision)
                p=-fdnew
                call My_dpotri_poQuery(H,dim,i)
                if(i==0) then
                    call syL2U(H,dim)
                    p=-matmul(H,fdnew)
                    phidnew=dot_product(fdnew,p)
                    a=1d0
                else!Hessian is not positive definite, use a as initial approximate inverse Hessian
                    p=-fdnew
                    phidnew=-dot_product(fdnew,fdnew)
                    if(-phidnew<QuasiNewtonTol) return
                    if(fnew==0d0) then
                        a=1d0
                    else
                        a=-fnew/phidnew
                    end if
                    s=x
                    y=fdnew
                    call StrongWolfe_fdwithf(c1,c2,f,fd,f_fd,x,a,p,fnew,phidnew,fdnew,dim)
                    phidnew=dot_product(fdnew,fdnew)
                    if(phidnew<QuasiNewtonTol) return
                    if(dot_product(p,p)*a*a<QuasiNewtonTol) then
                        if(QuasiNewtonWarning) then
                            write(*,'(1x,A84)')'BFGS warning: step length has converged, but gradient norm has not met accuracy goal'
                            write(*,'(1x,A56)')'A best estimation rather than exact solution is returned'
                            write(*,*)'Euclidean norm of gradient =',Sqrt(phidnew)
                            QuasiNewtonWarning=.false.
                        end if
                        return 
                    end if
                    s=x-s
                    y=fdnew-y
                    rho=1d0/dot_product(y,s)
                    U=-rho*vector_direct_product(y,s,dim,dim)
                    forall(i=1:dim)
                        U(i,i)=U(i,i)+1d0
                    end forall
                    H=a*U
                    H=matmul(transpose(U),H)+rho*vector_direct_product(s,s,dim,dim)
                    p=-matmul(H,fdnew)
                    phidnew=dot_product(fdnew,p)
                    a=1d0
                end if
            do iIteration=1,MaxQuasiNewtonIteration
                !Prepare
                s=x
                y=fdnew
                !Line search
                call StrongWolfe_fdwithf(c1,c2,f,fd,f_fd,x,a,p,fnew,phidnew,fdnew,dim)
                phidnew=dot_product(fdnew,fdnew)
                if(phidnew<QuasiNewtonTol) return
                if(dot_product(p,p)*a*a<QuasiNewtonTol) then
                    if(QuasiNewtonWarning) then
                        write(*,'(1x,A84)')'BFGS warning: step length has converged, but gradient norm has not met accuracy goal'
                        write(*,'(1x,A56)')'A best estimation rather than exact solution is returned'
                        write(*,*)'Euclidean norm of gradient =',Sqrt(phidnew)
                        QuasiNewtonWarning=.false.
                    end if
                    return 
                end if
                !Determine new direction
                if(mod(iIteration,freq)==0) then!Every freq steps compute exact Hessian
                    i=djacobi(fd_j,dim,dim,U,x,jacobiPrecision)
                    call My_dpotri_poQuery(U,dim,i)
                    if(i==0) then!Use exact Hessian if positive definite
                        call sycp(H,U,dim)
                        call syL2U(H,dim)
                        p=-matmul(H,fdnew)
                        phidnew=dot_product(fdnew,p)
                        a=1d0
                        cycle
                    end if
                end if
                !Otherwise use approximate Hessian
                s=x-s
                y=fdnew-y
                rho=1d0/dot_product(y,s)
                U=-rho*vector_direct_product(y,s,dim,dim)
                forall(i=1:dim)
                    U(i,i)=U(i,i)+1d0
                end forall
                H=matmul(transpose(U),matmul(H,U))+rho*vector_direct_product(s,s,dim,dim)
                p=-matmul(H,fdnew)
                phidnew=dot_product(fdnew,p)
                a=1d0
            end do
            if(iIteration>MaxQuasiNewtonIteration.and.QuasiNewtonWarning) then
                write(*,*)'Failed BFGS: max iteration exceeded!'
                write(*,*)'Euclidean norm of gradient =',Norm2(fdnew)
                QuasiNewtonWarning=.false.
            end if
        end subroutine BFGS_Strong_NumericalHessian_fdwithf

        !Memory usage = O( M * dim ). 2 < M < 31 is recommended
        !On input x is an initial guess of the minimum point of f(x), on exit x is the minimum point
        !Limited-memory Broyden–Fletcher–Goldfarb–Shanno (L-BFGS) quasi-Newton method, requiring Wolfe condition
        subroutine LBFGS(f,fd,x,dim,M)
            external::f,fd
            integer,intent(in)::dim,M
            real*8,dimension(dim),intent(inout)::x
            real*8,parameter::c1=1d-4,c2=0.9d0! 0 < c1 < c2 < 1: Wolfe constant
            integer::iIteration,i,recent
            real*8::a,fnew,phidnew
            real*8,dimension(dim)::p,fdnew,xold,fdold
            real*8,dimension(0:M)::rho,alpha
            real*8,dimension(dim,0:M)::s,y
            !Initialize
                call f(fnew,x,dim)
                call fd(fdnew,x,dim)
                p=-fdnew
                phidnew=-dot_product(fdnew,fdnew)
                if(-phidnew<QuasiNewtonTol) return
                if(fnew==0d0) then
                    a=1d0
                else
                    a=-fnew/phidnew
                end if
                xold=x
                fdold=fdnew
                !Initial approximate inverse Hessian = a
                call Wolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim)
                phidnew=dot_product(fdnew,fdnew)
                if(phidnew<QuasiNewtonTol) return
                if(dot_product(p,p)*a*a<QuasiNewtonTol) then
                    if(QuasiNewtonWarning) then
                        write(*,'(1x,A84)')'BFGS warning: step length has converged, but gradient norm has not met accuracy goal'
                        write(*,'(1x,A56)')'A best estimation rather than exact solution is returned'
                        write(*,*)'Euclidean norm of gradient =',Sqrt(phidnew)
                        QuasiNewtonWarning=.false.
                    end if
                    return 
                end if
                recent=0
                s(:,0)=x-xold
                y(:,0)=fdnew-fdold
                rho(0)=1d0/dot_product(y(:,0),s(:,0))
                do iIteration=1,M-1
                    xold=x
                    fdold=fdnew
                    !Determine new direction
                    p=fdnew
                    do i=recent,0,-1
                        alpha(i)=rho(i)*dot_product(s(:,i),p)
                        p=p-alpha(i)*y(:,i)
                    end do
                    p=p/rho(recent)/dot_product(y(:,recent),y(:,recent))
                    do i=0,recent
                        phidnew=rho(i)*dot_product(y(:,i),p)
                        p=p+(alpha(i)-phidnew)*s(:,i)
                    end do
                    p=-p
                    phidnew=dot_product(fdnew,p)
                    a=1d0
                    !Line search
                    call Wolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim)
                    recent=recent+1
                    s(:,recent)=x-xold
                    y(:,recent)=fdnew-fdold
                    rho(recent)=1d0/dot_product(y(:,recent),s(:,recent))
                end do
            do iIteration=1,MaxQuasiNewtonIteration
                xold=x
                fdold=fdnew
                !Determine new direction
                p=fdnew
                do i=recent,0,-1
                    alpha(i)=rho(i)*dot_product(s(:,i),p)
                    p=p-alpha(i)*y(:,i)
                end do
                do i=M-1,recent+1,-1
                    alpha(i)=rho(i)*dot_product(s(:,i),p)
                    p=p-alpha(i)*y(:,i)
                end do
                p=p/rho(recent)/dot_product(y(:,recent),y(:,recent))
                do i=recent+1,M-1
                    phidnew=rho(i)*dot_product(y(:,i),p)
                    p=p+(alpha(i)-phidnew)*s(:,i)
                end do
                do i=0,recent
                    phidnew=rho(i)*dot_product(y(:,i),p)
                    p=p+(alpha(i)-phidnew)*s(:,i)
                end do
                p=-p
                phidnew=dot_product(fdnew,p)
                a=1d0
                !Line search
                call Wolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim)
                phidnew=dot_product(fdnew,fdnew)
                if(phidnew<QuasiNewtonTol) return
                if(dot_product(p,p)*a*a<QuasiNewtonTol) then
                    if(QuasiNewtonWarning) then
                        write(*,'(1x,A86)')'L-BFGS warning: step length has converged, but gradient norm has not met accuracy goal'
                        write(*,'(1x,A56)')'A best estimation rather than exact solution is returned'
                        write(*,*)'Euclidean norm of gradient =',Sqrt(phidnew)
                        QuasiNewtonWarning=.false.
                    end if
                    return 
                end if
                recent=mod(recent+1,M)
                s(:,recent)=x-xold
                y(:,recent)=fdnew-fdold
                rho(recent)=1d0/dot_product(y(:,recent),s(:,recent))
            end do
            if(iIteration>MaxQuasiNewtonIteration.and.QuasiNewtonWarning) then
                write(*,*)'Failed L-BFGS: max iteration exceeded!'
                write(*,*)'Euclidean norm of gradient =',Norm2(fdnew)
                QuasiNewtonWarning=.false.
            end if
        end subroutine LBFGS
        !Strong Wolfe condition usually performs better
        subroutine LBFGS_Strong(f,fd,x,dim,M)
            external::f,fd
            integer,intent(in)::dim,M
            real*8,dimension(dim),intent(inout)::x
            real*8,parameter::c1=1d-4,c2=0.9d0! 0 < c1 < c2 < 1: Wolfe constant
            integer::iIteration,i,recent
            real*8::a,fnew,phidnew
            real*8,dimension(dim)::p,fdnew,xold,fdold
            real*8,dimension(0:M)::rho,alpha
            real*8,dimension(dim,0:M)::s,y
            !Initialize
                call f(fnew,x,dim)
                call fd(fdnew,x,dim)
                p=-fdnew
                phidnew=-dot_product(fdnew,fdnew)
                if(-phidnew<QuasiNewtonTol) return
                if(fnew==0d0) then
                    a=1d0
                else
                    a=-fnew/phidnew
                end if
                xold=x
                fdold=fdnew
                !Initial approximate inverse Hessian = a
                call StrongWolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim)
                phidnew=dot_product(fdnew,fdnew)
                if(phidnew<QuasiNewtonTol) return
                if(dot_product(p,p)*a*a<QuasiNewtonTol) then
                    if(QuasiNewtonWarning) then
                        write(*,'(1x,A84)')'BFGS warning: step length has converged, but gradient norm has not met accuracy goal'
                        write(*,'(1x,A56)')'A best estimation rather than exact solution is returned'
                        write(*,*)'Euclidean norm of gradient =',Sqrt(phidnew)
                        QuasiNewtonWarning=.false.
                    end if
                    return 
                end if
                recent=0
                s(:,0)=x-xold
                y(:,0)=fdnew-fdold
                rho(0)=1d0/dot_product(y(:,0),s(:,0))
                do iIteration=1,M-1
                    xold=x
                    fdold=fdnew
                    !Determine new direction
                    p=fdnew
                    do i=recent,0,-1
                        alpha(i)=rho(i)*dot_product(s(:,i),p)
                        p=p-alpha(i)*y(:,i)
                    end do
                    p=p/rho(recent)/dot_product(y(:,recent),y(:,recent))
                    do i=0,recent
                        phidnew=rho(i)*dot_product(y(:,i),p)
                        p=p+(alpha(i)-phidnew)*s(:,i)
                    end do
                    p=-p
                    phidnew=dot_product(fdnew,p)
                    a=1d0
                    !Line search
                    call StrongWolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim)
                    recent=recent+1
                    s(:,recent)=x-xold
                    y(:,recent)=fdnew-fdold
                    rho(recent)=1d0/dot_product(y(:,recent),s(:,recent))
                end do
            do iIteration=1,MaxQuasiNewtonIteration
                xold=x
                fdold=fdnew
                !Determine new direction
                p=fdnew
                do i=recent,0,-1
                    alpha(i)=rho(i)*dot_product(s(:,i),p)
                    p=p-alpha(i)*y(:,i)
                end do
                do i=M-1,recent+1,-1
                    alpha(i)=rho(i)*dot_product(s(:,i),p)
                    p=p-alpha(i)*y(:,i)
                end do
                p=p/rho(recent)/dot_product(y(:,recent),y(:,recent))
                do i=recent+1,M-1
                    phidnew=rho(i)*dot_product(y(:,i),p)
                    p=p+(alpha(i)-phidnew)*s(:,i)
                end do
                do i=0,recent
                    phidnew=rho(i)*dot_product(y(:,i),p)
                    p=p+(alpha(i)-phidnew)*s(:,i)
                end do
                p=-p
                phidnew=dot_product(fdnew,p)
                a=1d0
                !Line search
                call StrongWolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim)
                phidnew=dot_product(fdnew,fdnew)
                if(phidnew<QuasiNewtonTol) return
                if(dot_product(p,p)*a*a<QuasiNewtonTol) then
                    if(QuasiNewtonWarning) then
                        write(*,'(1x,A86)')'L-BFGS warning: step length has converged, but gradient norm has not met accuracy goal'
                        write(*,'(1x,A56)')'A best estimation rather than exact solution is returned'
                        write(*,*)'Euclidean norm of gradient =',Sqrt(phidnew)
                        QuasiNewtonWarning=.false.
                    end if
                    return 
                end if
                recent=mod(recent+1,M)
                s(:,recent)=x-xold
                y(:,recent)=fdnew-fdold
                rho(recent)=1d0/dot_product(y(:,recent),s(:,recent))
            end do
            if(iIteration>MaxQuasiNewtonIteration.and.QuasiNewtonWarning) then
                write(*,*)'Failed L-BFGS: max iteration exceeded!'
                write(*,*)'Euclidean norm of gradient =',Norm2(fdnew)
                QuasiNewtonWarning=.false.
            end if
        end subroutine LBFGS_Strong
        !When it is cheap to evaluate f' along with f
        subroutine LBFGS_Strong_fdwithf(f,fd,f_fd,x,dim,M)
            external::f,fd,f_fd
            integer,intent(in)::dim,M
            real*8,dimension(dim),intent(inout)::x
            real*8,parameter::c1=1d-4,c2=0.9d0! 0 < c1 < c2 < 1: Wolfe constant
            integer::iIteration,i,recent
            real*8::a,fnew,phidnew
            real*8,dimension(dim)::p,fdnew,xold,fdold
            real*8,dimension(0:M)::rho,alpha
            real*8,dimension(dim,0:M)::s,y
            !Initialize
                call f(fnew,x,dim)
                call fd(fdnew,x,dim)
                p=-fdnew
                phidnew=-dot_product(fdnew,fdnew)
                if(-phidnew<QuasiNewtonTol) return
                if(fnew==0d0) then
                    a=1d0
                else
                    a=-fnew/phidnew
                end if
                xold=x
                fdold=fdnew
                !Initial approximate inverse Hessian = a
                call StrongWolfe_fdwithf(c1,c2,f,fd,f_fd,x,a,p,fnew,phidnew,fdnew,dim)
                phidnew=dot_product(fdnew,fdnew)
                if(phidnew<QuasiNewtonTol) return
                if(dot_product(p,p)*a*a<QuasiNewtonTol) then
                    if(QuasiNewtonWarning) then
                        write(*,'(1x,A84)')'BFGS warning: step length has converged, but gradient norm has not met accuracy goal'
                        write(*,'(1x,A56)')'A best estimation rather than exact solution is returned'
                        write(*,*)'Euclidean norm of gradient =',Sqrt(phidnew)
                        QuasiNewtonWarning=.false.
                    end if
                    return 
                end if
                recent=0
                s(:,0)=x-xold
                y(:,0)=fdnew-fdold
                rho(0)=1d0/dot_product(y(:,0),s(:,0))
                do iIteration=1,M-1
                    xold=x
                    fdold=fdnew
                    !Determine new direction
                    p=fdnew
                    do i=recent,0,-1
                        alpha(i)=rho(i)*dot_product(s(:,i),p)
                        p=p-alpha(i)*y(:,i)
                    end do
                    p=p/rho(recent)/dot_product(y(:,recent),y(:,recent))
                    do i=0,recent
                        phidnew=rho(i)*dot_product(y(:,i),p)
                        p=p+(alpha(i)-phidnew)*s(:,i)
                    end do
                    p=-p
                    phidnew=dot_product(fdnew,p)
                    a=1d0
                    !Line search
                    call StrongWolfe_fdwithf(c1,c2,f,fd,f_fd,x,a,p,fnew,phidnew,fdnew,dim)
                    recent=recent+1
                    s(:,recent)=x-xold
                    y(:,recent)=fdnew-fdold
                    rho(recent)=1d0/dot_product(y(:,recent),s(:,recent))
                end do
            do iIteration=1,MaxQuasiNewtonIteration
                xold=x
                fdold=fdnew
                !Determine new direction
                p=fdnew
                do i=recent,0,-1
                    alpha(i)=rho(i)*dot_product(s(:,i),p)
                    p=p-alpha(i)*y(:,i)
                end do
                do i=M-1,recent+1,-1
                    alpha(i)=rho(i)*dot_product(s(:,i),p)
                    p=p-alpha(i)*y(:,i)
                end do
                p=p/rho(recent)/dot_product(y(:,recent),y(:,recent))
                do i=recent+1,M-1
                    phidnew=rho(i)*dot_product(y(:,i),p)
                    p=p+(alpha(i)-phidnew)*s(:,i)
                end do
                do i=0,recent
                    phidnew=rho(i)*dot_product(y(:,i),p)
                    p=p+(alpha(i)-phidnew)*s(:,i)
                end do
                p=-p
                phidnew=dot_product(fdnew,p)
                a=1d0
                !Line search
                call StrongWolfe_fdwithf(c1,c2,f,fd,f_fd,x,a,p,fnew,phidnew,fdnew,dim)
                phidnew=dot_product(fdnew,fdnew)
                if(phidnew<QuasiNewtonTol) return
                if(dot_product(p,p)*a*a<QuasiNewtonTol) then
                    if(QuasiNewtonWarning) then
                        write(*,'(1x,A86)')'L-BFGS warning: step length has converged, but gradient norm has not met accuracy goal'
                        write(*,'(1x,A56)')'A best estimation rather than exact solution is returned'
                        write(*,*)'Euclidean norm of gradient =',Sqrt(phidnew)
                        QuasiNewtonWarning=.false.
                    end if
                    return 
                end if
                recent=mod(recent+1,M)
                s(:,recent)=x-xold
                y(:,recent)=fdnew-fdold
                rho(recent)=1d0/dot_product(y(:,recent),s(:,recent))
            end do
            if(iIteration>MaxQuasiNewtonIteration.and.QuasiNewtonWarning) then
                write(*,*)'Failed L-BFGS: max iteration exceeded!'
                write(*,*)'Euclidean norm of gradient =',Norm2(fdnew)
                QuasiNewtonWarning=.false.
            end if
        end subroutine LBFGS_Strong_fdwithf

        !On input x is an initial guess of the minimum point of f(x), on exit x is the minimum point
        !Dai-Yun conjugate gradient method, requiring Wolfe condition 
        subroutine DYConjugateGradient(f,fd,x,dim)
            external::f,fd
            integer,intent(in)::dim
            real*8,dimension(dim),intent(inout)::x
            real*8,parameter::c1=1d-4,c2=0.45d0! 0 < c1 < c2 < 0.5: Wolfe constant
            integer::iIteration
            real*8::a,fnew,fold,phidnew,phidold
            real*8,dimension(dim)::p,fdnew,fdold
            !Initialize
                call f(fnew,x,dim)
                call fd(fdnew,x,dim)
                p=-fdnew
                phidnew=-dot_product(fdnew,fdnew)
                if(-phidnew<ConjugateGradientTol) return
                if(fnew==0d0) then
                    a=1d0
                else
                    a=-fnew/phidnew
                end if
            do iIteration=1,MaxConjugateGradientIteration
                !Prepare
                fold=fnew
                fdold=fdnew
                phidold=phidnew
                !Line search
                call Wolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim)
                phidnew=dot_product(fdnew,fdnew)
                if(phidnew<ConjugateGradientTol) return
                if(dot_product(p,p)*a*a<ConjugateGradientTol) then
                    if(ConjugateGradientWarning) then
                        write(*,'(1x,A98)')'Conjugate gradient warning: step length has converged, but gradient norm has not met accuracy goal'
                        write(*,'(1x,A56)')'A best estimation rather than exact solution is returned'
                        write(*,*)'Euclidean norm of gradient =',Sqrt(phidnew)
                        ConjugateGradientWarning=.false.
                    end if
                    return 
                end if
                !Determine new direction
                p=-fdnew+dot_product(fdnew,fdnew)/dot_product(fdnew-fdold,p)*p
                phidnew=dot_product(fdnew,p)
                if(phidnew>0d0) then!Ascent direction, reset to steepest descent direction
                    p=-fdnew
                    phidnew=-dot_product(fdnew,fdnew)
                end if
                a=a*phidold/phidnew
            end do
            if(iIteration>MaxConjugateGradientIteration.and.ConjugateGradientWarning) then
                write(*,*)'Failed conjugate gradient: max iteration exceeded!'
                write(*,*)'Euclidean norm of gradient =',Norm2(fdnew)
                ConjugateGradientWarning=.false.
            end if
        end subroutine DYConjugateGradient
        !To meet Nocedal's performance suggestion, Dai-Yun requires strong Wolfe condition
        subroutine DYConjugateGradient_Strong(f,fd,x,dim)
            external::f,fd
            integer,intent(in)::dim
            real*8,dimension(dim),intent(inout)::x
            real*8,parameter::c1=1d-4,c2=0.45d0! 0 < c1 < c2 < 0.5: Wolfe constant
            integer::iIteration
            real*8::a,fnew,fold,phidnew,phidold
            real*8,dimension(dim)::p,fdnew,fdold
            !Initialize
                call f(fnew,x,dim)
                call fd(fdnew,x,dim)
                p=-fdnew
                phidnew=-dot_product(fdnew,fdnew)
                if(-phidnew<ConjugateGradientTol) return
                if(fnew==0d0) then
                    a=1d0
                else
                    a=-fnew/phidnew
                end if
            do iIteration=1,MaxConjugateGradientIteration
                !Prepare
                fold=fnew
                fdold=fdnew
                phidold=phidnew
                !Line search
                call StrongWolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim)
                phidnew=dot_product(fdnew,fdnew)
                if(phidnew<ConjugateGradientTol) return
                if(dot_product(p,p)*a*a<ConjugateGradientTol) then
                    if(ConjugateGradientWarning) then
                        write(*,'(1x,A98)')'Conjugate gradient warning: step length has converged, but gradient norm has not met accuracy goal'
                        write(*,'(1x,A56)')'A best estimation rather than exact solution is returned'
                        write(*,*)'Euclidean norm of gradient =',Sqrt(phidnew)
                        ConjugateGradientWarning=.false.
                    end if
                    return 
                end if
                !Determine new direction
                p=-fdnew+dot_product(fdnew,fdnew)/dot_product(fdnew-fdold,p)*p
                phidnew=dot_product(fdnew,p)
                if(phidnew>0d0) then!Ascent direction, reset to steepest descent direction
                    p=-fdnew
                    phidnew=-dot_product(fdnew,fdnew)
                end if
                a=a*phidold/phidnew
            end do
            if(iIteration>MaxConjugateGradientIteration.and.ConjugateGradientWarning) then
                write(*,*)'Failed conjugate gradient: max iteration exceeded!'
                write(*,*)'Euclidean norm of gradient =',Norm2(fdnew)
                ConjugateGradientWarning=.false.
            end if
        end subroutine DYConjugateGradient_Strong
        !When it is cheap to evaluate f' along with f
        subroutine DYConjugateGradient_Strong_fdwithf(f,fd,f_fd,x,dim)
            external::f,fd,f_fd
            integer,intent(in)::dim
            real*8,dimension(dim),intent(inout)::x
            real*8,parameter::c1=1d-4,c2=0.45d0! 0 < c1 < c2 < 0.5: Wolfe constant
            integer::iIteration
            real*8::a,fnew,fold,phidnew,phidold
            real*8,dimension(dim)::p,fdnew,fdold
            !Initialize
                call f_fd(fnew,fdnew,x,dim)
                p=-fdnew
                phidnew=-dot_product(fdnew,fdnew)
                if(-phidnew<ConjugateGradientTol) return
                if(fnew==0d0) then
                    a=1d0
                else
                    a=-fnew/phidnew
                end if
            do iIteration=1,MaxConjugateGradientIteration
                !Prepare
                fold=fnew
                fdold=fdnew
                phidold=phidnew
                !Line search
                call StrongWolfe_fdwithf(c1,c2,f,fd,f_fd,x,a,p,fnew,phidnew,fdnew,dim)
                phidnew=dot_product(fdnew,fdnew)
                if(phidnew<ConjugateGradientTol) return
                if(dot_product(p,p)*a*a<ConjugateGradientTol) then
                    if(ConjugateGradientWarning) then
                        write(*,'(1x,A98)')'Conjugate gradient warning: step length has converged, but gradient norm has not met accuracy goal'
                        write(*,'(1x,A56)')'A best estimation rather than exact solution is returned'
                        write(*,*)'Euclidean norm of gradient =',Sqrt(phidnew)
                        ConjugateGradientWarning=.false.
                    end if
                    return 
                end if
                !Determine new direction
                p=-fdnew+dot_product(fdnew,fdnew)/dot_product(fdnew-fdold,p)*p
                phidnew=dot_product(fdnew,p)
                if(phidnew>0d0) then!Ascent direction, reset to steepest descent direction
                    p=-fdnew
                    phidnew=-dot_product(fdnew,fdnew)
                end if
                a=a*phidold/phidnew
            end do
            if(iIteration>MaxConjugateGradientIteration.and.ConjugateGradientWarning) then
                write(*,*)'Failed conjugate gradient: max iteration exceeded!'
                write(*,*)'Euclidean norm of gradient =',Norm2(fdnew)
                ConjugateGradientWarning=.false.
            end if
        end subroutine DYConjugateGradient_Strong_fdwithf

        !On input x is an initial guess of the minimum point of f(x), on exit x is the minimum point
        !Polak-Ribiere+ conjugate gradient method, requiring strong Wolfe condition
        subroutine PRConjugateGradient(f,fd,x,dim)
            external::f,fd
            integer,intent(in)::dim
            real*8,dimension(dim),intent(inout)::x
            real*8,parameter::c1=1d-4,c2=0.45d0! 0 < c1 < c2 < 0.5: Wolfe constant
            integer::iIteration
            real*8::a,fnew,fold,phidnew,phidold
            real*8,dimension(dim)::p,fdnew,fdold
            !Initialize
                call f(fnew,x,dim)
                call fd(fdnew,x,dim)
                p=-fdnew
                phidnew=-dot_product(fdnew,fdnew)
                if(-phidnew<ConjugateGradientTol) return
                if(fnew==0d0) then
                    a=1d0
                else
                    a=-fnew/phidnew
                end if
            do iIteration=1,MaxConjugateGradientIteration
                !Prepare
                fold=fnew
                fdold=fdnew
                phidold=phidnew
                !Line search
                call StrongWolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim)
                phidnew=dot_product(fdnew,fdnew)
                if(phidnew<ConjugateGradientTol) return
                if(dot_product(p,p)*a*a<ConjugateGradientTol) then
                    if(ConjugateGradientWarning) then
                        write(*,'(1x,A98)')'Conjugate gradient warning: step length has converged, but gradient norm has not met accuracy goal'
                        write(*,'(1x,A56)')'A best estimation rather than exact solution is returned'
                        write(*,*)'Euclidean norm of gradient =',Sqrt(phidnew)
                        ConjugateGradientWarning=.false.
                    end if
                    return 
                end if
                !Determine new direction
                p=-fdnew+dot_product(fdnew,fdnew-fdold)/dot_product(fdold,fdold)*p
                phidnew=dot_product(fdnew,p)
                if(phidnew>0d0) then!Ascent direction, reset to steepest descent direction
                    p=-fdnew
                    phidnew=-dot_product(fdnew,fdnew)
                end if
                a=a*phidold/phidnew
            end do
            if(iIteration>MaxConjugateGradientIteration.and.ConjugateGradientWarning) then
                write(*,*)'Failed conjugate gradient: max iteration exceeded!'
                write(*,*)'Euclidean norm of gradient =',Norm2(fdnew)
                ConjugateGradientWarning=.false.
            end if
        end subroutine PRConjugateGradient
        !When it is cheap to evaluate f' along with f
        subroutine PRConjugateGradient_fdwithf(f,fd,f_fd,x,dim)
            external::f,fd,f_fd
            integer,intent(in)::dim
            real*8,dimension(dim),intent(inout)::x
            real*8,parameter::c1=1d-4,c2=0.45d0! 0 < c1 < c2 < 0.5: Wolfe constant
            integer::iIteration
            real*8::a,fnew,fold,phidnew,phidold
            real*8,dimension(dim)::p,fdnew,fdold
            !Initialize
                call f_fd(fnew,fdnew,x,dim)
                p=-fdnew
                phidnew=-dot_product(fdnew,fdnew)
                if(-phidnew<ConjugateGradientTol) return
                if(fnew==0d0) then
                    a=1d0
                else
                    a=-fnew/phidnew
                end if
            do iIteration=1,MaxConjugateGradientIteration
                !Prepare
                fold=fnew
                fdold=fdnew
                phidold=phidnew
                !Line search
                call StrongWolfe_fdwithf(c1,c2,f,fd,f_fd,x,a,p,fnew,phidnew,fdnew,dim)
                phidnew=dot_product(fdnew,fdnew)
                if(phidnew<ConjugateGradientTol) return
                if(dot_product(p,p)*a*a<ConjugateGradientTol) then
                    if(ConjugateGradientWarning) then
                        write(*,'(1x,A98)')'Conjugate gradient warning: step length has converged, but gradient norm has not met accuracy goal'
                        write(*,'(1x,A56)')'A best estimation rather than exact solution is returned'
                        write(*,*)'Euclidean norm of gradient =',Sqrt(phidnew)
                        ConjugateGradientWarning=.false.
                    end if
                    return 
                end if
                !Determine new direction
                p=-fdnew+dot_product(fdnew,fdnew-fdold)/dot_product(fdold,fdold)*p
                phidnew=dot_product(fdnew,p)
                if(phidnew>0d0) then!Ascent direction, reset to steepest descent direction
                    p=-fdnew
                    phidnew=-dot_product(fdnew,fdnew)
                end if
                a=a*phidold/phidnew
            end do
            if(iIteration>MaxConjugateGradientIteration.and.ConjugateGradientWarning) then
                write(*,*)'Failed conjugate gradient: max iteration exceeded!'
                write(*,*)'Euclidean norm of gradient =',Norm2(fdnew)
                ConjugateGradientWarning=.false.
            end if
        end subroutine PRConjugateGradient_fdwithf
    !================ End =================
!------------------ End -------------------

!-------------- Trust region --------------
    !MKL trust-region nonlinear least square problem (trnlsp) solver wrapper
    !Solve f'(x) = 0 by minimizing merit function F(x) = f'(x)^2
    !    through trust-region method with model function m(p) = [ f'(x) + J(x) . p ]^2
    !    where f'(x) is M dimensional vector, x is N dimensional vector, J(x) is the Jacobian
    !This procedure can be interpreted in different ways other than f'(x) = 0:
    !    f'(x) can be viewed as Sqrt(weight)*residual terms,
    !        then trnlsp minimizes the square penalty (this is where its name comes from)
    !    When M = N, f'(x) can also be considered as the gradient of f(x), Jacobian = Hessian,
    !        but trnlsp doesn't necessarily optimize f(x) (unless f(x) = const * F(x) by coincidence)
    !        because the zero point of F(x) is merely the stationary point of f(x)
    
    !fd has the form of: subroutine fd(f'(x),x,M,N)
    !Jacobian has the form of: subroutine Jacobian(J(x),x,M,N)
    !Where f'(x) is M dimensional vector, x is N dimensional vector
    !On input x is an initial guess of f'(x) = 0, on exit x is the solution
    subroutine My_dtrnlsp(fd,Jacobian,x,M,N)
        external::fd,Jacobian
        integer,intent(in)::M,N
        real*8,dimension(N),intent(inout)::x
        !TotalIteration harvests the solver stops after how many interations
        !StopReason harvests why the solver has stopped:
        !    1,   max iteration exceeded
        !    2-6, MKLTrustRegionTol(StopReason-1) is met
        integer::TotalIteration,StopReason
        real*8::InitialResidual,FinalResidual,StepBound
        real*8,dimension(M)::fdx
        real*8,dimension(M,N)::J
        integer::RCI_request!Reverse communication interface parameter
        integer,dimension(6)::info!Results of input parameter checking
        type(handle_tr)::handle!Trust-region solver handle
        !Initialize
            fdx=0d0
            J=0d0
            StepBound=100d0
            if(dtrnlsp_init(handle,N,M,x,MKLTrustRegionTol,MaxMKLTrustRegionIteration,MaxMKLTrialStepIteration,StepBound)/=TR_SUCCESS) then
                call mkl_free_buffers
                return
            end if
            if(dtrnlsp_check(handle,N,M,J,fdx,MKLTrustRegionTol,info)/=TR_SUCCESS) then
                call mkl_free_buffers
                return
            else
                if(info(1)/=0.or.info(2)/=0.or.info(3)/=0.or.info(4)/=0) then
                    call mkl_free_buffers
                    return
                end if
            end if
            RCI_request=0
        do
            if (dtrnlsp_solve(handle,fdx,J,RCI_request)/=TR_SUCCESS) then
                call mkl_free_buffers
                return
            end if
            select case (RCI_request)
                case (-1,-2,-3,-4,-5,-6)
                    exit
                case (1)
                    call fd(fdx,x,M,N)
                case (2)
                    call Jacobian(J,x,M,N)
            end select
        end do
        !Clean up
            if (dtrnlsp_get(handle,TotalIteration,StopReason,InitialResidual,FinalResidual)/=TR_SUCCESS) then
                call mkl_free_buffers
                return
            end if
            if (dtrnlsp_delete(handle)/=TR_SUCCESS) then
                call mkl_free_buffers
                return
            end if
            call mkl_free_buffers
        if(StopReason/=3) then
            select case(StopReason)
                case(1)
                    if(trnlspWarning) then
                        write(*,*)'Failed trust region: max iteration exceeded!'
                        write(*,*)'Final residual =',FinalResidual
                        trnlspWarning=.false.
                    end if
                case(4)
                    if(trnlspWarning) then
                        write(*,*)'Failed trust region: singular Jacobian encountered!'
                        write(*,*)'Final residual =',FinalResidual
                        trnlspWarning=.false.
                    end if
                case default
                    if(trnlspWarning) then
                        write(*,'(1x,A87)')'Trust region warning: step length has converged, but residual has not met accuracy goal'
                        write(*,'(1x,A56)')'A best estimation rather than exact solution is returned'
                        write(*,*)'Final residual =',FinalResidual
                        trnlspWarning=.false.
                    end if
            end select
        end if
    end subroutine My_dtrnlsp
    
    !Cost warning: Jacobian is computed through central difference by calling djacobi
    !fd has the form of: subroutine fd(M,N,x,f'(x))
    !Where x is N dimensional vector, f'(x) is M dimensional vector
    !On input x is an initial guess of f'(x) = 0, on exit x is the solution
    subroutine My_dtrnlsp_NumericalJacobian(fd,x,M,N)
        external::fd
        integer,intent(in)::M,N
        real*8,dimension(N),intent(inout)::x
        !TotalIteration harvests the solver stops after how many interations
        !StopReason harvests why the solver has stopped:
        !    1,   max iteration exceeded
        !    2-6, MKLTrustRegionTol(StopReason-1) is met
        integer::TotalIteration,StopReason
        real*8::InitialResidual,FinalResidual,StepBound
        real*8,dimension(M)::fdx
        real*8,dimension(M,N)::J
        integer::RCI_request!Reverse communication interface parameter
        integer,dimension(6)::info!Results of input parameter checking
        type(handle_tr)::handle!Trust-region solver handle
        !Initialize
            fdx=0d0
            J=0d0
            StepBound=100d0
            if(dtrnlsp_init(handle,N,M,x,MKLTrustRegionTol,MaxMKLTrustRegionIteration,MaxMKLTrialStepIteration,StepBound)/=TR_SUCCESS) then
                call mkl_free_buffers
                return
            end if
            if(dtrnlsp_check(handle,N,M,J,fdx,MKLTrustRegionTol,info)/=TR_SUCCESS) then
                call mkl_free_buffers
                return
            else
                if(info(1)/=0.or.info(2)/=0.or.info(3)/=0.or.info(4)/=0) then
                    call mkl_free_buffers
                    return
                end if
            end if
            RCI_request=0
        do
            if (dtrnlsp_solve(handle,fdx,J,RCI_request)/=TR_SUCCESS) then
                call mkl_free_buffers
                return
            end if
            select case (RCI_request)
                case (-1,-2,-3,-4,-5,-6)
                    exit
                case (1)
                    call fd(M,N,x,fdx)
                case (2)
                    if(djacobi(fd,N,M,J,x,jacobiPrecision)/=TR_SUCCESS) then
                        call mkl_free_buffers
                        return
                    end if
            end select
        end do
        !Clean up
            if(dtrnlsp_get(handle,TotalIteration,StopReason,InitialResidual,FinalResidual)/=TR_SUCCESS) then
                call mkl_free_buffers
                return
            end if
            if(dtrnlsp_delete(handle)/=TR_SUCCESS) then
                call mkl_free_buffers
                return
            end if
            call mkl_free_buffers
        if(StopReason/=3) then
            select case(StopReason)
                case(1)
                    if(trnlspWarning) then
                        write(*,*)'Failed trust region: max iteration exceeded!'
                        write(*,*)'Final residual =',FinalResidual
                        trnlspWarning=.false.
                    end if
                case(4)
                    if(trnlspWarning) then
                        write(*,*)'Failed trust region: singular Jacobian encountered!'
                        write(*,*)'Final residual =',FinalResidual
                        trnlspWarning=.false.
                    end if
                case default
                    if(trnlspWarning) then
                        write(*,'(1x,A87)')'Trust region warning: step length has converged, but residual has not met accuracy goal'
                        write(*,'(1x,A56)')'A best estimation rather than exact solution is returned'
                        write(*,*)'Final residual =',FinalResidual
                        trnlspWarning=.false.
                    end if
            end select
        end if
    end subroutine My_dtrnlsp_NumericalJacobian
    
    !fd has the form of: subroutine fd(f'(x),x,M,N)
    !Jacobian has the form of: subroutine Jacobian(J(x),x,M,N)
    !Where f'(x) is M dimensional vector, x is N dimensional vector
    !Solve f'(x) = 0 with boundary condition low(i) <= x(i) <= up(i) for all i in [1,N]
    !low and up must satisfy low(i) < up(i) for all i in [1,N]
    !On input x is an initial guess of f'(x) = 0, on exit x is the solution
    subroutine My_dtrnlspbc(fd,Jacobian,x,low,up,M,N)
        external::fd,Jacobian
        integer,intent(in)::M,N
        real*8,dimension(N),intent(inout)::x
        real*8,dimension(N),intent(in)::low,up
        !TotalIteration harvests the solver stops after how many interations
        !StopReason harvests why the solver has stopped:
        !    1,   max iteration exceeded
        !    2-6, MKLTrustRegionTol(StopReason-1) is met
        integer::TotalIteration,StopReason
        real*8::InitialResidual,FinalResidual,StepBound
        real*8,dimension(M)::fdx
        real*8,dimension(M,N)::J
        integer::RCI_request!Reverse communication interface parameter
        integer,dimension(6)::info!Results of input parameter checking
        type(handle_tr)::handle!Trust-region solver handle
        !Initialize
            fdx=0d0
            J=0d0
            StepBound=100d0
            if(dtrnlspbc_init(handle,N,M,x,low,up,MKLTrustRegionTol,MaxMKLTrustRegionIteration,MaxMKLTrialStepIteration,StepBound)/=TR_SUCCESS) then
                call mkl_free_buffers
                return
            end if
            if(dtrnlspbc_check(handle,N,M,J,fdx,low,up,MKLTrustRegionTol,info)/=TR_SUCCESS) then
                call mkl_free_buffers
                return
            else
                if(info(1)/=0.or.info(2)/=0.or.info(3)/=0.or.info(4)/=0.or.info(5)/=0.or.info(6)/=0) then
                    call mkl_free_buffers
                    return
                end if
            end if
            RCI_request=0
        do
            if (dtrnlspbc_solve(handle,fdx,J,RCI_request)/=TR_SUCCESS) then
                call mkl_free_buffers
                return
            end if
            select case (RCI_request)
                case (-1,-2,-3,-4,-5,-6)
                    exit
                case (1)
                    call fd(fdx,x,M,N)
                case (2)
                    call Jacobian(J,x,M,N)
            end select
        end do
        !Clean up
            if (dtrnlspbc_get(handle,TotalIteration,StopReason,InitialResidual,FinalResidual)/=TR_SUCCESS) then
                call mkl_free_buffers
                return
            end if
            if (dtrnlspbc_delete(handle)/=TR_SUCCESS) then
                call mkl_free_buffers
                return
            end if
            call mkl_free_buffers
        do RCI_request=1,N
            if((x(RCI_request)<low(RCI_request).or.x(RCI_request)>up(RCI_request)).and.trnlspWarning) then
                write(*,*)'Failed trust region: boundary condition violated!'
                trnlspWarning=.false.
                exit
            end if
        end do
        if(StopReason/=3) then
            select case(StopReason)
                case(1)
                    if(trnlspWarning) then
                        write(*,*)'Failed trust region: max iteration exceeded!'
                        write(*,*)'Final residual =',FinalResidual
                        trnlspWarning=.false.
                    end if
                case(4)
                    if(trnlspWarning) then
                        write(*,*)'Failed trust region: singular Jacobian encountered!'
                        write(*,*)'Final residual =',FinalResidual
                        trnlspWarning=.false.
                    end if
                case default
                    if(trnlspWarning) then
                        write(*,'(1x,A87)')'Trust region warning: step length has converged, but residual has not met accuracy goal'
                        write(*,'(1x,A56)')'A best estimation rather than exact solution is returned'
                        write(*,*)'Final residual =',FinalResidual
                        trnlspWarning=.false.
                    end if
            end select
        end if
    end subroutine My_dtrnlspbc
    
    !Cost warning: Jacobian is computed through central difference by calling djacobi
    !fd has the form of: subroutine fd(M,N,x,f'(x))
    !Where f'(x) is M dimensional vector, x is N dimensional vector
    !Solve f'(x) = 0 with boundary condition low(i) <= x(i) <= up(i) for all i in [1,N]
    !low and up must satisfy low(i) < up(i) for all i in [1,N]
    !On input x is an initial guess of f'(x) = 0, on exit x is the solution
    subroutine My_dtrnlspbc_NumericalJacobian(fd,x,low,up,M,N)
        external::fd
        integer,intent(in)::M,N
        real*8,dimension(N),intent(inout)::x
        real*8,dimension(N),intent(in)::low,up
        !TotalIteration harvests the solver stops after how many interations
        !StopReason harvests why the solver has stopped:
        !    1,   max iteration exceeded
        !    2-6, MKLTrustRegionTol(StopReason-1) is met
        integer::TotalIteration,StopReason
        real*8::InitialResidual,FinalResidual,StepBound
        real*8,dimension(M)::fdx
        real*8,dimension(M,N)::J
        integer::RCI_request!Reverse communication interface parameter
        integer,dimension(6)::info!Results of input parameter checking
        type(handle_tr)::handle!Trust-region solver handle
        !Initialize
            fdx=0d0
            J=0d0
            StepBound=100d0
            if(dtrnlspbc_init(handle,N,M,x,low,up,MKLTrustRegionTol,MaxMKLTrustRegionIteration,MaxMKLTrialStepIteration,StepBound)/=TR_SUCCESS) then
                call mkl_free_buffers
                return
            end if
            if(dtrnlspbc_check(handle,N,M,J,fdx,low,up,MKLTrustRegionTol,info)/=TR_SUCCESS) then
                call mkl_free_buffers
                return
            else
                if(info(1)/=0.or.info(2)/=0.or.info(3)/=0.or.info(4)/=0.or.info(5)/=0.or.info(6)/=0) then
                    call mkl_free_buffers
                    return
                end if
            end if
            RCI_request=0
        do
            if (dtrnlspbc_solve(handle,fdx,J,RCI_request)/=TR_SUCCESS) then
                call mkl_free_buffers
                return
            end if
            select case (RCI_request)
                case (-1,-2,-3,-4,-5,-6)
                    exit
                case (1)
                    call fd(M,N,x,fdx)
                case (2)
                    if(djacobi(fd,N,M,J,x,jacobiPrecision)/=TR_SUCCESS) then
                        call mkl_free_buffers
                        return
                    end if
            end select
        end do
        !Clean up
            if (dtrnlspbc_get(handle,TotalIteration,StopReason,InitialResidual,FinalResidual)/=TR_SUCCESS) then
                call mkl_free_buffers
                return
            end if
            if (dtrnlspbc_delete(handle)/=TR_SUCCESS) then
                call mkl_free_buffers
                return
            end if
            call mkl_free_buffers
        do RCI_request=1,N
            if((x(RCI_request)<low(RCI_request).or.x(RCI_request)>up(RCI_request)).and.trnlspWarning) then
                write(*,*)'Failed trust region: boundary condition violated!'
                trnlspWarning=.false.
                exit
            end if
        end do
        if(StopReason/=3) then
            select case(StopReason)
                case(1)
                    if(trnlspWarning) then
                        write(*,*)'Failed trust region: max iteration exceeded!'
                        write(*,*)'Final residual =',FinalResidual
                        trnlspWarning=.false.
                    end if
                case(4)
                    if(trnlspWarning) then
                        write(*,*)'Failed trust region: singular Jacobian encountered!'
                        write(*,*)'Final residual =',FinalResidual
                        trnlspWarning=.false.
                    end if
                case default
                    if(trnlspWarning) then
                        write(*,'(1x,A87)')'Trust region warning: step length has converged, but residual has not met accuracy goal'
                        write(*,'(1x,A56)')'A best estimation rather than exact solution is returned'
                        write(*,*)'Final residual =',FinalResidual
                        trnlspWarning=.false.
                    end if
            end select
        end if
    end subroutine My_dtrnlspbc_NumericalJacobian
!------------------ End -------------------

!---------- Augmented Lagrangian ----------
    !Textbook is wrong: it claims that Lagrangian multiplier method
    !transforms a constrained optimization problem into an unconstrained one
    !However, L does not have a lower bound, it can be -infinity
    !when multiplier approaches ±infinity and constraint is not satisfied,
    !which makes all unconstrained minimizers fail
    !Lagrangian multiplier method actually turns a minimization problem into a saddle point problem,
    !which can only be solved by L'(x) = 0 rather than try decreasing f(x),
    !but the solution of L'(x) = 0 is not necessarily a minimum.
    !So Lagrangian multiplier method is only good when L has unique saddle point
    !This is why we have to turn to the augmented Lagrangian method
    
    !=========== Equality ===========
        !I will write someday (flag)
    !============= End ==============

    !========== Inequality ==========
        !I will write someday (flag)
    !============= End ==============
!------------------ End -------------------

!--------------- Heuristic ----------------
    !I will write someday (flag)
!------------------ End -------------------

end module NonlinearOptimization