!Nonlinear optimization routines
!
!Instruction:
!If you want to solve nonlinear equation(s), you may:
!    adopt MKL trust region solver
!    or build your own merit function then use other routines to search for its minimum
!This module mainly provides unconstrained local minimizer
!    only solvers in Augmented Lagrangian section are constrained
!    only solvers in Heuristic section are global
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
    !Suggestion:
    !    If dimensionality is low, adopt quasi-Newton (or Newton if Hessian is cheap and initial guess is close)
    !    If dimensionality is so high that O(dim^2) memory is unaffordable, adopt conjugate gradient or L-BFGS
    !Nomenclature:
    !    f = the target function to be minimized
    !    a = the line search step length
    !    p = the line search direction
    !    phi(a) = f( x + a * p ), so phi'(a) = f'( x + a * p ) . p
    !External procedure format:
    !    subroutine f(f(x),x,dim)
    !    subroutine fd(f'(x),x,dim)
    !    integer function f_fd(f(x),f'(x),x,dim)
    !    integer function fdd(f''(x),x,dim)
    !    dim dimensional vector x & f'(x), dim order matrix f''(x)
    !Required argument:
    !    external subroutine f & fd, integer dim, dim dimensional vector x
    !Common optional argument:
    !    f_fd: presence means evaluating f(x) & f'(x) together is cheaper than separately,
    !          strong Wolfe condition is thus applied, because Wolfe does not benefit from this
    !    Strong: (default = true) if true, use strong Wolfe condition instead of Wolfe condition
    !    Warning: (default = true) if false, all warnings will be suppressed
    !    MaxIteration: (default = 1000) max number of iterations to perform
    !    Tolerance: (default = 1d-15) convergence considered when || f'(x) ||_2 < Tolerance
    !               if search step < Tolerance before converging, terminate and throw a warning
    !    WolfeConst1 & WolfeConst2: 0 < WolfeConst1 < WolfeConst2 <  1  for Newton & quasi-Newton
    !                               0 < WolfeConst1 < WolfeConst2 < 0.5 for conjugate gradient
    !On input x is an initial guess, on exit x is a local minimum of f(x)

    !Newton-Raphson method, requiring Wolfe condition
    !Optional argument:
    !    fdd: presence means analytical Hessian is available, otherwise call djacobi for central difference Hessian
    subroutine NewtonRaphson(f,fd,x,dim,fdd,f_fd,Strong,Warning,MaxIteration,Tolerance,WolfeConst1,WolfeConst2)
        !Required argument
            external::f,fd
            integer,intent(in)::dim
            real*8,dimension(dim),intent(inout)::x
        !Optional argument
            integer,external,optional::fdd,f_fd
            logical,intent(in),optional::Strong,Warning
            integer,intent(in),optional::MaxIteration
            real*8,intent(in),optional::Tolerance,WolfeConst1,WolfeConst2
        logical::sw,warn,terminate
        integer::maxit,iIteration,info
        real*8::tol,c1,c2,a,fnew,phidnew,phidold
        real*8,dimension(dim)::p,fdnew
        real*8,dimension(dim,dim)::Hessian
        call Initialize()
        if(terminate) return
        if(present(fdd)) then!Analytical Hessian available
            if(present(f_fd)) then!Cheaper to evaluate f' along with f
                do iIteration=1,maxit!Main loop
                    phidold=phidnew!Prepare
                    call StrongWolfe_fdwithf(c1,c2,f,fd,f_fd,x,a,p,fnew,phidnew,fdnew,dim)!Line search
                    call After()!After search
                    if(terminate) return
                end do
            else
                if(sw) then!Use strong Wolfe condition instead of Wolfe condition
                    do iIteration=1,maxit!Main loop
                        phidold=phidnew!Prepare
                        call StrongWolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim)!Line search
                        call After()!After search
                        if(terminate) return
                    end do
                else
                    do iIteration=1,maxit!Main loop
                        phidold=phidnew!Prepare
                        call Wolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim)!Line search
                        call After()!After search
                        if(terminate) return
                    end do
                end if
            end if
        else
            if(present(f_fd)) then!Cheaper to evaluate f' along with f
                do iIteration=1,maxit!Main loop
                    phidold=phidnew!Prepare
                    call StrongWolfe_fdwithf(c1,c2,f,fd,f_fd,x,a,p,fnew,phidnew,fdnew,dim)!Line search
                    call After_NumericalHessian()!After search
                    if(terminate) return
                end do
            else
                if(sw) then!Use strong Wolfe condition instead of Wolfe condition
                    do iIteration=1,maxit!Main loop
                        phidold=phidnew!Prepare
                        call StrongWolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim)!Line search
                        call After_NumericalHessian()!After search
                        if(terminate) return
                    end do
                else
                    do iIteration=1,maxit!Main loop
                        phidold=phidnew!Prepare
                        call Wolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim)!Line search
                        call After_NumericalHessian()!After search
                        if(terminate) return
                    end do
                end if
            end if
        end if
        if(iIteration>maxit.and.warn) then
            write(*,*)'Failed Newton-Raphson: max iteration exceeded!'
            write(*,*)'Euclidean norm of gradient =',Norm2(fdnew)
        end if
        contains
            subroutine Initialize()!Parameters, initial f(x) & f'(x), initial direction & step length
                terminate=.false.
                !Set parameter according to optional argument
                    if(present(Strong)) then
                        sw=Strong
                    else
                        sw=.true.
                    end if
                    if(present(Warning)) then
                        warn=Warning
                    else
                        warn=.true.
                    end if
                    if(present(MaxIteration)) then
                        maxit=MaxIteration
                    else
                        maxit=1000
                    end if
                    if(present(Tolerance)) then!To save sqrt cost, tolerance is squared
                        tol=Tolerance*Tolerance
                    else
                        tol=1d-30
                    end if
                    if(present(WolfeConst1)) then
                        c1=WolfeConst1
                    else
                        c1=1d-4
                    end if
                    if(present(WolfeConst2)) then
                        c2=WolfeConst2
                    else
                        c2=0.9d0
                    end if
                if(present(f_fd)) then
                    info=f_fd(fnew,fdnew,x,dim)
                else
                    call f(fnew,x,dim)
                    call fd(fdnew,x,dim)
                end if
                if(present(fdd)) then
                    info=fdd(Hessian,x,dim)
                else
                    info=djacobi(fd_j,dim,dim,Hessian,x,1d-8)
                end if
                p=-fdnew
                call My_dposv(Hessian,p,dim,info)
                if(info==0) then
                    phidnew=dot_product(fdnew,p)
                    a=1d0
                else!Hessian is not positive definite, use steepest descent direction
                    p=-fdnew
                    phidnew=-dot_product(fdnew,fdnew)
                    if(-phidnew<tol) then
                        terminate=.true.
                        return
                    end if
                    if(fnew==0d0) then
                        a=1d0
                    else
                        a=-fnew/phidnew
                    end if
                end if
            end subroutine Initialize
            subroutine After()!Check convergence, determine new direction & step length
                !Check convergence
                phidnew=dot_product(fdnew,fdnew)
                if(phidnew<tol) then
                    terminate=.true.
                    return
                end if
                if(dot_product(p,p)*a*a<tol) then
                    if(warn) then
                        write(*,'(1x,A94)')'Newton-Raphson warning: step length has converged, but gradient norm has not met accuracy goal'
                        write(*,'(1x,A56)')'A best estimation rather than exact solution is returned'
                        write(*,*)'Euclidean norm of gradient =',Sqrt(phidnew)
                    end if
                    terminate=.true.
                    return
                end if
                !Determine new direction and step length
                p=-fdnew
                info=fdd(Hessian,x,dim)
                call My_dposv(Hessian,p,dim,info)
                if(info==0) then
                    phidnew=dot_product(fdnew,p)
                    a=1d0
                else!Hessian is not positive definite, use steepest descent direction
                    p=-fdnew
                    phidnew=-phidnew
                    a=a*phidold/phidnew
                end if
            end subroutine After
            subroutine After_NumericalHessian()!Check convergence, determine new direction & step length
                !Check convergence
                phidnew=dot_product(fdnew,fdnew)
                if(phidnew<tol) then
                    terminate=.true.
                    return
                end if
                if(dot_product(p,p)*a*a<tol) then
                    if(warn) then
                        write(*,'(1x,A94)')'Newton-Raphson warning: step length has converged, but gradient norm has not met accuracy goal'
                        write(*,'(1x,A56)')'A best estimation rather than exact solution is returned'
                        write(*,*)'Euclidean norm of gradient =',Sqrt(phidnew)
                    end if
                    terminate=.true.
                    return
                end if
                !Determine new direction and step length
                p=-fdnew
                info=djacobi(fd_j,dim,dim,Hessian,x,1d-8)
                call My_dposv(Hessian,p,dim,info)
                if(info==0) then
                    phidnew=dot_product(fdnew,p)
                    a=1d0
                else!Hessian is not positive definite, use steepest descent direction
                    p=-fdnew
                    phidnew=-phidnew
                    a=a*phidold/phidnew
                end if
            end subroutine After_NumericalHessian
            subroutine fd_j(M,N,x,fdx)!Reformat fd for djacobi
                integer,intent(in)::M,N
                real*8,dimension(N),intent(in)::x
                real*8,dimension(M),intent(out)::fdx
                call fd(fdx,x,N)
            end subroutine fd_j
    end subroutine NewtonRaphson

    !Broyden–Fletcher–Goldfarb–Shanno (BFGS) quasi-Newton method, requiring Wolfe condition
    !Optional argument:
    !    fdd: presence means analytical Hessian is available, otherwise call djacobi for central difference Hessian
    !    ExactStep: (default = 20) every how many steps compute exact hessian. [10,50] is recommended
    !               if ExactStep <= 0, exact Hessian will not be computed at all
    subroutine BFGS(f,fd,x,dim,fdd,ExactStep,f_fd,Strong,Warning,MaxIteration,Tolerance,WolfeConst1,WolfeConst2)
        !Required argument
            external::f,fd
            integer,intent(in)::dim
            real*8,dimension(dim),intent(inout)::x
        !Optional argument
            integer,external,optional::fdd,f_fd
            logical,intent(in),optional::Strong,Warning
            integer,intent(in),optional::ExactStep,MaxIteration
            real*8,intent(in),optional::Tolerance,WolfeConst1,WolfeConst2
        logical::sw,warn,terminate
        integer::freq,maxit,iIteration,i
        real*8::tol,c1,c2,a,fnew,phidnew,rho
        real*8,dimension(dim)::p,fdnew,s,y
        real*8,dimension(dim,dim)::U,H!Approximate inverse Hessian
        call Initialize()
        if(terminate) return
        if(freq>0) then!Exact Hessian will be computed
            if(present(fdd)) then!Analytical Hessian is available
                if(present(f_fd)) then!Cheaper to evaluate f' along with f
                    do iIteration=1,maxit!Main loop
                        s=x!Prepare
                        y=fdnew
                        call StrongWolfe_fdwithf(c1,c2,f,fd,f_fd,x,a,p,fnew,phidnew,fdnew,dim)!Line search
                        call After()!After search
                        if(terminate) return
                    end do
                else
                    if(sw) then!Use strong Wolfe condition instead of Wolfe condition
                        do iIteration=1,maxit!Main loop
                            s=x!Prepare
                            y=fdnew
                            call StrongWolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim)!Line search
                            call After()!After search
                            if(terminate) return
                        end do
                    else
                        do iIteration=1,maxit!Main loop
                            s=x!Prepare
                            y=fdnew
                            call Wolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim)!Line search
                            call After()!After search
                            if(terminate) return
                        end do
                    end if
                end if
            else
                if(present(f_fd)) then!Cheaper to evaluate f' along with f
                    do iIteration=1,maxit!Main loop
                        s=x!Prepare
                        y=fdnew
                        call StrongWolfe_fdwithf(c1,c2,f,fd,f_fd,x,a,p,fnew,phidnew,fdnew,dim)!Line search
                        call After_NumericalHessian()!After search
                        if(terminate) return
                    end do
                else
                    if(sw) then!Use strong Wolfe condition instead of Wolfe condition
                        do iIteration=1,maxit!Main loop
                            s=x!Prepare
                            y=fdnew
                            call StrongWolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim)!Line search
                            call After_NumericalHessian()!After search
                            if(terminate) return
                        end do
                    else
                        do iIteration=1,maxit!Main loop
                            s=x!Prepare
                            y=fdnew
                            call Wolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim)!Line search
                            call After_NumericalHessian()!After search
                            if(terminate) return
                        end do
                    end if
                end if
            end if
        else
            if(present(f_fd)) then!Cheaper to evaluate f' along with f
                do iIteration=1,maxit!Main loop
                    s=x!Prepare
                    y=fdnew
                    call StrongWolfe_fdwithf(c1,c2,f,fd,f_fd,x,a,p,fnew,phidnew,fdnew,dim)!Line search
                    call After_NoHessian()!After search
                    if(terminate) return
                end do
            else
                if(sw) then!Use strong Wolfe condition instead of Wolfe condition
                    do iIteration=1,maxit!Main loop
                        s=x!Prepare
                        y=fdnew
                        call StrongWolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim)!Line search
                        call After_NoHessian()!After search
                        if(terminate) return
                    end do
                else
                    do iIteration=1,maxit!Main loop
                        s=x!Prepare
                        y=fdnew
                        call Wolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim)!Line search
                        call After_NoHessian()!After search
                        if(terminate) return
                    end do
                end if
            end if
        end if
        if(iIteration>maxit.and.warn) then
            write(*,*)'Failed BFGS: max iteration exceeded!'
            write(*,*)'Euclidean norm of gradient =',Norm2(fdnew)
        end if
        contains
            subroutine Initialize()!Parameters, initial f(x) & f'(x), initial approximate inverse Hessian & direction & step length
                terminate=.false.
                !Set parameter according to optional argument
                    if(present(ExactStep)) then
                        freq=ExactStep
                    else
                        freq=20
                    end if
                    if(present(Strong)) then
                        sw=Strong
                    else
                        sw=.true.
                    end if
                    if(present(Warning)) then
                        warn=Warning
                    else
                        warn=.true.
                    end if
                    if(present(MaxIteration)) then
                        maxit=MaxIteration
                    else
                        maxit=1000
                    end if
                    if(present(Tolerance)) then!To save sqrt cost, tolerance is squared
                        tol=Tolerance*Tolerance
                    else
                        tol=1d-30
                    end if
                    if(present(WolfeConst1)) then
                        c1=WolfeConst1
                    else
                        c1=1d-4
                    end if
                    if(present(WolfeConst2)) then
                        c2=WolfeConst2
                    else
                        c2=0.9d0
                    end if
                if(present(f_fd)) then
                    i=f_fd(fnew,fdnew,x,dim)
                else
                    call f(fnew,x,dim)
                    call fd(fdnew,x,dim)
                end if
                if(freq>0) then
                    if(present(fdd)) then
                        i=fdd(H,x,dim)
                    else
                        i=djacobi(fd_j,dim,dim,H,x,1d-8)
                    end if
                    p=-fdnew
                    call My_dpotri(H,dim,i)
                    if(i==0) then
                        call syL2U(H,dim)
                        p=-matmul(H,fdnew)
                        phidnew=dot_product(fdnew,p)
                        a=1d0
                    end if
                end if
                if(freq<=0.or.i/=0) then!Hessian is either uncomputed or not positive definite, initial approximate inverse Hessian = a
                    p=-fdnew
                    phidnew=-dot_product(fdnew,fdnew)
                    if(-phidnew<tol) then
                        terminate=.true.
                        return
                    end if
                    if(fnew==0d0) then
                        a=1d0
                    else
                        a=-fnew/phidnew
                    end if
                    s=x
                    y=fdnew
                    if(present(f_fd)) then
                        call StrongWolfe_fdwithf(c1,c2,f,fd,f_fd,x,a,p,fnew,phidnew,fdnew,dim)
                    else
                        if(sw) then
                            call StrongWolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim)
                        else
                            call Wolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim)
                        end if
                    end if
                    phidnew=dot_product(fdnew,fdnew)
                    if(phidnew<tol) then
                        terminate=.true.
                        return
                    end if
                    if(dot_product(p,p)*a*a<tol) then
                        if(warn) then
                            write(*,'(1x,A84)')'BFGS warning: step length has converged, but gradient norm has not met accuracy goal'
                            write(*,'(1x,A56)')'A best estimation rather than exact solution is returned'
                            write(*,*)'Euclidean norm of gradient =',Sqrt(phidnew)
                        end if
                        terminate=.true.
                        return
                    end if
                    s=x-s
                    y=fdnew-y
                    rho=1d0/dot_product(y,s)
                    U=-rho*vector_direct_product(y,s,dim,dim)
                    forall(i=1:dim)
                        U(i,i)=U(i,i)+1d0
                    end forall
                    H=matmul(transpose(U),a*U)+rho*vector_direct_product(s,s,dim,dim)
                    p=-matmul(H,fdnew)
                    phidnew=dot_product(fdnew,p)
                    a=1d0
                end if
            end subroutine Initialize
            subroutine After()!Check convergence, update approximate inverse Hessian, determine new direction & step length
                !Check convergence
                phidnew=dot_product(fdnew,fdnew)
                if(phidnew<tol) then
                    terminate=.true.
                    return
                end if
                if(dot_product(p,p)*a*a<tol) then
                    if(warn) then
                        write(*,'(1x,A84)')'BFGS warning: step length has converged, but gradient norm has not met accuracy goal'
                        write(*,'(1x,A56)')'A best estimation rather than exact solution is returned'
                        write(*,*)'Euclidean norm of gradient =',Sqrt(phidnew)
                    end if
                    terminate=.true.
                    return
                end if
                !Determine new direction and step length, update approximate inverse Hessian
                i=mod(iIteration,freq)
                if(i==0) then!Every freq steps compute exact Hessian
                    i=fdd(U,x,dim)
                    call My_dpotri(U,dim,i)
                    if(i==0) then!Use exact Hessian if positive definite
                        call sycp(H,U,dim)
                        call syL2U(H,dim)
                        p=-matmul(H,fdnew)
                        phidnew=dot_product(fdnew,p)
                        a=1d0
                    end if
                end if
                if(i/=0) then!Exact Hessian is either uncomputed or not positive definite, update approximate Hessian
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
                end if
            end subroutine After
            subroutine After_NumericalHessian()!Check convergence, update approximate inverse Hessian, determine new direction & step length
                !Check convergence
                phidnew=dot_product(fdnew,fdnew)
                if(phidnew<tol) then
                    terminate=.true.
                    return
                end if
                if(dot_product(p,p)*a*a<tol) then
                    if(warn) then
                        write(*,'(1x,A84)')'BFGS warning: step length has converged, but gradient norm has not met accuracy goal'
                        write(*,'(1x,A56)')'A best estimation rather than exact solution is returned'
                        write(*,*)'Euclidean norm of gradient =',Sqrt(phidnew)
                    end if
                    terminate=.true.
                    return
                end if
                !Determine new direction and step length, update approximate inverse Hessian
                i=mod(iIteration,freq)
                if(i==0) then!Every freq steps compute exact Hessian
                    i=djacobi(fd_j,dim,dim,U,x,1d-8)
                    call My_dpotri(U,dim,i)
                    if(i==0) then!Use exact Hessian if positive definite
                        call sycp(H,U,dim)
                        call syL2U(H,dim)
                        p=-matmul(H,fdnew)
                        phidnew=dot_product(fdnew,p)
                        a=1d0
                    end if
                end if
                if(i/=0) then!Exact Hessian is either uncomputed or not positive definite, update approximate Hessian
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
                end if
            end subroutine After_NumericalHessian
            subroutine After_NoHessian()!Check convergence, update approximate inverse Hessian, determine new direction & step length
                !Check convergence
                phidnew=dot_product(fdnew,fdnew)
                if(phidnew<tol) then
                    terminate=.true.
                    return
                end if
                if(dot_product(p,p)*a*a<tol) then
                    if(warn) then
                        write(*,'(1x,A84)')'BFGS warning: step length has converged, but gradient norm has not met accuracy goal'
                        write(*,'(1x,A56)')'A best estimation rather than exact solution is returned'
                        write(*,*)'Euclidean norm of gradient =',Sqrt(phidnew)
                    end if
                    terminate=.true.
                    return
                end if
                !Determine new direction and step length, update approximate inverse Hessian
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
            end subroutine After_NoHessian
            subroutine fd_j(M,N,x,fdx)!Reformat for djacobi
                integer,intent(in)::M,N
                real*8,dimension(N),intent(in)::x
                real*8,dimension(M),intent(out)::fdx
                call fd(fdx,x,N)
            end subroutine fd_j
    end subroutine BFGS

    !Limited-memory Broyden–Fletcher–Goldfarb–Shanno (L-BFGS) quasi-Newton method, requiring Wolfe condition
    !Optional argument:
    !    Memory: (default = 10) memory usage = O( Memory * dim ). [3,30] is recommended
    subroutine LBFGS(f,fd,x,dim,Memory,f_fd,Strong,Warning,MaxIteration,Tolerance,WolfeConst1,WolfeConst2)
        !Required argument
            external::f,fd
            integer,intent(in)::dim
            real*8,dimension(dim),intent(inout)::x
        !Optional argument
            integer,external,optional::f_fd
            logical,intent(in),optional::Strong,Warning
            integer,intent(in),optional::Memory,MaxIteration
            real*8,intent(in),optional::Tolerance,WolfeConst1,WolfeConst2
        logical::sw,warn,terminate
        integer::M,maxit,iIteration,i,recent
        real*8::tol,c1,c2,a,fnew,phidnew
        real*8,dimension(dim)::p,fdnew,xold,fdold
        real*8,allocatable,dimension(:)::rho,alpha
        real*8,allocatable,dimension(:,:)::s,y
        call Initialize()
        if(terminate) return
        if(present(f_fd)) then!Cheaper to evaluate f' along with f
            do iIteration=1,maxit!Main loop
                call Before()!Before search
                call StrongWolfe_fdwithf(c1,c2,f,fd,f_fd,x,a,p,fnew,phidnew,fdnew,dim)!Line search
                call After()!After search
                if(terminate) return
            end do
        else
            if(sw) then!Use strong Wolfe condition instead of Wolfe condition
                do iIteration=1,maxit!Main loop
                    call Before()!Before search
                    call StrongWolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim)!Line search
                    call After()!After search
                    if(terminate) return
                end do
            else
                do iIteration=1,maxit!Main loop
                    call Before()!Before search
                    call Wolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim)!Line search
                    call After()!After search
                    if(terminate) return
                end do
            end if
        end if
        if(iIteration>maxit.and.warn) then
            write(*,'(1x,A38)')'Failed L-BFGS: max iteration exceeded!'
            write(*,*)'Euclidean norm of gradient =',Norm2(fdnew)
        end if
        !Clean up
            deallocate(rho)
            deallocate(alpha)
            deallocate(s)
            deallocate(y)
        contains
            subroutine Initialize()!Parameters, initial f(x) & f'(x), initial iteration history
                terminate=.false.
                !Set parameter according to optional argument
                    if(present(Memory)) then
                        M=max(0,Memory)
                    else
                        M=10
                    end if
                    if(present(Strong)) then
                        sw=Strong
                    else
                        sw=.true.
                    end if
                    if(present(Warning)) then
                        warn=Warning
                    else
                        warn=.true.
                    end if
                    if(present(MaxIteration)) then
                        maxit=MaxIteration
                    else
                        maxit=1000
                    end if
                    if(present(Tolerance)) then!To save sqrt cost, tolerance is squared
                        tol=Tolerance*Tolerance
                    else
                        tol=1d-30
                    end if
                    if(present(WolfeConst1)) then
                        c1=WolfeConst1
                    else
                        c1=1d-4
                    end if
                    if(present(WolfeConst2)) then
                        c2=WolfeConst2
                    else
                        c2=0.9d0
                    end if
                allocate(rho(0:M))
                allocate(alpha(0:M))
                allocate(s(dim,0:M))
                allocate(y(dim,0:M))
                if(present(f_fd)) then
                    i=f_fd(fnew,fdnew,x,dim)
                else
                    call f(fnew,x,dim)
                    call fd(fdnew,x,dim)
                end if
                p=-fdnew
                phidnew=-dot_product(fdnew,fdnew)
                if(-phidnew<tol) then
                    terminate=.true.
                    return
                end if
                if(fnew==0d0) then
                    a=1d0
                else
                    a=-fnew/phidnew
                end if
                xold=x
                fdold=fdnew
                !Initial approximate inverse Hessian = a
                if(present(f_fd)) then
                    call StrongWolfe_fdwithf(c1,c2,f,fd,f_fd,x,a,p,fnew,phidnew,fdnew,dim)
                else
                    if(sw) then
                        call StrongWolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim)
                    else
                        call Wolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim)
                    end if
                end if
                phidnew=dot_product(fdnew,fdnew)
                if(phidnew<tol) then
                    terminate=.true.
                    return
                end if
                if(dot_product(p,p)*a*a<tol) then
                    if(warn) then
                        write(*,'(1x,A84)')'BFGS warning: step length has converged, but gradient norm has not met accuracy goal'
                        write(*,'(1x,A56)')'A best estimation rather than exact solution is returned'
                        write(*,*)'Euclidean norm of gradient =',Sqrt(phidnew)
                    end if
                    terminate=.true.
                    return
                end if
                recent=0
                s(:,0)=x-xold
                y(:,0)=fdnew-fdold
                rho(0)=1d0/dot_product(y(:,0),s(:,0))
                do iIteration=1,M-1!Preiterate to get enough history
                    !Prepare
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
                    if(present(f_fd)) then!Line search
                        call StrongWolfe_fdwithf(c1,c2,f,fd,f_fd,x,a,p,fnew,phidnew,fdnew,dim)
                    else
                        if(sw) then
                            call StrongWolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim)
                        else
                            call Wolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim)
                        end if
                    end if
                    phidnew=dot_product(fdnew,fdnew)
                    if(phidnew<tol) then
                        terminate=.true.
                        return
                    end if
                    if(dot_product(p,p)*a*a<tol) then
                        if(warn) then
                            write(*,'(1x,A84)')'BFGS warning: step length has converged, but gradient norm has not met accuracy goal'
                            write(*,'(1x,A56)')'A best estimation rather than exact solution is returned'
                            write(*,*)'Euclidean norm of gradient =',Sqrt(phidnew)
                        end if
                        terminate=.true.
                        return
                    end if
                    recent=recent+1
                    s(:,recent)=x-xold
                    y(:,recent)=fdnew-fdold
                    rho(recent)=1d0/dot_product(y(:,recent),s(:,recent))
                end do
            end subroutine Initialize
            subroutine Before()!Prepare, determine new direction & step length
                xold=x!Prepare
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
            end subroutine Before
            subroutine After()!Check convergence, replace earliest iteration history with latest
                !Check convergence
                phidnew=dot_product(fdnew,fdnew)
                if(phidnew<tol) then
                    terminate=.true.
                    return
                end if
                if(dot_product(p,p)*a*a<tol) then
                    if(warn) then
                        write(*,'(1x,A86)')'L-BFGS warning: step length has converged, but gradient norm has not met accuracy goal'
                        write(*,'(1x,A56)')'A best estimation rather than exact solution is returned'
                        write(*,*)'Euclidean norm of gradient =',Sqrt(phidnew)
                    end if
                    terminate=.true.
                    return
                end if
                recent=mod(recent+1,M)
                s(:,recent)=x-xold
                y(:,recent)=fdnew-fdold
                rho(recent)=1d0/dot_product(y(:,recent),s(:,recent))
            end subroutine After
    end subroutine LBFGS

    !Conjugate gradient method, requiring either Wolfe or Strong Wolfe condition 
    !Optional argument:
    !    Method: (default = DY) which conjugate gradient method to use, currently available:
    !            DY (Dai-Yun), PR (Polak-Ribiere+)
    subroutine ConjugateGradient(f,fd,x,dim,Method,f_fd,Strong,Warning,MaxIteration,Tolerance,WolfeConst1,WolfeConst2)
        !Required argument
            external::f,fd
            integer,intent(in)::dim
            real*8,dimension(dim),intent(inout)::x
        !Optional argument
            character*2,intent(in),optional::Method
            integer,external,optional::f_fd
            logical,intent(in),optional::Strong,Warning
            integer,intent(in),optional::MaxIteration
            real*8,intent(in),optional::Tolerance,WolfeConst1,WolfeConst2
        logical::sw,warn,terminate
        character*2::type
        integer::maxit,iIteration,info
        real*8::tol,c1,c2,a,fnew,fold,phidnew,phidold
        real*8,dimension(dim)::p,fdnew,fdold
        call Initialize()
        if(terminate) return
        select case(type)
            case('DY')!Require Wolfe condition 
                if(present(f_fd)) then!Cheaper to evaluate f' along with f
                    do iIteration=1,maxit!Main loop
                        fold=fnew!Prepare
                        fdold=fdnew
                        phidold=phidnew
                        call StrongWolfe_fdwithf(c1,c2,f,fd,f_fd,x,a,p,fnew,phidnew,fdnew,dim)!Line search
                        call DY()!After search
                        if(terminate) return
                    end do
                else
                    if(sw) then!To meet Nocedal performance suggestion, Dai-Yun requires strong Wolfe condition
                        do iIteration=1,maxit!Main loop
                            fold=fnew!Prepare
                            fdold=fdnew
                            phidold=phidnew
                            call StrongWolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim)!Line search
                            call DY()!After search
                            if(terminate) return
                        end do
                    else
                        do iIteration=1,maxit!Main loop
                            fold=fnew!Prepare
                            fdold=fdnew
                            phidold=phidnew
                            call Wolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim)!Line search
                            call DY()!After search
                            if(terminate) return
                        end do
                    end if
                end if
            case('PR')!Require strong Wolfe condition 
                if(present(f_fd)) then!Cheaper to evaluate f' along with f
                    do iIteration=1,maxit!Main loop
                        fold=fnew!Prepare
                        fdold=fdnew
                        phidold=phidnew
                        call StrongWolfe_fdwithf(c1,c2,f,fd,f_fd,x,a,p,fnew,phidnew,fdnew,dim)!Line search
                        call PR()!After search
                        if(terminate) return
                    end do
                else
                    do iIteration=1,maxit!Main loop
                        fold=fnew!Prepare
                        fdold=fdnew
                        phidold=phidnew
                        call StrongWolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim)!Line search
                        call PR()!After search
                        if(terminate) return
                    end do
                end if
            case default!Throw a warning
                write(*,'(1x,A52,1x,A2)')'Program abort: unsupported conjugate gradient method',type
                stop
        end select
        if(iIteration>maxit.and.warn) then
            write(*,'(1x,A50)')'Failed conjugate gradient: max iteration exceeded!'
            write(*,*)'Euclidean norm of gradient =',Norm2(fdnew)
        end if
        contains
            subroutine Initialize()!Parameters, initial f(x) & f'(x), initial direction & step length
                terminate=.false.
                !Set parameter according to optional argument
                    if(present(Method)) then
                        type=Method
                    else
                        type='DY'
                    end if
                    if(present(Strong)) then
                        sw=Strong
                    else
                        sw=.true.
                    end if
                    if(present(Warning)) then
                        warn=Warning
                    else
                        warn=.true.
                    end if
                    if(present(MaxIteration)) then
                        maxit=MaxIteration
                    else
                        maxit=1000
                    end if
                    if(present(Tolerance)) then!To save sqrt cost, tolerance is squared
                        tol=Tolerance*Tolerance
                    else
                        tol=1d-30
                    end if
                    if(present(WolfeConst1)) then
                        c1=WolfeConst1
                    else
                        c1=1d-4
                    end if
                    if(present(WolfeConst2)) then
                        c2=WolfeConst2
                    else
                        c2=0.45d0
                    end if
                if(present(f_fd)) then
                    info=f_fd(fnew,fdnew,x,dim)
                else
                    call f(fnew,x,dim)
                    call fd(fdnew,x,dim)
                end if
                p=-fdnew
                phidnew=-dot_product(fdnew,fdnew)
                if(-phidnew<tol) then
                    terminate=.true.
                    return
                end if
                if(fnew==0d0) then
                    a=1d0
                else
                    a=-fnew/phidnew
                end if
            end subroutine Initialize
            subroutine DY()!Check convergence, determine new direction & step length
                !Check convergence
                phidnew=dot_product(fdnew,fdnew)
                if(phidnew<tol) then
                    terminate=.true.
                    return
                end if
                if(dot_product(p,p)*a*a<tol) then
                    if(warn) then
                        write(*,'(1x,A106)')'Dai-Yun conjugate gradient warning: step length has converged, but gradient norm has not met accuracy goal'
                        write(*,'(1x,A56)')'A best estimation rather than exact solution is returned'
                        write(*,*)'Euclidean norm of gradient =',Sqrt(phidnew)
                    end if
                    terminate=.true.
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
            end subroutine DY
            subroutine PR()!Check convergence, determine new direction & step length
                !Check convergence
                phidnew=dot_product(fdnew,fdnew)
                if(phidnew<tol) then
                    terminate=.true.
                    return
                end if
                if(dot_product(p,p)*a*a<tol) then
                    if(warn) then
                        write(*,'(1x,A113)')'Polak-Ribiere+ conjugate gradient warning: step length has converged, but gradient norm has not met accuracy goal'
                        write(*,'(1x,A56)')'A best estimation rather than exact solution is returned'
                        write(*,*)'Euclidean norm of gradient =',Sqrt(phidnew)
                    end if
                    terminate=.true.
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
            end subroutine PR
    end subroutine ConjugateGradient

    !=========== Line searcher ============
        !Some direction finders can only converge under strong Wolfe condition
        !Goldstein is only suitable for Newton, but not even for quasi-Newton
        !For Newton and quasi-Newton, the initial guess can always be a = 1, because their direction vector is well scaled
        !However, for methods whose direction is not determined by inverted (approximate) Hessian multiplying -gradient,
        !e.g., steepest descent and conjugate gradient, user has to come up with a good initial guess

        !Line search for a step length satisfying Wolfe condition
        !This routine is designed to minimize gradient computation
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
                    if(a<1d-37) then
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
                if(Abs(up-low)<1d-37) then
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
                        if(a<1d-37) return
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
                                if(a<1d-37) return
                            end do
                        end if
                    end if
                    if(a<1d-37) then
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
                if(Abs(up-low)<1d-37) return
            end do
        end subroutine StrongWolfeZoom
        
        !Almost same to StrongWolfe, but modified for cheap to evaluate f' along with f case
        subroutine StrongWolfe_fdwithf(c1,c2,f,fd,f_fd,x,a,p,fx,phid0,fdx,dim)
            real*8,intent(in)::c1,c2
            external::f,fd
            integer,external::f_fd
            integer,intent(in)::dim
            real*8,dimension(dim),intent(inout)::x
            real*8,intent(inout)::a
            real*8,dimension(dim),intent(in)::p
            real*8,intent(inout)::fx
            real*8,intent(in)::phid0
            real*8,dimension(dim),intent(out)::fdx
            real*8,parameter::increment=1.05d0! > 1
            integer::info
            real*8::c2_m_abs_phid0,fx0,atemp,aold,fold,phidnew,phidold
            real*8,dimension(dim)::x0
            !Initialize
                x0=x
                fx0=fx
                c2_m_abs_phid0=c2*Abs(phid0)
            !Check whether initial guess satisfies sufficient decrease condition
            x=x0+a*p
            info=f_fd(fx,fdx,x,dim)
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
                        info=f_fd(fx,fdx,x,dim)
                        phidnew=dot_product(fdx,p)
                        if(fx>=fold.or.phidnew<=0d0) then
                            x=x0
                            atemp=a
                            call StrongWolfeZoom_fdwithf(c1,c2,f_fd,x,a,p,fx0,phid0,aold,atemp,fold,fx,phidold,phidnew,fdx,dim)
                            fx=fx0
                            return
                        end if
                        if(a<1d-37) return
                    end do
                else!Search for larger a
                    do
                        aold=a
                        fold=fx
                        phidold=phidnew
                        a=aold*increment
                        x=x0+a*p
                        info=f_fd(fx,fdx,x,dim)
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
                                info=f_fd(fx,fdx,x,dim)
                                phidnew=dot_product(fdx,p)
                                if(fx>=fold.or.phidnew<=0d0) then
                                    x=x0
                                    atemp=a
                                    call StrongWolfeZoom_fdwithf(c1,c2,f_fd,x,a,p,fx0,phid0,aold,atemp,fold,fx,phidold,phidnew,fdx,dim)
                                    fx=fx0
                                    return
                                end if
                                if(a<1d-37) return
                            end do
                        end if
                    end if
                    if(a<1d-37) then
                        call fd(fdx,x,dim)
                        return
                    end if
                end do
            end if
        end subroutine StrongWolfe_fdwithf
        !Support StrongWolfe_fdwithf
        subroutine StrongWolfeZoom_fdwithf(c1,c2,f_fd,x,a,p,fx,phid0,low,up,flow,fup,phidlow,phidup,fdx,dim)
            real*8,intent(in)::c1,c2
            integer,external::f_fd
            integer,intent(in)::dim
            real*8,dimension(dim),intent(inout)::x
            real*8,intent(out)::a
            real*8,dimension(dim),intent(in)::p
            real*8,intent(inout)::fx
            real*8,intent(in)::phid0
            real*8,intent(inout)::low,up,flow,fup,phidlow,phidup
            real*8,dimension(dim),intent(out)::fdx
            integer::info
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
                info=f_fd(fx,fdx,x,dim)
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
                if(Abs(up-low)<1d-37) return
            end do
        end subroutine StrongWolfeZoom_fdwithf

        subroutine Goldstein()!Not implemented
            stop 'Program abort: Goldstein line search has not been implemented'
        end subroutine Goldstein
    !================ End =================
!------------------ End -------------------

!-------------- Trust region --------------
    !MKL trust-region nonlinear least square problem (trnlsp) solver wrapper
    !Solve f'(x) = 0 by minimizing merit function F(x) = f'(x)^2 through trust-region method
    !with model function m(p) = [ f'(x) + J(x) . p ]^2, where J(x) is the Jacobian
    !This procedure can be interpreted in different ways other than f'(x) = 0:
    !    f'(x) can be viewed as Sqrt(weight)*residual terms,
    !        then trnlsp minimizes the square penalty (this is where its name comes from)
    !    When M = N, f'(x) can also be considered as the gradient of f(x), Jacobian = Hessian,
    !        but trnlsp doesn't necessarily optimize f(x) (unless f(x) = const * F(x) by coincidence),
    !        it merely return a stationary point of f(x)
    !trnlspbc solves f'(x) = 0 subject to boundary condition low <= x <= up ( low < up )
    !External subroutine format:
    !    subroutine fd(f'(x),x,M,N)
    !    subroutine Jacobian(J(x),x,M,N)
    !IO format:
    !    M dimensional vector f'(x), N dimensional vectors x & low & up, M x N matrix J(x)
    !    On input x is an initial guess, on exit x is a solution of f'(x) = 0
    
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
                    if(djacobi(fd,N,M,J,x,1d-8)/=TR_SUCCESS) then
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
                    if(djacobi(fd,N,M,J,x,1d-8)/=TR_SUCCESS) then
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
    !Textbook is wrong: it claims Lagrangian multiplier method transforms constrained optimization into unconstrained one
    !However, L has no lower bound, approaching -infinity when multiplier diverges subject to violated constraint
    !Lagrangian multiplier method actually turns a minimization problem into a saddle point problem,
    !which cannot necessarily be solved through try decreasing f(x), making all unconstrained minimizers fail
    !Lagrangian multiplier method is only good when L has unique saddle point: we may simply minimize || L'(x) ||
    !In general case, we have to turn to the augmented Lagrangian method
    !Nomenclature:
    !    f = the target function to be minimized
    !    c = the constraint in standard form: c(x) = 0 for equality, c(x) >= 0 for inequality
    !External subroutine format:
    !    subroutine f(f(x),x,dim)
    !    subroutine fd(f'(x),x,dim)
    !    subroutine f_fd(f(x),f'(x),x,dim)
    !    subroutine fdd(f''(x),x,dim)
    !    subroutine c(c(x),x,dim)
    !IO format:
    !    dim dimensional vectors x & f'(x), dim order matrix f''(x)
    !    On input x is an initial guess, on exit x is a local minimum of f(x)

    subroutine AugmentedLagrangian(f)
        external::f
    end subroutine AugmentedLagrangian

    subroutine AugmentedLagrangianBC()!Not implemented
    end subroutine AugmentedLagrangianBC
!------------------ End -------------------

!--------------- Heuristic ----------------
    !Not implemented
!------------------ End -------------------

end module NonlinearOptimization