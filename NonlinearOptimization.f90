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
    use LinearAlgebra!For inverting exact Hessian
    use mkl_rci_type !For numerical exact Hessian & Trust region solver
    use mkl_rci      !For numerical exact Hessian & Trust region solver
    implicit none

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
    !    subroutine f & fd, dim dimensional vector x, integer dim
    !Common optional argument:
    !    f_fd: presence means evaluating f(x) & f'(x) together is cheaper than separately,
    !          strong Wolfe condition is thus applied, because Wolfe does not benefit from this
    !    Strong: (default = true) if true, use strong Wolfe condition instead of Wolfe condition
    !    Warning: (default = true) if false, all warnings will be suppressed
    !    MaxIteration: (default = 1000) max number of iterations to perform
    !    Precision: (default = 1d-15) convergence considered when || f'(x) ||_2 < Precision
    !    MinStepLength: (default = 1d-15) terminate if search step < MinStepLength before || f'(x) ||_2 converges
    !    WolfeConst1 & WolfeConst2: 0 < WolfeConst1 < WolfeConst2 <  1  for Newton & quasi-Newton
    !                               0 < WolfeConst1 < WolfeConst2 < 0.5 for conjugate gradient
    !On input x is an initial guess, on exit x is a local minimum of f(x)

    !Newton-Raphson method, requiring Wolfe condition
    !Optional argument: fdd: presence means analytical Hessian is available, otherwise call djacobi for central difference Hessian
    subroutine NewtonRaphson(f,fd,x,dim,fdd,f_fd,Strong,Warning,MaxIteration,Precision,MinStepLength,WolfeConst1,WolfeConst2)
        !Required argument
            external::f,fd
            integer,intent(in)::dim
            real*8,dimension(dim),intent(inout)::x
        !Optional argument
            integer,external,optional::fdd,f_fd
            logical,intent(in),optional::Strong,Warning
            integer,intent(in),optional::MaxIteration
            real*8,intent(in),optional::Precision,MinStepLength,WolfeConst1,WolfeConst2
        logical::sw,warn,terminate
        integer::maxit,iIteration,info
        real*8::tol,minstep,c1,c2,a,fnew,phidnew,phidold
        real*8,dimension(dim)::p,fdnew
        real*8,dimension(dim,dim)::Hessian
        !Initialize
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
                if(present(Precision)) then!To save sqrt cost, precision is squared
                    tol=Precision*Precision
                else
                    tol=1d-30
                end if
                if(present(MinStepLength)) then!To save sqrt cost, MinStepLength is squared
                    minstep=MinStepLength*MinStepLength
                else
                    minstep=1d-30
                end if
                if(present(WolfeConst1)) then
                    c1=max(1d-15,WolfeConst1)!Fail safe
                else
                    c1=1d-4
                end if
                if(present(WolfeConst2)) then
                    c2=min(1d0-1d-15,max(c1+1d-15,WolfeConst2))!Fail safe
                else
                    c2=0.9d0
                end if
            if(present(f_fd)) then!Initial f(x) & f'(x)
                info=f_fd(fnew,fdnew,x,dim)
            else
                call f(fnew,x,dim)
                call fd(fdnew,x,dim)
            end if
            !Initial direction & step length
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
                if(-phidnew<tol) return
                if(fnew==0d0) then
                    a=1d0
                else
                    a=-fnew/phidnew
                end if
            end if
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
            subroutine After()!Check convergence, determine new direction & step length
                !Check convergence
                phidnew=dot_product(fdnew,fdnew)
                if(phidnew<tol) then
                    terminate=.true.
                    return
                end if
                if(dot_product(p,p)*a*a<minstep) then
                    if(warn) then
                        write(*,'(1x,A94)')'Newton-Raphson warning: step length has converged, but gradient norm has not met accuracy goal'
                        write(*,'(1x,A56)')'A best estimation rather than exact solution is returned'
                        write(*,*)'Euclidean norm of gradient =',dSqrt(phidnew)
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
                if(dot_product(p,p)*a*a<minstep) then
                    if(warn) then
                        write(*,'(1x,A94)')'Newton-Raphson warning: step length has converged, but gradient norm has not met accuracy goal'
                        write(*,'(1x,A56)')'A best estimation rather than exact solution is returned'
                        write(*,*)'Euclidean norm of gradient =',dSqrt(phidnew)
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
    subroutine BFGS(f,fd,x,dim,fdd,ExactStep,f_fd,Strong,Warning,MaxIteration,Precision,MinStepLength,WolfeConst1,WolfeConst2)
        !Required argument
            external::f,fd
            integer,intent(in)::dim
            real*8,dimension(dim),intent(inout)::x
        !Optional argument
            integer,external,optional::fdd,f_fd
            logical,intent(in),optional::Strong,Warning
            integer,intent(in),optional::ExactStep,MaxIteration
            real*8,intent(in),optional::Precision,MinStepLength,WolfeConst1,WolfeConst2
        logical::sw,warn,terminate
        integer::freq,maxit,iIteration,i
        real*8::tol,minstep,c1,c2,a,fnew,phidnew,rho
        real*8,dimension(dim)::p,fdnew,s,y
        real*8,dimension(dim,dim)::U,H!Approximate inverse Hessian
        !Initialize
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
                if(present(Precision)) then!To save sqrt cost, precision is squared
                    tol=Precision*Precision
                else
                    tol=1d-30
                end if
                if(present(MinStepLength)) then!To save sqrt cost, MinStepLength is squared
                    minstep=MinStepLength*MinStepLength
                else
                    minstep=1d-30
                end if
                if(present(WolfeConst1)) then
                    c1=max(1d-15,WolfeConst1)!Fail safe
                else
                    c1=1d-4
                end if
                if(present(WolfeConst2)) then
                    c2=min(1d0-1d-15,max(c1+1d-15,WolfeConst2))!Fail safe
                else
                    c2=0.9d0
                end if
            if(present(f_fd)) then!Initial f(x) & f'(x)
                i=f_fd(fnew,fdnew,x,dim)
            else
                call f(fnew,x,dim)
                call fd(fdnew,x,dim)
            end if
            !Initial approximate inverse Hessian & direction & step length
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
                if(-phidnew<tol) return
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
                if(phidnew<tol) return
                if(dot_product(p,p)*a*a<minstep) then
                    if(warn) then
                        write(*,'(1x,A84)')'BFGS warning: step length has converged, but gradient norm has not met accuracy goal'
                        write(*,'(1x,A56)')'A best estimation rather than exact solution is returned'
                        write(*,*)'Euclidean norm of gradient =',dSqrt(phidnew)
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
                H=matmul(transpose(U),a*U)+rho*vector_direct_product(s,s,dim,dim)
                p=-matmul(H,fdnew)
                phidnew=dot_product(fdnew,p)
                a=1d0
            end if
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
            subroutine After()!Check convergence, update approximate inverse Hessian, determine new direction & step length
                !Check convergence
                phidnew=dot_product(fdnew,fdnew)
                if(phidnew<tol) then
                    terminate=.true.
                    return
                end if
                if(dot_product(p,p)*a*a<minstep) then
                    if(warn) then
                        write(*,'(1x,A84)')'BFGS warning: step length has converged, but gradient norm has not met accuracy goal'
                        write(*,'(1x,A56)')'A best estimation rather than exact solution is returned'
                        write(*,*)'Euclidean norm of gradient =',dSqrt(phidnew)
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
                if(dot_product(p,p)*a*a<minstep) then
                    if(warn) then
                        write(*,'(1x,A84)')'BFGS warning: step length has converged, but gradient norm has not met accuracy goal'
                        write(*,'(1x,A56)')'A best estimation rather than exact solution is returned'
                        write(*,*)'Euclidean norm of gradient =',dSqrt(phidnew)
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
                if(dot_product(p,p)*a*a<minstep) then
                    if(warn) then
                        write(*,'(1x,A84)')'BFGS warning: step length has converged, but gradient norm has not met accuracy goal'
                        write(*,'(1x,A56)')'A best estimation rather than exact solution is returned'
                        write(*,*)'Euclidean norm of gradient =',dSqrt(phidnew)
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
    !Optional argument: Memory: (default = 10) memory usage = O( Memory * dim ). [3,30] is recommended (must > 0)
    subroutine LBFGS(f,fd,x,dim,Memory,f_fd,Strong,Warning,MaxIteration,Precision,MinStepLength,WolfeConst1,WolfeConst2)
        !Required argument
            external::f,fd
            integer,intent(in)::dim
            real*8,dimension(dim),intent(inout)::x
        !Optional argument
            integer,external,optional::f_fd
            logical,intent(in),optional::Strong,Warning
            integer,intent(in),optional::Memory,MaxIteration
            real*8,intent(in),optional::Precision,MinStepLength,WolfeConst1,WolfeConst2
        logical::sw,warn,terminate
        integer::M,maxit,iIteration,i,recent
        real*8::tol,minstep,c1,c2,a,fnew,phidnew
        real*8,dimension(dim)::p,fdnew,xold,fdold
        real*8,allocatable,dimension(:)::rho,alpha
        real*8,allocatable,dimension(:,:)::s,y
        !Initialize
            terminate=.false.
            !Set parameter according to optional argument
                if(present(Memory)) then
                    M=max(1,Memory)!Fail safe
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
                if(present(Precision)) then!To save sqrt cost, precision is squared
                    tol=Precision*Precision
                else
                    tol=1d-30
                end if
                if(present(MinStepLength)) then!To save sqrt cost, MinStepLength is squared
                    minstep=MinStepLength*MinStepLength
                else
                    minstep=1d-30
                end if
                if(present(WolfeConst1)) then
                    c1=max(1d-15,WolfeConst1)!Fail safe
                else
                    c1=1d-4
                end if
                if(present(WolfeConst2)) then
                    c2=min(1d0-1d-15,max(c1+1d-15,WolfeConst2))!Fail safe
                else
                    c2=0.9d0
                end if
            allocate(rho(0:M))
            allocate(alpha(0:M))
            allocate(s(dim,0:M))
            allocate(y(dim,0:M))
            if(present(f_fd)) then!Initial f(x) & f'(x)
                i=f_fd(fnew,fdnew,x,dim)
            else
                call f(fnew,x,dim)
                call fd(fdnew,x,dim)
            end if
            !Initial iteration history
            p=-fdnew
            phidnew=-dot_product(fdnew,fdnew)
            if(-phidnew<tol) return
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
            if(phidnew<tol) return
            if(dot_product(p,p)*a*a<minstep) then
                if(warn) then
                    write(*,'(1x,A84)')'BFGS warning: step length has converged, but gradient norm has not met accuracy goal'
                    write(*,'(1x,A56)')'A best estimation rather than exact solution is returned'
                    write(*,*)'Euclidean norm of gradient =',dSqrt(phidnew)
                end if
                return
            end if
            recent=0
            s(:,0)=x-xold
            y(:,0)=fdnew-fdold
            rho(0)=1d0/dot_product(y(:,0),s(:,0))
            do iIteration=1,M-1!Preiterate to get enough history
                xold=x!Prepare
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
                if(phidnew<tol) return
                if(dot_product(p,p)*a*a<minstep) then
                    if(warn) then
                        write(*,'(1x,A84)')'BFGS warning: step length has converged, but gradient norm has not met accuracy goal'
                        write(*,'(1x,A56)')'A best estimation rather than exact solution is returned'
                        write(*,*)'Euclidean norm of gradient =',dSqrt(phidnew)
                    end if
                    return
                end if
                recent=recent+1
                s(:,recent)=x-xold
                y(:,recent)=fdnew-fdold
                rho(recent)=1d0/dot_product(y(:,recent),s(:,recent))
            end do
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
                if(dot_product(p,p)*a*a<minstep) then
                    if(warn) then
                        write(*,'(1x,A86)')'L-BFGS warning: step length has converged, but gradient norm has not met accuracy goal'
                        write(*,'(1x,A56)')'A best estimation rather than exact solution is returned'
                        write(*,*)'Euclidean norm of gradient =',dSqrt(phidnew)
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
    !Available methods: DY (Dai-Yun), PR (Polak-Ribiere+)
    !Optional argument: Method: (default = DY) which conjugate gradient method to use
    subroutine ConjugateGradient(f,fd,x,dim,Method,f_fd,Strong,Warning,MaxIteration,Precision,MinStepLength,WolfeConst1,WolfeConst2)
        !Required argument
            external::f,fd
            integer,intent(in)::dim
            real*8,dimension(dim),intent(inout)::x
        !Optional argument
            character*2,intent(in),optional::Method
            integer,external,optional::f_fd
            logical,intent(in),optional::Strong,Warning
            integer,intent(in),optional::MaxIteration
            real*8,intent(in),optional::Precision,MinStepLength,WolfeConst1,WolfeConst2
        logical::sw,warn,terminate
        character*2::type
        integer::maxit,iIteration,info
        real*8::tol,minstep,c1,c2,a,fnew,fold,phidnew,phidold
        real*8,dimension(dim)::p,fdnew,fdold
        !Initialize
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
                if(present(Precision)) then!To save sqrt cost, precision is squared
                    tol=Precision*Precision
                else
                    tol=1d-30
                end if
                if(present(MinStepLength)) then!To save sqrt cost, MinStepLength is squared
                    minstep=MinStepLength*MinStepLength
                else
                    minstep=1d-30
                end if
                if(present(WolfeConst1)) then
                    c1=max(1d-15,WolfeConst1)!Fail safe
                else
                    c1=1d-4
                end if
                if(present(WolfeConst2)) then
                    c2=min(1d0-1d-15,max(c1+1d-15,WolfeConst2))!Fail safe
                else
                    c2=0.45d0
                end if
            if(present(f_fd)) then!Initial f(x) & f'(x)
                info=f_fd(fnew,fdnew,x,dim)
            else
                call f(fnew,x,dim)
                call fd(fdnew,x,dim)
            end if
            !Initial direction & step length
            p=-fdnew
            phidnew=-dot_product(fdnew,fdnew)
            if(-phidnew<tol) return
            if(fnew==0d0) then
                a=1d0
            else
                a=-fnew/phidnew
            end if
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
            subroutine DY()!Check convergence, determine new direction & step length
                !Check convergence
                phidnew=dot_product(fdnew,fdnew)
                if(phidnew<tol) then
                    terminate=.true.
                    return
                end if
                if(dot_product(p,p)*a*a<minstep) then
                    if(warn) then
                        write(*,'(1x,A106)')'Dai-Yun conjugate gradient warning: step length has converged, but gradient norm has not met accuracy goal'
                        write(*,'(1x,A56)')'A best estimation rather than exact solution is returned'
                        write(*,*)'Euclidean norm of gradient =',dSqrt(phidnew)
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
                if(dot_product(p,p)*a*a<minstep) then
                    if(warn) then
                        write(*,'(1x,A113)')'Polak-Ribiere+ conjugate gradient warning: step length has converged, but gradient norm has not met accuracy goal'
                        write(*,'(1x,A56)')'A best estimation rather than exact solution is returned'
                        write(*,*)'Euclidean norm of gradient =',dSqrt(phidnew)
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
        !Goldstein is only suitable for Newton, but not even for quasi-Newton (so it has not been implemented)
        !For Newton and quasi-Newton, the initial guess can always be a = 1, because their direction vector is well scaled
        !However, for methods whose direction is not determined by inverted (approximate) Hessian multiplying -gradient,
        !e.g., steepest descent and conjugate gradient, user has to come up with a good initial guess
        !Input:  c1 & c2 are Wolfe constants (0<c1<c2<1), x is current x
        !        a is initial guess of a, p is current p, fx = f(x), phid0 = phi'(0)
        !Output: a harvests the step length satisfying certain condition, x = x + a * p, fx = f(x), fdx = f'(x)
        !Optional argument: Increment: (default = 1.05) each iteration change step length by how much time (must > 1)

        !Line search for a step length satisfying Wolfe condition
        !This routine is designed to minimize gradient computation
        subroutine Wolfe(c1,c2,f,fd,x,a,p,fx,phid0,fdx,dim,Increment)
            !Required argument
                real*8,intent(in)::c1,c2
                external::f,fd
                integer,intent(in)::dim
                real*8,dimension(dim),intent(inout)::x
                real*8,intent(inout)::a
                real*8,dimension(dim),intent(in)::p
                real*8,intent(inout)::fx
                real*8,intent(in)::phid0
                real*8,dimension(dim),intent(out)::fdx
            !Optional argument
                real*8,intent(in),optional::Increment
            real*8::incrmt,c2_m_abs_phid0,fx0,ftemp,atemp,aold,fold,phidx
            real*8,dimension(dim)::x0
            !Initialize
                if(present(Increment)) then
                    incrmt=max(1d0+1d-15,Increment)!Fail safe
                else
                    incrmt=1.05d0
                end if
                x0=x
                fx0=fx
                c2_m_abs_phid0=c2*Abs(phid0)
            !Check whether initial guess satisfies sufficient decrease condition
            x=x0+a*p
            call f(fx,x,dim)
            if(fx<=fx0+c1*a*phid0) then!Satisfied, search for larger a
                do
                    aold=a
                    fold=fx
                    a=aold*incrmt
                    x=x0+a*p
                    call f(fx,x,dim)
                    if(fx>fx0+c1*a*phid0) then
                        x=x0+aold*p
                        call fd(fdx,x,dim)
                        phidx=dot_product(fdx,p)
                        if(phidx>c2_m_abs_phid0) then
                            a=aold
                            fx=fold
                        else
                            atemp=a
                            ftemp=fx
                            call WolfeZoom(aold,atemp,fold,ftemp,phidx)
                        end if
                        return
                    end if
                end do
            else!Violated, search for smaller a
                do
                    aold=a
                    fold=fx
                    a=aold/incrmt
                    x=x0+a*p
                    call f(fx,x,dim)
                    if(fx<=fx0+c1*a*phid0) then
                        call fd(fdx,x,dim)
                        phidx=dot_product(fdx,p)
                        if(phidx<c2_m_abs_phid0) then
                            atemp=a
                            ftemp=fx
                            call WolfeZoom(atemp,aold,ftemp,fold,phidx)
                        end if
                        return
                    end if
                    if(a<1d-15) then
                        call fd(fdx,x,dim)
                        return
                    end if
                end do
            end if
            contains
                !Support Wolfe. low and up must satisfy:
                !    low < up
                !    low satisfies sufficient decrease condition, but up violates
                !    phi'(low) < 0
                subroutine WolfeZoom(low,up,flow,fup,phidlow)!Restore x and fx to input status before calling zoom
                    real*8,intent(inout)::low,up,flow,fup,phidlow
                    real*8::phidnew,phidlow_m_a
                    phidlow_m_a=phidlow*a!Initialize
                    do
                        !Updata a by quadratic interpolation
                            a=phidlow_m_a*a/2d0/(flow+phidlow_m_a-fup)
                            if(.not.(a>low.and.a<up)) a=(low+up)/2d0!Fail safe
                        x=x0+a*p
                        call f(fx,x,dim)
                        if(fx>fx0+c1*a*phid0) then
                            up=a
                            if(Abs(up-low)<1d-15) then
                                call fd(fdx,x,dim)
                                return
                            end if
                            fup=fx
                        else
                            call fd(fdx,x,dim)
                            phidnew=dot_product(fdx,p)
                            if(phidnew>c2_m_abs_phid0) return
                            low=a
                            if(Abs(up-low)<1d-15) return
                            flow=fx
                            phidlow=phidnew
                            phidlow_m_a=phidlow*a
                        end if
                    end do
                end subroutine WolfeZoom
        end subroutine Wolfe
        
        !Line search for a step length satisfying strong Wolfe condition
        subroutine StrongWolfe(c1,c2,f,fd,x,a,p,fx,phid0,fdx,dim,Increment)
            !Required argument
                real*8,intent(in)::c1,c2
                external::f,fd
                integer,intent(in)::dim
                real*8,dimension(dim),intent(inout)::x
                real*8,intent(inout)::a
                real*8,dimension(dim),intent(in)::p
                real*8,intent(inout)::fx
                real*8,intent(in)::phid0
                real*8,dimension(dim),intent(out)::fdx
            !Optional argument
                real*8,intent(in),optional::Increment
            real*8::incrmt,c2_m_abs_phid0,fx0,ftemp,atemp,aold,fold,phidnew,phidold
            real*8,dimension(dim)::x0
            !Initialize
                if(present(Increment)) then
                    incrmt=max(1d0+1d-15,Increment)!Fail safe
                else
                    incrmt=1.05d0
                end if
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
                        a=aold/incrmt
                        x=x0+a*p
                        call f(fx,x,dim)
                        call fd(fdx,x,dim)
                        phidnew=dot_product(fdx,p)
                        if(fx>=fold.or.phidnew<=0d0) then
                            atemp=a
                            ftemp=fx
                            call StrongWolfeZoom(aold,atemp,fold,ftemp,phidold,phidnew)
                            return
                        end if
                        if(a<1d-15) return
                    end do
                else!Search for larger a
                    do
                        aold=a
                        fold=fx
                        phidold=phidnew
                        a=aold*incrmt
                        x=x0+a*p
                        call f(fx,x,dim)
                        call fd(fdx,x,dim)
                        phidnew=dot_product(fdx,p)
                        if(fx>fx0+c1*a*phid0.or.fx>=fold) then
                            atemp=a
                            ftemp=fx
                            call StrongWolfeZoom(aold,atemp,fold,ftemp,phidold,phidnew)
                            return
                        end if
                        if(phidnew>0d0) then
                            if(Abs(phidnew)<=c2_m_abs_phid0) then
                                return
                            else
                                atemp=a
                                ftemp=fx
                                call StrongWolfeZoom(atemp,aold,ftemp,fold,phidnew,phidold)
                                fx=fx0
                            end if
                        end if
                    end do
                end if
            else!Violated, 1st search for smaller a satisfying sufficient decrease condition
                do
                    aold=a
                    fold=fx
                    a=aold/incrmt
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
                            atemp=a
                            ftemp=fx
                            call StrongWolfeZoom(atemp,aold,ftemp,fold,phidnew,phidold)
                            return
                        else!Search for such an a that phi(a) < phi(aold) & phid(a) > 0 is false
                            do
                                aold=a
                                fold=fx
                                phidold=phidnew
                                a=aold/incrmt
                                x=x0+a*p
                                call f(fx,x,dim)
                                call fd(fdx,x,dim)
                                phidnew=dot_product(fdx,p)
                                if(fx>=fold.or.phidnew<=0d0) then
                                    atemp=a
                                    ftemp=fx
                                    call StrongWolfeZoom(aold,atemp,fold,ftemp,phidold,phidnew)
                                    return
                                end if
                                if(a<1d-15) return
                            end do
                        end if
                    end if
                    if(a<1d-15) then
                        call fd(fdx,x,dim)
                        return
                    end if
                end do
            end if
            contains
                !Support StrongWolfe. low and up must satisfy:
                !    low satisfies sufficient decrease condition
                !    ( up - low ) * phi'(low) < 0
                !    [ low, up ] (or [ up, low ]) contains a step length satisfying strong Wolfe condition
                !        This means at least 1 of 3 following statements is true:
                !            up violates the sufficient decrease condition
                !            phi(up) >= phi(low)
                !            up < low & phi'(up) <= 0
                subroutine StrongWolfeZoom(low,up,flow,fup,phidlow,phidup)
                    real*8,intent(inout)::low,up,flow,fup,phidlow,phidup
                    real*8::phidnew,d1,d2
                    do while(Abs(up-low)>1d-15)
                        !Updata a by cubic interpolation
                            d1=phidlow+phidup-3d0*(flow-fup)/(low-up)
                            d2=up-low
                            if(d2>0d0) then
                                d2=dSqrt(d1*d1-phidlow*phidup)
                            else
                                d2=-dSqrt(d1*d1-phidlow*phidup)
                            end if
                            a=up-(up-low)*(phidup+d2-d1)/(phidup-phidlow+2d0*d2)
                            if(.not.(a>min(low,up).and.a<max(low,up))) a=(low+up)/2d0!Fail safe
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
                    end do
                end subroutine StrongWolfeZoom
        end subroutine StrongWolfe
        !When it is cheaper to evaluate f' along with f
        subroutine StrongWolfe_fdwithf(c1,c2,f,fd,f_fd,x,a,p,fx,phid0,fdx,dim,Increment)
            !Required argument
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
            !Optional argument
                real*8,intent(in),optional::Increment
            integer::info
            real*8::incrmt,c2_m_abs_phid0,fx0,ftemp,atemp,aold,fold,phidnew,phidold
            real*8,dimension(dim)::x0
            !Initialize
                if(present(Increment)) then
                    incrmt=max(1d0+1d-15,Increment)!Fail safe
                else
                    incrmt=1.05d0
                end if
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
                        a=aold/incrmt
                        x=x0+a*p
                        info=f_fd(fx,fdx,x,dim)
                        phidnew=dot_product(fdx,p)
                        if(fx>=fold.or.phidnew<=0d0) then
                            atemp=a
                            ftemp=fx
                            call StrongWolfeZoom_fdwithf(aold,atemp,fold,ftemp,phidold,phidnew)
                            return
                        end if
                        if(a<1d-15) return
                    end do
                else!Search for larger a
                    do
                        aold=a
                        fold=fx
                        phidold=phidnew
                        a=aold*incrmt
                        x=x0+a*p
                        info=f_fd(fx,fdx,x,dim)
                        phidnew=dot_product(fdx,p)
                        if(fx>fx0+c1*a*phid0.or.fx>=fold) then
                            atemp=a
                            ftemp=fx
                            call StrongWolfeZoom_fdwithf(aold,atemp,fold,ftemp,phidold,phidnew)
                            return
                        end if
                        if(phidnew>0d0) then
                            if(Abs(phidnew)<=c2_m_abs_phid0) return
                            atemp=a
                            ftemp=fx
                            call StrongWolfeZoom_fdwithf(atemp,aold,ftemp,fold,phidnew,phidold)
                            return
                        end if
                    end do
                end if
            else!Violated, 1st search for smaller a satisfying sufficient decrease condition
                do
                    aold=a
                    fold=fx
                    a=aold/incrmt
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
                            atemp=a
                            ftemp=fx
                            call StrongWolfeZoom_fdwithf(atemp,aold,ftemp,fold,phidnew,phidold)
                            return
                        else!Search for such an a that phi(a) < phi(aold) & phid(a) > 0 is false
                            do
                                aold=a
                                fold=fx
                                phidold=phidnew
                                a=aold/incrmt
                                x=x0+a*p
                                info=f_fd(fx,fdx,x,dim)
                                phidnew=dot_product(fdx,p)
                                if(fx>=fold.or.phidnew<=0d0) then
                                    atemp=a
                                    ftemp=fx
                                    call StrongWolfeZoom_fdwithf(aold,atemp,fold,ftemp,phidold,phidnew)
                                    return
                                end if
                                if(a<1d-15) return
                            end do
                        end if
                    end if
                    if(a<1d-15) then
                        call fd(fdx,x,dim)
                        return
                    end if
                end do
            end if
            contains
                subroutine StrongWolfeZoom_fdwithf(low,up,flow,fup,phidlow,phidup)
                    real*8,intent(inout)::low,up,flow,fup,phidlow,phidup
                    real*8::phidnew,d1,d2
                    do while(Abs(up-low)>1d-15)
                        !Updata a by cubic interpolation
                            d1=phidlow+phidup-3d0*(flow-fup)/(low-up)
                            d2=up-low
                            if(d2>0d0) then
                                d2=dSqrt(d1*d1-phidlow*phidup)
                            else
                                d2=-dSqrt(d1*d1-phidlow*phidup)
                            end if
                            a=up-(up-low)*(phidup+d2-d1)/(phidup-phidlow+2d0*d2)
                            if(.not.(a>min(low,up).and.a<max(low,up))) a=(low+up)/2d0!Fail safe
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
                    end do
                end subroutine StrongWolfeZoom_fdwithf
        end subroutine StrongWolfe_fdwithf
    !================ End =================
!------------------ End -------------------

!-------------- Trust region --------------
    !MKL trust-region nonlinear least square problem (trnlsp) solver wrapper
    !Solve f'(x) = 0 by minimizing merit function F(x) = f'(x)^2 through trust-region method
    !with model function m(p) = [ f'(x) + J(x) . p ]^2, where J(x) is the Jacobian
    !This algorithm can be interpreted in different ways other than f'(x) = 0:
    !    f'(x) can be viewed as Sqrt(weight)*residual terms,
    !        then trnlsp minimizes the square penalty (this is where its name comes from)
    !    When M = N, f'(x) can also be considered as the gradient of f(x), Jacobian = Hessian,
    !        but trnlsp doesn't necessarily optimize f(x) (unless f(x) = const * F(x) by coincidence),
    !        it merely return a stationary point of f(x)
    !External procedure format:
    !    subroutine fd(f'(x),x,M,N)
    !    integer function Jacobian(J(x),x,M,N)
    !    M dimensional vector f'(x), N dimensional vector x, M x N matrix J(x)
    !Required argument:
    !    subroutine fd, N dimensional vector x, integer M & N
    !Optional argument:
    !    Jacobian: presence means analytical Jacobian is available, otherwise call djacobi for central difference Jacobian
    !    low & up: presence means x subjects to boundary condition low <= x <= up (must both present & low < up)
    !    Warning: (default = true) if false, all warnings will be suppressed
    !    MaxIteration: (default = 1000) max number of iterations to perform
    !    MaxStepIteration: (default = 100) max number of iterations for determining each step length
    !    Precision: (default = 1d-15) convergence considered when || f'(x) ||_2 < Precision
    !    MinStepLength: (default = 1d-15) terminate if search step < MinStepLength before || f'(x) ||_2 converges
    !On input x is an initial guess, on exit x is a solution of f'(x) = 0
    subroutine TrustRegion(fd,x,M,N,Jacobian,low,up,Warning,MaxIteration,MaxStepIteration,Precision,MinStepLength)
        !Required argument
            external::fd
            integer,intent(in)::M,N
            real*8,dimension(N),intent(inout)::x
        !Optional argument
            integer,external,optional::Jacobian
            real*8,dimension(N),intent(in),optional::low,up
            logical,intent(in),optional::Warning
            integer,intent(in),optional::MaxIteration,MaxStepIteration
            real*8 ,intent(in),optional::Precision,MinStepLength
        !Reverse communication interface (RCI)
            integer::RCI_request!Recieve job request
            integer,dimension(6)::info!Results of input parameter checking
            type(handle_tr)::handle!Trust-region solver handle
        !Job control
            logical::warn
            integer::maxit=1000,maxstepit=1000
            !tol(1:5) contains the stopping criteria for solving f'(x) = 0:
            !    1, trust region radius < tol(1)
            !    2, || f'(x) ||_2 < tol(2)
            !    3, || Jacobian ||_1 < tol(3) 
            !    4, || s ||_2 < tol(4), where s is the trial step
            !    5, || f'(x) ||_2 - || f'(x) - Jacobian . s ||_2 < tol(5)
            !tol(6) is the precision of s calculation
            real*8,dimension(6)::tol
        !TotalIteration harvests the solver stops after how many interations
        !StopReason harvests why the solver has stopped:
        !    1,   max iteration exceeded
        !    2-6, tol(StopReason-1) is met
        integer::i,TotalIteration,StopReason
        real*8::InitialResidual,FinalResidual,StepBound
        real*8,dimension(M)::fdx
        real*8,dimension(M,N)::J
        !Initialize
            !Set parameter according to optional argument
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
                if(present(MaxStepIteration)) then
                    maxstepit=MaxStepIteration
                else
                    maxstepit=100
                end if
                tol=[1d-15,1d-15,1d-15,1d-15,1d-15,1d-15]
                if(present(Precision)) then
                    tol(2)=Precision
                    tol(5)=Precision
                end if
                if(present(MinStepLength)) then
                    tol(1)=MinStepLength
                    tol(4)=MinStepLength
                end if
            fdx=0d0
            J=0d0
            StepBound=100d0
            RCI_request=0
        if(present(low).and.present(up)) then
            if(dtrnlspbc_init(handle,N,M,x,low,up,tol,maxit,maxstepit,StepBound)/=TR_SUCCESS) then
                call mkl_free_buffers
                return
            end if
            if(dtrnlspbc_check(handle,N,M,J,fdx,low,up,tol,info)/=TR_SUCCESS) then
                call mkl_free_buffers
                return
            else
                if(info(1)/=0.or.info(2)/=0.or.info(3)/=0.or.info(4)/=0.or.info(5)/=0.or.info(6)/=0) then
                    call mkl_free_buffers
                    return
                end if
            end if
            if(present(Jacobian)) then
                do!Main loop
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
                            i=Jacobian(J,x,M,N)
                    end select
                end do
            else
                do!Main loop
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
                            if(djacobi(fd_j,N,M,J,x,1d-8)/=TR_SUCCESS) then
                                call mkl_free_buffers
                                return
                            end if
                    end select
                end do
            end if
            !Clean up
                if (dtrnlspbc_get(handle,TotalIteration,StopReason,InitialResidual,FinalResidual)/=TR_SUCCESS) then
                    call mkl_free_buffers
                    return
                end if
                if (dtrnlspbc_delete(handle)/=TR_SUCCESS) then
                    call mkl_free_buffers
                    return
                end if
            do i=1,N
                if((x(i)<low(i).or.x(i)>up(i)).and.warn) then
                    write(*,'(1x,A49)')'Failed trust region: boundary condition violated!'
                    exit
                end if
            end do
        else
            if(dtrnlsp_init(handle,N,M,x,tol,maxit,maxstepit,StepBound)/=TR_SUCCESS) then
                call mkl_free_buffers
                return
            end if
            if(dtrnlsp_check(handle,N,M,J,fdx,tol,info)/=TR_SUCCESS) then
                call mkl_free_buffers
                return
            else
                if(info(1)/=0.or.info(2)/=0.or.info(3)/=0.or.info(4)/=0) then
                    call mkl_free_buffers
                    return
                end if
            end if
            if(present(Jacobian)) then
                do!Main loop
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
                            i=Jacobian(J,x,M,N)
                    end select
                end do
            else
                do!Main loop
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
                            if(djacobi(fd_j,N,M,J,x,1d-8)/=TR_SUCCESS) then
                                call mkl_free_buffers
                                return
                            end if
                    end select
                end do
            end if
            !Clean up
                if (dtrnlsp_get(handle,TotalIteration,StopReason,InitialResidual,FinalResidual)/=TR_SUCCESS) then
                    call mkl_free_buffers
                    return
                end if
                if (dtrnlsp_delete(handle)/=TR_SUCCESS) then
                    call mkl_free_buffers
                    return
                end if
        end if
        call mkl_free_buffers
        if(StopReason/=3) then
            select case(StopReason)
                case(1)
                    if(warn) then
                        write(*,*)'Failed trust region: max iteration exceeded!'
                        write(*,*)'Final residual =',FinalResidual
                    end if
                case(4)
                    if(warn) then
                        write(*,*)'Failed trust region: singular Jacobian encountered!'
                        write(*,*)'Final residual =',FinalResidual
                    end if
                case default
                    if(warn) then
                        write(*,'(1x,A87)')'Trust region warning: step length has converged, but residual has not met accuracy goal'
                        write(*,'(1x,A56)')'A best estimation rather than exact solution is returned'
                        write(*,*)'Final residual =',FinalResidual
                    end if
            end select
        end if
        contains
            subroutine fd_j(M,N,x,fdx)!Reformat fd for djacobi
                integer,intent(in)::M,N
                real*8,dimension(N),intent(in)::x
                real*8,dimension(M),intent(out)::fdx
                call fd(fdx,x,N)
            end subroutine fd_j
    end subroutine TrustRegion
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