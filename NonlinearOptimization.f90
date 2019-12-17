!Nonlinear optimization routines
!
!Instruction:
!If you want to solve nonlinear equation(s), you may:
!    adopt MKL trust region solver
!    or build your own merit function then use other routines to search for its minimum
!This module mainly provides unconstrained local minimizer
!    only solvers in Equality constraint section are constrained
!    only solvers in Heuristic section are global
!Global optimization remains an unsolved problem, since it is NP-complete
!    Heuristic algorithm is general but naive: search entire space (pay exponentially)
!
!Reference: J. Nocedal, S. J. Wright, *Numerical Optimization 2nd edition* (Springer, 2006)
module NonlinearOptimization
    use LinearAlgebra
    use mkl_rci_type; use mkl_rci!For numerical exact Hessian & Trust region solver, disable them if you have no access to MKL
    implicit none

contains
!-------------- Line search --------------
    !Suggestion:
    !    If dimensionality is low, adopt quasi-Newton (or Newton if Hessian is cheap and initial guess is close)
    !    If dimensionality is so high that O(dim^2) memory is unaffordable, adopt conjugate gradient or L-BFGS
    !Nomenclature:
    !    f = the target function to be minimized
    !    a = the line search step length
    !    p = the line search direction
    !    phi(a) = f( x + a * p ), so phi'(a) = f'( x + a * p ) . p
    !External procedure:
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
    !    Increment: see line searcher
    !On input x is an initial guess, on exit x is a local minimum of f(x)

    !Newton-Raphson method, requiring Wolfe condition
    !Optional: fdd: presence means analytical Hessian is available, otherwise call djacobi for central difference Hessian
    subroutine NewtonRaphson(f,fd,x,dim,fdd,f_fd,Strong,Warning,MaxIteration,Precision,MinStepLength,WolfeConst1,WolfeConst2,Increment)
        !Required argument
            external::f,fd
            integer,intent(in)::dim
            real*8,dimension(dim),intent(inout)::x
        !Optional argument
            integer,external,optional::fdd,f_fd
            logical,intent(in),optional::Strong,Warning
            integer,intent(in),optional::MaxIteration
            real*8,intent(in),optional::Precision,MinStepLength,WolfeConst1,WolfeConst2,Increment
        logical::sw,warn,terminate
        integer::maxit,iIteration,info
        real*8::tol,minstep,c1,c2,a,fnew,phidnew,phidold
        real*8,dimension(dim)::p,fdnew
        real*8,dimension(dim,dim)::Hessian
        !Initialize
            terminate=.false.
            !Set parameter according to optional argument
                if(present(Strong)) then; sw=Strong
                    else; sw=.true.; end if
                if(present(Warning)) then; warn=Warning
                    else; warn=.true.; end if
                if(present(MaxIteration)) then; maxit=MaxIteration
                    else; maxit=1000; end if
                if(present(Precision)) then; tol=Precision*Precision!To save sqrt cost, precision is squared
                    else; tol=1d-30; end if
                if(present(MinStepLength)) then; minstep=MinStepLength*MinStepLength!To save sqrt cost, MinStepLength is squared
                    else; minstep=1d-30; end if
                if(present(WolfeConst1)) then; c1=max(1d-15,WolfeConst1)!Fail safe
                    else; c1=1d-4; end if
                if(present(WolfeConst2)) then; c2=min(1d0-1d-15,max(c1+1d-15,WolfeConst2))!Fail safe
                    else; c2=0.9d0; end if
            if(present(f_fd)) then!Initial f(x) & f'(x)
                info=f_fd(fnew,fdnew,x,dim)
            else
                call f(fnew,x,dim); call fd(fdnew,x,dim)
            end if
            !Initial direction & step length
            if(present(fdd)) then; info=fdd(Hessian,x,dim)
            else; info=djacobi(fd_j,dim,dim,Hessian,x,1d-8); end if
            p=-fdnew; call My_dposv(Hessian,p,dim,info)
            if(info==0) then
                phidnew=dot_product(fdnew,p); a=1d0
            else!Hessian is not positive definite, use steepest descent direction
                p=-fdnew; phidnew=-dot_product(fdnew,fdnew)
                if(-phidnew<tol) return
                if(fnew==0d0) then; a=1d0
                else; a=-fnew/phidnew; end if
            end if
        if(present(Increment)) then
            if(present(fdd)) then!Analytical Hessian available
                if(present(f_fd)) then!Cheaper to evaluate f' along with f
                    if(sw) then!Use strong Wolfe condition instead of Wolfe condition
                        do iIteration=1,maxit!Main loop
                            phidold=phidnew!Prepare
                            call StrongWolfe_fdwithf(c1,c2,f,fd,f_fd,x,a,p,fnew,phidnew,fdnew,dim,Increment=Increment)!Line search
                            call After()!After search
                            if(terminate) return
                        end do
                    else
                        do iIteration=1,maxit!Main loop
                            phidold=phidnew!Prepare
                            call Wolfe_fdwithf(c1,c2,f,fd,f_fd,x,a,p,fnew,phidnew,fdnew,dim,Increment=Increment)!Line search
                            call After()!After search
                            if(terminate) return
                        end do
                    end if
                else
                    if(sw) then!Use strong Wolfe condition instead of Wolfe condition
                        do iIteration=1,maxit!Main loop
                            phidold=phidnew!Prepare
                            call StrongWolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim,Increment=Increment)!Line search
                            call After()!After search
                            if(terminate) return
                        end do
                    else
                        do iIteration=1,maxit!Main loop
                            phidold=phidnew!Prepare
                            call Wolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim,Increment=Increment)!Line search
                            call After()!After search
                            if(terminate) return
                        end do
                    end if
                end if
            else
                if(present(f_fd)) then!Cheaper to evaluate f' along with f
                    if(sw) then!Use strong Wolfe condition instead of Wolfe condition
                        do iIteration=1,maxit!Main loop
                            phidold=phidnew!Prepare
                            call StrongWolfe_fdwithf(c1,c2,f,fd,f_fd,x,a,p,fnew,phidnew,fdnew,dim,Increment=Increment)!Line search
                            call After_NumericalHessian()!After search
                            if(terminate) return
                        end do
                    else
                        do iIteration=1,maxit!Main loop
                            phidold=phidnew!Prepare
                            call Wolfe_fdwithf(c1,c2,f,fd,f_fd,x,a,p,fnew,phidnew,fdnew,dim,Increment=Increment)!Line search
                            call After_NumericalHessian()!After search
                            if(terminate) return
                        end do
                    end if
                else
                    if(sw) then!Use strong Wolfe condition instead of Wolfe condition
                        do iIteration=1,maxit!Main loop
                            phidold=phidnew!Prepare
                            call StrongWolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim,Increment=Increment)!Line search
                            call After_NumericalHessian()!After search
                            if(terminate) return
                        end do
                    else
                        do iIteration=1,maxit!Main loop
                            phidold=phidnew!Prepare
                            call Wolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim,Increment=Increment)!Line search
                            call After_NumericalHessian()!After search
                            if(terminate) return
                        end do
                    end if
                end if
            end if
        else
            if(present(fdd)) then!Analytical Hessian available
                if(present(f_fd)) then!Cheaper to evaluate f' along with f
                    if(sw) then!Use strong Wolfe condition instead of Wolfe condition
                        do iIteration=1,maxit!Main loop
                            phidold=phidnew!Prepare
                            call StrongWolfe_fdwithf(c1,c2,f,fd,f_fd,x,a,p,fnew,phidnew,fdnew,dim)!Line search
                            call After()!After search
                            if(terminate) return
                        end do
                    else
                        do iIteration=1,maxit!Main loop
                            phidold=phidnew!Prepare
                            call Wolfe_fdwithf(c1,c2,f,fd,f_fd,x,a,p,fnew,phidnew,fdnew,dim)!Line search
                            call After()!After search
                            if(terminate) return
                        end do
                    end if
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
                    if(sw) then!Use strong Wolfe condition instead of Wolfe condition
                        do iIteration=1,maxit!Main loop
                            phidold=phidnew!Prepare
                            call StrongWolfe_fdwithf(c1,c2,f,fd,f_fd,x,a,p,fnew,phidnew,fdnew,dim)!Line search
                            call After_NumericalHessian()!After search
                            if(terminate) return
                        end do
                    else
                        do iIteration=1,maxit!Main loop
                            phidold=phidnew!Prepare
                            call Wolfe_fdwithf(c1,c2,f,fd,f_fd,x,a,p,fnew,phidnew,fdnew,dim)!Line search
                            call After_NumericalHessian()!After search
                            if(terminate) return
                        end do
                    end if
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
        end if
        if(iIteration>maxit.and.warn) then
            write(*,'(1x,A46)')'Failed Newton-Raphson: max iteration exceeded!'
            write(*,*)'Euclidean norm of gradient =',Norm2(fdnew)
        end if
        contains
        subroutine After()!Check convergence, determine new direction & step length
            !Check convergence
            phidnew=dot_product(fdnew,fdnew)
            if(phidnew<tol) then
                terminate=.true.; return
            end if
            if(dot_product(p,p)*a*a<minstep) then
                if(warn) then
                    write(*,'(1x,A94)')'Newton-Raphson warning: step length has converged, but gradient norm has not met accuracy goal'
                    write(*,'(1x,A56)')'A best estimation rather than exact solution is returned'
                    write(*,*)'Euclidean norm of gradient =',dSqrt(phidnew)
                end if
                terminate=.true.; return
            end if
            !Determine new direction and step length
            p=-fdnew; info=fdd(Hessian,x,dim); call My_dposv(Hessian,p,dim,info)
            if(info==0) then
                phidnew=dot_product(fdnew,p); a=1d0
            else!Hessian is not positive definite, use steepest descent direction
                p=-fdnew; phidnew=-phidnew; a=a*phidold/phidnew
            end if
        end subroutine After
        subroutine After_NumericalHessian()!Check convergence, determine new direction & step length
            !Check convergence
            phidnew=dot_product(fdnew,fdnew)
            if(phidnew<tol) then
                terminate=.true.; return
            end if
            if(dot_product(p,p)*a*a<minstep) then
                if(warn) then
                    write(*,'(1x,A94)')'Newton-Raphson warning: step length has converged, but gradient norm has not met accuracy goal'
                    write(*,'(1x,A56)')'A best estimation rather than exact solution is returned'
                    write(*,*)'Euclidean norm of gradient =',dSqrt(phidnew)
                end if
                terminate=.true.; return
            end if
            !Determine new direction and step length
            p=-fdnew; info=djacobi(fd_j,dim,dim,Hessian,x,1d-8); call My_dposv(Hessian,p,dim,info)
            if(info==0) then
                phidnew=dot_product(fdnew,p); a=1d0
            else!Hessian is not positive definite, use steepest descent direction
                p=-fdnew; phidnew=-phidnew; a=a*phidold/phidnew
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
    subroutine BFGS(f,fd,x,dim,fdd,ExactStep,f_fd,Strong,Warning,MaxIteration,Precision,MinStepLength,WolfeConst1,WolfeConst2,Increment)
        !Required argument
            external::f,fd
            integer,intent(in)::dim
            real*8,dimension(dim),intent(inout)::x
        !Optional argument
            integer,external,optional::fdd,f_fd
            logical,intent(in),optional::Strong,Warning
            integer,intent(in),optional::ExactStep,MaxIteration
            real*8,intent(in),optional::Precision,MinStepLength,WolfeConst1,WolfeConst2,Increment
        logical::sw,warn,terminate
        integer::freq,maxit,iIteration,i
        real*8::tol,minstep,c1,c2,a,fnew,phidnew,rho
        real*8,dimension(dim)::p,fdnew,s,y
        real*8,dimension(dim,dim)::U,H!Approximate inverse Hessian
        !Initialize
            terminate=.false.
            !Set parameter according to optional argument
                if(present(ExactStep)) then; freq=ExactStep
                    else; freq=20; end if
                if(present(Strong)) then; sw=Strong
                    else; sw=.true.; end if
                if(present(Warning)) then; warn=Warning
                    else; warn=.true.; end if
                if(present(MaxIteration)) then; maxit=MaxIteration
                    else; maxit=1000; end if
                if(present(Precision)) then; tol=Precision*Precision!To save sqrt cost, precision is squared
                    else; tol=1d-30; end if
                if(present(MinStepLength)) then; minstep=MinStepLength*MinStepLength!To save sqrt cost, MinStepLength is squared
                    else; minstep=1d-30; end if
                if(present(WolfeConst1)) then; c1=max(1d-15,WolfeConst1)!Fail safe
                    else; c1=1d-4; end if
                if(present(WolfeConst2)) then; c2=min(1d0-1d-15,max(c1+1d-15,WolfeConst2))!Fail safe
                    else; c2=0.9d0; end if
            if(present(f_fd)) then!Initial f(x) & f'(x)
                i=f_fd(fnew,fdnew,x,dim)
            else
                call f(fnew,x,dim); call fd(fdnew,x,dim)
            end if
            !Initial approximate inverse Hessian & direction & step length
            if(freq>0) then
                if(present(fdd)) then; i=fdd(H,x,dim)
                else; i=djacobi(fd_j,dim,dim,H,x,1d-8); end if
                p=-fdnew; call My_dpotri(H,dim,i)
                if(i==0) then
                    call syL2U(H,dim)
                    p=-matmul(H,fdnew); phidnew=dot_product(fdnew,p); a=1d0
                end if
            end if
            if(freq<=0.or.i/=0) then!Hessian is either uncomputed or not positive definite, initial approximate inverse Hessian = a
                p=-fdnew; phidnew=-dot_product(fdnew,fdnew)
                if(-phidnew<tol) return
                if(fnew==0d0) then; a=1d0
                else; a=-fnew/phidnew; end if
                s=x; y=fdnew
                if(present(Increment)) then
                    if(sw) then
                        call StrongWolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim,Increment=Increment)
                    else
                        call Wolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim,Increment=Increment)
                    end if
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
                s=x-s; y=fdnew-y; rho=1d0/dot_product(y,s)
                U=-rho*vector_direct_product(y,s,dim,dim)
                forall(i=1:dim); U(i,i)=U(i,i)+1d0; end forall
                H=matmul(transpose(U),a*U)+rho*vector_direct_product(s,s,dim,dim)
                p=-matmul(H,fdnew); phidnew=dot_product(fdnew,p); a=1d0
            end if
        if(present(Increment)) then
            if(freq>0) then!Exact Hessian will be computed
                if(present(fdd)) then!Analytical Hessian is available
                    if(present(f_fd)) then!Cheaper to evaluate f' along with f
                        if(sw) then!Use strong Wolfe condition instead of Wolfe condition
                            do iIteration=1,maxit!Main loop
                                s=x; y=fdnew!Prepare
                                call StrongWolfe_fdwithf(c1,c2,f,fd,f_fd,x,a,p,fnew,phidnew,fdnew,dim,Increment=Increment)!Line search
                                call After()!After search
                                if(terminate) return
                            end do
                        else
                            do iIteration=1,maxit!Main loop
                                s=x; y=fdnew!Prepare
                                call Wolfe_fdwithf(c1,c2,f,fd,f_fd,x,a,p,fnew,phidnew,fdnew,dim,Increment=Increment)!Line search
                                call After()!After search
                                if(terminate) return
                            end do
                        end if
                    else
                        if(sw) then!Use strong Wolfe condition instead of Wolfe condition
                            do iIteration=1,maxit!Main loop
                                s=x; y=fdnew!Prepare
                                call StrongWolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim,Increment=Increment)!Line search
                                call After()!After search
                                if(terminate) return
                            end do
                        else
                            do iIteration=1,maxit!Main loop
                                s=x; y=fdnew!Prepare
                                call Wolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim,Increment=Increment)!Line search
                                call After()!After search
                                if(terminate) return
                            end do
                        end if
                    end if
                else
                    if(present(f_fd)) then!Cheaper to evaluate f' along with f
                        if(sw) then!Use strong Wolfe condition instead of Wolfe condition
                            do iIteration=1,maxit!Main loop
                                s=x; y=fdnew!Prepare
                                call StrongWolfe_fdwithf(c1,c2,f,fd,f_fd,x,a,p,fnew,phidnew,fdnew,dim,Increment=Increment)!Line search
                                call After_NumericalHessian()!After search
                                if(terminate) return
                            end do
                        else
                            do iIteration=1,maxit!Main loop
                                s=x; y=fdnew!Prepare
                                call Wolfe_fdwithf(c1,c2,f,fd,f_fd,x,a,p,fnew,phidnew,fdnew,dim,Increment=Increment)!Line search
                                call After_NumericalHessian()!After search
                                if(terminate) return
                            end do
                        end if
                    else
                        if(sw) then!Use strong Wolfe condition instead of Wolfe condition
                            do iIteration=1,maxit!Main loop
                                s=x; y=fdnew!Prepare
                                call StrongWolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim,Increment=Increment)!Line search
                                call After_NumericalHessian()!After search
                                if(terminate) return
                            end do
                        else
                            do iIteration=1,maxit!Main loop
                                s=x; y=fdnew!Prepare
                                call Wolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim,Increment=Increment)!Line search
                                call After_NumericalHessian()!After search
                                if(terminate) return
                            end do
                        end if
                    end if
                end if
            else
                if(present(f_fd)) then!Cheaper to evaluate f' along with f
                    if(sw) then!Use strong Wolfe condition instead of Wolfe condition
                        do iIteration=1,maxit!Main loop
                            s=x; y=fdnew!Prepare
                            call StrongWolfe_fdwithf(c1,c2,f,fd,f_fd,x,a,p,fnew,phidnew,fdnew,dim,Increment=Increment)!Line search
                            call After_NoHessian()!After search
                            if(terminate) return
                        end do
                    else
                        do iIteration=1,maxit!Main loop
                            s=x; y=fdnew!Prepare
                            call Wolfe_fdwithf(c1,c2,f,fd,f_fd,x,a,p,fnew,phidnew,fdnew,dim,Increment=Increment)!Line search
                            call After_NoHessian()!After search
                            if(terminate) return
                        end do
                    end if
                else
                    if(sw) then!Use strong Wolfe condition instead of Wolfe condition
                        do iIteration=1,maxit!Main loop
                            s=x; y=fdnew!Prepare
                            call StrongWolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim,Increment=Increment)!Line search
                            call After_NoHessian()!After search
                            if(terminate) return
                        end do
                    else
                        do iIteration=1,maxit!Main loop
                            s=x; y=fdnew!Prepare
                            call Wolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim,Increment=Increment)!Line search
                            call After_NoHessian()!After search
                            if(terminate) return
                        end do
                    end if
                end if
            end if
        else
            if(freq>0) then!Exact Hessian will be computed
                if(present(fdd)) then!Analytical Hessian is available
                    if(present(f_fd)) then!Cheaper to evaluate f' along with f
                        if(sw) then!Use strong Wolfe condition instead of Wolfe condition
                            do iIteration=1,maxit!Main loop
                                s=x; y=fdnew!Prepare
                                call StrongWolfe_fdwithf(c1,c2,f,fd,f_fd,x,a,p,fnew,phidnew,fdnew,dim)!Line search
                                call After()!After search
                                if(terminate) return
                            end do
                        else
                            do iIteration=1,maxit!Main loop
                                s=x; y=fdnew!Prepare
                                call Wolfe_fdwithf(c1,c2,f,fd,f_fd,x,a,p,fnew,phidnew,fdnew,dim)!Line search
                                call After()!After search
                                if(terminate) return
                            end do
                        end if
                    else
                        if(sw) then!Use strong Wolfe condition instead of Wolfe condition
                            do iIteration=1,maxit!Main loop
                                s=x; y=fdnew!Prepare
                                call StrongWolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim)!Line search
                                call After()!After search
                                if(terminate) return
                            end do
                        else
                            do iIteration=1,maxit!Main loop
                                s=x; y=fdnew!Prepare
                                call Wolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim)!Line search
                                call After()!After search
                                if(terminate) return
                            end do
                        end if
                    end if
                else
                    if(present(f_fd)) then!Cheaper to evaluate f' along with f
                        if(sw) then!Use strong Wolfe condition instead of Wolfe condition
                            do iIteration=1,maxit!Main loop
                                s=x; y=fdnew!Prepare
                                call StrongWolfe_fdwithf(c1,c2,f,fd,f_fd,x,a,p,fnew,phidnew,fdnew,dim)!Line search
                                call After_NumericalHessian()!After search
                                if(terminate) return
                            end do
                        else
                            do iIteration=1,maxit!Main loop
                                s=x; y=fdnew!Prepare
                                call Wolfe_fdwithf(c1,c2,f,fd,f_fd,x,a,p,fnew,phidnew,fdnew,dim)!Line search
                                call After_NumericalHessian()!After search
                                if(terminate) return
                            end do
                        end if
                    else
                        if(sw) then!Use strong Wolfe condition instead of Wolfe condition
                            do iIteration=1,maxit!Main loop
                                s=x; y=fdnew!Prepare
                                call StrongWolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim)!Line search
                                call After_NumericalHessian()!After search
                                if(terminate) return
                            end do
                        else
                            do iIteration=1,maxit!Main loop
                                s=x; y=fdnew!Prepare
                                call Wolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim)!Line search
                                call After_NumericalHessian()!After search
                                if(terminate) return
                            end do
                        end if
                    end if
                end if
            else
                if(present(f_fd)) then!Cheaper to evaluate f' along with f
                    if(sw) then!Use strong Wolfe condition instead of Wolfe condition
                        do iIteration=1,maxit!Main loop
                            s=x; y=fdnew!Prepare
                            call StrongWolfe_fdwithf(c1,c2,f,fd,f_fd,x,a,p,fnew,phidnew,fdnew,dim)!Line search
                            call After_NoHessian()!After search
                            if(terminate) return
                        end do
                    else
                        do iIteration=1,maxit!Main loop
                            s=x; y=fdnew!Prepare
                            call Wolfe_fdwithf(c1,c2,f,fd,f_fd,x,a,p,fnew,phidnew,fdnew,dim)!Line search
                            call After_NoHessian()!After search
                            if(terminate) return
                        end do
                    end if
                else
                    if(sw) then!Use strong Wolfe condition instead of Wolfe condition
                        do iIteration=1,maxit!Main loop
                            s=x; y=fdnew!Prepare
                            call StrongWolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim)!Line search
                            call After_NoHessian()!After search
                            if(terminate) return
                        end do
                    else
                        do iIteration=1,maxit!Main loop
                            s=x; y=fdnew!Prepare
                            call Wolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim)!Line search
                            call After_NoHessian()!After search
                            if(terminate) return
                        end do
                    end if
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
                    terminate=.true.; return
                end if
                if(dot_product(p,p)*a*a<minstep) then
                    if(warn) then
                        write(*,'(1x,A84)')'BFGS warning: step length has converged, but gradient norm has not met accuracy goal'
                        write(*,'(1x,A56)')'A best estimation rather than exact solution is returned'
                        write(*,*)'Euclidean norm of gradient =',dSqrt(phidnew)
                    end if
                    terminate=.true.; return
                end if
                !Determine new direction and step length, update approximate inverse Hessian
                i=mod(iIteration,freq)
                if(i==0) then!Every freq steps compute exact Hessian
                    i=fdd(U,x,dim); call My_dpotri(U,dim,i)
                    if(i==0) then!Use exact Hessian if positive definite
                        call sycp(H,U,dim); call syL2U(H,dim)
                        p=-matmul(H,fdnew); phidnew=dot_product(fdnew,p); a=1d0
                    end if
                end if
                if(i/=0) then!Exact Hessian is either uncomputed or not positive definite, update approximate Hessian
                    s=x-s; y=fdnew-y; rho=1d0/dot_product(y,s)
                    U=-rho*vector_direct_product(y,s,dim,dim)
                    forall(i=1:dim); U(i,i)=U(i,i)+1d0; end forall
                    H=matmul(transpose(U),matmul(H,U))+rho*vector_direct_product(s,s,dim,dim)
                    p=-matmul(H,fdnew); phidnew=dot_product(fdnew,p); a=1d0
                end if
            end subroutine After
            subroutine After_NumericalHessian()!Check convergence, update approximate inverse Hessian, determine new direction & step length
                !Check convergence
                phidnew=dot_product(fdnew,fdnew)
                if(phidnew<tol) then
                    terminate=.true.; return
                end if
                if(dot_product(p,p)*a*a<minstep) then
                    if(warn) then
                        write(*,'(1x,A84)')'BFGS warning: step length has converged, but gradient norm has not met accuracy goal'
                        write(*,'(1x,A56)')'A best estimation rather than exact solution is returned'
                        write(*,*)'Euclidean norm of gradient =',dSqrt(phidnew)
                    end if
                    terminate=.true.; return
                end if
                !Determine new direction and step length, update approximate inverse Hessian
                i=mod(iIteration,freq)
                if(i==0) then!Every freq steps compute exact Hessian
                    i=djacobi(fd_j,dim,dim,U,x,1d-8)
                    call My_dpotri(U,dim,i)
                    if(i==0) then!Use exact Hessian if positive definite
                        call sycp(H,U,dim); call syL2U(H,dim)
                        p=-matmul(H,fdnew); phidnew=dot_product(fdnew,p); a=1d0
                    end if
                end if
                if(i/=0) then!Exact Hessian is either uncomputed or not positive definite, update approximate Hessian
                    s=x-s; y=fdnew-y; rho=1d0/dot_product(y,s)
                    U=-rho*vector_direct_product(y,s,dim,dim)
                    forall(i=1:dim); U(i,i)=U(i,i)+1d0; end forall
                    H=matmul(transpose(U),matmul(H,U))+rho*vector_direct_product(s,s,dim,dim)
                    p=-matmul(H,fdnew); phidnew=dot_product(fdnew,p); a=1d0
                end if
            end subroutine After_NumericalHessian
            subroutine After_NoHessian()!Check convergence, update approximate inverse Hessian, determine new direction & step length
                !Check convergence
                phidnew=dot_product(fdnew,fdnew)
                if(phidnew<tol) then
                    terminate=.true.; return
                end if
                if(dot_product(p,p)*a*a<minstep) then
                    if(warn) then
                        write(*,'(1x,A84)')'BFGS warning: step length has converged, but gradient norm has not met accuracy goal'
                        write(*,'(1x,A56)')'A best estimation rather than exact solution is returned'
                        write(*,*)'Euclidean norm of gradient =',dSqrt(phidnew)
                    end if
                    terminate=.true.; return
                end if
                !Determine new direction and step length, update approximate inverse Hessian
                s=x-s; y=fdnew-y; rho=1d0/dot_product(y,s)
                U=-rho*vector_direct_product(y,s,dim,dim)
                forall(i=1:dim); U(i,i)=U(i,i)+1d0; end forall
                H=matmul(transpose(U),matmul(H,U))+rho*vector_direct_product(s,s,dim,dim)
                p=-matmul(H,fdnew); phidnew=dot_product(fdnew,p); a=1d0
            end subroutine After_NoHessian
            subroutine fd_j(M,N,x,fdx)!Reformat for djacobi
                integer,intent(in)::M,N
                real*8,dimension(N),intent(in)::x
                real*8,dimension(M),intent(out)::fdx
                call fd(fdx,x,N)
            end subroutine fd_j
    end subroutine BFGS

    !Limited-memory Broyden–Fletcher–Goldfarb–Shanno (L-BFGS) quasi-Newton method, requiring Wolfe condition
    !Optional: Memory: (default = 10) memory usage = O( Memory * dim ). [3,30] is recommended (must > 0)
    subroutine LBFGS(f,fd,x,dim,Memory,f_fd,Strong,Warning,MaxIteration,Precision,MinStepLength,WolfeConst1,WolfeConst2,Increment)
        !Required argument
            external::f,fd
            integer,intent(in)::dim
            real*8,dimension(dim),intent(inout)::x
        !Optional argument
            integer,external,optional::f_fd
            logical,intent(in),optional::Strong,Warning
            integer,intent(in),optional::Memory,MaxIteration
            real*8,intent(in),optional::Precision,MinStepLength,WolfeConst1,WolfeConst2,Increment
        logical::sw,warn,terminate
        integer::mem,maxit,iIteration,i,recent
        real*8::tol,minstep,c1,c2,a,fnew,phidnew
        real*8,dimension(dim)::p,fdnew,xold,fdold
        real*8,allocatable,dimension(:)::rho,alpha
        real*8,allocatable,dimension(:,:)::s,y
        !Initialize
            terminate=.false.
            !Set parameter according to optional argument
                if(present(Memory)) then; mem=max(1,Memory)!Fail safe
                    else; mem=10; end if
                if(present(Strong)) then; sw=Strong
                    else; sw=.true.; end if
                if(present(Warning)) then; warn=Warning
                    else; warn=.true.; end if
                if(present(MaxIteration)) then; maxit=MaxIteration
                    else; maxit=1000; end if
                if(present(Precision)) then; tol=Precision*Precision!To save sqrt cost, precision is squared
                    else; tol=1d-30; end if
                if(present(MinStepLength)) then; minstep=MinStepLength*MinStepLength!To save sqrt cost, MinStepLength is squared
                    else; minstep=1d-30; end if
                if(present(WolfeConst1)) then; c1=max(1d-15,WolfeConst1)!Fail safe
                    else; c1=1d-4; end if
                if(present(WolfeConst2)) then; c2=min(1d0-1d-15,max(c1+1d-15,WolfeConst2))!Fail safe
                    else; c2=0.9d0; end if
            allocate(rho(0:mem)); allocate(alpha(0:mem)); allocate(s(dim,0:mem)); allocate(y(dim,0:mem))
            if(present(f_fd)) then!Initial f(x) & f'(x)
                i=f_fd(fnew,fdnew,x,dim)
            else
                call f(fnew,x,dim); call fd(fdnew,x,dim)
            end if
            !Initial iteration history
            p=-fdnew; phidnew=-dot_product(fdnew,fdnew)
            if(-phidnew<tol) return
            if(fnew==0d0) then; a=1d0
            else; a=-fnew/phidnew; end if
            xold=x; fdold=fdnew
            !Initial approximate inverse Hessian = a
            if(present(Increment)) then
                if(sw) then
                    call StrongWolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim,Increment=Increment)
                else
                    call Wolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim,Increment=Increment)
                end if
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
            s(:,0)=x-xold; y(:,0)=fdnew-fdold; rho(0)=1d0/dot_product(y(:,0),s(:,0))
            do iIteration=1,mem-1!Preiterate to get enough history
                xold=x; fdold=fdnew!Prepare
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
                p=-p; phidnew=dot_product(fdnew,p); a=1d0
                if(present(Increment)) then!Line search
                    if(sw) then
                        call StrongWolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim,Increment=Increment)
                    else
                        call Wolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim,Increment=Increment)
                    end if
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
                s(:,recent)=x-xold; y(:,recent)=fdnew-fdold; rho(recent)=1d0/dot_product(y(:,recent),s(:,recent))
            end do
        if(present(Increment)) then
            if(present(f_fd)) then!Cheaper to evaluate f' along with f
                if(sw) then!Use strong Wolfe condition instead of Wolfe condition
                    do iIteration=1,maxit!Main loop
                        call Before()!Before search
                        call StrongWolfe_fdwithf(c1,c2,f,fd,f_fd,x,a,p,fnew,phidnew,fdnew,dim,Increment=Increment)!Line search
                        call After()!After search
                        if(terminate) return
                    end do
                else
                    do iIteration=1,maxit!Main loop
                        call Before()!Before search
                        call Wolfe_fdwithf(c1,c2,f,fd,f_fd,x,a,p,fnew,phidnew,fdnew,dim,Increment=Increment)!Line search
                        call After()!After search
                        if(terminate) return
                    end do
                end if
            else
                if(sw) then!Use strong Wolfe condition instead of Wolfe condition
                    do iIteration=1,maxit!Main loop
                        call Before()!Before search
                        call StrongWolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim,Increment=Increment)!Line search
                        call After()!After search
                        if(terminate) return
                    end do
                else
                    do iIteration=1,maxit!Main loop
                        call Before()!Before search
                        call Wolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim,Increment=Increment)!Line search
                        call After()!After search
                        if(terminate) return
                    end do
                end if
            end if
        else
            if(present(f_fd)) then!Cheaper to evaluate f' along with f
                if(sw) then!Use strong Wolfe condition instead of Wolfe condition
                    do iIteration=1,maxit!Main loop
                        call Before()!Before search
                        call StrongWolfe_fdwithf(c1,c2,f,fd,f_fd,x,a,p,fnew,phidnew,fdnew,dim)!Line search
                        call After()!After search
                        if(terminate) return
                    end do
                else
                    do iIteration=1,maxit!Main loop
                        call Before()!Before search
                        call Wolfe_fdwithf(c1,c2,f,fd,f_fd,x,a,p,fnew,phidnew,fdnew,dim)!Line search
                        call After()!After search
                        if(terminate) return
                    end do
                end if
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
        end if
        if(iIteration>maxit.and.warn) then
            write(*,'(1x,A38)')'Failed L-BFGS: max iteration exceeded!'
            write(*,*)'Euclidean norm of gradient =',Norm2(fdnew)
        end if
        deallocate(rho); deallocate(alpha); deallocate(s); deallocate(y)!Clean up
        contains
            subroutine Before()!Prepare, determine new direction & step length
                xold=x; fdold=fdnew!Prepare
                !Determine new direction
                p=fdnew
                do i=recent,0,-1
                    alpha(i)=rho(i)*dot_product(s(:,i),p)
                    p=p-alpha(i)*y(:,i)
                end do
                do i=mem-1,recent+1,-1
                    alpha(i)=rho(i)*dot_product(s(:,i),p)
                    p=p-alpha(i)*y(:,i)
                end do
                p=p/rho(recent)/dot_product(y(:,recent),y(:,recent))
                do i=recent+1,mem-1
                    phidnew=rho(i)*dot_product(y(:,i),p)
                    p=p+(alpha(i)-phidnew)*s(:,i)
                end do
                do i=0,recent
                    phidnew=rho(i)*dot_product(y(:,i),p)
                    p=p+(alpha(i)-phidnew)*s(:,i)
                end do
                p=-p; phidnew=dot_product(fdnew,p); a=1d0
            end subroutine Before
            subroutine After()!Check convergence, replace earliest iteration history with latest
                !Check convergence
                phidnew=dot_product(fdnew,fdnew)
                if(phidnew<tol) then
                    terminate=.true.; return
                end if
                if(dot_product(p,p)*a*a<minstep) then
                    if(warn) then
                        write(*,'(1x,A86)')'L-BFGS warning: step length has converged, but gradient norm has not met accuracy goal'
                        write(*,'(1x,A56)')'A best estimation rather than exact solution is returned'
                        write(*,*)'Euclidean norm of gradient =',dSqrt(phidnew)
                    end if
                    terminate=.true.; return
                end if
                recent=mod(recent+1,mem)
                s(:,recent)=x-xold; y(:,recent)=fdnew-fdold; rho(recent)=1d0/dot_product(y(:,recent),s(:,recent))
            end subroutine After
    end subroutine LBFGS

    !Conjugate gradient method, requiring either Wolfe or Strong Wolfe condition
    !Available methods: DY (Dai-Yun), PR (Polak-Ribiere+)
    !Optional: Method: (default = DY) which conjugate gradient method to use
    subroutine ConjugateGradient(f,fd,x,dim,Method,f_fd,Strong,Warning,MaxIteration,Precision,MinStepLength,WolfeConst1,WolfeConst2,Increment)
        !Required argument
            external::f,fd
            integer,intent(in)::dim
            real*8,dimension(dim),intent(inout)::x
        !Optional argument
            character*32,intent(in),optional::Method
            integer,external,optional::f_fd
            logical,intent(in),optional::Strong,Warning
            integer,intent(in),optional::MaxIteration
            real*8,intent(in),optional::Precision,MinStepLength,WolfeConst1,WolfeConst2,Increment
        logical::sw,warn,terminate
        character*32::type
        integer::maxit,iIteration,info
        real*8::tol,minstep,c1,c2,a,fnew,fold,phidnew,phidold
        real*8,dimension(dim)::p,fdnew,fdold
        !Initialize
            terminate=.false.
            !Set parameter according to optional argument
                if(present(Method)) then; type=Method
                    else; type='DY'; end if
                if(present(Strong)) then; sw=Strong
                    else; sw=.true.; end if
                if(present(Warning)) then; warn=Warning
                    else; warn=.true.; end if
                if(present(MaxIteration)) then; maxit=MaxIteration
                    else; maxit=1000; end if
                if(present(Precision)) then; tol=Precision*Precision!To save sqrt cost, precision is squared
                    else; tol=1d-30; end if
                if(present(MinStepLength)) then;  minstep=MinStepLength*MinStepLength!To save sqrt cost, MinStepLength is squared
                    else; minstep=1d-30; end if
                if(present(WolfeConst1)) then; c1=max(1d-15,WolfeConst1)!Fail safe
                    else; c1=1d-4; end if
                if(present(WolfeConst2)) then; c2=min(1d0-1d-15,max(c1+1d-15,WolfeConst2))!Fail safe
                    else; c2=0.45d0; end if
            if(present(f_fd)) then!Initial f(x) & f'(x)
                info=f_fd(fnew,fdnew,x,dim)
            else
                call f(fnew,x,dim); call fd(fdnew,x,dim)
            end if
            !Initial direction & step length
            p=-fdnew; phidnew=-dot_product(fdnew,fdnew)
            if(-phidnew<tol) return
            if(fnew==0d0) then; a=1d0
            else; a=-fnew/phidnew; end if
        select case(type)
            case('DY')!Require Wolfe condition
                if(present(Increment)) then
                    if(present(f_fd)) then!Cheaper to evaluate f' along with f
                        if(sw) then!Use strong Wolfe condition instead of Wolfe condition
                            do iIteration=1,maxit!Main loop
                                fold=fnew; fdold=fdnew; phidold=phidnew!Prepare
                                call StrongWolfe_fdwithf(c1,c2,f,fd,f_fd,x,a,p,fnew,phidnew,fdnew,dim,Increment=Increment)!Line search
                                call DY()!After search
                                if(terminate) return
                            end do
                        else
                            do iIteration=1,maxit!Main loop
                                fold=fnew; fdold=fdnew; phidold=phidnew!Prepare
                                call Wolfe_fdwithf(c1,c2,f,fd,f_fd,x,a,p,fnew,phidnew,fdnew,dim,Increment=Increment)!Line search
                                call DY()!After search
                                if(terminate) return
                            end do
                        end if
                    else
                        if(sw) then!To meet Nocedal performance suggestion, Dai-Yun requires strong Wolfe condition
                            do iIteration=1,maxit!Main loop
                                fold=fnew; fdold=fdnew; phidold=phidnew!Prepare
                                call StrongWolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim,Increment=Increment)!Line search
                                call DY()!After search
                                if(terminate) return
                            end do
                        else
                            do iIteration=1,maxit!Main loop
                                fold=fnew; fdold=fdnew; phidold=phidnew!Prepare
                                call Wolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim,Increment=Increment)!Line search
                                call DY()!After search
                                if(terminate) return
                            end do
                        end if
                    end if
                else
                    if(present(f_fd)) then!Cheaper to evaluate f' along with f
                        if(sw) then!Use strong Wolfe condition instead of Wolfe condition
                            do iIteration=1,maxit!Main loop
                                fold=fnew; fdold=fdnew; phidold=phidnew!Prepare
                                call StrongWolfe_fdwithf(c1,c2,f,fd,f_fd,x,a,p,fnew,phidnew,fdnew,dim)!Line search
                                call DY()!After search
                                if(terminate) return
                            end do
                        else
                            do iIteration=1,maxit!Main loop
                                fold=fnew; fdold=fdnew; phidold=phidnew!Prepare
                                call Wolfe_fdwithf(c1,c2,f,fd,f_fd,x,a,p,fnew,phidnew,fdnew,dim)!Line search
                                call DY()!After search
                                if(terminate) return
                            end do
                        end if
                    else
                        if(sw) then!To meet Nocedal performance suggestion, Dai-Yun requires strong Wolfe condition
                            do iIteration=1,maxit!Main loop
                                fold=fnew; fdold=fdnew; phidold=phidnew!Prepare
                                call StrongWolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim)!Line search
                                call DY()!After search
                                if(terminate) return
                            end do
                        else
                            do iIteration=1,maxit!Main loop
                                fold=fnew; fdold=fdnew; phidold=phidnew!Prepare
                                call Wolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim)!Line search
                                call DY()!After search
                                if(terminate) return
                            end do
                        end if
                    end if
                end if
            case('PR')!Require strong Wolfe condition
                if(present(Increment)) then
                    if(present(f_fd)) then!Cheaper to evaluate f' along with f
                        do iIteration=1,maxit!Main loop
                            fold=fnew; fdold=fdnew; phidold=phidnew!Prepare
                            call StrongWolfe_fdwithf(c1,c2,f,fd,f_fd,x,a,p,fnew,phidnew,fdnew,dim,Increment=Increment)!Line search
                            call PR()!After search
                            if(terminate) return
                        end do
                    else
                        do iIteration=1,maxit!Main loop
                            fold=fnew; fdold=fdnew; phidold=phidnew!Prepare
                            call StrongWolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim,Increment=Increment)!Line search
                            call PR()!After search
                            if(terminate) return
                        end do
                    end if
                else
                    if(present(f_fd)) then!Cheaper to evaluate f' along with f
                        do iIteration=1,maxit!Main loop
                            fold=fnew; fdold=fdnew; phidold=phidnew!Prepare
                            call StrongWolfe_fdwithf(c1,c2,f,fd,f_fd,x,a,p,fnew,phidnew,fdnew,dim)!Line search
                            call PR()!After search
                            if(terminate) return
                        end do
                    else
                        do iIteration=1,maxit!Main loop
                            fold=fnew; fdold=fdnew; phidold=phidnew!Prepare
                            call StrongWolfe(c1,c2,f,fd,x,a,p,fnew,phidnew,fdnew,dim)!Line search
                            call PR()!After search
                            if(terminate) return
                        end do
                    end if
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
                terminate=.true.; return
            end if
            if(dot_product(p,p)*a*a<minstep) then
                if(warn) then
                    write(*,'(1x,A106)')'Dai-Yun conjugate gradient warning: step length has converged, but gradient norm has not met accuracy goal'
                    write(*,'(1x,A56)')'A best estimation rather than exact solution is returned'
                    write(*,*)'Euclidean norm of gradient =',dSqrt(phidnew)
                end if
                terminate=.true.; return
            end if
            !Determine new direction
            p=-fdnew+dot_product(fdnew,fdnew)/dot_product(fdnew-fdold,p)*p
            phidnew=dot_product(fdnew,p)
            if(phidnew>0d0) then!Ascent direction, reset to steepest descent direction
                p=-fdnew; phidnew=-dot_product(fdnew,fdnew)
            end if
            a=a*phidold/phidnew
        end subroutine DY
        subroutine PR()!Check convergence, determine new direction & step length
            !Check convergence
            phidnew=dot_product(fdnew,fdnew)
            if(phidnew<tol) then
                terminate=.true.; return
            end if
            if(dot_product(p,p)*a*a<minstep) then
                if(warn) then
                    write(*,'(1x,A113)')'Polak-Ribiere+ conjugate gradient warning: step length has converged, but gradient norm has not met accuracy goal'
                    write(*,'(1x,A56)')'A best estimation rather than exact solution is returned'
                    write(*,*)'Euclidean norm of gradient =',dSqrt(phidnew)
                end if
                terminate=.true.; return
            end if
            !Determine new direction
            p=-fdnew+dot_product(fdnew,fdnew-fdold)/dot_product(fdold,fdold)*p
            phidnew=dot_product(fdnew,p)
            if(phidnew>0d0) then!Ascent direction, reset to steepest descent direction
                p=-fdnew; phidnew=-dot_product(fdnew,fdnew)
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
        !Optional: Increment: (default = 1.05) each iteration change a by how much time (must > 1)

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
            if(present(Increment)) then; incrmt=max(1d0+1d-15,Increment)!Fail safe
            else; incrmt=1.05d0; end if
            x0=x; fx0=fx; c2_m_abs_phid0=c2*dAbs(phid0)
            !Check whether initial guess satisfies sufficient decrease condition
            x=x0+a*p; call f(fx,x,dim)
            if(fx<=fx0+c1*a*phid0) then!Satisfied, search for larger a
                do
                    aold=a; fold=fx
                    a=aold*incrmt; x=x0+a*p; call f(fx,x,dim)
                    if(fx>fx0+c1*a*phid0) then
                        x=x0+aold*p
                        call fd(fdx,x,dim)
                        phidx=dot_product(fdx,p)
                        if(phidx>c2_m_abs_phid0) then
                            a=aold; fx=fold
                        else
                            atemp=a; ftemp=fx
                            call zoom(aold,atemp,fold,ftemp,phidx)
                        end if
                        return
                    end if
                end do
            else!Violated, search for smaller a
                do
                    aold=a; fold=fx
                    a=aold/incrmt; x=x0+a*p; call f(fx,x,dim)
                    if(fx<=fx0+c1*a*phid0) then
                        call fd(fdx,x,dim)
                        phidx=dot_product(fdx,p)
                        if(phidx<c2_m_abs_phid0) then
                            atemp=a; ftemp=fx
                            call zoom(atemp,aold,ftemp,fold,phidx)
                        end if
                        return
                    end if
                    if(a<1d-15) then
                        call fd(fdx,x,dim); return
                    end if
                end do
            end if
            contains
            !low & up must satisfy:
            !    low < up
            !    low satisfies sufficient decrease condition, but up violates
            !    phi'(low) < 0
            subroutine zoom(low,up,flow,fup,phidlow)!Restore x and fx to input status before calling zoom
                real*8,intent(inout)::low,up,flow,fup,phidlow
                real*8::phidnew,phidlow_m_a
                phidlow_m_a=phidlow*a!Initialize
                do
                    !Updata a by quadratic interpolation
                        a=phidlow_m_a*a/2d0/(flow+phidlow_m_a-fup)
                        if(.not.(a>low.and.a<up)) a=(low+up)/2d0!Fail safe
                    x=x0+a*p; call f(fx,x,dim)
                    if(fx>fx0+c1*a*phid0) then
                        up=a
                        if(up-low<1d-15.or.(up-low)/max(dAbs(low),dAbs(up))<1d-15) then
                            call fd(fdx,x,dim); return
                        end if
                        fup=fx
                    else
                        call fd(fdx,x,dim); phidnew=dot_product(fdx,p)
                        if(phidnew>c2_m_abs_phid0) return
                        low=a
                        if(up-low<1d-15.or.(up-low)/max(dAbs(low),dAbs(up))<1d-15) return
                        flow=fx; phidlow=phidnew; phidlow_m_a=phidlow*a
                    end if
                end do
            end subroutine zoom
        end subroutine Wolfe
        !When it is cheaper to evaluate f' along with f
        subroutine Wolfe_fdwithf(c1,c2,f,fd,f_fd,x,a,p,fx,phid0,fdx,dim,Increment)!CURRENTLY NO BETTER THAN Wolfe
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
            real*8::incrmt,c2_m_abs_phid0,fx0,ftemp,atemp,aold,fold,phidx
            real*8,dimension(dim)::x0
            !Initialize
            if(present(Increment)) then; incrmt=max(1d0+1d-15,Increment)!Fail safe
            else; incrmt=1.05d0; end if
            x0=x; fx0=fx; c2_m_abs_phid0=c2*dAbs(phid0)
            !Check whether initial guess satisfies sufficient decrease condition
            x=x0+a*p; call f(fx,x,dim)
            if(fx<=fx0+c1*a*phid0) then!Satisfied, search for larger a
                do
                    aold=a; fold=fx
                    a=aold*incrmt; x=x0+a*p
                    call f(fx,x,dim)
                    if(fx>fx0+c1*a*phid0) then
                        x=x0+aold*p; call fd(fdx,x,dim)
                        phidx=dot_product(fdx,p)
                        if(phidx>c2_m_abs_phid0) then
                            a=aold; fx=fold
                        else
                            atemp=a; ftemp=fx
                            call zoom(aold,atemp,fold,ftemp,phidx)
                        end if
                        return
                    end if
                end do
            else!Violated, search for smaller a
                do
                    aold=a; fold=fx
                    a=aold/incrmt; x=x0+a*p; call f(fx,x,dim)
                    if(fx<=fx0+c1*a*phid0) then
                        call fd(fdx,x,dim)
                        phidx=dot_product(fdx,p)
                        if(phidx<c2_m_abs_phid0) then
                            atemp=a; ftemp=fx
                            call zoom(atemp,aold,ftemp,fold,phidx)
                        end if
                        return
                    end if
                    if(a<1d-15) then
                        call fd(fdx,x,dim); return
                    end if
                end do
            end if
            contains
            !low & up must satisfy:
            !    low < up
            !    low satisfies sufficient decrease condition, but up violates
            !    phi'(low) < 0
            subroutine zoom(low,up,flow,fup,phidlow)!Restore x and fx to input status before calling zoom
                real*8,intent(inout)::low,up,flow,fup,phidlow
                real*8::phidnew,phidlow_m_a
                phidlow_m_a=phidlow*a!Initialize
                do
                    !Updata a by quadratic interpolation
                        a=phidlow_m_a*a/2d0/(flow+phidlow_m_a-fup)
                        if(.not.(a>low.and.a<up)) a=(low+up)/2d0!Fail safe
                    x=x0+a*p; call f(fx,x,dim)
                    if(fx>fx0+c1*a*phid0) then
                        up=a
                        if(up-low<1d-15.or.(up-low)/max(dAbs(low),dAbs(up))<1d-15) then
                            call fd(fdx,x,dim); return
                        end if
                        fup=fx
                    else
                        call fd(fdx,x,dim); phidnew=dot_product(fdx,p)
                        if(phidnew>c2_m_abs_phid0) return
                        low=a
                        if(up-low<1d-15.or.(up-low)/max(dAbs(low),dAbs(up))<1d-15) return
                        flow=fx; phidlow=phidnew; phidlow_m_a=phidlow*a
                    end if
                end do
            end subroutine zoom
        end subroutine Wolfe_fdwithf
        
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
            if(present(Increment)) then; incrmt=max(1d0+1d-15,Increment)!Fail safe
            else; incrmt=1.05d0; end if
            x0=x; fx0=fx; c2_m_abs_phid0=c2*dAbs(phid0)
            !Check whether initial guess satisfies sufficient decrease condition
            x=x0+a*p; call f(fx,x,dim)
            if(fx<=fx0+c1*a*phid0) then!Satisfied, try to search for larger a
                call fd(fdx,x,dim)
                phidnew=dot_product(fdx,p)
                if(phidnew>0d0) then!Curve is heading up
                    if(dAbs(phidnew)<=c2_m_abs_phid0) return
                    do!Else have to search for smaller a that phi(a) < phi(aold) & phid(a) > 0 is false
                        aold=a; fold=fx; phidold=phidnew
                        a=aold/incrmt; x=x0+a*p; call f(fx,x,dim); call fd(fdx,x,dim); phidnew=dot_product(fdx,p)
                        if(fx>=fold.or.phidnew<=0d0) then
                            atemp=a; ftemp=fx
                            call zoom(aold,atemp,fold,ftemp,phidold,phidnew)
                            return
                        end if
                        if(a<1d-15) return
                    end do
                else!Search for larger a
                    do
                        aold=a; fold=fx; phidold=phidnew
                        a=aold*incrmt; x=x0+a*p; call f(fx,x,dim); call fd(fdx,x,dim); phidnew=dot_product(fdx,p)
                        if(fx>fx0+c1*a*phid0.or.fx>=fold) then
                            atemp=a; ftemp=fx
                            call zoom(aold,atemp,fold,ftemp,phidold,phidnew)
                            return
                        end if
                        if(phidnew>0d0) then
                            if(dAbs(phidnew)<=c2_m_abs_phid0) then; return
                            else
                                atemp=a; ftemp=fx
                                call zoom(atemp,aold,ftemp,fold,phidnew,phidold)
                                fx=fx0
                            end if
                        end if
                    end do
                end if
            else!Violated, 1st search for smaller a satisfying sufficient decrease condition
                do
                    aold=a; fold=fx
                    a=aold/incrmt; x=x0+a*p; call f(fx,x,dim)
                    if(fx<=fx0+c1*a*phid0) then!Found, then look at slope
                        call fd(fdx,x,dim)
                        phidnew=dot_product(fdx,p)
                        if(dAbs(phidnew)<=c2_m_abs_phid0) return
                        if(phidnew<0d0) then!Within [a, aold]
                            x=x0+aold*p; call fd(fdx,x,dim); phidold=dot_product(fdx,p)
                            atemp=a; ftemp=fx
                            call zoom(atemp,aold,ftemp,fold,phidnew,phidold)
                            return
                        else!Search for such an a that phi(a) < phi(aold) & phid(a) > 0 is false
                            do
                                aold=a; fold=fx; phidold=phidnew
                                a=aold/incrmt; x=x0+a*p; call f(fx,x,dim); call fd(fdx,x,dim); phidnew=dot_product(fdx,p)
                                if(fx>=fold.or.phidnew<=0d0) then
                                    atemp=a; ftemp=fx
                                    call zoom(aold,atemp,fold,ftemp,phidold,phidnew)
                                    return
                                end if
                                if(a<1d-15) return
                            end do
                        end if
                    end if
                    if(a<1d-15) then
                        call fd(fdx,x,dim); return
                    end if
                end do
            end if
            contains
            !low & up must satisfy:
            !    low satisfies sufficient decrease condition
            !    ( up - low ) * phi'(low) < 0
            !    [ low, up ] (or [ up, low ]) contains a step length satisfying strong Wolfe condition
            !        This means at least 1 of 3 following statements is true:
            !            up violates the sufficient decrease condition
            !            phi(up) >= phi(low)
            !            up < low & phi'(up) <= 0
            subroutine zoom(low,up,flow,fup,phidlow,phidup)
                real*8,intent(inout)::low,up,flow,fup,phidlow,phidup
                real*8::phidnew,d1,d2
                do
                    !Updata a by cubic interpolation
                        d1=phidlow+phidup-3d0*(flow-fup)/(low-up); d2=up-low
                        if(d2>0d0) then; d2=dSqrt(d1*d1-phidlow*phidup)
                        else; d2=-dSqrt(d1*d1-phidlow*phidup); end if
                        a=up-(up-low)*(phidup+d2-d1)/(phidup-phidlow+2d0*d2)
                        if(.not.(a>min(low,up).and.a<max(low,up))) a=(low+up)/2d0!Fail safe
                    x=x0+a*p; call f(fx,x,dim); call fd(fdx,x,dim); phidnew=dot_product(fdx,p)
                    if(fx>fx0+c1*a*phid0.or.fx>=flow) then
                        up=a; fup=fx; phidup=phidnew
                    else
                        if(dAbs(phidnew)<=c2_m_abs_phid0) return
                        if(phidnew*(up-low)>=0d0) then
                            up=low; fup=flow; phidup=phidlow
                        end if
                        low=a; flow=fx; phidlow=phidnew
                    end if
                    if(dAbs(up-low)<1d-15.or.dAbs(up-low)/max(dAbs(low),dAbs(up))<1d-15) return
                end do
            end subroutine zoom
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
                if(present(Increment)) then; incrmt=max(1d0+1d-15,Increment)!Fail safe
                else; incrmt=1.05d0; end if
                x0=x; fx0=fx; c2_m_abs_phid0=c2*dAbs(phid0)
            !Check whether initial guess satisfies sufficient decrease condition
            x=x0+a*p; info=f_fd(fx,fdx,x,dim)
            if(fx<=fx0+c1*a*phid0) then!Satisfied, try to search for larger a
                phidnew=dot_product(fdx,p)
                if(phidnew>0d0) then!Curve is heading up
                    if(dAbs(phidnew)<=c2_m_abs_phid0) return
                    do!Else have to search for smaller a that phi(a) < phi(aold) & phid(a) > 0 is false
                        aold=a; fold=fx; phidold=phidnew
                        a=aold/incrmt; x=x0+a*p; info=f_fd(fx,fdx,x,dim); phidnew=dot_product(fdx,p)
                        if(fx>=fold.or.phidnew<=0d0) then
                            atemp=a; ftemp=fx
                            call zoom(aold,atemp,fold,ftemp,phidold,phidnew)
                            return
                        end if
                        if(a<1d-15) return
                    end do
                else!Search for larger a
                    do
                        aold=a; fold=fx; phidold=phidnew
                        a=aold*incrmt; x=x0+a*p; info=f_fd(fx,fdx,x,dim); phidnew=dot_product(fdx,p)
                        if(fx>fx0+c1*a*phid0.or.fx>=fold) then
                            atemp=a; ftemp=fx
                            call zoom(aold,atemp,fold,ftemp,phidold,phidnew)
                            return
                        end if
                        if(phidnew>0d0) then
                            if(dAbs(phidnew)<=c2_m_abs_phid0) return
                            atemp=a; ftemp=fx
                            call zoom(atemp,aold,ftemp,fold,phidnew,phidold)
                            return
                        end if
                    end do
                end if
            else!Violated, 1st search for smaller a satisfying sufficient decrease condition
				do
                    aold=a; fold=fx
                    a=aold/incrmt; x=x0+a*p; call f(fx,x,dim)
                    if(fx<=fx0+c1*a*phid0) then!Found, then look at slope
                        call fd(fdx,x,dim); phidnew=dot_product(fdx,p)
                        if(dAbs(phidnew)<=c2_m_abs_phid0) return
                        if(phidnew<0d0) then!Within [a, aold]
                            x=x0+aold*p; call fd(fdx,x,dim); phidold=dot_product(fdx,p)
                            atemp=a; ftemp=fx
                            call zoom(atemp,aold,ftemp,fold,phidnew,phidold)
                            return
                        else!Search for such an a that phi(a) < phi(aold) & phid(a) > 0 is false
                            do
                                aold=a; fold=fx; phidold=phidnew
                                a=aold/incrmt; x=x0+a*p; info=f_fd(fx,fdx,x,dim); phidnew=dot_product(fdx,p)
                                if(fx>=fold.or.phidnew<=0d0) then
                                    atemp=a; ftemp=fx
                                    call zoom(aold,atemp,fold,ftemp,phidold,phidnew)
                                    return
                                end if
                                if(a<1d-15) return
                            end do
                        end if
                    end if
                    if(a<1d-15) then
                        call fd(fdx,x,dim); return
                    end if
                end do
            end if
            contains
            !low & up must satisfy:
            !    low satisfies sufficient decrease condition
            !    ( up - low ) * phi'(low) < 0
            !    [ low, up ] (or [ up, low ]) contains a step length satisfying strong Wolfe condition
            !        This means at least 1 of 3 following statements is true:
            !            up violates the sufficient decrease condition
            !            phi(up) >= phi(low)
            !            up < low & phi'(up) <= 0
            subroutine zoom(low,up,flow,fup,phidlow,phidup)
                real*8,intent(inout)::low,up,flow,fup,phidlow,phidup
                real*8::phidnew,d1,d2
				do
                    !Updata a by cubic interpolation
                        d1=phidlow+phidup-3d0*(flow-fup)/(low-up); d2=up-low
                        if(d2>0d0) then; d2=dSqrt(d1*d1-phidlow*phidup)
                        else; d2=-dSqrt(d1*d1-phidlow*phidup); end if
                        a=up-(up-low)*(phidup+d2-d1)/(phidup-phidlow+2d0*d2)
                        if(.not.(a>min(low,up).and.a<max(low,up))) a=(low+up)/2d0!Fail safe
                    x=x0+a*p; info=f_fd(fx,fdx,x,dim); phidnew=dot_product(fdx,p)
                    if(fx>fx0+c1*a*phid0.or.fx>=flow) then
                        up=a; fup=fx; phidup=phidnew
                    else
                        if(dAbs(phidnew)<=c2_m_abs_phid0) return
                        if(phidnew*(up-low)>=0d0) then
                            up=low; fup=flow; phidup=phidlow
                        end if
                        low=a; flow=fx; phidlow=phidnew
                    end if
                    if(dAbs(up-low)<1d-15.or.dAbs(up-low)/max(dAbs(low),dAbs(up))<1d-15) return
                end do
            end subroutine zoom
        end subroutine StrongWolfe_fdwithf
    !================ End =================
!------------------ End ------------------

!------------- Trust region --------------
    !MKL trust-region nonlinear least square problem (trnlsp) solver wrapper
    !Solve f'(x) = 0 by minimizing merit function F(x) = f'(x)^2 through trust-region method
    !with model function m(p) = [ f'(x) + J(x) . p ]^2, where J(x) is the Jacobian
    !This algorithm can be interpreted in different ways other than f'(x) = 0:
    !    f'(x) can be viewed as Sqrt(weight)*residual terms,
    !        then trnlsp minimizes the square penalty (this is where its name comes from)
    !    When M = N, f'(x) can also be considered as the gradient of f(x), Jacobian = Hessian,
    !        but trnlsp doesn't necessarily optimize f(x) (unless f(x) = const * F(x) by coincidence),
    !        it merely return a stationary point of f(x)
    !External procedure:
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
                if(present(Warning)) then; warn=Warning
                    else; warn=.true.; end if
                if(present(MaxIteration)) then; maxit=MaxIteration
                    else; maxit=1000; end if
                if(present(MaxStepIteration)) then; maxstepit=MaxStepIteration
                    else; maxstepit=100; end if
                tol=[1d-15,1d-15,1d-15,1d-15,1d-15,1d-15]
                if(present(Precision)) tol(2)=Precision
                if(present(MinStepLength)) then
                    tol(1)=MinStepLength
                    tol(4)=MinStepLength
                    tol(5)=MinStepLength
                end if
            fdx=0d0; J=0d0; StepBound=100d0; RCI_request=0
        if(present(low).and.present(up)) then
            if(dtrnlspbc_init(handle,N,M,x,low,up,tol,maxit,maxstepit,StepBound)/=TR_SUCCESS) then
                write(*,*)'Trust region abort: invalid initialization'
                call mkl_free_buffers; return
            end if
            if(dtrnlspbc_check(handle,N,M,J,fdx,low,up,tol,info)/=TR_SUCCESS) then
                write(*,*)'Trust region abort: check failed'
                call mkl_free_buffers; return
            else
                if(info(1)/=0.or.info(2)/=0.or.info(3)/=0.or.info(4)/=0.or.info(5)/=0.or.info(6)/=0) then
                    write(*,*)'Trust region abort: check was not passed, the information is:'
                    write(*,*)info
                    call mkl_free_buffers; return
                end if
            end if
            if(present(Jacobian)) then
                do!Main loop
                    if (dtrnlspbc_solve(handle,fdx,J,RCI_request)/=TR_SUCCESS) then
                        call mkl_free_buffers; return
                    end if
                    select case (RCI_request)
                        case (-1,-2,-3,-4,-5,-6); exit
                        case (1); call fd(fdx,x,M,N)
                        case (2); i=Jacobian(J,x,M,N)
                    end select
                end do
            else
                do!Main loop
                    if (dtrnlspbc_solve(handle,fdx,J,RCI_request)/=TR_SUCCESS) then
                        call mkl_free_buffers; return
                    end if
                    select case (RCI_request)
                        case (-1,-2,-3,-4,-5,-6); exit
                        case (1); call fd(fdx,x,M,N)
                        case (2)
                            if(djacobi(fd_j,N,M,J,x,1d-8)/=TR_SUCCESS) then
                                call mkl_free_buffers; return
                            end if
                    end select
                end do
            end if
            !Clean up
                if (dtrnlspbc_get(handle,TotalIteration,StopReason,InitialResidual,FinalResidual)/=TR_SUCCESS) then
                    call mkl_free_buffers; return
                end if
                if (dtrnlspbc_delete(handle)/=TR_SUCCESS) then
                    call mkl_free_buffers; return
                end if
            do i=1,N
                if((x(i)<low(i).or.x(i)>up(i)).and.warn) then
                    write(*,'(1x,A49)')'Failed trust region: boundary condition violated!'
                    exit
                end if
            end do
        else
            if(dtrnlsp_init(handle,N,M,x,tol,maxit,maxstepit,StepBound)/=TR_SUCCESS) then
                write(*,*)'Trust region abort: invalid initialization'
                call mkl_free_buffers; return
            end if
            if(dtrnlsp_check(handle,N,M,J,fdx,tol,info)/=TR_SUCCESS) then
                write(*,*)'Trust region abort: check failed'
                call mkl_free_buffers; return
            else
                if(info(1)/=0.or.info(2)/=0.or.info(3)/=0.or.info(4)/=0) then
                    write(*,*)'Trust region abort: check was not passed, the information is:'
                    write(*,*)info
                    call mkl_free_buffers; return
                end if
            end if
            if(present(Jacobian)) then
                do!Main loop
                    if (dtrnlsp_solve(handle,fdx,J,RCI_request)/=TR_SUCCESS) then
                        call mkl_free_buffers; return
                    end if
                    select case (RCI_request)
                        case (-1,-2,-3,-4,-5,-6); exit
                        case (1); call fd(fdx,x,M,N)
                        case (2); i=Jacobian(J,x,M,N)
                    end select
                end do
            else
                do!Main loop
                    if (dtrnlsp_solve(handle,fdx,J,RCI_request)/=TR_SUCCESS) then
                        call mkl_free_buffers; return
                    end if
                    select case (RCI_request)
                        case (-1,-2,-3,-4,-5,-6); exit
                        case (1); call fd(fdx,x,M,N)
                        case (2)
                            if(djacobi(fd_j,N,M,J,x,1d-8)/=TR_SUCCESS) then
                                call mkl_free_buffers; return
                            end if
                    end select
                end do
            end if
            !Clean up
                if (dtrnlsp_get(handle,TotalIteration,StopReason,InitialResidual,FinalResidual)/=TR_SUCCESS) then
                    call mkl_free_buffers; return
                end if
                if (dtrnlsp_delete(handle)/=TR_SUCCESS) then
                    call mkl_free_buffers; return
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
                call fd(fdx,x,M,N)
            end subroutine fd_j
    end subroutine TrustRegion
!------------------ End ------------------

!---------- Equality constraint ----------
    !Lagrangian multiplier method is a classical way to treat equality constraint:
    !    L = f - lamda . c
    !where f is the target function to be minimized, lamda is Lagrangian multiplier, c is equality constraint c(x) = 0
    !Textbook is wrong: it claims Lagrangian multiplier method transforms constrained optimization into unconstrained one
    !However, L has no lower bound, since lamda . c can approach infinity when c != 0 and lamda diverges
    !Lagrangian multiplier method actually turns a minimization problem into a saddle point problem,
    !which cannot necessarily be solved through decreasing L, deteriorating all unconstrained minimizers
    !Lagrangian multiplier method is numerically feasible only when at least 1 of the following statements is true:
    !    L has unique saddle point
    !    The initial guess is sufficiently close to the exact solution
    !under which circumstance we may simply minimize || L'(x) ||
    !In general case, we have to turn to the augmented Lagrangian method:
    !    Augmented Lagrangian = f - lamda . c + miu / 2 * c . c
    !where miu is constraint violation penalty strength
    !    miu should >= 1 because c ~ ( lamda - lamda_true ) / miu
    !Suggestion:
    !    Augmented Lagrangian has ill conditioned Hessian when miu is too large, deteriorating performance of line searchers
    !    so do not take too much interations nor push accuracy to double precision limit
    !External procedure:
    !    subroutine f(f(x),x,N)
    !    subroutine fd(f'(x),x,N)
    !    integer function f_fd(f(x),f'(x),x,N)
    !    integer function fdd(f''(x),x,N)
    !    subroutine c(c(x),x,M,N)
    !    subroutine cd(c'(x),x,M,N)
    !    integer function cdd(c''(x),x,M,N)
    !    N dimensional vector x & f'(x), N order matrix f''(x),
    !    M dimensional vector c, N x M matrix c'(x), N x N x M 3rd-order tensor c''(x)
    !Common required argument:
    !    subroutine f & fd & c & cd, N dimensional vector x, integer N & M
    !On input x is an initial guess, on exit x is a local minimum of f(x) subject to constraint c(x)

    !Lagrangian multiplier method for equality constraint
    !Please read the instruction above to make sure this is really feasible for your problem
    !Additional required argument:
    !    lamda: on input is an initial guess of Lagranguan multiplier, on exit is the solution
    !Optional argument:
    !    Warning: (default = true) if false, all warnings will be suppressed
    !    MaxIteration: (default = 1000) max number of iterations to perform
    !    Precision: (default = 1d-15) convergence considered when || L'(x) ||_2 < Precision
    subroutine LagrangianMultiplier(fd,fdd,c,cd,cdd,x,lamda,N,M,&
        Warning,MaxIteration,Precision)
        !Required argument
            external::fd,c,cd; integer,external::fdd,cdd
            integer,intent(in)::N,M
            real*8,dimension(N),intent(inout)::x; real*8,dimension(M),intent(inout)::lamda
        !Optional argument
            logical,intent(in),optional::Warning
            integer,intent(in),optional::MaxIteration
            real*8,intent(in),optional::Precision
        !Job control
            logical::warn; integer::maxit; real*8::tol
        integer::iIteration,i,dim
        real*8,dimension(N+M)::minusLd; real*8,dimension(N+M,N+M)::Ldd
        real*8,dimension(M)::cx; real*8,dimension(N,M)::cdx; real*8,dimension(N,N,M)::cddx
        !Set parameter according to optional argument
            if(present(Warning)) then; warn=Warning
                else; warn=.true.; end if
            if(present(MaxIteration)) then; maxit=MaxIteration
                else; maxit=1000; end if
            if(present(Precision)) then; tol=Precision*Precision!To save sqrt cost, precision is squared
                else; tol=1d-30; end if
        dim=N+M
        do iIteration=1,maxit
            !Construct -Ld and Ldd
                call fd(minusLd(1:N),x,N); call c(cx,x,M,N); call cd(cdx,x,M,N)
                minusLd(1:N)=matmul(cdx,lamda)-minusLd(1:N); minusLd(N+1:dim)=cx
                if(dot_product(minusLd,minusLd)<tol) return!Converged
                i=fdd(Ldd(1:N,1:N),x,N); i=cdd(cddx,x,M,N)
                forall(i=1:N)
                    Ldd(i,1:N)=Ldd(i,1:N)-matmul(cddx(i,:,:),lamda)
                end forall
                Ldd(N+1:dim,1:N)=-transpose(cdx); Ldd(N+1:dim,N+1:dim)=0d0
            !Newton iteration
                call My_dsysv(Ldd,minusLd,dim)
                x=x+minusLd(1:N); lamda=lamda+minusLd(N+1:dim)
        end do
        if(iIteration>maxit.and.warn) then
            write(*,'(1x,A53)')'Failed Lagrangian multiplier: max iteration exceeded!'
            call fd(minusLd(1:N),x,N); call c(cx,x,M,N); call cd(cdx,x,M,N)
            minusLd(1:N)=matmul(cdx,lamda)-minusLd(1:N); minusLd(N+1:dim)=cx
            write(*,*)'Euclidean norm of Lagrangian gradient =',Norm2(minusLd)
        end if
    end subroutine LagrangianMultiplier

    !Augmented Lagrangian method for equality constraint
    !Optional argument:
    !    UnconstrainedSolver: (default = BFGS) specify the unconstraind solver to use, every line searcher is available
    !    lamda0: (default = 0) initial guess of lamda
    !    miu0: (default = 1) initial miu (must >= 1)
    !    All line search optional arguments are also optional here and will be passed to line searchers,
	!    some of them also control augmented Lagrangian behaviour: Warning, MaxIteration, Precision, Increment
    !        MaxIteration: max number of augmented Lagrangian iterations to perform
    !        Precision: convergence considered when || c(x) ||_2 < Precision
	!        Increment: each iteration change miu by how much time
    subroutine AugmentedLagrangian(f,fd,c,cd,x,N,M,UnconstrainedSolver,lamda0,miu0,&
        f_fd,Strong,Warning,MaxIteration,Precision,MinStepLength,WolfeConst1,WolfeConst2,Increment,fdd,cdd,ExactStep,Memory,Method)
        !Required argument
            external::f,fd,c,cd
            integer,intent(in)::N,M
            real*8,dimension(N),intent(inout)::x
        !Optional argument
            character*32,intent(in),optional::UnconstrainedSolver
            real*8,dimension(M),intent(in),optional::lamda0
            integer,external,optional::f_fd,fdd,cdd
            logical,intent(in),optional::Strong,Warning
            integer,intent(in),optional::MaxIteration,ExactStep,Memory
            real*8,intent(in),optional::miu0,Precision,MinStepLength,WolfeConst1,WolfeConst2,Increment
            character*32,intent(in),optional::Method
        !Job control
            character*32::solver,type
            logical::sw,warn
            integer::maxit,freq,mem
            real*8::tol,minstep,c1,c2,incrmt
        integer::iIteration,i
        real*8::tolsq,miu
        real*8,dimension(M)::lamda,cx
        real*8,dimension(N,M)::cdx
        real*8,dimension(N,N,M)::cddx
        real*8,dimension(N,N)::Lddxtemp
        !Set parameter according to optional argument
            !Augmented Lagrangian optional argument
                if(present(UnconstrainedSolver)) then; solver=UnconstrainedSolver
                    else; solver='BFGS'; end if
                if(present(lamda0)) then; lamda=lamda0
                    else; lamda=0d0; end if
                if(present(miu0)) then; miu=max(1d0,miu0)
                    else; miu=1d0; end if
            !Common line search optional argument
                if(present(Strong)) then; sw=Strong
                    else; sw=.true.; end if
                if(present(Warning)) then; warn=Warning
                    else; warn=.true.; end if
                if(present(MaxIteration)) then; maxit=MaxIteration
                    else; maxit=1000; end if
                if(present(Precision)) then; tol=Precision
                    else; tol=1d-15; end if
                if(present(MinStepLength)) then; minstep=MinStepLength
                    else;  minstep=1d-15; end if
                if(present(WolfeConst1)) then; c1=max(1d-15,WolfeConst1)!Fail safe
                    else; c1=1d-4; end if
                if(present(WolfeConst2)) then
                    c2=min(1d0-1d-15,max(c1+1d-15,WolfeConst2))!Fail safe
                else
                    if(present(UnconstrainedSolver)) then
                        if(UnconstrainedSolver=='ConjugateGradient') then; c2=0.45d0
                            else; c2=0.9d0; end if
                    else
                        c2=0.9d0
                    end if
                end if
                if(present(Increment)) then; incrmt=Increment
                    else; incrmt=1.05d0; end if
            !Solver specific optional argument
                if(present(ExactStep)) then; freq=ExactStep
                    else; freq=20; end if
                if(present(Memory)) then; mem=max(1,Memory)!Fail safe
                    else; mem=10; end if
                if(present(Method)) then; type=Method
                    else; type='DY'; end if
		tolsq=tol*tol!To save sqrt cost, precision is squared
        select case(solver)
            case('NewtonRaphson')
                if(present(fdd).and.present(cdd)) then
                    if(present(f_fd)) then
                        do iIteration=1,maxit
                            call NewtonRaphson(L,Ld,x,N,f_fd=L_Ld_fdwithf,fdd=Ldd,&
                            Strong=sw,Warning=warn,MaxIteration=maxit,Precision=tol,MinStepLength=minstep,WolfeConst1=c1,WolfeConst2=c2,Increment=incrmt)
                            call c(cx,x,M,N)
                            if(dot_product(cx,cx)<tolsq) exit
                            lamda=lamda-miu*cx; miu=miu*incrmt
                        end do
                    else
                        do iIteration=1,maxit
                            call NewtonRaphson(L,Ld,x,N,f_fd=L_Ld,fdd=Ldd,&
                            Strong=sw,Warning=warn,MaxIteration=maxit,Precision=tol,MinStepLength=minstep,WolfeConst1=c1,WolfeConst2=c2,Increment=incrmt)
                            call c(cx,x,M,N)
                            if(dot_product(cx,cx)<tolsq) exit
                            lamda=lamda-miu*cx; miu=miu*incrmt
                        end do
                    end if
                else
                    if(present(f_fd)) then
                        do iIteration=1,maxit
                            call NewtonRaphson(L,Ld,x,N,f_fd=L_Ld_fdwithf,&
                            Strong=sw,Warning=warn,MaxIteration=maxit,Precision=tol,MinStepLength=minstep,WolfeConst1=c1,WolfeConst2=c2,Increment=incrmt)
                            call c(cx,x,M,N)
                            if(dot_product(cx,cx)<tolsq) exit
                            lamda=lamda-miu*cx; miu=miu*incrmt
                        end do
                    else
						do iIteration=1,maxit
                            call NewtonRaphson(L,Ld,x,N,f_fd=L_Ld,&
                            Strong=sw,Warning=warn,MaxIteration=maxit,Precision=tol,MinStepLength=minstep,WolfeConst1=c1,WolfeConst2=c2,Increment=incrmt)
                            call c(cx,x,M,N)
                            if(dot_product(cx,cx)<tolsq) exit
                            lamda=lamda-miu*cx; miu=miu*incrmt
                        end do
                    end if
                end if
            case('BFGS')
                if(present(fdd).and.present(cdd)) then
                    if(present(f_fd)) then
                        do iIteration=1,maxit
                            call BFGS(L,Ld,x,N,f_fd=L_Ld_fdwithf,fdd=Ldd,ExactStep=freq,&
                            Strong=sw,Warning=warn,MaxIteration=maxit,Precision=tol,MinStepLength=minstep,WolfeConst1=c1,WolfeConst2=c2,Increment=incrmt)
                            call c(cx,x,M,N)
                            if(dot_product(cx,cx)<tolsq) exit
                            lamda=lamda-miu*cx; miu=miu*incrmt
                        end do
                    else
                        do iIteration=1,maxit
                            call BFGS(L,Ld,x,N,f_fd=L_Ld,fdd=Ldd,ExactStep=freq,&
                            Strong=sw,Warning=warn,MaxIteration=maxit,Precision=tol,MinStepLength=minstep,WolfeConst1=c1,WolfeConst2=c2,Increment=incrmt)
                            call c(cx,x,M,N)
                            if(dot_product(cx,cx)<tolsq) exit
                            lamda=lamda-miu*cx; miu=miu*incrmt
                        end do
                    end if
                else
                    if(present(f_fd)) then
                        do iIteration=1,maxit
                            call BFGS(L,Ld,x,N,f_fd=L_Ld_fdwithf,ExactStep=freq,&
                            Strong=sw,Warning=warn,MaxIteration=maxit,Precision=tol,MinStepLength=minstep,WolfeConst1=c1,WolfeConst2=c2,Increment=incrmt)
                            call c(cx,x,M,N)
                            if(dot_product(cx,cx)<tolsq) exit
                            lamda=lamda-miu*cx; miu=miu*incrmt
                        end do
                    else
                        do iIteration=1,maxit
                            call BFGS(L,Ld,x,N,f_fd=L_Ld,ExactStep=freq,&
                            Strong=sw,Warning=warn,MaxIteration=maxit,Precision=tol,MinStepLength=minstep,WolfeConst1=c1,WolfeConst2=c2,Increment=incrmt)
                            call c(cx,x,M,N)
                            if(dot_product(cx,cx)<tolsq) exit
                            lamda=lamda-miu*cx; miu=miu*incrmt
                        end do
                    end if
                end if
            case('LBFGS')
                if(present(f_fd)) then
                    do iIteration=1,maxit
                        call LBFGS(L,Ld,x,N,f_fd=L_Ld_fdwithf,Memory=mem,&
                        Strong=sw,Warning=warn,MaxIteration=maxit,Precision=tol,MinStepLength=minstep,WolfeConst1=c1,WolfeConst2=c2,Increment=incrmt)
                        call c(cx,x,M,N)
                        if(dot_product(cx,cx)<tolsq) exit
                        lamda=lamda-miu*cx; miu=miu*incrmt
                    end do
                else
                    do iIteration=1,maxit
                        call LBFGS(L,Ld,x,N,f_fd=L_Ld,Memory=mem,&
                        Strong=sw,Warning=warn,MaxIteration=maxit,Precision=tol,MinStepLength=minstep,WolfeConst1=c1,WolfeConst2=c2,Increment=incrmt)
                        call c(cx,x,M,N)
                        if(dot_product(cx,cx)<tolsq) exit
                        lamda=lamda-miu*cx; miu=miu*incrmt
                    end do
                end if
            case('ConjugateGradient')
                if(present(f_fd)) then
                    do iIteration=1,maxit
                        call ConjugateGradient(L,Ld,x,N,f_fd=L_Ld_fdwithf,Method=type,&
                        Strong=sw,Warning=warn,MaxIteration=maxit,Precision=tol,MinStepLength=minstep,WolfeConst1=c1,WolfeConst2=c2,Increment=incrmt)
                        call c(cx,x,M,N)
                        if(dot_product(cx,cx)<tolsq) exit
                        lamda=lamda-miu*cx; miu=miu*incrmt
                    end do
                else
                    do iIteration=1,maxit
                        call ConjugateGradient(L,Ld,x,N,f_fd=L_Ld,Method=type,&
                        Strong=sw,Warning=warn,MaxIteration=maxit,Precision=tol,MinStepLength=minstep,WolfeConst1=c1,WolfeConst2=c2,Increment=incrmt)
                        call c(cx,x,M,N)
                        if(dot_product(cx,cx)<tolsq) exit
                        lamda=lamda-miu*cx; miu=miu*incrmt
                    end do
                end if
            case default!Throw a warning
                write(*,'(1x,A47,1x,A2)')'Program abort: unsupported unconstrained solver',solver
                stop
        end select
        if(iIteration>maxit.and.warn) then
            write(*,'(1x,A52)')'Failed augmented Lagrangian: max iteration exceeded!'
            write(*,*)'Euclidean norm of constraint violation =',Norm2(cx)
        end if
        contains
            subroutine L(Lx,x,N)
                integer,intent(in)::N
                real*8,dimension(N),intent(in)::x
                real*8,intent(out)::Lx
                call f(Lx,x,N); call c(cx,x,M,N)
                Lx=Lx-dot_product(lamda,cx)+miu/2d0*dot_product(cx,cx)
            end subroutine L
            subroutine Ld(Ldx,x,N)
                integer,intent(in)::N
                real*8,dimension(N),intent(in)::x
                real*8,dimension(N),intent(out)::Ldx
                call fd(Ldx,x,N); call c(cx,x,M,N); call cd(cdx,x,M,N)
                Ldx=Ldx+matmul(cdx,miu*cx-lamda)
            end subroutine Ld
            integer function L_Ld(Lx,Ldx,x,N)!Compute c & cd together is cheaper
                integer,intent(in)::N
                real*8,dimension(N),intent(in)::x
                real*8,intent(out)::Lx
                real*8,dimension(N),intent(out)::Ldx
                call f(Lx,x,N); call c(cx,x,M,N)
                Lx=Lx-dot_product(lamda,cx)+miu/2d0*dot_product(cx,cx)
                call fd(Ldx,x,N); call cd(cdx,x,M,N)
                Ldx=Ldx+matmul(cdx,miu*cx-lamda)
                L_Ld=0!return 0
            end function L_Ld
            integer function L_Ld_fdwithf(Lx,Ldx,x,N)!When f_fd is available
                integer,intent(in)::N
                real*8,dimension(N),intent(in)::x
                real*8,intent(out)::Lx
                real*8,dimension(N),intent(out)::Ldx
                i=f_fd(Lx,Ldx,x,N); call c(cx,x,M,N)
                Lx=Lx-dot_product(lamda,cx)+miu/2d0*dot_product(cx,cx)
                call cd(cdx,x,M,N)
                Ldx=Ldx+matmul(cdx,miu*cx-lamda)
                L_Ld_fdwithf=0!return 0
            end function L_Ld_fdwithf
            integer function Ldd(Lddx,x,N)!When fdd & cdd is available
                integer,intent(in)::N
                real*8,dimension(N),intent(in)::x
                real*8,dimension(N,N),intent(out)::Lddx
                i=fdd(Lddx,x,N); i=cdd(cddx,x,M,N); call c(cx,x,M,N); call cd(cdx,x,M,N)
                cx=miu*cx-lamda
                forall(i=1:N)
                    Lddxtemp(:,i)=matmul(cddx(i,:,:),cx)
                end forall
                Lddx=Lddx+Lddxtemp+matmul(cdx,transpose(cdx))
                Ldd=0!return 0
            end function Ldd
    end subroutine AugmentedLagrangian
!------------------ End ------------------

!--------------- Heuristic ---------------
    !Not implemented, because I cannot figure out a way to implement these algorithms as black box
!------------------ End ------------------

end module NonlinearOptimization