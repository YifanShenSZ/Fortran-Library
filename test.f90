!A test program on Fortran-Library
program main
    use General
    use Mathematics
    use LinearAlgebra
    use GeometryTransformation
    use NonlinearOptimization
    use Nonadiabatic
    implicit none
    integer::i,j,k,l,freq,memory,dim,M,N
    integer,dimension(10)::indicesort
    real*8::difference,DoubleTemp,dbtp
    real*8,dimension(3)::avec,bvec,cvec,mass
    real*8,dimension(4)::quaternion
    real*8,dimension(9)::H2Ogeom,H2Ogeom1,H2Ogeom2
    real*8,dimension(10)::eigval,x,low,up
    real*8,dimension(10,10)::A,B,C,eigvec
    real*8,dimension(10,10,10)::A3,B3,C3
    real*8,dimension(10,10,10,10)::A4
    complex*16,dimension(10)::zeigval
    complex*16,dimension(10,10)::zA,zB,zeigvec

    call BetterRandomSeed()
    write(*,*)'This is a test program on Fortran-Library'
    write(*,*)'Correct routines should print close to 0'
    write(*,*)

write(*,*)'Testing all linear algebra routines...'
    write(*,*)
    write(*,*)'Testing vector operation...'
        write(*,*)
        write(*,*)'cross_product, triple_product'
            call random_number(avec)
            call random_number(bvec)
            call random_number(cvec)
            write(*,*)triple_product(avec,bvec,cvec)-dot_product(cross_product(avec,bvec),cvec)
    write(*,*)
    write(*,*)'Testing matrix operation...'
        write(*,*)
        write(*,*)'matmul_dgemm, matmul_dsymm'
            do j=1,10
                do i=j,10
                    call BetterRandomNumber(A(i,j))
                end do
            end do
            call syL2U(A,10)
            do j=1,10
                do i=j,10
                    call BetterRandomNumber(B(i,j))
                end do
            end do
            call syL2U(B,10)
            write(*,*)My_dlange('M',matmul(A,B)-matmul_dgemm(A,B,10,10,10),10,10)
            write(*,*)My_dlange('M',matmul(A,B)-matmul_dsymm(A,B,10,10),10,10)
        write(*,*)
        write(*,*)'symatmulasy, asymatmulsy'
            do j=1,10
                do i=j,10
                    call BetterRandomNumber(A(i,j))
                end do
            end do
            call syL2U(A,10)
            do i=1,10
                do j=1,10
                    B(i,j)=dFactorial2(i)-dFactorial2(j)
                end do
            end do
            write(*,*)My_dlange('M',matmul(A,B)-symatmulasy(A,B,10),10,10)
            write(*,*)My_dlange('M',matmul(B,A)-asymatmulsy(B,A,10),10,10)
    write(*,*)
    write(*,*)'Testing high order tensor operation...'
        write(*,*)
        write(*,*)'sy3matmulsy'
            do j=1,10
                do i=j,10
                    call random_number(A3(:,i,j))
                end do
            end do
            call sy3L2U(A3,10,10)
            do j=1,10
                do i=j,10
                    call BetterRandomNumber(B(i,j))
                end do
            end do
            call syL2U(B,10)
            B3=sy3matmulsy(A3,B,10,10)
            forall(i=1:10)
                C3(i,:,:)=matmul(A3(i,:,:),B)
            end forall
            DoubleTemp=0d0
            do j=1,10
                do i=1,10
                    DoubleTemp=DoubleTemp+dot_product(B3(:,i,j)-C3(:,i,j),B3(:,i,j)-C3(:,i,j))
                end do
            end do
            write(*,*)DoubleTemp
        write(*,*)
        write(*,*)'sy3matdotmul'
            do j=1,10
                do i=j,10
                    call random_number(A3(:,i,j))
                end do
            end do
            call sy3L2U(A3,10,10)
            do j=1,10
                do i=j,10
                    call random_number(B3(:,i,j))
                end do
            end do
            call sy3L2U(B3,10,10)
            A=0d0
            do i=1,10
                do j=1,10
                    do k=1,10
                        A(i,j)=A(i,j)+dot_product(A3(:,i,k),B3(:,k,j))
                    end do
                end do
            end do
            write(*,*)My_dlange('M',sy3matdotmul(A3,B3,10,10)-A,10,10)
        write(*,*)
        write(*,*)'sy4matdotmulsy3'
            do j=1,10
                do i=j,10
                    call random_number(A4(:,:,i,j))
                end do
            end do
            call sy4L2U(A4,10,10,10)
            do j=1,10
                do i=j,10
                    call random_number(B3(:,i,j))
                end do
            end do
            call sy3L2U(B3,10,10)
            A3=0d0
            do i=1,10
                do j=1,10
                    do k=1,10
                        forall(l=1:10)
                            A3(l,i,j)=A3(l,i,j)+dot_product(A4(l,:,i,k),B3(:,k,j))
                        end forall
                    end do
                end do
            end do
            B3=sy4matdotmulsy3(A4,B3,10,10,10)
            DoubleTemp=0d0
            do j=1,10
                do i=1,10
                    DoubleTemp=DoubleTemp+dot_product(A3(:,i,j)-B3(:,i,j),A3(:,i,j)-B3(:,i,j))
                end do
            end do
            write(*,*)DoubleTemp
        write(*,*)
        write(*,*)'sy3UnitaryTransformation'
            do j=1,10
                do i=j,10
                    call random_number(A3(:,i,j))
                end do
            end do
            call sy3L2U(A3,10,10)
            do j=1,10
                do i=j,10
                    call BetterRandomNumber(B(i,j))
                end do
            end do
            call My_dsyev('V',B,eigval,10)
            B3=sy3UnitaryTransformation(A3,B,10,10)
            forall(i=1:10)
                C3(i,:,:)=matmul(transpose(B),matmul(A3(i,:,:),B))
            end forall
            DoubleTemp=0d0
            do j=1,10
                do i=j,10
                    DoubleTemp=DoubleTemp+dot_product(B3(:,i,j)-C3(:,i,j),B3(:,i,j)-C3(:,i,j))
                end do
            end do
            write(*,*)DoubleTemp
    write(*,*)
    write(*,*)'Testing linear solver & eigensystem...'
        do j=1,10
            do i=1,10
                call BetterRandomNumber(A(i,j))
            end do
        end do
        A=matmul(transpose(A),A)
        do i=1,10
            call BetterRandomNumber(DoubleTemp)
            A(i,i)=A(i,i)+DoubleTemp+1d-6
        end do
        call random_number(low)
        B=A
        x=low
        call My_dgesv(B,x,10)
        write(*,*)
        write(*,*)'dgesv'
            B=A
            up=low
            call My_dgesv(B,up,10)
            write(*,*)norm2(up-x)
        write(*,*)
        write(*,*)'dgetri'
            B=A
            call My_dgetri(B,10)
            write(*,*)norm2(matmul(B,low)-x)
        write(*,*)
        write(*,*)'dgeev'
            B=A
            call My_dgeev('V',B,eigval,up,eigvec,10)
            up=matmul(eigvec,matmul(transpose(eigvec),low)/eigval)
            write(*,*)norm2(up-x)
        write(*,*)
        write(*,*)'dsysv'
            B=A
            up=low
            call My_dsysv(B,up,10)
            write(*,*)norm2(up-x)
        write(*,*)
        write(*,*)'dsytri'
            B=A
            call My_dsytri(B,10)
            call syL2U(B,10)
            write(*,*)norm2(matmul(B,low)-x)
        write(*,*)
        write(*,*)'dsyev'
            eigvec=A
            call My_dsyev('V',eigvec,eigval,10)
            up=matmul(eigvec,matmul(transpose(eigvec),low)/eigval)
            write(*,*)norm2(up-x)
        write(*,*)
        write(*,*)'dposv'
            B=A
            up=low
            call My_dposv(B,up,10)
            write(*,*)norm2(up-x)
        write(*,*)
        write(*,*)'dpotri'
            B=A
            call My_dpotri(B,10)
            call syL2U(B,10)
            write(*,*)norm2(matmul(B,low)-x)
        write(*,*)
        write(*,*)'zgeev, zheev'
            do j=1,10
                do i=1,10
                    call BetterRandomNumber(DoubleTemp)
                    call BetterRandomNumber(dbtp)
                    zA(i,j)=cmplx(DoubleTemp,dbtp)
                end do
            end do
            zA=matmul(conjg(transpose(zA)),zA)
            zB=zA
            call My_zgeev('V',zB,zeigval,zeigvec,10)
            low=abs(zeigval)
            forall(i=1:10)
                indicesort(i)=i
            end forall
            call dQuickSort(low,1,10,indicesort,10)
            forall(i=1:10)
                zB(:,i)=zeigvec(:,i)
            end forall
            call My_zheev('V',zA,eigval,10)
            write(*,*)norm2(low-eigval)
            write(*,*)z2normge(matmul(zB,zA),10,10)-1d0
        write(*,*)
    write(*,*)
    write(*,*)'Testing matrix norm...'
        write(*,*)
        write(*,*)'dlange, dlansy'
        do i=1,10
            do j=1,10
                A(i,j)=dFactorial2(i)*dFactorial2(j)
            end do
        end do
        write(*,*)My_dlange('M',A,10,10)-My_dlansy('M',A,10)
        write(*,*)My_dlange('F',A,10,10)-My_dlansy('F',A,10)
        write(*,*)My_dlange('1',A,10,10)-My_dlansy('1',A,10)
        write(*,*)My_dlange('I',A,10,10)-My_dlansy('I',A,10)
        write(*,*)
        write(*,*)'d2normge'
        eigvec=matmul(transpose(A),A)
        call My_dsyev('N',eigvec,eigval,10)
        write(*,*)d2normge(A,10,10)-dSqrt(maxval(eigval))
    write(*,*)
write(*,*)'Linear algebra routines test passed'
write(*,*)

write(*,*)'Testing all nonlinear-optimization solvers...'
    write(*,*)
    freq=10
    memory=10
    dim=10
    M=10
    N=10
    low=-1d0
    up=1d0
    write(*,*)'Newton'
        call random_number(x)
        NewtonRaphsonWarning=.true.
        call NewtonRaphson(f,fd,fdd,x,dim)
        write(*,*)norm2(x)
    write(*,*)
    write(*,*)'Newton_S'
        call random_number(x)
        NewtonRaphsonWarning=.true.
        call NewtonRaphson_Strong(f,fd,fdd,x,dim)
        write(*,*)norm2(x)
    write(*,*)
    write(*,*)'Newton_S_fdwithf'
        call random_number(x)
        NewtonRaphsonWarning=.true.
        call NewtonRaphson_Strong_fdwithf(f,fd,f_fd,fdd,x,dim)
        write(*,*)norm2(x)
    write(*,*)
    write(*,*)'BFGS'
       call random_number(x)
       QuasiNewtonWarning=.true.
       call BFGS(f,fd,fdd,x,dim,freq)
       write(*,*)norm2(x)
    write(*,*)
    write(*,*)'BFGS_cheap'
       call random_number(x)
       QuasiNewtonWarning=.true.
       call BFGS_cheap(f,fd,x,dim)
        write(*,*)norm2(x)
    write(*,*)
    write(*,*)'BFGS_NH'
        call random_number(x)
        QuasiNewtonWarning=.true.
        call BFGS_NumericalHessian(f,fd,fd_j,x,dim,freq)
        write(*,*)norm2(x)
    write(*,*)
    write(*,*)'BFGS_S'
        call random_number(x)
        QuasiNewtonWarning=.true.
        call BFGS_Strong(f,fd,fdd,x,dim,freq)
        write(*,*)norm2(x)
    write(*,*)
    write(*,*)'BFGS_S_cheap'
        call random_number(x)
        QuasiNewtonWarning=.true.
        call BFGS_Strong_cheap(f,fd,x,dim)
        write(*,*)norm2(x)
    write(*,*)
    write(*,*)'BFGS_S_NH'
        call random_number(x)
        QuasiNewtonWarning=.true.
        call BFGS_Strong_NumericalHessian(f,fd,fd_j,x,dim,freq)
        write(*,*)norm2(x)
    write(*,*)
    write(*,*)'BFGS_S_fdwithf'
        call random_number(x)
        QuasiNewtonWarning=.true.
        call BFGS_Strong_fdwithf(f,fd,f_fd,fdd,x,dim,freq)
        write(*,*)norm2(x)
    write(*,*)
    write(*,*)'BFGS_S_cheap_fdwithf'
        call random_number(x)
        QuasiNewtonWarning=.true.
        call BFGS_Strong_cheap_fdwithf(f,fd,f_fd,x,dim)
        write(*,*)norm2(x)
    write(*,*)
    write(*,*)'BFGS_S_NH_fdwithf'
        call random_number(x)
        QuasiNewtonWarning=.true.
        call BFGS_Strong_NumericalHessian_fdwithf(f,fd,f_fd,fd_j,x,dim,freq)
        write(*,*)norm2(x)
    write(*,*)
    write(*,*)'LBFGS'
        call random_number(x)
        QuasiNewtonWarning=.true.
        call LBFGS(f,fd,x,dim,memory)
        write(*,*)norm2(x)
    write(*,*)
    write(*,*)'LBFGS_S'
        call random_number(x)
        QuasiNewtonWarning=.true.
        call LBFGS_Strong(f,fd,x,dim,memory)
        write(*,*)norm2(x)
    write(*,*)
    write(*,*)'LBFGS_S_fdwithf'
        call random_number(x)
        QuasiNewtonWarning=.true.
        call LBFGS_Strong_fdwithf(f,fd,f_fd,x,dim,memory)
        write(*,*)norm2(x)
    write(*,*)
    write(*,*)'DY'
        call random_number(x)
        ConjugateGradientWarning=.true.
        call DYConjugateGradient(f,fd,x,dim)
        write(*,*)norm2(x)
    write(*,*)
    write(*,*)'DY_S'
        call random_number(x)
        ConjugateGradientWarning=.true.
        call DYConjugateGradient_Strong(f,fd,x,dim)
        write(*,*)norm2(x)
    write(*,*)
    write(*,*)'DY_S_fdwithf'
        call random_number(x)
        ConjugateGradientWarning=.true.
        call DYConjugateGradient_Strong_fdwithf(f,fd,f_fd,x,dim)
        write(*,*)norm2(x)
    write(*,*)
    write(*,*)'PR'
        call random_number(x)
        ConjugateGradientWarning=.true.
        call PRConjugateGradient(f,fd,x,dim)
        write(*,*)norm2(x)
    write(*,*)
    write(*,*)'PR_fdwithf'
        call random_number(x)
        ConjugateGradientWarning=.true.
        call PRConjugateGradient_fdwithf(f,fd,f_fd,x,dim)
        write(*,*)norm2(x)
    write(*,*)
    write(*,*)'dtrnlsp'
        call random_number(x)
        trnlspWarning=.true.
        call my_dtrnlsp(fd_tr,fdd_tr,x,M,N)
        write(*,*)norm2(x)
    write(*,*)
    write(*,*)'dtrnlsp_NJ'
        call random_number(x)
        trnlspWarning=.true.
        call my_dtrnlsp_NumericalJacobian(fd_j,x,M,N)
        write(*,*)norm2(x)
    write(*,*)
    write(*,*)'dtrnlspbc'
        call random_number(x)
        trnlspWarning=.true.
        call my_dtrnlspbc(fd_tr,fdd_tr,x,low,up,M,N)
        write(*,*)norm2(x)
    write(*,*)
    write(*,*)'dtrnlspbc_NJ'
        call random_number(x)
        trnlspWarning=.true.
        call My_dtrnlspbc_NumericalJacobian(fd_j,x,low,up,M,N)
        write(*,*)norm2(x)
    write(*,*)
write(*,*)'Nonlinear-optimization solvers test passed'
write(*,*)

write(*,*)'Testing all geometry transformation routines...'
    write(*,*)
    write(*,*)'Standardize geometry'
        H2Ogeom(1:3)=[0d0,0d0,1d0]
        H2Ogeom(4:6)=[0d0,-1d0,0d0]
        H2Ogeom(7:9)=[0d0,1d0,0d0]
        mass=[18d0,1d0,1d0]
        H2Ogeom1=H2Ogeom
        call StandardizeGeometry(H2Ogeom1,mass,3)
        quaternion=BetterRandomUnitQuaternion()
        H2Ogeom2=H2Ogeom
        do i=1,3
            H2Ogeom2(3*i-2:3*i)=Rotate(quaternion,H2Ogeom2(3*i-2:3*i))
        end do
        call StandardizeGeometry2Reference(H2Ogeom2,H2Ogeom1,difference,mass,3)
        write(*,*)difference
    write(*,*)
write(*,*)'Geometry transformation routines test passed'
write(*,*)

write(*,*)'Testing all nonadiabatic routines...'
    call InitializePhaseFixing(10)
    write(*,*)
    write(*,*)'dFixdHPhase'
        call random_number(A3)
        B3=A3
        do i=1,10
            DoubleTemp=BetterRandomNumber_minus1to1()
            DoubleTemp=DoubleTemp/dAbs(DoubleTemp)
            B3(:,:,i)=DoubleTemp*B3(:,:,i)
            B3(:,i,:)=DoubleTemp*B3(:,i,:)
        end do
        call dFixdHPhase(B3,A3,difference,10,10)
        write(*,*)difference
    write(*,*)
    write(*,*)'dFixHPhaseBydH'
        call random_number(A)
        B=A
        call random_number(A3)
        B3=A3
        do i=1,10
            DoubleTemp=BetterRandomNumber_minus1to1()
            DoubleTemp=DoubleTemp/dAbs(DoubleTemp)
            B(:,i)=DoubleTemp*B(:,i)
            B(i,:)=DoubleTemp*B(i,:)
            B3(:,:,i)=DoubleTemp*B3(:,:,i)
            B3(:,i,:)=DoubleTemp*B3(:,i,:)
        end do
        call dFixHPhaseBydH(B,B3,A3,difference,10,10)
        write(*,*)My_dlansy('M',A-B,10)
    write(*,*)
    write(*,*)'dAssignBasisPhaseBydH'
        A=UnitMatrix(10)
        call random_number(A3)
        B3=A3
        do i=1,10
            DoubleTemp=BetterRandomNumber_minus1to1()
            DoubleTemp=DoubleTemp/dAbs(DoubleTemp)
            B(:,i)=DoubleTemp*A(:,i)
            B3(:,:,i)=DoubleTemp*B3(:,:,i)
            B3(:,i,:)=DoubleTemp*B3(:,i,:)
        end do
        call dAssignBasisPhaseBydH(B,B3,A3,difference,10,10)
        if(B(1,1)==A(1,1)) then
            write(*,*)'B = A',My_dlange('M',A-B,10,10)
        else
            write(*,*)'B =-A',My_dlange('M',A+B,10,10)
        end if
    write(*,*)
    write(*,*)'dFixHPhase_AssignBasisPhaseBydH'
        call random_number(A)
        B=A
        C=UnitMatrix(10)
        call random_number(A3)
        B3=A3
        do i=1,10
            DoubleTemp=BetterRandomNumber_minus1to1()
            DoubleTemp=DoubleTemp/dAbs(DoubleTemp)
            eigvec(:,i)=DoubleTemp*C(:,i)
            B(:,i)=DoubleTemp*B(:,i)
            B(i,:)=DoubleTemp*B(i,:)
            B3(:,:,i)=DoubleTemp*B3(:,:,i)
            B3(:,i,:)=DoubleTemp*B3(:,i,:)
        end do
        call dFixHPhase_AssignBasisPhaseBydH(B,eigvec,B3,A3,difference,10,10)
        if(eigvec(1,1)==C(1,1)) then
            write(*,*)'B = A',My_dlansy('M',A-B,10),My_dlange('M',C-eigvec,10,10)
        else
            write(*,*)'B =-A',My_dlansy('M',A-B,10),My_dlange('M',C+eigvec,10,10)
        end if
    write(*,*)
write(*,*)'Nonadiabatic routines test passed'
write(*,*)

contains
!Routines for testing nonlinear-optimization solvers
    subroutine f(fx,x,dim)
        real*8,intent(out)::fx
        integer,intent(in)::dim
        real*8,dimension(dim),intent(in)::x
        integer::i
        fx=0d0
        do i=1,dim
            fx=fx+x(i)**4
        end do
    end subroutine f
    
    subroutine fd(fdx,x,dim)
        integer,intent(in)::dim
        real*8,dimension(dim),intent(out)::fdx
        real*8,dimension(dim),intent(in)::x
        integer::i
        forall(i=1:dim)
            fdx(i)=4d0*x(i)**3
        end forall
    end subroutine fd

    subroutine f_fd(fx,fdx,x,dim)
        integer,intent(in)::dim
        real*8,intent(out)::fx
        real*8,dimension(dim),intent(out)::fdx
        real*8,dimension(dim),intent(in)::x
        integer::i
        fx=0d0
        do i=1,dim
            fx=fx+x(i)**4
            fdx(i)=4d0*x(i)**3
        end do
    end subroutine f_fd

    subroutine fdd(fddx,x,N)
        integer,intent(in)::N
        real*8,dimension(N,N),intent(out)::fddx
        real*8,dimension(N),intent(in)::x
        integer::i
        fddx=0d0
        forall(i=1:dim)
            fddx(i,i)=12d0*x(i)*x(i)
        end forall
    end subroutine fdd

    subroutine fd_tr(fdx,x,M,N)
        integer,intent(in)::M,N
        real*8,dimension(M),intent(out)::fdx
        real*8,dimension(N),intent(in)::x
        integer::i
        forall(i=1:dim)
            fdx(i)=4d0*x(i)**3
        end forall
    end subroutine fd_tr

    subroutine fdd_tr(fddx,x,M,N)
        integer,intent(in)::M,N
        real*8,dimension(M,N),intent(out)::fddx
        real*8,dimension(N),intent(in)::x
        integer::i
        fddx=0d0
        forall(i=1:dim)
            fddx(i,i)=12d0*x(i)*x(i)
        end forall
    end subroutine fdd_tr

    subroutine fd_j(M,N,x,fdx)
        integer,intent(in)::M,N
        real*8,dimension(M),intent(out)::fdx
        real*8,dimension(N),intent(in)::x
        integer::i
        forall(i=1:dim)
            fdx(i)=4d0*x(i)**3
        end forall
    end subroutine fd_j

end program main