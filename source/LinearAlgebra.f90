!Vector & matrix & high order tensor operation, BLAS & LAPACK routine wrapper
!
!Instruction:
!Naming convention (following LAPACK):
!    (overloaded) s = real*4, d = real*8, c = complex*8, z = complex*16
!    ge = general matrix, sy = real symmetric matrix, asy = anti symmetric matrix
!    he = Hermitian matrix, ahe = anti Hermitian matrix
!    po = real symmetric or Hermitian positive definite matrix
!Only use lower triangle of sy & he & po, strictly lower triangle of asy & ahe, otherwise specified
module LinearAlgebra
    implicit none

!Overload
    !--------------- Vector ----------------
        interface cross_product
            module procedure scross_product, dcross_product
        end interface cross_product

        interface triple_product
            module procedure striple_product, dtriple_product
        end interface triple_product
    !----------------- End -----------------

    !--------------- Matrix ----------------
        interface syL2U
            module procedure ssyL2U, dsyL2U
        end interface syL2U
    !----------------- End -----------------
    
    !------------ Linear solver ------------
        interface My_gesv
            module procedure My_dgesv, My_dgesvM
        end interface My_gesv

        interface My_sysv
            module procedure My_dsysv, My_dsysvM
        end interface My_sysv

        interface My_posv
            module procedure My_dposv, My_dposvM
        end interface My_posv

        interface My_getri
            module procedure My_dgetri, My_zgetri
        end interface My_getri

        interface My_potri
            module procedure My_spotri, My_dpotri
        end interface My_potri

        interface GeneralizedInverseTranspose
            module procedure sGeneralizedInverseTranspose, dGeneralizedInverseTranspose
        end interface GeneralizedInverseTranspose
    !----------------- End -----------------

    !------------- Eigensystem -------------
        interface My_geev
            module procedure My_dgeev, My_zgeev
        end interface My_geev
    !----------------- End -----------------

    !------------- Matrix norm -------------
        interface norm2ge
            module procedure dnorm2ge, znorm2ge
        end interface norm2ge

        interface My_lange
            module procedure My_dlange, My_zlange
        end interface My_lange
    !----------------- End -----------------

contains
!--------------- Vector ----------------
    !cross_product(a,b) = a x b
    !float a, b
    function scross_product(a, b)
        real*4,dimension(3),intent(in)::a,b
        real*4,dimension(3)::scross_product
        scross_product(1)=a(2)*b(3)-a(3)*b(2)
        scross_product(2)=a(3)*b(1)-a(1)*b(3)
        scross_product(3)=a(1)*b(2)-a(2)*b(1)
    end function scross_product
    !double a, b
    function dcross_product(a, b)
        real*8,dimension(3),intent(in)::a,b
        real*8,dimension(3)::dcross_product
        dcross_product(1)=a(2)*b(3)-a(3)*b(2)
        dcross_product(2)=a(3)*b(1)-a(1)*b(3)
        dcross_product(3)=a(1)*b(2)-a(2)*b(1)
    end function dcross_product
    
    !triple_product(a,b,c) = ( a x b ) . c
    !float a, b, c
    real*4 function striple_product(a, b, c)
        real*4,dimension(3),intent(in)::a,b,c
        striple_product=c(1)*(a(2)*b(3)-a(3)*b(2))+c(2)*(a(3)*b(1)-a(1)*b(3))+c(3)*(a(1)*b(2)-a(2)*b(1))
    end function striple_product
    !double a, b, c
    real*8 function dtriple_product(a, b, c)
        real*8,dimension(3),intent(in)::a,b,c
        dtriple_product=c(1)*(a(2)*b(3)-a(3)*b(2))+c(2)*(a(3)*b(1)-a(1)*b(3))+c(3)*(a(1)*b(2)-a(2)*b(1))
    end function dtriple_product
    
    !M dimensional vector a, N dimensional vector b, vector_direct_product(a,b) = a b
    function vector_direct_product(a, b, M, N)
        integer,intent(in)::M,N
        real*8,dimension(M),intent(in)::a
        real*8,dimension(N),intent(in)::b
        real*8,dimension(M,N)::vector_direct_product
        integer::i
        forall(i=1:N)
            vector_direct_product(:,i)=a*b(i)
        end forall
    end function vector_direct_product

    !N dimensional vector a, vector_direct_square(a) = a a, which is symmetric
    function vector_direct_square(a, N)
        integer,intent(in)::N
        real*8,dimension(N),intent(in)::a
        real*8,dimension(N,N)::vector_direct_square
        integer::i,j
        forall(i=1:N,j=1:N,i>=j)
            vector_direct_square(i,j)=a(i)*a(j)
        end forall
    end function vector_direct_square

    !M dimensional vector a, N dimensional vector b, return [ a^T, b^T ]^T
    function vector_direct_sum(a, b, M, N)
        integer,intent(in)::M,N
        real*8,dimension(M),intent(in)::a
        real*8,dimension(N),intent(in)::b
        real*8,dimension(M+N)::vector_direct_sum
        vector_direct_sum(1:M)=a
        vector_direct_sum(M+1:M+N)=b
    end function vector_direct_sum
!----------------- End -----------------

!--------------- Matrix ----------------
    !N order matrix A, return det(A)
    real*8 function determinant(A, N)
        integer,intent(in)::N
        real*8,dimension(N,N),intent(in)::A
        integer::i; integer,dimension(N)::ipiv; real*8::sign
        real*8,dimension(N,N)::Acopy
        Acopy=A; call dgetrf(N,N,Acopy,N,ipiv,i)
        if(ipiv(1)==1) then; sign=1d0; else; sign=-1d0; end if
        determinant=Acopy(1,1)
        do i=2,N
            if(ipiv(i)/=i) sign=-sign
            determinant=determinant*Acopy(i,i)
        end do
        determinant=determinant*sign
    end function determinant

    !N order matrix A, return Tr(A)
    real*8 function Trace(A, N)
        integer,intent(in)::N
        real*8,dimension(N,N),intent(in)::A
        integer::i
        Trace=0d0
        do i=1,N; Trace=Trace+A(i,i); end do
    end function Trace

    !N order matrix A, return the main diagonal vector of A
    function DiagVector(A, N)
        integer,intent(in)::N
        real*8,dimension(N,N),intent(in)::A
        real*8,dimension(N)::DiagVector
        integer::i
        forall(i=1:N); DiagVector(i)=A(i,i); end forall
    end function DiagVector

    !M x K matrix A, K x N matrix B, return A . B
    function matmul_dgemm(A, B, M, K, N)
        integer,intent(in)::M,K,N
        real*8,dimension(M,K),intent(in)::A
        real*8,dimension(K,N),intent(in)::B
        real*8,dimension(M,N)::matmul_dgemm
        call dgemm('N','N',M,N,K,1d0,A,M,B,K,0d0,matmul_dgemm,M)
    end function

    !M order real symmetric matrix A, M x N matrix B, return A . B
    function matmul_dsymm(A, B, M, N)
        integer,intent(in)::M,N
        real*8,dimension(M,M),intent(in)::A
        real*8,dimension(M,N),intent(in)::B
        real*8,dimension(M,N)::matmul_dsymm
        call dsymm('L','L',M,N,1d0,A,M,B,M,0d0,matmul_dsymm,M)
    end function

    !M x N matrix A, N dimensional vector x, return A . x
    function mvmul_dgemv(A, x, M, N)
        integer,intent(in)::M,N
        real*8,dimension(M,N),intent(in)::A
        real*8,dimension(N),intent(in)::x
        real*8,dimension(N)::mvmul_dgemv
        call dgemv('N',M,N,1d0,A,M,x,1,0d0,mvmul_dgemv,1)
    end function mvmul_dgemv

    !MA x NA matrix A, MB x NB matrix B, ( A direct_product B )_ijkl = A_ij * B_kl
    function matrix_direct_product(A, B, MA, NA, MB, NB)
        integer,intent(in)::MA,NA,MB,NB
        real*8,dimension(MA,NA),intent(in)::A
        real*8,dimension(MB,NB),intent(in)::B
        real*8,dimension(MA,NA,MB,NB)::matrix_direct_product
        integer::i,j
        forall(i=1:MB,j=1:NB)
            matrix_direct_product(:,:,i,j)=A*B(i,j)
        end forall
    end function matrix_direct_product

    !M order matrix A, N order matrix B, A direct_sum B = M + N order block diagonal matrix with A and B as diagonal blocks
    function matrix_direct_sum(A, B, M, N)
        integer,intent(in)::M,N
        real*8,dimension(M,M),intent(in)::A
        real*8,dimension(N,N),intent(in)::B
        real*8,dimension(M+N,M+N)::matrix_direct_sum
        matrix_direct_sum(1:M,1:M)=A
        matrix_direct_sum(M+1:M+N,1:M)=0d0
        matrix_direct_sum(1:M,M+1:M+N)=0d0
        matrix_direct_sum(M+1:M+N,M+1:M+N)=B
    end function matrix_direct_sum

    !N order blank matrix A, N order real symmetric matrix B, perform A = B
    subroutine sycp(A, B, N)
        integer,intent(in)::N
        real*8,dimension(N,N),intent(out)::A
        real*8,dimension(N,N),intent(in)::B
        integer::i,j
        forall(i=1:N,j=1:N,i>=j)
            A(i,j)=B(i,j)
        end forall
    end subroutine sycp

    !N order matrix A, strictly upper triangle is blank, copy strictly lower triangle elements to strictly upper triangle
    !float A
    subroutine ssyL2U(A, N)
        integer,intent(in)::N
        real*4,dimension(N,N),intent(inout)::A
        integer::i,j
        forall(i=1:N-1,j=2:N,i<j); A(i,j)=A(j,i); end forall
    end subroutine ssyL2U
    !double A
    subroutine dsyL2U(A, N)
        integer,intent(in)::N
        real*8,dimension(N,N),intent(inout)::A
        integer::i,j
        forall(i=1:N-1,j=2:N,i<j); A(i,j)=A(j,i); end forall
    end subroutine dsyL2U

    !N order real symmetric matrix A, N order anti symmetric matrix B, return A . B
    function symatmulasy(A, B, N)
        integer,intent(in)::N
        real*8,dimension(N,N),intent(in)::A,B
        real*8,dimension(N,N)::symatmulasy
        integer::i,j,k
        symatmulasy=0d0
        do j=1,N
            do i=1,j
                do k=1,i
                    symatmulasy(i,j)=symatmulasy(i,j)-A(i,k)*B(j,k)
                end do
                do k=i+1,j-1
                    symatmulasy(i,j)=symatmulasy(i,j)-A(k,i)*B(j,k)
                end do
                do k=j+1,N
                    symatmulasy(i,j)=symatmulasy(i,j)+A(k,i)*B(k,j)
                end do
            end do
            do i=j+1,N
                do k=1,j-1
                    symatmulasy(i,j)=symatmulasy(i,j)-A(i,k)*B(j,k)
                end do
                do k=j+1,i
                    symatmulasy(i,j)=symatmulasy(i,j)+A(i,k)*B(k,j)
                end do
                do k=i+1,N
                    symatmulasy(i,j)=symatmulasy(i,j)+A(k,i)*B(k,j)
                end do
            end do
        end do
    end function symatmulasy

    !N order anti symmetric matrix A, N order real symmetric matrix B, return A . B
    function asymatmulsy(A, B, N)
        integer,intent(in)::N
        real*8,dimension(N,N),intent(in)::A,B
        real*8,dimension(N,N)::asymatmulsy
        integer::i,j,k
        asymatmulsy=0d0
        do j=1,N
            do i=1,j
                do k=1,i-1
                    asymatmulsy(i,j)=asymatmulsy(i,j)+A(i,k)*B(j,k)
                end do
                do k=i+1,j
                    asymatmulsy(i,j)=asymatmulsy(i,j)-A(k,i)*B(j,k)
                end do
                do k=j+1,N
                    asymatmulsy(i,j)=asymatmulsy(i,j)-A(k,i)*B(k,j)
                end do
            end do
            do i=j+1,N
                do k=1,j
                    asymatmulsy(i,j)=asymatmulsy(i,j)+A(i,k)*B(j,k)
                end do
                do k=j+1,i-1
                    asymatmulsy(i,j)=asymatmulsy(i,j)+A(i,k)*B(k,j)
                end do
                do k=i+1,N
                    asymatmulsy(i,j)=asymatmulsy(i,j)-A(k,i)*B(k,j)
                end do
            end do
        end do
    end function asymatmulsy
!----------------- End -----------------

!---------- High order tensor ----------
    !dim x N x N 3rd-order tensor A, Trace3(i) = Tr[A(i,:,:)]
    function Trace3(A, dim, N)
        integer,intent(in)::dim,N
        real*8,dimension(dim,N,N),intent(in)::A
        real*8,dimension(dim)::Trace3
        integer::i
        Trace3=0d0
        do i=1,N
            Trace3=Trace3+A(:,i,i)
        end do
    end function Trace3

    !dim1 x dim2 x N x N 4th-order tensor A, Trace4(i,j) = Tr[A(i,j,:,:)]
    function Trace4(A, dim1, dim2, N)
        integer,intent(in)::dim1,dim2,N
        real*8,dimension(dim1,dim2,N,N),intent(in)::A
        real*8,dimension(dim1,dim2)::Trace4
        integer::i
        Trace4=0d0
        do i=1,N
            Trace4=Trace4+A(:,:,i,i)
        end do
    end function Trace4

    !dim x M x N 3rd-order tensor A, transpose on M x N
    function Transpose3(A, dim, M, N)
        integer,intent(in)::dim,M,N
        real*8,dimension(dim,M,N),intent(in)::A
        real*8,dimension(dim,N,M)::Transpose3
        integer::i
        forall(i=1:dim)
            Transpose3(i,:,:)=Transpose(A(i,:,:))
        end forall
    end function Transpose3

    !dim1 x dim2 x M x N 4rd-order tensor A, transpose on M x N
    function Transpose4(A, dim1, dim2, M, N)
        integer,intent(in)::dim1,dim2,M,N
        real*8,dimension(dim1,dim2,M,N),intent(in)::A
        real*8,dimension(dim1,dim2,N,M)::Transpose4
        integer::i,j
        forall(i=1:dim1,j=1:dim2)
            Transpose4(i,j,:,:)=Transpose(A(i,j,:,:))
        end forall
    end function Transpose4

    !dim x N x N 3rd-order tensor A, strictly upper triangle of N x N is blank
    !Copy strictly lower triangle elements to strictly upper triangle
    subroutine sy3L2U(A, dim, N)
        integer,intent(in)::dim,N
        real*8,dimension(dim,N,N),intent(inout)::A
        integer::i,j
        forall(i=1:N-1,j=2:N,i<j)
            A(:,i,j)=A(:,j,i)
        end forall
    end subroutine sy3L2U

    !dim1 x dim2 x N x N 3rd-order tensor A, strictly upper triangle of N x N is blank
    !Copy strictly lower triangle elements to strictly upper triangle
    subroutine sy4L2U(A, dim1, dim2, N)
        integer,intent(in)::dim1,dim2,N
        real*8,dimension(dim1,dim2,N,N),intent(inout)::A
        integer::i,j
        forall(i=1:N-1,j=2:N,i<j)
            A(:,:,i,j)=A(:,:,j,i)
        end forall
    end subroutine sy4L2U

    !dim x N x N 3rd-order tensor A, real symmetric in N x N. N order real symmetric matrix B
    !Matrix multiplication on N x N
    function sy3matmulsy(A, B, dim, N)
        integer,intent(in)::dim,N
        real*8,dimension(dim,N,N),intent(in)::A
        real*8,dimension(N,N),intent(in)::B
        real*8,dimension(dim,N,N)::sy3matmulsy
        integer::i,j,k,l
        sy3matmulsy=0d0
        do j=1,N
            do i=1,j
                do k=1,i
                    sy3matmulsy(:,i,j)=sy3matmulsy(:,i,j)+A(:,i,k)*B(j,k)
                end do
                do k=i+1,j
                    sy3matmulsy(:,i,j)=sy3matmulsy(:,i,j)+A(:,k,i)*B(j,k)
                end do
                do k=j+1,N
                    sy3matmulsy(:,i,j)=sy3matmulsy(:,i,j)+A(:,k,i)*B(k,j)
                end do
            end do
            do i=j+1,N
                do k=1,j
                    sy3matmulsy(:,i,j)=sy3matmulsy(:,i,j)+A(:,i,k)*B(j,k)
                end do
                do k=j+1,i
                    sy3matmulsy(:,i,j)=sy3matmulsy(:,i,j)+A(:,i,k)*B(k,j)
                end do
                do k=i+1,N
                    sy3matmulsy(:,i,j)=sy3matmulsy(:,i,j)+A(:,k,i)*B(k,j)
                end do
            end do
        end do
    end function sy3matmulsy

    !dim x N x N 3rd-order tensor A & B, real symmetric in N x N
    !Vector dot product on dim, matrix multiplication on N x N
    function sy3matdotmul(A, B, dim, N)
        integer,intent(in)::dim,N
        real*8,dimension(dim,N,N),intent(in)::A,B
        real*8,dimension(N,N)::sy3matdotmul
        integer::i,j,k
        sy3matdotmul=0d0
        do j=1,N
            do i=1,j
                do k=1,i
                    sy3matdotmul(i,j)=sy3matdotmul(i,j)+dot_product(A(:,i,k),B(:,j,k))
                end do
                do k=i+1,j
                    sy3matdotmul(i,j)=sy3matdotmul(i,j)+dot_product(A(:,k,i),B(:,j,k))
                end do
                do k=j+1,N
                    sy3matdotmul(i,j)=sy3matdotmul(i,j)+dot_product(A(:,k,i),B(:,k,j))
                end do
            end do
            do i=j+1,N
                do k=1,j
                    sy3matdotmul(i,j)=sy3matdotmul(i,j)+dot_product(A(:,i,k),B(:,j,k))
                end do
                do k=j+1,i
                    sy3matdotmul(i,j)=sy3matdotmul(i,j)+dot_product(A(:,i,k),B(:,k,j))
                end do
                do k=i+1,N
                    sy3matdotmul(i,j)=sy3matdotmul(i,j)+dot_product(A(:,k,i),B(:,k,j))
                end do
            end do
        end do
    end function sy3matdotmul

    !dim1 x dim2 x N x N 4th-order tensor A, dim2 x N x N 3rd-order tensor B, real symmetric in N x N
    !Vector dot product on dim2, matrix multiplication on N x N
    function sy4matdotmulsy3(A, B, dim1, dim2, N)
        integer,intent(in)::dim1,dim2,N
        real*8,dimension(dim1,dim2,N,N),intent(in)::A
        real*8,dimension(dim2,N,N),intent(in)::B
        real*8,dimension(dim1,N,N)::sy4matdotmulsy3
        integer::i,j,k,l
        sy4matdotmulsy3=0d0
        do j=1,N
            do i=1,j
                do k=1,i
                    forall(l=1:dim1)
                        sy4matdotmulsy3(l,i,j)=sy4matdotmulsy3(l,i,j)+dot_product(A(l,:,i,k),B(:,j,k))
                    end forall
                end do
                do k=i+1,j
                    forall(l=1:dim1)
                        sy4matdotmulsy3(l,i,j)=sy4matdotmulsy3(l,i,j)+dot_product(A(l,:,k,i),B(:,j,k))
                    end forall
                end do
                do k=j+1,N
                    forall(l=1:dim1)
                        sy4matdotmulsy3(l,i,j)=sy4matdotmulsy3(l,i,j)+dot_product(A(l,:,k,i),B(:,k,j))
                    end forall
                end do
            end do
            do i=j+1,N
                do k=1,j
                    forall(l=1:dim1)
                        sy4matdotmulsy3(l,i,j)=sy4matdotmulsy3(l,i,j)+dot_product(A(l,:,i,k),B(:,j,k))
                    end forall
                end do
                do k=j+1,i
                    forall(l=1:dim1)
                        sy4matdotmulsy3(l,i,j)=sy4matdotmulsy3(l,i,j)+dot_product(A(l,:,i,k),B(:,k,j))
                    end forall
                end do
                do k=i+1,N
                    forall(l=1:dim1)
                        sy4matdotmulsy3(l,i,j)=sy4matdotmulsy3(l,i,j)+dot_product(A(l,:,k,i),B(:,k,j))
                    end forall
                end do
            end do
        end do
    end function sy4matdotmulsy3

    !dim x N x N 3rd-order tensor A, anti symmetric in N x N. N order real symmetric matrix B
    !Matrix multiplication on N x N
    function asy3matmulsy(A, B, dim, N)
        integer,intent(in)::dim,N
        real*8,dimension(dim,N,N),intent(in)::A
        real*8,dimension(N,N),intent(in)::B
        real*8,dimension(dim,N,N)::asy3matmulsy
        integer::i,j,k,l
        asy3matmulsy=0d0
        do j=1,N
            do i=1,j
                do k=1,i-1
                    asy3matmulsy(:,i,j)=asy3matmulsy(:,i,j)+A(:,i,k)*B(j,k)
                end do
                do k=i+1,j
                    asy3matmulsy(:,i,j)=asy3matmulsy(:,i,j)-A(:,k,i)*B(j,k)
                end do
                do k=j+1,N
                    asy3matmulsy(:,i,j)=asy3matmulsy(:,i,j)-A(:,k,i)*B(k,j)
                end do
            end do
            do i=j+1,N
                do k=1,j
                    asy3matmulsy(:,i,j)=asy3matmulsy(:,i,j)+A(:,i,k)*B(j,k)
                end do
                do k=j+1,i-1
                    asy3matmulsy(:,i,j)=asy3matmulsy(:,i,j)+A(:,i,k)*B(k,j)
                end do
                do k=i+1,N
                    asy3matmulsy(:,i,j)=asy3matmulsy(:,i,j)-A(:,k,i)*B(k,j)
                end do
            end do
        end do
    end function asy3matmulsy

    !dim1 x N x N 3rd-order tensor A, anti symmetric in N x N
    !dim2 x N x N 3rd-order tensor B, real symmetric in N x N
    !Vector direct product on dim1 and dim2, matrix multiplication on N x N
    function asy3matdirectmulsy3(A, B, dim1, dim2, N)
        integer,intent(in)::dim1,dim2,N
        real*8,dimension(dim1,N,N),intent(in)::A
        real*8,dimension(dim2,N,N),intent(in)::B
        real*8,dimension(dim1,dim2,N,N)::asy3matdirectmulsy3
        integer::i,j,k
        asy3matdirectmulsy3=0d0
        do j=1,N
            do i=1,j
                do k=1,i-1
                    asy3matdirectmulsy3(:,:,i,j)=asy3matdirectmulsy3(:,:,i,j)+vector_direct_product(A(:,i,k),B(:,j,k),dim1,dim2)
                end do
                do k=i+1,j
                    asy3matdirectmulsy3(:,:,i,j)=asy3matdirectmulsy3(:,:,i,j)-vector_direct_product(A(:,k,i),B(:,j,k),dim1,dim2)
                end do
                do k=j+1,N
                    asy3matdirectmulsy3(:,:,i,j)=asy3matdirectmulsy3(:,:,i,j)-vector_direct_product(A(:,k,i),B(:,k,j),dim1,dim2)
                end do
            end do
            do i=j+1,N
                do k=1,j
                    asy3matdirectmulsy3(:,:,i,j)=asy3matdirectmulsy3(:,:,i,j)+vector_direct_product(A(:,i,k),B(:,j,k),dim1,dim2)
                end do
                do k=j+1,i-1
                    asy3matdirectmulsy3(:,:,i,j)=asy3matdirectmulsy3(:,:,i,j)+vector_direct_product(A(:,i,k),B(:,k,j),dim1,dim2)
                end do
                do k=i+1,N
                    asy3matdirectmulsy3(:,:,i,j)=asy3matdirectmulsy3(:,:,i,j)-vector_direct_product(A(:,k,i),B(:,k,j),dim1,dim2)
                end do
            end do
        end do
    end function asy3matdirectmulsy3

    !dim x N x N 3rd-order tensor A, real symmetric in N x N. N order unitary matrix U
    !Return U^T . A . U on N x N of A
    function sy3UnitaryTransformation(A, U, dim, N)
        integer,intent(in)::dim,N
        real*8,dimension(dim,N,N),intent(in)::A
        real*8,dimension(N,N),intent(in)::U
        real*8,dimension(dim,N,N)::sy3UnitaryTransformation
        integer::i,j,k,l
        forall(i=1:N,j=1:N,i>=j)
            sy3UnitaryTransformation(:,i,j)=0d0
        end forall
        do j=1,N
            do i=j,N
                do k=1,N
                    sy3UnitaryTransformation(:,i,j)=sy3UnitaryTransformation(:,i,j)+U(k,i)*U(k,j)*A(:,k,k)
                    do l=k+1,N
                        sy3UnitaryTransformation(:,i,j)=sy3UnitaryTransformation(:,i,j)+(U(k,i)*U(l,j)+U(l,i)*U(k,j))*A(:,l,k)
                    end do
                end do
            end do
        end do
    end function sy3UnitaryTransformation

    !dim1 x dim2 x N x N 4th-order tensor A, real symmetric in N x N. N order unitary matrix U
    !Return U^T . A . U on N x N of A
    function sy4UnitaryTransformation(A, U, dim1, dim2, N)
        integer,intent(in)::dim1,dim2,N
        real*8,dimension(dim1,dim2,N,N),intent(in)::A
        real*8,dimension(N,N),intent(in)::U
        real*8,dimension(dim1,dim2,N,N)::sy4UnitaryTransformation
        integer::i,j,k,l
        forall(i=1:N,j=1:N,i>=j)
            sy4UnitaryTransformation(:,:,i,j)=0d0
        end forall
        do j=1,N
            do i=j,N
                do k=1,N
                    sy4UnitaryTransformation(:,:,i,j)=sy4UnitaryTransformation(:,:,i,j)+U(k,i)*U(k,j)*A(:,:,k,k)
                    do l=k+1,N
                        sy4UnitaryTransformation(:,:,i,j)=sy4UnitaryTransformation(:,:,i,j)+(U(k,i)*U(l,j)+U(l,i)*U(k,j))*A(:,:,l,k)
                    end do
                end do
            end do
        end do
    end function sy4UnitaryTransformation
!----------------- End -----------------

!------------ Linear solver ------------
    !N order matrix A
    !ge method: Doolittle LU decomposition
    !sy method: Bunch-Kaufman L . D . L^T decomposition (D is 1x1 and 2x2 block diagonal)
    !po method: Cholesky L . L^T decomposition

    !=========== Solve ===========
        !Solve the linear system A . x = b
        !b harvests the solution x
        
        !A harvests the Doolittle LU decomposition
        !N dimensional vector b
        subroutine My_dgesv(A, b, N)
            integer,intent(in)::N
            real*8,dimension(N,N),intent(inout)::A
            real*8,dimension(N),intent(inout)::b
            integer::info
            integer,dimension(N)::ipiv
            call dgesv(N,1,A,N,ipiv,b,N,info)
        end subroutine My_dgesv
        !N x M matrix b
        subroutine My_dgesvM(A, b, N, M)
            integer,intent(in)::N,M
            real*8,dimension(N,N),intent(inout)::A
            real*8,dimension(N,M),intent(inout)::b
            integer::info
            integer,dimension(N)::ipiv
            call dgesv(N,M,A,N,ipiv,b,N,info)
        end subroutine My_dgesvM
    
        !A harvests the Bunch-Kaufman L . D . L^T decomposition (D is 1x1 and 2x2 block diagonal)
        !N dimensional vector b
        subroutine My_dsysv(A, b, N)
            integer,intent(in)::N
            real*8,dimension(N,N),intent(inout)::A
            real*8,dimension(N),intent(inout)::b
            integer::info
            integer,dimension(N)::ipiv
            real*8,dimension(N)::work
            call dsysv('L',N,1,A,N,ipiv,b,N,work,N,info)
        end subroutine My_dsysv
        !N x M matrix b
        subroutine My_dsysvM(A, b, N, M)
            integer,intent(in)::N,M
            real*8,dimension(N,N),intent(inout)::A
            real*8,dimension(N,M),intent(inout)::b
            integer::info
            integer,dimension(N)::ipiv
            real*8,dimension(N)::work
            call dsysv('L',N,M,A,N,ipiv,b,N,work,N,info)
        end subroutine My_dsysvM
    
        !Optional: info: If A is po, info returns 0; else, the info-th leading minor of A <= 0 and solving failed
        !A harvests the Cholesky L . L^T decomposition. A will be overwritten even fail
        !b will not be overwritten if fail
        !N dimensional vector b
        subroutine My_dposv(A, b, N, info)
            integer,intent(in)::N
            real*8,dimension(N,N),intent(inout)::A
            real*8,dimension(N),intent(inout)::b
            integer,intent(out),optional::info
            integer::temp
            if(present(info)) then
                call dposv('L',N,1,A,N,b,N,info)
            else
                call dposv('L',N,1,A,N,b,N,temp)
            end if
        end subroutine My_dposv
        !N x M matrix b
        subroutine My_dposvM(A, b, N, M, info)
            integer,intent(in)::N,M
            real*8,dimension(N,N),intent(inout)::A
            real*8,dimension(N,M),intent(inout)::b
            integer,intent(out),optional::info
            integer::temp
            if(present(info)) then
                call dposv('L',N,M,A,N,b,N,info)
            else
                call dposv('L',N,M,A,N,b,N,temp)
            end if
        end subroutine My_dposvM
    !============ End ============

    !========== Inverse ==========
        !Inverse A
        !A harvests its inverse

        subroutine My_dgetri(A, N)
            integer,intent(in)::N
            real*8,dimension(N,N),intent(inout)::A
            integer::info
            integer,dimension(N)::ipiv
            real*8,dimension(N)::work
            call dgetrf(N,N,A,N,ipiv,info)
            call dgetri(N,A,N,ipiv,work,N,info)
        end subroutine My_dgetri
        subroutine My_zgetri(A, N)
            integer,intent(in)::N
            complex*16,dimension(N,N),intent(inout)::A
            integer::info
            integer,dimension(N)::ipiv
            complex*16,dimension(N)::work
            call dgetrf(N,N,A,N,ipiv,info)
            call dgetri(N,A,N,ipiv,work,N,info)
        end subroutine My_zgetri

        subroutine My_dsytri(A, N)
            integer,intent(in)::N
            real*8,dimension(N,N),intent(inout)::A
            integer::info
            integer,dimension(N)::ipiv
            real*8,dimension(N)::work
            call dsytrf('L',N,A,N,ipiv,work,N,info)
            call dsytri('L',N,A,N,ipiv,work,info)
        end subroutine My_dsytri

        !Optional: info: if A is po, info returns 0; else, the info-th leading minor of A <= 0 and inversing failed
        !A will be overwritten even fail
        !float A
        subroutine My_spotri(A, N, info)
            integer,intent(in)::N
            real*4,dimension(N,N),intent(inout)::A
            integer,intent(out),optional::info
            integer::temp
            if(present(info)) then
                call spotrf('L',N,A,N,info)
                if(info/=0) return
                call spotri('L',N,A,N,info)
            else
                call spotrf('L',N,A,N,temp)
                if(temp/=0) return
                call spotri('L',N,A,N,temp)
            end if
        end subroutine My_spotri
        !double A
        subroutine My_dpotri(A, N, info)
            integer,intent(in)::N
            real*8,dimension(N,N),intent(inout)::A
            integer,intent(out),optional::info
            integer::temp
            if(present(info)) then
                call dpotrf('L',N,A,N,info)
                if(info/=0) return
                call dpotri('L',N,A,N,info)
            else
                call dpotrf('L',N,A,N,temp)
                if(temp/=0) return
                call dpotri('L',N,A,N,temp)
            end if
        end subroutine My_dpotri

        !M x N real matrix A. A harvests the transpose of its generalized inverse
        !float A
        subroutine sGeneralizedInverseTranspose(A, M, N)
            integer,intent(in)::M,N
            real*4,dimension(M,N),intent(inout)::A
            real*4,dimension(M,M)::AAT
            AAT=matmul(A,transpose(A))
            call My_potri(AAT,M)
            call syL2U(AAT,M)
            A=matmul(AAT,A)
        end subroutine sGeneralizedInverseTranspose
        !double A
        subroutine dGeneralizedInverseTranspose(A, M, N)
            integer,intent(in)::M,N
            real*8,dimension(M,N),intent(inout)::A
            real*8,dimension(M,M)::AAT
            AAT=matmul(A,transpose(A))
            call My_potri(AAT,M)
            call syL2U(AAT,M)
            A=matmul(AAT,A)
        end subroutine dGeneralizedInverseTranspose
    !============ End ============
!----------------- End -----------------

!------------- Eigensystem -------------
    !N order matrix A
    !Compute the eigen system of A
    !jobtype: 'N' eigenvalues only, 'V' eigenvectors as well
    !ge method: QR decomposition
    !sy/he method: unitary transform to real symmetric tridiagonal, then
    !    if 'N', Pal-Walker-Kahan QR; else 'V', implicit QR

    !A will be overwritten
    !eigvec harvests the normalized eigenvectors
    !eigvalr harvests the real part of the eigenvalues,
    !eigvali harvests the imaginary part of the eigenvalues,
    !Only right eigen: A . eigvec(:,j) = [ eigvalr(j) + i * eigvali(j) ] * eigvec(:,j)
    subroutine My_dgeev(jobtype, A,  eigvalr,eigvali, eigvec, N)
        character,intent(in)::jobtype
        integer,intent(in)::N
        real*8,dimension(N,N),intent(in)::A
        real*8,dimension(N),intent(out)::eigvalr,eigvali
        real*8,dimension(N,N),intent(out)::eigvec
        integer::info
        real*8,dimension(5*N)::work
        real*8,dimension(1,N)::vl
        call dgeev('N',jobtype,N,A,N,eigvalr,eigvali,vl,1,eigvec,N,work,5*N,info)
    end subroutine My_dgeev
    !eigval harvests the eigenvalues
    !Only right eigen: A . eigvec(:,j) = eigval(j) * eigvec(:,j)
    subroutine My_zgeev(jobtype, A, eigval, eigvec, N)
        character,intent(in)::jobtype
        integer,intent(in)::N
        complex*16,dimension(N,N),intent(in)::A
        complex*16,dimension(N),intent(inout)::eigval
        complex*16,dimension(N,N),intent(inout)::eigvec
        integer::lwork,info
        real*8,dimension(2*n)::rwork
        complex*16,dimension(3*N)::work
        complex*16,dimension(1,N)::vl
        call zgeev('N',jobtype,N,A,N,eigval,vl,1,eigvec,N,work,3*N,rwork,info)
    end subroutine My_zgeev
    
    !eigval harvests the eigenvalues in ascending order, A harvests the normalized eigenvectors
    !A will be overwritten even for 'N' job
    subroutine My_dsyev(jobtype, A, eigval, N)
        integer,intent(in)::N
        character,intent(in)::jobtype
        real*8,dimension(N,N),intent(inout)::A
        real*8,dimension(N),intent(out)::eigval
        integer::info
        real*8,dimension(3*N)::work
        call dsyev(jobtype,'L',N,A,N,eigval,work,3*N,info)
    end subroutine My_dsyev

    !gtype: 1, A . eigvec = S . eigvec . diag(eigval)
    !       2, A . S . eigvec = eigvec . diag(eigval)
    !N order real symmetric positive definite matrix S
    !Optional: info: return 0 if normal termination; else, the info-th leading minor of S <= 0 and diagonalization failed
    !eigval harvests the eigenvalues in ascending order, A harvests the eigenvectors normalized under S metric
    !S harvests the Cholesky L . L^T decomposition. S will be overwritten even fail
    !A will be overwritten even for 'N' job
    subroutine My_dsygv(gtype, jobtype, A, S, eigval, N, info)
        integer,intent(in)::gtype,N
        character,intent(in)::jobtype
        real*8,dimension(N,N),intent(inout)::A,S
        real*8,dimension(N),intent(out)::eigval
        integer,intent(out),optional::info
        integer::temp
        real*8,dimension(3*N)::work
        if(present(info)) then
            call dsygv(gtype,jobtype,'L',N,A,N,S,N,eigval,work,3*N,info)
            if(info>0) info=info-N
        else
            call dsygv(gtype,jobtype,'L',N,A,N,S,N,eigval,work,3*N,temp)
        end if
    end subroutine My_dsygv

    !eigval harvests the eigenvalues in ascending order, A harvests the normalized eigenvectors
    !A will be overwritten even for 'N' job
    subroutine My_zheev(jobtype, A, eigval, N)
        integer,intent(in)::N
        character,intent(in)::jobtype
        complex*16,dimension(N,N),intent(inout)::A
        real*8,dimension(N),intent(out)::eigval
        integer::info
        real*8,dimension(3*N-2)::rwork
        complex*16,dimension(2*N)::work
        call zheev(jobtype,'L',N,A,N,eigval,work,2*N,rwork,info)
    end subroutine My_zheev
!----------------- End -----------------

!------------- Matrix norm -------------
    !============= 2-norm ==============
        !M x N matrix A, return 2 norm of A
    
        !Double A
        real*8 function dnorm2ge(A, M, N)
            integer,intent(in)::M,N
            real*8,dimension(M,N),intent(in)::A
            real*8,dimension(N)::eigval
            real*8,dimension(N,N)::ATA
            ATA=matmul(transpose(A),A)
            call My_dsyev('N',ATA,eigval,N)
            dnorm2ge=Sqrt(maxval(eigval))
        end function dnorm2ge
        !Complex A
        real*8 function znorm2ge(A, M, N)
            integer,intent(in)::M,N
            complex*16,dimension(M,N),intent(in)::A
            real*8,dimension(N)::eigval
            complex*16,dimension(N,N)::ATA
            ATA=matmul(conjg(transpose(A)),A)
            call My_zheev('N',ATA,eigval,N)
            znorm2ge=Sqrt(maxval(eigval))
        end function znorm2ge
    !=============== End ===============

    !========== F-norm square ==========
        !Return Frobenius norm square of M x N matrix A

        real*8 function dgeFrobeniusSquare(A, M, N)
            integer,intent(in)::M,N
            real*8,dimension(M,N),intent(in)::A
            integer::i
            dgeFrobeniusSquare=dot_product(A(:,1),A(:,1))
            do i=2,N
                dgeFrobeniusSquare=dgeFrobeniusSquare+dot_product(A(:,i),A(:,i))
            end do
        end function dgeFrobeniusSquare

        real*8 function dsyFrobeniusSquare(A, N)
            integer,intent(in)::N
            real*8,dimension(N,N),intent(in)::A
            integer::i
            real*8::temp
            dsyFrobeniusSquare=A(1,1)*A(1,1)
            temp=dot_product(A(2:N,1),A(2:N,1))
            do i=2,N-1
                dsyFrobeniusSquare=dsyFrobeniusSquare+A(i,i)*A(i,i)
                temp=temp+dot_product(A(i+1:N,i),A(i+1:N,i))
            end do
            dsyFrobeniusSquare=dsyFrobeniusSquare+A(N,N)*A(N,N)+2d0*temp
        end function dsyFrobeniusSquare
    !=============== End ===============
    
    !=========== Other norms ===========
        !M x N matrix A
        !jobtype: 'M', return max(abs(A(i,j))) (note this is not a subordinate norm)
        !         'F', return Frobenius norm = Sqrt(sum of element squares) = Sqrt[Tr(A^T.A)]
        !              also called Euclidean norm, note this is not a subordinate norm
        !         '1', return 1 norm
        !         'I', return infinity norm

        !Double A
        real*8 function My_dlange(jobtype, A, M, N)
            character,intent(in)::jobtype
            integer,intent(in)::M,N
            real*8,dimension(M,N),intent(in)::A
            real*8,external::dlange
            real*8,dimension(M)::work
            My_dlange=dlange(jobtype,M,N,A,M,work)
        end function My_dlange
        !Complex A
        real*8 function My_zlange(jobtype, A, M, N)
            character,intent(in)::jobtype
            integer,intent(in)::M,N
            complex*16,dimension(M,N)::A
            real*8,external::zlange
            real*8,dimension(M)::work
            My_zlange=zlange(jobtype,M,N,A,M,work)
        end function My_zlange

        real*8 function My_dlansy(jobtype, A, N)
            character,intent(in)::jobtype
            integer,intent(in)::N
            real*8,dimension(N,N),intent(in)::A
            real*8,external::dlansy
            real*8,dimension(N)::work
            My_dlansy=dlansy(jobtype,'L',N,A,N,work)
        end function My_dlansy
    !=============== End ===============
!----------------- End -----------------

end module LinearAlgebra