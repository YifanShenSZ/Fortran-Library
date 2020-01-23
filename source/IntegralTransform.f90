!Integral transform routine, currently only Fourier transform available
!
!Instruction:
!Fourier transform and its inverse is defined by: (This is unitary transformation definition)
!    F(k) = Integrate[ exp(-ikx) * f(x), { x, -Infinity, Infinity }] / Sqrt(2pi)
!    f(x) = Integrate[ exp( ikx) * F(k), { k, -Infinity, Infinity }] / Sqrt(2pi)
!This has to be carried out through obtaining enough data points then numerically integrating
!For discrete data points subject to periodic boundary condition, we may define discrete Fourier transform:
!    F(w_j) = Sum[ f(t_i) * exp(-i w_j t_i), { i, 0, N-1 }]
!    f(t_i) = Sum[ F(w_j) * exp( i w_j t_i), { j, 0, N-1 }]
!This definition is an interpolation. Complexity = O(N^2), for N = 2^n can be reduced to O(NlogN)
!The input data set consists of N uniformly spaced data points {x_i, f(x_i)}, x_i+1 - x_i = dx
!Periodic boundary condition is assumed as f(x_I) = f(x_[I mod N]) for all I < 0 or I >= N
!From periodic boundary condition, there can only be N different frequencies: 2pi/Ndx * {0, 1, ..., N-1}
module IntegralTransform
    use Mathematics
    use mkl_dfti!For FFT, disable it if you have no access to MKL
    implicit none

contains
!Fourier transform lx points (x,psy) to lk points (k,phi). x must be uniformly spaced
!Input: x, psy, k. Output: phi
subroutine dFourierTransform(x,psy,lx,k,phi,lk)
    integer,intent(in)::lx,lk
    real*8,dimension(lx),intent(in)::x
    complex*16,dimension(lx),intent(in)::psy
    real*8,dimension(lk),intent(in)::k
    complex*16,dimension(lk),intent(out)::phi
    integer::i,j
    do j=1,lk
        phi(j)=exp(-ci*k(j)*x(1))*psy(1)
        do i=2,lx
            phi(j)=phi(j)+exp(-ci*k(j)*x(i))*psy(i)
        end do
    end do
    phi=phi*(x(2)-x(1))/sqrt2pi
end subroutine dFourierTransform

!Inverse Fourier transform lk points (k,phi) to lx points (x,psy). k must be uniformly spaced
!Input: k, phi, x. Output: psy
subroutine dInverseFourierTransform(k,phi,lk,x,psy,lx)
    integer,intent(in)::lk,lx
    real*8,dimension(lk),intent(in)::k
    complex*16,dimension(lk),intent(in)::phi
    real*8,dimension(lx),intent(in)::x
    complex*16,dimension(lx),intent(out)::psy
    integer::i,j
    do j=1,lx
        psy(j)=exp(ci*k(1)*x(j))*phi(1)
        do i=2,lk
            psy(j)=psy(j)+exp(ci*k(i)*x(j))*phi(i)
        end do
    end do
    psy=psy*(k(2)-k(1))/sqrt2pi
end subroutine dInverseFourierTransform

!Fast fourier transform 2^n data points
!On input psy contains {psy(x_i)}, on exit psy harvests {phi(k_j)}
!Although no explicit output, {k_j} = 2pi / (2^n * dx) * {0, 1, ..., 2^n-1}
subroutine FFT(psy,n)
    integer,intent(in)::n
    complex*16,dimension(2**n),intent(inout)::psy
    type(DFTI_DESCRIPTOR),pointer::handle
    integer::status
    !Step1: create descriptor, 4 arguments are:
    !    1, handle (some mkl handle storing everything)
    !    2, precision of the transform, available: DFTI_SINGLE, DFTI_DOUBLE
    !    3, dimension of the transform (this is a 1D wrapper, you have to organize your own loop for higher dimension)
    !    4, length of the data set
    status=DftiCreateDescriptor(handle,DFTI_DOUBLE,DFTI_COMPLEX,1,2**n)
    status=DftiCommitDescriptor(handle)!Step2: commit the handle
    !Step3: compute the transformation. For backward, use DftiComputeBackward
    status=DftiComputeForward(handle,psy)
    status=DftiFreeDescriptor(handle)!Step4: free space
end subroutine FFT

end module IntegralTransform