!Statistical routine: data analysis, probability distribution
module Statistics
    use Mathematics; use LinearAlgebra
    implicit none

contains
!Sample variance of data
real*8 function Variance(data, N)
    integer,intent(in)::N
    real*8,dimension(N),intent(in)::data
    Variance=(sum(data*data)-sum(data)**2/dble(N))/dble(N-1)
end function Variance

!R^2 (coefficient of determination) for the prediction of data
real*8 function RSquare(prediction, data, N)
    integer,intent(in)::N
    real*8,dimension(N),intent(in)::prediction,data
    real*8,dimension(N)::dev
    dev=prediction-data; RSquare=1d0-sum(dev*dev)/(sum(data*data)-sum(data)**2/dble(N))
end function RSquare

!N dimensional normal distribution
real*8 function NormalDistribution(x, average, covariance, N)
    integer,intent(in)::N
    real*8,dimension(N),intent(in)::x,average
    real*8,dimension(N,N),intent(in)::covariance
    real*8,dimension(N)::displacement
    real*8,dimension(N,N)::covinv
    displacement=x-average; covinv=covariance; call My_dsytri(covinv,N)
    NormalDistribution=dExp(-0.5d0*dot_product(displacement,matmul(covinv,displacement)))&
                       /sqrt2pi**N/dSqrt(determinant(covariance,N))
end function NormalDistribution

end module Statistics