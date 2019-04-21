!General basic derived type, general basic routine
module General
    implicit none

!Derived type
    !Analog to double 2nd pointer in C
    !Example: type(d2PArray),allocatable,dimension(:)::A
    !         A(i).Array(j) is the same to A[i][j] in C (forget about starting from 0 or 1)
    type d2PArray
        real*8,allocatable,dimension(:)::Array
    end type d2PArray

contains
!Show date hour minute second
subroutine ShowTime()
    integer::value(1:8)
    call date_and_time(values=value)
    write(*,*)value(3),'d',value(5),':',value(6),':',value(7)
end subroutine ShowTime

! x_input = x_output * 10^i
subroutine dScientificNotation(x,i)
    real*8,intent(inout)::x
    integer,intent(out)::i
    i=0
    do
        if(x<1d0) then
            x=x*10d0
            i=i-1
        else
            exit
        end if
    end do
    do
        if(x>=10d0) then
            x=x/10d0
            i=i+1
        else
            exit
        end if
    end do
end subroutine dScientificNotation

!---------- Tensor object ----------
    !Return a N dimensional array filled with 1d0
    function ones(N)
        integer,intent(in)::N
        real*8,dimension(N)::ones
        ones=1d0
    end function ones
    
    !Return a N order unit matrix
    function UnitMatrix(N)
        integer,intent(in)::N
        real*8,dimension(N,N)::UnitMatrix
        integer::i
        UnitMatrix=0d0
        forall(i=1:N)
            UnitMatrix(i,i)=1d0
        end forall
    end function UnitMatrix
    
    !Input:  N dimensional vector x
    !Output: N order diagonal matrix with main diagonal vector = x
    function diag(x,N)
        integer,intent(in)::N
        real*8,dimension(N),intent(in)::x
        real*8,dimension(N,N)::diag
        integer::i
        diag=0d0
        forall(i=1:N)
            diag(i,i)=x(i)
        end forall
    end function diag
!--------------- End ---------------

!---------- Random number ----------
    !A better version of random_seed
    subroutine BetterRandomSeed()
        integer ii,nn,value(1:8)
        integer,allocatable :: seed(:)
        double precision flagd
        call random_seed(size=nn)
        allocate(seed(nn))
        call date_and_time(values=value)
        seed = value(8)+37*(/(ii-1,ii=1,nn)/)
        call random_seed(put=seed)
        deallocate(seed)
        do ii=1,value(6)*3600+value(7)*60+value(8)
            call random_number(flagd)
        end do
    end subroutine BetterRandomSeed

    !Make random_number more random
    subroutine BetterRandomNumber(harvest)
        implicit none
        integer ii,value(1:8)
        real*8::harvest
        call date_and_time(values=value)
        do ii=1,max(value(8)/100,2)
          call random_number(harvest)
        enddo
    end subroutine BetterRandomNumber

    !Generate a random number in [-1,1)
    real*8 function random_number_minus1to1()
        call random_number(random_number_minus1to1)
        random_number_minus1to1=2d0*random_number_minus1to1-1d0
    end function
    !More random
    real*8 function BetterRandomNumber_minus1to1()
        call BetterRandomNumber(BetterRandomNumber_minus1to1)
        BetterRandomNumber_minus1to1=2d0*BetterRandomNumber_minus1to1-1d0
    end function

    !Gaussian random number generator using box-muller method
    !Ref: http://en.wikipedia.org/wiki/gaussian_random_variable
    real*8 function GaussianRandomNumber(mean,sigma)
        real*8,intent(in)::mean,sigma
        real*8::r1,r2
        call random_number(r1)
        call random_number(r2)
        GaussianRandomNumber=mean+sigma*dsqrt(-2d0*dlog(r1))*dcos(6.283185307179586d0*r2)
    end function GaussianRandomNumber
    !More random
    real*8 function BetterGaussianRandomNumber(mean,sigma)
        real*8,intent(in)::mean,sigma
        real*8::r1,r2
        call BetterRandomNumber(r1)
        call BetterRandomNumber(r2)
        BetterGaussianRandomNumber=mean+sigma*dsqrt(-2d0*dlog(r1))*dcos(6.283185307179586d0*r2)
    end function BetterGaussianRandomNumber

    !Return a random unit quaternion
    function RandomUnitQuaternion()
        real*8,dimension(4)::RandomUnitQuaternion
        real*8,dimension(2)::r1,r2
        do while(.true.)
            call random_number(r1)
            r1=r1*2d0-1d0
            if(dot_product(r1,r1)<1d0) exit
        end do
        do while(.true.)
            call random_number(r2)
            r2=r2*2d0-1d0
            if(dot_product(r2,r2)<1d0) exit
        end do
        RandomUnitQuaternion(1:2)=r1
        RandomUnitQuaternion(3:4)=r2*Sqrt((1d0-dot_product(r1,r1))/dot_product(r2,r2))
    end function RandomUnitQuaternion
    !More random
    function BetterRandomUnitQuaternion()
        real*8,dimension(4)::BetterRandomUnitQuaternion
        real*8,dimension(2)::r1,r2
        do while(.true.)
            call BetterRandomNumber(r1(1))
            call BetterRandomNumber(r1(2))
            r1=r1*2d0-1d0
            if(dot_product(r1,r1)<1d0) exit
        end do
        do while(.true.)
            call BetterRandomNumber(r2(1))
            call BetterRandomNumber(r2(2))
            r2=r2*2d0-1d0
            if(dot_product(r2,r2)<1d0) exit
        end do
        BetterRandomUnitQuaternion(1:2)=r1
        BetterRandomUnitQuaternion(3:4)=r2*Sqrt((1d0-dot_product(r1,r1))/dot_product(r2,r2))
    end function BetterRandomUnitQuaternion
!--------------- End ---------------

!------------- Sorting -------------
    !Inputs: N order vector item & indices, item contains the elements to be sorted, indices(i)=i
    !        sort from item(first) to item(last), inclusive
    !Returns: item (modified) contains elements in ascending order
    !         indices(i) is the original index of sorted item(i)
    recursive subroutine dQuickSort(item,first,last,indices,N)
        integer,intent(in)::N
        real*8,dimension(N),intent(inout)::item
        integer,intent(in)::first,last
        integer,dimension(N),intent(inout)::indices
        integer::mid
        if (First<Last) then
           call dSplit(item,first,last,mid,indices,N)
           call dQuickSort(item,first,mid-1,indices,N)
           call dQuickSort(item,mid+1,last,indices,N)
        end if
        contains
            subroutine dSplit(item,low,high,mid,indices,N)
                integer,intent(in)::N
                real*8,dimension(N),intent(inout)::Item
                integer,intent(in)::low,high
                integer,intent(inout)::mid
                integer,dimension(N),intent(inout)::indices
                integer::left,right,iPivot,iSwap
                real*8::pivot,swap
                left=low
                right=high
                pivot=item(low)
                iPivot=indices(low)
                do while ( left < right )
                    do while ( left < right .and. item(right) >= pivot )
                        right = right - 1
                    end do
                    do while ( left < right .and. item(left) <= pivot )
                        left = left + 1
                    end do
                    if (left < right) then
                        swap        = item(left)
                        item(left)  = item(right)
                        item(right) = swap
                        iSwap          = indices(left)
                        indices(left)  = indices(right)
                        indices(right) = iSwap
                    end if
                end do
                item(low)   = item(right)
                item(right) = pivot
                mid         = right
                indices(low)   = indices(right)
                indices(right) = iPivot
            end subroutine dSplit
    end subroutine dQuickSort
    
    !Inputs: N order vector item & indices, item contains the elements to be sorted, indices(i)=i
    !        sort from item(first) to item(last), inclusive
    !Returns: item (modified) contains elements in ascending order
    !         indices(i) is the original index of sorted item(i)
    recursive subroutine iQuickSort(item,first,last,indices,N)
        integer,intent(in)::N
        integer,dimension(N),intent(inout)::item
        integer,intent(in)::first,last
        integer,dimension(N),intent(inout)::indices
        integer::mid
        if (First<Last) then
           call iSplit(item,first,last,mid,indices,N)
           call iQuickSort(item,first,mid-1,indices,N)
           call iQuickSort(item,mid+1,last,indices,N)
        end if
        contains
            subroutine iSplit(item,low,high,mid,indices,N)
                integer,intent(in)::N
                integer,dimension(N),intent(inout)::Item
                integer,intent(in)::low,high
                integer,intent(inout)::mid
                integer,dimension(N),intent(inout)::indices
                integer::left,right,iPivot,iSwap,pivot,swap
                left=low
                right=high
                pivot=item(low)
                iPivot=indices(low)
                do while ( left < right )
                    do while ( left < right .and. item(right) >= pivot )
                        right = right - 1
                    end do
                    do while ( left < right .and. item(left) <= pivot )
                        left = left + 1
                    end do
                    if (left < right) then
                        swap        = item(left)
                        item(left)  = item(right)
                        item(right) = swap
                        iSwap          = indices(left)
                        indices(left)  = indices(right)
                        indices(right) = iSwap
                    end if
                end do
                item(low)   = item(right)
                item(right) = pivot
                mid         = right
                indices(low)   = indices(right)
                indices(right) = iPivot
            end subroutine iSplit
    end subroutine iQuickSort
!--------------- End ---------------

end module General