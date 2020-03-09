!General basic derived type, general basic routine
module General
    implicit none

!Derived type
    !Analog to 2nd pointer in C
    type i2PArray
        integer,allocatable,dimension(:)::Array
    end type i2PArray

    type d2PArray
        real*8,allocatable,dimension(:)::Array
    end type d2PArray

    type d2PMatrix
        real*8,allocatable,dimension(:,:)::Matrix
    end type d2PMatrix

!Overload
    interface QuickSort
        module procedure dQuickSort, iQuickSort
    end interface QuickSort

    interface MergeSort
        module procedure dMergeSort, iMergeSort
    end interface MergeSort

contains
subroutine ShowTime()!Show date hour minute second
    integer,dimension(8)::time
    call date_and_time(values=time)
    write(*,'(1x,I4,1x,A4,1x,I2,1x,A5,1x,I2,1x,A3,1x,I2,A1,I2,A1,I2)')time(1),'year',time(2),'month',time(3),'day',time(5),':',time(6),':',time(7)
end subroutine ShowTime

subroutine dScientificNotation(x, i)! x_input = x_output * 10^i
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
    !Return an N dimensional array filled with 1d0
    function ones(N)
        integer,intent(in)::N
        real*8,dimension(N)::ones
        ones=1d0
    end function ones
    
    !Return an N order unit matrix
    function UnitMatrix(N)
        integer,intent(in)::N
        real*8,dimension(N,N)::UnitMatrix
        integer::i
        UnitMatrix=0d0; forall(i=1:N); UnitMatrix(i,i)=1d0; end forall
    end function UnitMatrix
    
    !Input:  N dimensional vector x
    !Output: N order diagonal matrix with main diagonal vector = x
    function diag(x, N)
        integer,intent(in)::N
        real*8,dimension(N),intent(in)::x
        real*8,dimension(N,N)::diag
        integer::i
        diag=0d0; forall(i=1:N); diag(i,i)=x(i); end forall
    end function diag
!--------------- End ---------------

!---------- Random number ----------
    !A better version of random_seed
    subroutine BetterRandomSeed()
        integer::i,j
        integer,dimension(8)::time
        integer,allocatable,dimension(:)::seed
        real*8::flagd
        call random_seed(size=j)
        allocate(seed(j))
        call date_and_time(values=time)
        seed=time(8)+37d0*(/(i-1,i=1,j)/)
        call random_seed(put=seed)
        deallocate(seed)
        do i=1,time(6)*3600d0+time(7)*60d0+time(8)
            call random_number(flagd)
        end do
    end subroutine BetterRandomSeed

    !Make random_number more random
    real*8 function BetterRandomNumber()
        integer::i
        integer,dimension(8)::time
        real*8::harvest
        call date_and_time(values=time)
        do i=1,max(time(8)/100,2)
          call random_number(harvest)
        end do
        BetterRandomNumber=harvest
    end function BetterRandomNumber

    !Generate a random number in [-1,1)
    real*8 function random_number_minus1to1()
        call random_number(random_number_minus1to1)
        random_number_minus1to1=2d0*random_number_minus1to1-1d0
    end function
    !More random
    real*8 function BetterRandomNumber_minus1to1()
        BetterRandomNumber_minus1to1=BetterRandomNumber()
        BetterRandomNumber_minus1to1=2d0*BetterRandomNumber_minus1to1-1d0
    end function

    !Gaussian random number generator using box-muller method
    !Ref: http://en.wikipedia.org/wiki/gaussian_random_variable
    real*8 function GaussianRandomNumber(mean, sigma)
        real*8,intent(in)::mean,sigma
        real*8::r1,r2
        call random_number(r1)
        call random_number(r2)
        GaussianRandomNumber=mean+sigma*dsqrt(-2d0*dlog(r1))*dcos(6.283185307179586d0*r2)
    end function GaussianRandomNumber
    !More random
    real*8 function BetterGaussianRandomNumber(mean, sigma)
        real*8,intent(in)::mean,sigma
        real*8::r1,r2
        r1=BetterRandomNumber()
        r2=BetterRandomNumber()
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
            r1(1)=BetterRandomNumber()
            r1(2)=BetterRandomNumber()
            r1=r1*2d0-1d0
            if(dot_product(r1,r1)<1d0) exit
        end do
        do while(.true.)
            r2(1)=BetterRandomNumber()
            r2(2)=BetterRandomNumber()
            r2=r2*2d0-1d0
            if(dot_product(r2,r2)<1d0) exit
        end do
        BetterRandomUnitQuaternion(1:2)=r1
        BetterRandomUnitQuaternion(3:4)=r2*Sqrt((1d0-dot_product(r1,r1))/dot_product(r2,r2))
    end function BetterRandomUnitQuaternion
!--------------- End ---------------

!------------- Sorting -------------
    !Two classical O(NlogN) sorting methods are implemented:
    !    Quick sort is faster, but unstable
    !    Merge sort is stable, but slower

    !Inputs: N dimensional vector item & indices, item contains the elements to be sorted, indices(i)=i
    !        sort from item(first) to item(last), inclusive
    !Returns: item (modified) contains elements in ascending order
    !         indices(i) is the original index of sorted item(i)
    !Double item
    recursive subroutine dQuickSort(item, first, last, indices, N)
        integer,intent(in)::N,first,last
        real*8,dimension(N),intent(inout)::item
        integer,dimension(N),intent(inout)::indices
        integer::mid
        if (First<Last) then
           call split(item,first,last,mid,indices,N)
           call dQuickSort(item,first,mid-1,indices,N); call dQuickSort(item,mid+1,last,indices,N)
        end if
        contains
        subroutine split(item,low,high,mid,indices,N)
            integer,intent(in)::N,low,high
            real*8,dimension(N),intent(inout)::Item
            integer,intent(inout)::mid
            integer,dimension(N),intent(inout)::indices
            integer::left,right,iPivot,iSwap
            real*8::pivot,swap
            left=low; right=high; pivot=item(low); iPivot=indices(low)
            do while(left<right)
                do while(left<right.and.item(right)>=pivot)
                    right=right-1
                end do
                do while(left<right.and.item(left)<=pivot)
                    left=left+1
                end do
                if(left<right) then
                    swap=item(left); item(left)=item(right); item(right)=swap
                    iSwap=indices(left); indices(left)=indices(right); indices(right)=iSwap
                end if
            end do
            item(low)=item(right); item(right)=pivot; mid=right
            indices(low)=indices(right); indices(right)=iPivot
        end subroutine split
    end subroutine dQuickSort
    !Integer item
    recursive subroutine iQuickSort(item, first, last, indices, N)
        integer,intent(in)::N,first,last
        integer,dimension(N),intent(inout)::item,indices
        integer::mid
        if (First<Last) then
           call split(item,first,last,mid,indices,N)
           call iQuickSort(item,first,mid-1,indices,N); call iQuickSort(item,mid+1,last,indices,N)
        end if
        contains
        subroutine split(item,low,high,mid,indices,N)
            integer,intent(in)::N,low,high
            integer,dimension(N),intent(inout)::item,indices
            integer,intent(inout)::mid
            integer::left,right,iPivot,iSwap,pivot,swap
            left=low; right=high; pivot=item(low); iPivot=indices(low)
            do while(left<right)
                do while(left<right.and.item(right)>=pivot)
                    right=right-1
                end do
                do while(left<right.and.item(left)<=pivot)
                    left=left+1
                end do
                if(left<right) then
                    swap=item(left); item(left)=item(right); item(right)=swap
                    iSwap=indices(left); indices(left)=indices(right); indices(right)=iSwap
                end if
            end do
            item(low)=item(right); item(right)=pivot; mid=right
            indices(low)=indices(right); indices(right)=iPivot
        end subroutine split
    end subroutine iQuickSort
    
    !Inputs: N dimensional vector item, item contains the elements to be sorted
    !        sort from item(first) to item(last), inclusive
    !Returns: item (modified) contains elements in ascending order
    !         NRP harvests the number of reverse pairs
    !Double item
    subroutine dMergeSort(item, first, last, NRP, N)
        integer,intent(in)::N,first,last
        real*8,dimension(N),intent(inout)::item
        integer,intent(out)::NRP
        integer::length
        real*8,allocatable,dimension(:)::work
        NRP=0; length=last-first+1; allocate(work((length+1)/2))
        call sort(item(first:last),length,work)
        contains
        recursive subroutine sort(A, N, work)
            integer,intent(in)::N
            real*8,dimension(N),intent(inout)::A
            real*8,dimension((N+1)/2),intent(out)::work
            integer::NA,NB
            real*8::swap
            if(N<2) return
            if(N==2) then
                if(A(1)>A(2)) then
                    swap=A(1); A(1)=A(2); A(2)=swap; NRP=NRP+1
                end if
                return
            end if
            NA=(N+1)/2; NB=N-NA
            call sort(A(1:NA),NA,work); call sort(A(NA+1:N),NB,work)
            if(A(NA)>A(NA+1)) then
                work(1:NA)=A(1:NA); call merge(work(1:NA),NA,A(NA+1:N),NB,A,N)
            end if
        end subroutine sort
        subroutine merge(A, NA, B, NB, C, NC)
            integer,intent(in):: NA,NB,NC!Normal usage: NA + NB = NC
            real*8,dimension(NA),intent(inout)::A
            real*8,dimension(NB),intent(in)::B
            real*8,dimension(NC),intent(inout)::C!B overlays C(NA+1:NC)
            integer::i,j,k
            i=1; j=1; k=1
            do while(i<=NA.and.j<=NB)
                if(A(i)<=B(j)) then
                    C(k)=A(i); i=i+1
                else
                    C(k)=B(j); j=j+1; NRP=NRP+NA-i+1
                end if
                k=k+1
            end do
            do while(i<=NA)!B overlays C(NA+1:NC), no need to account for case i reaches NA earlier
               C(k)=A(i); i=i+1; k=k+1
            end do
        end subroutine merge
    end subroutine dMergeSort
    !Integer item
    subroutine iMergeSort(item, first, last, NRP, N)
        integer,intent(in)::N,first,last
        integer,dimension(N),intent(inout)::item
        integer,intent(out)::NRP
        integer::length
        integer,allocatable,dimension(:)::work
        NRP=0; length=last-first+1; allocate(work((length+1)/2))
        call sort(item(first:last),length,work)
        contains
        recursive subroutine sort(A, N, work)
            integer,intent(in)::N
            integer,dimension(N),intent(inout)::A
            integer,dimension((N+1)/2),intent(out)::work
            integer::NA,NB,swap
            if(N<2) return
            if(N==2) then
                if(A(1)>A(2)) then
                    swap=A(1); A(1)=A(2); A(2)=swap; NRP=NRP+1
                end if
                return
            end if
            NA=(N+1)/2; NB=N-NA
            call sort(A(1:NA),NA,work); call sort(A(NA+1:N),NB,work)
            if(A(NA)>A(NA+1)) then
                work(1:NA)=A(1:NA); call merge(work(1:NA),NA,A(NA+1:N),NB,A,N)
            end if
        end subroutine sort
        subroutine merge(A, NA, B, NB, C, NC)
            integer,intent(in):: NA,NB,NC!Normal usage: NA + NB = NC
            integer,dimension(NA),intent(inout)::A
            integer,dimension(NB),intent(in)::B
            integer,dimension(NC),intent(inout)::C!B overlays C(NA+1:NC)
            integer::i,j,k
            i=1; j=1; k=1
            do while(i<=NA.and.j<=NB)
                if(A(i)<=B(j)) then
                    C(k)=A(i); i=i+1
                else
                    C(k)=B(j); j=j+1; NRP=NRP+NA-i+1
                end if
                k=k+1
            end do
            do while(i<=NA)!B overlays C(NA+1:NC), no need to account for case i reaches NA earlier
               C(k)=A(i); i=i+1; k=k+1
            end do
        end subroutine merge
    end subroutine iMergeSort
!--------------- End ---------------

end module General