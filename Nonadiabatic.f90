!Common routines applied in nonadiabatic chemistry
module Nonadiabatic
    implicit none

!Global variable
    !Phase fixing
        integer::Nonadiabatic_NPhasePossibilities,Nonadiabatic_NPhaseDifferencePossibilities
        real*8,allocatable,dimension(:,:)::Nonadiabatic_PhasePossibility,Nonadiabatic_PhaseDifferencePossibility

contains
!Input:  N dimensional ascendingly sorted array energy
!Output: degenerate harvests whether exists almost degenerate energy levels (energy difference < threshold)
subroutine CheckDegeneracy(degenerate,threshold,energy,N)
    logical,intent(out)::degenerate
    real*8,intent(in)::threshold
    integer,intent(in)::N
    real*8,dimension(N),intent(in)::energy
    integer::i
    degenerate=.false.
    do i=1,N-1
        if(energy(i+1)-energy(i)<threshold) then
            degenerate=.true.
            exit
        end if  
    end do
end subroutine CheckDegeneracy

!A, eigval, eigvec satisfy: A . eigvec(:,i) = eigval(i) * eigvec(:,i) for all 1 <= i <= N
!dA = ▽A in A representation is dim x N x N 3rd-order tensor
!Return dim x N x N 3rd-order tensor M satisfying: ▽eigvec = eigvec M (matrix multiplication on N x N)
!Note that M is anti hermitian, so only fill in strictly lower triangle
function deigvec_ByKnowneigval_dA(eigval,dA,dim,N)
    integer,intent(in)::dim,N
    real*8,dimension(N),intent(in)::eigval
    real*8,dimension(dim,N,N),intent(in)::dA
    real*8,dimension(dim,N,N)::deigvec_ByKnowneigval_dA
    integer::i,j
    forall(i=2:N,j=1:N-1,i>j)
        deigvec_ByKnowneigval_dA(:,i,j)=dA(:,i,j)/(eigval(j)-eigval(i))
    end forall
end function deigvec_ByKnowneigval_dA

!--------------------- Phase fixing ---------------------
    !Eigenvector has indeterminate phase, consequently any inner product involving two
    !different states also does not have determinate phase. Sometimes we need to fix it

    !Generate the permutation list of phase and phase difference according to number of states (N)
    subroutine InitializePhaseFixing(N)
        integer,intent(in)::N
        integer::i,j
        !Basis have 2^N possible phases
            Nonadiabatic_NPhasePossibilities=ishft(1,N)-1!Unchanged case is excluded
            if(allocated(Nonadiabatic_PhasePossibility)) deallocate(Nonadiabatic_PhasePossibility)
            allocate(Nonadiabatic_PhasePossibility(N,Nonadiabatic_NPhasePossibilities))
            Nonadiabatic_PhasePossibility(1,1)=-1d0
            Nonadiabatic_PhasePossibility(2:N,1)=1d0
            do i=2,Nonadiabatic_NPhasePossibilities
                Nonadiabatic_PhasePossibility(:,i)=Nonadiabatic_PhasePossibility(:,i-1)
                j=1
                do while(Nonadiabatic_PhasePossibility(j,i)==-1d0)
                    Nonadiabatic_PhasePossibility(j,i)=1d0
                    j=j+1
                end do
                Nonadiabatic_PhasePossibility(j,i)=-1d0
            end do
        !Basis have 2^(N-1) possible phase difference, so we arbitrarily assign 1st basis has phase = 1
            Nonadiabatic_NPhaseDifferencePossibilities=ishft(1,N-1)-1!Unchanged case is excluded
            if(allocated(Nonadiabatic_PhaseDifferencePossibility)) deallocate(Nonadiabatic_PhaseDifferencePossibility)
            allocate(Nonadiabatic_PhaseDifferencePossibility(N,Nonadiabatic_NPhaseDifferencePossibilities))
            Nonadiabatic_PhaseDifferencePossibility(1,1)=1d0
            Nonadiabatic_PhaseDifferencePossibility(2,1)=-1d0
            Nonadiabatic_PhaseDifferencePossibility(3:N,1)=1d0
            do i=2,Nonadiabatic_NPhaseDifferencePossibilities
                Nonadiabatic_PhaseDifferencePossibility(:,i)=Nonadiabatic_PhaseDifferencePossibility(:,i-1)
                j=2
                do while(Nonadiabatic_PhaseDifferencePossibility(j,i)==-1d0)
                    Nonadiabatic_PhaseDifferencePossibility(j,i)=1d0
                    j=j+1
                end do
                Nonadiabatic_PhaseDifferencePossibility(j,i)=-1d0
            end do
    end subroutine InitializePhaseFixing

    !dim x N x N 3-order tensor dH
    !Fix off-diagonal element phase of dH by minimizing its difference from reference
    !difference harvest the minimum || dH - dH_Ref ||_F^2
    subroutine dFixdHPhase(dH,dH_Ref,difference,dim,N)
        integer,intent(in)::dim,N
        real*8,dimension(dim,N,N),intent(inout)::dH
        real*8,dimension(dim,N,N),intent(in)::dH_Ref
        real*8,intent(out)::difference
        integer::indexmin,i,istate,jstate
        real*8::temp
        real*8,dimension(dim)::d
        !Initialize: the difference for the input phase
            indexmin=0
            difference=0d0
            do istate=1,N
                do jstate=istate+1,N
                    d=dH(:,jstate,istate)-dH_Ref(:,jstate,istate)
                    difference=difference+dot_product(d,d)
                end do
            end do
        do i=1,Nonadiabatic_NPhaseDifferencePossibilities!Try out all possibilities
            temp=0d0
            do jstate=2,N!1st basis is assigned to have phase = 1
                d=Nonadiabatic_PhaseDifferencePossibility(jstate,i)*dH(:,jstate,1)-dH_Ref(:,jstate,1)
                temp=temp+dot_product(d,d)
            end do
            do istate=2,N
                do jstate=istate+1,N
                    d=Nonadiabatic_PhaseDifferencePossibility(istate,i)*Nonadiabatic_PhaseDifferencePossibility(jstate,i)*dH(:,jstate,istate)-dH_Ref(:,jstate,istate)
                    temp=temp+dot_product(d,d)
                end do
            end do
            if(temp<difference) then
                indexmin=i
                difference=temp
            end if
        end do
        if(indexmin>0) then!Fix off-diagonals phase = the one with smallest difference
            forall(jstate=2:N)!1st basis is assigned to have phase = 1
                dH(:,jstate,1)=Nonadiabatic_PhaseDifferencePossibility(jstate,indexmin)*dH(:,jstate,1)
            end forall
            forall(istate=2:N-1,jstate=3:N,jstate>istate)
                dH(:,jstate,istate)=Nonadiabatic_PhaseDifferencePossibility(istate,indexmin)*Nonadiabatic_PhaseDifferencePossibility(jstate,indexmin)*dH(:,jstate,istate)
            end forall
        end if
        difference=difference*2d0!Above only counted strictly lower triangle
        do istate=1,N!Diagonals
            d=dH(:,istate,istate)-dH_Ref(:,istate,istate)
            difference=difference+dot_product(d,d)
        end do
    end subroutine dFixdHPhase
    
    !N order matrix H, dim x N x N 3-order tensor dH
    !Fix off-diagonal element phase of dH by minimizing its difference from reference
    !Phases of H off-diagonal elements is then fixed accordingly
    !difference harvest the minimum || dH - dH_Ref ||_F^2
    subroutine dFixHPhaseBydH(H,dH,dH_Ref,difference,dim,N)
        integer,intent(in)::dim,N
        real*8,dimension(N,N),intent(inout)::H
        real*8,dimension(dim,N,N),intent(inout)::dH
        real*8,dimension(dim,N,N),intent(in)::dH_Ref
        real*8,intent(out)::difference
        integer::indexmin,i,istate,jstate
        real*8::temp
        real*8,dimension(dim)::d
        !Initialize: the difference for the input phase
            indexmin=0
            difference=0d0
            do istate=1,N
                do jstate=istate+1,N
                    d=dH(:,jstate,istate)-dH_Ref(:,jstate,istate)
                    difference=difference+dot_product(d,d)
                end do
            end do
        do i=1,Nonadiabatic_NPhaseDifferencePossibilities!Try out all possibilities
            temp=0d0
            do jstate=2,N!1st basis is assigned to have phase = 1
                d=Nonadiabatic_PhaseDifferencePossibility(jstate,i)*dH(:,jstate,1)-dH_Ref(:,jstate,1)
                temp=temp+dot_product(d,d)
            end do
            do istate=2,N
                do jstate=istate+1,N
                    d=Nonadiabatic_PhaseDifferencePossibility(istate,i)*Nonadiabatic_PhaseDifferencePossibility(jstate,i)*dH(:,jstate,istate)-dH_Ref(:,jstate,istate)
                    temp=temp+dot_product(d,d)
                end do
            end do
            if(temp<difference) then
                indexmin=i
                difference=temp
            end if
        end do
        if(indexmin>0) then!Fix off-diagonals phase = the one with smallest difference
            forall(jstate=2:N)!1st basis is assigned to have phase = 1
                 H(  jstate,1)=Nonadiabatic_PhaseDifferencePossibility(jstate,indexmin)* H(  jstate,1)
                dH(:,jstate,1)=Nonadiabatic_PhaseDifferencePossibility(jstate,indexmin)*dH(:,jstate,1)
            end forall
            forall(istate=2:N-1,jstate=3:N,jstate>istate)
                 H(  jstate,istate)=Nonadiabatic_PhaseDifferencePossibility(istate,indexmin)*Nonadiabatic_PhaseDifferencePossibility(jstate,indexmin)* H(  jstate,istate)
                dH(:,jstate,istate)=Nonadiabatic_PhaseDifferencePossibility(istate,indexmin)*Nonadiabatic_PhaseDifferencePossibility(jstate,indexmin)*dH(:,jstate,istate)
            end forall
        end if
        difference=difference*2d0!Above only counted strictly lower triangle
        do istate=1,N!Diagonals
            d=dH(:,istate,istate)-dH_Ref(:,istate,istate)
            difference=difference+dot_product(d,d)
        end do
    end subroutine dFixHPhaseBydH
    
    !N order matrix phi, dim x N x N 3-order tensor dH
    !Fix off-diagonal element phase of dH by minimizing its difference from reference
    !Phase of basis phi is then assigned accordingly
    !Warning: 1st basis is arbitrarily assigned phase = 1 because dH can only determine phase difference
    !difference harvest the minimum || dH - dH_Ref ||_F^2
    subroutine dAssignBasisPhaseBydH(phi,dH,dH_Ref,difference,dim,N)
        integer,intent(in)::dim,N
        real*8,dimension(N,N),intent(inout)::phi
        real*8,dimension(dim,N,N),intent(inout)::dH
        real*8,dimension(dim,N,N),intent(in)::dH_Ref
        real*8,intent(out)::difference
        integer::indexmin,i,istate,jstate
        real*8::temp
        real*8,dimension(dim)::d
        !Initialize: the difference for the input phase
            indexmin=0
            difference=0d0
            do istate=1,N
                do jstate=istate+1,N
                    d=dH(:,jstate,istate)-dH_Ref(:,jstate,istate)
                    difference=difference+dot_product(d,d)
                end do
            end do
        do i=1,Nonadiabatic_NPhaseDifferencePossibilities!Try out all possibilities
            temp=0d0
            do jstate=2,N!1st basis is assigned to have phase = 1
                d=Nonadiabatic_PhaseDifferencePossibility(jstate,i)*dH(:,jstate,1)-dH_Ref(:,jstate,1)
                temp=temp+dot_product(d,d)
            end do
            do istate=2,N
                do jstate=istate+1,N
                    d=Nonadiabatic_PhaseDifferencePossibility(istate,i)*Nonadiabatic_PhaseDifferencePossibility(jstate,i)*dH(:,jstate,istate)-dH_Ref(:,jstate,istate)
                    temp=temp+dot_product(d,d)
                end do
            end do
            if(temp<difference) then
                indexmin=i
                difference=temp
            end if
        end do
        if(indexmin>0) then!Assign phase = the one with smallest difference
            forall(jstate=2:N)!1st basis is assigned to have phase = 1
                phi(:,jstate)=Nonadiabatic_PhaseDifferencePossibility(jstate,indexmin)*phi(:,jstate)
                dH(:,jstate,1)=Nonadiabatic_PhaseDifferencePossibility(jstate,indexmin)*dH(:,jstate,1)
            end forall
            forall(istate=2:N-1,jstate=3:N,jstate>istate)
                dH(:,jstate,istate)=Nonadiabatic_PhaseDifferencePossibility(istate,indexmin)*Nonadiabatic_PhaseDifferencePossibility(jstate,indexmin)*dH(:,jstate,istate)
            end forall
        end if
        difference=difference*2d0!Above only counted strictly lower triangle
        do istate=1,N!Diagonals
            d=dH(:,istate,istate)-dH_Ref(:,istate,istate)
            difference=difference+dot_product(d,d)
        end do
    end subroutine dAssignBasisPhaseBydH
    
    !N order matrix H & phi, dim x N x N 3-order tensor dH
    !Fix off-diagonal element phase of dH by minimizing its difference from reference
    !Phases of H off-diagonal elements is then fixed accordingly
    !Phase of basis phi is then assigned accordingly
    !Warning: 1st basis is arbitrarily assigned phase = 1 because dH can only determine phase difference
    !difference harvest the minimum || dH - dH_Ref ||_F^2
    subroutine dFixHPhase_AssignBasisPhaseBydH(H,phi,dH,dH_Ref,difference,dim,N)
        integer,intent(in)::dim,N
        real*8,dimension(N,N),intent(inout)::H,phi
        real*8,dimension(dim,N,N),intent(inout)::dH
        real*8,dimension(dim,N,N),intent(in)::dH_Ref
        real*8,intent(out)::difference
        integer::indexmin,i,istate,jstate
        real*8::temp
        real*8,dimension(dim)::d
        !Initialize: the difference for the input phase
            indexmin=0
            difference=0d0
            do istate=1,N
                do jstate=istate+1,N
                    d=dH(:,jstate,istate)-dH_Ref(:,jstate,istate)
                    difference=difference+dot_product(d,d)
                end do
            end do
        do i=1,Nonadiabatic_NPhaseDifferencePossibilities!Try out all possibilities
            temp=0d0
            do jstate=2,N!1st basis is assigned to have phase = 1
                d=Nonadiabatic_PhaseDifferencePossibility(jstate,i)*dH(:,jstate,1)-dH_Ref(:,jstate,1)
                temp=temp+dot_product(d,d)
            end do
            do istate=2,N
                do jstate=istate+1,N
                    d=Nonadiabatic_PhaseDifferencePossibility(istate,i)*Nonadiabatic_PhaseDifferencePossibility(jstate,i)*dH(:,jstate,istate)-dH_Ref(:,jstate,istate)
                    temp=temp+dot_product(d,d)
                end do
            end do
            if(temp<difference) then
                indexmin=i
                difference=temp
            end if
        end do
        if(indexmin>0) then!Fix phase = the one with smallest difference
            forall(jstate=2:N)!1st basis is assigned to have phase = 1
                phi(:,jstate) =Nonadiabatic_PhaseDifferencePossibility(jstate,indexmin)*phi(:,jstate)
                 H(  jstate,1)=Nonadiabatic_PhaseDifferencePossibility(jstate,indexmin)* H(  jstate,1)
                dH(:,jstate,1)=Nonadiabatic_PhaseDifferencePossibility(jstate,indexmin)*dH(:,jstate,1)
            end forall
            forall(istate=2:N-1,jstate=3:N,jstate>istate)
                 H(  jstate,istate)=Nonadiabatic_PhaseDifferencePossibility(istate,indexmin)*Nonadiabatic_PhaseDifferencePossibility(jstate,indexmin)* H(  jstate,istate)
                dH(:,jstate,istate)=Nonadiabatic_PhaseDifferencePossibility(istate,indexmin)*Nonadiabatic_PhaseDifferencePossibility(jstate,indexmin)*dH(:,jstate,istate)
            end forall
        end if
        difference=difference*2d0!Above only counted strictly lower triangle
        do istate=1,N!Diagonals
            d=dH(:,istate,istate)-dH_Ref(:,istate,istate)
            difference=difference+dot_product(d,d)
        end do
    end subroutine dFixHPhase_AssignBasisPhaseBydH
!------------------------- End --------------------------

end module Nonadiabatic