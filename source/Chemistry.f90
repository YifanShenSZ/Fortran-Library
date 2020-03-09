!Common routines applied in chemistry
module Chemistry
    use Mathematics
    implicit none

!Global variable
    !Phase fixing
        integer::Chemistry_NPhasePossibilities,Chemistry_NPhaseDifferencePossibilities
        real*8,allocatable,dimension(:,:)::Chemistry_PhasePossibility,Chemistry_PhaseDifferencePossibility

contains
integer function Symbol2Number(element)!Return element number based on input element symbol
    character*2,intent(in)::element
    select case(element)
        case('H' ); Symbol2Number=1
        case('He'); Symbol2Number=2
        case('Li'); Symbol2Number=3
        case('Be'); Symbol2Number=4
        case('B' ); Symbol2Number=5
        case('C' ); Symbol2Number=6
        case('N' ); Symbol2Number=7
        case('O' ); Symbol2Number=8
        case('F' ); Symbol2Number=9
        case('Ne'); Symbol2Number=10
        case('Na'); Symbol2Number=11
        case('Mg'); Symbol2Number=12
        case('Al'); Symbol2Number=13
        case('Si'); Symbol2Number=14
        case('P' ); Symbol2Number=15
        case('S' ); Symbol2Number=16
        case('Cl'); Symbol2Number=17
        case('Ar'); Symbol2Number=18
        case('K' ); Symbol2Number=19
        case('Ca'); Symbol2Number=20
        case('Sc'); Symbol2Number=21
        case('Ti'); Symbol2Number=22
        case('V' ); Symbol2Number=23
        case('Cr'); Symbol2Number=24
        case('Mn'); Symbol2Number=25
        case('Fe'); Symbol2Number=26
        case('Co'); Symbol2Number=27
        case('Ni'); Symbol2Number=28
        case('Cu'); Symbol2Number=29
        case('Zn'); Symbol2Number=30
        case('Ga'); Symbol2Number=31
        case('Ge'); Symbol2Number=32
        case('As'); Symbol2Number=33
        case('Se'); Symbol2Number=34
        case('Br'); Symbol2Number=35
        case('Kr'); Symbol2Number=36
        case default
            write(*,'(1x,A42,1x,A2)')'Program abort: unsupported element symbol:',element
            write(*,'(1x,A37)')'Current version only supports H to Kr'
            stop
    end select
end function Symbol2Number
character*2 function Number2Symbol(element)!Vice versa
    integer,intent(in)::element
    select case(element)
        case(1 ); Number2Symbol='H'
        case(2 ); Number2Symbol='He'
        case(3 ); Number2Symbol='Li'
        case(4 ); Number2Symbol='Be'
        case(5 ); Number2Symbol='B'
        case(6 ); Number2Symbol='C'
        case(7 ); Number2Symbol='N'
        case(8 ); Number2Symbol='O'
        case(9 ); Number2Symbol='F'
        case(10); Number2Symbol='Ne'
        case(11); Number2Symbol='Na'
        case(12); Number2Symbol='Mg'
        case(13); Number2Symbol='Al'
        case(14); Number2Symbol='Si'
        case(15); Number2Symbol='P'
        case(16); Number2Symbol='S'
        case(17); Number2Symbol='Cl'
        case(18); Number2Symbol='Ar'
        case(19); Number2Symbol='K'
        case(20); Number2Symbol='Ca'
        case(21); Number2Symbol='Sc'
        case(22); Number2Symbol='Ti'
        case(23); Number2Symbol='V'
        case(24); Number2Symbol='Cr'
        case(25); Number2Symbol='Mn'
        case(26); Number2Symbol='Fe'
        case(27); Number2Symbol='Co'
        case(28); Number2Symbol='Ni'
        case(29); Number2Symbol='Cu'
        case(30); Number2Symbol='Zn'
        case(31); Number2Symbol='Ga'
        case(32); Number2Symbol='Ge'
        case(33); Number2Symbol='As'
        case(34); Number2Symbol='Se'
        case(35); Number2Symbol='Br'
        case(36); Number2Symbol='Kr'
        case default
            write(*,'(1x,A42,1x,I3)')'Program abort: unsupported element number:',element
            write(*,'(1x,A37)')'Current version only supports 1 to 36'
            stop
    end select
end function Number2Symbol

!Write molecular structure and vibration to an Avogadro format file
!Required: symbol: element symbol of each atom
!          structure: Cartesian coordinate of each atom in Angstrom
!          freq: vibrational angular frequencies in cm^-1
!          mode: Cartesian normal modes
!Optional: file: (default = 'avogadro.log') output file name
subroutine Avogadro_Vibration(NAtoms, symbol, structure, vibdim, freq, mode, file)
    integer,intent(in)::NAtoms,vibdim
    character*2,dimension(NAtoms),intent(in)::symbol
    real*8,dimension(3,NAtoms),intent(in)::structure
    real*8,dimension(vibdim),intent(in)::freq
    real*8,dimension(3,NAtoms,vibdim),intent(in)::mode
    character(*),optional,intent(in)::file
    integer::i; integer,dimension(NAtoms)::number
    if(present(file)) then; open(unit=99,file=file,status='replace')
    else; open(unit=99,file='avogadro.log',status='replace'); end if
        write(99,'(A29)')'---------- Comment ----------'
        write(99,'(A76)')'    Open this file in Avogadro to visualize the molecule and the normal mode'
        write(99,'(A65)')'    Only the standard orientation and normal modes are meaningful'
        write(99,'(A83)')'    Other lines are meant to cheat Avogadro to consider this file as a Gaussian log'
        write(99,'(A29)')'------------ End ------------'
        write(99,*)
        write(99,'(A36)')'Gaussian, Inc.  All Rights Reserved.'
        write(99,'(A16)')' # freq hf/3-21g'
        write(99,'(A29)')' Charge =  0 Multiplicity = 1'
        write(99,*)
        write(99,'(A71)')'                         Standard orientation:                         '
        write(99,'(A70)')' ---------------------------------------------------------------------'
        write(99,'(A66)')' Center     Atomic      Atomic             Coordinates (Angstroms)'
        write(99,'(A67)')' Number     Number       Type             X           Y           Z'
        write(99,'(A70)')' ---------------------------------------------------------------------'
        do i=1,NAtoms
            number(i)=Symbol2Number(symbol(i))
            write(99,'(I7,I11,I12,4x,3F12.6)')i,number(i),0,structure(:,i)
        end do
        write(99,'(A70)')' ---------------------------------------------------------------------'
        do i=0,vibdim/3-1; call Print3Columns(3*i+1); end do
        select case(mod(vibdim,3))
        case(1); call PrintColumn(vibdim)
        case(2); call Print2Columns(vibdim-1)
        end select
    close(99)
    contains
    subroutine PrintColumn(start)
        integer,intent(in)::start
        integer::i
        write(99,'(I23)')start
        write(99,'(A15,F12.4)')' Frequencies --',freq(start)
        write(99,'(A15,F12.4)')' Red. masses --',0d0
        write(99,'(A15,F12.4)')' Frc consts  --',0d0
        write(99,'(A15,F12.4)')' IR Inten    --',0d0
        write(99,'(A15,F12.4)')' Raman Activ --',0d0
        write(99,'(A15,F12.4)')' Depolar (P) --',0d0
        write(99,'(A15,F12.4)')' Depolar (U) --',0d0
        write(99,'(A31)')'  Atom  AN      X      Y      Z'
        do i=1,NAtoms
            write(99,'(I6,I4,2x,3F7.2)')i,number(i),mode(:,i,start)
        end do
    end subroutine PrintColumn
    subroutine Print2Columns(start)
        integer,intent(in)::start
        integer::i
        write(99,'(2I23)')start,start+1
        write(99,'(A15,F12.4,F23.4)')' Frequencies --',freq(start),freq(start+1)
        write(99,'(A15,F12.4,F23.4)')' Red. masses --',0d0,0d0
        write(99,'(A15,F12.4,F23.4)')' Frc consts  --',0d0,0d0
        write(99,'(A15,F12.4,F23.4)')' IR Inten    --',0d0,0d0
        write(99,'(A15,F12.4,F23.4)')' Raman Activ --',0d0,0d0
        write(99,'(A15,F12.4,F23.4)')' Depolar (P) --',0d0,0d0
        write(99,'(A15,F12.4,F23.4)')' Depolar (U) --',0d0,0d0
        write(99,'(A54)')'  Atom  AN      X      Y      Z        X      Y      Z'
        do i=1,NAtoms
            write(99,'(I6,I4,2x,3F7.2,2x,3F7.2)')i,number(i),mode(:,i,start),mode(:,i,start+1)
        end do
    end subroutine Print2Columns
    subroutine Print3Columns(start)
        integer,intent(in)::start
        integer::i
        write(99,'(3I23)')start,start+1,start+2
        write(99,'(A15,F12.4,2F23.4)')' Frequencies --',freq(start),freq(start+1),freq(start+2)
        write(99,'(A15,F12.4,2F23.4)')' Red. masses --',0d0,0d0,0d0
        write(99,'(A15,F12.4,2F23.4)')' Frc consts  --',0d0,0d0,0d0
        write(99,'(A15,F12.4,2F23.4)')' IR Inten    --',0d0,0d0,0d0
        write(99,'(A15,F12.4,2F23.4)')' Raman Activ --',0d0,0d0,0d0
        write(99,'(A15,F12.4,2F23.4)')' Depolar (P) --',0d0,0d0,0d0
        write(99,'(A15,F12.4,2F23.4)')' Depolar (U) --',0d0,0d0,0d0
        write(99,'(A77)')'  Atom  AN      X      Y      Z        X      Y      Z        X      Y      Z'
        do i=1,NAtoms
            write(99,'(I6,I4,2x,3F7.2,2x,3F7.2,2x,3F7.2)')i,number(i),mode(:,i,start),mode(:,i,start+1),mode(:,i,start+2)
        end do
    end subroutine Print3Columns
end subroutine Avogadro_Vibration

!Input:  N dimensional ascendingly sorted array energy
!Output: degenerate harvests whether exists almost degenerate energy levels (energy difference < threshold)
subroutine CheckDegeneracy(degenerate, threshold, energy, N)
    logical,intent(out)::degenerate
    real*8,intent(in)::threshold
    integer,intent(in)::N
    real*8,dimension(N),intent(in)::energy
    integer::i
    degenerate=.false.
    do i=1,N-1
        if(energy(i+1)-energy(i)<threshold) then
            degenerate=.true.; exit
        end if  
    end do
end subroutine CheckDegeneracy

!At conical intersection there is a gauge degree of freedom, conical intersection adapted coordinate is
!gauge g . h = 0, where g & h are force difference & interstate coupling between intersected states
!Note this gauge does not determine adiabatic states uniquely: 8 possibilities in total
!    The transformation rotation angle can differ by arbitrary integer times of pi / 4,
!        so there are 4 possibilities (differ by pi is only a total phase change)
!    The state ordering could be exchanged, introducing 2 times of possibilities
!This can also be viewed as: we only know g & h will be along 2 lines, so there are 8 different assignments
!Reference: D. R. Yarkony, J. Chem. Phys. 112, 2111 (2000)
!Required: grad1 & grad2: energy gradient on 1st & 2nd intersected potential energy surfaces
!          h: interstate coupling between the intersected states
!          dim: integer specifying the dimension of grad1 & grad2 & h
!Optional: phi1 & phi2: wavefunction of 1st & 2nd intersected states
!          gref & href: reference g & h to uniquely determine gh orthogonalization as the 1 with smallest difference to reference out of 8
!On exit grad1, grad2, h (and optionally phi1, phi2) will be gauged
subroutine ghOrthogonalization(grad1, grad2, h, dim, phi1, phi2, gref, href)
    !Required argument:
        integer,intent(in)::dim
        real*8,dimension(dim),intent(inout)::grad1,grad2,h
    !Optional argument:
        real*8,dimension(:),intent(inout),optional::phi1,phi2
		real*8,dimension(dim),intent(in),optional::gref,href
    logical::exchange,exchangemin; integer::i
    real*8::theta,sinsqtheta,cossqtheta,sin2theta,thetamin,differencemin,difference,differencex
    real*8,dimension(dim)::g,dh11,dh21,dh22,dh11min,dh21min,dh22min,gtemp,vectemp
    real*8,allocatable,dimension(:)::phitemp
    g=(grad2-grad1)/2d0; sinsqtheta=dot_product(g,h)
    if(present(gref).and.present(href)) then
        if(dAbs(sinsqtheta)<1d-15) then
            thetamin=0d0!Try principle value
            dh21min=h
            vectemp=g-gref; difference=dot_product(vectemp,vectemp)
            vectemp=g+gref; differencex=dot_product(vectemp,vectemp)
            if(differencex<difference) then
                exchangemin=.true.; vectemp=h-href; differencemin=differencex+dot_product(vectemp,vectemp)
                dh11min=grad2; dh22min=grad1
            else
                exchangemin=.false.; vectemp=h-href; differencemin=difference+dot_product(vectemp,vectemp)
                dh11min=grad1; dh22min=grad2
            end if
        else
            theta=dot_product(g,g)-dot_product(h,h)
            if(dAbs(theta)<1d-15) then; thetamin=pid8
            else; thetamin=atan(2d0*sinsqtheta/theta)/4d0; end if
            sinsqtheta=sin(thetamin); cossqtheta=cos(thetamin)!Try principle value
            sin2theta=2d0*sinsqtheta*cossqtheta; sinsqtheta=sinsqtheta*sinsqtheta; cossqtheta=cossqtheta*cossqtheta
            dh11min=cossqtheta*grad1+sinsqtheta*grad2-sin2theta*h
            dh22min=sinsqtheta*grad1+cossqtheta*grad2+sin2theta*h
            dh21min=(cossqtheta-sinsqtheta)*h-sin2theta*g
            gtemp=(dh22min-dh11min)/2d0
            vectemp=gtemp-gref; difference=dot_product(vectemp,vectemp)
            vectemp=gtemp+gref; differencex=dot_product(vectemp,vectemp)
            if(differencex<difference) then
                exchangemin=.true.; vectemp=dh21min-href; differencemin=differencex+dot_product(vectemp,vectemp)
                vectemp=dh11min; dh11min=dh22min; dh22min=vectemp
            else
                exchangemin=.false.; vectemp=dh21min-href; differencemin=difference+dot_product(vectemp,vectemp)
            end if
        end if
        theta=thetamin; do i=1,3!Try 3 remaining solutions
            theta=theta+pid4
            sinsqtheta=sin(theta); cossqtheta=cos(theta)
            sin2theta=2d0*sinsqtheta*cossqtheta; sinsqtheta=sinsqtheta*sinsqtheta; cossqtheta=cossqtheta*cossqtheta
            dh11=cossqtheta*grad1+sinsqtheta*grad2-sin2theta*h
            dh22=sinsqtheta*grad1+cossqtheta*grad2+sin2theta*h
            dh21=(cossqtheta-sinsqtheta)*h-sin2theta*g
            gtemp=(dh22-dh11)/2d0
            vectemp=gtemp-gref; difference=dot_product(vectemp,vectemp)
            vectemp=gtemp+gref; differencex=dot_product(vectemp,vectemp)
            if(differencex<difference) then; exchange=.true.; vectemp=dh21-href; difference=differencex+dot_product(vectemp,vectemp)
            else; exchange=.false.; vectemp=dh21-href; difference=difference+dot_product(vectemp,vectemp); end if
            if(difference<differencemin) then
                differencemin=difference; exchangemin=exchange; thetamin=theta; dh21min=dh21
                if(exchange) then; dh11min=dh22; dh22min=dh11
                else; dh11min=dh11; dh22min=dh22; end if
            end if
        end do
        grad1=dh11min; grad2=dh22min; h=dh21min
        if(present(phi1).and.present(phi2)) then!Also gauge wavefunctions
            if(size(phi1)==size(phi2)) then
                sinsqtheta=sin(thetamin); cossqtheta=cos(thetamin)
                allocate(phitemp(size(phi1))); phitemp=phi1
                if(exchangemin) then
                    phi1=sinsqtheta*phitemp+cossqtheta*phi2
                    phi2=cossqtheta*phitemp-sinsqtheta*phi2
                else
                    phi1=cossqtheta*phitemp-sinsqtheta*phi2
                    phi2=sinsqtheta*phitemp+cossqtheta*phi2
                end if
                deallocate(phitemp)
            else
                write(*,'(1x,A89)')'gh orthogonolization warning: inconsistent size of wavefunctions, they will not be gauged'
            end if
        end if
    else
        if(dAbs(sinsqtheta)<1d-15) return
        theta=dot_product(g,g)-dot_product(h,h)
        if(dAbs(theta)<1d-15) then; theta=pid8
        else; theta=atan(2d0*sinsqtheta/theta)/4d0; end if
        sinsqtheta=dSin(theta); cossqtheta=dCos(theta)
        if(present(phi1).and.present(phi2)) then!Also gauge wavefunctions
            if(size(phi1)==size(phi2)) then
                allocate(phitemp(size(phi1)))
                phitemp=phi1
                phi1=cossqtheta*phitemp-sinsqtheta*phi2
                phi2=sinsqtheta*phitemp+cossqtheta*phi2
                deallocate(phitemp)
            else
                write(*,'(1x,A89)')'gh orthogonolization warning: inconsistent size of wavefunctions, they will not be gauged'
            end if
        end if
        sin2theta=2d0*sinsqtheta*cossqtheta; sinsqtheta=sinsqtheta*sinsqtheta; cossqtheta=cossqtheta*cossqtheta
        dh11=grad1
        grad1=cossqtheta*dh11+sinsqtheta*grad2-sin2theta*h
        grad2=sinsqtheta*dh11+cossqtheta*grad2+sin2theta*h
        h=(cossqtheta-sinsqtheta)*h-sin2theta*g
    end if
end subroutine ghOrthogonalization

!----------- Phase fixing -----------
    !Eigenvector has indeterminate phase, consequently any inner product involving two
    !different states also does not have determinate phase. Sometimes we need to fix it

    !Generate the permutation list of phase and phase difference according to number of states (N)
    subroutine InitializePhaseFixing(N)
        integer,intent(in)::N
        integer::i,j
        !Basis have 2^N possible phases
            Chemistry_NPhasePossibilities=ishft(1,N)-1!Unchanged case is excluded
            if(allocated(Chemistry_PhasePossibility)) deallocate(Chemistry_PhasePossibility)
            allocate(Chemistry_PhasePossibility(N,Chemistry_NPhasePossibilities))
            Chemistry_PhasePossibility(1,1)=-1d0
            Chemistry_PhasePossibility(2:N,1)=1d0
            do i=2,Chemistry_NPhasePossibilities
                Chemistry_PhasePossibility(:,i)=Chemistry_PhasePossibility(:,i-1)
                j=1
                do while(Chemistry_PhasePossibility(j,i)==-1d0)
                    Chemistry_PhasePossibility(j,i)=1d0
                    j=j+1
                end do
                Chemistry_PhasePossibility(j,i)=-1d0
            end do
        !Basis have 2^(N-1) possible phase difference, so we arbitrarily assign 1st basis has phase = 1
            Chemistry_NPhaseDifferencePossibilities=ishft(1,N-1)-1!Unchanged case is excluded
            if(allocated(Chemistry_PhaseDifferencePossibility)) deallocate(Chemistry_PhaseDifferencePossibility)
            allocate(Chemistry_PhaseDifferencePossibility(N,Chemistry_NPhaseDifferencePossibilities))
            Chemistry_PhaseDifferencePossibility(1,1)=1d0
            Chemistry_PhaseDifferencePossibility(2,1)=-1d0
            Chemistry_PhaseDifferencePossibility(3:N,1)=1d0
            do i=2,Chemistry_NPhaseDifferencePossibilities
                Chemistry_PhaseDifferencePossibility(:,i)=Chemistry_PhaseDifferencePossibility(:,i-1)
                j=2
                do while(Chemistry_PhaseDifferencePossibility(j,i)==-1d0)
                    Chemistry_PhaseDifferencePossibility(j,i)=1d0
                    j=j+1
                end do
                Chemistry_PhaseDifferencePossibility(j,i)=-1d0
            end do
    end subroutine InitializePhaseFixing

    !dim x N x N 3-order tensor dH
    !Fix off-diagonal element phase of dH by minimizing its difference from reference
    !difference harvest the minimum || dH - dH_Ref ||_F^2
    subroutine dFixdHPhase(dH, dH_Ref, difference, dim, N)
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
        do i=1,Chemistry_NPhaseDifferencePossibilities!Try out all possibilities
            temp=0d0
            do jstate=2,N!1st basis is assigned to have phase = 1
                d=Chemistry_PhaseDifferencePossibility(jstate,i)*dH(:,jstate,1)-dH_Ref(:,jstate,1)
                temp=temp+dot_product(d,d)
            end do
            do istate=2,N
                do jstate=istate+1,N
                    d=Chemistry_PhaseDifferencePossibility(istate,i)*Chemistry_PhaseDifferencePossibility(jstate,i)*dH(:,jstate,istate)-dH_Ref(:,jstate,istate)
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
                dH(:,jstate,1)=Chemistry_PhaseDifferencePossibility(jstate,indexmin)*dH(:,jstate,1)
            end forall
            forall(istate=2:N-1,jstate=3:N,jstate>istate)
                dH(:,jstate,istate)=Chemistry_PhaseDifferencePossibility(istate,indexmin)*Chemistry_PhaseDifferencePossibility(jstate,indexmin)*dH(:,jstate,istate)
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
    subroutine dFixHPhaseBydH(H, dH, dH_Ref, difference, dim, N)
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
        do i=1,Chemistry_NPhaseDifferencePossibilities!Try out all possibilities
            temp=0d0
            do jstate=2,N!1st basis is assigned to have phase = 1
                d=Chemistry_PhaseDifferencePossibility(jstate,i)*dH(:,jstate,1)-dH_Ref(:,jstate,1)
                temp=temp+dot_product(d,d)
            end do
            do istate=2,N
                do jstate=istate+1,N
                    d=Chemistry_PhaseDifferencePossibility(istate,i)*Chemistry_PhaseDifferencePossibility(jstate,i)*dH(:,jstate,istate)-dH_Ref(:,jstate,istate)
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
                 H(  jstate,1)=Chemistry_PhaseDifferencePossibility(jstate,indexmin)* H(  jstate,1)
                dH(:,jstate,1)=Chemistry_PhaseDifferencePossibility(jstate,indexmin)*dH(:,jstate,1)
            end forall
            forall(istate=2:N-1,jstate=3:N,jstate>istate)
                 H(  jstate,istate)=Chemistry_PhaseDifferencePossibility(istate,indexmin)*Chemistry_PhaseDifferencePossibility(jstate,indexmin)* H(  jstate,istate)
                dH(:,jstate,istate)=Chemistry_PhaseDifferencePossibility(istate,indexmin)*Chemistry_PhaseDifferencePossibility(jstate,indexmin)*dH(:,jstate,istate)
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
    subroutine dAssignBasisPhaseBydH(phi, dH, dH_Ref, difference, dim, N)
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
        do i=1,Chemistry_NPhaseDifferencePossibilities!Try out all possibilities
            temp=0d0
            do jstate=2,N!1st basis is assigned to have phase = 1
                d=Chemistry_PhaseDifferencePossibility(jstate,i)*dH(:,jstate,1)-dH_Ref(:,jstate,1)
                temp=temp+dot_product(d,d)
            end do
            do istate=2,N
                do jstate=istate+1,N
                    d=Chemistry_PhaseDifferencePossibility(istate,i)*Chemistry_PhaseDifferencePossibility(jstate,i)*dH(:,jstate,istate)-dH_Ref(:,jstate,istate)
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
                phi(:,jstate)=Chemistry_PhaseDifferencePossibility(jstate,indexmin)*phi(:,jstate)
                dH(:,jstate,1)=Chemistry_PhaseDifferencePossibility(jstate,indexmin)*dH(:,jstate,1)
            end forall
            forall(istate=2:N-1,jstate=3:N,jstate>istate)
                dH(:,jstate,istate)=Chemistry_PhaseDifferencePossibility(istate,indexmin)*Chemistry_PhaseDifferencePossibility(jstate,indexmin)*dH(:,jstate,istate)
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
    subroutine dFixHPhase_AssignBasisPhaseBydH(H, phi, dH, dH_Ref, difference, dim, N)
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
        do i=1,Chemistry_NPhaseDifferencePossibilities!Try out all possibilities
            temp=0d0
            do jstate=2,N!1st basis is assigned to have phase = 1
                d=Chemistry_PhaseDifferencePossibility(jstate,i)*dH(:,jstate,1)-dH_Ref(:,jstate,1)
                temp=temp+dot_product(d,d)
            end do
            do istate=2,N
                do jstate=istate+1,N
                    d=Chemistry_PhaseDifferencePossibility(istate,i)*Chemistry_PhaseDifferencePossibility(jstate,i)*dH(:,jstate,istate)-dH_Ref(:,jstate,istate)
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
                phi(:,jstate) =Chemistry_PhaseDifferencePossibility(jstate,indexmin)*phi(:,jstate)
                 H(  jstate,1)=Chemistry_PhaseDifferencePossibility(jstate,indexmin)* H(  jstate,1)
                dH(:,jstate,1)=Chemistry_PhaseDifferencePossibility(jstate,indexmin)*dH(:,jstate,1)
            end forall
            forall(istate=2:N-1,jstate=3:N,jstate>istate)
                 H(  jstate,istate)=Chemistry_PhaseDifferencePossibility(istate,indexmin)*Chemistry_PhaseDifferencePossibility(jstate,indexmin)* H(  jstate,istate)
                dH(:,jstate,istate)=Chemistry_PhaseDifferencePossibility(istate,indexmin)*Chemistry_PhaseDifferencePossibility(jstate,indexmin)*dH(:,jstate,istate)
            end forall
        end if
        difference=difference*2d0!Above only counted strictly lower triangle
        do istate=1,N!Diagonals
            d=dH(:,istate,istate)-dH_Ref(:,istate,istate)
            difference=difference+dot_product(d,d)
        end do
    end subroutine dFixHPhase_AssignBasisPhaseBydH
!--------------- End ----------------

!A, eigval, eigvec satisfy: A . eigvec(:,i) = eigval(i) * eigvec(:,i) for all 1 <= i <= N
!dA = ▽A in A representation is dim x N x N 3rd-order tensor
!Return dim x N x N 3rd-order tensor M satisfying: ▽eigvec = eigvec M (matrix multiplication on N x N)
!Note that M is anti hermitian, so only fill in strictly lower triangle
function deigvec_ByKnowneigval_dA(eigval, dA, dim, N)
    integer,intent(in)::dim,N
    real*8,dimension(N),intent(in)::eigval
    real*8,dimension(dim,N,N),intent(in)::dA
    real*8,dimension(dim,N,N)::deigvec_ByKnowneigval_dA
    integer::i,j
    forall(i=2:N,j=1:N-1,i>j)
        deigvec_ByKnowneigval_dA(:,i,j)=dA(:,i,j)/(eigval(j)-eigval(i))
    end forall
end function deigvec_ByKnowneigval_dA

!---------- Python special ----------
    subroutine py_ghOrthogonalization(grad1, grad2, h, dim)
        integer,intent(in)::dim
        real*8,dimension(dim),intent(inout)::grad1,grad2,h
        real*8::theta,sinsqtheta,cossqtheta,sin2theta
        real*8,dimension(dim)::g,dh11
        g=(grad2-grad1)/2d0; sinsqtheta=dot_product(g,h)
        if(dAbs(sinsqtheta)<1d-15) return
        theta=dot_product(g,g)-dot_product(h,h)
        if(dAbs(theta)<1d-15) then; theta=pid8
        else; theta=atan(2d0*sinsqtheta/theta)/4d0; end if
        sinsqtheta=dSin(theta); cossqtheta=dCos(theta)
        sin2theta=2d0*sinsqtheta*cossqtheta; sinsqtheta=sinsqtheta*sinsqtheta; cossqtheta=cossqtheta*cossqtheta
        dh11=grad1
        grad1=cossqtheta*dh11+sinsqtheta*grad2-sin2theta*h
        grad2=sinsqtheta*dh11+cossqtheta*grad2+sin2theta*h
        h=(cossqtheta-sinsqtheta)*h-sin2theta*g
    end subroutine py_ghOrthogonalization
!--------------- End ----------------

end module Chemistry