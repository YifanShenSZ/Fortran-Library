!Common geometry transformation applied in molecule computation:
!    Standard geometry
!    Cartesian <-> internal coodinate
!    Normal mode and vibrational frequency
module GeometryTransformation
    use General
    use LinearAlgebra
    use NonlinearOptimization
    implicit none

!Derived type
    !Example: type(InternalCoordinateDefinition),allocatable,dimension(:)::IntCDef
    !         IntCDef(iIntC) stands for the iIntC-th internal coordinate
    !         IntCDef(iIntC).NMotions is the number of motions involved in this internal coordinate
    !         IntCDef(iIntC).motion(iMotion) stands for its iMotion-th motion
    !         IntCDef(iIntC).motion(iMotion).type is its type
    !         IntCDef(iIntC).motion(iMotion).coeff is its normalized linear combination coefficient
    !         IntCDef(iIntC).motion(iMotion).atom(i) is its i-th involved atom
    type InvolvedMotion
        character*10::type!Currently only support stretching, bending, torsion
        real*8::coeff
        !For stretching, the motion coordinate is the bond length atom1_atom2
        !For bending,   the motion coordinate is bond angle atom1_atom2_atom3, range [0,pi]
        !For torsion,   the motion coordinate is dihedral angle atom1_atom2_atom3_atom4, range (-pi,pi]
        !               n_abc (the normal vector of plane abc) is a unit vector along r_ba x r_bc,
        !               Dihedral angle has same sign to n_123 x n_234 . r_23
        integer,allocatable,dimension(:)::atom
    end type InvolvedMotion
    type InternalCoordinateDefinition
        integer::NMotions
        type(InvolvedMotion),allocatable,dimension(:)::motion
    end type InternalCoordinateDefinition

!Global variable
    type(InternalCoordinateDefinition),allocatable,dimension(:)::GeometryTransformation_IntCDef!short for INTernal Coordinate DEFinition

contains
!Standardize a geometry (and optionally gradient) (optionally to a reference)
!Here we define a standard geometry as:
!    the centre of mass at origin
!    the rotational principle axes along xyz axes, with the smallest corresponding to x axis, 2nd to y, 3rd to z
!Note this definition does not determine geometry uniquely:
!    axes may take different positive direction as long as forming right-hand system (4 possibilities)
!Required argument: 3NAtoms order vector geom (geom[3i-2:3i] corresponds to [x,y,z] of i-th atom),
!                   NAtoms order vector mass (mass[i] is the mass of i-th atom), NAtoms, NStates
!Optional argument:
!    reference: uniquely define the standard geometry by the 1 with smallest difference to reference out of 4
!    difference: harvest || output geom - reference ||_2^2
!    grad: 3NAtoms order vector gradient to transform to standard coordinate
!    nadgrad: 3NAtoms x NStates x NStates 3-rd order tensor to transform to standard coordinate
subroutine StandardizeGeometry(geom,mass,NAtoms,NStates,reference,difference,grad,nadgrad)
    !Required argument
        integer,intent(in)::NAtoms,NStates
        real*8,dimension(3*NAtoms),intent(inout)::geom
        real*8,dimension(NAtoms),intent(in)::mass
    !Optional argument
        real*8,dimension(3*NAtoms),intent(in),optional::reference
        real*8,intent(out),optional::difference
        real*8,dimension(3*NAtoms),intent(inout),optional::grad
        real*8,dimension(3*NAtoms,NStates,NStates),intent(inout),optional::nadgrad
    integer::i,indicemin,istate,jstate
    real*8::SystemMass,temp
    real*8,dimension(3)::com!Short for Centre Of Mass
    real*8,dimension(3,3)::UT,moi!Short for Momentum Of Inertia
    real*8,dimension(3*NAtoms,4)::r!To store 4 legal geometries
    !Compute current centre of mass
        com=0d0
        SystemMass=0d0
        do i=1,NAtoms
            com=com+mass(i)*geom(3*i-2:3*i)
            SystemMass=SystemMass+mass(i)
        end do
        com=com/SystemMass
    !Shift centre of mass to origin, then get momentum of inertia in centre of mass frame
        moi=0d0
        do i=1,NAtoms
            geom(3*i-2:3*i)=geom(3*i-2:3*i)-com
            moi=moi+mass(i)*(dot_product(geom(3*i-2:3*i),geom(3*i-2:3*i))*UnitMatrix(3)-vector_direct_product(geom(3*i-2:3*i),geom(3*i-2:3*i),3,3))
        end do
    !Diagonalize momentum of inertia
        call My_dsyev('V',moi,com,3)
        if(triple_product(moi(:,1),moi(:,2),moi(:,3))<0d0) moi(:,1)=-moi(:,1)
        UT=transpose(moi)
    !Transform geometry (and optionally gradient) to principle axes frame
    if(present(reference)) then
        !The positive direction has 4 legal choices, determine it by comparing difference to the reference geometry
        forall(i=1:NAtoms)
            r(3*i-2:3*i,1)=matmul(UT,geom(3*i-2:3*i))
            !Flip x and y positive directions
            r(3*i-2:3*i-1,2)=-r(3*i-2:3*i-1,1)
            r(3*i,2)=r(3*i,1)
            !Flip x and z positive directions
            r(3*i-2,3)=-r(3*i-2,1)
            r(3*i,3)=-r(3*i,1)
            r(3*i-1,3)=r(3*i-1,1)
            !Flip y and z positive directions
            r(3*i-1:3*i,4)=-r(3*i-1:3*i,1)
            r(3*i-2,4)=r(3*i-2,1)
        end forall
        indicemin=1!Which one has smallest difference?
        temp=dot_product(r(:,1)-reference,r(:,1)-reference)!The smallest difference
        do i=2,4
            SystemMass=dot_product(r(:,i)-reference,r(:,i)-reference)
            if(SystemMass<temp) then
                temp=SystemMass
                indicemin=i
            end if
        end do
        if(present(difference)) difference=temp!Harvest || output geom - reference ||_2^2
        geom=r(:,indicemin)!Determine the geometry
        select case(indicemin)!Determine the principle axes
            case(2)
                moi(:,1:2)=-moi(:,1:2)
            case(3)
                moi(:,1)=-moi(:,1)
                moi(:,3)=-moi(:,3)
            case(4)
                moi(:,2:3)=-moi(:,2:3)
        end select
    else
        forall(i=1:NAtoms)
            geom(3*i-2:3*i)=matmul(UT,geom(3*i-2:3*i))
        end forall
    end if
    if(present(grad)) then
        forall(i=1:NAtoms)
            grad(3*i-2:3*i)=matmul(moi,grad(3*i-2:3*i))
        end forall
    end if
    if(present(nadgrad)) then
        forall(i=1:NAtoms,istate=1:NStates,jstate=1:NStates,istate>=jstate)
            nadgrad(3*i-2:3*i,istate,jstate)=matmul(moi,nadgrad(3*i-2:3*i,istate,jstate))
        end forall
    end if
end subroutine StandardizeGeometry

!---------- Cartesian <-> Internal ----------
    !An interal coordinate is the linear combination of several translationally and rotationally invariant displacements,
    !    but only displacements under same unit can be combined, i.e., you must treat length and angle separately,
    !    unless appropriate metric tensor is applied
    !It is OK to define more than 3NAtoms-6 (or 3NAtoms-5 for linear molecule) internal coordinates,
    !    but only 3NAtoms-6 (or 3NAtoms-5 for linear molecule) partial derivatives are independent
    !!Although the transformation from Cartesian coordinate to internal coordinate is not necessarily linear,
    !    for infinitesimal displacement it is linear, corresponding to a matrix form: dq = B dr
    !    where dq is internal coordinate differentiation, dr is Cartesian coordinate differentiation,
    !    B is Jacobian(q,r) (historically called Wilson B matrix)
    !r is a 3NAtoms order vector with r[3*i-2:3*i] corresponding to the coordinate of i-th atom
    !Nomenclature:
    !    cartdim & intdim: Cartesian & internal space dimension
    !    cartgrad & intgrad: Cartesian & internal gradient
    !    cartnadgrad & intnadgrad: Cartesian & internal nonadiabatic gradient (cartdim & intdim x NStates x NStates 3rd-order tensor)

    !This is a demo on how to define internal coordinate: read Columbus7 intcfl file
    subroutine demo_DefineInternalCoordinate(intdim)
        integer,intent(out)::intdim!Return the dimension of internal space
        character*10,allocatable,dimension(:)::MotionType
        character*24::CharTemp24
        integer::i,j,k,NLines
        integer,allocatable,dimension(:)::KLine
        real*8::DbTemp
        open(unit=99,file='intcfl',status='old')
            !Get how many lines are the definition for internal coordinates & how many internal coordinates there are
                NLines=0
                intdim=0
                read(99,*)!First line is always 'TEXAS'
                do while(.true.)
                    read(99,'(A24)')CharTemp24
                    if (index(CharTemp24,'STRE')==0.and.index(CharTemp24,'BEND')==0.and.index(CharTemp24,'TORS')==0) exit
                    NLines=NLines+1
                    if(scan(CharTemp24,'K')==1) intdim=intdim+1
                end do
                rewind 99
            !Get whether a line is the start of a new internal coordinate, and what type of motion a line stands for
                allocate(KLine(intdim+1))
                KLine(intdim+1)=NLines+1
                allocate(MotionType(NLines))
                j=1
                read(99,*)!First line is always 'TEXAS'
                do i=1,NLines
                    read(99,'(A24)')CharTemp24
                    if(scan(CharTemp24,'K')==1) then
                        KLine(j)=i
                        j=j+1
                    end if
                    if(index(CharTemp24,'STRE')>0) then
                        MotionType(i)='stretching'
                    else if(index(CharTemp24,'BEND')>0) then
                        MotionType(i)='bending'
                    else if(index(CharTemp24,'TORS')>0) then
                        MotionType(i)='torsion'
                    end if
                end do
                rewind 99
            !Finally read Columbus internal coordinate definition
                allocate(GeometryTransformation_IntCDef(intdim))
                k=1!Counter for line
                read(99,*)!First line is always 'TEXAS'
                do i=1,intdim
                    GeometryTransformation_IntCDef(i).NMotions=KLine(i+1)-KLine(i)
                    allocate(GeometryTransformation_IntCDef(i).motion(GeometryTransformation_IntCDef(i).NMotions))
                    if(GeometryTransformation_IntCDef(i).NMotions==1) then
                        GeometryTransformation_IntCDef(i).motion(1).type=MotionType(k)
                        GeometryTransformation_IntCDef(i).motion(1).coeff=1d0!Only 1 motion
                        select case(MotionType(k))
                            case('stretching')
                                allocate(GeometryTransformation_IntCDef(i).motion(1).atom(2))
                                read(99,'(A28,I5,1x,I9)')CharTemp24,GeometryTransformation_IntCDef(i).motion(1).atom
                            case('bending')
                                allocate(GeometryTransformation_IntCDef(i).motion(1).atom(3))
                                read(99,'(A28,I5,1x,I9,1x,I9)')CharTemp24,GeometryTransformation_IntCDef(i).motion(1).atom(1),&
                                    GeometryTransformation_IntCDef(i).motion(1).atom(3),GeometryTransformation_IntCDef(i).motion(1).atom(2)
                            case('torsion')
                                allocate(GeometryTransformation_IntCDef(i).motion(1).atom(4))
                                read(99,'(A28,I5,1x,I9,1x,I9,1x,I9)')CharTemp24,GeometryTransformation_IntCDef(i).motion(1).atom
                            case default!Throw a warning
                                write(*,'(1x,A51,1x,A10)')'Program abort: unsupported internal coordinate type',MotionType(k)
                                stop
                        end select
                        k=k+1
                    else
                        DbTemp=0d0
                        do j=1,GeometryTransformation_IntCDef(i).NMotions
                            GeometryTransformation_IntCDef(i).motion(j).type=MotionType(k)
                            select case(MotionType(k))
                                case('stretching')
                                    allocate(GeometryTransformation_IntCDef(i).motion(j).atom(2))
                                    read(99,'(A10,F10.7,8x,I6,1x,I9)')CharTemp24,&
                                        GeometryTransformation_IntCDef(i).motion(j).coeff,GeometryTransformation_IntCDef(i).motion(j).atom
                                case('bending')
                                    allocate(GeometryTransformation_IntCDef(i).motion(j).atom(3))
                                    read(99,'(A10,F10.7,8x,I6,1x,I9,1x,I9)')CharTemp24,GeometryTransformation_IntCDef(i).motion(j).coeff,&
                                        GeometryTransformation_IntCDef(i).motion(j).atom(1),GeometryTransformation_IntCDef(i).motion(j).atom(3),GeometryTransformation_IntCDef(i).motion(j).atom(2)
                                case('torsion')
                                    allocate(GeometryTransformation_IntCDef(i).motion(j).atom(4))
                                    read(99,'(A10,F10.7,8x,I6,1x,I9,1x,I9,1x,I9)')CharTemp24,&
                                        GeometryTransformation_IntCDef(i).motion(j).coeff,GeometryTransformation_IntCDef(i).motion(j).atom
                                case default!Throw a warning
                                    write(*,'(1x,A51,1x,A10)')'Program abort: unsupported internal coordinate type',MotionType(k)
                                    stop
                            end select
                            k=k+1
                            DbTemp=DbTemp+GeometryTransformation_IntCDef(i).motion(j).coeff*GeometryTransformation_IntCDef(i).motion(j).coeff
                        end do
                        DbTemp=Sqrt(DbTemp)
                        forall(j=1:GeometryTransformation_IntCDef(i).NMotions)
                            GeometryTransformation_IntCDef(i).motion(j).coeff=GeometryTransformation_IntCDef(i).motion(j).coeff/DbTemp
                        end forall
                    end if
                end do
        close(99)
        !Clean up
            deallocate(MotionType)
            deallocate(KLine)
    end subroutine demo_DefineInternalCoordinate

    !========== Cartesian -> Internal ==========
        !Transform geometry (and optionally gradient) from Cartesian coordinate to internal coordinate
        !Required: r, cartdim, q, intdim, NStates
        !Optional: cartgrad, intgrad, cartnadgrad, intnadgrad
        !r & cartgrad & cartnadgrad are the input Cartesian space value, q & intgrad & intnadgrad harvest corresponding internal space value
        subroutine Cartesian2Internal(r,cartdim,q,intdim,NStates,cartgrad,intgrad,cartnadgrad,intnadgrad)
            !Required argument
                integer,intent(in)::cartdim,intdim,NStates
                real*8,dimension(cartdim),intent(in)::r
                real*8,dimension(intdim),intent(out)::q
            !Optional argument
                real*8,dimension(cartdim),intent(in),optional::cartgrad
                real*8,dimension(intdim),intent(out),optional::intgrad
                real*8,dimension(cartdim,NStates,NStates),intent(in),optional::cartnadgrad
                real*8,dimension(intdim,NStates,NStates),intent(out),optional::intnadgrad
            integer::i,j
            real*8,dimension(intdim,cartdim)::B
            call WilsonBMatrixAndInternalCoordinateq(B,q,r,intdim,cartdim)
            if(present(cartgrad).and.present(intgrad)) then
                call dGeneralizedInverseTranspose(B,intdim,cartdim)
                intgrad=matmul(B,cartgrad)
                if(present(cartnadgrad).and.present(intnadgrad)) then
                    forall(i=1:NStates,j=1:NStates)
                        intnadgrad(:,i,j)=matmul(B,cartnadgrad(:,i,j))
                    end forall
                end if
            else if(present(cartnadgrad).and.present(intnadgrad)) then
                call dGeneralizedInverseTranspose(B,intdim,cartdim)
                forall(i=1:NStates,j=1:NStates)
                    intnadgrad(:,i,j)=matmul(B,cartnadgrad(:,i,j))
                end forall
            end if
        end subroutine Cartesian2Internal
        
        !Generate Wilson B matrix and internal coordinate q from Cartesian coordinate r and internal coordinate definition
        !The definition is a global variable, so only take r as input
        subroutine WilsonBMatrixAndInternalCoordinateq(B,q,r,intdim,cartdim)
            integer,intent(in)::intdim,cartdim
            real*8,dimension(intdim,cartdim),intent(out)::B
            real*8,dimension(intdim),intent(out)::q
            real*8,dimension(cartdim),intent(in)::r
            integer::iIntC,iMotion
            real*8::qMotion
            real*8,dimension(cartdim)::BRowVector
            B=0d0
            q=0d0
            do iIntC=1,intdim
                do iMotion=1,GeometryTransformation_IntCDef(iIntC).NMotions
                    select case(GeometryTransformation_IntCDef(iIntC).motion(iMotion).type)
                        case('stretching')
                            call bAndStretching(BRowVector,qMotion,r,GeometryTransformation_IntCDef(iIntC).motion(iMotion),cartdim)
                        case('bending')
                            call bAndBending(BRowVector,qMotion,r,GeometryTransformation_IntCDef(iIntC).motion(iMotion),cartdim)
                        case('torsion')
                            call bAndTorsion(BRowVector,qMotion,r,GeometryTransformation_IntCDef(iIntC).motion(iMotion),cartdim)
                        case default!Throw a warning
                            write(*,'(1x,A51,1x,A10)')'Program abort: unsupported internal coordinate type',GeometryTransformation_IntCDef(iIntC).motion(iMotion).type
                            stop
                    end select
                    B(iIntC,:)=B(iIntC,:)+GeometryTransformation_IntCDef(iIntC).motion(iMotion).coeff*BRowVector
                    q(iIntC)=q(iIntC)+GeometryTransformation_IntCDef(iIntC).motion(iMotion).coeff*qMotion
                end do
            end do
            contains
                !Generate the transformation vector b from dr to dq: b . dr = dq
                !Transform from Cartesian coordinate r to a certain motion coordinate q
                !Internal coordinate is the linear combination of several motions,
                !so b contributes (but not necessarily equals) to one row of Wilson B matrix
                ! ( i-th row vector of B ) . dr =  d( i-th internal coordinate )
                !For stretching, q = bond length
                subroutine bAndStretching(b,q,r,motion,cartdim)
                    integer,intent(in)::cartdim
                    real*8,dimension(cartdim),intent(out)::b
                    real*8,intent(out)::q
                    real*8,dimension(cartdim),intent(in)::r
                    type(InvolvedMotion),intent(in)::motion
                    real*8,dimension(3)::runit12
                    b=0d0
                    runit12=r(3*motion.atom(2)-2:3*motion.atom(2))-r(3*motion.atom(1)-2:3*motion.atom(1))
                    q=Norm2(runit12)
                    runit12=runit12/q
                    b(3*motion.atom(1)-2:3*motion.atom(1))=-runit12
                    b(3*motion.atom(2)-2:3*motion.atom(2))=runit12
                end subroutine bAndStretching
                !For bending, q = bond angle
                subroutine bAndBending(b,q,r,motion,cartdim)
                    integer,intent(in)::cartdim
                    real*8,dimension(cartdim),intent(out)::b
                    real*8,intent(out)::q
                    real*8,dimension(cartdim),intent(in)::r
                    type(InvolvedMotion),intent(in)::motion
                    real*8::r21,r23,costheta,sintheta
                    real*8,dimension(3)::runit21,runit23
                    b=0d0
                    runit21=r(3*motion.atom(1)-2:3*motion.atom(1))-r(3*motion.atom(2)-2:3*motion.atom(2))
                        r21=Norm2(runit21)
                        runit21=runit21/r21
                    runit23=r(3*motion.atom(3)-2:3*motion.atom(3))-r(3*motion.atom(2)-2:3*motion.atom(2))
                        r23=Norm2(runit23)
                        runit23=runit23/r23
                    costheta=dot_product(runit21,runit23)
                    sintheta=dSqrt(1d0-costheta*costheta)
                    b(3*motion.atom(1)-2:3*motion.atom(1))=(costheta*runit21-runit23)/(sintheta*r21)
                    b(3*motion.atom(3)-2:3*motion.atom(3))=(costheta*runit23-runit21)/(sintheta*r23)
                    b(3*motion.atom(2)-2:3*motion.atom(2))=-b(3*motion.atom(1)-2:3*motion.atom(1))-b(3*motion.atom(3)-2:3*motion.atom(3))
                    q=acos(costheta)
                end subroutine bAndBending
                !For torsion, q = dihedral angle
                subroutine bAndTorsion(b,q,r,motion,cartdim)
                    integer,intent(in)::cartdim
                    real*8,dimension(cartdim),intent(out)::b
                    real*8,intent(out)::q
                    real*8,dimension(cartdim),intent(in)::r
                    type(InvolvedMotion),intent(in)::motion
                    real*8::r21,r23,r43,costheta1,sintheta1,costheta2,sintheta2
                    real*8,dimension(3)::runit23,n123,n234
                    b=0d0!Initialize
                    !Prepare
                    n123=r(3*motion.atom(1)-2:3*motion.atom(1))-r(3*motion.atom(2)-2:3*motion.atom(2))
                        r21=Norm2(n123)
                        n123=n123/r21
                    runit23=r(3*motion.atom(3)-2:3*motion.atom(3))-r(3*motion.atom(2)-2:3*motion.atom(2))
                        r23=Norm2(runit23)
                        runit23=runit23/r23
                    n234=r(3*motion.atom(3)-2:3*motion.atom(3))-r(3*motion.atom(4)-2:3*motion.atom(4))
                        r43=Norm2(n234)
                        n234=n234/r43
                    costheta1=dot_product(n123,runit23)
                    sintheta1=dSqrt(1d0-costheta1*costheta1)
                    n123=cross_product(n123,runit23)
                        n123=n123/sintheta1
                    costheta2=dot_product(runit23,n234)
                    sintheta2=dSqrt(1d0-costheta2*costheta2)
                    n234=cross_product(runit23,n234)
                        n234=n234/sintheta2
                    !Output
                    b(3*motion.atom(1)-2:3*motion.atom(1))=n123/(r21*sintheta1)
                    b(3*motion.atom(4)-2:3*motion.atom(4))=-n234/(r43*sintheta2)
                    b(3*motion.atom(2)-2:3*motion.atom(2))=(r21*costheta1-r23)/(r21*r23*sintheta1)*n123+costheta2/(r23*sintheta2)*n234
                    b(3*motion.atom(3)-2:3*motion.atom(3))=(r23-r43*costheta2)/(r23*r43*sintheta2)*n234-costheta1/(r23*sintheta1)*n123
                    q=acos(dot_product(n123,n234))
                    if(triple_product(n123,n234,runit23)<0d0) q=-q
                end subroutine bAndTorsion
        end subroutine WilsonBMatrixAndInternalCoordinateq

        !If you want internal coordinate only, you may adopt convenient functions below

        !Generate internal coordinate q from Cartesian coordinate r and internal coordinate definition
        !The definition is a global variable, so only take r as input
        function InternalCoordinateq(r,intdim,cartdim)
            integer,intent(in)::intdim,cartdim
            real*8,dimension(cartdim),intent(in)::r
            real*8,dimension(intdim)::InternalCoordinateq
            integer::iIntC,iMotion
            InternalCoordinateq=0d0
            do iIntC=1,intdim
                do iMotion=1,GeometryTransformation_IntCDef(iIntC).NMotions
                    select case(GeometryTransformation_IntCDef(iIntC).motion(iMotion).type)
                        case('stretching')
                            InternalCoordinateq(iIntC)=InternalCoordinateq(iIntC)+&
                                GeometryTransformation_IntCDef(iIntC).motion(iMotion).coeff*stretching(r,GeometryTransformation_IntCDef(iIntC).motion(iMotion).atom,cartdim)
                        case('bending')
                            InternalCoordinateq(iIntC)=InternalCoordinateq(iIntC)+&
                                GeometryTransformation_IntCDef(iIntC).motion(iMotion).coeff*bending(r,GeometryTransformation_IntCDef(iIntC).motion(iMotion).atom,cartdim)
                        case('torsion')
                            InternalCoordinateq(iIntC)=InternalCoordinateq(iIntC)+&
                                GeometryTransformation_IntCDef(iIntC).motion(iMotion).coeff*torsion(r,GeometryTransformation_IntCDef(iIntC).motion(iMotion).atom,cartdim)
                        case default!Throw a warning
                            write(*,'(1x,A51,1x,A10)')'Program abort: unsupported internal coordinate type',GeometryTransformation_IntCDef(iIntC).motion(iMotion).type
                            stop
                    end select
                end do
            end do
        end function InternalCoordinateq
        
        !Transform from Cartesian coordinate r to a certain motion coordinate q, atomlist defines which atoms are involved
        !For stretching, q = bond length atom1_atom2
        real*8 function stretching(r,atomlist,cartdim)
            integer,intent(in)::cartdim
            real*8,dimension(cartdim),intent(in)::r
            integer,dimension(2),intent(in)::atomlist
            real*8,dimension(3)::r12
            r12=r(3*atomlist(2)-2:3*atomlist(2))-r(3*atomlist(1)-2:3*atomlist(1))
            stretching=Norm2(r12)
        end function stretching
        !For bending, q = bond angle atom1_atom2_atom3, range [0,pi]
        real*8 function bending(r,atomlist,cartdim)
            integer,intent(in)::cartdim
            real*8,dimension(cartdim),intent(in)::r
            integer,dimension(3),intent(in)::atomlist
            real*8,dimension(3)::runit21,runit23
            runit21=r(3*atomlist(1)-2:3*atomlist(1))-r(3*atomlist(2)-2:3*atomlist(2))
                runit21=runit21/Norm2(runit21)
            runit23=r(3*atomlist(3)-2:3*atomlist(3))-r(3*atomlist(2)-2:3*atomlist(2))
                runit23=runit23/Norm2(runit23)
            bending=acos(dot_product(runit21,runit23))
        end function bending
        !For torsion, q = dihedral angle atom1_atom2_atom3_atom4, range (-pi,pi]
        !    n_abc (the normal vector of plane abc) is a unit vector along r_ba x r_bc,
        !    Dihedral angle has same sign to n_123 x n_234 . r_23
        real*8 function torsion(r,atomlist,cartdim)
            integer,intent(in)::cartdim
            real*8,dimension(cartdim),intent(in)::r
            integer,dimension(4),intent(in)::atomlist
            real*8,dimension(3)::r21,r23,r43
            r21=r(3*atomlist(1)-2:3*atomlist(1))-r(3*atomlist(2)-2:3*atomlist(2))
            r23=r(3*atomlist(3)-2:3*atomlist(3))-r(3*atomlist(2)-2:3*atomlist(2))
            r43=r(3*atomlist(3)-2:3*atomlist(3))-r(3*atomlist(4)-2:3*atomlist(4))
            r21=cross_product(r21,r23)!r21 = n123, temporarily
                r21=r21/Norm2(r21)
            r43=cross_product(r23,r43)!r43 = n234, temporarily
                r43=r43/Norm2(r43)
            torsion=acos(dot_product(r21,r43))
            if(triple_product(r21,r43,r23)<0d0) torsion=-torsion
        end function torsion
    !=================== End ===================

    !========== Cartesian <- Internal ==========
        !Generate Cartesian coordinate r from internal coordinate q and internal coordinate definition
        !The definition is a global variable, so only take r as input
        !Optional argument:
        !    mass: if present, will standardize r (otherwise r varies by arbitrary translation & rotation)
        !    r0: initial guess of r; if mass is also present, will standardize with r0 as reference
        function CartesianCoordinater(q,cartdim,intdim,mass,r0)
            !Required argument
                integer,intent(in)::cartdim,intdim
                real*8,dimension(intdim)::q
            !Optional argument
                real*8,dimension(cartdim/3),intent(in),optional::mass
                real*8,dimension(cartdim),intent(in),optional::r0
            real*8,dimension(cartdim)::CartesianCoordinater
            if(present(r0)) then!Initial guess
                CartesianCoordinater=r0
            else
                call random_number(CartesianCoordinater)
            end if
            call TrustRegion(Residue,CartesianCoordinater,intdim,cartdim)
            if(present(mass)) then
                if(present(r0)) then
                    call StandardizeGeometry(CartesianCoordinater,mass,cartdim/3,1,reference=r0)
                else
                    call StandardizeGeometry(CartesianCoordinater,mass,cartdim/3,1)
                end if
            end if
            contains
                subroutine Residue(res,r,intdim,cartdim)
                    integer,intent(in)::intdim,cartdim
                    real*8,dimension(intdim),intent(out)::res
                    real*8,dimension(cartdim),intent(in)::r
                    res=InternalCoordinateq(r,intdim,cartdim)-q
                end subroutine Residue
        end function CartesianCoordinater

        !Transform geometry (and optionally gradient) from Cartesian coordinate to internal coordinate
        !Required: q, intdim, r, cartdim, NStates
        !Optional: mass & r0 (same as above), intgrad, cartgrad, intnadgrad, cartnadgrad
        !q & intgrad & intnadgrad are the input internal space value, r & cartgrad & cartnadgrad harvest corresponding Cartesian space value
        subroutine Internal2Cartesian(q,intdim,r,cartdim,NStates,mass,r0,intgrad,cartgrad,intnadgrad,cartnadgrad)
            !Required argument
                integer,intent(in)::intdim,cartdim,NStates
                real*8,dimension(intdim),intent(in)::q
                real*8,dimension(cartdim),intent(out)::r
            !Optional argument
                real*8,dimension(cartdim/3),intent(in),optional::mass
                real*8,dimension(cartdim),intent(in),optional::r0
                real*8,dimension(intdim),intent(in),optional::intgrad
                real*8,dimension(cartdim),intent(out),optional::cartgrad
                real*8,dimension(intdim,NStates,NStates),intent(in),optional::intnadgrad
                real*8,dimension(cartdim,NStates,NStates),intent(out),optional::cartnadgrad
            integer::i,j
            real*8,dimension(intdim)::qtemp
            real*8,dimension(intdim,cartdim)::B
            if(present(mass)) then
                if(present(r0)) then
                    r=CartesianCoordinater(q,cartdim,intdim,mass,r0)
                else
                    r=CartesianCoordinater(q,cartdim,intdim,mass)
                end if
            else
                if(present(r0)) then
                    r=CartesianCoordinater(q,cartdim,intdim,r0)
                else
                    r=CartesianCoordinater(q,cartdim,intdim)
                end if
            end if
            call WilsonBMatrixAndInternalCoordinateq(B,qtemp,r,intdim,cartdim)
            if(present(intgrad).and.present(cartgrad)) cartgrad=matmul(transpose(B),intgrad)
            if(present(intnadgrad).and.present(cartnadgrad)) then
                forall(i=1:NStates,j=1:NStates)
                    cartnadgrad(:,i,j)=matmul(transpose(B),intnadgrad(:,i,j))
                end forall
            end if
        end subroutine Internal2Cartesian
    !=================== End ===================
!------------------- End --------------------

!--------------- Normal mode ----------------
	!Normal mode is the eigenvectors of Hessian:
	!    In cartesian coordinate, it is the usual eigenvector
	!    In  internal coordinate, it is the generalized eigenvector of G under Hessian metric
	!G is built from mass and Wilson B matrix, for details see Wilson GF method in: (note that Wilson calls Hessian by F)
	!E. B. Wilson, J. C. Decius, P. C. Cross, Molecular viobrations: the theory of infrared and Raman vibrational spectra (Dover, 1980)

	!Input:  3NAtoms order real symmetric matrix H = Hessian in Cartesian coordinate (will be overwritten)
	!Output: vibdim order vector freq = vibrational angular frequencies (negative if imaginary)
	!        vibdim order matrix mode = normal modes in input frame contained in each column
	!Lowest 3 * NAtoms - vibdim modes are considered translation and rotation thus ruled out
	subroutine VibrationAnalysis(freq,mode,vibdim,H,mass,NAtoms)
		integer,intent(in)::vibdim,NAtoms
		real*8,dimension(vibdim),intent(out)::freq
		real*8,dimension(vibdim,vibdim),intent(out)::mode
        real*8,dimension(3*NAtoms,3*NAtoms),intent(inout)::H
		real*8,dimension(NAtoms),intent(in)::mass
		integer::i
		integer,dimension(3*NAtoms)::indices
		real*8,dimension(NAtoms)::sqrtmass
		real*8,dimension(3*NAtoms)::freqall,freqabs
		!Obtain freq^2 and normal modes
		    sqrtmass=dSqrt(mass)
		    forall(i=1:NAtoms)
		    	H(:,3*i-2:3*i)=H(:,3*i-2:3*i)/sqrtmass(i)
		    	H(3*i-2:3*i,:)=H(3*i-2:3*i,:)/sqrtmass(i)
		    end forall
			call My_dsyev('V',H,freqall,3*NAtoms)
		!Rule out 3 * NAtoms - vibdim translations and rotations
		    freqabs=dAbs(freqall)
		    forall(i=1:3*NAtoms)
		        indices(i)=i
		    end forall
			call dQuickSort(freqabs,1,3*NAtoms,indices,3*NAtoms)
			freqabs=freqall
		    forall(i=1:vibdim)
				freqall(i)=freqabs(indices(i+3*NAtoms-vibdim))
				mode(:,i)=H(:,indices(i+3*NAtoms-vibdim))
            end forall
		do i=1,vibdim!freq^2 -> freq
			if(freqall(i)<0d0) then
                freq(i)=-dSqrt(-freqall(i))
            else
                freq(i)=dSqrt(freqall(i))
            end if
		end do
		!Sort freq ascendingly, then sort normal modes accordingly
		    forall(i=1:vibdim)
		        indices(i)=i
			end forall
			call dQuickSort(freq,1,vibdim,indices(1:vibdim),vibdim)
			H(:,1:vibdim)=mode
			forall(i=1:vibdim)
                mode(:,i)=H(:,indices(i))
            end forall
	end subroutine VibrationAnalysis

    !Use Wilson GF method to obtain normal mode and vibrational frequency from Hessian in internal coordinate
    !Input:  intdim order real symmetric matrix H = Hessian in internal coordinate
    !             intdim x 3*NAtoms matrix B      = Wilson B matrix
    !              NAtoms order array mass        = mass of each atom
    !Output: freq = vibrational angular frequencies (negative if imaginary)
    !         H   = normal modes in input frame
    subroutine WilsonGFMethod(freq,H,intdim,B,mass,NAtoms)
        integer,intent(in)::intdim,NAtoms
        real*8,dimension(intdim),intent(out)::freq
        real*8,dimension(intdim,intdim),intent(inout)::H
        real*8,dimension(intdim,3*NAtoms),intent(inout)::B
        real*8,dimension(NAtoms),intent(in)::mass
        integer::i
        integer,dimension(intdim)::indices
        real*8,dimension(intdim)::freqtemp
        real*8,dimension(intdim,intdim)::GF
        real*8,dimension(intdim,3*NAtoms)::Btemp
        !GF method: obtain freq^2 and normal modes
            forall(i=1:NAtoms)
                Btemp(:,3*i-2:3*i)=B(:,3*i-2:3*i)/mass(i)
            end forall
            call syL2U(H,intdim)
            GF=matmul(matmul(Btemp,transpose(B)),H)
            call My_dgeev('V',GF,freq,freqtemp,H,intdim)!Hessian is not necessarily positive definite, so do not call sygv
        do i=1,intdim!freq^2 -> freq
            if(freq(i)<0d0) then
                freq(i)=-dSqrt(-freq(i))
            else
                freq(i)=dSqrt(freq(i))
            end if
        end do
        !Sort freq ascendingly, then sort normal modes accordingly
            forall(i=1:intdim)
                indices(i)=i
            end forall
            call dQuickSort(freq,1,intdim,indices,intdim)
            GF=H
            forall(i=1:intdim)
                H(:,i)=GF(:,indices(i))
            end forall
	end subroutine WilsonGFMethod
!------------------- End --------------------

end module GeometryTransformation