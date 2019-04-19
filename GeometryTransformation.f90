!Common geometry transformation applied in molecule computation
module GeometryTransformation
    use General
    use LinearAlgebra
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
!----------- Standardize geometry -----------
    !Here we define a standard geometry as:
    !    the centre of mass is at origin
    !    the rotational principle axes are along xyz axes, with the smallest corresponding to x axis, 2nd to y, 3rd to z
    !Note this definition does not determine geometry uniquely:
    !    axes may take different positive direction as long as forming right-hand system (4 possibilities)

    !Standardize a geometry
    subroutine StandardizeGeometry(geom,mass,NAtoms)
        integer,intent(in)::NAtoms
        real*8,dimension(3*NAtoms),intent(inout)::geom
        real*8,dimension(NAtoms),intent(in)::mass
        integer::i
        real*8::SystemMass
        real*8,dimension(3)::com!Short for Centre Of Mass
        real*8,dimension(3,3)::moi!Short for Momentum Of Inertia
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
        moi=transpose(moi)
        !Transform geometry and gradients to principle axes frame
        forall(i=1:NAtoms)
            geom(3*i-2:3*i)=matmul(moi,geom(3*i-2:3*i))
        end forall
    end subroutine StandardizeGeometry
    !Also transform gradient
    subroutine StandardizeGeometry_Gradient(geom,grad,mass,NAtoms)
        integer,intent(in)::NAtoms
        real*8,dimension(3*NAtoms),intent(inout)::geom,grad
        real*8,dimension(NAtoms),intent(in)::mass
        integer::i
        real*8::SystemMass
        real*8,dimension(3)::com!Short for Centre Of Mass
        real*8,dimension(3,3)::UT,moi!Short for Momentum Of Inertia
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
        !Transform geometry and gradients to principle axes frame
        forall(i=1:NAtoms)
            geom(3*i-2:3*i)=matmul(UT,geom(3*i-2:3*i))
            grad(3*i-2:3*i)=matmul(moi,grad(3*i-2:3*i))
        end forall
    end subroutine StandardizeGeometry_Gradient
    !In nonadiabatic process, we need the gradient of NStates order matrix
    subroutine StandardizeGeometry_NonadabaticGradient(geom,grad,mass,NAtoms,NStates)
        integer,intent(in)::NAtoms,NStates
        real*8,dimension(3*NAtoms),intent(inout)::geom
        real*8,dimension(3*NAtoms,NStates,NStates),intent(inout)::grad
        real*8,dimension(NAtoms),intent(in)::mass
        integer::i,istate,jstate
        real*8::SystemMass
        real*8,dimension(3)::com!Short for Centre Of Mass
        real*8,dimension(3,3)::UT,moi!Short for Momentum Of Inertia
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
        !Transform geometry and gradients to principle axes frame
        forall(i=1:NAtoms)
            geom(3*i-2:3*i)=matmul(UT,geom(3*i-2:3*i))
            forall(istate=1:NStates,jstate=1:NStates,istate>=jstate)
                grad(3*i-2:3*i,istate,jstate)=matmul(moi,grad(3*i-2:3*i,istate,jstate))
            end forall
        end forall
    end subroutine StandardizeGeometry_NonadabaticGradient

    !StandardizeGeometry a point to a reference: choose the 1 with smallest difference to reference out of 4
    !difference = 2 norm square of difference
    subroutine StandardizeGeometry2Reference(geom,reference,difference,mass,NAtoms)
        integer,intent(in)::NAtoms
        real*8,dimension(3*NAtoms),intent(inout)::geom
        real*8,dimension(3*NAtoms),intent(in)::reference
        real*8,intent(out)::difference
        real*8,dimension(NAtoms),intent(in)::mass
        integer::i,indicemin
        real*8::SystemMass
        real*8,dimension(3)::com!Short for Centre Of Mass
        real*8,dimension(3,3)::moi!Short for Momentum Of Inertia
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
        moi=transpose(moi)
        !The positive direction has 4 legal choices, determine it by comparing difference to the reference geometry
        do i=1,NAtoms
            r(3*i-2:3*i,1)=matmul(moi,geom(3*i-2:3*i))
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
        end do
        indicemin=1!Which one has smallest difference?
        difference=dot_product(r(:,1)-reference,r(:,1)-reference)!The smallest difference
        do i=2,4
            SystemMass=dot_product(r(:,i)-reference,r(:,i)-reference)
            if(SystemMass<difference) then
                difference=SystemMass
                indicemin=i
            end if
        end do
        geom=r(:,indicemin)!Transform geometry to principle axes frame
    end subroutine StandardizeGeometry2Reference
    !Also transform gradient
    subroutine StandardizeGeometry2Reference_Gradient(geom,grad,reference,difference,mass,NAtoms)
        integer,intent(in)::NAtoms
        real*8,dimension(3*NAtoms),intent(inout)::geom,grad
        real*8,dimension(3*NAtoms),intent(in)::reference
        real*8,intent(out)::difference
        real*8,dimension(NAtoms),intent(in)::mass
        integer::i,indicemin
        real*8::SystemMass
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
        !The positive direction has 4 legal choices, determine it by comparing difference to the reference geometry
        do i=1,NAtoms
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
        end do
        indicemin=1!Which one has smallest difference?
        difference=dot_product(r(:,1)-reference,r(:,1)-reference)!The smallest difference
        do i=2,4
            SystemMass=dot_product(r(:,i)-reference,r(:,i)-reference)
            if(SystemMass<difference) then
                difference=SystemMass
                indicemin=i
            end if
        end do
        geom=r(:,indicemin)!Transform geometry to principle axes frame
        !Transform gradients to principle axes frame
        select case(indicemin)
            case(2)
                moi(:,1:2)=-moi(:,1:2)
            case(3)
                moi(:,1)=-moi(:,1)
                moi(:,3)=-moi(:,3)
            case(4)
                moi(:,2:3)=-moi(:,2:3)
            case default
        end select
        forall(i=1:NAtoms)
            grad(3*i-2:3*i)=matmul(moi,grad(3*i-2:3*i))
        end forall
    end subroutine StandardizeGeometry2Reference_Gradient
    !In nonadiabatic process, we need the gradient of NStates order matrix
    subroutine StandardizeGeometry2Reference_NonadiabaticGradient(geom,grad,reference,difference,mass,NAtoms,NStates)
        integer,intent(in)::NAtoms,NStates
        real*8,dimension(3*NAtoms),intent(inout)::geom
        real*8,dimension(3*NAtoms,NStates,NStates),intent(inout)::grad
        real*8,dimension(3*NAtoms),intent(in)::reference
        real*8,intent(out)::difference
        real*8,dimension(NAtoms),intent(in)::mass
        integer::i,istate,jstate
        real*8::SystemMass
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
        !The positive direction has 4 legal choices, determine it by comparing difference to the reference geometry
        do i=1,NAtoms
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
        end do
        istate=1!Which one has smallest difference?
        difference=dot_product(r(:,1)-reference,r(:,1)-reference)!The smallest difference
        do i=2,4
            SystemMass=dot_product(r(:,i)-reference,r(:,i)-reference)
            if(SystemMass<difference) then
                difference=SystemMass
                istate=i
            end if
        end do
        geom=r(:,istate)!Transform geometry to principle axes frame
        !Transform gradients to principle axes frame
        select case(istate)
            case(2)
                moi(:,1:2)=-moi(:,1:2)
            case(3)
                moi(:,1)=-moi(:,1)
                moi(:,3)=-moi(:,3)
            case(4)
                moi(:,2:3)=-moi(:,2:3)
            case default
        end select
        forall(i=1:NAtoms,istate=1:NStates,jstate=1:NStates,istate>=jstate)
            grad(3*i-2:3*i,istate,jstate)=matmul(moi,grad(3*i-2:3*i,istate,jstate))
        end forall
    end subroutine StandardizeGeometry2Reference_NonadiabaticGradient
!------------------- End --------------------

!---------- Cartesian <-> Internal ----------
    !An interal coordinate is the linear combination of several translationally and rotationally invariant displacements,
    !    but only displacements under same unit can be combined, i.e., you must treat length and angle separately,
    !    unless appropriate metric tensor is applied
    !It is OK to define more than 3N-6 (or 3N-5 for linear molecule) internal coordinates,
    !    but only 3N-6 (or 3N-5 for linear molecule) partial derivatives are independent
    !!Although the transformation from Cartesian coordinate to internal coordinate is not necessarily linear,
    !    for infinitesimal displacement it is linear, corresponding to a matrix form: dq = B dr
    !    where dq is internal coordinate differentiation, dr is Cartesian coordinate differentiation,
    !    B is Jacobian(q,r) (historically called Wilson B matrix)

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
        !-------------- Transform geometry only --------------
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
                                    GeometryTransformation_IntCDef(iIntC).motion(iMotion).coeff*stretching(r,GeometryTransformation_IntCDef(iIntC).motion(iMotion),cartdim)
                            case('bending')
                                InternalCoordinateq(iIntC)=InternalCoordinateq(iIntC)+&
                                    GeometryTransformation_IntCDef(iIntC).motion(iMotion).coeff*bending(r,GeometryTransformation_IntCDef(iIntC).motion(iMotion),cartdim)
                            case('torsion')
                                InternalCoordinateq(iIntC)=InternalCoordinateq(iIntC)+&
                                    GeometryTransformation_IntCDef(iIntC).motion(iMotion).coeff*torsion(r,GeometryTransformation_IntCDef(iIntC).motion(iMotion),cartdim)
                            case default!Throw a warning
                                write(*,'(1x,A51,1x,A10)')'Program abort: unsupported internal coordinate type',GeometryTransformation_IntCDef(iIntC).motion(iMotion).type
                                stop
                        end select
                    end do
                end do
            end function InternalCoordinateq
            
            !Transform from Cartesian coordinate r to a certain motion coordinate q
            !For stretching, q = bond length
            real*8 function stretching(r,motion,cartdim)
                integer,intent(in)::cartdim
                real*8,dimension(cartdim),intent(in)::r
                type(InvolvedMotion),intent(in)::motion
                real*8,dimension(3)::r12
                r12=r(3*motion.atom(2)-2:3*motion.atom(2))-r(3*motion.atom(1)-2:3*motion.atom(1))
                stretching=Norm2(r12)
            end function stretching
            !For bending, q = bond angle
            real*8 function bending(r,motion,cartdim)
                integer,intent(in)::cartdim
                real*8,dimension(cartdim),intent(in)::r
                type(InvolvedMotion),intent(in)::motion
                real*8,dimension(3)::runit21,runit23
                runit21=r(3*motion.atom(1)-2:3*motion.atom(1))-r(3*motion.atom(2)-2:3*motion.atom(2))
                    runit21=runit21/Norm2(runit21)
                runit23=r(3*motion.atom(3)-2:3*motion.atom(3))-r(3*motion.atom(2)-2:3*motion.atom(2))
                    runit23=runit23/Norm2(runit23)
                bending=acos(dot_product(runit21,runit23))
            end function bending
            !For torsion, q = dihedral angle
            real*8 function torsion(r,motion,cartdim)
                integer,intent(in)::cartdim
                real*8,dimension(cartdim),intent(in)::r
                type(InvolvedMotion),intent(in)::motion
                real*8,dimension(3)::r21,r23,r43
                r21=r(3*motion.atom(1)-2:3*motion.atom(1))-r(3*motion.atom(2)-2:3*motion.atom(2))
                r23=r(3*motion.atom(3)-2:3*motion.atom(3))-r(3*motion.atom(2)-2:3*motion.atom(2))
                r43=r(3*motion.atom(3)-2:3*motion.atom(3))-r(3*motion.atom(4)-2:3*motion.atom(4))
                r21=cross_product(r21,r23)!r21 = n123, temporarily
                    r21=r21/Norm2(r21)
                r43=cross_product(r23,r43)!r43 = n234, temporarily
                    r43=r43/Norm2(r43)
                torsion=acos(dot_product(r21,r43))
                if(triple_product(r21,r43,r23)<0d0) torsion=-torsion
            end function torsion
        !------------------------ End ------------------------
            
        !---------- Transform geometry and gradient ----------
            !Transform geometry and gradient from Cartesian coordinate to internal coordinate
            subroutine Cartesian2Internal(r,cartgrad,cartdim,q,intgrad,intdim)
                integer,intent(in)::cartdim,intdim
                real*8,dimension(cartdim),intent(in)::r,cartgrad
                real*8,dimension(intdim),intent(out)::q,intgrad
                real*8,dimension(intdim,cartdim)::B
                call WilsonBMatrixAndInternalCoordinateq(B,q,r,intdim,cartdim)
                call Cartesian2InternalGradient(cartgrad,cartdim,intgrad,intdim,B)
            end subroutine Cartesian2Internal
            !In nonadiabatic process, we need the gradient of NStates order matrix
            subroutine Cartesian2Internal_Nonadiabatic(r,cartgrad,cartdim,q,intgrad,intdim,NStates)
                integer,intent(in)::cartdim,intdim,NStates
                real*8,dimension(cartdim),intent(in)::r
                real*8,dimension(cartdim,NStates,NStates),intent(in)::cartgrad
                real*8,dimension(intdim),intent(out)::q
                real*8,dimension(intdim,NStates,NStates),intent(out)::intgrad
                real*8,dimension(intdim,cartdim)::B
                call WilsonBMatrixAndInternalCoordinateq(B,q,r,intdim,cartdim)
                call Cartesian2InternalNonadiabaticGradient(cartgrad,cartdim,intgrad,intdim,B,NStates)
            end subroutine Cartesian2Internal_Nonadiabatic
            
            !========== Build Wilson B matrix and transform geometry ==========
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
                end subroutine WilsonBMatrixAndInternalCoordinateq
                
                !Support WilsonBMatrixAndInternalCoordinateq
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
                    sintheta=Sqrt(1d0-costheta*costheta)
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
                    sintheta1=Sqrt(1d0-costheta1*costheta1)
                    n123=cross_product(n123,runit23)
                        n123=n123/sintheta1
                    costheta2=dot_product(runit23,n234)
                    sintheta2=Sqrt(1d0-costheta2*costheta2)
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
            !============================== End ===============================
            
            !Generate internal gradient from Wilson B matrix and Cartesian gradient
            subroutine Cartesian2InternalGradient(cartgrad,cartdim,intgrad,intdim,B)
                integer,intent(in)::cartdim,intdim
                real*8,dimension(cartdim),intent(in)::cartgrad
                real*8,dimension(intdim),intent(out)::intgrad
                real*8,dimension(intdim,cartdim),intent(in)::B
                real*8,dimension(intdim,cartdim)::GT!G = the generalized inversion of B
                GT=B
                call dGeneralizedInverseTranspose(GT,intdim,cartdim)
                intgrad=matmul(GT,cartgrad)
            end subroutine Cartesian2InternalGradient
            !In nonadiabatic process, we need the gradient of NStates order matrix
            subroutine Cartesian2InternalNonadiabaticGradient(cartgrad,cartdim,intgrad,intdim,B,NStates)
                integer,intent(in)::cartdim,intdim,NStates
                real*8,dimension(cartdim,NStates,NStates),intent(in)::cartgrad
                real*8,dimension(intdim,NStates,NStates),intent(out)::intgrad
                real*8,dimension(intdim,cartdim),intent(in)::B
                integer::i,j
                real*8,dimension(intdim,cartdim)::GT!G = the generalized inversion of B
                GT=B
                call dGeneralizedInverseTranspose(GT,intdim,cartdim)
                forall(i=1:NStates,j=1:NStates)
                    intgrad(:,i,j)=matmul(GT,cartgrad(:,i,j))
                end forall
            end subroutine Cartesian2InternalNonadiabaticGradient
        !------------------------ End ------------------------
    !=================== End ===================

    !========== Cartesian <- Internal ==========
    !I will write someday (flag)
    !=================== End ===================
!------------------- End --------------------

!Use Wilson GF method to analyze vibration from Hessian in internal coordinate
!Input:  intdim order real symmetric matrix H = Hessian in internal coordinate
!        intdim x cartdim matrix B = Wilson B matrix
!        NAtoms order array mass = mass of each atom
!Output: freqr =   real    part of vibrational angular frequencies
!        freqi = imaginary part of vibrational angular frequencies
!        H = normal coordinates in input frame
subroutine VibrationAnalysis(freqr,freqi,H,intdim,B,cartdim,mass,NAtoms)
    integer,intent(in)::intdim,cartdim,NAtoms
    real*8,dimension(intdim),intent(out)::freqr,freqi
    real*8,dimension(intdim,intdim),intent(inout)::H
    real*8,dimension(intdim,cartdim),intent(inout)::B
    real*8,dimension(NAtoms),intent(in)::mass
    integer::i
    real*8,dimension(intdim,intdim)::GF
    real*8,dimension(intdim,cartdim)::Btemp
    forall(i=1:NAtoms)
        Btemp(:,3*i-2:3*i)=B(:,3*i-2:3*i)/mass(i)
    end forall
    call syL2U(H,intdim)
    GF=matmul(matmul(Btemp,transpose(B)),H)
    call My_dgeev('V',GF,freqr,freqi,H,intdim)
end subroutine VibrationAnalysis

end module GeometryTransformation