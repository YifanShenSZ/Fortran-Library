!Common geometry transformation applied in molecule computation:
!    Standardize geometry
!    Assimilate geometry
!    Cartesian <-> internal coodinate
!    Vibrational frequency and normal mode
!
!Reference: E. B. Wilson, J. C. Decius, P. C. Cross, *Molecular viobrations: the theory of infrared and Raman vibrational spectra* (Dover, 1980)
module GeometryTransformation
    use General; use Mathematics; use LinearAlgebra
    use NonlinearOptimization
    implicit none

!Derived type
    !Example: type(IntCoordDef),allocatable,dimension(:)::IntCoordDef
    !         IntCoordDef(iIntCoord) stands for the iIntCoord-th internal coordinate
    !         IntCoordDef(iIntCoord).NMotions is the number of motions involved in this internal coordinate
    !         IntCoordDef(iIntCoord).motion(iMotion) stands for its iMotion-th motion
    !         IntCoordDef(iIntCoord).motion(iMotion).type is its type
    !         IntCoordDef(iIntCoord).motion(iMotion).coeff is its normalized linear combination coefficient
    !         IntCoordDef(iIntCoord).motion(iMotion).atom(i) is its i-th involved atom
    type InvolvedMotion
        !Currently only support stretching, bending, torsion, OutOfPlane
        !stretching: the motion coordinate is bond length atom1_atom2
        !bending   : the motion coordinate is bond angle atom1_atom2_atom3, range [0,pi]
        !            derivative encounters singularity at pi
        !torsion   : the motion coordinate is dihedral angle atom1_atom2_atom3_atom4, range (-pi,pi]
        !            specifically, angle between plane 123 and plane 234
        !            dihedral angle has same sign to n_123 x n_234 . r_23
        !            where n_abc (the normal vector of plane abc) is the unit vector along r_ab x r_bc
        !            dihedral angle value encounters discontinuity at pi
        !OutOfPlane: the motion coordinate is out of plane angle atom1_atom2_atom3_atom4, range [-pi/2,pi/2]
        !            specifically, bond 12 out of plane 234
        character*10::type
        real*8::coeff
        integer,allocatable,dimension(:)::atom
    end type InvolvedMotion
    type IntCoordDef!short for INTernal COORDinate DEFinition
        integer::NMotions
        type(InvolvedMotion),allocatable,dimension(:)::motion
    end type IntCoordDef

!GeometryTransformation module only variable
    type(IntCoordDef),allocatable,dimension(:)::GeometryTransformation_IntCoordDef!short for INTernal COORDinate DEFinition

contains
!Standardize a geometry (and optionally gradient) (optionally to a reference)
!Here we define a standard geometry as: (also called standard orientaion)
!    the centre of mass at origin
!    the rotational principle axes along xyz axes, with the smallest corresponding to x axis, 2nd to y, 3rd to z
!Note this definition does not determine geometry uniquely:
!    axes may take different positive direction as long as forming right-hand system (4 possibilities)
!Required argument: 
!    geom: geom(:,i) corresponds to [x,y,z] of i-th atom
!    mass: mass(i) is the mass of i-th atom
!    NAtoms, NStates
!Optional argument:
!    ref : uniquely define the standard geometry by the one with smallest difference to ref out of 4
!    diff: harvest mass^2 weighted || geom - ref ||_F^2
!    grad: gradient (of multiple electronic states) to transform to standard coordinate
subroutine StandardizeGeometry(geom,mass,NAtoms,NStates,ref,diff,grad)
    !Required argument
        integer,intent(in)::NAtoms,NStates
        real*8,dimension(3,NAtoms),intent(inout)::geom
        real*8,dimension(NAtoms),intent(in)::mass
    !Optional argument
        real*8,dimension(3,NAtoms),intent(in),optional::ref
        real*8,intent(out),optional::diff
        real*8,dimension(3,NAtoms,NStates,NStates),intent(inout),optional::grad
    integer::indicemin,i,istate,jstate
    real*8::mindiff,dbletemp
    real*8,dimension(3)::com!centre of mass
    real*8,dimension(3,3)::UT,moi!moment of inertia
    real*8,dimension(3,NAtoms,4)::r!To store 4 legal geometries
    !Shift centre of mass to origin, then get momentum of inertia in centre of mass frame
        com=matmul(geom,mass)/sum(mass); moi=0d0
        do i=1,NAtoms
            geom(:,i)=geom(:,i)-com
            moi=moi+mass(i)*(dot_product(geom(:,i),geom(:,i))*UnitMatrix(3)-vector_direct_product(geom(:,i),geom(:,i),3,3))
        end do
    !Diagonalize momentum of inertia
        call My_dsyev('V',moi,com,3)
        if(triple_product(moi(:,1),moi(:,2),moi(:,3))<0d0) moi(:,1)=-moi(:,1)!Should be rotation rather than reflection
        UT=transpose(moi)
    !Transform geometry (and optionally gradient) to principle axes frame
    if(present(ref)) then
        !The positive direction has 4 legal choices, determine it by comparing difference to the reference
        forall(i=1:NAtoms)
            r(:,i,1)=matmul(UT,geom(:,i))
            !Flip x and y positive directions
            r(1:2,i,2)=-r(1:2,i,1)
            r(  3,i,2)= r(  3,i,1)
            !Flip x and z positive directions
            r(1,i,3)=-r(1,i,1)
            r(2,i,3)= r(2,i,1)
            r(3,i,3)=-r(3,i,1)
            !Flip y and z positive directions
            r(  1,i,4)= r(  1,i,1)
            r(2:3,i,4)=-r(2:3,i,1)
        end forall
        !Which one has smallest difference?
        indicemin=1; mindiff=difference(r(:,:,1))
        do i=2,4
            dbletemp=difference(r(:,:,i))
            if(dbletemp<mindiff) then
                indicemin=i; mindiff=dbletemp
            end if
        end do
        if(present(diff)) diff=mindiff!Harvest mass^2 weighted || geom - ref ||_2^2
        geom=r(:,:,indicemin)!Determine the unique geometry
        select case(indicemin)!Determine the principle axes
        case(2); moi(:,1:2)=-moi(:,1:2)
        case(3); moi(:,1)=-moi(:,1); moi(:,3)=-moi(:,3)
        case(4); moi(:,2:3)=-moi(:,2:3)
        end select
        if(present(grad)) UT=transpose(moi)!Update U^T accordingly
    else
        forall(i=1:NAtoms); geom(:,i)=matmul(UT,geom(:,i)); end forall
    end if
    if(present(grad)) then
        forall(i=1:NAtoms,istate=1:NStates,jstate=1:NStates,istate>=jstate)
            grad(:,i,istate,jstate)=matmul(UT,grad(:,i,istate,jstate))
        end forall
    end if
    contains
    real*8 function difference(geom)
        real*8,dimension(3,NAtoms),intent(in)::geom
        integer::i; real*8,dimension(3,NAtoms)::temp
        temp=geom-ref
        forall(i=1:NAtoms); temp(:,i)=temp(:,i)*mass(i); end forall
        difference=dgeFrobeniusSquare(temp,3,NAtoms)
    end function difference
end subroutine StandardizeGeometry

!Assimilate a geometry to a reference
!Required argument: 
!    geom: geom(:,i) corresponds to [x,y,z] of i-th atom
!    ref : the reference
!    mass: mass(i) is the mass of i-th atom
!    NAtoms
!Optional argument:
!    diff: harvest mass^2 weighted || geom - ref ||_F^2
!    init: initial value of rotation, see rotation section
subroutine AssimilateGeometry(geom,ref,mass,NAtoms,diff,init)
    !Required argument
        integer,intent(in)::NAtoms
        real*8,dimension(3,NAtoms),intent(inout)::geom
        real*8,dimension(3,NAtoms),intent(in)::ref
        real*8,dimension(NAtoms),intent(in)::mass
    !Optional argument
        real*8,intent(out),optional::diff
        real*8,dimension(3),intent(in),optional::init
    !Rotation
        real*8::diffmin,difftemp,sintheta
        real*8,dimension(3)::qind
        real*8,dimension(4)::q,qmin
        real*8,dimension(3,NAtoms)::geommin,geomtemp
    integer::i!Work variable
    !Shift the centre of mass of the input geometry to the reference
        qind=matmul(ref-geom,mass)/sum(mass)
        forall(i=1:NAtoms); geom(:,i)=geom(:,i)+qind; end forall
    !Rotate geometry to minimize mass^2 weighted || geom - ref ||_F^2
        !Rotation is done by unit quaternion q=[cos(alpha/2),sin(alpha/2)*axis]
        !Here we let the independent variables be alpha/2, theta, phi
        !    where axis=[sin(theta)cos(phi),sin(theta)sin(phi),cos(theta)]
        if(present(init)) then; qind=init!Initial value = user input
        else!or random value
            !No rotation
            qmin=[1d0,0d0,0d0,0d0]
            diffmin=difference(geom)
            !Randomly sample 1000000 orientations
            do i=1,1000000
                q=RandomUnitQuaternion()
                geomtemp=geom; call Rotate(q,geomtemp,NAtoms); difftemp=difference(geomtemp)
                if(difftemp<diffmin) then; qmin=q; diffmin=difftemp; end if
            end do
            !Try [0,360) degree with step length = 1 degree along x, y, z
                !0 degree along x
                    qind=[0d0,pid2,0d0]
                    q(1)=dCos(qind(1)); sintheta=dSin(qind(2))
                    q(2:4)=dSin(qind(1))*[sintheta*dCos(qind(3)),sintheta*dSin(qind(3)),dCos(qind(2))]
                    geomtemp=geom; call Rotate(q,geomtemp,NAtoms); difftemp=difference(geomtemp)
                    if(difftemp<diffmin) then; qmin=q; diffmin=difftemp; end if
                    do i=1,359
                        qind(1)=qind(1)+DegInRad/2d0
                        q(1)=dCos(qind(1)); sintheta=dSin(qind(2))
                        q(2:4)=dSin(qind(1))*[sintheta*dCos(qind(3)),sintheta*dSin(qind(3)),dCos(qind(2))]
                        geomtemp=geom; call Rotate(q,geomtemp,NAtoms); difftemp=difference(geomtemp)
                        if(difftemp<diffmin) then; qmin=q; diffmin=difftemp; end if
                    end do
                !0 degree along y
                    qind=[0d0,pid2,pid2]
                    q(1)=dCos(qind(1)); sintheta=dSin(qind(2))
                    q(2:4)=dSin(qind(1))*[sintheta*dCos(qind(3)),sintheta*dSin(qind(3)),dCos(qind(2))]
                    geomtemp=geom; call Rotate(q,geomtemp,NAtoms); difftemp=difference(geomtemp)
                    if(difftemp<diffmin) then; qmin=q; diffmin=difftemp; end if
                    do i=1,359
                        qind(1)=qind(1)+DegInRad/2d0
                        q(1)=dCos(qind(1)); sintheta=dSin(qind(2))
                        q(2:4)=dSin(qind(1))*[sintheta*dCos(qind(3)),sintheta*dSin(qind(3)),dCos(qind(2))]
                        geomtemp=geom; call Rotate(q,geomtemp,NAtoms); difftemp=difference(geomtemp)
                        if(difftemp<diffmin) then; qmin=q; diffmin=difftemp; end if
                    end do
                !0 degree along z
                    qind=0d0
                    q(1)=dCos(qind(1)); sintheta=dSin(qind(2))
                    q(2:4)=dSin(qind(1))*[sintheta*dCos(qind(3)),sintheta*dSin(qind(3)),dCos(qind(2))]
                    geomtemp=geom; call Rotate(q,geomtemp,NAtoms); difftemp=difference(geomtemp)
                    if(difftemp<diffmin) then; qmin=q; diffmin=difftemp; end if
                    do i=1,359
                        qind(1)=qind(1)+DegInRad/2d0
                        q(1)=dCos(qind(1)); sintheta=dSin(qind(2))
                        q(2:4)=dSin(qind(1))*[sintheta*dCos(qind(3)),sintheta*dSin(qind(3)),dCos(qind(2))]
                        geomtemp=geom; call Rotate(q,geomtemp,NAtoms); difftemp=difference(geomtemp)
                        if(difftemp<diffmin) then; qmin=q; diffmin=difftemp; end if
                    end do
            !Use the one with smallest difference
            if(qmin(1)==1d0) then!No rotation, axis is arbitrary
                qind=[0d0,pid2,pid2]!Arbitrarily let axis = y
            else
                qind(1)=dACos(qmin(1))
                q(1:3)=qmin(2:4)/dSqrt(1d0-qmin(1)*qmin(1))!Save axis
                qind(2)=dACos(q(3))
                q(1:2)=q(1:2)/dSqrt(1d0-q(3)*q(3))
                qind(3)=dACos(q(1)); if(q(2)<0d0) qind(3)=-qind(3)
            end if
        end if
        !Search for an optimal rotation
        call TrustRegion(Residue,qind,3*NAtoms,3,Warning=.false.)
        !Perform the rotation
        q(1)=dCos(qind(1)); sintheta=dSin(qind(2))
        q(2:4)=dSin(qind(1))*[sintheta*dCos(qind(3)),sintheta*dSin(qind(3)),dCos(qind(2))]
        call Rotate(q,geom,NAtoms)
    if(present(diff)) diff=difference(geom)
    contains
    real*8 function difference(geom)
        real*8,dimension(3,NAtoms),intent(in)::geom
        integer::i; real*8,dimension(3,NAtoms)::temp
        temp=geom-ref
        forall(i=1:NAtoms); temp(:,i)=temp(:,i)*mass(i); end forall
        difference=dgeFrobeniusSquare(temp,3,NAtoms)
    end function difference
    subroutine Residue(res,qind,cartdim,qdim)
        integer,intent(in)::cartdim,qdim
        real*8,dimension(cartdim),intent(out)::res
        real*8,dimension(qdim),intent(in)::qind
        real*8::sintheta; real*8,dimension(4)::q
        real*8,dimension(3,NAtoms)::geomtemp
        q(1)=dCos(qind(1)); sintheta=dSin(qind(2))
        q(2:4)=dSin(qind(1))*[sintheta*dCos(qind(3)),sintheta*dSin(qind(3)),dCos(qind(2))]
        geomtemp=geom; call Rotate(q,geomtemp,NAtoms); geomtemp=geomtemp-ref
        forall(i=1:NAtoms); geomtemp(:,i)=geomtemp(:,i)*mass(i); end forall
        res=reshape(geomtemp,[cartdim])
    end subroutine Residue
end subroutine AssimilateGeometry

!---------- Cartesian <-> Internal ----------
    !An interal coordinate is the linear combination of several translationally and rotationally invariant displacements
    !    but only displacements under same unit can be combined, i.e. you must treat length and angle separately
    !    unless appropriate metric tensor is applied
    !It is OK to define more than 3NAtoms-6 (or 3NAtoms-5 for linear molecule) internal coordinates,
    !    but only 3NAtoms-6 (or 3NAtoms-5 for linear molecule) partial derivatives are independent
    !Although the transformation from Cartesian coordinate to internal coordinate is not necessarily linear
    !    for infinitesimal displacement it is linear, corresponding to a matrix form: dq = B . dr
    !    where dq is internal coordinate differentiation, dr is Cartesian coordinate differentiation
    !    B is Jacobian(q,r) (historically called Wilson B matrix)
    !r is a 3NAtoms order vector with r[3*i-2:3*i] corresponding to the coordinate of i-th atom
    !Nomenclature:
    !    cartdim & intdim: Cartesian & internal space dimensionality
    !    cartgrad & intgrad: Cartesian & internal coordinate gradient (cartdim & intdim x NStates x NStates 3rd-order tensor)

    !Define internal coordinate, return the internal space dimensionality
    !Input:  format: internal coordinate definition format (Available: Columbus7, default)
    !        (optional) file: (default = 'intcfl' for Columbus7, 'IntCoordDef' for default) internal coordinate definition file name
    !Output: the internal space dimensionality
    !        also set the module-wide variable GeometryTransformation_IntCoordDef
    !        which will be refered by all routines in this section
    !See InvolvedMotion in 'Derived type' section for available types and ordering of atoms
    integer function DefineInternalCoordinate(format,file)
        character(*),intent(in)::format
        character(*),optional,intent(in)::file
        integer::intdim
        if(allocated(GeometryTransformation_IntCoordDef)) deallocate(GeometryTransformation_IntCoordDef)
        select case(format)
        case('Columbus7'); call Columbus7()
        case default; call default()
        end select
        DefineInternalCoordinate=intdim
        contains
        !First line is always 'TEXAS'
        !New internal coordinate line starts with 'K'
        subroutine Columbus7()
            integer::NDef
            integer,allocatable,dimension(:)::NewLine
            character*10,allocatable,dimension(:)::MotionType
            character*24::chartemp; integer::i,j,k; real*8::dbletemp
            if(present(file)) then; open(unit=99,file=file,status='old')
            else; open(unit=99,file='intcfl',status='old'); end if
                !The number of motion definition lines & internal coordinates
                    NDef=0; intdim=0; read(99,*)
                    do
                        read(99,'(A24)',iostat=i)chartemp
                        if(i/=0&!End of file or no definition
                        .or.(index(chartemp,'STRE')==0.and.index(chartemp,'BEND')==0&
                        .and.index(chartemp,'TORS')==0.and.index(chartemp,'OUT' )==0)) exit
                        NDef=NDef+1
                        if(scan(chartemp,'K')==1) intdim=intdim+1
                    end do; rewind 99
                !New internal coordinate lines & motions of line
                    allocate(NewLine(intdim+1)); NewLine(intdim+1)=NDef+1
                    allocate(MotionType(NDef))
                    k=1; read(99,*)
                    do i=1,NDef
                        read(99,'(A24)')chartemp
                        if(scan(chartemp,'K')==1) then; NewLine(k)=i; k=k+1; end if
                        if(index(chartemp,'STRE')>0) then; MotionType(i)='stretching'
                        else if(index(chartemp,'BEND')>0) then; MotionType(i)='bending'
                        else if(index(chartemp,'TORS')>0) then; MotionType(i)='torsion'
                        else if(index(chartemp,'OUT')>0) then; MotionType(i)='OutOfPlane'; end if
                    end do; rewind 99
                !Finally read internal coordinate definition. Linear combinations are normalized
                    allocate(GeometryTransformation_IntCoordDef(intdim))
                    k=1; read(99,*)
                    do i=1,intdim
                        GeometryTransformation_IntCoordDef(i).NMotions=NewLine(i+1)-NewLine(i)
                        allocate(GeometryTransformation_IntCoordDef(i).motion(GeometryTransformation_IntCoordDef(i).NMotions))
                        if(GeometryTransformation_IntCoordDef(i).NMotions==1) then
                            GeometryTransformation_IntCoordDef(i).motion(1).type=MotionType(k)
                            GeometryTransformation_IntCoordDef(i).motion(1).coeff=1d0
                            select case(MotionType(k))
                            case('stretching')
                                allocate(GeometryTransformation_IntCoordDef(i).motion(1).atom(2))
                                read(99,'(A28,I5,1x,I9)')chartemp,&
                                GeometryTransformation_IntCoordDef(i).motion(1).atom
                            case('bending')
                                allocate(GeometryTransformation_IntCoordDef(i).motion(1).atom(3))
                                read(99,'(A28,I6,1x,I9,1x,I9)')chartemp,&
                                GeometryTransformation_IntCoordDef(i).motion(1).atom(1),&
                                GeometryTransformation_IntCoordDef(i).motion(1).atom(3),&
                                GeometryTransformation_IntCoordDef(i).motion(1).atom(2)
                            case('torsion')
                                allocate(GeometryTransformation_IntCoordDef(i).motion(1).atom(4))
                                read(99,'(A28,I6,1x,I9,1x,I9,1x,I9)')chartemp,&
                                GeometryTransformation_IntCoordDef(i).motion(1).atom
                            case('OutOfPlane')
                                allocate(GeometryTransformation_IntCoordDef(i).motion(1).atom(4))
                                read(99,'(A28,I6,1x,I9,1x,I9,1x,I9)')chartemp,&
                                GeometryTransformation_IntCoordDef(i).motion(1).atom(1),&
                                GeometryTransformation_IntCoordDef(i).motion(1).atom(4),&
                                GeometryTransformation_IntCoordDef(i).motion(1).atom(2),&
                                GeometryTransformation_IntCoordDef(i).motion(1).atom(3)
                            case default; write(*,*)'Program abort: unsupported internal coordinate type '//trim(adjustl(MotionType(k))); stop
                            end select
                            k=k+1
                        else
                            dbletemp=0d0
                            do j=1,GeometryTransformation_IntCoordDef(i).NMotions
                                GeometryTransformation_IntCoordDef(i).motion(j).type=MotionType(k)
                                select case(MotionType(k))
                                case('stretching')
                                    allocate(GeometryTransformation_IntCoordDef(i).motion(j).atom(2))
                                    read(99,'(A10,F10.7,8x,I5,1x,I9)')chartemp,&
                                    GeometryTransformation_IntCoordDef(i).motion(j).coeff,&
                                    GeometryTransformation_IntCoordDef(i).motion(j).atom
                                case('bending')
                                    allocate(GeometryTransformation_IntCoordDef(i).motion(j).atom(3))
                                    read(99,'(A10,F10.7,8x,I6,1x,I9,1x,I9)')chartemp,&
                                    GeometryTransformation_IntCoordDef(i).motion(j).coeff,&
                                    GeometryTransformation_IntCoordDef(i).motion(j).atom(1),&
                                    GeometryTransformation_IntCoordDef(i).motion(j).atom(3),&
                                    GeometryTransformation_IntCoordDef(i).motion(j).atom(2)
                                case('torsion')
                                    allocate(GeometryTransformation_IntCoordDef(i).motion(j).atom(4))
                                    read(99,'(A10,F10.7,8x,I6,1x,I9,1x,I9,1x,I9)')chartemp,&
                                    GeometryTransformation_IntCoordDef(i).motion(j).coeff,&
                                    GeometryTransformation_IntCoordDef(i).motion(j).atom
                                case('OutOfPlane')
                                    allocate(GeometryTransformation_IntCoordDef(i).motion(j).atom(4))
                                    read(99,'(A10,F10.7,8x,I6,1x,I9,1x,I9,1x,I9)')chartemp,&
                                    GeometryTransformation_IntCoordDef(i).motion(j).coeff,&
                                    GeometryTransformation_IntCoordDef(i).motion(j).atom(1),&
                                    GeometryTransformation_IntCoordDef(i).motion(j).atom(4),&
                                    GeometryTransformation_IntCoordDef(i).motion(j).atom(2),&
                                    GeometryTransformation_IntCoordDef(i).motion(j).atom(3)
                                case default; write(*,*)'Program abort: unsupported internal coordinate type '//trim(adjustl(MotionType(k))); stop
                                end select
                                k=k+1
                                dbletemp=dbletemp+GeometryTransformation_IntCoordDef(i).motion(j).coeff*GeometryTransformation_IntCoordDef(i).motion(j).coeff
                            end do
                            dbletemp=Sqrt(dbletemp)
                            forall(j=1:GeometryTransformation_IntCoordDef(i).NMotions)
                                GeometryTransformation_IntCoordDef(i).motion(j).coeff=&
                                GeometryTransformation_IntCoordDef(i).motion(j).coeff/dbletemp
                            end forall
                        end if
                    end do
            close(99)
            deallocate(NewLine); deallocate(MotionType)!Clean up
        end subroutine Columbus7
        !First 6 spaces of a line are reserved to indicate the start of new internal coordinate
        !Example:
        ! coor |   coeff   |    type     |      atom
        !--------------------------------------------------
        !     1    1.000000    stretching     1     2          # Comment
        !          1.000000    stretching     1     3
        !     2    1.000000    stretching     1     2
        !         -1.000000    stretching     1     3
        !     3    1.000000       bending     2     1     3
        subroutine default()
            integer::NDef
            integer,allocatable,dimension(:)::NewLine
            character*10,allocatable,dimension(:)::MotionType
            character*10::chartemp; integer::i,j,k,l; real*8::dbletemp
            if(present(file)) then; open(unit=99,file=file,status='old')
            else; open(unit=99,file='IntCoordDef',status='old'); end if
                !The number of motion definition lines & internal coordinates
                    NDef=0
                    do
                        read(99,'(I6)',iostat=i)j
                        if(i/=0) exit!End of file
                        NDef=NDef+1
                        if(j>0) intdim=j
                    end do; rewind 99
                !New internal coordinate lines & motions of line
                    allocate(NewLine(intdim+1)); NewLine(intdim+1)=NDef+1
                    allocate(MotionType(NDef))
                    k=1
                    do i=1,NDef
                        read(99,'(I6)',advance='no')j
                        if(j>0) then; NewLine(k)=i; k=k+1; end if
                        read(99,*)dbletemp,MotionType(i)
                    end do; rewind 99
                !Finally read internal coordinate definition. Linear combinations are normalized
                    allocate(GeometryTransformation_IntCoordDef(intdim))
                    k=1
                    do i=1,intdim
                        GeometryTransformation_IntCoordDef(i).NMotions=NewLine(i+1)-NewLine(i)
                        allocate(GeometryTransformation_IntCoordDef(i).motion(GeometryTransformation_IntCoordDef(i).NMotions))
                        dbletemp=0d0
                        do j=1,GeometryTransformation_IntCoordDef(i).NMotions
                            read(99,'(I6)',advance='no')l
                            GeometryTransformation_IntCoordDef(i).motion(j).type=MotionType(k)
                            select case(MotionType(k))
                            case('stretching')
                                allocate(GeometryTransformation_IntCoordDef(i).motion(j).atom(2))
                                read(99,*)GeometryTransformation_IntCoordDef(i).motion(j).coeff,chartemp,&
                                          GeometryTransformation_IntCoordDef(i).motion(j).atom
                            case('bending')
                                allocate(GeometryTransformation_IntCoordDef(i).motion(j).atom(3))
                                read(99,*)GeometryTransformation_IntCoordDef(i).motion(j).coeff,chartemp,&
                                          GeometryTransformation_IntCoordDef(i).motion(j).atom
                            case('torsion')
                                allocate(GeometryTransformation_IntCoordDef(i).motion(j).atom(4))
                                read(99,*)GeometryTransformation_IntCoordDef(i).motion(j).coeff,chartemp,&
                                          GeometryTransformation_IntCoordDef(i).motion(j).atom
                            case('OutOfPlane')
                                allocate(GeometryTransformation_IntCoordDef(i).motion(j).atom(4))
                                read(99,*)GeometryTransformation_IntCoordDef(i).motion(j).coeff,chartemp,&
                                          GeometryTransformation_IntCoordDef(i).motion(j).atom
                            case default; write(*,*)'Program abort: unsupported internal coordinate type '//trim(adjustl(MotionType(k))); stop
                            end select
                            k=k+1
                            dbletemp=dbletemp+GeometryTransformation_IntCoordDef(i).motion(j).coeff*GeometryTransformation_IntCoordDef(i).motion(j).coeff
                        end do
                        dbletemp=Sqrt(dbletemp)
                        forall(j=1:GeometryTransformation_IntCoordDef(i).NMotions)
                            GeometryTransformation_IntCoordDef(i).motion(j).coeff=&
                            GeometryTransformation_IntCoordDef(i).motion(j).coeff/dbletemp
                        end forall
                    end do
            close(99)
            deallocate(NewLine); deallocate(MotionType)!Clean up
        end subroutine default
    end function DefineInternalCoordinate

    !========== Cartesian -> Internal ==========
        !Convert r to q
        subroutine InternalCoordinateq(r,q,cartdim,intdim)
            integer,intent(in)::cartdim,intdim
            real*8,dimension(cartdim),intent(in)::r
            real*8,dimension(intdim),intent(out)::q
            integer::iIntC,iMotion
            q=0d0
            do iIntC=1,intdim
                do iMotion=1,GeometryTransformation_IntCoordDef(iIntC).NMotions
                    select case(GeometryTransformation_IntCoordDef(iIntC).motion(iMotion).type)
                    case('stretching')
                        q(iIntC)=q(iIntC)&
                            +GeometryTransformation_IntCoordDef(iIntC).motion(iMotion).coeff&
                            *stretching(r,GeometryTransformation_IntCoordDef(iIntC).motion(iMotion).atom,cartdim)
                    case('bending')
                        q(iIntC)=q(iIntC)&
                            +GeometryTransformation_IntCoordDef(iIntC).motion(iMotion).coeff&
                            *bending(r,GeometryTransformation_IntCoordDef(iIntC).motion(iMotion).atom,cartdim)
                    case('torsion')
                        q(iIntC)=q(iIntC)&
                            +GeometryTransformation_IntCoordDef(iIntC).motion(iMotion).coeff&
                            *torsion(r,GeometryTransformation_IntCoordDef(iIntC).motion(iMotion).atom,cartdim)
                    case('OutOfPlane')
                        q(iIntC)=q(iIntC)&
                            +GeometryTransformation_IntCoordDef(iIntC).motion(iMotion).coeff&
                            *OutOfPlane(r,GeometryTransformation_IntCoordDef(iIntC).motion(iMotion).atom,cartdim)
                    end select
                end do
            end do
            contains
            !Transform from Cartesian coordinate r to a certain motion coordinate q, atom defines which atoms are involved
            !For stretching, q = bond length
            real*8 function stretching(r,atom,cartdim)
                integer,intent(in)::cartdim
                real*8,dimension(cartdim),intent(in)::r
                integer,dimension(2),intent(in)::atom
                real*8,dimension(3)::r12
                r12=r(3*atom(2)-2:3*atom(2))-r(3*atom(1)-2:3*atom(1))
                stretching=Norm2(r12)
            end function stretching
            !For bending, q = bond angle
            real*8 function bending(r,atom,cartdim)
                integer,intent(in)::cartdim
                real*8,dimension(cartdim),intent(in)::r
                integer,dimension(3),intent(in)::atom
                real*8,dimension(3)::runit21,runit23
                runit21=r(3*atom(1)-2:3*atom(1))-r(3*atom(2)-2:3*atom(2))
                    runit21=runit21/Norm2(runit21)
                runit23=r(3*atom(3)-2:3*atom(3))-r(3*atom(2)-2:3*atom(2))
                    runit23=runit23/Norm2(runit23)
                bending=acos(dot_product(runit21,runit23))
            end function bending
            !For torsion, q = dihedral angle
            real*8 function torsion(r,atom,cartdim)
                integer,intent(in)::cartdim
                real*8,dimension(cartdim),intent(in)::r
                integer,dimension(4),intent(in)::atom
                real*8,dimension(3)::r12,r23,r34,n123,n234
                r12=r(3*atom(2)-2:3*atom(2))-r(3*atom(1)-2:3*atom(1))
                r23=r(3*atom(3)-2:3*atom(3))-r(3*atom(2)-2:3*atom(2))
                r34=r(3*atom(4)-2:3*atom(4))-r(3*atom(3)-2:3*atom(3))
                n123=cross_product(r12,r23); n123=n123/Norm2(n123)
                n234=cross_product(r23,r34); n234=n234/Norm2(n234)
                torsion=acos(dot_product(n123,n234))
                if(triple_product(n123,n234,r23)<0d0) torsion=-torsion
            end function torsion
            !For out of plane, q = out of plane angle
            real*8 function OutOfPlane(r,atom,cartdim)
                integer,intent(in)::cartdim
                real*8,dimension(cartdim),intent(in)::r
                integer,dimension(4),intent(in)::atom
                real*8,dimension(3)::r21,r23,r24
                r21=r(3*atom(1)-2:3*atom(1))-r(3*atom(2)-2:3*atom(2))
                r23=r(3*atom(3)-2:3*atom(3))-r(3*atom(2)-2:3*atom(2))
                r24=r(3*atom(4)-2:3*atom(4))-r(3*atom(2)-2:3*atom(2))
                r23=cross_product(r23,r24)
                OutOfPlane=asin(dot_product(r23/norm2(r23),r21/norm2(r21)))
            end function OutOfPlane
        end subroutine InternalCoordinateq

        !Convert r & cartgrad to q & intgrad
        subroutine Cartesian2Internal(r,cartgrad,q,intgrad,cartdim,intdim,NStates)
            integer,intent(in)::cartdim,intdim,NStates
            real*8,dimension(cartdim),intent(in)::r
            real*8,dimension(cartdim,NStates,NStates),intent(in)::cartgrad
            real*8,dimension(intdim),intent(out)::q
            real*8,dimension(intdim,NStates,NStates),intent(out)::intgrad
            integer::i,j; real*8,dimension(intdim,cartdim)::B
            call WilsonBMatrixAndInternalCoordinateq(r,B,q,cartdim,intdim)
            call dGeneralizedInverseTranspose(B,intdim,cartdim)
            forall(i=1:NStates,j=1:NStates); intgrad(:,i,j)=matmul(B,cartgrad(:,i,j)); end forall
        end subroutine Cartesian2Internal
        
        !From r, generate B & q
        subroutine WilsonBMatrixAndInternalCoordinateq(r,B,q,cartdim,intdim)
            integer,intent(in)::cartdim,intdim
            real*8,dimension(cartdim),intent(in)::r
            real*8,dimension(intdim,cartdim),intent(out)::B
            real*8,dimension(intdim),intent(out)::q
            integer::iIntC,iMotion; real*8::qMotion; real*8,dimension(cartdim)::BRowVector
            B=0d0; q=0d0
            do iIntC=1,intdim
                do iMotion=1,GeometryTransformation_IntCoordDef(iIntC).NMotions
                    select case(GeometryTransformation_IntCoordDef(iIntC).motion(iMotion).type)
                    case('stretching'); call bAndStretching(BRowVector,qMotion,r,GeometryTransformation_IntCoordDef(iIntC).motion(iMotion).atom,cartdim)
                    case('bending')   ; call bAndBending   (BRowVector,qMotion,r,GeometryTransformation_IntCoordDef(iIntC).motion(iMotion).atom,cartdim)
                    case('torsion')   ; call bAndTorsion   (BRowVector,qMotion,r,GeometryTransformation_IntCoordDef(iIntC).motion(iMotion).atom,cartdim)
                    case('OutOfPlane'); call bAndOutOfPlane(BRowVector,qMotion,r,GeometryTransformation_IntCoordDef(iIntC).motion(iMotion).atom,cartdim)
                    end select
                    B(iIntC,:)=B(iIntC,:)+GeometryTransformation_IntCoordDef(iIntC).motion(iMotion).coeff*BRowVector
                    q(iIntC)  =q(iIntC)  +GeometryTransformation_IntCoordDef(iIntC).motion(iMotion).coeff*qMotion
                end do
            end do
            contains
            !Generate the transformation vector b from dr to dq: b . dr = dq
            !Transform from Cartesian coordinate r to a certain motion coordinate q
            !Internal coordinate is the linear combination of several motions,
            !so b contributes (but not necessarily equals) to one row of Wilson B matrix
            ! d( i-th internal coordinate ) = ( i-th row vector of B ) . dr
            !For stretching, q = bond length
            subroutine bAndStretching(b,q,r,atom,cartdim)
                integer,intent(in)::cartdim
                real*8,dimension(cartdim),intent(out)::b
                real*8,intent(out)::q
                real*8,dimension(cartdim),intent(in)::r
                integer,dimension(2),intent(in)::atom
                real*8,dimension(3)::runit12
                b=0d0!Initialize
                runit12=r(3*atom(2)-2:3*atom(2))-r(3*atom(1)-2:3*atom(1))
                q=Norm2(runit12)
                runit12=runit12/q
                b(3*atom(1)-2:3*atom(1))=-runit12
                b(3*atom(2)-2:3*atom(2))=runit12
            end subroutine bAndStretching
            !For bending, q = bond angle
            subroutine bAndBending(b,q,r,atom,cartdim)
                integer,intent(in)::cartdim
                real*8,dimension(cartdim),intent(out)::b
                real*8,intent(out)::q
                real*8,dimension(cartdim),intent(in)::r
                integer,dimension(3),intent(in)::atom
                real*8::r21,r23,costheta,sintheta
                real*8,dimension(3)::runit21,runit23
                b=0d0!Initialize
                !Prepare
                runit21=r(3*atom(1)-2:3*atom(1))-r(3*atom(2)-2:3*atom(2))
                    r21=Norm2(runit21); runit21=runit21/r21
                runit23=r(3*atom(3)-2:3*atom(3))-r(3*atom(2)-2:3*atom(2))
                    r23=Norm2(runit23); runit23=runit23/r23
                costheta=dot_product(runit21,runit23); sintheta=dSqrt(1d0-costheta*costheta)
                !Output
                b(3*atom(1)-2:3*atom(1))=(costheta*runit21-runit23)/(sintheta*r21)
                b(3*atom(3)-2:3*atom(3))=(costheta*runit23-runit21)/(sintheta*r23)
                b(3*atom(2)-2:3*atom(2))=-b(3*atom(1)-2:3*atom(1))-b(3*atom(3)-2:3*atom(3))
                q=acos(costheta)
            end subroutine bAndBending
            !For torsion, q = dihedral angle
            subroutine bAndTorsion(b,q,r,atom,cartdim)
                integer,intent(in)::cartdim
                real*8,dimension(cartdim),intent(out)::b
                real*8,intent(out)::q
                real*8,dimension(cartdim),intent(in)::r
                integer,dimension(4),intent(in)::atom
                real*8::r12,r23,r34,sin123,cos123,sin234,cos234
                real*8,dimension(3)::runit12,runit23,runit34,n123,n234
                b=0d0!Initialize
                !Prepare
                runit12=r(3*atom(2)-2:3*atom(2))-r(3*atom(1)-2:3*atom(1))
                r12=Norm2(runit12); runit12=runit12/r12
                runit23=r(3*atom(3)-2:3*atom(3))-r(3*atom(2)-2:3*atom(2))
                r23=Norm2(runit23); runit23=runit23/r23
                runit34=r(3*atom(4)-2:3*atom(4))-r(3*atom(3)-2:3*atom(3))
                r34=Norm2(runit34); runit34=runit34/r34
                cos123=-dot_product(runit12,runit23); sin123=dSqrt(1d0-cos123*cos123)
                n123=cross_product(runit12,runit23)/sin123
                cos234=-dot_product(runit23,runit34); sin234=dSqrt(1d0-cos234*cos234)
                n234=cross_product(runit23,runit34)/sin234
                !Output
                b(3*atom(1)-2:3*atom(1))=-n123/(r12*sin123)
                b(3*atom(2)-2:3*atom(2))=(r23-r12*cos123)/(r12*r23*sin123)*n123-cos234/(r23*sin234)*n234
                b(3*atom(3)-2:3*atom(3))=(r34*cos234-r23)/(r23*r34*sin234)*n234+cos123/(r23*sin123)*n123
                b(3*atom(4)-2:3*atom(4))= n234/(r34*sin234)
                q=acos(dot_product(n123,n234))
                if(triple_product(n123,n234,runit23)<0d0) q=-q
            end subroutine bAndTorsion
            !For out of plane, q = out of plane angle
            subroutine bAndOutOfPlane(b,q,r,atom,cartdim)
                integer,intent(in)::cartdim
                real*8,dimension(cartdim),intent(out)::b
                real*8,intent(out)::q
                real*8,dimension(cartdim),intent(in)::r
                integer,dimension(4),intent(in)::atom
                real*8::r21,r23,r24,sin324,cos324,sin324sq,sintheta,costheta,tantheta
                real*8,dimension(3)::runit21,runit23,runit24
                b=0d0!Initialize
                !Prepare
                runit21=r(3*atom(1)-2:3*atom(1))-r(3*atom(2)-2:3*atom(2))
                r21=Norm2(runit21); runit21=runit21/r21
                runit23=r(3*atom(3)-2:3*atom(3))-r(3*atom(2)-2:3*atom(2))
                r23=Norm2(runit23); runit23=runit23/r23
                runit24=r(3*atom(4)-2:3*atom(4))-r(3*atom(2)-2:3*atom(2))
                r24=Norm2(runit24); runit24=runit24/r24
                cos324=dot_product(runit23,runit24)
                sin324=dSqrt(1d0-cos324*cos324); sin324sq=sin324*sin324
                sintheta=triple_product(runit23,runit24,runit21)/sin324
                costheta=dSqrt(1d0-sintheta*sintheta); tantheta=sintheta/costheta
                !Output
                b(3*atom(1)-2:3*atom(1))=(cross_product(runit23,runit24)/costheta/sin324-tantheta*runit21)/r21
                b(3*atom(3)-2:3*atom(3))=(cross_product(runit24,runit21)/costheta/sin324-tantheta/sin324sq*(runit23-cos324*runit24))/r23
                b(3*atom(4)-2:3*atom(4))=(cross_product(runit21,runit23)/costheta/sin324-tantheta/sin324sq*(runit24-cos324*runit23))/r24
                b(3*atom(2)-2:3*atom(2))=-b(3*atom(1)-2:3*atom(1))-b(3*atom(3)-2:3*atom(3))-b(3*atom(4)-2:3*atom(4))
                q=asin(sintheta)
            end subroutine bAndOutOfPlane
        end subroutine WilsonBMatrixAndInternalCoordinateq
    !=================== End ===================

    !========== Cartesian <- Internal ==========
        !Convert q to r
        !Please note that r may vary with arbitrary translation & rotation
        !Optional argument: r0: (default = random) initial guess of r
        subroutine CartesianCoordinater(q,r,intdim,cartdim,r0)
            !Required argument
                integer,intent(in)::intdim,cartdim
                real*8,dimension(intdim),intent(in)::q
                real*8,dimension(cartdim),intent(out)::r
            !Optional argument
                real*8,dimension(cartdim),intent(in),optional::r0
            if(present(r0)) then; r=r0!Initial guess
            else; call random_number(r); end if
            call TrustRegion(Residue,r,cartdim,cartdim,Jacobian=Jacobian,Precision=1d-10)
            contains
            subroutine Residue(res,r,dim,cartdim)
                integer,intent(in)::dim,cartdim
                real*8,dimension(dim),intent(out)::res
                real*8,dimension(cartdim),intent(in)::r
                call InternalCoordinateq(r,res(1:intdim),cartdim,intdim)
                res(1:intdim)=res(1:intdim)-q
                res(intdim+1:dim)=0d0
            end subroutine Residue
            integer function Jacobian(Jacob,r,dim,cartdim)
                integer,intent(in)::dim,cartdim
                real*8,dimension(dim,cartdim),intent(out)::Jacob
                real*8,dimension(cartdim),intent(in)::r
                real*8,dimension(intdim)::qtemp
                call WilsonBMatrixAndInternalCoordinateq(r,Jacob(1:intdim,:),qtemp,cartdim,intdim)
                Jacob(intdim+1:dim,:)=0d0
                Jacobian=0!Return 0
            end function Jacobian
        end subroutine CartesianCoordinater

        !Convert q & intgrad to r & cartgrad
        !Please note that r may vary with arbitrary translation & rotation
        !so cartgrad varies with same rotation to r
        !Optional argument: r0: (default = random) initial guess of r
        subroutine Internal2Cartesian(q,intgrad,r,cartgrad,intdim,cartdim,NStates,r0)
            !Required argument
                integer,intent(in)::intdim,cartdim,NStates
                real*8,dimension(intdim),intent(in)::q
                real*8,dimension(intdim,NStates,NStates),intent(in)::intgrad
                real*8,dimension(cartdim),intent(out)::r
                real*8,dimension(cartdim,NStates,NStates),intent(out)::cartgrad
            !Optional argument
                real*8,dimension(cartdim),intent(in),optional::r0
            integer::i,j; real*8,dimension(intdim)::qtemp; real*8,dimension(intdim,cartdim)::B
            if(present(r0)) then; call CartesianCoordinater(q,r,intdim,cartdim,r0=r0)
            else; call CartesianCoordinater(q,r,intdim,cartdim); end if
            call WilsonBMatrixAndInternalCoordinateq(r,B,qtemp,cartdim,intdim)
            forall(i=1:NStates,j=1:NStates); cartgrad(:,i,j)=matmul(transpose(B),intgrad(:,i,j)); end forall
        end subroutine Internal2Cartesian
    !=================== End ===================
!------------------- End --------------------

!--------------- Normal mode ----------------
	!Normal mode is the mass weighted eigenvector of Hessian:
	!    In Cartesian coordinate, it is the usual eigenvector
	!    In  internal coordinate, it is the generalized eigenvector of G under Hessian metric
	!G is built from mass and Wilson B matrix, for details see Wilson GF method in reference
    !For visualization, the output normal mode is not mass weighted

    !Obtain normal mode and vibrational frequency from Hessian in Cartesian coordinate
    !Note that normal modes are normalized by mass weight metric
    !Input:   H  : Cartesian coordinate Hessian (will be overwritten)
    !        mass: mass of each atom
    !Output: freq: vibrational angular frequencies (negative if imaginary)
    !        mode: normal modes contained in each column
    !Lowest 3 * NAtoms - vibdim modes are considered translation and rotation thus ruled out
    subroutine AnalyzeVibration(H,mass,freq,mode,NAtoms,vibdim)
        integer,intent(in)::NAtoms,vibdim
        real*8,dimension(3*NAtoms,3*NAtoms),intent(inout)::H
		real*8,dimension(NAtoms),intent(in)::mass
		real*8,dimension(vibdim),intent(out)::freq
		real*8,dimension(3*NAtoms,vibdim),intent(out)::mode
		integer::i; integer,dimension(3*NAtoms)::indice
		real*8,dimension(NAtoms)::sqrtmass; real*8,dimension(3*NAtoms)::freqall,freqabs
		!Obtain freq^2 and normal modes
		sqrtmass=dSqrt(mass)
		forall(i=1:NAtoms)
			H(:,3*i-2:3*i)=H(:,3*i-2:3*i)/sqrtmass(i)
			H(3*i-2:3*i,:)=H(3*i-2:3*i,:)/sqrtmass(i)
		end forall
        call My_dsyev('V',H,freqall,3*NAtoms)
        forall(i=1:3*NAtoms); H(3*i-2:3*i,:)=H(3*i-2:3*i,:)/sqrtmass(i); end forall
		!Rule out 3 * NAtoms - vibdim translations and rotations
		freqabs=dAbs(freqall)
		forall(i=1:3*NAtoms); indice(i)=i; end forall
		call dQuickSort(freqabs,1,3*NAtoms,indice,3*NAtoms)
		freqabs=freqall
		forall(i=1:vibdim)
			freqall(i)=freqabs(indice(i+3*NAtoms-vibdim))
			mode(:,i)=H(:,indice(i+3*NAtoms-vibdim))
        end forall
		do i=1,vibdim!freq^2 -> freq
			if(freqall(i)<0d0) then; freq(i)=-dSqrt(-freqall(i))
            else; freq(i)=dSqrt(freqall(i)); end if
		end do
		!Sort freq ascendingly, then sort normal modes accordingly
		forall(i=1:vibdim); indice(i)=i; end forall
		call dQuickSort(freq,1,vibdim,indice(1:vibdim),vibdim)
		H(:,1:vibdim)=mode; forall(i=1:vibdim); mode(:,i)=H(:,indice(i)); end forall
	end subroutine AnalyzeVibration

    !Use Wilson GF method to obtain normal mode and vibrational frequency from Hessian in internal coordinate
    !Note that normal modes are normalized by Hessian metric
    !Input:     H    : internal coordinate Hessian
    !           B    : Wilson B matrix
    !          mass  : mass of each atom
    !Output:   freq  : vibrational angular frequencies (negative if imaginary)
    !        intmode : internal coordinate normal modes contained in each column (Wilson L matrix)
    !          Linv  : Wilson L^-1 matrix
    !        cartmode: Cartesian coordinate normal modes contained in each column
    subroutine WilsonGFMethod(H,B,mass,freq,intmode,Linv,cartmode,intdim,NAtoms)
        integer,intent(in)::intdim,NAtoms
        real*8,dimension(intdim,intdim),intent(in)::H
        real*8,dimension(intdim,3*NAtoms),intent(in)::B
        real*8,dimension(NAtoms),intent(in)::mass
        real*8,dimension(intdim),intent(out)::freq
        real*8,dimension(intdim,intdim),intent(out)::intmode,Linv
        real*8,dimension(3*NAtoms,intdim),intent(out)::cartmode
        integer::i; integer,dimension(intdim)::indice
        real*8,dimension(intdim)::freqtemp!For imaginary frequency case
        real*8,dimension(intdim,intdim)::Htemp
        real*8,dimension(intdim,3*NAtoms)::Btemp
        !Try solving G . H . l = l . w^2 in generalized eigenvalue manner
        !LAPACK will normalized l by l(:,i)^T . H . l(:,j) = delta_ij,
        !but the true solution is L(:,i) . H . L(:,j) = w^2
        !This is why I call l raw normal mode. w is freq
        Htemp=H!Save Hessian, since it may still be needed if not positive definite
        forall(i=1:NAtoms); Btemp(:,3*i-2:3*i)=B(:,3*i-2:3*i)/mass(i); end forall
        intmode=matmul(Btemp,transpose(B)); call My_dsygv(2,'V',intmode,Htemp,freq,intdim,info=i)
        if(i==0) then!H is positive definite, freq^2 and raw normal modes are normally obtained: convert w^2 to w and l to L
            Linv=matmul(transpose(intmode),H)
            forall(i=1:intdim)
                freq(i)=dSqrt(freq(i))
                intmode(:,i)=intmode(:,i)*freq(i); Linv(i,:)=Linv(i,:)/freq(i)
            end forall
        else!resolve (G . H) . l = l . w^2 by general eigensolver then convert w^2 to w and l to L
            Htemp=matmul(intmode,H); call My_dgeev('V',Htemp,freq,freqtemp,intmode,intdim)
            forall(i=1:intdim)!Raw normal modes -> normal modes
                intmode(:,i)=intmode(:,i)*dSqrt(freq(i)/dot_product(intmode(:,i),matmul(H,intmode(:,i))))
            end forall
            do i=1,intdim!freq^2 -> freq
                if(freq(i)<0d0) then; freq(i)=-dSqrt(-freq(i))
                else; freq(i)=dSqrt(freq(i)); end if
            end do
            forall(i=1:intdim); indice(i)=i; end forall!Sort freq ascendingly, then sort normal modes accordingly
            call dQuickSort(freq,1,intdim,indice,intdim)
            Htemp=intmode; forall(i=1:intdim); intmode(:,i)=Htemp(:,indice(i)); end forall
            Linv=intmode; call My_dgetri(Linv,intdim)
        end if
        !Convert internal coordinate normal mode to Cartesian coordinate normal mode
        !L . dQ = dq = B . dr, where Q denotes internal coordinate normal mode
        !So L contains internal coordinate normal mode in each column
        !   B^g . L contains Cartesian coordinate normal mode in each column
        Btemp=B
        call dGeneralizedInverseTranspose(Btemp,intdim,3*NAtoms)
        cartmode=matmul(transpose(Btemp),intmode)
    end subroutine WilsonGFMethod
!------------------- End --------------------

end module GeometryTransformation