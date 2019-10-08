!Clustering routines
!
!For impatient:
!    For crazy amount of data, use K-means
!    For large but not insane, use GMM
!
!Instruction: (notation see corresponding subroutine)
!    K-means:
!        pros: low O(N dim K) cost
!        cons: Each cluster must be a sphere
!              Might fall into local minimum on large data set. In expectation, O(logK) worse than global minimum
!    Gaussian mixture model (GMM):
!        pros: O(N dim^2 K) cost, low for N
!        cons: High cost for dim
!              Each cluster must be an ellipse
!              Might fall into local maximum on large data set
module Clustering
    use General; use Mathematics; use LinearAlgebra
    implicit none

contains
!Input:  dim x N matrix data containing N dim dimensional data points
!        N order vector weight tells how many times each point repeats (equivalently the weight of each point)
!        K specifies the number of clusters
!        (optional) if initialized = true, use user input value of centre as initial guess
!                   else initialize centre by K-means++
!Output: dim x K matrix centre harvests K cluster centres
!        N order integer vector ascription telling the ascription of each data point
subroutine Kmeans(N,dim,data,weight,K,centre,ascription,initialized)
    integer,intent(in)::N,dim,K
    real*8,dimension(dim,N),intent(in)::data
    real*8,dimension(N),intent(in)::weight
    real*8,dimension(dim,K),intent(inout)::centre
    integer,dimension(N),intent(out)::ascription
    logical,intent(in),optional::initialized
    integer::i,j; integer,dimension(N)::ascriptionold
    real*8::mindsq,dbletemp; real*8,dimension(dim)::vectemp!Expectation work space
    real*8,dimension(K)::population!Maximization work space
    logical::geninit; integer::icentre; real*8,dimension(N)::dsq!K-means++ work space
    geninit=.true.; if(present(initialized)) geninit=.not.initialized
    if(geninit) then!Choose initial centre by K-means++
        centre(:,1)=data(:,ceiling(dble(N)*BetterRandomNumber()))!Randomly choose a point as 1st centre
        ascription=1
        do icentre=2,K!Find remaining K-1 centres
            dbletemp=0d0!Store sum of distance square
            do i=1,N
                vectemp=data(:,i)-centre(:,ascription(i)); dsq(i)=weight(i)*dot_product(vectemp,vectemp)
                dbletemp=dbletemp+dsq(i)
            end do
            dsq=dsq/dbletemp!Renormalize to be probability
            dbletemp=BetterRandomNumber()!Sum of distance square is no longer useful, use the variable as work space
            do i=1,N!Randomly determine the new centre
                if(dbletemp<dsq(i)) exit
                dbletemp=dbletemp-dsq(i)
            end do
            centre(:,icentre)=data(:,i)
            do i=1,N!Update ascription
                ascription(i)=1
                vectemp=data(:,i)-centre(:,1); mindsq=dot_product(vectemp,vectemp)
                do j=2,icentre!Only icentre centres have been decided
                    vectemp=data(:,i)-centre(:,j); dbletemp=dot_product(vectemp,vectemp)
                    if(dbletemp<mindsq) then
                        ascription(i)=j; mindsq=dbletemp
                    end if
                end do
            end do
        end do
    else!K-means++ will get intial ascription ready, so we also prepare it here
        do i=1,N
            ascription(i)=1
            vectemp=data(:,i)-centre(:,1); mindsq=dot_product(vectemp,vectemp)
            do j=2,K
                vectemp=data(:,i)-centre(:,j); dbletemp=dot_product(vectemp,vectemp)
                if(dbletemp<mindsq) then
                    ascription(i)=j; mindsq=dbletemp
                end if
            end do
        end do
    end if
    ascriptionold=ascription!Pre loop
    do!Main loop
        population=0d0; centre=0d0!Maximization: update centre
        do i=1,N
            j=ascription(i)
            population(j)=population(j)+weight(i); centre(:,j)=centre(:,j)+weight(i)*data(:,i)
        end do
        forall(i=1:K)
            centre(:,i)=centre(:,i)/population(i)
        end forall
        do i=1,N!Expectation: update ascription
            ascription(i)=1
            vectemp=data(:,i)-centre(:,1); mindsq=dot_product(vectemp,vectemp)
            do j=2,K
                vectemp=data(:,i)-centre(:,j); dbletemp=dot_product(vectemp,vectemp)
                if(dbletemp<mindsq) then
                    ascription(i)=j; mindsq=dbletemp
                end if
            end do
        end do
        do i=1,N!Converge when ascription no longer changes
            if(ascription(i)/=ascriptionold(i)) exit
        end do
        if(i>N) return
        ascriptionold(i:N)=ascription(i:N)!Get ready for next loop
    end do
end subroutine Kmeans

!Input:  dim x N matrix data containing N dim dimensional data points
!        N order vector weight tells how many times each point repeats (equivalently the weight of each point)
!        K specifies the number of clusters
!        Optional: initialized: (default = false) if true, use user input value of centre as initial guess
!                               else initialize centre by K-means
!                  Precision: (default = 1d-15) convergence considered when responsibility change < Precision
!Output: K order vector population harvests the population of each gaussian
!        dim x K matrix centre harvests K cluster centres
!        dim x dim x K 3rd-order tensor covariance harvests K cluster covariance matrices
!        N x K order matrix responsibility telling the responsibility of each data point
!Fractional weight is interpreted as infinite data set. On exit, weight is normalized
subroutine GaussianMixtureModel(N,dim,data,weight,K,population,centre,covariance,responsibility,initialized,Precision)
    integer,intent(in)::N,dim,K
    real*8,dimension(dim,N),intent(in)::data
    real*8,dimension(N),intent(inout)::weight
    real*8,dimension(K),intent(out)::population
    real*8,dimension(dim,K),intent(inout)::centre
    real*8,dimension(dim,dim,K),intent(out)::covariance
    real*8,dimension(N,K),intent(out)::responsibility
    logical,intent(in),optional::initialized
    real*8,intent(in),optional::Precision
    integer::icentre,i,j
    real*8::tol,dbletemp
    real*8,dimension(dim)::vectemp
    real*8,dimension(dim,dim)::matrixtemp
    real*8,dimension(N,K)::responsibilityold!Responsibility old
    !Initialization work space
    logical::geninit
    integer,dimension(N)::ascription
    real*8::mindsq
    if(present(Precision)) then
        tol=Precision
    else
        tol=1d-15
    end if
    geninit=.true.
    if(present(initialized)) geninit=.not.initialized
    if(geninit) then!Choose initial centre by K-means
        call Kmeans(N,dim,data,weight,K,centre,ascription)
    else!K-means will get intial ascription ready, so we also prepare it here
        do i=1,N
            ascription(i)=1
            vectemp=data(:,i)-centre(:,1); mindsq=dot_product(vectemp,vectemp)
            do j=2,K
                vectemp=data(:,i)-centre(:,j); dbletemp=dot_product(vectemp,vectemp)
                if(dbletemp<mindsq) then
                    ascription(i)=j; mindsq=dbletemp
                end if
            end do
        end do
    end if
    !Pre loop: centre has done, initialize responsibility & population & covariance based on ascription
        responsibilityold=0d0; forall(i=1:N); responsibilityold(i,ascription(i))=1d0; end forall
        population=0d0; forall(i=1:dim,j=1:dim,icentre=1:K,i>=j); covariance(i,j,icentre)=0d0; end forall
        do i=1,N; if(dble(floor(weight(i)))/=weight(i)) exit; end do
        if(i>N) then!All weight are integer, which indicates a finite data set
            do i=1,N
                j=ascription(i)
                population(j)=population(j)+weight(i)
                vectemp=data(:,i)-centre(:,j)
                covariance(:,:,j)=covariance(:,:,j)+weight(i)*vector_direct_square(vectemp,dim)
            end do
            forall(icentre=1:K)
                forall(i=1:dim)!Sample variance
                    covariance(i,i,icentre)=covariance(i,i,icentre)/(population(icentre)-1d0)
                end forall
                forall(i=1:dim,j=1:dim,i>j)!Sample covariance
                    covariance(i,j,icentre)=covariance(i,j,icentre)/population(icentre)
                end forall
            end forall
            dbletemp=sum(weight); population=population/dbletemp; weight=weight/dbletemp
        else!Infinite sample variance does not need that -1
            weight=weight/sum(weight)
            do i=1,N
                j=ascription(i)
                population(j)=population(j)+weight(i)
                vectemp=data(:,i)-centre(:,j)
                covariance(:,:,j)=covariance(:,:,j)+weight(i)*vector_direct_square(vectemp,dim)
            end do
            forall(i=1:dim,j=1:dim,icentre=1:K,i>=j)
                covariance(i,j,icentre)=covariance(i,j,icentre)/population(icentre)
            end forall
        end if
        dbletemp=sqrt2pi**N!Preparation for normal distribution calculation
    do!Main loop
        !Expectation
        do icentre=1,K!1st, prepare for normal distribution calculation
            forall(i=1:dim,j=1:dim,i>=j)
                matrixtemp(i,j)=covariance(i,j,icentre)
            end forall
            call syL2U(matrixtemp,dim)
            !Store the pre-exponential coefficient in population
            population(icentre)=population(icentre)/dbletemp/dSqrt(determinant(matrixtemp,dim))
            call My_dsytri(covariance(:,:,icentre),dim)!Invert covariance matrices
            call syL2U(covariance(:,:,icentre),dim)
        end do
        forall(i=1:K,j=1:N)!2nd, compute the probability density on each data point contributed by each gaussian
            responsibility(j,i)=population(i)&!Store in responsibility
                *dExp(-0.5d0*dot_product(data(:,j)-centre(:,i),matmul(covariance(:,:,i),data(:,j)-centre(:,i))))
        end forall
        forall(i=1:N)!Finally, update responsibility
            responsibility(i,:)=responsibility(i,:)/sum(responsibility(i,:))
        end forall
        if(My_dlange('M',responsibility-responsibilityold,N,K)<tol) exit!Check convergence
        !Maximization
        forall(i=1:N)!Prepare
            responsibilityold(i,:)=weight(i)*responsibility(i,:)
        end forall
        forall(i=1:K)!Update population
            population(i)=sum(responsibilityold(:,i))
        end forall
        centre=matmul(data,responsibilityold)!Update centre & covariance
        do icentre=1,K
            forall(i=1:dim,j=1:dim,i>=j)
                covariance(i,j,icentre)=0d0
            end forall
            do j=1,N
                vectemp=data(:,j)-centre(:,icentre)
                covariance(:,:,icentre)=vector_direct_square(vectemp,dim)*responsibilityold(j,icentre)
            end do
        end do
        forall(icentre=1:K)
            centre(:,icentre)=centre(:,icentre)/population(icentre)
            forall(i=1:dim,j=1:dim,i>=j)
                covariance(i,j,icentre)=covariance(i,j,icentre)/population(icentre)
            end forall
        end forall
        responsibilityold=responsibility!Get ready for next loop
    end do
    !Exit after expectation, do a maximization accordingly
    forall(i=1:N)!Prepare
        responsibilityold(i,:)=weight(i)*responsibility(i,:)
    end forall
    forall(i=1:K)!Update population
        population(i)=sum(responsibilityold(:,i))
    end forall
    centre=matmul(data,responsibilityold)!Update centre & covariance
    do icentre=1,K
        forall(i=1:dim,j=1:dim,i>=j)
            covariance(i,j,icentre)=0d0
        end forall
        do j=1,N
            vectemp=data(:,j)-centre(:,icentre)
            covariance(:,:,icentre)=vector_direct_square(vectemp,dim)*responsibilityold(j,icentre)
        end do
    end do
    forall(icentre=1:K)
        centre(:,icentre)=centre(:,icentre)/population(icentre)
        forall(i=1:dim,j=1:dim,i>=j)
            covariance(i,j,icentre)=covariance(i,j,icentre)/population(icentre)
        end forall
    end forall
end subroutine GaussianMixtureModel

end module Clustering