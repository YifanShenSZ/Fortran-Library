!Mathematical and physical constant, mathematical routine
module Mathematics
    implicit none

!Constant
    !Mathematical constant
        real*8,parameter::pi=3.141592653589793d0,enature=2.718281828459045d0,&
            GoldenRatio=0.6180339887498949d0,EulerGamma=0.57721566490153286060651209d0
        complex*16,parameter::ci=(0d0,1d0)
        !Frequently used form (m = multiply, d = divide, inv = inverse)
        real*8,parameter::Sqrt2=1.4142135623730951d0,Sqrt3=1.7320508075688772d0,invSqrt2=0.7071067811865475d0,&
            pim2=6.283185307179586d0,pim4=12.566370614359172d0,&
            pid2=1.5707963267948966d0,pid4=0.7853981633974483d0,pid8=0.39269908169872414d0,&
            pisqd3=3.2898681336964524d0,&
            sqrtpi=1.7724538509055159d0,sqrt2pi=2.506628274631d0
    !Physical constant (If has unit, then in atomic unit)
        real*8,parameter::NAvogadro=6.02214076d23
    !Unit conversion (Multiplying xIny converts x to y)
        real*8,parameter::DegInRad=0.017453292519943295,&!degree to radius
            !Common unit to atomic unit
            AMUInAU=1822.888486192d0,&!Mass
            AInAU=1.8897261339212517d0,&!Length
            kJmolInAU=0.00038087967507991464d0,cm_1InAu=4.556335830019422d-6,&!Energy
            fsInAU=41.341373336561354d0,&!Time
            KInAU=3.166813539739535d-6,&!Temperature
            barInAU=3.39882737736419d-9!Pressure

!Overload
    interface RK4
        module procedure dRK4_t, dRK4, zRK4
    end interface RK4

contains
!----------------- Special function -----------------
    real*8 function Gaussian(x, miu, sigma)
        real*8,intent(in)::x,miu,sigma
        real*8::temp
        temp=(x-miu)/sigma; Gaussian=exp(-temp*temp/2d0)/(sqrt2pi*sigma)
    end function Gaussian
    real*8 function dGaussiandmiu(x, miu, sigma)
        real*8,intent(in)::x,miu,sigma
        real*8::temp
        temp=(x-miu)/sigma
        dGaussiandmiu=2d0/(sqrt2pi*sigma*sigma)*temp*exp(-temp*temp/2d0)
    end function dGaussiandmiu
    real*8 function dGaussiandsigma(x, miu, sigma)
        real*8,intent(in)::x,miu,sigma
        real*8::temp
        temp=(x-miu)/sigma; temp=temp*temp
        dGaussiandsigma=2d0/(sqrt2pi*sigma*sigma)*(temp-1d0)*exp(-temp/2d0)
    end function dGaussiandsigma

    real*8 function Lorentzian(x, miu, sigma)
        real*8,intent(in)::x,miu,sigma
        real*8::temp
        temp=(x-miu)/sigma; Lorentzian=1d0/(pi*sigma)/(1d0+temp*temp)
    end function Lorentzian
    real*8 function dLorentziandmiu(x, miu, sigma)
        real*8,intent(in)::x,miu,sigma
        real*8::temp,temp1
        temp=(x-miu)/sigma; temp1=1d0+temp*temp
        dLorentziandmiu=2d0/(pi*sigma*sigma)*temp/(temp1*temp1)
    end function dLorentziandmiu
    real*8 function dLorentziandsigma(x, miu, sigma)
        real*8,intent(in)::x,miu,sigma
        real*8::temp,temp1,temp2
        temp=(x-miu)/sigma; temp=temp*temp
        temp1=1d0+temp
        temp2=pi*sigma*sigma
        dLorentziandsigma=1d0/temp2/temp1*(2d0*temp/temp1-1d0)
    end function dLorentziandsigma

    real*8 function inverse_erfc(x)
        real*8::x
        logical::flag; real*8::temp
        flag=.false.
        if(x>1d0) then
            x=2d0-x; flag=.true.
        end if
        if(x<0.1d0) then
            temp=-0.4515827052894548d0-2d0*log(x)
            inverse_erfc=Sqrt((temp-log(temp))/2d0)
        else if(x<0.25d0) then
            temp=4.995264535887506d0*(x-0.15)
            inverse_erfc=1.0179024648320276d0-temp/2d0+temp**2*0.2544756162080069d0&
                -temp**3*0.21435423798518619d0+temp**4*0.20605638309556853d0
        else 
            temp=(1d0-x)*1.7724538509055159d0
            inverse_erfc=temp/2d0+temp**3/24d0+temp**5*7d0/960d0+temp**7*127d0/80640d0&
                +temp**9*4369d0/11612160d0
        end if
        if(flag) then
            x=2d0-x; inverse_erfc=-inverse_erfc
        end if
    end function inverse_erfc

    !========== Gaussian integral ==========
        !Integrate[ 1 / Sqrt(2pi) / sigma * Exp(-0.5 * (x / sigma)^2) * x^i, {x, -Infinity, Infinity}]
        real*8 function GaussianIntegral(i, sigma)
            integer,intent(in)::i
            real*8,intent(in)::sigma
            GaussianIntegral=0d0
            if(mod(i,2)==0) then
                GaussianIntegral=dFactorial2(i-1)*sigma**i
            end if
        end function GaussianIntegral

        !Divide sigma^i to make GaussianIntegral dimensionless
        real*8 function GaussianIntegraldsig(i)
            integer,intent(in)::i
            GaussianIntegraldsig=0d0
            if(mod(i,2)==0) then
                GaussianIntegraldsig=dFactorial2(i-1)
            end if
        end function GaussianIntegraldsig

        !Extend to two variable x & p case x^i * p^j
        real*8 function BinaryGaussianIntegral(i, j, sigmax, sigmap, rho)
            integer,intent(in)::i,j
            real*8,intent(in)::sigmax,sigmap,rho
            integer::k,minimum,maximum
            BinaryGaussianIntegral=0d0
            if(mod(i+j,2)==0) then
                minimum=min(i,j)
                maximum=max(i,j)
                do k=0,minimum/2
                    BinaryGaussianIntegral=BinaryGaussianIntegral&
                        +rho**(minimum-2*k)*dCombination(minimum,2*k)*dFactorial2(2*k-1)*dPermutation(maximum,minimum-2*k)*dFactorial2(maximum-minimum+2*k-1)
                end do
                BinaryGaussianIntegral=BinaryGaussianIntegral*sigmax**i*sigmap**j
            end if
        end function BinaryGaussianIntegral

        real*8 function BinaryGaussianIntegraldsig(i, j, rho)
            integer,intent(in)::i,j
            real*8,intent(in)::rho
            integer::k,minimum,maximum
            BinaryGaussianIntegraldsig=0d0
            if(mod(i+j,2)==0) then
                minimum=min(i,j)
                maximum=max(i,j)
                do k=0,minimum/2
                    BinaryGaussianIntegraldsig=BinaryGaussianIntegraldsig&
                        +rho**(minimum-2*k)*dCombination(minimum,2*k)*dFactorial2(2*k-1)*dPermutation(maximum,minimum-2*k)*dFactorial2(maximum-minimum+2*k-1)
                end do
            end if
        end function BinaryGaussianIntegraldsig
    !================= End =================

    !=========== Gamma function ============
        real*8 function lnGamma(x)
            real*8,intent(in)::x
            real*8::x1,x2,y
            real*8,dimension(5)::r4=[0.279195317918525d0,0.4917317610505968d0,&
                0.0692910599291889d0,3.350343815022304d0,6.012459259764103d0]
            real*8,dimension(9)::r1=[-2.66685511495d0,-24.4387534237d0,-21.9698958928d0, &
                11.1667541262d0,3.13060547623d0,0.607771387771d0,&
                11.9400905721d0,31.4690115749d0,15.2346874070d0],&
                r2=[-78.3359299449d0,-142.046296688d0,137.519416416d0,&
                78.6994924154d0,4.16438922228d0,47.0668766060d0,&
                313.399215894d0,263.505074721d0,43.3400022514d0],&
                r3=[-2.12159572323d5,2.30661510616d5,2.74647644705d4,&
                -4.02621119975d4,-2.29660729780d3,-1.16328495004d5, &
                -1.46025937511d5,-2.42357409629d4,-5.70691009324d2]
            lngamma=0d0
            if (x<1.5d0) then
                if (x<0.5d0) then
                    lngamma=-log(x)
                    y=x+1d0
                else
                    y=x
                    x1=x-1d0
                end if
                lngamma=lngamma+x1*((((r1(5)*y+r1(4))*y+r1(3))*y+r1(2))*y+r1(1))/((((y+r1(9))*y+r1(8))*y+r1(7))*y+r1(6))
            else if(x<4d0) then
                y=x-2d0
                lngamma=y*((((r2(5)*x+r2(4))*x+r2(3))*x+r2(2))*x+r2(1))/((((x+r2(9))*x+r2(8))*x+r2(7))*x+r2(6))
            else if(x<12d0) then
                lngamma = ((((r3(5)*x+r3(4))*x+r3(3))*x+r3(2))*x+r3(1))/((((x+r3(9))*x+r3(8))*x+r3(7))*x+r3(6))
            else
                y=log(x)
                lngamma=x*(y-1d0)-0.5d0*y+9.18938533204673d-1
                if (x<5.1d5) then
                    x1=1.0d0/x
                    x2=x1**2
                    lngamma=lngamma+x1*((r4(3)*x2+r4(2))*x2+r4(1))/((x2+r4(5))*x2+r4(4))
                end if
            end if
        end function lnGamma
    
        !Regularized incomplete gamma function
        real*8 function gamma_regularized_inc_lower(p, x)
            real*8,intent(in)::p,x
            integer::i
            real*8::a,an,arg,b,dif,factor,g,gin,rn,term
            real*8,dimension(6)::pn
            g=lngamma(p)
            arg=p*log(x)-x-g
            factor=exp(arg)
            if(x<1d0.or.x<p) then
                gin=1d0
                term=1d0
                rn=p
                do while(term>1d-15)
                    rn=rn+1d0
                    term=term*x/rn
                    gin=gin+term
                end do
                gamma_regularized_inc_lower=gin*factor/p
            else
                a=1d0-p
                b=a+x+1d0
                term=0d0
                pn(1)=1d0
                pn(2)=x
                pn(3)=x+1d0
                pn(4)=x*b
                gin=pn(3)/pn(4)
                do
                    a=a+1d0
                    b=b+2d0
                    term=term+1d0
                    an=a*term
                    forall(i=1:2)
                        pn(i+4)=b*pn(i+2)-an*pn(i)
                    end forall
                    if(pn(6)/=0d0) then
                        rn=pn(5)/pn(6)
                        dif=abs(gin-rn)
                        if(dif<1d-15.and.dif<1d-15*rn) then
                            gamma_regularized_inc_lower = 1.0d0 - factor * gin
                            exit
                        end if
                        gin=rn
                    end if
                    do i=1,4
                        pn(i)=pn(i+2)
                    end do
                    if (1d37<abs(pn(5))) then
                        forall(i=1:4)
                            pn(i)=pn(i)/1d37
                        end forall
                    end if
                end do
            end if
        end function gamma_regularized_inc_lower

        !Incomplete gamma function
        real*8 function gamma_inc(p, x)
            real*8 ,intent(in)::p,x
            real*8::temp,temp2,temp3,temp6
            if(p<1d-15) then!10 significant figures
                if(x>18.78741462316355d0) then
                    gamma_inc=0d0
                else if(x>13.29d0) then
                    temp=x-15.99d0
                    temp2=temp*temp
                    temp3=temp*temp2
                    temp6=temp3*temp3
                    gamma_inc=-eulergamma-log(x)&
                        +3.3491791984585926d0+temp*0.06253907982075221d0-temp2*(0.0019555649204058842d0+temp2*3.823889781974916d-6)&
                        +temp3*(8.153164493904012d-5+temp2*1.9125482223137583d-7)&
                        +temp6*(-9.95754526846618d-9+temp*5.323631012903137d-10-temp2*2.8955510452797626d-11&
                        +temp3*1.5900561610283544d-12-temp*temp3*8.753766072858468d-14)
                else if(x>9.29d0) then
                    temp=x-11.26d0
                    temp2=temp*temp
                    temp3=temp*temp2
                    temp6=temp3*temp3
                    gamma_inc=-eulergamma-log(x)&
                        +2.998473344163554d0+temp*0.08880880303082772d0-temp2*(0.003942980690850256d0+temp2*1.548920196156643d-5)&
                        +temp3*(0.00023325998949853596d0+temp2*1.0909454673108656d-6)&
                        +temp6*(-7.915055845502563d-8+temp*5.798241768265512d-9-temp2*4.222086908296498d-10&
                        +temp3*3.01783879825018d-11-temp*temp3*2.096958494268975d-12)
                else if(x>6.24d0) then
                    temp=x-7.74d0
                    temp2=temp*temp
                    temp3=temp*temp2
                    temp6=temp3*temp3
                    gamma_inc=-eulergamma-log(x)&
                        +2.623667685122219d0+temp*0.12914275561045493d0-temp2*(0.008314449873086318+temp2*6.614403525922297d-5)&
                        +temp3*(7.067770869351106d-4+temp2*6.368169476536038d-6)&
                        +temp6*(-6.075635427151326d-7+temp*5.6129848279670643d-8-temp2*4.951311575317511d-9&
                        +temp3*4.1372418821110844d-10-temp*temp3*3.2617274753582414d-11)
                else if(x>3.82d0) then
                    temp=x-5.01d0
                    temp2=temp*temp
                    temp3=temp*temp2
                    temp6=temp3*temp3
                    gamma_inc=-eulergamma-log(x)&
                        +2.1897864802195994d0+temp*0.19826928077719455d0-temp2*(0.01912159455797797d0+temp2*0.00029220534183235477d0)&
                        +temp3*(0.002322537422676534d0+temp2*3.5563555405274084d-5)&
                        +temp6*(-4.066098341347632d-6+temp*4.3146411654336176d-7-temp2*4.2331758731517724d-8&
                        +temp3*3.841319120017832d-9-temp*temp3*3.231267642291151d-10)
                else if(x>1.85d0) then
                    temp=x-2.81d0
                    temp2=temp*temp
                    temp3=temp*temp2
                    temp6=temp3*temp3
                    gamma_inc=-eulergamma-log(x)&
                        +1.6270397285101232d0+temp*0.33444662192442226d0-temp2*(0.048797442977232874d0+temp2*0.0012441706978081175)&
                        +temp3*(0.008006214680697388d0+temp2*0.00017566843000357574d0)&
                        +temp6*(-2.233890154694154d-5+temp*2.5630586828795826d-6-temp2*2.667249013198523d-7&
                        +temp3*2.533095435935187d-8-temp*temp3*2.2088886053649257d-9)
                else if(x>0.723d0) then
                    temp=x-1.25d0
                    temp2=temp*temp
                    temp3=temp*temp2
                    temp6=temp3*temp3
                    gamma_inc=-eulergamma-log(x)&
                        +0.9467725887416528d0+temp*0.570796162511848d0-temp2*(0.11371654626066321d0+temp2*0.003918751159257355d0)&
                        +temp3*(0.022448185090995054d0+temp2*0.000597968762856782d0)&
                        +temp6*(-8.030717872653526d-5+temp*9.590827815628838d-6-temp2*1.0289604856217327d-6&
                        +temp3*1.0008090251673041d-7-temp*temp3*8.895816642345844d-9)
                else
                    temp2=x*x
                    temp3=x*temp2
                    temp6=temp3*temp3
                    gamma_inc=-eulergamma-log(x)&
                        +x-temp2*(0.25d0+temp2/96d0)+temp3*(0.05555555555555555d0+temp2/600d0)&
                        +temp6*(-2.314814814814815d-4+x/35280d0-temp2/322560d0+temp3/3265920d0-x*temp3/36288d3)
                end if
            else
                gamma_inc=(1d0-gamma_regularized_inc_lower(p,x))*gamma(p)
            end if
        end function gamma_inc
    !================= End =================
!----------------------- End ------------------------

!------------------ Combinatorics -------------------
    !Exact factorial for 0 <= N <= 20, double precision for 21 <= N <= 40, 8 significant figures for N >= 41
    !Optional: Warning: (default = true) if true, throw a warning when using 8 significant figures approximation
    real*8 function dFactorial(N, Warning)
        integer,intent(in)::N
        logical,intent(in),optional::Warning
        integer::i; real*8::temp
        select case(N)
            case(0); dFactorial=1d0
            case(1); dFactorial=1d0
            case(2); dFactorial=2d0
            case(3); dFactorial=6d0
            case(4); dFactorial=24d0
            case(5); dFactorial=120d0
            case(6); dFactorial=720d0
            case(7); dFactorial=5040d0
            case(8); dFactorial=40320d0
            case(9); dFactorial=362880d0
            case(10); dFactorial=3628800d0
            case(11); dFactorial=39916800d0
            case(12); dFactorial=479001600d0
            case(13); dFactorial=6227020800d0
            case(14); dFactorial=87178291200d0
            case(15); dFactorial=1307674368d3
            case(16); dFactorial=20922789888d3
            case(17); dFactorial=355687428096d3
            case(18); dFactorial=6402373705728d3
            case(19); dFactorial=121645100408832d3
            case(20); dFactorial=243290200817664d4
            case(21); dFactorial=5109094217170944d4
            case(22); dFactorial=112400072777760768d4
            case(23); dFactorial=2585201673888497664d4
            case(24); dFactorial=62044840173323943936d4
            case(25); dFactorial=15511210043330985984d6
            case(26); dFactorial=403291461126605635584d6
            case(27); dFactorial=10888869450418352160768d6
            case(28); dFactorial=304888344611713860501504d6
            case(29); dFactorial=8841761993739701954543616d6
            case(30); dFactorial=26525285981219105863630848d7
            case(31); dFactorial=822283865417792281772556288d7
            case(32); dFactorial=26313083693369353016721801216d7
            case(33); dFactorial=868331761881188649551819440128d7
            case(34); dFactorial=29523279903960414084761860964352d7
            case(35); dFactorial=103331479663861449296666513375232d8
            case(36); dFactorial=3719933267899012174679994481508352d8
            case(37); dFactorial=137637530912263450463159795815809024d8
            case(38); dFactorial=5230226174666011117600072241000742912d8
            case(39); dFactorial=203978820811974433586402817399028973568d8
            case(40); dFactorial=815915283247897734345611269596115894272d9
            case default
                if((.not.present(Warning)).or.Warning) write(*,'(1x,A22)')'Not accurate factorial'
                temp = (9d0 * N)**3.141592653589793d0
                dFactorial = Sqrt(6.283185307179586d0 * N) &
                           * (N / enature)**N &
                           * Exp(1d0 / 12d0 / N - Log(9d0 * N) / (temp - 1d0 / temp))
        end select
    end function dFactorial

    !Exact double factorial for -1 <= N <= 30, double precision for 31 <= N <= 40, depend on factorial for N >= 41
    real*8 function dFactorial2(N)
        integer,intent(in)::N
        integer::i
        select case(N)
            case(-1); dFactorial2=1d0
            case(0); dFactorial2=1d0
            case(1); dFactorial2=1d0
            case(2); dFactorial2=2d0
            case(3); dFactorial2=3d0
            case(4); dFactorial2=8d0
            case(5); dFactorial2=15d0
            case(6); dFactorial2=48d0
            case(7); dFactorial2=105d0
            case(8); dFactorial2=384d0
            case(9); dFactorial2=945d0
            case(10); dFactorial2=3840d0
            case(11); dFactorial2=10395d0
            case(12); dFactorial2=46080d0
            case(13); dFactorial2=135135d0
            case(14); dFactorial2=645120d0
            case(15); dFactorial2=2027025d0
            case(16); dFactorial2=10321920d0
            case(17); dFactorial2=34459425d0
            case(18); dFactorial2=185794560d0
            case(19); dFactorial2=654729075d0
            case(20); dFactorial2=3715891200d0
            case(21); dFactorial2=13749310575d0
            case(22); dFactorial2=81749606400d0
            case(23); dFactorial2=316234143225d0
            case(24); dFactorial2=1961990553600d0
            case(25); dFactorial2=7905853580625d0
            case(26); dFactorial2=51011754393600d0
            case(27); dFactorial2=213458046676875d0
            case(28); dFactorial2=14283291230208d2
            case(29); dFactorial2=6190283353629375d0
            case(30); dFactorial2=42849873690624d3
            case(31); dFactorial2=191898783962510625d0
            case(32); dFactorial2=1371195958099968000d0
            case(33); dFactorial2=6332659870762850625d0
            case(34); dFactorial2=46620662575398912000d0
            case(35); dFactorial2=221643095476699771875d0
            case(36); dFactorial2=1678343852714360832000d0
            case(37); dFactorial2=8200794532637891559375d0
            case(38); dFactorial2=63777066403145711616000d0
            case(39); dFactorial2=319830986772877770815625d0
            case(40); dFactorial2=2551082656125828464640000d0
            case default
                if(mod(n,2)==1) then
                    dFactorial2 = dFactorial(N + 1) / (2**((N + 1) / 2) * dFactorial((N + 1) / 2))
                else
                    dFactorial2 = 2**(N / 2) * dFactorial(N / 2)
                end if
        end select
    end function dFactorial2

    !Exact permutation for 0 <= M <= 10 or N < 10, depend on factorial otherwise
    real*8 function dPermutation(M, N)
        integer,intent(in)::M,N
        integer::i
        if(M<2.or.N==0) then
            dPermutation=1d0
        else if(N==1) then
            dPermutation=M
        else if(N==M.or.N==M-1) then
            dPermutation=dFactorial(M)
        else
            select case(M)
                case(4); dPermutation=12d0
                case(5)
                    select case(N)
                        case(2); dPermutation=20d0
                        case(3); dPermutation=60d0
                    end select
                case(6)
                    select case(N)
                        case(2); dPermutation=30d0
                        case(3); dPermutation=120d0
                        case(4); dPermutation=360d0
                    end select
                case(7)
                    select case(N)
                        case(2); dPermutation=42d0
                        case(3); dPermutation=210d0
                        case(4); dPermutation=840d0
                        case(5); dPermutation=2520d0
                    end select
                case(8)
                    select case(N)
                        case(2); dPermutation=56d0
                        case(3); dPermutation=336d0
                        case(4); dPermutation=1680d0
                        case(5); dPermutation=6720d0
                        case(6); dPermutation=20160d0
                    end select
                case(9)
                    select case(N)
                        case(2); dPermutation=72d0
                        case(3); dPermutation=504d0
                        case(4); dPermutation=3024d0
                        case(5); dPermutation=15120d0
                        case(6); dPermutation=60480d0
                        case(7); dPermutation=181440d0
                    end select
                case(10)
                    select case(N)
                        case(2); dPermutation=90d0
                        case(3); dPermutation=720d0
                        case(4); dPermutation=5040d0
                        case(5); dPermutation=30240d0
                        case(6); dPermutation=151200d0
                        case(7); dPermutation=604800d0
                        case(8); dPermutation=1814400d0
                    end select
                case default
                    if(M<20.and.N<10) then
                        dPermutation=M*(M-1)
                        do i=M-2,M-N+1,-1
                            dPermutation=dPermutation*i
                        end do
                    else
                        dPermutation=dFactorial(M)/dFactorial(M-N)
                    end if
            end select
        end if
    end function dPermutation

    !Exact combination for 0 <= M <= 11 or (N = 0 or N = 1 or N = M - 1 or N = M), depend on permutation and factorial otherwise
    real*8 function dCombination(M, N)
        integer,intent(in)::M,N
        integer::ntemp
        if(M<2.or.N==0.or.N==M) then
            dCombination=1d0
        else if(N==1.or.N==(M-1)) then
            dCombination=M
        else 
            select case(M)
                case(4); dCombination=6d0
                case(5); dCombination=10d0
                case default
                    if (N < M / 2) then
                        ntemp=M-N
                    else
                        ntemp=N
                    end if
                    select case(M)
                        case(6)
                            select case(ntemp)
                                case(4); dCombination=15d0
                                case(3); dCombination=20d0
                            end select
                        case(7)
                            select case(ntemp)
                                case(5); dCombination=21d0
                                case(4); dCombination=35d0
                            end select
                        case(8)
                            select case(ntemp)
                                case(6); dCombination=28d0
                                case(5); dCombination=56d0
                                case(4); dCombination=70d0
                            end select
                        case(9)
                            select case(ntemp)
                                case(7); dCombination=36d0
                                case(6); dCombination=84d0
                                case(5); dCombination=126d0
                            end select
                        case(10)
                            select case(ntemp)
                                case(8); dCombination=45d0
                                case(7); dCombination=120d0
                                case(6); dCombination=210d0
                                case(5); dCombination=252d0
                            end select
                        case(11)
                            select case(ntemp)
                                case(9); dCombination=55d0
                                case(8); dCombination=165d0
                                case(7); dCombination=330d0
                                case(6); dCombination=462d0
                            end select
                        case default
                            dCombination=dPermutation(M,N)/dFactorial(N)
                    end select
            end select
        end if
    end function dCombination

    !Exact factorial for N <= 20, 8 bytes integer cannot represent N >= 21
    integer*8 function iFactorial(N, Warning)
        integer,intent(in)::N
        logical,intent(in),optional::Warning
        select case(N)
            case(0); iFactorial=1 
            case(1); iFactorial=1 
            case(2); iFactorial=2 
            case(3); iFactorial=6 
            case(4); iFactorial=24 
            case(5); iFactorial=120 
            case(6); iFactorial=720 
            case(7); iFactorial=5040 
            case(8); iFactorial=40320 
            case(9); iFactorial=362880 
            case(10); iFactorial=3628800 
            case(11); iFactorial=39916800 
            case(12); iFactorial=479001600 
            case(13); iFactorial=6227020800 
            case(14); iFactorial=87178291200 
            case(15); iFactorial=1307674368000
            case(16); iFactorial=20922789888000
            case(17); iFactorial=355687428096000
            case(18); iFactorial=6402373705728000
            case(19); iFactorial=121645100408832000
            case(20); iFactorial=2432902008176640000
            case default; write(*,*)'Failed integer factorial: 8 bytes integer upper limit exceeded!'
        end select
    end function iFactorial

    !Exact double factorial for N <= 33, 8 bytes integer cannot represent N >= 34
    integer*8 function iFactorial2(N)
        integer,intent(in)::N
        select case(N)
            case(-1); iFactorial2=1 
            case(0); iFactorial2=1 
            case(1); iFactorial2=1 
            case(2); iFactorial2=2 
            case(3); iFactorial2=3 
            case(4); iFactorial2=8 
            case(5); iFactorial2=15 
            case(6); iFactorial2=48 
            case(7); iFactorial2=105 
            case(8); iFactorial2=384 
            case(9); iFactorial2=945 
            case(10); iFactorial2=3840 
            case(11); iFactorial2=10395 
            case(12); iFactorial2=46080 
            case(13); iFactorial2=135135 
            case(14); iFactorial2=645120 
            case(15); iFactorial2=2027025 
            case(16); iFactorial2=10321920 
            case(17); iFactorial2=34459425 
            case(18); iFactorial2=185794560 
            case(19); iFactorial2=654729075 
            case(20); iFactorial2=3715891200 
            case(21); iFactorial2=13749310575 
            case(22); iFactorial2=81749606400 
            case(23); iFactorial2=316234143225 
            case(24); iFactorial2=1961990553600 
            case(25); iFactorial2=7905853580625 
            case(26); iFactorial2=51011754393600 
            case(27); iFactorial2=213458046676875 
            case(28); iFactorial2=1428329123020800 
            case(29); iFactorial2=6190283353629375 
            case(30); iFactorial2=42849873690624000 
            case(31); iFactorial2=191898783962510625 
            case(32); iFactorial2=1371195958099968000 
            case(33); iFactorial2=6332659870762850625 
            case default
                write(*,'(1x,A69)')'Failed integer double factorial: 8 bytes integer upper limit exceeded!'
        end select
    end function iFactorial2

    !Exact permutation for 0 <= M <= 10 or N < 10, depend on factorial otherwise
    integer*8 function iPermutation(M, N)
        integer,intent(in)::M,N
        integer::i
        if(M<2.or.N==0) then
            iPermutation=1 
        else if(N==1) then
            iPermutation=M
        else if(N==M.or.N==M-1) then
            iPermutation=iFactorial(M)
        else
            select case(M)
                case(4); iPermutation=12 
                case(5)
                    select case(N)
                        case(2); iPermutation=20 
                        case(3); iPermutation=60 
                    end select
                case(6)
                    select case(N)
                        case(2); iPermutation=30 
                        case(3); iPermutation=120 
                        case(4); iPermutation=360 
                    end select
                case(7)
                    select case(N)
                        case(2); iPermutation=42 
                        case(3); iPermutation=210 
                        case(4); iPermutation=840 
                        case(5); iPermutation=2520 
                    end select
                case(8)
                    select case(N)
                        case(2); iPermutation=56 
                        case(3); iPermutation=336 
                        case(4); iPermutation=1680 
                        case(5); iPermutation=6720 
                        case(6); iPermutation=20160 
                    end select
                case(9)
                    select case(N)
                        case(2); iPermutation=72 
                        case(3); iPermutation=504 
                        case(4); iPermutation=3024 
                        case(5); iPermutation=15120 
                        case(6); iPermutation=60480 
                        case(7); iPermutation=181440 
                    end select
                case(10)
                    select case(N)
                        case(2); iPermutation=90 
                        case(3); iPermutation=720 
                        case(4); iPermutation=5040 
                        case(5); iPermutation=30240 
                        case(6); iPermutation=151200 
                        case(7); iPermutation=604800 
                        case(8); iPermutation=1814400 
                    end select
                case default
                    if(M<20.and.N<10) then
                        iPermutation=M*(M-1)
                        do i=M-2,M-N+1,-1
                            iPermutation=iPermutation*i
                        end do
                    else
                        iPermutation=iFactorial(M)/iFactorial(M-N)
                    end if
            end select
        end if
    end function iPermutation

    !Exact combination for 0 <= M <= 11 or (N = 0 or N = 1 or N = M - 1 or N = M), depend on permutation and factorial otherwise
    integer*8 function iCombination(M, N)
        integer,intent(in)::M,N
        integer::ntemp
        if(M<2.or.N==0.or.N==M) then
            iCombination=1 
        else if(N==1.or.N==(M-1)) then
            iCombination=M
        else 
            select case(M)
                case(4)
                    iCombination=6 
                case(5)
                    iCombination=10 
                case default
                    if (N < M / 2) then
                        ntemp = M - N
                    else
                        ntemp=N
                    end if
                    select case(M)
                        case(6)
                            select case(ntemp)
                                case(4); iCombination=15 
                                case(3); iCombination=20 
                            end select
                        case(7)
                            select case(ntemp)
                                case(5); iCombination=21 
                                case(4); iCombination=35 
                            end select
                        case(8)
                            select case(ntemp)
                                case(6); iCombination=28 
                                case(5); iCombination=56 
                                case(4); iCombination=70 
                            end select
                        case(9)
                            select case(ntemp)
                                case(7); iCombination=36 
                                case(6); iCombination=84 
                                case(5); iCombination=126 
                            end select
                        case(10)
                            select case(ntemp)
                                case(8); iCombination=45 
                                case(7); iCombination=120 
                                case(6); iCombination=210 
                                case(5); iCombination=252 
                            end select
                        case(11)
                            select case(ntemp)
                                case(9); iCombination=55 
                                case(8); iCombination=165 
                                case(7); iCombination=330 
                                case(6); iCombination=462 
                            end select
                        case default
                            iCombination=iPermutation(M,N)/iFactorial(N)
                    end select
            end select
        end if
    end function iCombination
!----------------------- End ------------------------

!-------------------- Quaternion --------------------
    !Here a quaternion is represented by a 4 dimensional vector

    !Quaternion multiplication a * b
    function quamul(a, b)
        real*8,dimension(4),intent(in)::a,b
        real*8,dimension(4)::quamul
        quamul(1)=a(1)*b(1)-a(2)*b(2)-a(3)*b(3)-a(4)*b(4)
        quamul(2)=a(2)*b(1)+a(1)*b(2)+a(4)*b(3)-a(3)*b(4)
        quamul(3)=a(3)*b(1)+a(1)*b(3)+a(2)*b(4)-a(4)*b(2)
        quamul(4)=a(4)*b(1)+a(1)*b(4)+a(3)*b(2)-a(2)*b(3)
    end function quamul
    
    !Rotate each column of 3 x NAtoms matrix r by unit quaternion q
    !q = [ cos(alpha/2), sin(alpha/2) * axis ]
    subroutine Rotate(q, r, NAtoms)
        real*8,dimension(4),intent(in)::q
        integer,intent(in)::NAtoms
        real*8,dimension(3,NAtoms),intent(inout)::r
        integer::i; real*8,dimension(4)::qstar,qtemp,qresult
        qstar(1)=q(1); qstar(2:4)=-q(2:4); qtemp(1)=0d0
        do i=1,NAtoms
            qtemp(2:4)=r(:,i)
            qresult=quamul(quamul(qstar,qtemp),q)
            r(:,i)=qresult(2:4)
        end do
    end subroutine Rotate
!----------------------- End ------------------------

!---------- Ordinary differential equation ----------
    !The ordinary differential equation is provided by a subroutine to calculate du/dt(t, u)
    !From tau and u(tau), calculate u(tau + dt)

    !Input:  tau, old = u(tau), f, dt, dim
    !Output: new harvests u(dt)
    !For du/dt(t, u), subroutine f(du/dt, t, u, dim)
    !For du/dt(u), subroutine f(du/dt, u, dim)
    !where du/dt and u are dim dimentional vectors

    !Runge Kutta 4 order
    !For du/dt(t, u), double old, new
    subroutine dRK4_t(t, old, new, f, dt, dim)
        real*8, intent(in)::t
        integer, intent(in)::dim
        real*8, dimension(dim), intent(in)::old
        real*8, dimension(dim), intent(out)::new
        external::f
        real*8,intent(in)::dt
        real*8::dtd2
        real*8,dimension(dim)::k1, k2, k3, k4
        dtd2 = dt / 2d0
        call f(k1, t       , old            , dim)
        call f(k2, t + dtd2, old + k1 * dtd2, dim)
        call f(k3, t + dtd2, old + k2 * dtd2, dim)
        call f(k4, t + dt  , old + k3 * dt  , dim)
        new = old + dt / 6d0 * (k1 + 2d0 * k2 + 2d0 * k3 + k4)
    end subroutine dRK4_t
    !For du/dt(u), double old, new
    subroutine dRK4(old, new, f, dt, dim)
        integer,intent(in)::dim
        real*8,dimension(dim),intent(in)::old
        real*8,dimension(dim),intent(out)::new
        external::f
        real*8,intent(in)::dt
        real*8::dtd2
        real*8,dimension(dim)::k1,k2,k3,k4
        dtd2=dt/2d0
        call f(k1,old,dim)
        call f(k2,old+k1*dtd2,dim)
        call f(k3,old+k2*dtd2,dim)
        call f(k4,old+k3*dt,dim)
        new=old+dt/6d0*(k1+2d0*k2+2d0*k3+k4)
    end subroutine dRK4
    !For du/dt(u), complex old, new
    subroutine zRK4(old, new, f, dt, dim)
        integer,intent(in)::dim
        complex*16,dimension(dim),intent(in)::old
        complex*16,dimension(dim),intent(out)::new
        external::f
        real*8,intent(in)::dt
        real*8::dtd2
        complex*16,dimension(dim)::k1,k2,k3,k4
        dtd2=dt/2d0
        call f(k1,old,dim)
        call f(k2,old+k1*dtd2,dim)
        call f(k3,old+k2*dtd2,dim)
        call f(k4,old+k3*dt,dim)
        new=old+dt/6d0*(k1+2d0*k2+2d0*k3+k4)
    end subroutine zRK4

    !Perdiction-correction 2 order, predictor = Euler, corrector = backward Euler
    subroutine dPredictCorrect2(old, new, f, dt, dim, Warning, MaxIteration, Precision)
        !Required argument
            integer,intent(in)::dim
            real*8,dimension(dim),intent(in)::old
            real*8,dimension(dim),intent(inout)::new
            external::f
            real*8,intent(in)::dt
        !Optional argument
            logical,intent(in),optional::Warning
            integer,intent(in),optional::MaxIteration
            real*8 ,intent(in),optional::Precision
        logical::warn
        integer::maxit,i
        real*8::tol,dtd2,absdev,reldev
        real*8,dimension(dim)::k,olditer,kiter
        !Set parameter according to optional argument
            if(present(Warning)) then
                warn=Warning
            else
                warn=.true.
            end if
            if(present(MaxIteration)) then
                maxit=MaxIteration
            else
                maxit=20
            end if
            if(present(Precision)) then
                tol=precision
            else
                tol=1d-15
            end if
        dtd2=dt/2d0
        call f(k,old,dim)
        olditer=old+dt*k
        call f(kiter,olditer,dim)
        olditer=old+dtd2*(k+kiter)
        do i=1,maxit
            call f(kiter,olditer,dim)
            new=old+dtd2*(k+kiter)
            absdev=maxval(Abs(new-olditer))
            reldev=maxval(Abs((new-olditer)/new))
            if(absdev<tol.or.reldev<tol) exit
            olditer=new
        end do
        if(i>maxit.and.warn) then
            write(*,'(1x,A53)')'Failed perdiction-correction: max iteration exceeded!'
            write(*,*)'absdev =',absdev
            write(*,*)'reldev =',reldev
        end if
    end subroutine dPredictCorrect2
!----------------------- End ------------------------

!------------------- Integration --------------------
    !Integrate[f(x), {x, low, up}]
    !Optional argument:
    !    MinIteration: (default = 4) numerical integration may fraudulently converge with too few grids
    !                  at least 2^MinIteration grids will be used
    !                  MinIteration >= 4, because Romberg series require at least 16 grids
    !    MaxIteration (default = 25), AbsolutePrecision (default = 1d-15), RelativePrecision (default = 1d-15)
    real*8 function dRomberg(f, low, up, Warning, MinIteration, MaxIteration, AbsolutePrecision, RelativePrecision)
        !Required argument
            real*8,external::f
            real*8,intent(in)::low,up
        !Optional argument
            logical,intent(in),optional::Warning
            integer,intent(in),optional::MinIteration,MaxIteration
            real*8 ,intent(in),optional::AbsolutePrecision,RelativePrecision
        integer::minit,maxit,i,j,grids
        real*8::abstol,reltol,absdev,reldev,t,told,s,sold,c,cold,rold,dx,dxm,temp,start
        !Set parameter according to optional argument
            if(present(MinIteration)) then; minit=max(4,MinIteration)
                else; minit=4; end if
            if(present(MaxIteration)) then; maxit=MaxIteration
                else; maxit=25; end if
            if(present(AbsolutePrecision)) then; abstol=Abs(AbsolutePrecision)
                else; abstol=1d-15; end if
            if(present(RelativePrecision)) then; reltol=Abs(RelativePrecision)
                else; reltol=1d-15; end if
        dx=(up-low)/2d0!t_0
        t=(f(low)+f(up))*dx
        told=t!s_0
        t=told/2d0+dx*f(low+dx); s=(4d0*t-told)/3d0
        dx=dx/2d0
        told=t; sold=s!c_0
        t=told/2d0+dx*(f(low+dx)+f(up-dx)); s=(4d0*t-told)/3d0; c=(16d0*s-sold)/15d0
        dx=dx/2d0
        told=t; sold=s; cold=c!r_0
        temp=dx*3d0
        t=told/2d0+dx*(f(low+dx)+f(low+temp)+f(up-temp)+f(up-dx)); s=(4d0*t-told)/3d0; c=(16d0*s-sold)/15d0; dromberg=(64d0*c-cold)/63d0
        grids=8!Preiteration
        do i=4,minit
            temp=0d0; dxm=dx; dx=dx/2d0; grids=grids*2
            told=t; sold=s; cold=c; rold=dromberg
            start=low-dx
            do j=1,grids-1,2
                start=start+dxm; temp=temp+f(start)
            end do
            t=told/2d0+dx*temp; s=(4d0*t-told)/3d0; c=(16d0*s-sold)/15d0; dromberg=(64d0*c-cold)/63d0
        end do
        do i=minit+1,maxit!Main loop
            temp=0d0; dxm=dx; dx=dx/2d0; grids=grids*2
            told=t; sold=s; cold=c; rold=dromberg
            start=low-dx
            do j=1,grids-1,2
                start=start+dxm; temp=temp+f(start)
            end do
            t=told/2d0+dx*temp; s=(4d0*t-told)/3d0; c=(16d0*s-sold)/15d0; dromberg=(64d0*c-cold)/63d0
            absdev=abs(dromberg-rold); reldev=absdev/abs(dromberg)
            if(absdev<abstol.or.reldev<reltol) exit
        end do
        if(i>maxit.and.((.not.present(Warning)).or.Warning)) then
            write(*,'(1x,A46)')'Failed Romberg integration: max divide exceed!'
            write(*,*)'Absolute deviation =',absdev
            write(*,*)'Relative deviation =',reldev
        end if
    end function dRomberg
!----------------------- End ------------------------

end module Mathematics