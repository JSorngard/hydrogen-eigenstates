program eigensplines
    implicit none

    !Variables
    integer,parameter::k=4,pts=30,neqs=pts-2*(k-1),nsplines=pts-k,n=3,l=1,nqpts=k+2,out_unit=20
    real*8,parameter::rb=5.29177211E-11,lstkntptSI=10*rb,lstkntpt=10.,hplanck=6.62607E-34,me=9.109383E-31
    real*8,parameter::pi=3.14159265359,hbar=hplanck/(2*pi),epsilonnaught=8.854E-12,Q=1.6E-19
    integer,parameter::Z=1,LWORK=936
    integer::posi,i,j,ii,jj,row,collumn,left,indiff,iii,jjj,INFO
    real*8,dimension(nqpts)::abscissas,weights
    real*8::r,temp,xm,xr,s,dx,a,b,x,b1,b2,db1,db2,temp1,temp2,temp3
    real*8,external::bspline,dbspline,glquad,func,legendrep,sphereV,sphereA,const,splinesquared,dsplinesquared
    real*8,dimension(pts)::kntpts
    real*8,dimension(nsplines,nsplines)::preLHS,preRHS
    real*8,dimension(nsplines-2,nsplines-2):: LHS,RHS,VL,VR
    real*8,dimension(nsplines-2)::ALPHAR,ALPHAI,BETA
    complex*16,dimension(nsplines-2)::eigens
    integer,dimension(nsplines-2)::ipiv
    real*8,dimension(LWORK)::WORK
    real*8,dimension(k,nqpts/2)::splinestore1,splinestore2,dsplinestore1,dsplinestore2
    real*8,dimension(k,nqpts,nsplines-k+1)::splinebox,dsplinebox
    real*8,dimension(k,nqpts)::bigsplinestore,bigdsplinestore
    real*8,dimension(nqpts/2)::poswghts,negwghts
    real*8,dimension(k)::sumstore,indices
    real*8,dimension(k,k)::intstore
    logical,external::isinf
    integer,external::splinestart,splineend
    real*8,dimension(3)::array1
    integer,dimension(3)::array2
    complex*8::enrg

    !Interfaces
    !Interface
    !    real*8 function spline(x,inkntpts,inn,ini,ink)
    !        real*8,intent(in)::x
    !        integer,optional::inn,ink,ini
    !        real*8,optional,dimension(:)::inkntpts
    !    end function spline
    !end interface

    !Initialization of variables
    call zeroes(preLHS,neqs,nsplines)
    
    !Generate the abscissas and weights for future use of Gaussian quadrature
    call gauleg(real(-1,kind=8),real(1,kind=8),abscissas,weights,nqpts)
    poswghts=weights((size(weights)/2+1):size(weights))
    negwghts=weights(1:size(weights)/2)

    !Generate a knot sequence with equal spacing from 0 to lstkntpt.
    do i=1,k
        kntpts(i)=real(0,kind=8)
    enddo
    do i=k+1,pts-k+1
        kntpts(i)=kntpts(i-1)+lstkntpt/(neqs-1)
    enddo
    do i=pts-k+2,pts
        kntpts(i)=kntpts(i-1)
    enddo

    !Loopa genom alla närliggande par av knutpunkter.                   -check
    !när (abs(i-j).le.(k-1)) så har integralen ett värde                -check
    !beräkna bsplvb och dbsplvb i alla abscissor mellan knutpunkterna.  -check
    !summera ihop värdena multiplicerat med vikterna.
    call zeroes(dsplinestore1,k,nqpts/2)
    call zeroes(dsplinestore2,k,nqpts/2)
    call zeroes(splinestore1,k,nqpts/2)
    call zeroes(splinestore2,k,nqpts/2)
    call zeroes(LHS,nsplines-2,nsplines-2)
    call zeroes(RHS,nsplines-2,nsplines-2)
    call zeroes(intstore,k,k)
    !write(*,"(4d11.3)") dsplinestore1(:,:)
    do left=k,nsplines !loop genom alla par av knutpunkter för integration
        a=kntpts(left) !denna loop tar fram värdet för bsplinesen i alla relevanta punkter
        b=kntpts(left+1)
        xm=0.5*(b+a)
        xr=0.5*(b-a)
        dx=0
        do j=(nqpts/2+1),nqpts
            dx=xr*abscissas(j)
            call getsplines(kntpts,pts,k,xm+dx,bigsplinestore(:,j))
            call getsplines(kntpts,pts,k,xm-dx,bigsplinestore(:,(nqpts+1-j)))
            call getdsplines(kntpts,pts,k,xm+dx,bigdsplinestore(:,j),1)
            call getdsplines(kntpts,pts,k,xm-dx,bigdsplinestore(:,(nqpts+1-j)),1)
            !splinebox(:,:,left-k+1)=bigsplinestore
            !dsplinebox(:,:,left-k+1)=bigdsplinestore
            do ii=1,k
                do jj=1,k
                    iii=left-k+ii !the i-index of the bspline at position ii of the vector returned by bsplvb
                    jjj=left-k+jj
                    s=0
                    if (abs(iii-jjj)<=k-1) then !if the integral would be nonzero
                        do i=1,nqpts !integration loop
                            dx=xr*abscissas(i)
                            x=dx+xm
                            b1=bigsplinestore(ii,i)*bigsplinestore(jj,i)
                            db1=bigdsplinestore(ii,i)*bigdsplinestore(jj,i)
                            temp1=real(0.5,kind=8)*db1
                            temp2=(real(l,kind=8)*(real(l,kind=8)+real(1,kind=8))/x**real(2,kind=8)-real(Z,kind=8)/x)*b1
                            s=s+xr*weights(i)*(temp1+temp2)
                            temp3=temp3+xr*weights(i)*b1*b2
                            !Vi måste få med rätt omskalat x!
                            !Fler knutpunkter?
                        enddo
                    endif
                    preLHS(iii,jjj)=preLHS(iii,jjj)+s
                    preLHS(jjj,iii)=preLHS(jjj,iii)+s
                    preRHS(iii,jjj)=preRHS(iii,jjj)+temp3
                    preRHS(jjj,iii)=preRHS(jjj,iii)+temp3
                enddo
            enddo
            !call getsplines(kntpts,pts,k,xm+dx,splinestore1(:,(j-nqpts/2)))
            !call getsplines(kntpts,pts,k,xm-dx,splinestore2(:,(j-nqpts/2)))
            !call getdsplines(kntpts,pts,k,xm+dx,dsplinestore1(:,(j-nqpts/2)),1)
            !call getdsplines(kntpts,pts,k,xm-dx,dsplinestore2(:,(j-nqpts/2)),1)
            !splinebox(:,nqpts+1-j,left-k+1)=splinestore2(:,j-nqpts/2)
            !splinebox(:,j,left-k+1)=splinestore1(:,j-nqpts/2)
            !dsplinebox(:,nqpts+1-j,left-k+1)=dsplinestore2(:,j-nqpts/2)
            !dsplinebox(:,j,left-k+1)=dsplinestore1(:,j-nqpts/2)
        enddo
    enddo

    LHS=preLHS(2:nsplines-1,2:nsplines-1)
    RHS=preRHS(2:nsplines-1,2:nsplines-1)
    !call wbmatrix(LHS,nsplines-2,nsplines-2)
    !call wbmatrix(RHS,nsplines-2,nsplines-2)
    call dggev("V","N",nsplines-2,LHS,nsplines-2,RHS,nsplines-2,ALPHAR,ALPHAI,BETA,VL,nsplines-2,VR,nsplines-2,WORK,LWORK,INFO)
    
    do i=1,size(ipiv)
        ipiv(i)=i
    enddo

    do i=1,size(ALPHAR)
        eigens(i)=complex(ALPHAR(i),ALPHAI(i)/BETA(i))
    enddo
    call qsort(1,size(eigens),eigens,ipiv,size(eigens))

    write(*,*) "The eigenvalues are:"
    do j=1,min(size(ALPHAR),5)
        write(*,"(f11.3)") real(eigens(j))
    enddo
    write(*,*) "       ."
    write(*,*) "       ."
    write(*,*) "       ."
    
    !write(*,*) 5*8*nsplines

    !open(unit=out_unit,file="dspl.txt",action="write",status="replace")
    !do i=0,100
    !    x=kntpts(k+1)*real(i,kind=8)/100
    !    call dbsplvb(kntpts,pts,k,1,x,k,dsplinestore1(:,3),1)
    !    write(out_unit,*) x,",",dsplinestore1(1,3)
    !enddo
    !close(out_unit)
    !open(unit=out_unit,file="spl.txt",action="write",status="replace")
    !do i=0,100
    !    x=kntpts(k+1)*real(i,kind=8)/100
    !    call bsplvb(kntpts,k,1,x,k,dsplinestore1(:,3))
    !    write(out_unit,*) x,",",dsplinestore1(1,3)
    !enddo
    !close(out_unit)
end program eigensplines

!function lhs(kntpts,pts,j,i,k,l,x)
!    implicit none
!    real*8,intent(in) :: x
!    integer,intent(in)::pts,j,i,k,l
!    real*8,intent(in),dimension(pts)::kntpts
!    real*8 :: lhs,a,b,c
!    real*8::dsplinesquared,splinesquared
!    real*8,external::dbspline,bspline
!
!    dsplinesquared=dbspline(kntpts,pts,j,k,k,x,1)*dbspline(kntpts,pts,i,k,k,x,1)
!    splinesquared=bspline(kntpts,pts,j,k,k,x)*bspline(kntpts,pts,i,k,k,x)
!
!    lhs=0.5*dsplinesquared+(l*(l+1)/x**2-1/x)*splinesquared    
!end function lhs

!function const(x)
!    implicit none
!    real*8,intent(in)::x
!    real*8::const
!    const=real(1,kind=8)
!end function const

!function func(x)
!    implicit none
!    real*8,intent(in)::x
!    real*8::func
!    func=x**2
!end function func

!function sphereV(r)
!    implicit none
!    real*8,intent(in)::r
!    real*8::sphereV
!    sphereV=4*3.14159265359*r**3/3
!end function sphereV

!function sphereA(r)
!    implicit none
!    real*8,intent(in)::r
!    real*8::sphereA
!    sphereA=4*3.14159265359*r**2
!end function sphereA

!function spline(x,inkntpts,inn,ini,ink)
!    implicit none
!    real*8,intent(in) :: x
!    integer,optional::inn,ink,ini
!    real*8,optional,dimension(:)::inkntpts
!    integer,save::nn,kk,ii
!    real*8,dimension(:),allocatable,save::kntptss
!    real*8 :: spline
!    real*8,external::bspline
!
!    if (present(inn)) then
!        nn=inn
!        if (allocated(kntptss)) then
!            deallocate(kntptss)
!        endif
!        allocate(kntptss(nn))
!    endif
!    if(present(inkntpts)) then
!        kntptss=inkntpts
!    endif
!    if (present(ink)) then
!        kk=ink
!    endif
!    if (present(ini)) then
!        ii=ini
!    endif
!
!    spline=bspline(kntptss,nn,ii,kk,kk,x)    
!end function spline

!function test(x)
!    implicit none
!    real*8,intent(in)::x
!    real*8::test
!    integer,parameter::nn=30,kk=4,neqsy=nn-2*(kk-1)
!    real*8,dimension(nn)::kntptss
!    real*8,external::bspline
!    real*8,parameter::rbz=5.29177211E-11,lstkntptz=10*rbz
!    integer::ii
!    do ii=kk+1,nn-kk+1
!        kntptss(ii)=kntptss(ii-1)+lstkntptz/(neqsy-1)
!    enddo
!    do ii=nn-kk+2,nn
!        kntptss(ii)=kntptss(ii-1)
!    enddo
!    test=bspline(kntptss,nn,1,kk,kk,x)
!    !write(*,*) x,test
!end function test