program eigensplines
    implicit none

    !Variables
    integer,parameter::k=7,pts=50,neqs=pts-2*(k-1),nsplines=pts-k,n=3,l=0,nqpts=k+1,out_unit=20
    real*8,parameter::rb=5.29177211E-11,lstkntptSI=10*rb,lstkntpt=50.,hplanck=6.62607E-34,me=9.109383E-31
    real*8,parameter::pi=3.14159265359,hbar=hplanck/(2*pi),epsilonnaught=8.854E-12,Q=1.6E-19
    integer,parameter::Z=1
    integer::posi,i,j,ii,jj,row,collumn,left,indiff,iii,jjj,INFO,LWORK,resolution
    real*8,dimension(nqpts)::abscissas,weights,dx
    real*8::r,temp,xm,xr,s,a,b,x,b1,b2,db1,db2,temp1,temp2,temp3,kinetic,centrifugal,potential
    real*8,external::bspline,dbspline,glquad,resum_splines
    real*8,dimension(pts)::kntpts
    real*8,dimension(nsplines,nsplines)::preLHS,preRHS
    real*8,dimension(nsplines-2,nsplines-2):: LHS,RHS,VL,VR
    real*8,dimension(nsplines,nsplines-2):: nVR,nVL
    real*8,dimension(nsplines-2)::ALPHAR,ALPHAI,BETA
    complex*16,dimension(nsplines-2)::eigens
    integer,dimension(nsplines-2)::ipiv
    real*8,dimension(:),allocatable::WORK
    real*8,dimension(k,nqpts/2)::splinestore1,splinestore2,dsplinestore1,dsplinestore2
    real*8,dimension(k,nqpts,nsplines-k+1)::splinebox,dsplinebox
    real*8,dimension(k,nqpts)::bigsplinestore,bigdsplinestore
    real*8,dimension(k)::sumstore,indices
    real*8,dimension(k,1)::splinestore,dsplinestore
    real*8,dimension(k,1)::teststore
    real*8,dimension(k,k)::intstorel,intstorer,matstore,dmatstore
    logical,external::isinf
    integer,external::splinestart,splineend
    real*8,dimension(3)::array1
    integer,dimension(3)::array2
    complex*8::enrg

    !Generate the abscissas and weights for future use of Gaussian quadrature
    call gauleg(-1.d0,1.d0,abscissas,weights,nqpts)

    !Generate a knot sequence with equal spacing from 0 to lstkntpt.
    do i=1,k
        kntpts(i)=0.d0
    enddo
    do i=k+1,pts-k+1
        kntpts(i)=kntpts(i-1)+lstkntpt/(neqs-1)
    enddo
    do i=pts-k+2,pts
        kntpts(i)=kntpts(i-1)
    enddo
    !write(*,"(f11.3)") kntpts

    !Zeroing memory
    LHS=0.d0
    RHS=0.d0
    preLHS=0.d0
    preRHS=0.d0
    
    do left=k,nsplines !loop genom alla par av knutpunkter för integration
        intstorel=0.d0
        intstorer=0.d0
        matstore=0.d0
        dmatstore=0.d0
        splinestore=0.d0
        dsplinestore=0.d0
        a=kntpts(left)
        b=kntpts(left+1)
        xm=0.5d0*(b+a)
        xr=0.5d0*(b-a)
        dx=xr*abscissas
        do i=1,nqpts !integration loop
            x=dx(i)+xm
            call getsplines(kntpts,pts,k,x,splinestore(:,1))
            call getdsplines(kntpts,pts,k,x,dsplinestore(:,1),1)

            matstore=matmul(splinestore,transpose(splinestore))
            intstorer=intstorer+xr*weights(i)*matstore
           
            dmatstore=matmul(dsplinestore,transpose(dsplinestore))
            dmatstore=dmatstore*.5d0
            centrifugal=.5d0*real(l,kind=8)*(real(l,kind=8)+1.d0)/(x**2.d0)
            potential=-real(Z,kind=8)/x
            matstore=matstore*(centrifugal+potential)
            intstorel=intstorel+xr*weights(i)*(dmatstore+matstore)
            !matstore=matstore*potential
            !intstorel=intstorel+xr*weights(i)*(matstore)
        enddo
        preLHS(left-k+1:left,left-k+1:left)=preLHS(left-k+1:left,left-k+1:left)+intstorel
        preRHS(left-k+1:left,left-k+1:left)=preRHS(left-k+1:left,left-k+1:left)+intstorer
    enddo


    LHS=preLHS(2:nsplines-1,2:nsplines-1)
    RHS=preRHS(2:nsplines-1,2:nsplines-1)
    
    !write(*,*) RHS(1,:)
    !call wbmatrix(LHS,nsplines-2,nsplines-2)
    !call wbmatrix(RHS,nsplines-2,nsplines-2)
    
    !Zeroing memory for lapack
    ALPHAR=0.d0
    ALPHAI=0.d0
    BETA=0.d0
    VR=0.d0
    VL=0.d0

    !Determine optimal size of working emmory
    call dggev('N','V',nsplines-2,LHS,nsplines-2,RHS,nsplines-2,ALPHAR,ALPHAI,BETA,VL,nsplines-2,VR,nsplines-2,array1,-1,INFO)
    LWORK=array1(1)
    allocate(WORK(LWORK))

    !Zeroing memory for lapack
    ALPHAR=0.d0
    ALPHAI=0.d0
    BETA=0.d0
    VR=0.d0
    VL=0.d0
    WORK=0.d0

    !Compute eigenvalues
    call dggev('N','V',nsplines-2,LHS,nsplines-2,RHS,nsplines-2,ALPHAR,ALPHAI,BETA,VL,nsplines-2,VR,nsplines-2,WORK,LWORK,INFO)
    !write(*,*) "Lapack response code: ",INFO
    do i=1,size(ALPHAR)
        eigens(i)=complex(ALPHAR(i)/BETA(i),ALPHAI(i)/BETA(i))
    enddo

    do i=1,size(ipiv)
        ipiv(i)=i
    enddo
    call qsort(1,size(eigens),eigens,ipiv,size(eigens))
    
    VR=VR(:,ipiv)
    do i=1,nsplines-2
        nVR(:,i)=(/0.d0,VR(:,i),0.d0/)
    enddo

    !write(*,*) "The sorted eigenvalues are:"
    !call summarizearray(real(eigens),size(real(eigens)),5)
    write(*,*) "The negative eigenvalues are:"
    do i=1,size(eigens)
        if (real(eigens(i),kind=8)<0) then
            write(*,*) real(eigens(i),kind=8)
        endif
    enddo
    
    resolution=1000
    open(unit=out_unit,file="wavefunc.txt",action="write",status="replace")
    do i=0,resolution
        x=lstkntpt*real(i,kind=8)/resolution
        
        write(out_unit,*) x,",",resum_splines(kntpts,pts,k,nVR(:,1),size(nVR(:,1)),x)**2.d0
    enddo
    close(out_unit)
    !open(unit=out_unit,file="spl.txt",action="write",status="replace")
    !do i=0,100
    !    x=kntpts(k+1)*real(i,kind=8)/100
    !    call bsplvb(kntpts,k,1,x,k,dsplinestore1(:,3))
    !    write(out_unit,*) x,",",dsplinestore1(1,3)
    !enddo
    !close(out_unit)
end program eigensplines

function resum_splines(kntpts,pts,k,expcoeffs,coeffs,x)
    implicit none
    integer,intent(in) :: pts,coeffs,k
    real*8,intent(in) :: x
    real*8,dimension(pts),intent(in) :: kntpts
    real*8,dimension(coeffs),intent(in) :: expcoeffs
    real*8,dimension(k) :: splines
    real*8 :: resum_splines
    integer :: i,left

    !Determine which knotpoint is the one directly to the left of x
    do i=1,pts
        if (kntpts(i)>x) then
            left=i-1
            exit
        endif
    enddo

    call getsplines(kntpts,pts,k,x,splines)

    resum_splines=sum(expcoeffs(left-k+1:left)*splines)
end function resum_splines