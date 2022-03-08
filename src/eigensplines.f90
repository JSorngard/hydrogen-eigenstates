program eigensplines
    use iso_fortran_env
    implicit none

    !Variables
    integer,  parameter :: l=0,nn=1
    integer,  parameter :: k=7,pts=50,neqs=pts-2*(k-1),nsplines=pts-k,n=3,nqpts=k+1,out_unit=20
    real(real64), parameter :: rb=5.29177211E-11,lstkntptSI=10*rb,lstkntpt=50.,hplanck=6.62607E-34,me=9.109383E-31
    real(real64), parameter :: pi=3.14159265359,hbar=hplanck/(2*pi),epsilonnaught=8.854E-12,Q=1.6E-19
    integer,  parameter :: Z=1
    logical,  parameter :: iprint=.false.

    integer :: i,left,resolution
    real(real64), dimension(nqpts) :: abscissas,weights
    real(real64), external :: bspline,dbspline,glquad
    real(real64), dimension(pts) :: kntpts
    real(real64), dimension(nsplines-2,nsplines-2) :: LHS, RHS
    real(real64), dimension(nsplines,nsplines-2) :: nVR
    real(real64), dimension(nsplines-2) :: ALPHAR, ALPHAI, BETA
    real(real64), dimension(nsplines-2, nsplines-2) :: VL, VR
    complex*16, dimension(nsplines-2) :: eigens
    integer, dimension(nsplines-2) :: ipiv
    logical, external :: isinf
    integer, external :: splinestart,splineend
    real(real64) :: x

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

    if(iprint) then
        write(*,*) "Generating equation..."
    end if

    call build_equation(l,kntpts,pts,k,LHS,RHS,abscissas,weights,nsplines,nqpts)

    if(iprint) then
        write(*,*) "Solving..."
    end if

    call solve_eigensystem(LHS,RHS,nsplines,ALPHAR,ALPHAI,BETA,VL,VR,eigens)

    if(iprint) then
        write(*,*) "Sorting eigenvalues..."
    end if
    do i=1,size(ipiv)
        ipiv(i)=i
    enddo
    call qsort(1,size(eigens),eigens,ipiv,size(eigens))
    
    VR=VR(:,ipiv)
    do i=1,nsplines-2
        nVR(:,i)=(/0.d0,VR(:,i),0.d0/)
    enddo

    write(*,*) "The negative eigenvalues are:"
    do i=1,size(eigens)
        if (real(eigens(i),kind=8)<0) then
            write(*,*) real(eigens(i),kind=8)
        endif
    enddo
    
    if(iprint) then
        write(*,*) "Writing wavefunction to file..."
    end if
    resolution=1000
    open(unit=out_unit,file="wavefunc.txt",action="write",status="replace")
    do i=0,resolution
        x=lstkntpt*real(i,kind=8)/resolution
        
        write(out_unit,*) x,",",resum_splines(kntpts,pts,k,nVR(:,nn),size(nVR(:,nn)),x)**2.d0
    enddo
    close(out_unit)

contains

    subroutine solve_eigensystem(LHS,RHS,nsplines,ALPHAR,ALPHAI,BETA,VL,VR,eigens)
        integer,    intent(in)                                              :: nsplines
        real(real64),   intent(in), dimension(nsplines-2, nsplines-2)       :: LHS, RHS
        real(real64),   intent(out),   dimension(nsplines-2)                :: ALPHAR, ALPHAI, BETA
        real(real64),   intent(out),   dimension(nsplines-2, nsplines-2)    :: VL, VR
        complex*16, intent(out),   dimension(nsplines-2)                    :: eigens

        real(real64), dimension(3) :: array1
        integer :: INFO, LWORK
        real(real64), dimension(:), allocatable :: WORK

        !Zeroing memory for lapack
        ALPHAR=0.d0
        ALPHAI=0.d0
        BETA=0.d0
        VR=0.d0
        VL=0.d0

        !Determine optimal size of working memory
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
        
        if(INFO/=0) then
            write(*,*) "WARNING: LAPACK response code: ",INFO," when solving eigensystem"
        end if

        !Work out eigenvalues from lapack response
        do i=1,size(ALPHAR)
            eigens(i)=complex(ALPHAR(i)/BETA(i),ALPHAI(i)/BETA(i))
        enddo
    end subroutine solve_eigensystem

    subroutine build_equation(l,kntpts,pts,k,LHS,RHS,abscissas,weights,nsplines,nqpts)
        integer,  intent(in)                                            :: l, k, pts
        integer,  intent(in)                                            :: nsplines, nqpts
        real(real64), intent(in),    dimension(pts)                     :: kntpts
        real(real64), intent(in),    dimension(nqpts)                   :: abscissas, weights
        real(real64), intent(inout), dimension(nsplines-2,nsplines-2)   :: LHS, RHS
        
        real(real64), dimension(nqpts) :: dx
        real(real64), dimension(nsplines,nsplines) :: preLHS,preRHS
        real(real64), dimension(k,1) :: splinestore, dsplinestore
        real(real64), dimension(k,k) :: intstorel, intstorer, matstore, dmatstore
        real(real64) :: xm, xr, a, b, x, centrifugal, potential
        integer :: left


        !Zeroing memory
        LHS=0.d0
        RHS=0.d0
        preLHS=0.d0
        preRHS=0.d0

        do left=k,nsplines !loop through all pairs of knot points for integration
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
            do i=1,k+1 !integration loop
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
            enddo
            preLHS(left-k+1:left,left-k+1:left)=preLHS(left-k+1:left,left-k+1:left)+intstorel
            preRHS(left-k+1:left,left-k+1:left)=preRHS(left-k+1:left,left-k+1:left)+intstorer
        enddo


        LHS=preLHS(2:nsplines-1,2:nsplines-1)
        RHS=preRHS(2:nsplines-1,2:nsplines-1)
    end subroutine build_equation

    function resum_splines(kntpts,pts,k,expcoeffs,coeffs,x)
        integer,  parameter :: DP=kind(1.d0)
        integer,intent(in) :: pts,coeffs,k
        real(DP),intent(in) :: x
        real(DP),dimension(pts),intent(in) :: kntpts
        real(DP),dimension(coeffs),intent(in) :: expcoeffs
        real(DP),dimension(k) :: splines
        real(DP) :: resum_splines
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

end program eigensplines

