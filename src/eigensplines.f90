program eigensplines
    use iso_fortran_env
    implicit none

    !Physical constants
    real(real64), parameter :: rb=5.29177211E-11, hplanck=6.62607E-34, me=9.109383E-31
    real(real64), parameter :: pi=4*atan(1.d0), hbar=hplanck/(2*pi), epsilonnaught=8.854E-12, Q=1.6E-19
    integer,  parameter :: Z=1
    
    !Variables, can be changed to change the system
    integer,  parameter :: l=0, nn=1
    integer,  parameter :: k=7, pts=50, out_unit=20
    integer,  parameter :: neqs=pts-2*(k-1), nsplines=pts-k, n=3, nqpts=k+1, syssize=nsplines-2
    real(real64), parameter :: lstkntpt=50., lstkntptSI=10*rb
    logical,  parameter :: iprint=.false.

    !Intermediate variables used for computation
    integer :: i, resolution
    real(real64), dimension(nqpts) :: abscissas, weights
    real(real64), external :: bspline, dbspline, glquad
    real(real64), dimension(pts) :: kntpts
    real(real64), dimension(syssize, syssize) :: LHS, RHS
    real(real64), dimension(nsplines, syssize) :: nVR
    real(real64), dimension(syssize) :: ALPHAR, ALPHAI, BETA
    real(real64), dimension(syssize, syssize) :: VL, VR
    complex(real64), dimension(syssize) :: eigens
    integer, dimension(syssize) :: ipiv
    logical, external :: isinf
    integer, external :: splinestart, splineend
    real(real64) :: x

    !Generate the abscissas and weights for future use of Gaussian quadrature
    call gauleg(-1.d0, 1.d0, abscissas, weights, nqpts)

    !Generate a knot sequence with equal spacing from 0 to lstkntpt.
    do i = 1, k
        kntpts(i) = 0.d0
    enddo
    do i = k+1, pts-k+1
        kntpts(i) = kntpts(i-1) + lstkntpt/(neqs-1)
    enddo
    do i = pts-k+2, pts
        kntpts(i) = kntpts(i-1)
    enddo

    if(iprint) then
        write(*,*) "Generating equation..."
    end if

    call build_equation(l, kntpts, pts, k, LHS, RHS, abscissas, weights, nsplines, nqpts)

    if(iprint) then
        write(*,*) "Solving..."
    end if

    call solve_eigensystem(LHS, RHS, syssize, ALPHAR, ALPHAI, BETA, VL, VR, eigens)

    !Sort the eigenvalues by energy
    if(iprint) then
        write(*,*) "Sorting eigenvalues..."
    end if
    do i = 1, syssize
        ipiv(i) = i
    enddo
    call qsort(1, syssize, eigens, ipiv, syssize)
    !Reorder the array of eigen vectors in the same way
    VR = VR(:,ipiv)
    do i = 1, syssize
        nVR(:,i) = (/0.d0, VR(:,i), 0.d0/)
    enddo

    write(*,*) "The negative eigenvalues are:"
    do i = 1, syssize
        if (dble(eigens(i)) < 0) then
            write(*,*) dble(eigens(i))
        endif
    enddo
    
    if(iprint) then
        write(*,*) "Writing wavefunction to file..."
    end if
    resolution = 1000
    open(unit=out_unit, file="wavefunc.txt", action="write", status="replace")
    do i = 0, resolution
        x = lstkntpt*dble(i)/resolution
        
        write(out_unit,*) x,",",resum_splines(kntpts,pts,k,nVR(:,nn),size(nVR(:,nn)),x)**2.d0
    enddo
    close(out_unit)

contains

    subroutine solve_eigensystem(LHS,RHS,syssize,ALPHAR,ALPHAI,BETA,VL,VR,eigens)
        integer,         intent(in)                                 :: syssize
        real(real64),    intent(in),  dimension(syssize, syssize)   :: LHS, RHS
        real(real64),    intent(out), dimension(syssize)            :: ALPHAR, ALPHAI, BETA
        real(real64),    intent(out), dimension(syssize, syssize)   :: VL, VR
        complex(real64), intent(out), dimension(syssize)            :: eigens

        real(real64), dimension(3) :: array1
        integer :: INFO, LWORK
        real(real64), dimension(:), allocatable :: WORK

        !Zeroing memory for lapack
        ALPHAR = 0.d0
        ALPHAI = 0.d0
        BETA   = 0.d0
        VR     = 0.d0
        VL     = 0.d0

        !Determine optimal size of working memory
        call dggev('N','V',syssize,LHS,syssize,RHS,syssize,ALPHAR,ALPHAI,BETA,VL,syssize,VR,syssize,array1,-1,INFO)
        LWORK=array1(1)
        allocate(WORK(LWORK))

        !Zeroing memory for lapack
        ALPHAR = 0.d0
        ALPHAI = 0.d0
        BETA   = 0.d0
        VR     = 0.d0
        VL     = 0.d0
        WORK   = 0.d0

        !Compute eigenvalues
        call dggev('N','V',syssize,LHS,syssize,RHS,syssize,ALPHAR,ALPHAI,BETA,VL,syssize,VR,syssize,WORK,LWORK,INFO)
        
        if(INFO /= 0) then
            write(*,*) "WARNING: LAPACK response code: ",INFO," when solving the system. System size: ",syssize
            stop
        end if

        !Work out eigenvalues from lapack response
        do i = 1, syssize
            eigens(i) = complex(ALPHAR(i)/BETA(i),ALPHAI(i)/BETA(i))
        enddo
    end subroutine solve_eigensystem


    subroutine build_equation(l,kntpts,pts,k,LHS,RHS,abscissas,weights,nsplines,nqpts)
        integer,      intent(in)                                        :: l, k, pts
        integer,      intent(in)                                        :: nsplines, nqpts
        real(real64), intent(in),    dimension(pts)                     :: kntpts
        real(real64), intent(in),    dimension(nqpts)                   :: abscissas, weights
        real(real64), intent(inout), dimension(nsplines-2,nsplines-2)   :: LHS, RHS
        
        real(real64), dimension(nqpts) :: dx
        real(real64), dimension(nsplines,nsplines) :: preLHS,preRHS
        real(real64), dimension(k) :: splinestore, dsplinestore
        real(real64), dimension(k,k) :: intstorel, intstorer, matstore, dmatstore
        real(real64) :: xm, xr, a, b, x, centrifugal, potential
        integer :: left

        LHS    = 0.d0
        RHS    = 0.d0
        preLHS = 0.d0
        preRHS = 0.d0

        do left=k,nsplines !loop through all pairs of knot points for integration
            intstorel    = 0.d0
            intstorer    = 0.d0
            matstore     = 0.d0
            dmatstore    = 0.d0
            splinestore  = 0.d0
            dsplinestore = 0.d0

            a = kntpts(left)
            b = kntpts(left+1)
            
            xm = 0.5d0*(b+a)
            xr = 0.5d0*(b-a)
            dx = xr*abscissas

            do i=1,k+1 !integration loop
                x=dx(i)+xm
                call getsplines(kntpts,pts,k,x,splinestore)
                call getdsplines(kntpts,pts,k,x,dsplinestore,1)

                call outer_product(splinestore,splinestore,matstore,k,k)
                intstorer=intstorer+xr*weights(i)*matstore
               
                call outer_product(dsplinestore,dsplinestore,dmatstore,k,k)
                dmatstore=dmatstore*.5d0
                centrifugal=.5d0*dble(l)*(dble(l)+1.d0)/(x**2.d0)
                potential=-dble(Z)/x
                matstore=matstore*(centrifugal+potential)
                intstorel=intstorel+xr*weights(i)*(dmatstore+matstore)
            enddo
            preLHS(left-k+1:left,left-k+1:left)=preLHS(left-k+1:left,left-k+1:left)+intstorel
            preRHS(left-k+1:left,left-k+1:left)=preRHS(left-k+1:left,left-k+1:left)+intstorer
        enddo


        LHS=preLHS(2:nsplines-1,2:nsplines-1)
        RHS=preRHS(2:nsplines-1,2:nsplines-1)
    end subroutine build_equation


    subroutine outer_product(vec1,vec2,result,n,m)
        integer,      intent(in)                    :: n, m
        real(real64), intent(in),  dimension(m)     :: vec1
        real(real64), intent(in),  dimension(n)     :: vec2
        real(real64), intent(out), dimension(m,n)   :: result
        result = 0.d0
        call dger(m,n,1.d0,vec1,1,vec2,1,result,m)
    end subroutine outer_product


    function resum_splines(kntpts,pts,k,expcoeffs,coeffs,x)
        integer,intent(in) :: pts,coeffs,k
        real(real64),intent(in) :: x
        real(real64),dimension(pts),intent(in) :: kntpts
        real(real64),dimension(coeffs),intent(in) :: expcoeffs
        real(real64),dimension(k) :: splines
        real(real64) :: resum_splines
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

