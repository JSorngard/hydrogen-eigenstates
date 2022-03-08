subroutine getsplines(kntpts,n,k,x,out)
    use iso_fortran_env
    implicit none
    real(real64), intent(in)                    :: x
    integer,      intent(in)                    :: k , n
    real(real64), intent(in),    dimension(n)   :: kntpts
    real(real64), intent(inout), dimension(k)   :: out

    integer :: i, left
    logical, external :: isinf

    left = 0

    if (size(kntpts) < 2*k + 1) then
        left = -1
    elseif (n < 2*k + 1) then
        left = -2
    elseif(k < 1)then
        left = -3
    elseif(x < kntpts(k))then
        left = -4
    elseif(size(out) /= k)then
        left = -5
    endif

    if (left < 0) then
        write(*,*) "getsplines error: the input is invalid in argument",-left
    endif
    
    do i = k, n-k+2 !Determine which knotpoint is the one directly to the left of x
        if (kntpts(i) > x) then
            left = i - 1
            exit
        endif
    enddo

    if (x < kntpts(1) .or. x > kntpts(size(kntpts))) then
        do i = 1, k
            out(i) = 0.d0
        enddo
    endif
    if (x > kntpts(size(kntpts))) then
        write(*,*) "getsplines error: requested splines outside the domain (right)"
        write(*,*) "at x=",x
    elseif (x < kntpts(1)) then
        write(*,*) "getsplines error: requested splines outside the domain (left)"
        write(*,*) "at x=",x
    endif

    call bsplvb(kntpts, k, 1, x, left, out)
    
    do i = 1, size(out)
        if(isinf(out(i))) then
            write(*,*) i,"getsplines error: results from bsplvb are infinity"
            exit
        endif
    enddo
end subroutine getsplines