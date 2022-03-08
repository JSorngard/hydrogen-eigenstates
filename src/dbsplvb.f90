recursive subroutine dbsplvb(kntpts,pts,order,index,x,left,out,dorder)
    !knutpunkter,antal knytpunkter, ordning, 1, x, index hos knutpunkten till vänster om x, rum för output
    use iso_fortran_env
    implicit none
    real(real64), intent(in)                        :: x
    integer,      intent(in)                        :: order, index, left, dorder, pts
    real(real64), intent(inout), dimension(order)   :: out
    
    real(real64), dimension(pts) :: kntpts
    real(real64), dimension(order) :: bsplines, losplines, hisplines, a, b, splinea, splineb
    integer :: i, mini, maxi

    if (left < order) then
        write(*,*) "dsplvb error: left is smaller than the index of the first physical point"
        write(*,*) "left is:",left,"the first physical point is at:",order
        write(*,*) "x=",x
    endif

    !Get all non-zero splines
    call bsplvb(kntpts, order-1, index, x, left, bsplines)

    bsplines(order) = 0.d0
    
    losplines = (/0.d0, bsplines(1:(order-1))/)

    !Add the fact that the b-spline of index left+order+1 is 0
    hisplines = (/bsplines(2:(order)), 0.d0/)

    !Calculate the mixing coefficients of the b-splines
    maxi = min(left, pts-order)
    mini = maxi - order + 1
    a = kntpts((mini+order-1):(maxi+order-1)) - kntpts((mini):(maxi))
    b = kntpts((mini+order):(maxi+order)) - kntpts((mini+1):(maxi+1))
    !Replace all elements of a and b with 1/a and 1/b respectively. 
    !Unless they contain a zero element, in which case that element is not edited
    if (any(a == 0.d0)) then
        do i = 1, size(a)
            if (a(i) /= 0.d0) then
                a(i) = 1.d0/a(i)
            endif
        enddo
    else
        a = 1.d0/a
    endif

    if (any(b == 0.d0)) then
        do i = 1, size(b)
            if (b(i) /= 0.d0) then
                b(i) = 1.d0/b(i)
            endif
        enddo
    else
        b = 1.d0/b
    endif
    if (dorder == 1) then
        out = (order - 1)*(a*losplines - b*bsplines)
    else
        call dbsplvb(kntpts, pts, order-1, index, x, left, losplines, dorder-1)
        losplines=(/0.d0, losplines(1:(order-2))/)
        call dbsplvb(kntpts, pts, order-1, index, x, left, bsplines, dorder-1)
        bsplines(order) = 0.d0
        bsplines(order-1) = 0.d0
        out = (order - 1)*(a*losplines - b*bsplines)
    endif
end subroutine dbsplvb