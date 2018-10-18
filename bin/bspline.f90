!Returns the value of B-spline B(i,k;x) for knot points kntpts.
recursive function BSpline(kntpts,n,i,k,origk,x) result(result)
    implicit none

    !Variable declaration
    integer,intent(in)::n
    real*8, dimension(n), intent(in):: kntpts
    integer,intent(in):: i,k,origk
    real*8,intent(in):: x
    real*8:: a,b,splinea,splineb
    real*8:: result
    
    !Executable statements
    if (i>n-k .or. i<1) then
        write(*,*) "Index must be between 1 and ",n-k,", but is ",i
    endif

    if (k==1) then
        if ((kntpts(i) <= x .and. x < kntpts(i+1)).and.(kntpts(i)<kntpts(i+1))) then
            result=1
        else
            result=0
        endif
        return
    endif
    !If we are at the last point on the interval for the last B-spline we return its limit instead of its value.
    if (i==N-origk) then
        if (x==kntpts(N)) then
            result=1
            return
        endif
    endif

    !If k!=1 we calculate mixing coefficients (in cases with division by zero, the result is set to zero)
    a=kntpts(i+k-1)-kntpts(i)
    if (a/=0) then
        a=(x-kntpts(i))/a
    else
        a=0
    endif

    b=kntpts(i+k)-kntpts(i+1)
    if (b/=0) then
        b=(kntpts(i+k)-x)/b
    else
        b=0
    endif

    splinea=0
    splineb=0
    if (a/=0) then
        splinea=Bspline(kntpts,n,i,k-1,origk,x)
    endif
    if (b/=0) then
        splineb=BSpline(kntpts,n,i+1,k-1,origk,x)
    endif
    
    result=a*splinea+b*splineb
end function BSpline