recursive function legendrep(n,x) result(P)
    implicit none
    real*8,intent(in)::x
    real*8::P
    integer,intent(in)::n
    if (n==0) then
        P=1
    elseif (n==1) then
        P=x
    else
        P=((2*n-3)*x*legendrep(n-1,x)-(n-1)*legendrep(n-2,x))/n
    endif
end function legendrep