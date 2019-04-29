!Calculates the first or n:th derivative of the BSpline 
recursive function dBSpline(kntpts,n,i,k,origk,x,o) result(res)
    implicit none

    integer,intent(in)::n
    real*8, dimension(n),intent(in):: kntpts
    integer,intent(in)::i,k,origk
    integer,intent(in),optional::o
    real*8,intent(in)::x
    real*8::a,b,splinea,splineb,res
    integer::order
    real*8::BSpline
    
    if (present(o)) then
        order=o
    else
        order=1
    endif
    
    !a=kntpts(i+k)-kntpts(i)
    a=kntpts(i+k-1)-kntpts(i)
    if (a/=0) then
        !a=k/a
        a=(k-1)/a
    else
        a=0
    endif

    !b=kntpts(i+k+1)-kntpts(i+1)
    b=kntpts(i+k)-kntpts(i+1)
    if (b/=0) then
        !b=-k/b
        b=-(k-1)/b
    else
        b=0
    endif
    
    splinea=0
    splineb=0
    if (order==1) then
        if (a/=0) then
            splinea=BSpline(kntpts,n,i,k-1,origk,x)
        endif
        if (b/=0) then
            splineb=BSpline(kntpts,n,i+1,k-1,origk,x)
        endif
        res=a*splinea+b*splineb    
    else
        if (a/=0) then
            splinea=dBSpline(kntpts,n,i,k-1,origk,x,order-1)
        endif
        if (b/=0) then
            splineb=dBSpline(kntpts,n,i+1,k-1,origk,x,order-1)
        endif
        res=a*splinea+b*splineb
    endif
end function dBSpline