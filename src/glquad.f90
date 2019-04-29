function glquad(abscissas,weights,n,a,b,func)
    implicit none

    !Variable declaration
    integer,intent(in)::n
    real*8,intent(in),dimension(n)::abscissas,weights
    real*8,intent(in)::a,b
    real*8::xm,xr,s,dx,glquad
    integer::i

    !Interfaces
    interface
        real*8 function func(x)
            real*8,intent(in)::x
        end function func
    end interface

    !Initialize outout
    glquad=0

    !Calculate domain scalings
    xm=0.5*(b+a)
    xr=0.5*(b-a)
    s=0
    dx=0

    !Perform the quadrature sum
    !do i=1,n/2
    do i=n/2+1,n
        dx=xr*abscissas(i)
        s=s+weights(i)*(func(xm+dx)+func(xm-dx))
    enddo
    glquad=s*xr
    if (isnan(glquad)) then
        write(*,*) xm,xr, "WARNING: MEMORY CORRUPTION SUPPRESSION NOT WORKING"
    endif
end function glquad