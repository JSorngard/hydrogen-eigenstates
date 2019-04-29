function sglquad(abscissas,weights,n,a,b,func,kntpts,pts,rowindex,columnindex,k,l)
    implicit none

    !Variable declaration
    integer,intent(in)::n,rowindex,columnindex,k,pts,l
    real*8,intent(in),dimension(n)::abscissas,weights
    real*8,intent(in)::a,b
    real*8,intent(in),dimension(pts)::kntpts
    real*8::xm,xr,s,dx,sglquad
    integer::i

    !Interfaces
    interface
        real*8 function func(kntpts,pts,j,i,k,l,x)
            real*8,intent(in)::x
            integer,intent(in)::pts
            real*8,intent(in),dimension(pts)::kntpts
            integer,intent(in)::j,i,k,l
        end function func
    end interface

    !Initialize outout
    sglquad=0

    !Calculate domain scalings
    xm=0.5*(b+a)
    xr=0.5*(b-a)
    s=0
    dx=0
    !write(*,*) func(kntpts,pts,rowindex,columnindex,k,l,xr*abscissas(1))
    !Perform the quadrature sum
    do i=1,n/2
        dx=xr*abscissas(i)
        s=s+weights(i)*(func(kntpts,pts,rowindex,columnindex,k,l,xm+dx)+func(kntpts,pts,rowindex,columnindex,k,l,xm-dx))        
    enddo
    sglquad=s*xr
end function sglquad