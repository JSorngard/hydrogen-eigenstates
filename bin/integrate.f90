pure function integrate(x,y) result(r)
    !Calculates the integral of an array y with respect to x using the trapezoid
    !approximation. Note that the mesh spacing of x does not have to be uniform.
    real*8,intent(in)::x(:)           !Variable x
    real*8,intent(in)::y(size(x))     !Function y(x)
    real*8::r                          !integral int(y(x))dx

    !Integrate using trapezoidal rule
    associate(n=>size(x))
        r=sum((y(1+1:n-0)+y(1+0:n-1))*(x(1+1:n-0)-x(1+0:n-1)))/2
    end associate
end function