function isinf(x)
    implicit none
    real*8,intent(in) :: x
    logical:: isinf
    real*8::dbl_prec_var
    real*8,parameter::infinity=huge(dbl_prec_var)

    isinf=((real(x,kind=8)>infinity).or.(real(x,kind=8))<-infinity)
end function isinf