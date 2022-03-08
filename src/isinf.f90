function isinf(x)
    use iso_fortran_env
    implicit none
    
    real(real64), intent(in) :: x
    
    logical :: isinf
    real(real64) :: dbl_prec_var
    real(real64), parameter :: infinity=huge(dbl_prec_var)

    isinf=((dble(x)>infinity).or.(dble(x))<-infinity)
end function isinf