subroutine summarizearray(a,n,k)
    implicit none
    integer,intent(in) :: n,k
    integer :: i
    real*8,dimension(n),intent(in) :: a

    if (n<1) then
        write(*,*) "summarizearray error: invalid array size"
    endif
    if (k<0) then
        write(*,*) "summarizearray error: invalid trim"
    endif
    if (n>2*k+3) then
        do i=1,min(n,k)
            write(*,"(f11.3)") a(i)
        enddo

        write(*,*) "       ."
        write(*,*) "       ."
        write(*,*) "       ."

        do i=max(n-k+1,1),n
            write(*,"(f11.3)") a(i)
        enddo 
    else
        write(*,"(f11.3)") a
    endif
end subroutine summarizearray