!Returns inte index of the first instance of the specified element in vector v
function frstinstind(v,n,tar)
    implicit none
    integer,intent(in) :: n,tar
    integer,dimension(n),intent(in) :: v
    integer :: frstinstind,i
    do i=1,n
        if (v(i)==tar) then
            frstinstind=i
        else
            frstinstind=0
        endif
    enddo
end function frstinstind