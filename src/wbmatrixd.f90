subroutine wbmatrixd(A,n,m)
    implicit none
    integer,intent(in)::n,m
    integer::i,j
    real*8,dimension(n,m),intent(in)::A
    write(*,*) "Structure:"
    do i=1,n
        do j=1,m
            if (A(i,j)/=0) then
                write(*,"(i1)",advance="no") 1
            else
                write(*,"(i1)",advance="no") 0
            endif
        enddo
        write(*,*)
    enddo
    write(*,*) "Value of bands"
    do i=1,n
        do j=1,m
            if (A(i,j)/=0) then
                write(*,"(d11.3)",advance="no") A(i,j)
            endif
        enddo
        write(*,*)
    enddo
end subroutine wbmatrixd