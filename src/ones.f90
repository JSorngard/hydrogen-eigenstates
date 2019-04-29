subroutine ones(matrix,i,j)
    implicit none
    integer,intent(in)::i,j
    real*8,intent(inout),dimension(i,j)::matrix

    matrix(:,:)=real(1,kind=8)
end subroutine ones