function splinestart(i,k)
    implicit none
    integer,intent(in)::i,k
    integer::splinestart
    splinestart=max(k,i)
end function splinestart

function splineend(N,i,k)
    implicit none
    integer,intent(in) :: N,i,k
    integer :: splineend
    splineend=min(N-k+1,k+i)
end function splineend