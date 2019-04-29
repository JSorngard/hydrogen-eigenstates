!Returns the index of the bspline at position i of the vector returned by bsplvb
function getbsplineindex(left,N,i,k)
    implicit none
    integer,intent(in) :: left,N,i,k
    integer :: getbsplineindex
    
    getbsplineindex=0
    
    !Input verification
    if (left<k.or.left>N-k+1) then
        getbsplineindex=-1
    endif
    if (N<2*k+1) then
        getbsplineindex=-2
    endif
    if (i<left-k+1) then
        getbsplineindex=-3
    endif
    if (k<1) then
        getbsplineindex=-4
    endif

    getbsplineindex=left-k-i

    if (getbsplineindex<0) then
        write(*,*) "getbsplineindex error: the input is invalid at position", -getbsplineindex
    endif
    
end function getbsplineindex