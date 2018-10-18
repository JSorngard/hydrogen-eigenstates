!Given the leftmost point of a pair of adjacent knotpoints, returns the position
!of bspline i in the vector returned by bsplvb, 0 if it is not present, or -j if
!the j:th input is invalid
function getdbsplvbindex(N,i,k,left)
    implicit none
    integer,intent(in) :: N,i,k,left
    integer,dimension(k)::indices
    integer :: getdbsplvbindex,j
    integer,external::frstinstind
    
    !Initialize output variable
    getdbsplvbindex=0
    
    !Verification of input
    if (N<2*k+1) then
        getdbsplvbindex=-1
    endif
    if (i<1 .or.(i>N-k)) then
        getdbsplvbindex=-2
    endif
    if (k<1) then
        getdbsplvbindex=-3
    endif
    if ((left<k).or.(left>(N-k+1))) then
        getdbsplvbindex=-4
    endif

    !Computation of index
    if (getdbsplvbindex>=0) then
        do j=1,k
            indices(j)=left-k+j
        enddo
        getdbsplvbindex=frstinstind(indices,size(indices),i)
    endif

    if (getdbsplvbindex<0) then
        write(*,*) "getbsplvb error: the input is invalid at position", -getdbsplvbindex
    elseif (getdbsplvbindex==0) then
        write(*,*) "getbsplvb error: the requested bspline has no value between the given knotpoints"
    endif

end function getdbsplvbindex