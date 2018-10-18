subroutine getdsplines(kntpts,n,k,x,out,o)
    implicit none
    real*8,intent(in) :: x
    integer,intent(in)::k,n
    integer,intent(in),optional::o
    real*8,dimension(n),intent(in)::kntpts
    real*8,dimension(k),intent(inout)::out
    integer::i,left,order
    logical,external::isinf


    left=0

    if (present(o)) then
        order=o
    else
        order=1
    endif

    if (size(kntpts)<2*k+1) then
        left=-1
    elseif (n<2*k+1) then
        left=-2
    elseif(k<1)then
        left=-3
    elseif(x<kntpts(k))then
        left=-4
    elseif(size(out)/=k)then
        left=-5
    elseif(order<1)then
        left=-6
    endif

    if (left<0) then
        write(*,*) "getdsplines error: the input is invalid in argument",-left
    endif

    do i=k,n-k+2 !Determine which knotpoint is the one directly to the left of x
        if (kntpts(i)>x) then
            left=i-1
            exit
        endif
    enddo

    if (i>n-k+1) then
        write(*,*) "getdsplines error: the knotpoint sequence has some problem (no point in the sequence is to the right of x)"
    endif

    call dbsplvb(kntpts,n,k,1,x,left,out,order)
    
    if (x<kntpts(k).or.x>kntpts(size(kntpts))) then
        do i=1,order
            out(i)=0
        enddo
    endif
    if (x<kntpts(k)) then
        write(*,*) "getdsplines error: requested splines outside the domain (left)"
        write(*,*) "at x=",x
    elseif (x>kntpts(size(kntpts))) then
        write(*,*) "getdsplines error: requested splines outside the domain (right)"
        write(*,*) "at x=",x
    endif


     do i=1,size(out)
        if(isinf(out(i))) then
            write(*,*) i,"getdsplines error: results from dbsplvb are infinity"
            exit
        endif
    enddo


end subroutine getdsplines