!left: the index of the knotpoint directly to the left
!grid: the array of knotpoints
!npg: the number of elements in the knotpoint array
!db_matrix: the matrix to which the derivatives of all bsplines are written
!pos: the array containing all x-positions where the b-splines are to be evaluated
!npa: the number of elements in the pos array
subroutine get_dbsplines(left,grid,npg,db_matrix,pos,npa)
implicit none

integer,intent(in)::npg,npa,left			!length grid, length abcs, matrix dimensions and left
real*8,intent(in),dimension(npa)::pos
real*8,intent(in),dimension(npg)::grid
real*8,intent(out),dimension(npa,npa-1)::db_matrix

real*8,dimension(npa-2)::bsp_tmp
real*8,dimension(npa,npa-2)::btmp_matrix
real*8,dimension(npa)::zero_vec,t1,t2
integer::k

do k=1, npa
	call bsplvb(grid,npa-2,1,pos(k),left,bsp_tmp)
	btmp_matrix(k,:) = bsp_tmp
end do

zero_vec = 0.d0

do k=1, npa-1
	if( k == 1)then
		t1 = zero_vec
		t2 = btmp_matrix(:,k)/(grid(left+k)-grid(left-npa+k+2))
	else if( k == npa-1)then
		t1 = btmp_matrix(:,k-1)/(grid(left+k-1)-grid(left-npa+k+1))
		t2 = zero_vec
	else
		t1 = btmp_matrix(:,k-1)/(grid(left+k-1)-grid(left-npa+k+1))
		t2 = btmp_matrix(:,k)/(grid(left+k)-grid(left-npa+k+2))
	end if
	
	db_matrix(:,k) = (npa-2)*(t1-t2) 
end do


end subroutine get_dbsplines
