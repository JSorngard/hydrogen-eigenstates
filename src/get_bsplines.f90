!left: the index of the knotpoint directly to the left
!grid: the array of knotpoints
!npg: the number of elements in the knotpoint array
!b_matrix: the matrix to which all bsplines are written
!pos: the array containing all x-positions where the b-splines are to be evaluated
!npa: the number of elements in the pos array
subroutine get_bsplines(left,grid,npg,b_matrix,pos,npa)
implicit none

integer,intent(in)::npg,npa,left			!length grid, length abcs, matrix dimensions and left
real*8,intent(in),dimension(npg)::grid
real*8,intent(out),dimension(npa,npa-1)::b_matrix
real*8,intent(in),dimension(npa)::pos

real*8,dimension(npa-1)::bsp_tmp
integer::k


do k=1, npa
	call bsplvb(grid,npa-1,1,pos(k),left,bsp_tmp)
	b_matrix(k,:) = bsp_tmp
	
end do

end subroutine get_bsplines
