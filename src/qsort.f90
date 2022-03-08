!(start,end,array1,array2,size)
!start: startindex of the vector to be sorted (most often 1)
!end: last index of the vector to be sorted (most often the size)
!array1: the vector to be sorted (min to max)
!array2: the vector that will contain the index switches so that other matrices
!that are disordered in the same way can be sorted instantly
!It must however contain the range start:end when the subroutine is called
!Usage of array2: sorted_matrix=matrix(:,array2)
!size: the size of the input vectors (array1 and array2)
!Author: Jimmy Vinbladh
recursive subroutine qsort(left,right,array,index_vec,len_arr)

implicit none

integer,intent(in)::len_arr
integer,intent(in)::left,right
integer,dimension(len_arr),intent(inout)::index_vec
complex*16,dimension(len_arr), intent(inout)::array

integer::bp,bl,br

if (right .gt. left)then
    call quick_sort(left,right,array,index_vec,len_arr,bp)
    br = bp-1
    bl = bp+1
    call qsort(left,br,array,index_vec,len_arr)
    call qsort(bl,right,array,index_vec,len_arr)
end if

end subroutine qsort

subroutine quick_sort(leftt,rightt,array,index_vec,len_arr,bp)
implicit none

integer,intent(in)::len_arr
integer,intent(in)::leftt,rightt
integer,dimension(len_arr),intent(inout)::index_vec
complex*16,dimension(len_arr), intent(inout)::array
integer,intent(out)::bp

integer::p,org_right,temp_ind, piv_ind,left,right
complex*16::temp_el,piv

left = leftt
right = rightt


p = floor((left+right)/2.d0)
piv         = array(p)
array(p)    = array(right)
piv_ind     = index_vec(p)
index_vec(p)= index_vec(right)
org_right   = right

do while (left .lt. right)
    right = right-1
	do while (real(array(left)) .le. real(piv) .and. left .lt. right)
        left = left+1
	end do
	do while (real(array(right)) .gt. real(piv) .and. right .gt. left)
        right = right-1
	end do
	if (real(array(left)) .gt. real(array(right)))then
        temp_el         = array(right)
        temp_ind        = index_vec(right)
        array(right)    = array(left)
        index_vec(right)= index_vec(left)
        array(left)     = temp_el
        index_vec(left) = temp_ind
	end if
end do


if (real(array(right)) .gt. real(piv)) then
    array(org_right)    = array(right)
    index_vec(org_right)= index_vec(right)
    array(right)        = piv
    index_vec(right)    = piv_ind
    bp                  = right
else
    array(org_right)    = array(right+1)
    index_vec(org_right)= index_vec(right+1)
    array(right+1)      = piv
    index_vec(right+1)  = piv_ind
    bp                  = right+1
end if

end subroutine quick_sort


