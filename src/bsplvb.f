*     Author: C. De. Boor
      SUBROUTINE bsplvb(t,jhigh,index,x,left,biatx)
*      INCLUDE 'standa.typ'
*      INCLUDE 'spldim.def'
*      include 'ndim_inc'
      PARAMETER (JMAX=100)
      integer index,jhigh,left,i,j,jp1
      real*8 t,x,biatx,deltal,deltar,saved,term
      DIMENSION biatx(jhigh),t(left+jhigh),deltal(jmax),deltar(jmax)
      SAVE deltal,deltar
      DATA j/1/
*      write(6,*) ' jmax=',jmax
      GO TO (10,20),index
 10   j = 1
      biatx(1) = 1.d0
      IF (j .GE. jhigh) GO TO 99

 20   CONTINUE
	 jp1 = j + 1
	 deltar(j) = t(left+j) - x
	 deltal(j) = x - t(left+1-j)
*         write(6,'(1pd12.4,2(i5,1pd14.6))')
*     :   x,left+j,t(left+j),left+1-j,t(left+1-j)
*         write(6,'(i3,1p3d12.4)') j,deltal(j),deltar(j),
*     :   abs(deltal(j)-deltar(j))
	 saved = 0.d0
	 DO i = 1,j
*         write(6,'(2i3,1p3d12.4)') i,j,deltal(jp1-1),deltar(i),
*     :   abs(deltal(jp1-1)-deltar(i))

	     term = biatx(i)/(deltar(i) + deltal(jp1-i))
	     biatx(i) = saved + deltar(i)*term
	     saved = deltal(jp1-i)*term
	 END DO
	 biatx(jp1) = saved
	 j = jp1
	 IF (j .LT. jhigh) GO TO 20
 99   RETURN
      END
