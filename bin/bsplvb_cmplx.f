      SUBROUTINE bsplvb_cmplx(t,jhigh,index,x,left,biatx)
      IMPLICIT None
      integer jhigh,index,left
      complex*16 t(left+jhigh),biatx(jhigh)
      complex*16 x
*local
      integer i,j,jp1
      complex*16 saved,term
      complex*16 deltal(100),deltar(100)
      SAVE deltal,deltar
      DATA j/1/
      

      GO TO (10,20),index

 10   j = 1      
      biatx(1) = 1.d0            
      IF (j .GE. jhigh) GO TO 99

 20   CONTINUE
	 jp1 = j + 1
	 deltar(j) = t(left+j) - x
	 deltal(j) = x - t(left+1-j)         
	 saved = dcmplx(0.d0,0.d0)
	 DO i = 1,j
	     term = biatx(i)/(deltar(i) + deltal(jp1-i))
	     biatx(i) = saved + deltar(i)*term
	     saved = deltal(jp1-i)*term

	 END DO
	 biatx(jp1) = saved
	 j = jp1
	 IF (j .LT. jhigh) GO TO 20
 99   RETURN
      END
