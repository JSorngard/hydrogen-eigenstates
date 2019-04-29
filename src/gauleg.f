      SUBROUTINE GAULEG(X1,X2,XA,WA,NP)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DATA ZERO,ONE,TWO,THREE,FOUR,FIVE,SIX,EPS/
     : 0.0D0,1.0D0,2.0D0,3.0D0,4.0D0,5.0D0,6.0D0,3.D-15/
      DIMENSION XA(NP),WA(NP)
      M=(NP+1)/2
      XM=0.5D0*(X2+X1)
      XL=0.5D0*(X2-X1)


      DO 12 I=1,M
        Z=DCOS(3.14159265358979D0*(I-0.25D0)/(NP+0.5D0))

 1      CONTINUE
        P1=ONE
        P2=ZERO
        DO 11 J=1,NP
          P3=P2
          P2=P1
          P1=(DBLE(2*J-1)*Z*P2-DBLE(J-1)*P3)/DBLE(J)
 11     CONTINUE
        PP=DBLE(NP)*(Z*P1-P2)/(Z*Z-ONE)
        Z1=Z
        Z=Z1-P1/PP
        IF(ABS(Z-Z1).GT.EPS) GOTO 1
        XA(I)=XM-XL*Z
        XA(NP+1-I)=XM+XL*Z
        WA(I)=TWO*XL/((ONE-Z*Z)*PP*PP)
        WA(NP+1-I)=WA(I)

 12   CONTINUE
      RETURN
      END
