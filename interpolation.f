      module interpolation
      private
      public fival,cfival,FFC,FFR4
      contains
!************************************************************************
!*     REAL 4-point lagrange interpolation routine.
!*     interpolates thr FUNCTION value fival at point r from an
!*     array of points stored in fdis(ndm). this array is assumed
!*     to be defined such that the first element fdis(1) CONTAINS
!*     the FUNCTION value at r=xv(1) and xv(2 .. ndm) are monotonically
!*     increasing.
!************************************************************************
      FUNCTION cfival(r,xv,fdis,ndm,alpha)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 cfival,fdis(ndm),y1,y2,y3,y4
      DIMENSION xv(ndm)
      IF(r.GT.xv(ndm)) go to 9
      DO 5 k=1,ndm-2
 5    IF(r.LT.xv(k)) go to 6
      k=ndm-2
 6    nst=MAX(k-1,1)
      x1=xv(nst)
      x2=xv(nst+1)
      x3=xv(nst+2)
      x4=xv(nst+3)
      y1=fdis(nst+0)
      y2=fdis(nst+1)
      y3=fdis(nst+2)
      y4=fdis(nst+3)
      pii1=(x1-x2)*(x1-x3)*(x1-x4)
      pii2=(x2-x1)*(x2-x3)*(x2-x4)
      pii3=(x3-x1)*(x3-x2)*(x3-x4)
      pii4=(x4-x1)*(x4-x2)*(x4-x3)
      xd1=r-x1
      xd2=r-x2
      xd3=r-x3
      xd4=r-x4
      pi1=xd2*xd3*xd4
      pi2=xd1*xd3*xd4
      pi3=xd1*xd2*xd4
      pi4=xd1*xd2*xd3
      cfival=y1*pi1/pii1+y2*pi2/pii2+y3*pi3/pii3+y4*pi4/pii4
      RETURN
 9    cfival=fdis(ndm) * EXP(alpha*(xv(ndm)-r))
      RETURN
      END


c-----------------------------------------------------------------------

************************************************************************
*     REAL 4-point lagrange interpolation routine.
*     interpolates thr FUNCTION value fival at point r from an
*     array of points stored in fdis(ndm). this array is assumed
*     to be defined such that the first element fdis(1) CONTAINS
*     the FUNCTION value at r=xv(1) and xv(2 .. ndm) are monotonically
*     increasing.
************************************************************************
      FUNCTION fival(r,xv,fdis,ndm,alpha)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 fdis(ndm),y1,y2,y3,y4
      DIMENSION xv(ndm)
      IF(r.GT.xv(ndm)) go to 9
      DO 5 k=1,ndm-2
 5    IF(r.LT.xv(k)) go to 6
      k=ndm-2
 6    nst=MAX(k-1,1)
      x1=xv(nst)
      x2=xv(nst+1)
      x3=xv(nst+2)
      x4=xv(nst+3)
      y1=fdis(nst+0)
      y2=fdis(nst+1)
      y3=fdis(nst+2)
      y4=fdis(nst+3)
      pii1=(x1-x2)*(x1-x3)*(x1-x4)
      pii2=(x2-x1)*(x2-x3)*(x2-x4)
      pii3=(x3-x1)*(x3-x2)*(x3-x4)
      pii4=(x4-x1)*(x4-x2)*(x4-x3)
      xd1=r-x1
      xd2=r-x2
      xd3=r-x3
      xd4=r-x4
      pi1=xd2*xd3*xd4
      pi2=xd1*xd3*xd4
      pi3=xd1*xd2*xd4
      pi4=xd1*xd2*xd3
      fival=y1*pi1/pii1+y2*pi2/pii2+y3*pi3/pii3+y4*pi4/pii4
      RETURN
 9    fival=fdis(ndm) * EXP(alpha*(xv(ndm)-r))
      RETURN
      END

c-----------------------------------------------------------------------

      FUNCTION FFC(PP,F,N)
      COMPLEX*16 FFC,F(N)
      REAL*8 PP
      PARAMETER(X=.16666666666667)
      I=PP
      IF(I.LE.0) GO TO 2
      IF(I.GE.N-2) GO TO 4
    1 P=PP-I
      P1=P-1.
      P2=P-2.
      Q=P+1.
      FFC=(-P2*F(I)+Q*F(I+3))*(P*P1*X)+(P1*F(I+1)-P*F(I+2))*(Q*P2*.5)
      RETURN
    2 IF(I.LT.0) GO TO 3
      I=1
      GO TO 1
    3 FFC=F(1)
      RETURN
    4 IF(I.GT.N-2) GO TO 5
      I=N-3
      GO TO 1
    5 FFC=F(N)
      RETURN
      END




      FUNCTION FFR4(Y,F,N)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 F(N),P,P1,P2,Q,X,FFR4
      REAL*8 Y
      PARAMETER(X=.16666666666667)
      P=Y
      I=P
      IF(I.LE.0) GO TO 2
      IF(I.GE.N-2) GO TO 4
    1 P=P-I
      P1=P-1.
      P2=P-2.
      Q=P+1.
      FFR4=(-P2*F(I)+Q*F(I+3))*P*P1*X+(P1*F(I+1)-P*F(I+2))*Q*P2*.5
      RETURN
    2 IF(I.LT.0) GO TO 3
      I=1
      GO TO 1
    3 FFR4=F(1)
      RETURN
    4 IF(I.GT.N-2) GO TO 5
      I=N-3
      GO TO 1
    5 FFR4=F(N)
      RETURN
      END


      end module
