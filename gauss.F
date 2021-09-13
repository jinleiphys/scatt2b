C>     this version of gauleg is embedded in a module
C>     and includes the standard transformation trns
C>     tol is set to 1E-12, which should work with
C>     REAL = double precision


      MODULE gauss
       USE precision
       PRIVATE

       PUBLIC gauleg,trns,print_gauss,trnsmin,gauleg_4_parts,gauleg_6_parts

       CONTAINS

       ! divided the whole interval into 4 parts
       subroutine gauleg_4_parts(N4,x1,x2,X,W)
       implicit none
       INTEGER N, N4
       REAL(dpreal) x1,x2,X(N4),W(N4)
       real(dpreal) x_step

       N=N4/4
       x_step=(x2-x1)/4.0_dpreal

       call gauleg(N,x1,x1+x_step,X(1:N),W(1:N))
       call gauleg(N,x1+x_step,x1+2.*x_step,X(N+1:2*N),W(N+1:2*N))
       call gauleg(N,x1+2.*x_step,x1+3.*x_step,X(2*N+1:3*N),W(2*N+1:3*N))
       call gauleg(N,x1+3.*x_step,x2,X(3*N+1:N4),W(3*N+1:N4))

       end subroutine


       ! divided the whole interval into 6 parts
       subroutine gauleg_6_parts(N6,x1,x2,X,W)
       implicit none
       INTEGER N, N6
       REAL(dpreal) x1,x2,X(N6),W(N6)
       real(dpreal) x_step

       N=N6/6
       x_step=(x2-x1)/6.0_dpreal

       call gauleg(N,x1,x1+x_step,X(1:N),W(1:N))
       call gauleg(N,x1+1.*x_step,x1+2.*x_step,X(1*N+1:2*N),W(1*N+1:2*N))
       call gauleg(N,x1+2.*x_step,x1+3.*x_step,X(2*N+1:3*N),W(2*N+1:3*N))
       call gauleg(N,x1+3.*x_step,x1+4.*x_step,X(3*N+1:4*N),W(3*N+1:4*N))
       call gauleg(N,x1+4.*x_step,x1+5.*x_step,X(4*N+1:5*N),W(4*N+1:5*N))
       call gauleg(N,x1+5.*x_step,x2,X(5*N+1:N6),W(5*N+1:N6))

       end subroutine



!    this routine only print out info on the version number


C>     this subroutine calculates standard gauss Legendre points
C>     between x1 and x2 (usually -1.0_dpreal and  1.0_dpreal)
C>     N is the number of mesh points required.
C>     The grid and the weights are stored in the arrays X and W
C>     @param[in] x1 lower boundary
C>     @param[in] x2 upper boundary
C>     @param[in] N number of grid points
C>     @param[out] X grid points
C>     @param[out] W integration weights
       SUBROUTINE gauleg(N,x1,x2,X,W)
        IMPLICIT NONE
        INTEGER N
        REAL(dpreal) x1,x2,X(N),W(N)
        REAL(dpreal) z1,z,xm,xl,pp,p3,p2,p1,pi,tol
        INTEGER m,i,j

        pi=acos(-1.0)
        tol=1.E-12

        m=(n+1)/2
        xm=0.5*(x2+x1)
        xl=0.5*(x2-x1)

        DO 10 i=1,m
         z=cos(pi*(i-0.25)/(N+0.5))

 20      CONTINUE
         p1=1.0E0
         p2=0.0E0
         DO 30 j=1,N
          p3=p2
          p2=p1
          p1=((2*j-1)*z*p2-(j-1)*p3)/j
 30      CONTINUE
         pp=N*(z*p1-p2)/(z*z-1.0E0)
         z1=z
         z=z1-p1/pp
         IF( abs(z1-z) .GT. tol) GOTO 20 ! Scheifenende

         X(i) = xm - xl*z
         X(n+1-i) = xm + xl*z
         W(i) = 2.E0*xl/((1.0-z*z)*pp*pp)
         W(n+1-i) = W(i)
 10     CONTINUE
       END SUBROUTINE gauleg





C     the subsroutine trns uses gauleg to obtain Gauss-Legendre grid points
C     and performs a hyperbolic transformation on the points and weights
C     resulting new points and weights that are more even spread in the intervals
C     The grid has NP1/2 points between 0 and P1, NP1/2 points between P1 and P2
C     and NP2 points between P2 and P3
C     The complete grid runs from 0 to P3 and has NP=NP1+NP2 grid points
C     It is your responsibility to make sure that NP=NP1+NP2
C     P1,P2,P3 defines the intervals
C     NP1,NP2,NP=NP1+NP2 the number of mesh points.
C     grid and weights are stored in XP and AP on exit


       SUBROUTINE TRNS(NP1,NP2,NP,P1,P2,P3,XP,AP)
       IMPLICIT NONE
C     ===============
C
C     TRNS BELEGT DIE FELDER XP UND AP MIT TRANSFORMIERTEN
C     GAUSS-LEGENDRE-PUNKTEN UND GEWICHTEN
C
C     NP1 PUNKTE WERDEN UEBER DIE HYPERBOLISCHE TRANSFORMATION
C
C     X --> (1.+X) / (1./P1-(1./P1-2./P2)*X)
C
C     AUF DAS INTERVALL (0.;P2) ABGEBILDET, WOBEI
C     NP1/2 PUNKTE IN (0.;P1) UND
C     NP1/2 PUNKTE IN (P1;P2) LIEGEN
C
C     NP2 PUNKTE WERDEN UEBER DIE LINEARE TRANSFORMATION
C
C     X --> (P3+P2)/2. + (P3-P2)/2.*X
C
C     AUF DAS INTERVALL (P2;P3) ABGEBILDET
C
C     NP = NP1 + NP2
C
       REAL(dpreal) XP1(NP1),AP1(NP1),XP2(NP2),AP2(NP2)
       REAL(dpreal) XP(NP),AP(NP)
       REAL(dpreal) P1,P2,P3
       REAL(dpreal),PARAMETER :: eins=1.0_dpreal
       INTEGER NP1,NP2,NP
       REAL(dpreal) XX,X,A,DELPH
       INTEGER I

       CALL gauleg(NP1,-eins,eins,XP1,AP1)

       DO 1 I=1,NP1
        X=XP1(I)
        A=AP1(I)
        XX=1.0_dpreal/P1-(1.0_dpreal/P1-2.0_dpreal/P2)*X
        XP1(I)=(1.0_dpreal+X) / XX
 1      AP1(I)=(2.0_dpreal/P1-2.0_dpreal/P2)*A / XX**2
C
       IF(NP2 .NE. 0) THEN

        CALL gauleg(NP2,-eins,eins,XP2,AP2)


        DO 2 I=1,NP2
         X=XP2(I)
         A=AP2(I)
         DELPH=(P3-P2)/2.0_dpreal
         XP2(I)=(P3+P2)/2.0_dpreal + DELPH*X
 2       AP2(I)=DELPH*A
       ENDIF
C
       DO 3 I=1,NP1
        XP(I)=XP1(I)
 3      AP(I)=AP1(I)
C
       IF(NP2 .NE. 0) THEN
        DO 4 I=1,NP2
         XP(I+NP1)=XP2(I)
 4       AP(I+NP1)=AP2(I)
       ENDIF
C
       RETURN
      END SUBROUTINE trns

C  same as TRNS except that the first interval starts not at zero but at PMIN

      SUBROUTINE TRNSMIN(NP1,NP2,NP,PMIN,P1,P2,P3,XP,AP)
       IMPLICIT NONE
       REAL(dpreal) XP1(NP1),AP1(NP1),XP2(NP2),AP2(NP2)
       REAL(dpreal) XP(NP),AP(NP)
       REAL(dpreal) P1,P2,P3,PMIN
       REAL(dpreal),PARAMETER :: eins=1.0_dpreal
       INTEGER NP1,NP2,NP,I
       REAL(dpreal) AA,BB,CC,X,A,XX1,XX2,DELPH

C     =============
C
C     TRNS BELEGT DIE FELDER XP UND AP MIT TRANSFORMIERTEN
C     GAUSS-LEGENDRE-PUNKTEN UND GEWICHTEN
C
C     NP1 PUNKTE WERDEN UEBER DIE HYPERBOLISCHE TRANSFORMATION
C
C     X --> (A+B*X) / (1.+C*X)
C     AA=P1 ; BB=( P1*(PMIN+P2)-2.*PMIN*P2 )/( P2-PMIN ) ;
C     CC=( 2.*P1-PMIN-P2 )/( P2-PMIN )
C
C     AUF DAS INTERVALL (PMIN;P2) ABGEBILDET, WOBEI
C     NP1/2 PUNKTE IN (PMIN;P1) UND
C     NP1/2 PUNKTE IN (P1;P2) LIEGEN
C
C     NP2 PUNKTE WERDEN UEBER DIE LINEARE TRANSFORMATION
C
C     X --> (P3+P2)/2. + (P3-P2)/2.*X
C
C     AUF DAS INTERVALL (P2;P3) ABGEBILDET
C
C     NP = NP1 + NP2
C
C
      AA=P1
      BB=( P1*(PMIN+P2)-2.*PMIN*P2 )/( P2-PMIN )
      CC=( 2.*P1-PMIN-P2 )/( P2-PMIN )
      CALL gauleg(NP1,-eins,eins,XP1,AP1)
      DO 1 I=1,NP1
      X=XP1(I)
      A=AP1(I)
      XX1=AA+BB*X
      XX2=1.+CC*X
      XP1(I)=XX1/XX2
    1 AP1(I)=( BB*XX2-XX1*CC )*A/(XX2*XX2)
      IF(NP2 .NE. 0) THEN
      CALL gauleg(NP2,-eins,eins,XP2,AP2)
      DO 2 I=1,NP2
      X=XP2(I)
      A=AP2(I)
      DELPH=(P3-P2)/2.
      XP2(I)=(P3+P2)/2. + DELPH*X
    2 AP2(I)=DELPH*A
      ENDIF
C
      DO 3 I=1,NP1
      XP(I)=XP1(I)
    3 AP(I)=AP1(I)
C                QQ
      IF(NP2 .NE. 0) THEN
      DO 4 I=1,NP2
      XP(I+NP1)=XP2(I)
    4 AP(I+NP1)=AP2(I)
      ENDIF
C
      RETURN
      END SUBROUTINE TRNSMIN






      END MODULE gauss
