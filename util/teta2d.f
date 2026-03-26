      DOUBLE PRECISION FUNCTION TETA2D( X , Y )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCUL DE TETA 2-EME COORDONNEE CYLINDRIQUE EN FONCTION DE X , Y
C ----- COORDONNEES CARTESIENNES 2D
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS        MARS 1982
C2345X7..............................................................012
      DOUBLE PRECISION PI,X,Y,DATAN
      DATA             PI/ 3.14159265358979312D0 /

      IF( X .LT. 0D0 ) GOTO 30
      IF( X .EQ. 0D0 ) GOTO 20

C     X > 0
C     =====
      IF( Y .LT. 0D0 ) GOTO 16
C
C     X > 0 , Y>= 0
C     -------------
      TETA2D = DATAN( Y / X )
      GOTO 9999
C
C     X > 0 , Y < 0
C     -------------
   16 TETA2D = 2.D0 * PI + DATAN( Y / X )
      GOTO 9999
C
C     X = 0
C     =====
   20 IF( Y .LT. 0D0 ) GOTO 26
      IF( Y .EQ. 0D0 ) GOTO 24

C     X = 0 , Y > 0
C     -------------
      TETA2D = 0.5D0 * PI
      GOTO 9999
C
C     X = Y + 0
C     ---------
   24 TETA2D = 0.D0
      GOTO 9999
C
C     X = 0 , Y < 0
C     -------------
   26 TETA2D = 1.5D0 * PI
      GOTO 9999
C
C     X <= 0
C     ======
   30 TETA2D = PI + DATAN( Y / X )

 9999 RETURN
      END
