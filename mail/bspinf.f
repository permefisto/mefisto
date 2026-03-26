      SUBROUTINE BSPINF( KDEGRE , LU , XYZ ,
     %                   B , A , FACM , MATRIX ,
     %                   S , IERR )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    GENERER LES COEFFICIENTS DES POLYNOMES SUR CHAQUE INTERVALLE
C -----    D'UNE B-SPLINE FERMEE UNIFORME D'INTERPOLATION
C
C ENTREES:
C --------
C KDEGRE : DEGRE DES POLYNOMES DE LA B-SPLINE
C LU     : NOMBRE DE POINTS D'INTERPOLATION
C XYZ    : 3 COORDONNEES DES POINTS DE CONTROLE
C
C AUXILIAIRES :
C -------------
C B      : VALEURS DES B(J,M)(TI)
C A      : VALEURS INTERMEDIAIRES
C FACM   : 1/M!
C MATRIX   : LA MATRICE DES BJ,M(RI) A INVERSER
C
C SORTIES:
C --------
C S      : LES COEFFICIENTS DES POLYNOMES SUR CHAQUE INTERVALLE
C IERR   : =0 SI PAS D'ERREUR
C          >0 SI ERREUR RENCONTREE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS       JUIN  1990
C2345X7..............................................................012
      include"./incl/gsmenu.inc"
      REAL     B(-KDEGRE:1,0:KDEGRE),
     %         A(-KDEGRE:0,0:KDEGRE),
     %         XYZ(1:3,0:LU-1),
     %         FACM(0:KDEGRE),
     %         S(0:KDEGRE,0:LU-1,1:3)
      DOUBLE PRECISION MATRIX(LU,LU+3), SS
C
C     LE TABLEAU DES 1/M!
      FACM(0) = 1.0
      FACM(1) = 1.0
      SS      = 1.D0
      DO 10 M=2,KDEGRE
        SS = SS * M
        FACM(M) = REAL( 1.D0 / SS )
 10   CONTINUE
C
C     LA VALEUR DE U0
      KS2 = KDEGRE / 2
      IF( MOD(KDEGRE,2) .EQ. 0 ) THEN
C        KDEGRE=2 KS2
         IMPAIR = 0
         U0     = KS2 + 0.5
         NOU0   = KS2
      ELSE
C        KDEGRE=2 KS2 + 1
         IMPAIR = 1
         U0     = KS2 + 1
         NOU0   = KS2 + 1
      ENDIF
C
C     CALCUL DES BJM(U0) POUR J=-KDEGRE,...,0
C     =======================================
      LB = (KDEGRE+2) * (KDEGRE+1)
      CALL AZEROR( LB , B )
      B(0,0) = 1.0
      DO 50 M=1,KDEGRE
         DO 48 J=-M,0
C           LE VRAI INDICE DE B
            JB = NOU0 + J
C           EVALUATION DES FRACTIONS DE LA FORMULE
            B(J,M) = ( ( U0     - JB ) * B(J,M-1) +
     %                 ( JB+M+1 - U0 ) * B(J+1,M-1) ) / M
 48      CONTINUE
 50   CONTINUE
C
C     CONSTRUCTION DE LA MATRICE MATRIX ET DES 3 SECONDS MEMBRES
C     ========================================================
C     MISE A ZERO DE LA MATRICE
      CALL AZEROD( LU*LU , MATRIX )
C
      DO 60 I=0,LU-1
C
C        RANGEMENT DE B(RI) DANS LA LIGNE I DE LA MATRICE MATRIX
         DO 57 J=-KDEGRE,0
            M = 1 + NOU0 + I + J
            IF( M .LE.  0 ) M = M + LU
            IF( M .GT. LU ) M = M - LU
            MATRIX(1+I,M) = B(J,KDEGRE)
 57      CONTINUE
C
C        GENERATION DES 3 SECONDS MEMBRES
         DO 59 NC=1,3
            MATRIX(1+I,LU+NC) = XYZ(NC,I)
 59      CONTINUE
 60   CONTINUE
C
C     INVERSION DU SYSTEME LINEAIRE
C     =============================
      CALL GAUSPT( LU , 3 , MATRIX , I )
      IF( I .NE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'PROBLEME DANS BSPINF'
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF
C
C     CALCUL DES BJM(0) POUR J=-KDEGRE,...,0   (T(I)=I)
C     ======================================
      IF( IMPAIR .EQ. 0 ) THEN
C        SI KDEGRE EST IMPAIR LES B SONT DEJA CALCULES
         CALL AZEROR( LB , B )
         B(0,0) = 1.0
         DO 80 M=1,KDEGRE
            DO 70 J=-M,0
C              EVALUATION DES FRACTIONS DE LA FORMULE
               B(J,M) = ( ( J+M+1 )  * B(J+1,M-1)
     %                -        J     * B(J,M-1)  ) / M
 70         CONTINUE
 80      CONTINUE
      ENDIF
C
C     GENERATION DES POLYNOMES DE LA B-SPLINE
C     =======================================
C     LE TABLEAU A SERA COMPLETE AU FUR ET A MESURE QUE I CROIT
      DO 90 J=-KDEGRE,0
         A(J,0) = 0
 90   CONTINUE
C
      DO 1000 I=0,LU-1
C
C        CALCULS POUR CHAQUE COMPOSANTE DES POINTS
         DO 900 NC=1,3
C
C           INITIALISATION
            DO 120 J=-KDEGRE,0
               M = 1 + J + I
               IF( M .LE. 0 ) M = M + LU
               A(J,0) = REAL( MATRIX( M , LU + NC ) )
 120        CONTINUE
C
C           LE CALCUL PROPREMENT DIT
            DO 240 M=1,KDEGRE
               DO 230 J=-KDEGRE+M,0
                  A(J,M) = A(J,M-1) - A(J-1,M-1)
 230           CONTINUE
 240        CONTINUE
C
C           LES COEFFICIENTS S(M,I,NC) DU POLYNOME
            DO 320 M=0,KDEGRE
               SS = 0D0
               DO 310 J=-KDEGRE+M,0
                  SS = SS + A(J,M) * B(J,KDEGRE-M)
 310           CONTINUE
               S(M,I,NC) = REAL( FACM(M) * SS )
 320        CONTINUE
C
 900     CONTINUE
 1000 CONTINUE
      END
