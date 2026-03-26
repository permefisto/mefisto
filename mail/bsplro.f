      SUBROUTINE BSPLRO( NBCOMP , KDEGRE , LR , L , NBNOBR , NBINBS ,
     %                   R , POIDS , T , B , XYZ , A , S , FACM ,
     %                   NO )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    GENERER LES COEFFICIENTS DES POLYNOMES SUR CHAQUE INTERVALLE
C -----    D'UNE B-SPLINE RATIONNELLE OU NON OUVERTE NON UNIFORME
C
C ENTREES:
C --------
C NBCOMP : NOMBRE DE COMPOSANTES (3 BSPLINE POLYNOMIALE,4 RATIONNELLE)
C KDEGRE : DEGRE DES POLYNOMES DE LA B-SPLINE
C LR     : NOMBRE-1 D'EXTREMITES DIFFERENTES DANS T
C L      : NOMBRE DE POINTS-1 DE CONTROLE (R) DE LA LIGNE
C NBNOBR : NOMBRE DE NOEUDS T
C R      : LES ABSCISSES PARAMETRE IDENTIFIES
C POIDS  : POIDS DES POINTS DE CONTROLE
C XYZ    : 3 COORDONNEES DES POINTS DE CONTROLE
C
C AUXILIAIRES :
C -------------
C T      : VALEURS DU PARAMETRES AUX NOEUDS
C B      : VALEURS DES B(J,M)(RI)
C A      : VALEURS INTERMEDIAIRES
C FACM   : 1/M!
C NO     : NUMERO DU DERNIER NOEUD DE CHAQUE POINT
C
C SORTIES:
C --------
C S      : LES COEFFICIENTS DES POLYNOMES SUR CHAQUE INTERVALLE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS       JUIN  1990
C2345X7..............................................................012
      INTEGER  NO(0:LR)
      REAL     R (0:LR),POIDS(0:L-1),
     %         T(0:NBNOBR-1),
     %         B(-KDEGRE:KDEGRE,0:KDEGRE),
     %         A(-KDEGRE:0,0:KDEGRE),
     %         XYZ(1:3,0:L-1),
     %         FACM(0:KDEGRE),
     %         S(0:KDEGRE,0:NBINBS-1,1:NBCOMP)
      DOUBLE PRECISION F
C
C     LE TABLEAU DES 1/M!
      FACM(0) = 1.
      FACM(1) = 1.
      F       = 1.D0
      DO 10 M=2,KDEGRE
        F = F * M
        FACM(M) = REAL( 1.D0 / F )
 10   CONTINUE
C
C     LA BOUCLE SUR LES INTERVALLES R(I)-R(I+1)
C     =========================================
      DO 1000 I=0,LR-1
C        LE NUMERO DU DERNIER NOEUD DU POINT I
         NOI = NO( I )
C
C        CALCUL DES BJM(RI)
         DO 45 J=-KDEGRE,KDEGRE
            B(J,0) = 0
 45      CONTINUE
         B(0,0) = 1.0
C
         DO 50 M=1,KDEGRE
            DO 48 J=-KDEGRE,KDEGRE-M
C              LE VRAI INDICE DE B
               JA = NOI + J
C              EVALUATION DES FRACTIONS DE LA FORMULE
               IF( T(JA+M) .NE. T(JA) ) THEN
                  U1 = ( R(I) - T(JA) ) / ( T(JA+M) - T(JA) )
               ELSE
                  U1 = 0
               ENDIF
               KM = JA + M + 1
               IF( T(KM) .NE. T(JA+1) ) THEN
                  U2 = ( T(KM) - R(I) ) / ( T(KM) - T(JA+1) )
               ELSE
                  U2 = 0
               ENDIF
               B(J,M) = U1 * B(J,M-1) + U2 * B(J+1,M-1)
 48         CONTINUE
 50      CONTINUE
C
C        CALCULS POUR CHAQUE COMPOSANTE DES POINTS
         DO 900 NC=1,NBCOMP
C
C           INITIALISATION
            IF ( NBCOMP .EQ. 4 ) THEN
C              RATIONNELLE
               DO 100 J=-KDEGRE,0
                  JA = NOI + J
                  IF( NC .LE. 3 ) THEN
                     A(J,0) = POIDS(JA) * XYZ( NC , JA )
                  ELSE
                     A(J,0) = POIDS(JA)
                  ENDIF
 100           CONTINUE
            ELSE
C              POLYNOMIALE
               DO 120 J=-KDEGRE,0
                  A(J,0) = XYZ( NC , NOI + J )
 120           CONTINUE
            ENDIF
C
C           LE CALCUL PROPREMENT DIT
            DO 240 M=1,KDEGRE
               KM = KDEGRE - M + 1
               DO 230 J=-KDEGRE+M,0
C                 LE VRAI INDICE DE A
                  JA = NOI + J
                  IF( T(JA+KM) .NE. T(JA) ) THEN
                     A(J,M) = KM*(A(J,M-1)-A(J-1,M-1))/(T(JA+KM)-T(JA))
                  ELSE
                     A(J,M) = 0
                  ENDIF
 230           CONTINUE
 240        CONTINUE
C
C           LES COEFFICIENTS S(M,I,NC) DU POLYNOME
            DO 320 M=0,KDEGRE
               U1 = 0
               DO 310 J=-KDEGRE+M,0
                  U1 = U1 + A(J,M) * B(J,KDEGRE-M)
 310           CONTINUE
               S(M,I,NC) = FACM(M) * U1
 320        CONTINUE
C
 900     CONTINUE
 1000 CONTINUE
      END
