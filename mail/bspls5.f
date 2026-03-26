      SUBROUTINE BSPLS5( I , KDEGRE , LR , R , LT , T ,
     %                   NBPOIN , POINT ,
     %                   B , A , S , FACM ,
     %                   NORT )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    GENERER LES COEFFICIENTS DES 3 POLYNOMES SUR L'INTERVALLE
C -----    [R(I),R(I+1)[  D'UNE LIGNE B-SPLINE
C
C ENTREES:
C --------
C I      : NUMERO DE 0 A LR-1 DE L'INTERVALLE EN X TRAITE
C KDEGRE : DEGRE DES POLYNOMES DE LA B-SPLINE
C LR     : NOMBRE-1 D'INTERVALLES DE LA LIGNE B-SPLINE
C LT     : NOMBRE DE POINTS-1 DE CONTROLE (T) DE LA LIGNE
C R      : EXTREMITES DES INTERVALLES
C T      : VALEURS DU PARAMETRE POUR CALCULER LES B(J,M)
C NBPOIN : NOMBRE TOTAL DE POINTS DE CONTROLE DE LA SURFACE
C POINT  : TABLEAU DE TOUS LES POINTS DE CONTROLE DE LA SURFACE
C          OBTENU APRES INVERSION DU SYSTEME LINEAIRE
C
C AUXILIAIRES :
C -------------
C B      : VALEURS DES B(J,M)(R(I))
C A      : VALEURS INTERMEDIAIRES
C FACM   : 1/M!
C NORT   : NUMERO DU DERNIER NOEUD T DE CHAQUE POINT R
C
C SORTIES:
C --------
C R      : LE PARAMETRE DES POINTS ASSOCIES A T
C S      : LES COEFFICIENTS DES POLYNOMES SUR CHAQUE INTERVALLE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS       JUIN  1990
C2345X7..............................................................012
      INTEGER  NORT(0:LR)
      REAL     R(0:LR),
     %         T(0:LT+KDEGRE),
     %         POINT(1:NBPOIN,1:3),
     %         B(-KDEGRE:KDEGRE,0:KDEGRE),
     %         A(-KDEGRE:0,0:KDEGRE),
     %         FACM(0:KDEGRE),
     %         S(0:KDEGRE,1:3)
C
C     GENERATION DU POLYNOME DE LA B-SPLINE
C     =====================================
C        LA VALEUR DU PARAMETRE
         SS  = R(I)
C        LE NUMERO DU DERNIER NOEUD DU POINT I
         NOI = NORT( I )
C
C        CALCUL DES BJM(RI)
         DO 85 J=-KDEGRE,KDEGRE
            B(J,0) = 0
 85      CONTINUE
         B(0,0) = 1.0
C
         DO 90 M=1,KDEGRE
            DO 88 J=-KDEGRE,KDEGRE-M
C              LE VRAI INDICE DE B
               JB = NOI + J
C              EVALUATION DES FRACTIONS DE LA FORMULE
               IF( T(JB+M) .NE. T(JB) ) THEN
                  U1 = ( SS - T(JB) ) / ( T(JB+M) - T(JB) )
               ELSE
                  U1 = 0
               ENDIF
               IF( T(JB+M+1) .NE. T(JB+1) ) THEN
                  U2 = ( T(JB+M+1) - SS )
     %               / ( T(JB+M+1) - T(JB+1) )
               ELSE
                  U2 = 0
               ENDIF
               B(J,M) = U1 * B(J,M-1) + U2 * B(J+1,M-1)
 88         CONTINUE
 90      CONTINUE
C
C        CALCULS POUR CHAQUE COMPOSANTE DES POINTS
         DO 900 NC=1,3
C
C           INITIALISATION
            DO 120 J=-KDEGRE,0
               A(J,0) =POINT( NOI + J , NC )
 120        CONTINUE
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
               SS = 0
               DO 310 J=-KDEGRE+M,0
                  SS = SS + A(J,M) * B(J,KDEGRE-M)
 310           CONTINUE
               S(M,NC) = FACM(M) * SS
 320        CONTINUE
C
 900     CONTINUE
      END
