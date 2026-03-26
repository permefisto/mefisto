      SUBROUTINE BSPLS4( I  , DEGREX , LR , R , LT , T ,
     %                   NORTY , DEGREY , NBPX , NBPY , POINT ,
     %                   B  , A , FACM , NORXTX , S )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    GENERER LES COEFFICIENTS DES 3 POLYNOMES SUR L'INTERVALLE
C -----    [R(I),R(I+1)[  D'UNE LIGNE B-SPLINE
C
C ENTREES:
C --------
C I      : NUMERO DE 0 A LR-1 DE L'INTERVALLE EN X TRAITE
C DEGREX : DEGRE DES POLYNOMES DE LA B-SPLINE
C LR     : NOMBRE-1 D'INTERVALLES DE LA LIGNE B-SPLINE
C R      : EXTREMITES DES INTERVALLES
C LT     : NOMBRE DE POINTS-1 DE CONTROLE (T) DE LA LIGNE
C T      : VALEURS DU PARAMETRE POUR CALCULER LES B(J,M)
C NORTY  : NUMERO DU DERNIER NOEUD EN Y DE RY(JY)
C DEGREY : DEGRE EN Y DES POLYNOMES
C NBPX   : NOMBRE DE POINTS EN X
C NBPY   : NOMBRE DE POINTS EN Y
C POINT  : TABLEAU DE TOUS LES POINTS DE CONTROLE DE LA SURFACE
C          OBTENU APRES INVERSION DU SYSTEME LINEAIRE
C
C AUXILIAIRES :
C -------------
C B      : VALEURS DES B(J,M)(R(I))
C AY     : VALEURS INTERMEDIAIRES
C FACM   : 1/M!
C NORXTX : NUMERO DU DERNIER NOEUD T DE CHAQUE POINT R
C
C SORTIES:
C --------
C S      : LES COEFFICIENTS DES POLYNOMES EN X SUR LES -DEGREY NOEUDS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS       JUIN  1990
C2345X7..............................................................012
      INTEGER  NORXTX(0:LR),DEGREX,DEGREY
      REAL     R(0:LR),
     %         T(0:LT+DEGREX),
     %         B(-DEGREX:1,0:DEGREX),
     %         A(-DEGREX:0,0:DEGREX),
     %         FACM(0:DEGREX),
     %         S(0:DEGREX,1:3,-DEGREY:0)
      DOUBLE   PRECISION POINT(1:NBPX,1:NBPY,1:3),SS
C
C     GENERATION DU POLYNOME DE LA B-SPLINE
C     =====================================
C        LA VALEUR DU PARAMETRE
         SS  = R(I)
C        LE NUMERO DU DERNIER NOEUD DU POINT I
         NORTX = NORXTX( I )
C
C        CALCUL DES BJM(RI)
C        ------------------
         CALL AZEROR( (DEGREX+2)*(DEGREX+1) , B )
         B(0,0) = 1.0
         DO 20 M=1,DEGREX
            DO 10 J=-M,0
C              LE VRAI INDICE DE B
               JB = NORTX + J
C              EVALUATION DES FRACTIONS DE LA FORMULE
               IF( T(JB+M) .NE. T(JB) ) THEN
                  U1 = REAL( ( SS - T(JB) ) / ( T(JB+M) - T(JB) ) )
               ELSE
                  U1 = 0
               ENDIF
               KB = JB + M + 1
               IF( T(KB) .NE. T(JB+1) ) THEN
                  U2 = REAL( ( T(KB) - SS )
     %                     / ( T(KB) - T(JB+1) ) )
               ELSE
                  U2 = 0
               ENDIF
               B(J,M) = U1 * B(J,M-1) + U2 * B(J+1,M-1)
 10         CONTINUE
 20      CONTINUE
C
C        LES POLYNOMES POUR IY=JY-DEGREY A JY
C        ------------------------------------
         DO 300 IY=-DEGREY,0
C
C           CALCULS POUR CHAQUE COMPOSANTE DES POINTS
            DO 200 NC=1,3
C
C              INITIALISATION
               DO 120 J=-DEGREX,0
                  A(J,0) = REAL(POINT( 1 +NORTX +J, 1 +NORTY +IY , NC ))
 120           CONTINUE
C
C              LE CALCUL PROPREMENT DIT
               DO 140 M=1,DEGREX
                  KB = DEGREX - M + 1
                  DO 130 J=-DEGREX+M,0
C                    LE VRAI INDICE DE A
                     JB = NORTX + J
                     IF( T(JB+KB) .NE. T(JB) ) THEN
                        A(J,M) = KB * (A(J,M-1)-A(J-1,M-1))
     %                              / (T(JB+KB)-T(JB))
                     ELSE
                        A(J,M) = 0
                     ENDIF
 130              CONTINUE
 140           CONTINUE
C
C              LES COEFFICIENTS S(M,I,NC) DU POLYNOME
               DO 160 M=0,DEGREX
                  SS = 0
                  DO 150 J=-DEGREX+M,0
                     SS = SS + A(J,M) * B(J,DEGREX-M)
 150              CONTINUE
                  S(M,NC,IY) = REAL( FACM(M) * SS )
 160           CONTINUE
 200        CONTINUE
 300     CONTINUE
C
      RETURN
      END
