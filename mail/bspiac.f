      SUBROUTINE BSPIAC( LIGFER , KDEGRE , NBPCBL , LR ,
     %                   NBSOLI , XYZPC , R , U , S , XYZSOM ,
     %                   XYZTGS , NBTGS )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LES 3 COORDONNEES DES NBSOLI-1 SOMMETS DE LA LIGNE
C -----    ET LES 3 COMPOSANTES DES 2 TANGENTES DES ARETES
C          D'UNE B-SPLINE POLYNOMIALE D'INTERPOLATION
C
C ENTREES:
C --------
C LIGFER : 1 LIGNE FERMEE , 0 SI LIGNE OUVERTE
C KDEGRE : DEGRE DES POLYNOMES DE LA B-SPLINE
C NBPCBL : NOMBRE DE POINTS D'INTERPOLATION
C LR     : NOMBRE DE POINTS-1 DE CONTROLE (R) DE LA LIGNE
C NBSOLI : NOMBRE DE SOMMETS DE LA LIGNE
C XYZPC  : LES COORDONNEES DES L+1 POINTS DE CONTROLE
C R      : LES ABSCISSES PARAMETRE AYANT POUR IMAGE LES POINTS CONTROLE
C U      : LES ABSCISSES PARAMETRE DES NOEUDS D'INTERPOLATION
C S      : LES COEFFICIENTS DES POLYNOMES SUR CHAQUE INTERVALLE
C NBTGS  : NOMBRE DE TANGENTES DE LA LIGNE
C
C SORTIES:
C --------
C XYZSOM : 3 COORDONNEES DES POINTS DE CONTROLE
C XYZTGS : 3 COORDONNEES DES TANGENTES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS       MAI   1996
C2345X7..............................................................012
      REAL     R(0:LR),U(0:NBPCBL),
     %         XYZPC(1:3,1:NBPCBL),
     %         S(0:KDEGRE,0:LR-1,1:3),
     %         XYZSOM(1:3,1:NBSOLI),
     %         XYZTGS(1:3,1:NBTGS)
      REAL     TG(3)
C
      NBARLI = NBSOLI - 1
      LU     = NBPCBL - 1
C
C     LE PARAMETRE DU 1-ER POINT
      XYZSOM(3,1) = U(0)
C
C     SOMMETS EQUIDISTANTS DANS CHAQUE INTERVALLE DU PARAMETRE
C     ========================================================
      NBA = NBARLI / LU
      N   = 1
      DO 20 J=1,LU
         PAS = ( U(J) - U(J-1) ) / NBA
         DO 10 I=1,NBA
            N = N + 1
            XYZSOM(3,N) = XYZSOM(3,N-1) + PAS
 10      CONTINUE
 20   CONTINUE
C
C     LE PARAMETRE DU DERNIER POINT
      XYZSOM(3,NBSOLI) = U(LU)
C
C     LA LONGUEUR DU PREMIER INTERVALLE
      H00 = XYZSOM(3,2) - XYZSOM(3,1)
C
C     PASSAGE AUX 3 COORDONNEES DES SOMMETS DE LA LIGNE DANS R**3
C     ===========================================================
      NUTG = 0
      H1   = 0
      I    = 0
      DO 100 N=1,NBARLI
C
C        LA VALEUR DU PARAMETRE
         PAS = XYZSOM(3,N)
         H0  = H1
         H1  = XYZSOM(3,N+1) - XYZSOM(3,N)
C
 50      IF( PAS .GE. R(I+1) ) THEN
C
C           PASSAGE A L'INTERVALLE SUIVANT DE R
            I = I + 1
            GOTO 50
C
         ELSE
C
C           PAS EST DANS L'INTERVALLE R(I) R(I+1)
            RR = PAS - R(I)
C
C           CALCUL DES 3 COORDONNEES DU SOMMET 1 DE L'ARETE N
            DO 80 J=1,3
               A = S(KDEGRE,I,J)
               DO 70 M=KDEGRE-1,0,-1
                  A = A * RR + S(M,I,J)
 70            CONTINUE
               XYZSOM(J,N) = A
 80         CONTINUE
C
            IF( NBTGS .GT. 0 ) THEN
C              CALCUL DES 3 COORDONNEES DE LA TANGENTE DU SOMMET 1 DE L'ARETE N
               NUTG = NUTG + 1
               DO 90 J=1,3
                  A = KDEGRE * S(KDEGRE,I,J)
                  DO 95 M=KDEGRE-1,1,-1
                     A = A * RR + M * S(M,I,J)
 95               CONTINUE
                  TG(J) = A
                  XYZTGS(J,NUTG) = A * H1
 90            CONTINUE
C
               IF( N .GT. 1 ) THEN
C                 STOCKAGE DE LA TANGENTE DU SOMMET 2 DE L'ARETE PRECEDENTE
                  DO 69 J=1,3
                     XYZTGS(J,NUTG-1) = -TG(J) * H0
 69               CONTINUE
               ENDIF
C
C              INCREMENTATION DE LA TANGENTE A CALCULER ENSUITE
               NUTG = NUTG + 1
            ENDIF
C
         ENDIF
C
 100  CONTINUE
C
      IF ( LIGFER .EQ. 0 ) THEN
C
C        LIGNE OUVERTE
C        -------------
C        CALCUL DES COORDONNEES DE LA TANGENTE RELATIVE AU DERNIER SOMMET
         IF ( NBTGS .GT. 0 ) THEN
C           LA VALEUR DU PARAMETRE AU POINT NBSOLI DE LA B-SPLINE
C           XYZSOM(3,NBSOLI) EST DANS L'INTERVALLE R(I) R(I+1)
            I  = LR - 1
            RR = XYZSOM(3,NBSOLI) - R(I)
            DO 300 J=1,3
               A = KDEGRE * S(KDEGRE,I,J)
               DO 305 M=KDEGRE-1,1,-1
                  A = A * RR + M * S(M,I,J)
 305           CONTINUE
               XYZTGS(J,NBTGS) = -A * H1
 300        CONTINUE
         ENDIF
C
C        LE PREMIER ET DERNIER SOMMET DE LA LIGNE
         DO 200 J=1,3
            XYZSOM( J ,      1 ) = XYZPC( J ,    1   )
            XYZSOM( J , NBSOLI ) = XYZPC( J , NBPCBL )
 200     CONTINUE
C
      ELSE
C
C        LIGNE FERMEE
C        ------------
         IF( NBTGS .GT. 0 ) THEN
            DO 310 J=1,3
               XYZTGS(J,NUTG) = -(XYZTGS(J,1) / H00) * H1
 310        CONTINUE
         ENDIF
C
C        DERNIER=PREMIER SOMMET
         DO 320 J=1,3
            XYZSOM( J , NBSOLI ) = XYZSOM( J , 1 )
            XYZTGS( J , NBTGS  ) = XYZTGS( J , 1 )
 320     CONTINUE
C
      ENDIF
      END
