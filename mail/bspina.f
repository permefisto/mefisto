      SUBROUTINE BSPINA( LIGFER , KDEGRE , NBPCBL , LR ,
     %                   NBSOLI , RAGEBL , XYZPC , R , S , XYZSOM ,
     %                   XYZTGS , NBTGS )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LES 3 COORDONNEES DES NBSOLI-1 SOMMETS DE LA LIGNE
C -----    ET LES 2 TANGENTES DE CHAQUE ARETE D'UNE COURBE B-SPLINE
C          POLYNOMIALE D'INTERPOLATION
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
C S      : LES COEFFICIENTS DES POLYNOMES SUR CHAQUE INTERVALLE
C NBTGS  : NOMBRE DE TANGENTES DE LA LIGNE (2 FOIS LE NOMBRE D'ARETES)
C
C SORTIES:
C --------
C XYZSOM : 3 COORDONNEES DES POINTS DE CONTROLE
C XYZTGS : 3 COORDONNEES DES TANGENTES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS       MAI   1996
C2345X7..............................................................012
      REAL     R(0:LR),
     %         XYZPC(1:3,1:NBPCBL),
     %         S(0:KDEGRE,0:LR-1,1:3),
     %         XYZSOM(1:3,1:NBSOLI),
     %         XYZTGS(1:3,1:NBTGS)
      REAL     TG(3)
C
C     LE PARAMETRE DU 1-ER ET DERNIER POINT
      IF( LIGFER .EQ. 1 .AND. MOD(KDEGRE,2) .EQ. 0 ) THEN
         U0 = KDEGRE / 2 + 0.5
      ELSE
         U0 = KDEGRE + 1
      ENDIF
      IF( LIGFER .EQ. 0 ) THEN
         U0 = R(0)
      ENDIF
      XYZSOM(3,1) = U0
C
C     PROGRESSION GEOMETRIQUE OU NON ?
      NBS1 = NBSOLI - 1
      IF( ABS(RAGEBL-1.) .LT. 1.E-3 ) THEN
C
C        SOMMETS EQUIDISTANTS DANS L'ESPACE DU PARAMETRE
C        ===============================================
         RAGEBL = 1.0
         PAS    = ( R(LR) - R(0) ) / NBS1
         DO 20 I=2,NBSOLI
            XYZSOM(3,I) = XYZSOM(3,I-1) + PAS
 20      CONTINUE
C
      ELSE
C
C        SOMMETS EN PROGRESSION GEOMETRIQUE DANS L'ESPACE DU PARAMETRE
C        =============================================================
         PAS = (R(LR)-R(0)) * (1.-RAGEBL) / (1.-RAGEBL**NBS1)
         DO 40 I=2,NBSOLI
            XYZSOM(3,I) = XYZSOM(3,I-1) + PAS
            PAS         = PAS * RAGEBL
 40      CONTINUE
      ENDIF
C     ICI: LA VALEUR DU PARAMETRE DE CHAQUE SOMMET EST RANGE DANS XYZSOM(3,*)!
C
      PASINI = XYZSOM(3,2) - XYZSOM(3,1)
      PASFIN = XYZSOM(3,NBSOLI) - XYZSOM(3,NBS1)
C
C     PASSAGE AUX 3 COORDONNEES DES SOMMETS ET DES TANGENTES DE LA LIGNE DANS R*
C     ==========================================================================
      NUTG   = 0
      H1     = 0.0
      I      = 0
      TRANSL = 0.0
      DO 100 N=1,NBS1
C
C        LA LONGUEUR DU SOUS INTERVALLE N DU PARAMETRE
         H0  = H1
         H1  = XYZSOM(3,N+1) - XYZSOM(3,N)
C        LA VALEUR DU PARAMETRE AU POINT N DE LA B-SPLINE
 50      PAS = XYZSOM(3,N) + TRANSL
         IF( PAS .GE. R(I+1) ) THEN
C           PASSAGE A L'INTERVALLE SUIVANT DE R
            I = I + 1
            IF( LIGFER .EQ. 1 ) THEN
C              LA LIGNE EST FERMEE
               IF( I .GE. LR ) THEN
C                 RETOUR AU DEBUT DU PARAMETRE
                  I      = 0
                  TRANSL = R(0) - R(LR)
               ENDIF
            ENDIF
            GOTO 50
C
         ELSE
C
C           PAS EST DANS L'INTERVALLE R(I) R(I+1)
C           ATTENTION A LA DIFFERENCE PAS - R(I) POUR LE CALCUL DU POLYNOME!
            RR = PAS - R(I)
C
            IF ( NBTGS .GT. 0 ) THEN
C
C              CALCUL DES 3 COMPOSANTES DE LA TG DU SOMMET 1 DE L'ARETE N
               NUTG = NUTG + 1
               DO 68 J=1,3
C                 CALCUL DE LA COMPOSANTE J DE LA DERIVEE PAR RAPPORT A R
                  A = KDEGRE * S(KDEGRE,I,J)
                  DO 64 M=KDEGRE-1,1,-1
                     A = A * RR + M * S(M,I,J)
 64               CONTINUE
                  TG(J) = A
                  XYZTGS(J,NUTG) = A * H1
 68            CONTINUE
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
C           LES 3 COORDONNEES DU SOMMET N DE LA B-SPLINE
            DO 80 J=1,3
               A = S(KDEGRE,I,J)
               DO 70 M=KDEGRE-1,0,-1
                  A = A * RR + S(M,I,J)
 70            CONTINUE
               XYZSOM(J,N) = A
 80         CONTINUE
C
         ENDIF
 100  CONTINUE
C
      IF ( LIGFER .EQ. 0 ) THEN
C
C        LIGNE OUVERTE
C        -------------
         IF( NBTGS .GT. 0 ) THEN
C           LE CALCUL DE LA DERNIERE TANGENTE DE LA B-SPLINE
C           LA VALEUR DU PARAMETRE AU POINT NBSOLI DE LA B-SPLINE
            PAS = XYZSOM(3,NBSOLI)
C           PAS EST DANS L'INTERVALLE R(I) R(I+1)
C           ATTENTION A LA DIFFERENCE PAS - R(I) POUR LE CALCUL DU POLYNOME!
            RR = PAS - R(I)
            DO 270 J=1,3
               A = KDEGRE * S(KDEGRE,I,J)
               DO 260 M=KDEGRE-1,1,-1
                  A = A * RR + M * S(M,I,J)
 260           CONTINUE
               XYZTGS(J,NUTG) = - A * PASFIN
 270        CONTINUE
         ENDIF
C
C        LE PREMIER ET DERNIER SOMMET DE LA LIGNE OUVERTE
         DO 280 J=1,3
            XYZSOM( J,      1 ) = XYZPC( J,    1   )
            XYZSOM( J, NBSOLI ) = XYZPC( J, NBPCBL )
 280     CONTINUE
C
      ELSE
C
C        LIGNE FERMEE
C        ------------
         IF( NBTGS .GT. 0 ) THEN
            DO 310 J=1,3
               XYZTGS(J,NUTG) = -(XYZTGS(J,1) / PASINI) * PASFIN
 310        CONTINUE
            IF( NUTG .NE. NBTGS ) CALL XVPAUSE
         ENDIF
C
C        LE DERNIER EST LE MEME QUE LE PREMIER SOMMET
         DO 320 J=1,3
            XYZSOM( J, NBSOLI ) = XYZSOM( J, 1 )
 320     CONTINUE
C
      ENDIF
      END
