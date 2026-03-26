         SUBROUTINE SUEXQ2( NBS1,   NBS2,   NBA1,   NBA2,   NBTGS,
     S                      MNXYT1, MNNTG1, MNXYT2, MNNTG2,
     S                      MNXYT3, MNNTG3, MNXYT4, MNNTG4,
     S                      COSO,   COTG,   STCARR, TGPURE )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     MAILLER UN QUADRANGLE C1 DEFINI PAR SES BORDS
C -----     PAR LA METHODE DE COOK-GORDON-HALL AVEC OU SANS TANGENTES
C
C ENTREES :
C ---------
C NBS1    : NOMBRE DE POINTS SUR LES LIGNES "HORIZONTALES"
C NBS2    : NOMBRE DE POINTS SUR LES LIGNES "VERTICALES"
C NBA1    : NOMBRE D'ARETES  SUR LES LIGNES "HORIZONTALES"
C NBA2    : NOMBRE D'ARETES  SUR LES LIGNES "VERTICALES"
C NBTGS   : EXISTENCE (>0) OU NON (=0) DE TANGENTES
C MNXYT1234: ADRESSE MCN DE LA PREMIERE COMPOSANTE DES TANGENTES DU COTE 1234
C MNNTG1234: ADRESSE MCN DES 2 NUMEROS DES TANGENTES DES ARETES  DU COTE 1234
C COSO    : COORDONNEES DES SOMMETS DES 4 COTES (BORD) DU MAILLAGE
C
C SORTIES :
C ---------
C STCARR  : COORDONNEES DES NBS1 x NBS2 SOMMETS DU CARRE UNITE
C COSO    : COORDONNEES DES NBS1 x NBS2 SOMMETS DU MAILLAGE
C COTG    : 3 COMPOSANTES DES 8 TANGENTES DES (NBA1)(NBA2)QUADRANGLES
C
C TABLEAUX AUXILIAIRES:
C ---------------------
C STCARR  : COORDONNEES DES NBS1 x NBS2 SOMMETS DU CARRE UNITE
C TGPURE  : SI NBTGS>0 LES 2 DERIVEES PURES EN TOUS LES SOMMETS DU CARRE UNITE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY     ANALYSE NUMERIQUE UPMC PARIS   SEPTEMBRE 1988
C MODIFS : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS   SEPTEMBRE 1996
C23456---------------------------------------------------------------012
      COMMON / EPSSSS / EPZERO,EPSXYZ
      REAL              STCARR(2,NBS1,NBS2)
      REAL              COSO(3,NBS1,NBS2)
      REAL              COTG(3,8,NBA1,NBA2)
      REAL              TGPURE(3,2,NBS1,NBS2)
C
C     CALCUL DES 2 COORDONNEES DES SOMMETS DES 4 COTES (BORD) DU CARRE UNITE
C     ======================================================================
      RL1 = 0.
      RL2 = 0.
      DO 1 I = 2, NBS1
C        ** LE COTE 1
         RL1 = RL1 + DIST2P( COSO(1,I-1,1), COSO(1,I,1) )
         STCARR(1,I,1) = RL1
         STCARR(2,I,1) = 0.0
C        ** LE COTE 3
         RL2 = RL2 + DIST2P( COSO(1,I-1,NBS2), COSO(1,I,NBS2) )
         STCARR(1,I,NBS2) = RL2
         STCARR(2,I,NBS2) = 1.0
 1    CONTINUE
C
C     CALCUL DE LA LONGUEUR DES COTES 1 ET 3 : NORMALISATION A LA LONGUEUR 1
      DO 2 I = 2, NBS1-1
C        ** LE COTE 1
         STCARR(1,I,1) = STCARR(1,I,1) / RL1
C        L'ORDONNEE VAUT 0.0
C        ** LE COTE 3
         STCARR(1,I,NBS2) = STCARR(1,I,NBS2) / RL2
C        L'ORDONNEE VAUT 1.0
 2    CONTINUE
C
      RL1 = 0.
      RL2 = 0.
      DO 3 J = 2, NBS2
C        ** LE COTE 2
         RL1 = RL1 + DIST2P( COSO(1,NBS1,J-1), COSO(1,NBS1,J) )
         STCARR(1,NBS1,J) = 1.0
         STCARR(2,NBS1,J) = RL1
C        ** LE COTE 4
         RL2 = RL2 + DIST2P( COSO(1,1,J-1), COSO(1,1,J) )
         STCARR(1,1,J) = 0.0
         STCARR(2,1,J) = RL2
 3    CONTINUE
C
C     CALCUL DE LA LONGUEUR DES COTES 2 ET 4 : NORMALISATION
      DO 4 J = 2, NBS2-1
C        ** LE COTE 2
         STCARR(2,NBS1,J) = STCARR(2,NBS1,J) / RL1
C        L'ABSCISSE VAUT 1.0
C        ** LE COTE 4
         STCARR(2,1,J) = STCARR(2,1,J) / RL2
C        L'ABSCISSE VAUT 0.0
 4    CONTINUE
C
C     LES COORDONNEES DES 4 SOMMETS DU CARRE UNITE
      STCARR(1,1,1) = 0.0
      STCARR(2,1,1) = 0.0
      STCARR(1,NBS1,1) = 1.0
      STCARR(2,NBS1,1) = 0.0
      STCARR(1,1,NBS2) = 0.0
      STCARR(2,1,NBS2) = 1.0
      STCARR(1,NBS1,NBS2) = 1.0
      STCARR(2,NBS1,NBS2) = 1.0
C
C     CALCUL DES 2 COORDONNEES DES SOMMETS INTERNES DU CARRE UNITE
C     ============================================================
      DO 8 J = 2, NBA2
C        LES COORDONNEES DES POINTS  M2 ET M4 DES COTES 2 ET 4
         Y2 = STCARR(2,NBS1,J)
         Y4 = STCARR(2,1,J)
         DO 6 I = 2, NBA1
C           LES COORDONNEES DES POINTS  M1 ET M3 DES COTES 1 ET 3
            X1 = STCARR(1,I,1)
            X3 = STCARR(1,I,NBS2)
C           INTERSECTION DES DROITES M1M3 ET M2M4
            X3X1 = X3 - X1
            Y2Y4 = Y2 - Y4
            IF( ABS(X3X1) .LT. EPSXYZ * (ABS(X1)+ABS(X3)) ) THEN
               XM = X1
               YM = Y4 + Y2Y4 * XM
            ELSE
               IF( ABS(Y2Y4) .LT. EPSXYZ * (ABS(Y2)+ABS(Y4)) ) THEN
                 YM = Y2
                 XM = X1 + X3X1 * YM
               ELSE
                 YM = ( Y4 + Y2Y4 * X1 ) / ( 1 - Y2Y4 * X3X1 )
                 XM = X1 + X3X1 * YM
               ENDIF
            ENDIF
            STCARR(1,I,J) = XM
            STCARR(2,I,J) = YM
 6       CONTINUE
 8    CONTINUE
C
      IF( NBTGS .LE. 0 ) THEN
C
C        PAS DE TANGENTES:
C        CALCUL SEUL DES 3 COORDONNEES DES SOMMETS INTERNES
C        ==================================================
         CALL QUSTIN( NBS1, NBS2, STCARR, COSO )
C
      ELSE
C
C        EXISTENCE DE TANGENTES:
C        CALCUL DES Ds/Du ET Ds/Dv AUX SOMMETS INTERNES DU CARRE UNITE
C        CALCUL DES 3 COORDONNEES DES SOMMETS INTERNES DU QUADRANGLE A MAILLER
C        CALCUL DES DS/DU ET DS/DV AUX 4 SOMMETS DES QUADRANGLES DU QUADRANGLE
C        =====================================================================
         CALL QUSTTG( NBS1,   NBS2,   NBA1,   NBA2,
     S                MNXYT1, MNNTG1, MNXYT2, MNNTG2,
     S                MNXYT3, MNNTG3, MNXYT4, MNNTG4,
     S                COSO,   COTG,   STCARR, TGPURE )
      ENDIF
      END
