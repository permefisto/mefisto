      SUBROUTINE SQUADR( NUCOTE, NBTGS,  NBSOCT, MNSOCT,
     %                   NBARLI, MNXYTG, MNNTGL, MNSOFA, MNFASU,
     %                   MNCARR, XY, MN2DER, XYZDER )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    GENERER LE TABLEAU DES 3 COORDONNEES D'UNE QUADRANGULATION-
C -----    TRIANGULATION D'UN QUADRANGLE COURBE AVEC OU SANS TANGENTES
C
C ENTREES :
C ---------
C NUCOTE : NUCOTE(3) NUMERO DE 1 A 4 DU COTE 3 PARMI LES 4 LIGNES...
C          SI NUCOTE(3)=-4 LE COTE 3 EST LA LIGNE 4 A PARCOURIR
C                          EN SENS INVERSE DE SON RANGEMENT...
C          POUR LE SENS C1:S1S2 C2:S2S3 C3:S3S4 C4:S4S1
C
C                S4----<-------S3
C                |             |
C                |             |
C               \/             /\
C                |             |
C                |             |
C                S1---->-------S2
C
C NBTGS  : =0 SI PAS DE TANGENTES STOCKEES POUR LES 4 LIGNES COTES
C          >0 SINON
C NBSOCT : NOMBRE DE SOMMETS DES 4 COTES DU QUADRANGLE
C MNSOCT : ADRESSE MCN DU TABLEAU XYZSOMMET DES 4 LIGNES
C NBARLI : NOMBRE D'ARETES DE LA LIGNE  (0 SI PROBLEME RENCONTRE)
C MNXYTG : ADRESSE MCN DU TABLEAU DES 3 COMPOSANTES DES TANGENTES
C MNNTGL : ADRESSE MCN DU TABLEAU DES 2 NUMEROS DES TANGENTES DES
C          NBARLI ARETES DE LA LIGNE    (0 SI PROBLEME RENCONTRE)
C          TABLEAU A DETRUIRE EN FIN D'UTILISATION
C MNSOFA : ADRESSE MCN DU TMS 'XYZSOMMET' DE LA SURFACE
C          CF '~/TD/D/A___XYZSOMMET'
C MNFASU : ADRESSE MCN DU TMS 'NSEF' DES NUMEROS DES SOMMETS DES EF
C          CF ~/TD/D/A___NSEF
C
C SORTIES:
C --------
C MNCARR : ADRESSE MCN DU TABLEAU XY DES 2 COORDONNEES DES SOMMETS
C          DU MAILLAGE DU CARRE UNITE
C XY     : 2 COORDONNEES DES SOMMETS DE LA TRIA-QUADRANGULATION
C          SUR LE CARRE UNITE DE REFERENCE
C MN2DER : ADRESSE MCN DU TABLEAU XYZDER DES 3 COMPOSANTES
C          DES 2 DERIVEES ds/du ET ds/dv EN CHAQUE SOMMET DU QUADRANGLE
C XYZDER : 3 COMPOSANTES DES 2 DERIVEES ds/du ET ds/dv PURES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS     OCTOBRE 1996
C2345X7..............................................................012
      PARAMETER       (ITERMX=5)
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/pp.inc"
      COMMON           MCN(MOTMCN)
      REAL            RMCN(1)
      EQUIVALENCE     (MCN(1),RMCN(1))
      REAL             XY(2,*),XYZDER(3,2,*)
C
      INTEGER          NUCOTE(4),NBSOCT(4),NBSTCT(4),MNSOCT(4),MNCUCT(4)
      INTEGER          NBARLI(4),MNXYST(4),
     %                 MNXYTG(4),MNNTGL(4)
      REAL             RLONGC(4)
C
C     FONCTION INTERNE
      XYZ(K,N) =  RMCN(MN+3*N+K)
C
C     LE NUMERO UTILISATEUR DE CHACUN DES 4 COTES
      N1   = ABS( NUCOTE(1) )
      N2   = ABS( NUCOTE(2) )
      N3   = ABS( NUCOTE(3) )
      N4   = ABS( NUCOTE(4) )
C
C     LE NOMBRE DE SOMMETS
C     NBS1   : NOMBRE DE SOMMETS DU COTE 1 ET 3 DU QUADRANGLE
C     NBS2   : NOMBRE DE SOMMETS DU COTE 2 DU QUADRANGLE   NBS2 > NBS4
C     NBS4   : NOMBRE DE SOMMETS DU COTE 4 DU QUADRANGLE
      NBS1 = NBSOCT( N1 )
      NBS2 = NBSOCT( N2 )
C     NBS3 = NBSOCT( N3 )
      NBS4 = NBSOCT( N4 )
      NBS  = NUSOTQ( NBS1, NBS2, NBS4, NBS1, NBS2 )
C
C     LE NUMERO DU SOMMET HAUT GAUCHE DU QUADRANGLE
      NUSTHG = NUSOTQ( NBS1, NBS2, NBS4, 1, NBS4 )
C
C     LE DERNIER SOMMET D'UN QUADRANGLE D'UNE LIGNE
      NDSQ = NBS1 - NBS2 + NBS4
C
C     CREATION DES 2 COORDONNEES DES SOMMETS INTERNES DU MAILLAGE DU CARRE UNITE
C     --------------------------------------------------------------------------
      DO 50 JJ=2,NBS2-1
         DO 40 II=2,NBS1-1
C           RECHERCHE DES COORDONNEES DU SOMMET (II,JJ)
C           DANS LE CARRE DE REFERENCE
            IF( JJ .LT. NBS4 ) THEN
C              LIGNE REGULIERE DE QT
               N = II + (JJ-1) * NBS1
               XY(1,N) = FLOAT(II-1) / (NBS1-1)
               IF( II .LE. NDSQ ) THEN
                  XY(2,N) = FLOAT(JJ-1) / (NBS4-1)
               ELSE
                  XY(2,N) = FLOAT(JJ-1) / (NBS2-NBS1+II-1)
               ENDIF
            ELSE
C              LIGNE IRREGULIERE DE TRIANGLES
               IF( II .LE. NDSQ+JJ-NBS4 ) GOTO 40
C              DANS LE TRIANGLE TOPOLOGIQUE SUPERIEUR DROIT
               N = NUSOTQ( NBS1, NBS2, NBS4, II, JJ )
               XY(1,N) = FLOAT(II-1) / (NBS1-1)
               XY(2,N) = FLOAT(JJ-1) / (NBS2-1-(NBS1-II))
            ENDIF
 40      CONTINUE
 50   CONTINUE
C
C     ATTENTION: LA LOGIQUE DE RANGEMENT DES SOMMETS FRONTALIERS
C     DANS LE SENS DIRECT CHANGE A NOUVEAU EN REVENANT AU
C     RANGEMENT PAR COTES CHANGE A NOUVEAU EN REVENANT AU
C     RANGEMENT PAR COTES INITIAL
C     COTE 1 S1->S2  COTE 2 S2->S3  COTE 3 S4->S3 COTE 4 S1->S4
C
C                S4---->-------S3
C                |             |
C                |             |
C               /\             /\
C                |             |
C                |             |
C                S1---->-------S2
C
      NUCOTE(3) = -NUCOTE(3)
      NUCOTE(4) = -NUCOTE(4)
C
C     CALCUL DES COORDONNEES CURVILIGNES HOMOGENES DES LIGNES
C     DES 4 COTES DU QUADRANGLE COURBE
C     -------------------------------------------------------
      CALL CRCUCT( 4, NUCOTE, MNSOCT, RLONGC, MNCUCT )
C
C     CREATION DES 2 COORDONNEES DES SOMMETS DES 4 COTES DU MAILLAGE DU CARRE UN
C     --------------------------------------------------------------------------
C     CALCUL DES ABSCISSES CURVILIGNES DES SOMMETS SUR CHAQUE COTE DU CARRE
C     LE COTE 1
      MN = MNCUCT(1) - 1
      DO 60 II=1,NBS1
         XY(1,II) = RMCN(MN+II)
         XY(2,II) = 0.0
 60   CONTINUE
C
C     LE COTE 2
      MN = MNCUCT(2) - 1
      NS0C2 = NBS1
      DO 70 JJ=2,NBS2
         NS1C2 = NUSOTQ( NBS1, NBS2, NBS4, NBS1, JJ )
         XY(1,NS1C2) = 1.0
         XY(2,NS1C2) = RMCN(MN+JJ)
         NS0C2 = NS1C2
 70   CONTINUE
C
C     LE COTE 3
      MN = MNCUCT(3) - 1
      NS0C3 = NUSOTQ( NBS1, NBS2, NBS4, 1, NBS4 )
      XY(1,NS0C3) = 0.0
      XY(2,NS0C3) = 1.0
      DO 80 II=2,NBS1
         IF( II .LE. NDSQ ) THEN
            JJ = NBS4
         ELSE
            JJ = NBS2 - NBS1 + II
         ENDIF
         NS1C3 = NUSOTQ( NBS1, NBS2, NBS4, II, JJ )
         XY(1,NS1C3) = RMCN(MN+II)
         XY(2,NS1C3) = 1.0
         NS0C3 = NS1C3
 80   CONTINUE
C
C     LE COTE 4
      MN = MNCUCT(4) - 1
      NS0C4 = 1
      DO 90 JJ=2,NBS4
         NS1C4 = NUSOTQ( NBS1, NBS2, NBS4, 1, JJ )
         XY(1,NS1C4) = 0.0
         XY(2,NS1C4) = RMCN(MN+JJ)
         NS0C4 = NS1C4
 90   CONTINUE
C
C     BARYCENTRAGES DES SOMMETS INTERNES AU CARRE POUR PRENDRE EN COMPTE
C     LA REPARTITION DES SOMMETS SUR LES 4 COTES DU CARRE
C     ==================================================================
      DO 150 K=1,ITERMX
         DO 130 JJ=2,NBS2-1
            DO 120 II=2,NBS1-1
               IF( II .LT. NDSQ ) THEN
C                 BARYCENTRE DU N S O E DES COLONNES DE QUADRANGLES
                  IF( JJ .GE. NBS4 ) GOTO 120
                  N = II + (JJ-1) * NBS1
                  NS0C1 = N-1
                  NS1C1 = N+1
                  NS0C2 = N-NBS1
                  NS1C2 = N+NBS1
               ELSE
                  IF( JJ .GE. NBS2-NBS1+II ) GOTO 120
C                 DANS LES COLONNES DE TRIANGLES
C                 BARYCENTRE = NO SO NE SE
                  N = NUSOTQ( NBS1, NBS2, NBS4, II, JJ )
C                 LES 4 POINTS VOISINS DU POINT N
                  NS0C1 = NUSOTQ( NBS1, NBS2, NBS4, II-1, JJ   )
                  NS1C1 = NUSOTQ( NBS1, NBS2, NBS4, II-1, JJ-1 )
                  NS0C2 = NUSOTQ( NBS1, NBS2, NBS4, II+1, JJ   )
                  NS1C2 = NUSOTQ( NBS1, NBS2, NBS4, II+1, JJ+1 )
               ENDIF
               IF( II .EQ. NDSQ ) NS1C1 = NS0C1
C
C              LE BARYCENTRE DES POINTS
               XY(1,N) = ( XY(1,NS0C1) + XY(1,NS1C1)
     %                   + XY(1,NS0C2) + XY(1,NS1C2) ) * 0.25
               XY(2,N) = ( XY(2,NS0C1) + XY(2,NS1C1)
     %                   + XY(2,NS0C2) + XY(2,NS1C2) ) * 0.25
 120        CONTINUE
 130     CONTINUE
 150  CONTINUE
C
C     CALCUL DU POINT IMAGE SUR LE QUADRANGLE COURBE
C     ==============================================
      IF( NBTGS .LE. 0 ) THEN
C
C        PAS DE TANGENTES AUX ARETES DES EF
C        ----------------------------------
         DO 300 JJ=2,NBS2-1
            DO 290 II=2,NBS1-1
C              RECHERCHE DES COORDONNEES DU SOMMET (II,JJ)
C              DANS LE CARRE DE REFERENCE
               IF( JJ .LE. NBS4 ) THEN
C                 LIGNE REGULIERE DE QT
                  N = II + (JJ-1) * NBS1
               ELSE
C                 LIGNE IRREGULIERE DE TRIANGLES
                  IF( JJ .GE. NBS2 - NBS1 + II ) GOTO 290
C                 DANS LE TRIANGLE TOPOLOGIQUE SUPERIEUR DROIT
                  N = NUSOTQ( NBS1, NBS2, NBS4, II, JJ )
               ENDIF
               X0 = XY(1,N)
               Y0 = XY(2,N)
C
C              RECHERCHE DES 2 SOMMETS ENCADRANTS X0 SUR LE COTE 1
               NS1C1 = 0
               NS0C1 = 1
               DO 210 I=2,NBS1
                  NS1C1 = I
                  IF( X0 .LE. XY(1,NS1C1) ) GOTO 215
                  NS0C1 = NS1C1
 210           CONTINUE
C
C              RECHERCHE DES 2 SOMMETS ENCADRANTS Y0 SUR LE COTE 2
 215           NS0C2 = NBS1
               NS1C2 = 0
               DO 220 J=2,NBS2
                  NS1C2 = NUSOTQ( NBS1, NBS2, NBS4, NBS1, J )
                  IF( Y0 .LE. XY(2,NS1C2) ) GOTO 225
                  NS0C2 = NS1C2
 220           CONTINUE
C
C              RECHERCHE DES 2 SOMMETS ENCADRANTS X0 SUR LE COTE 3
 225           NS0C3 = NUSTHG
               NS1C3 = 0
               DO 230 I=2,NBS1
                  IF( I .LE. NDSQ ) THEN
                     J = NBS4
                  ELSE
                     J = NBS4 + I - NDSQ
                  ENDIF
                  NS1C3 = NUSOTQ( NBS1, NBS2, NBS4, I, J )
                  IF( X0 .LE. XY(1,NS1C3) ) GOTO 235
                  NS0C3 = NS1C3
 230           CONTINUE
C
C              RECHERCHE DES 2 SOMMETS ENCADRANTS Y0 SUR LE COTE 4
 235           NS0C4 = 1
               NS1C4 = 0
               DO 240 J=2,NBS4
                  NS1C4 = NUSOTQ( NBS1, NBS2, NBS4, 1, J )
                  IF( Y0 .LE. XY(2,NS1C4) ) GOTO 245
                  NS0C4 = NS1C4
 240           CONTINUE
C
C              LA TRANSFORMATION DU SOMMET N DU CARRE SUR LE QUADRANGLE
 245           X1 = 1.0 - X0
               Y1 = 1.0 - Y0
               MN = MNSOFA+WYZSOM-4
               DO 270 K=1,3
                  RMCN(MN+3*N+K) =
     %              Y1 * ( (X0-XY(1,NS0C1)) * XYZ(K,NS1C1) +
     %                     (XY(1,NS1C1)-X0) * XYZ(K,NS0C1) )
     %                   / ( XY(1,NS1C1)-XY(1,NS0C1) )
     %           +  X0 * ( (Y0-XY(2,NS0C2)) * XYZ(K,NS1C2) +
     %                     (XY(2,NS1C2)-Y0) * XYZ(K,NS0C2) )
     %                   / ( XY(2,NS1C2)-XY(2,NS0C2) )
     %           +  Y0 * ( (X0-XY(1,NS0C3)) * XYZ(K,NS1C3) +
     %                     (XY(1,NS1C3)-X0) * XYZ(K,NS0C3) )
     %                   / ( XY(1,NS1C3)-XY(1,NS0C3) )
     %           +  X1 * ( (Y0-XY(2,NS0C4)) * XYZ(K,NS1C4) +
     %                     (XY(2,NS1C4)-Y0) * XYZ(K,NS0C4) )
     %                   / ( XY(2,NS1C4)-XY(2,NS0C4) )
     %          - ( X1*Y1*XYZ(K,1)   + X0*Y1*XYZ(K,NBS1) +
     %              X0*Y0*XYZ(K,NBS) + X1*Y0*XYZ(K,NUSTHG) )
 270           CONTINUE
 290        CONTINUE
 300     CONTINUE
C
      ELSE
C
C        PRESENCE DE TANGENTES AUX ARETES DES EF
C        CALCUL DES TABLEAUX DES COORDONNEES DES SOMMETS DES 4 COTES
C        DU NOMBRE D'ARETES ET DES TANGENTES DES 4 COTES DU QUADRANGLE COURBE
C        --------------------------------------------------------------------
         CALL CRTGCT( 4, NUCOTE, MNSOCT, MNXYST, NBARLI, MNNTGL )
         NBSTCT(1) = NBSOCT( N1 )
         NBSTCT(2) = NBSOCT( N2 )
         NBSTCT(3) = NBSOCT( N3 )
         NBSTCT(4) = NBSOCT( N4 )
C
         DO 350 JJ=1,NBS2
            DO 320 II=1,NBS1
C
C              RECHERCHE DES COORDONNEES (X0,Y0) DU SOMMET (II,JJ)
C              DANS LE CARRE DE REFERENCE
               IF( JJ .LE. NBS4 ) THEN
C                 LIGNE REGULIERE DE QT
                  N = II + (JJ-1) * NBS1
               ELSE
C                 LIGNE IRREGULIERE DE TRIANGLES
                  IF( JJ .GT. NBS2 - NBS1 + II ) GOTO 320
C                 DANS LE TRIANGLE TOPOLOGIQUE SUPERIEUR DROIT
                  N = NUSOTQ( NBS1, NBS2, NBS4, II, JJ )
               ENDIF
               X0 = XY(1,N)
               Y0 = XY(2,N)
C
C              CALCUL DES 3 COORDONNEES DU SOMMET SUR LE QUADRANGLE COURBE
C              CALCUL DES DERIVEES  ds/du ET ds/dv PURES
               CALL QU1STG( X0,     Y0,   NBSTCT,
     S              RMCN(MNCUCT(1)), RMCN(MNCUCT(2)),
     S              RMCN(MNCUCT(3)), RMCN(MNCUCT(4)),
     S              MNXYST(N1), MNXYST(N2), MNXYST(N3), MNXYST(N4),
     S              MNXYTG(N1), MNXYTG(N2), MNXYTG(N3), MNXYTG(N4),
     S              MNNTGL(N1), MNNTGL(N2), MNNTGL(N3), MNNTGL(N4),
     S              RMCN(MNSOFA+WYZSOM-3+3*N), XYZDER(1,1,N) )
 320        CONTINUE
 350     CONTINUE
C
C        LA TRANSFORMATION DES TANGENTES CALCULEES EN CHAQUE SOMMET
C        DU MAILLAGE DU CARRE UNITE EN LES 6 OU 8 TANGENTES SUR
C        CHAQUE TRIANGLE OU QUADRANGLE DU MAILLAGE DU CARRE UNITE
C        ----------------------------------------------------------
         CALL QTDETG( MNCARR, MN2DER, MNSOFA, MNFASU )
      ENDIF
      END
