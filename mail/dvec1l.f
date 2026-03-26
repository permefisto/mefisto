      SUBROUTINE DVEC1L( NBSOMT, NBSOMM, MXSOMM, PXYD, AREMAX,
     %                   COSINU, NOANC1, NOANC2,
     %                   MXSTEC, NBSTEC, NOSTEC, N1ST, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TROUVER LA SUITE DES SOMMETS FORMANT L'ENVELOPPE CONVEXE
C -----    DES NBSOMT PREMIERS POINTS DU TABLEAU PXYD
C          LES SOMMETS ALIGNES DE LA PREMIERE ET DERNIERE ARETES SONT
C          CONSERVES DANS CETTE ENVELOPPE CONVEXE
C
C ENTREES:
C --------
C NBSOMT : NOMBRE DE SOMMETS A TRIANGULER
C NBSOMM : NOMBRE DE SOMMETS ACTUELS DANS PXYD
C PXYD   : TABLEAU DES COORDONNEES 2D DES POINTS
C AREMAX : LONGUEUR DE LA PLUS LONGUE ARETE DE LA FRONTIERE
C MXSTEC : NOMBRE MAXIMAL D'ARETES DE L'ENVELOPPE CONVEXE
C
C TABLEAUX AUXILIAIRES :
C ----------------------
C COSINU : REEL2 (1:NBSOMT)
C NOANC1 : ENTIER(1:NBSOMT)
C NOANC2 : ENTIER(1:NBSOMT)
C
C SORTIES:
C --------
C NBSOMM : NOMBRE DE SOMMETS ACTUELS DANS PXYD
C NBSTEC : NOMBRE DE SOMMETS DE L'ENVELOPPE CONVEXE
C NOSTEC : NOSTEC(1,NA) NUMERO DU SOMMET DE L'ARETE NA DE L'ENVELOPPE
C          NOSTEC(2,NA) NUMERO DE L'ARETE SUIVANTE DANS L'ENVELOPPE
C          NA=1 EST LA PREMIERE ARETE DE L'ENVELOPPE PARCOURUE SELON
C          LE SENS DES AIGUILLES D'UNE MONTRE
C          LE CHAINAGE EST CIRCULAIRE
C N1ST   : NUMERO PXYD DU POINT D'ORDONNEE MINIMALE ET
C          A EGALITE D'ABSCISSE MINIMALE DE L'ENVELOPPE CONVEXE
C IERR   : 0 SI PAS D'ERREUR
C          1 SATURATION DES SOMMETS
C          5 SATURATION DES ARETES DU TABLEAU NOSTEC
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC       JUIN 1994
C....................................................................012
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      LOGICAL           TRATRI
      COMMON / DV2DCO / TRATRI
C     TRACE OU NON DES TRIANGLES GENERES DANS LA TRIANGULATION
      INTEGER           NOANC1(NBSOMT),
     %                  NOANC2(NBSOMT),
     %                  NOSTEC(1:2,1:MXSTEC)
      DOUBLE PRECISION  PXYD(3,MXSOMM),
     %                  COSINU(NBSOMT),
     %                  XYB(2)
      DOUBLE PRECISION  COS0,COS1,X,Y,XP1,YP1,PV,
     %                  D0,D1,D2,X1X,Y1Y,X2X,Y2Y,EPSARE
C
C     RECHERCHE DU POINT N1ST DANS PXYD D'ORDONNEE MINIMALE
C     ET A EGALITE D'ABSCISSE MINIMALE. INITIALISATION NOANC1
C     =======================================================
      YP1    = PXYD(2,1)
      N1ST   = 1
      EPSARE = AREMAX * 1D-8
      DO 10 N=2,NBSOMT
         Y = PXYD(2,N)
         IF( Y .GT. YP1 + EPSARE ) GOTO 10
         IF( Y .LT. YP1 - EPSARE ) THEN
            YP1  = Y
            N1ST = N
         ELSE
            IF( PXYD(1,N) .LT. PXYD(1,N1ST) ) THEN
               N1ST = N
            ENDIF
         ENDIF
 10   CONTINUE
C
C     CALCUL DE -COSINUS DES ANGLES (P1X,P1PN) INITIALISATION NOANC1
C     ==============================================================
      XP1 = PXYD(1,N1ST)
      YP1 = PXYD(2,N1ST)
      XYB(1) = XP1
      XYB(2) = YP1
      DO 20 N=1,NBSOMT
         IF( N1ST .NE. N ) THEN
            X = XP1 - PXYD(1,N)
            Y = PXYD(2,N) - YP1
            COSINU(N) = X / SQRT( X * X + Y * Y )
C           LA CONTRIBUTION DU POINT AU BARYCENTRE
            XYB(1) = XYB(1) + PXYD(1,N)
            XYB(2) = XYB(2) + PXYD(2,N)
         ELSE
            COSINU(N) = -1.01D0
         ENDIF
         NOANC1(N) = N
 20   CONTINUE
C
C     LES 2 COORDONNEES DU BARYCENTRE DES POINTS DE L'ENVELOPPE
      XYB(1) = XYB(1) / NBSOMT
      XYB(2) = XYB(2) / NBSOMT
C
C     TRI DES -COSINUS SELON LEUR VALEUR CROISSANTE
C     =============================================
      CALL TRITRD( NBSOMT, COSINU, NOANC1 )
C
C     PERMUTATION DES SOMMETS ALIGNES AVEC LE POINT 1 SUR LA PREMIERE ARETE
C     SELON LA DISTANCE CROISSANTE AU SOMMET 1
C     =====================================================================
      NP0  = NOANC1(2)
      COS0 = COSINU(2)
      DO 27 NBSA1 = 3, NBSOMT
         NP1  = NOANC1( NBSA1 )
         COS1 = COSINU( NBSA1 )
CCC  07/03/97         IF( ABS(COS1-COS0) .GT. 1D-4 ) GOTO 30
         IF( ABS(COS1-COS0) .GT. 1D-7 ) GOTO 30
 27   CONTINUE
 30   NBSA1 = NBSA1 - 1
C
C     ICI LES POINTS 1 A NBSA1 SONT ALIGNES
C     TRI SELON LES DISTANCES CROISSANTES AU POINT N1ST
      COS0 = COSINU(2)
      DO 32 I=1,NBSA1
         NP0 = NOANC1( I )
C        LA DISTANCE DU POINT I AU POINT N1ST
         COSINU(I)  = (PXYD(1,NP0)-XP1)**2 + (PXYD(2,NP0)-YP1)**2
 32   CONTINUE
      CALL TRITRD( NBSA1, COSINU, NOANC1 )
      DO 34 I=1,NBSA1
C        LE COSINUS EST REGENERE
         COSINU(I) = COS0
 34   CONTINUE
C
C     VERIFICATION DU BON ORDRE PAR LE CALCUL DE L'ANGLE (PI PI+1, PI B)
      DO 40 I=2,NBSA1-1
         NP0 = NOANC1(I)
         NP1 = NOANC1(I+1)
C        ATTENTION SINUS( PI PI+1, PI B )
         X = PXYD(1,NP1) - PXYD(1,NP0)
         Y = PXYD(2,NP1) - PXYD(2,NP0)
         X1X = XYB(1) - PXYD(1,NP0)
         Y1Y = XYB(2) - PXYD(2,NP0)
         COS0 = ( X * Y1Y - Y * X1X ) /
     %           SQRT( (X**2+Y**2) * (X1X**2+Y1Y**2) )
         IF( COS0 .LT. 0 ) THEN
C           POINT MAL PLACE : PERMUTATION DES 2 SOMMETS
            NOANC1( I   ) = NP1
            NOANC1( I+1 ) = NP0
         ENDIF
 40   CONTINUE
C
C     PERMUTATION DES DERNIERS SOMMETS ALIGNES AVEC LE POINT 1 SUR LA DERNIERE A
C     SELON LA DISTANCE DECROISSANTE AU SOMMET 1
C     ==========================================================================
      COS0 = COSINU( NBSOMT )
      DO 43 NBSA2=NBSOMT-1, 1, -1
         IF( ABS( COS0 - COSINU(NBSA2) ) .GT. 1D-4 ) GOTO 44
 43   CONTINUE
 44   NBSA2 = NBSA2+1
C     ICI LES SOMMETS NOANC1(NBSOMT) A NOANC1(NBSA2) SONT ALIGNES AVEC 1
C
C     ICI LES POINTS N1ST ET  NBSOMT A NBSA2 SONT ALIGNES
C     TRI SELON LES DISTANCES CROISSANTES AU POINT N1ST
      COS0 = COSINU(NBSOMT)
      DO 51 I=NBSOMT,NBSA2,-1
         NP0 = NOANC1( I )
C        LA DISTANCE DU POINT I AU POINT N1ST
         COSINU(I) = -( (PXYD(1,NP0)-XP1)**2 + (PXYD(2,NP0)-YP1)**2 )
 51   CONTINUE
      CALL TRITRD( NBSOMT-NBSA2+1, COSINU(NBSA2), NOANC1(NBSA2) )
      DO 52 I=NBSA2,NBSOMT
C        LE COSINUS EST REGENERE
         COSINU(I) = COS0
 52   CONTINUE
C
C     VERIFICATION DU BON ORDRE PAR LE CALCUL DE L'ANGLE (PI PI+1, PI B)
      DO 54 I=NBSOMT-1,NBSA2,-1
         NP0 = NOANC1(I)
         NP1 = NOANC1(I+1)
C        ATTENTION SINUS( PI PI+1, PI B )
         X = PXYD(1,NP1) - PXYD(1,NP0)
         Y = PXYD(2,NP1) - PXYD(2,NP0)
         X1X = XYB(1) - PXYD(1,NP0)
         Y1Y = XYB(2) - PXYD(2,NP0)
         COS0 = ( X * Y1Y - Y * X1X ) /
     %           SQRT( (X**2+Y**2) * (X1X**2+Y1Y**2) )
         IF( COS0 .LT. 0 ) THEN
C           POINT MAL PLACE : PERMUTATION DES 2 SOMMETS
            NOANC1( I   ) = NP1
            NOANC1( I+1 ) = NP0
         ENDIF
 54   CONTINUE
C
C     EN CAS D'EGALITE DE L'ANGLE, TRI SELON LA DISTANCE A P1
C     ET SUPPRESSION DES POINTS INTERNES ALIGNES
C     =======================================================
      NOANC2(NBSA1) = NOANC1(NBSA1)
      NBSENV = NBSA1
      NP0    = NOANC1( NBSA1+1 )
      COS0   = COSINU( NBSA1+1 )
      DO 55 N=NBSA1+2,NBSA2
C        LE POINT SELON LE TRI DES ANGLES
         NP1  = NOANC1( N )
         COS1 = COSINU( N )
         IF( COS1 .GT. COS0 ) THEN
C           LE POINT N N'EST PAS ALIGNE AVEC N1ST ET NP0
            NBSENV = NBSENV + 1
            NOANC2(NBSENV) = NP0
            COS0   = COS1
            NP0    = NP1
         ELSE
C           LE POINT N EST ALIGNE AVEC N1ST ET NP0
            D0  = (PXYD(1,NP0)-XP1)**2 + (PXYD(2,NP0)-YP1)**2
            D1  = (PXYD(1,NP1)-XP1)**2 + (PXYD(2,NP1)-YP1)**2
            IF( D0 .LT. D1 ) THEN
C              LE PLUS ELOIGNE DE N1ST EST NP1
               NP0 = NP1
            ENDIF
         ENDIF
 55   CONTINUE
C     LE DERNIER EST AJOUTE
      NBSENV = NBSENV + 1
      NOANC2(NBSENV) = NP0
C
C     GESTION DU CONTOUR FERME DE L'ENVELOPPE CONVEXE ACTUELLE
C     ========================================================
C     CHAINAGE CIRCULAIRE (SENS DES AIGUILLES D'UNE MONTRE)
C     L'ARETE 1 DE L'ENVELOPPE (SOMMET, ARETE SUIVANTE)
      NOSTEC(1,1) = NOANC1(1)
      NOSTEC(2,1) = 2
C     L'ARETE 2 DE L'ENVELOPPE (SOMMET, ARETE SUIVANTE)
      NOSTEC(1,2) = NOANC2(NBSA1+1)
      NOSTEC(2,2) = 3
C     L'ARETE 3 DE L'ENVELOPPE (SOMMET, ARETE SUIVANTE)
      NOSTEC(1,3) = NOANC2(NBSA1)
      NOSTEC(2,3) = 1
C
C     CHAINAGE DES ARETES VIDES DE L'ENVELOPPE
C     ========================================
C     N1STEC NUMERO DE LA PREMIERE ARETE VIDE DESNOSTEC
      N1STEC = 4
      DO 57 I=N1STEC,MXSTEC
C        CHAINAGE SUR LE SUIVANT
         NOSTEC(2,I) = I+1
 57   CONTINUE
      NOSTEC(2,MXSTEC) = 0
      NOS = NBSA1+2
C
C     =======================================================
C     LA BOUCLE DE L'ENVELOPPE CONVEXE DES SOMMETS 4 A NBSENV
C     =======================================================
 60   IF( NOS .LE. NBSENV ) THEN
C
C        LE VRAI NUMERO DU SOMMET NOS A TRAITER MAINTENANT
         NP = NOANC2(NOS)
C
C        LA PREMIERE ARETE DE L'ENVELOPPE
         J1  = 1
C
C        L'ARETE J1 DE L'ENVELOPPE EST ELLE VISIBLE DE NP=NOANC1(NOS)?
C        =============================================================
         NBA = 0
         X   = PXYD(1,NP)
         Y   = PXYD(2,NP)
C
 70      J2  = NOSTEC(2,J1)
C        LES 2 SOMMETS DE L'ARETE J1
         NS1 = NOSTEC(1,J1)
         NS2 = NOSTEC(1,J2)
         X1  = REAL( PXYD(1,NS1) )
         Y1  = REAL( PXYD(2,NS1) )
         X2  = REAL( PXYD(1,NS2) )
         Y2  = REAL( PXYD(2,NS2) )
         X1X = X1 - X
         Y1Y = Y1 - Y
         D1  = X1X * X1X + Y1Y * Y1Y
         X2X = X2 - X
         Y2Y = Y2 - Y
         D2  = X2X * X2X + Y2Y * Y2Y
C        NP NS1  PRODUIT VECTORIEL NP NS2
         PV  = ( X1X * Y2Y - Y1Y * X2X ) / SQRT( D1 * D2 )
         IF( PV .GT. -1D-5 ) THEN
C
C           NP EST IL ALIGNE AVEC NS1 NS2 ?
            IF( PV .LE. 1D-4 ) THEN
C              OUI : NP EST ALIGNE AVEC NS1 NS2
C              CALCUL DU SIGNE DU COSINUS(NS1 NP, NS1 NS2)
               COS1 = X1X * (X1-X2) + Y1Y * (Y1-Y2)
               IF( COS1 .LT. 0 ) GOTO 72
C              ICI L'ANGLE (NS1 NP , NS1 NS2) EST POSITIF TRES PETIT
            ENDIF
C
C           AJOUT DU POINT NP A L'ENVELOPPE
            IF( NBA .EQ. 0 ) THEN
C              SUPPRESSION D'UNE ARETE ET AJOUT DE 2 AUTRES
               IF( N1STEC .LE. 0 ) THEN
C                 SATURATION DES ARETES DE L'ENVELOPPE
                  NBLGRC(NRERR) = 1
                  KERR(1) = 'SATURATION DE L''ENVELOPPE DES ARETES'
                  CALL LEREUR
                  IERR = 5
                  GOTO 9999
               ENDIF
C
C              LA NOUVELLE ARETE A AJOUTER EST CREEE
               N1     = N1STEC
               N1STEC = NOSTEC(2,N1STEC)
               NOSTEC(1,N1) = NP
               NOSTEC(2,N1) = J2
C              MISE A JOUR DE L'ARETE J1 DE L'ENVELOPPE
               NOSTEC(2,J1) = N1
C              POUR LA SUITE
               J0  = N1
               J1  = J2
               NBA = 1
C
            ELSE
C
C              MODIFICATION DE L'ARETE
               NOSTEC(2,J0) = J2
C              J1 EST VIDE ET RENDUE A NOSTEC
               NOSTEC(2,J1) = N1STEC
               N1STEC = J1
C              POUR LA SUITE
               J1 = J2
            ENDIF
            GOTO 70
C
C           FIN DU REPETER JUSQU'A ARETE NON SEPARANTE
         ENDIF
C
C        TRACE DE L'ENVELOPPE CONVEXE
 72      IF( TRATRI ) THEN
            CALL EFFACE
            CALL DVTRCF( 1, 2, NOSTEC, PXYD )
         ENDIF
C        FIN TANT QUE IL EXISTE UN POINT A AJOUTER
         NOS = NOS + 1
         GOTO 60
      ENDIF
C
C     CALCUL DE LA LONGUEUR MOYENNE DES ARETES DE L'ENVELOPPE
C     =======================================================
      D0     = 0
      NBSTEC = 0
      NP1    = NOSTEC(1,1)
      N1     = NOSTEC(2,1)
C
 75   NBSTEC = NBSTEC + 1
      NP0    = NP1
      NP1    = NOSTEC(1,N1)
      D0     = D0 + SQRT( (PXYD(1,NP1)-PXYD(1,NP0))**2 +
     %                    (PXYD(2,NP1)-PXYD(2,NP0))**2  )
      N1 = NOSTEC(2,N1)
      IF( N1 .NE. 1 ) GOTO 75
C
      D0 = MAX( D0 / NBSTEC , DBLE(AREMAX) )
C
C     PARCOURS ET DECOUPAGE DES ARETES TROP LONGUES
C     =============================================
      N0   = 1
      N1   = NOSTEC(2,1)
      IF( NBSA2 .LT. NBSOMT ) THEN
C        LA PREMIERE ARETE EST EVITEE
         N0 = N1
         N1 = NOSTEC(2,N0)
      ENDIF
C
 80   NP0  = NOSTEC(1,N0)
      NP1  = NOSTEC(1,N1)
      D1   = SQRT( (PXYD(1,NP1)-PXYD(1,NP0))**2 +
     %             (PXYD(2,NP1)-PXYD(2,NP0))**2  )
C     LE NOMBRE D'ARETES MOYENNES DANS CETTE ARETE
      N = INT( D1 / D0 )
      IF( N .GE. 2 ) THEN
C
C        N-1 POINTS SONT INTERCALES SUR CETTE ARETE
         X = (PXYD(1,NP1)-PXYD(1,NP0)) / N
         Y = (PXYD(2,NP1)-PXYD(2,NP0)) / N
         DO 82 I=1,N-1
C
C           CONSTRUCTION DU POINT
            IF( NBSOMM .GE. MXSOMM ) THEN
               NBLGRC(NRERR) = 1
               KERR(1) = 'SATURATION DES SOMMETS'
               CALL LEREUR
               IERR = 1
               GOTO 9999
            ENDIF
            NBSOMM = NBSOMM + 1
            PXYD(1,NBSOMM) = PXYD(1,NP0) + I * X
            PXYD(2,NBSOMM) = PXYD(2,NP0) + I * Y
            PXYD(3,NBSOMM) = AREMAX*3
C
C           LA NOUVELLE ARETE A AJOUTER DANS LE CHAINAGE DE L'ENVELOPPE
            N2     = N1STEC
            N1STEC = NOSTEC(2,N1STEC)
C           LE SOMMET
            NOSTEC(1,N2) = NBSOMM
C           LE CHAINAGE SUR L'ARETE SUIVANTE
            NOSTEC(2,N2) = N1
C           LA PRECEDENTE ARETE POINTE SUR LA NOUVELLE
            NOSTEC(2,N0) = N2
C           LA SUIVANTE DE LA PRECEDENTE EST LA NOUVELLE
            N0 = N2
 82      CONTINUE
      ENDIF
C
C     PASSAGE A L'ARETE SUIVANTE
      IF( N1 .NE. 1 ) THEN
         N0 = N1
         N1 = NOSTEC(2,N1)
C        S'IL EXISTE DES SOMMETS SUR LA DERNIERE ARETE PAS D'AJOUT
         IF( N1 .EQ. 1 .AND. NBSA1 .GT. 2 ) GOTO 86
         GOTO 80
      ENDIF
C
C     AJOUT DES SOMMETS DE LA DERNIERE ARETE ALIGNES AVEC LE POINT 1
C     ==============================================================
C     NOANC1(NBSA2+1 A NBSOMT ) DOIVENT ETRE AJOUTES A L'ENVELOPPE CONVEXE
 86   I = 1
      DO 85 J = NBSOMT, NBSA2+1, -1
         IF( N1STEC .LE. 0 ) THEN
C           SATURATION DES ARETES DE L'ENVELOPPE
            NBLGRC(NRERR) = 1
            KERR(1) = 'SATURATION DE L''ENVELOPPE DES ARETES'
            CALL LEREUR
            IERR = 5
            GOTO 9999
         ENDIF
C        LA NOUVELLE ARETE A AJOUTER
         N1     = N1STEC
         N1STEC = NOSTEC(2,N1STEC)
         NOSTEC(1,N1) = NOANC1( J )
         NOSTEC(2,N1) = NOSTEC(2,I)
         NOSTEC(2,I) = N1
         I           = N1
 85   CONTINUE
C
C     AJOUT DES SOMMETS DE L'ARETE 1 DE L'ENVELOPPE CONVEXE
C     =====================================================
C     LE NOMBRE DE SOMMETS DE L'ENVELOPPE CONVEXE
      NBSTEC = 1
      NOS    = 1
 90   IF( NOSTEC(2,NOS) .NE. 1 ) THEN
         NBSTEC = NBSTEC + 1
         NOS    = NOSTEC(2,NOS)
         GOTO 90
      ENDIF
C
      DO 95 J = NBSA1-1, 2, -1
         IF( N1STEC .LE. 0 ) THEN
C           SATURATION DES ARETES DE L'ENVELOPPE
            NBLGRC(NRERR) = 1
            KERR(1) = 'SATURATION DE L''ENVELOPPE DES ARETES'
            CALL LEREUR
            IERR = 5
            GOTO 9999
         ENDIF
C        LA NOUVELLE ARETE A AJOUTER
         N1     = N1STEC
         N1STEC = NOSTEC(2,N1STEC)
         NOSTEC(1,N1)  = NOANC1( J )
         NOSTEC(2,N1)  = NOSTEC(2,NOS)
         NOSTEC(2,NOS) = N1
         NOS           = N1
 95   CONTINUE
      NBSTEC = NBSTEC + NBSA1 - 2
C
C     VERIFICATION DU BON ORDRE PAR LE CALCUL DE L'ANGLE (PI PI+1, PI B)
C     QUI NE DOIT PAS ETRE SUPERIEUR A 180 DEGRES
C     ATTENTION L'ENVELOPPE CONVEXE EST PARCOURUE DANS LE SENS INDIRECT
C     ENVELOPPE PARCOURUE SELON LE SENS DES AIGUILLES D'UNE MONTRE
C     ==================================================================
      N0 = 1
      DO 200 I=1,NBSTEC
C        L'ARETE SUIVANTE
         N1  = NOSTEC(2,N0)
C        LE NUMERO DU SOMMET DE L'ARETE N1
         NP0 = NOSTEC(1,N1)
C        L'ARETE SUIVANTE
         N2  = NOSTEC(2,N1)
C        LE NUMERO DU SOMMET DE L'ARETE N2
         NP1 = NOSTEC(1,N2)
C        ATTENTION SINUS( PI PI+1, PI B )
         X = PXYD(1,NP1) - PXYD(1,NP0)
         Y = PXYD(2,NP1) - PXYD(2,NP0)
         X1X = XYB(1) - PXYD(1,NP0)
         Y1Y = XYB(2) - PXYD(2,NP0)
         COS0 = ( X * Y1Y - Y * X1X ) /
     %           SQRT( (X**2+Y**2) * (X1X**2+Y1Y**2) )
         IF( COS0 .GT. 0 ) THEN
C           POINT MAL PLACE : PERMUTATION DES 2 SOMMETS
            WRITE(IMPRIM,*) 'DVEC1L: SOMMETS ',NP0,' ',NP1,
     %                      ' A PERMUTER. REVOIR LES SEUILS'
            NOSTEC(2,N1) = NOSTEC(2,N2)
            NOSTEC(2,N2) = N1
            NOSTEC(2,N0) = N2
            N1 = N2
         ENDIF
C
C        PASSAGE A L'ARETE SUIVANTE
         N0 = N1
 200  CONTINUE
C
C     TRACE DE L'ENVELOPPE CONVEXE FINALE
      IF( TRATRI ) THEN
         CALL EFFACE
         CALL DVTRCF( 1, 2, NOSTEC, PXYD )
      ENDIF
C
 9999 RETURN
      END
