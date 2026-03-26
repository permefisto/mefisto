      SUBROUTINE TRDVAN( MXSOMM, NBSOMT, MXTRIA, PXYD , AREMAX,
     %                   N1TRVI, NOTRIA, CETRIA,
     %                   COSINU, NOANCP, MXARCF, CF,
     %                   N1ST,   N2ST,   IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    GENERER LES TRIANGLES DELAUNAY DES SOMMETS
C -----    1 A NBSOMT DU TABLEAU PXYD A L'AIDE D'UN TRI DES ANGLES
C
C ENTREES:
C --------
C MXSOMM : NOMBRE MAXIMAL DE SOMMETS PERMIS POUR LA TRIANGULATION
C NBSOMT : NOMBRE DE SOMMETS A TRIANGULER
C MXTRIA : NOMBRE MAXIMAL DE TRIANGLES DECLARABLES
C PXYD   : TABLEAU DES COORDONNEES 2D DES POINTS
C AREMAX : LONGUEUR DE LA PLUS LONGUE ARETE DE LA FRONTIERE
C
C ENTREES ET SORTIES :
C --------------------
C N1TRVI : NUMERO DU 1 PREMIER TRIANGLE VIDE DANS LE TABLEAU NOTRIA
C          LE CHAINAGE DES TRIANGLES VIDES SE FAIT SUR NOTRIA(4,.)
C NOTRIA : LISTE DES TRIANGLES
C                 ------- ------- ------- -------- -------- --------
C  PAR TRIANGLE : SOMMET1 SOMMET2 SOMMET3 TR_VOIS1 TR_VOIS2 TR_VOIS3
C                 ------- ------- ------- -------- -------- --------
C                 SOMMET    EST LE NUMERO DU SOMMET
C                 TR_VOIS i EST LE NUMERO DANS NOTRIA DU TRIANGLE
C                               ADJACENT PAR L'ARETE i
C
C CETRIA : COORDONNEES DU CENTRE DU CERCLE CIRCONSCRIT ET
C          CARRE DU RAYON
C                         ------- ------- --------
C          PAR TRIANGLE : XCENTRE YCENTRE RAYON**2
C                         ------- ------- --------
C          TABLEAU REEL2(3,MXTRIA)
C
C MXARCF : NOMBRE MAXIMAL D'ARETES DU CF
C
C SORTIES:
C --------
C N1ST   : NUMERO PXYD DU PREMIER SOMMET MAILLE
C N2ST   : NUMERO PXYD DU SECOND  SOMMET MAILLE
C IERR   : 0 SI PAS D'ERREUR
C          1 SATURATION DES SOMMETS
C          3 SATURATION DES TRIANGLES
C          4 SATURATION DES ARETES DE LA PILE
C          5 SATURATION DES ARETES DU CF
C          9 TOUS LES POINTS SONT ALIGNES
C
C TABLEAUX AUXILIAIRES :
C ----------------------
C COSINU : REEL2(1:NBSOMT)
C NOANCP : ENTIER(1:NBSOMT)
C CF     : ENTIER(1:3,1:MXARCF)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC    JANVIER 1994
C....................................................................012
      DOUBLE PRECISION  EPSMIN,      EPSSIN
      PARAMETER        (EPSMIN=1D-8, EPSSIN=1D-9)
      include"./incl/gsmenu.inc"
      include"./incl/trvari.inc"
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      INTEGER           NOTRIA(6,MXTRIA),
     %                  NOANCP(NBSOMT),
     %                  CF(1:3,1:MXARCF)
      DOUBLE PRECISION  PXYD(3,MXSOMM),
     %                  CETRIA(3,MXTRIA),
     %                  COSINU(NBSOMT)
      DOUBLE PRECISION  X,Y,XP1,YP1,PV,D0,D1,D2,X1X,Y1Y,X2X,Y2Y,
     %                  CO,EPSARE
C
C     RECHERCHE DU POINT N1ST DANS PXYD D'ORDONNEE MINIMALE
C     ET A EGALITE D'ABSCISSE MINIMALE. INITIALISATION NOANCP
C     =======================================================
      YP1    = PXYD(2,1)
      N1ST   = 1
      EPSARE = AREMAX * EPSMIN
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
C     CALCUL DU COSINUS DES ANGLES (P1X,P1PN) INITIALISATION NOANCP
C     =============================================================
      XP1 = PXYD(1,N1ST)
      YP1 = PXYD(2,N1ST)
      DO 20 N=1,NBSOMT
         IF( N1ST .NE. N ) THEN
            X = PXYD(1,N) - XP1
            Y = PXYD(2,N) - YP1
            COSINU(N) = X / SQRT( X * X + Y * Y )
            IF( ABS(COSINU(N)) .LT. EPSSIN ) COSINU(N) = 0D0
         ELSE
            COSINU(N) = 1.01D0
         ENDIF
         NOANCP(N) = N
 20   CONTINUE
C
C     TRI DES COSINUS SELON LEUR VALEUR DECROISSANTE
C     ==============================================
      CALL TRITRD( NBSOMT, COSINU, NOANCP )
C     MISE SOUS LA FORME DECROISSANTE
      N2 = NBSOMT / 2
      DO 22 N=1,N2
         N1 = NBSOMT + 1 - N
         X          = COSINU(N )
         COSINU(N ) = COSINU(N1)
         COSINU(N1) = X
         NAR        = NOANCP(N )
         NOANCP(N ) = NOANCP(N1)
         NOANCP(N1) = NAR
 22   CONTINUE
C
C     EN CAS D'EGALITE DE L'ANGLE, TRI SELON LA DISTANCE A P1
C     =======================================================
      J0 = 2
      X  = COSINU(2)
      DO 27 N=3,NBSOMT
         Y  = COSINU( N )
         IF( ABS(X-Y) .LE. EPSSIN ) THEN
C
C           LE POINT N EST ALIGNE AVEC N1ST ET J0
C           TRI SELON LA DISTANCE AU POINT P1
            DO 23 K=N,J0+1,-1
               NS1 = NOANCP(K)
               NS0 = NOANCP(K-1)
               D1  = (PXYD(1,NS1)-XP1)**2 + (PXYD(2,NS1)-YP1)**2
               D0  = (PXYD(1,NS0)-XP1)**2 + (PXYD(2,NS0)-YP1)**2
               IF( D0 .GT. D1 ) THEN
C                 PERMUTATION DES POINTS K-1 ET K
C                 LE PLUS ELOIGNE EST NUMEROTE EN DERNIER
                  NOANCP(K  ) = NS0
                  NOANCP(K-1) = NS1
               ELSE
                  GOTO 26
               ENDIF
 23         CONTINUE
         ELSE
C
C           LE POINT N N'EST PAS ALIGNE AVEC N1ST ET J0
C           ALIGNEMENT EN Y DES POINTS INTERMEDIAIRES
            CALL ALIGPD( J0, N-1, NOANCP, PXYD, EPSARE )
C
            J0 = N
         ENDIF
 26      X = Y
 27   CONTINUE
C
C     ALIGNEMENT EN Y DES POINTS INTERMEDIAIRES
      CALL ALIGPD( J0, NBSOMT, NOANCP, PXYD, EPSARE )
C
C     RECHERCHE DU PREMIER POINT N NON ALIGNE AVEC 1 2
C     ================================================
      N2ST = NOANCP(2)
      X2   = REAL( PXYD(1,N2ST) )
      Y2   = REAL( PXYD(2,N2ST) )
      DO 28 N=3,NBSOMT
C        NP N1ST  PRODUIT VECTORIEL NP N2ST
         NP = NOANCP(N)
         X  = REAL( PXYD(1,NP) )
         Y  = REAL( PXYD(2,NP) )
         X1X = XP1 - X
         Y1Y = YP1 - Y
         D1  = X1X * X1X + Y1Y * Y1Y
         X2X = X2 - X
         Y2Y = Y2 - Y
         D2  = X2X * X2X + Y2Y * Y2Y
         PV  = ( X1X * Y2Y - Y1Y * X2X ) / SQRT(D1*D2)
         IF( ABS(PV) .GT. EPSSIN ) GOTO 30
 28   CONTINUE
C
C     ICI TOUS LES POINTS SONT ALIGNES
      NBLGRC(NRERR) = 1
      KERR(1) = 'SURFACE REDUITE A UN SEGMENT DE DROITE'
      CALL LEREUR
      IERR = 9
      GOTO 9999
C
C     LES N-2 PREMIERS TRIANGLES
C     ==========================
 30   DO 40 I=1,N-2
C        TRIANGLE I DE SOMMET I,I+1,N APRES TRI
         NS1 = NOANCP(I)
         NS2 = NOANCP(I+1)
         NOTRIA(1,I) = NS1
         NOTRIA(2,I) = NS2
         NOTRIA(3,I) = NP
         NOTRIA(4,I) = 0
         NOTRIA(5,I) = I+1
         NOTRIA(6,I) = I-1
C        CALCUL DU CENTRE DU CERCLE CIRCONSCRIT ET DU CARRE DU RAYON
         IER = 1
         CALL CENCED( PXYD(1,NS1), PXYD(1,NS2), PXYD(1,NP),
     %                CETRIA(1,I), IER )
         IF( IER .NE. 0 ) THEN
            WRITE(IMPRIM,*)'TRDVAN:TRIANGLE ',I,
     %      ' DEGENERE SOMMET1=',NS1,' ',NS2,' ',NP
         ENDIF
CCCC        TRACE DU TRIANGLE
CCC         NCOUL = MOD(I,NDCOUL-N1COUL)+N1COUL
CCC         CALL DVTRTR( PXYD, NOTRIA, I, NCOUL, NCNOIR )
 40   CONTINUE
C
C     CORRECTIONS
C     LE TRIANGLE OPPOSE PAR L'ARETE 2 DU DERNIER TRIANGLE
      NOTRIA(5,N-2) = 0
C
C     CHAINAGE DES TRIANGLES VIDES
C     ============================
      N1TRVI = N-1
C     CHAINAGE DES TRIANGLES VIDES SUIVANTS
      DO 45 I = N1TRVI , MXTRIA
C        LE NUMERO DU SOMMET 1 = 0 => TRIANGLE VIDE
         NOTRIA(1,I) = 0
C        LE TRIANGLE VIDE SUIVANT
         NOTRIA(4,I) = I + 1
 45   CONTINUE
C     LA FIN DU CHAINAGE DES SUIVANTS VIDES
      NOTRIA(4,MXTRIA) = 0
C
C     GESTION DU CONTOUR FERME DE L'ENVELOPPE CONVEXE ACTUELLE
C     ========================================================
C     CHAINAGE CIRCULAIRE (SENS DES AIGUILLES D'UNE MONTRE)
      NS1 = NOANCP(1)
      DO 50 I=1,N
C        L'ARETE I DU CF (SOMMET,TRIANGLE OPPOSE, SUIVANTE)
         CF(1,I) = NS1
         CF(3,I) = I+1
         NS1 = NOANCP( N + 1 - I )
 50   CONTINUE
      CF(3,N) = 1
C
C     LE TRIANGLE OPPOSE A L'ARETE
      CF(2,1) = 1
      CF(2,2) = N-2
      DO 55 I=1,N-2
         CF(2,I+2) = N-1-I
 55   CONTINUE
      N0S = N
C
C     CHAINAGE DES ARETES VIDES DU CF
C     ===============================
C     N1ARCF NUMERO DE LA PREMIERE ARETE VIDE DES CF
      N1ARCF = N+1
      DO 57 I=N1ARCF,MXARCF
C        CHAINAGE SUR LE SUIVANT
         CF(3,I) = I+1
 57   CONTINUE
      CF(3,MXARCF) = 0
C
C     ===================================================
C     LA BOUCLE DE TRIANGULATION DES SOMMETS N+1 A NBSOMT
C     ===================================================
      NOS = N+1
 60   IF( NOS .LE. NBSOMT ) THEN
C
C        LE VRAI NUMERO DU SOMMET NOS
         NP = NOANCP(NOS)
C        SON COSINUS
         CO = COSINU(NOS)
C
C        LE TRIANGLE QUI PRECEDE
         NT0 = 0
         NBT = 0
C        LA PREMIERE ARETE DU CF
         J1 = 1
C
C        RECHERCHE DU POINT SUIVANT NON SUR LA DROITE 1-NP
         DO 62 N=NOS+1,NBSOMT
            IF( ABS(COSINU(N)-CO) .GT. EPSSIN ) GOTO 63
 62      CONTINUE
         N = NBSOMT + 1
C
 63      N = N - 1
C        N EST LE DERNIER POINT SUR LA DROITE N1ST-NP
         IF( N .GT. NOS ) THEN
C
C           IL EXISTE PLUSIEURS POINTS ALIGNES AVEC N1ST
C           ============================================
C           GENERATION DU TRIANGLE CF(1)-NP
C           L'ARETE SUIVANTE
            J2 = CF(3,J1)
C           LES 2 SOMMETS DE L'ARETE J1 DU CF
            NS1  = CF(1,J1)
C           LE TRIANGLE ADJACENT
            NTOP = CF(2,J1)
            NS2  = CF(1,J2)
C
            IF( NOS-1 .EQ. N0S ) THEN
C
C              LE TRIANGLE NS1-NS2-NP EST IL DEGENERE ?
C              ----------------------------------------
               X   = REAL( PXYD(1,NP) )
               Y   = REAL( PXYD(2,NP) )
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
C              NP NS1  PRODUIT VECTORIEL NP NS2
               PV  = ( X1X * Y2Y - Y1Y * X2X ) / SQRT(D1*D2)
               IF( PV .LT. EPSSIN ) THEN
C                 TRIANGLE DEGENERE
                  J0  = J1
                  J1  = J2
                  J2  = CF(3,J2)
                  NS1 = NS2
                  NS2 = CF(1,J2)
                  NT0 = CF(2,J1)
                  NOS = NOS - 1
                  GOTO 64
               ENDIF
            ENDIF
C
C           UN TRIANGLE NS1-NS2-NP A FORMER
            CALL CRETRD( NS1, NS2, NP, NTOP, 0, NT0,
     %                   PXYD, N1TRVI, NOTRIA, CETRIA,
     %                   NT )
            IF( NT .LE. 0 ) THEN
C              SATURATION DES TRIANGLES
               IERR = 3
               GOTO 9999
            ENDIF
            IF( CETRIA(3,NT) .GE. 1D28 ) THEN
               WRITE(IMPRIM,*)'TRDVAN:TRIANGLE ',NT,
     %       ' DEGENERE SOMMET2=',(NOTRIA(K,NT),K=1,3)
            ENDIF
            NBT = NBT + 1
CCCC           TRACE DU TRIANGLE
CCC            NCOUL = MOD(I,NDCOUL-N1COUL)+N1COUL
CCC            CALL DVTRTR( PXYD, NOTRIA, NT, NCOUL, NCNOIR )
C
C           MISE A JOUR DU CF
            CF(2,J1) = NT
            J0  = J1
            NS1 = NP
            NT0 = NT
C
 64         DO 65 I=NOS+1,N
C
C              GENERATION DES TRIANGLES ARETE DE LA DROITE N1ST-NP  NS2
C              UN TRIANGLE NS1-NS2-NP A FORMER
               NP = NOANCP(I)
               CALL  CRETRD( NS1, NS2, NP, NT0, 0, 0,
     %                       PXYD, N1TRVI, NOTRIA, CETRIA,
     %                       NT )
               IF( NT .LE. 0 ) THEN
C                 SATURATION DES TRIANGLES
                  IERR = 3
                  GOTO 9999
               ENDIF
               IF( CETRIA(3,NT) .GE. 1D28 ) THEN
                  WRITE(IMPRIM,*)'TRDVAN:TRIANGLE ',NT,
     %          ' DEGENERE SOMMET3=',(NOTRIA(K,NT),K=1,3)
               ENDIF
CCCC              TRACE DU TRIANGLE
CCC               NCOUL = MOD(I,NDCOUL-N1COUL)+N1COUL
CCC               CALL DVTRTR( PXYD, NOTRIA, NT, NCOUL, NCNOIR )
C
C              MISE A JOUR DU CF
               IF( N1ARCF .LE. 0 ) THEN
C                 SATURATION DES ARETES DU CF
                  NBLGRC(NRERR) = 1
                  KERR(1) = 'SATURATION DU CF DES ARETES'
                  CALL LEREUR
                  IERR = 5
                  GOTO 9999
               ENDIF
C              LA NOUVELLE ARETE DU CF A AJOUTER
               J1     = N1ARCF
               N1ARCF = CF(3,N1ARCF)
               CF(1,J1) = NS1
               CF(2,J1) = NT
C              CHAINAGE DU PRECEDENT
               CF(3,J0 ) = J1
C
C              POUR LA SUITE
               J0  = J1
               NS1 = NP
               NT0 = NT
               NBT = NBT + 1
 65         CONTINUE
C
C           L'ARETE NP-NS2 DU CF
            IF( N1ARCF .LE. 0 ) THEN
C              SATURATION DES ARETES DU CF
               NBLGRC(NRERR) = 1
               KERR(1) = 'SATURATION DU CF DES ARETES'
               CALL LEREUR
               IERR = 5
               GOTO 9999
            ENDIF
C           LA NOUVELLE ARETE A AJOUTER
            J1     = N1ARCF
            N1ARCF = CF(3,N1ARCF)
            CF(1,J1) = NP
            CF(2,J1) = NT
C           CHAINAGE DU PRECEDENT
            CF(3,J0 ) = J1
            CF(3,J1 ) = J2
C
C           LE POINT EN COURS DE TRAITEMENT
            NOS = N
            J0  = J1
            J1  = J2
         ENDIF
C
C        L'ARETE J1 DU CF EST ELLE VISIBLE DE NP=NOANCP(NOS)?
C        ====================================================
         X  = PXYD(1,NP)
         Y  = PXYD(2,NP)
C
 70      J2 = CF(3,J1)
C        LES 2 SOMMETS DE L'ARETE J1
         NS1 = CF(1,J1)
         NS2 = CF(1,J2)
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
         PV  = ( X1X * Y2Y - Y1Y * X2X ) / SQRT(D1*D2)
         IF( PV .GT. EPSSIN ) THEN
C
C           INITIALISATION DU NOUVEAU TRIANGLE NT=NS1-NS2-NP
            NBT  = NBT + 1
            NTOP = CF(2,J1)
            CALL CRETRD( NS1, NS2, NP, NTOP, 0, NT0,
     %                   PXYD, N1TRVI, NOTRIA, CETRIA,
     %                   NT )
            IF( NT .LE. 0 ) THEN
C              SATURATION DES TRIANGLES
               IERR = 3
               GOTO 9999
            ENDIF
            IF( CETRIA(3,NT) .GE. 1D28 ) THEN
               WRITE(IMPRIM,*)'TRDVAN:TRIANGLE ',NT,
     %       ' DEGENERE SOMMET4=',(NOTRIA(K,NT),K=1,3)
            ENDIF
CCCC           TRACE DU TRIANGLE NT
CCC            NCOUL = MOD(NT,NDCOUL-N1COUL)+N1COUL
CCC            CALL DVTRTR( PXYD, NOTRIA, NT, NCOUL, NCNOIR )
C
            NT0 = NT
C
C           MISE A JOUR DU CF
            IF( NBT .EQ. 1 ) THEN
C              SUPPRESSION D'UNE ARETE ET AJOUT DE 2 AUTRES
               IF( N1ARCF .LE. 0 ) THEN
C                 SATURATION DES ARETES DU CF
                  NBLGRC(NRERR) = 1
                  KERR(1) = 'SATURATION DU CF DES ARETES'
                  CALL LEREUR
                  IERR = 5
                  GOTO 9999
               ENDIF
C              LA NOUVELLE ARETE A AJOUTER
               N1     = N1ARCF
               N1ARCF = CF(3,N1ARCF)
               CF(1,N1) = NP
               CF(2,N1) = NT
               CF(3,N1) = J2
C              MISE A JOUR DE L'ARETE J1 DU CF
               CF(2,J1) = NT
               CF(3,J1) = N1
C
C              POUR LA SUITE
               J0 = N1
               J1 = J2
C
            ELSE
C
C              SUPPRESSION D'UNE ARETE DANS LE CF
               CF(2,J0) = NT
               CF(3,J0) = J2
C              J1 EST VIDE ET RENDUE AU CF
               CF(3,J1) = N1ARCF
               N1ARCF = J1
C              POUR LA SUITE
               J1 = J2
            ENDIF
            GOTO 70
C           FIN DU REPETER JUSQU'A ARETE NON SEPARANTE
         ENDIF
C
C        VERIFICATION QUE LE POINT NP EST BIEN LE SOMMET D'AU
C        MOINS UN TRIANGLE
         IF( NBT .LE. 0 ) THEN
C           NP APPARTIENT A AUCUN TRIANGLE
            IF( CF(3,J2) .EQ. 1 ) THEN
C              TOUT LE CF A ETE PARCOURU SANS TROUVER UN ANGLE NON PLAT
               WRITE(IMPRIM,*) 'TRDVAN:ALGO A REVOIR'
               CALL XVPAUSE
            ENDIF
C           JUSQU'A LA PREMIERE ARETE DU CF FAISANT UN ANGLE NON NUL
            J0 = J1
            J1 = J2
            GOTO 70
         ENDIF
C
C        FIN TANT QUE IL EXISTE UN POINT A AJOUTER
         NOS = NOS + 1
         GOTO 60
      ENDIF
C
 9999 RETURN
      END
