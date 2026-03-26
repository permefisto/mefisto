       SUBROUTINE ROTAN2( NBANGL, ANGLES, PT1,    PT2,
     %                    FERME,  NBSOLI, SOMLIG, NBSOAX, NOSOAX,
     %                    NUTYML, NBARLI, NSARLI,
     %                    NBTGLI, XYZTGL, NUTGLI,
     %                    NBSOSU, SOMSUR, NUTYMS, NSFASU,
     %                    NBTGSU, XYZTGS, NUTGSU,
     %                    IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    GENERATION DU MAILLAGE DE LA SURFACE PAR ROTATION D'UNE LIGNE
C -----
C
C ENTREES :
C ---------
C NBANGL  : LE NOMBRE DES ANGLES DE ROTATION
C ANGLES  : LES NBANGL ANGLES DE ROTATION( EN DEGRES )
C PT1 PT2 : LES 3 COORDONNEES DES 2 POINTS DE DEFINITION DE L'AXE
C FERME   : VRAI SI LA SOMME DES ANGLES VAUT 360 DEGRES, FAUX SINON
C NBSOLI  : LE NOMBRE DE SOMMETS DE LA LIGNE
C SOMLIG  : LES 3 COORDONNEES DES NBSOLI SOMMETS DE LA LIGNE
C NBSOAX  : LE NOMBRE DE SOMMETS DE LA LIGNE SUR L'AXE
C NOSOAX  : +LE NUMERO DU SOMMET S'IL N'EST PAS SUR L'AXE PT1-PT2
C           -LE NUMERO DU SOMMET S'IL   EST     SUR L'AXE PT1-PT2
C NUTYML  : NUMERO TYPE DU MAILLAGE DE LA LIGNE
C           0 MAILLAGE NON STRUCTURE => TABLEAU NSARLI     NECESSAIRE
C           2 MAILLAGE     STRUCTURE => TABLEAU NSARLI NON NECESSAIRE
C NBARLI  : LE NOMBRE D'ARETES DE LA LIGNE
C NSARLI  : LE NUMERO DES 2 SOMMETS DE CHAQUE ARETE DE LA LIGNE (1 A NBSOLI)
C NBTGLI  : LE NOMBRE DE TANGENTES STOCKEES DE LA LIGNE
C XYZTGL  : LES 3 COMPOSANTES DES TANGENTES DE LA LIGNE
C NUTGLI  : LES 2 NUMEROS DES TANGENTES DES NBARLI ARETES DE LA LIGNE
C NBSOSU  : LE NOMBRE DE SOMMETS DE LA SURFACE MAILLEE
C NUTYMS  : LE NUMERO TYPE DU MAILLAGE DE LA SURFACE
C           0 MAILLAGE NON STRUCTURE => TABLEAU NSFASU     GENERE
C           4 MAILLAGE     STRUCTURE => TABLEAU NSFASU NON GENERE
C
C SORTIES :
C ---------
C SOMSUR  : LES 3 COORDONNEES DES NBSOSU SOMMETS DE LA SURFACE MAILLEE
C NSFASU  : LES 4 SOMMETS DES FACES DE LA SURFACE
C NBTGSU  : LE NOMBRE DE TANGENTES DE LA SURFACE MAILLEE
C XYZTGS  : LES 3 COMPOSANTES DES VECTEURS TANGENTS
C NUTGSU  : LES 8 NUMEROS DES TANGENTES DE CHACUN DES EF DE LA SURFACE
C IERR    : 1 SI TOUTE LA LIGNE SUR L'AXE PT1-PT2 => SORTIE SANS FACES
C           2 SI  UNE ARETE EST SUR L'AXE PT1-PT2 => SORTIE SANS FACES
C           0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS   NOVEMBRE 1996
C23456---------------------------------------------------------------012
      include"./incl/gsmenu.inc"
      COMMON / EPSSSS / EPZERO,EPSXYZ
      REAL              ANGLES(NBANGL), PT1(3), PT2(3),
     %                  SOMLIG(3,NBSOLI), SOMSUR(3,NBSOSU),
     %                  XYZTGL(3,NBTGLI),
     %                  XYZTGS(3,NBTGSU)
      INTEGER           NOSOAX(NBSOLI),
     %                  NSARLI(2,NBARLI),
     %                  NUTGLI(2,NBARLI),
     %                  NSFASU(4,NBARLI,NBANGL),
     %                  NUTGSU(8,NBARLI,NBANGL)
C
      REAL              XYZ(3), XYZL(3)
      LOGICAL           FERME
      DOUBLE PRECISION  D2D3(3,3)
C
C     PI / 180 => PIS180
      PIS180 = ATAN(1.) / 45.0
C
C     CHANGEMENT DE REPERE PT1: PT1-PT2, PT1-1-ER POINT DE LA LIGNE
C                                            NON COLINEAIRE A L'AXE
      COSMIN = 2
      N      = 0
      DO 5 I=1,NBSOLI
         IF( NOSOAX(I) .LT. 0 ) GOTO 5
C        RECHERCHE DU COSINUS MINIMAL EN VALEUR ABSOLUE
         COSI = ABS( COS3PT( PT1, PT2, SOMLIG(1,I) ) )
         IF( COSI .LT. COSMIN ) THEN
C           MEILLEUR MINIMUM
            COSMIN = COSI
            N      = I
            IF( COSMIN .LT. 0.2 ) GOTO 8
         ENDIF
 5    CONTINUE
      IF( N .EQ. 0 ) GOTO 9900
C
 8    CALL DF3D2D( PT1 , PT2 , SOMLIG(1,N), D2D3 , IERR )
      IF( IERR .LT. 0 ) THEN
         NBLGRC(NRERR) = 2
         KERR(1) = 'AXE ET LIGNE PRESQUE ALIGNES'
         KERR(2) = 'RESULTATS A VERIFIER'
         CALL LEREUR
      ENDIF
      IF( IERR .GT. 0 ) GOTO 9900
C
C     LA LIGNE N'EST PAS COLINEAIRE A PT1-PT2
C     LES COORDONNEES DES SOMMETS DE LA 1-ERE LIGNE DANS LA SURFACE
      CALL TRTATA( SOMLIG , SOMSUR , 3*NBSOLI )
C
C     LE TABLEAU NOSOAX EST MODIFIE
C     NOSOAX(SOMMET NON SUR L'AXE)=+SON NUMERO DE 1 A NBSOLI-NBSOAX
C     NOSOAX(SOMMET     SUR L'AXE)=-SON NUMERO DE 1 A NBSOAX
      ND1  = 0
      NBSS = 0
      DO 30 N=1,NBSOLI
         IF( NOSOAX(N) .LT. 0 ) THEN
            IF( ND1 .EQ. 1 ) THEN
C              2 POINTS SUR L'AXE SUCCESSIFS +> INTERDIT
               NBLGRC(NRERR) = 2
               KERR(1) = '2 SOMMETS SUCCESSIFS DE LA LIGNE SUR L''AXE'
               KERR(2) = 'SITUATION INTERDITE'
               CALL LEREUR
               RETURN
            ENDIF
            NBSS = NBSS + 1
            NOSOAX(N) = -NBSS
            ND1 = 1
         ELSE
            NOSOAX(N) = N - NBSS
            ND1 = 0
         ENDIF
 30   CONTINUE
C
C     LES NBANGL ROTATIONS
      A    = 0
      COSA = 1
      SINA = 0
      NBSS   = NBSOLI
      NBSOTR = NBSOLI - NBSOAX
      ND1    = 0
      ND2    = NBSOLI
      NBTGSU = 0
C
      DO 100 N=1,NBANGL
C
C        L'ANGLE N ACTUEL DE ROTATION
         A0    = A
         COSA0 = COSA
         SINA0 = SINA
C
C        L'ANGLE SUIVANT
         A = A + ANGLES(N)
C
C        LA LONGUEUR DE L'ARC DE CERCLE
         R = ANGLES(N) * PIS180
C
C        LE COSINUS ET SINUS
         ARAD = A * PIS180
         COSA = COS( ARAD )
         SINA = SIN( ARAD )
C
C        CAS PARTICULIER DE LA SURFACE FERMEE ET DU DERNIER SECTEUR
C        POUR LEQUEL LA DERNIERE LIGNE EST EN FAIT LA PREMIERE
         IF( FERME .AND. N.EQ.NBANGL ) GOTO 60
C
C        LES COORDONNEES DES SOMMETS DE LA LIGNE DE ROTATION
C        ===================================================
         DO 50 I=1,NBSOLI
            IF( NOSOAX(I) .GT. 0 ) THEN
C              UN SOMMET A AJOUTER CAR NON SUR L'AXE
               NBSS = NBSS + 1
C              LE SOMMET DE LA LIGNE EST RAMENE DANS LE NOUVEAU REPERE
               CALL CH3D3D( PT1, D2D3, SOMLIG(1,I), SOMSUR(1,NBSS) )
C              LA ROTATION AUTOUR DE L'AXE PT1-PT2 CAD X
C              DANS LE NOUVEAU REPERE L'AXE DE ROTATION EST OX
               S2 = SOMSUR(2,NBSS)
               S3 = SOMSUR(3,NBSS)
               XYZ(1) = SOMSUR(1,NBSS)
               XYZ(2) = COSA * S2 - SINA * S3
               XYZ(3) = SINA * S2 + COSA * S3
C              RETOUR AU REPERE INITIAL
               CALL CH3D3R( PT1, D2D3, XYZ, SOMSUR(1,NBSS) )
            ENDIF
 50      CONTINUE
C
C        LES TANGENTES DES FACES ENGENDREES PAR CETTE ROTATION
C        =====================================================
 60      DO 80 I=1,NBARLI
C
C           CALCUL DES 8 TANGENTES DE L'EF (I,N)
            CALL AZEROI( 8, NUTGSU(1,I,N) )
C
c           LES 4 TANGENTES A LA LIGNE GENERATRICE POUR CET EF (I,N)
            DO 70 K=1,2
C
C              ROTATION DES TANGENTES DE LA LIGNE
C              ----------------------------------
C              LE NUMERO DE LA TANGENTE K DE L'ARETE I DE LA LIGNE GENERATRICE
               IF( NBTGLI .GT. 0 ) THEN
C                 IL EXISTE DES TANGENTES : LE NUMERO DE LA TANGENTE
                  NT = NUTGLI(K,I)
               ELSE
C                 PAS DE TANGENTE
                  NT = 0
               ENDIF
               IF( NT .NE. 0 ) THEN
C
C                 LE SIGNE DU NUMERO DE LA TANGENTE
                  IF( NT .LT. 0 ) THEN
                     LESIGN = -1
                     NT     = -NT
                  ELSE
                     LESIGN = 1
                  ENDIF
C
C                 LE VECTEUR TANGENT DE LA LIGNE EST RAMENE DANS LE NOUVEAU REPE
                  CALL CV3D3D( D2D3, XYZTGL(1,NT), XYZL )
                  S2 = XYZL(2)
                  S3 = XYZL(3)
C
C                 LA ROTATION AUTOUR DE L'AXE PT1-PT2 D'UN ANGLE A0
C                 DANS LE NOUVEAU REPERE L'AXE DE ROTATION EST OX
                  XYZ(1) = XYZL(1)
                  XYZ(2) = COSA0 * S2 - SINA0 * S3
                  XYZ(3) = SINA0 * S2 + COSA0 * S3
C
C                 RETOUR AU REPERE INITIAL DE CETTE NOUVELLE TANGENTE
                  NBTGSU = NBTGSU + 1
                  CALL CV3D3R( D2D3, XYZ, XYZTGS(1,NBTGSU) )
                  IF( K .EQ. 1 ) THEN
                     NUTEF = 1
                  ELSE
                     NUTEF = 4
                  ENDIF
                  NUTGSU( NUTEF, I, N ) = LESIGN * NBTGSU
C
C                 LA ROTATION AUTOUR DE L'AXE PT1-PT2 D'UN ANGLE A
                  XYZ(1) = XYZL(1)
                  XYZ(2) = COSA * S2 - SINA * S3
                  XYZ(3) = SINA * S2 + COSA * S3
C
C                 RETOUR AU REPERE INITIAL DE CETTE NOUVELLE TANGENTE
                  NBTGSU = NBTGSU + 1
                  CALL CV3D3R( D2D3, XYZ, XYZTGS(1,NBTGSU) )
                  IF( K .EQ. 1 ) THEN
                     NUTEF = 8
                  ELSE
                     NUTEF = 5
                  ENDIF
                  NUTGSU( NUTEF, I, N ) = LESIGN * NBTGSU
               ENDIF
C
 70         CONTINUE
C
C           LES AUTRES TANGENTES A LA ROTATION
C           ----------------------------------
            IF( NOSOAX(I) .GT. 0 ) THEN
C              LE SOMMET DE LA LIGNE EST RAMENE DANS LE NOUVEAU REPERE
               CALL CH3D3D( PT1, D2D3, SOMLIG(1,I), XYZL )
               S2 = XYZL(2)
               S3 = XYZL(3)
C
C              LE VECTEUR TANGENT A LA ROTATION D'ANGLE A0
               XYZ(1) = 0
               XYZ(2) = (-SINA0 * S2 - COSA0 * S3) * R
               XYZ(3) = ( COSA0 * S2 - SINA0 * S3) * R
C              RETOUR AU REPERE INITIAL
               NBTGSU = NBTGSU + 1
               CALL CV3D3R( D2D3, XYZ, XYZTGS(1,NBTGSU) )
               NUTGSU( 2, I, N ) = NBTGSU
C
C              MOINS LE VECTEUR TANGENT DU A LA ROTATION D'ANGLE A
               XYZ(1) = 0
               XYZ(2) = ( SINA * S2 + COSA * S3) * R
               XYZ(3) = (-COSA * S2 + SINA * S3) * R
C              RETOUR AU REPERE INITIAL
               NBTGSU = NBTGSU + 1
               CALL CV3D3R( D2D3, XYZ, XYZTGS(1,NBTGSU) )
               NUTGSU( 7, I, N ) = NBTGSU
C
            ELSE
C
C              LE SOMMET 1 ET 4 DU QUADRANGLE FORME LE SOMMET 1 DU TRIANGLE
C              LA TG 8 DU QUADRANGLE EST LA TANGENTE 2 DU TRIANGLE
               NUTGSU( 2, I, N ) = NUTGSU( 8, I, N )
               NUTGSU( 7, I, N ) = 0
               NUTGSU( 8, I, N ) = 0
C
            ENDIF
C
            IF( NOSOAX(I+1) .GT. 0 ) THEN
C              LE SOMMET DE LA LIGNE EST RAMENE DANS LE NOUVEAU REPERE
               CALL CH3D3D( PT1, D2D3, SOMLIG(1,I+1), XYZL )
               S2 = XYZL(2)
               S3 = XYZL(3)
C
C              LE VECTEUR TANGENT A LA ROTATION D'ANGLE A0
               XYZ(1) = 0
               XYZ(2) = (-SINA0 * S2 - COSA0 * S3) * R
               XYZ(3) = ( COSA0 * S2 - SINA0 * S3) * R
C              RETOUR AU REPERE INITIAL
               NBTGSU = NBTGSU + 1
               CALL CV3D3R( D2D3, XYZ, XYZTGS(1,NBTGSU) )
               NUTGSU( 3, I, N ) = NBTGSU
C
C              MOINS LE VECTEUR TANGENT DU A LA ROTATION D'ANGLE A
               XYZ(1) = 0
               XYZ(2) = ( SINA * S2 + COSA * S3) * R
               XYZ(3) = (-COSA * S2 + SINA * S3) * R
C              RETOUR AU REPERE INITIAL
               NBTGSU = NBTGSU + 1
               CALL CV3D3R( D2D3, XYZ, XYZTGS(1,NBTGSU) )
               NUTGSU( 6, I, N ) = NBTGSU
C
            ELSE
C
C              LE SOMMET 2 ET 3 DU QUADRANGLE FORME LE SOMMET 2 DU TRIANGLE
C              LA TG 5 DU QUADRANGLE EST LA TANGENTE 3 DU TRIANGLE
C              LA TG 7 DU QUADRANGLE EST LA TANGENTE 5 DU TRIANGLE
C              LA TG 8 DU QUADRANGLE EST LA TANGENTE 6 DU TRIANGLE
               NUTGSU( 3, I, N ) = NUTGSU( 5, I, N )
               NUTGSU( 5, I, N ) = NUTGSU( 7, I, N )
               NUTGSU( 6, I, N ) = NUTGSU( 8, I, N )
               NUTGSU( 7, I, N ) = 0
               NUTGSU( 8, I, N ) = 0
C
            ENDIF
C
 80      CONTINUE
C
         IF( NUTYMS .EQ. 0 ) THEN
C
C           LA SURFACE N'EST PAS UN QUADRANGLE STRUCTURE
C           CAR IL EXISTE AU MOINS UN POINT SUR L'AXE OU BIEN
C           LA SURFACE SE REFERME PAR ROTATION DE 360 DEGRES
            DO 90 I=1,NBARLI
C              LE NUMERO DES 2 SOMMETS DE CETTE ARETE
               IF( NUTYML .EQ. 2 ) THEN
                  NS1 = I
                  NS2 = I + 1
               ELSE
                  NS1 = NSARLI(1,I)
                  NS2 = NSARLI(2,I)
               ENDIF
C              LE NUMERO DES 2 AUTRES SOMMETS DU QUADRANGLE 1=>4 2=>3
               NS4 = NOSOAX( NS1 )
               NS3 = NOSOAX( NS2 )
C
               IF( N .EQ. 1 ) THEN
C                 LA PREMIERE TRANCHE
                  IF( NS4 .LT. 0 ) THEN
                     IF( NS3 .LT. 0 ) THEN
C                       ARETE SUR L'AXE INTERDITE
                        NBLGRC(NRERR) = 1
                        KERR(1) = 'ARETE SUR L''AXE INTERDITE'
                        CALL LEREUR
                        IERR = 2
                        RETURN
                     ELSE
C                       TRIANGLE S1 S2 S3 AVEC S1 SUR L'AXE
C                       NS1 = NS1
C                       NS2 = NS2
                        NS3 = NOSOAX(NS2) + NBSOLI
                        NS4 = 0
                     ENDIF
                  ELSE
                     IF( NS3 .LT. 0 ) THEN
C                       TRIANGLE S1 S2 S4 AVEC S2 SUR L'AXE
C                       NS1 = NS1
C                       NS2 = NS2
                        NS3 = NOSOAX(NS1) + NBSOLI
                        NS4 = 0
                     ELSE
C                       QUADRANGLE 12 34
C                       NS1 = NS1
C                       NS2 = NS2
                        NS3 = NOSOAX(NS2) + NBSOLI
                        NS4 = NOSOAX(NS1) + NBSOLI
                     ENDIF
                  ENDIF
               ELSE
C                 LES TRANCHES SUIVANTES
                  IF( NS4 .LT. 0 ) THEN
                     IF( NS3 .GT. 0 ) THEN
C                       TRIANGLE S1 S2 S3 AVEC S1 SUR L'AXE
C                       NS1 = NS1
                        NS2 = NS3 + ND1
                        NS3 = NS3 + ND2
                        NS4 = 0
                     ENDIF
                  ELSE
                     IF( NS3 .LT. 0 ) THEN
C                       TRIANGLE S1 S2 S4 AVEC S2 SUR L'AXE
                        NS1 = NS4 + ND1
C                       NS2 = NS2
                        NS3 = NS4 + ND2
                        NS4 = 0
                     ELSE
C                       QUADRANGLE 12 34
                        NS1 = NS4 + ND1
                        NS2 = NS3 + ND1
                        NS3 = NS3 + ND2
                        NS4 = NS4 + ND2
                     ENDIF
                  ENDIF
               ENDIF
C
               NSFASU(1,I,N) = NS1
               NSFASU(2,I,N) = NS2
               NSFASU(3,I,N) = NS3
               NSFASU(4,I,N) = NS4
 90         CONTINUE
C           LE DECALAGE POUR LA TRANCHE SUIVANTE
            ND1 = ND2
            ND2 = ND2 + NBSOTR
         ENDIF
 100  CONTINUE
C
C     EN CAS DE FERMETURE REMISE A JOUR DE LA DERNIERE LIGNE=PREMIERE LIGNE
C     =====================================================================
      IF( FERME ) THEN
         DO 200 I=1,NBARLI
            IF( NOSOAX(NSFASU(2,I,1)) .LT. 0 ) THEN
C              LE 3-EME SOMMET DE LA DERNIERE LIGNE EST LE 1-EME DE LA PREMIERE
               NSFASU(3,I,NBANGL) = NSFASU(1,I,1)
            ELSE
C              LE 3-EME SOMMET DE LA DERNIERE LIGNE EST LE 2-EME DE LA PREMIERE
               NSFASU(3,I,NBANGL) = NSFASU(2,I,1)
            ENDIF
C           LE 4-EME SOMMET DE LA DERNIERE LIGNE EST LE 1-EME DE LA PREMIERE
            IF( NSFASU(4,I,NBANGL) .NE. 0 ) THEN
               NSFASU(4,I,NBANGL) = NSFASU(1,I,1)
            ENDIF
 200     CONTINUE
      ENDIF
      RETURN
C
C     ERREUR
 9900 NBLGRC(NRERR) = 1
      KERR(1) = 'AXE ET LIGNE ALIGNES'
      CALL LEREUR
      IERR = 2
      END
