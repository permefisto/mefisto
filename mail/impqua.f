      SUBROUTINE IMPQUA( NUTYOB, NMOBJT, MNNSEF, MNXYZS,
     %                   NBEFMQ, QUAMIN, SURVOLEF )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    AFFICHER L'HISTOGRAMME DES QUALITES DES EF
C -----    LA QUALITE MOYENNE, MINIMALE ET L'ECART TYPE A 1
C          DES QUALITES DES EF DU MAILLAGE D'UNE SURFACE OU D'UN VOLUME

C ENTREES:
C --------
C NUTYOB : NUMERO DU TYPE DE PLSVO
C          =1: point, 2: ligne, 3: surface, 4: volume, 5: objet
C NMOBJT : NOM DU PLSV
C MNNSEF : ADRESSE MCN DU TABLEAU MS DU MAILLAGE
C MNXYZS : ADRESSE MCN DU TABLEAU MS DES SOMMETS

C SORTIES:
C --------
C NBEFMQ : NOMBRE D'EF DE QUALITE INFERIEURE A LA QUALITE MINIMALE QUEFMI
C QUAMIN : QUALITE MINIMALE DES EF DU MAILLAGE
C SURVOLEF : SURFACE des EF ou VOLUME des EF SELON LE TYPE LSV
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS       AVRIL 1991
C....................................................................012
ccc   PARAMETER        (QUEFMI=0.1)
ccc   PARAMETER        (QUEFMI=0.08)
      PARAMETER        (QUEFMI=0.03)
C                       QUALITE MINIMALE SOUHAITEE
      include"./incl/ntmnlt.inc"

C     LE NUMERO D'UNITE DU CLAVIER, FENETRE DES AFFICHAGES ET
C     LE PARAMETRE DE NIVEAU D'INTERACTIVITE
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)

      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      include"./incl/gsmenu.inc"
      include"./incl/langue.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/pp.inc"
      COMMON             MCN(MOTMCN)
      REAL              RMCN(MOTMCN)
      EQUIVALENCE      (MCN(1),RMCN(1))
      CHARACTER*(*)     NMOBJT
      CHARACTER*10      NMTYOB
      CHARACTER*128     KNOM
      CHARACTER*36      TITRE
      INTEGER           NOSOEL(64), NBEFQU(-3:10)
      INTEGER           NBTYEF(0:9)
      REAL              XYZEFS(6,2)

C     PROTECTION CONTRE LES ADRESSES ERRONNEES
      NBEFMQ = 0
      IF( MNNSEF .LE. 0 .OR. MNXYZS .LE. 0 ) RETURN

      WRITE(IMPRIM,*)
      WRITE(IMPRIM,19000)
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,*) 'impqua: la QUALITE des EF du MAILLAGE'
      ELSE
         WRITE(IMPRIM,*) 'impqua: FINITE ELEMENT QUALITIES of the MESH'
      ENDIF

C     DIMENSION DE L'ESPACE DES COORDONNEES
C     PAS D'IMPRESSION POUR LA DIMENSION 6
      NBCOOR = MCN( MNXYZS+WBCOOR )
      IF( NBCOOR .NE. 3 ) RETURN

C     TRACE SEULEMENT DE LA QUALITE DES SURFACES ET DES VOLUMES
      IF( NUTYOB .NE. 3 .AND. NUTYOB .NE. 4 ) RETURN

C     VALIDITE DE LA SURFACE ou VOLUME
      IF( NUTYOB .EQ. 3 ) THEN

C        OUVERTURE DU TMS DE LA SURFACE
         CALL LXLXOU( NTSURF, NMOBJT, NTSF, MNSF   )
         IF( NTSF .LE. 0 ) THEN
            IF( LANGAG .EQ. 0 ) THEN
               PRINT*,'impqua: SURFACE ',NMOBJT,' INCONNUE'
            ELSE
               PRINT*,'impqua: SURFACE ',NMOBJT,' UNKNOWN'
            ENDIF
            RETURN
         ENDIF
         CALL LXTSOU( NTSF,  'NSEF', NTNSEF, MNNSEF )
         IF( NTNSEF .LE. 0 ) THEN
            IF( LANGAG .EQ. 0 ) THEN
               PRINT*,'impqua: SURFACE ',NMOBJT,' TMS NSEF INCONNU'
            ELSE
               PRINT*,'impqua: SURFACE ',NMOBJT,' TMS NSEF UNKNOWN'
            ENDIF
            RETURN
         ENDIF
         CALL LXTSOU( NTSF,  'XYZSOMMET', NTXYZS, MNXYZS )
         IF( NTXYZS .LE. 0 ) THEN
            IF( LANGAG .EQ. 0 ) THEN
               PRINT*,'impqua: SURFACE ',NMOBJT,' TMS XYZSOMMET INCONNU'
            ELSE
               PRINT*,'impqua: SURFACE ',NMOBJT,' TMS XYZSOMMET UNKNOWN'
            ENDIF
            RETURN
         ENDIF

      ELSE

C        OUVERTURE DU TMS DU VOLUME
         CALL LXLXOU( NTVOLU, NMOBJT, NTSF, MNSF   )
         IF( NTSF .LE. 0 ) THEN
            IF( LANGAG .EQ. 0 ) THEN
               PRINT*,'impqua: VOLUME ',NMOBJT,' INCONNUE'
            ELSE
               PRINT*,'impqua: VOLUME ',NMOBJT,' UNKNOWN'
            ENDIF
            RETURN
         ENDIF
         CALL LXTSOU( NTSF,  'NSEF', NTNSEF, MNNSEF )
         IF( NTNSEF .LE. 0 ) THEN
            IF( LANGAG .EQ. 0 ) THEN
               PRINT*,'impqua: VOLUME ',NMOBJT,' TMS NSEF INCONNU'
            ELSE
               PRINT*,'impqua: VOLUME ',NMOBJT,' TMS NSEF UNKNOWN'
            ENDIF
            RETURN
         ENDIF
         CALL LXTSOU( NTSF,  'XYZSOMMET', NTXYZS, MNXYZS )
         IF( NTXYZS .LE. 0 ) THEN
            IF( LANGAG .EQ. 0 ) THEN
               PRINT*,'impqua: VOLUME ',NMOBJT,' TMS XYZSOMMET INCONNU'
            ELSE
               PRINT*,'impqua: VOLUME ',NMOBJT,' TMS XYZSOMMET UNKNOWN'
            ENDIF
            RETURN
         ENDIF

      ENDIF

C     SI NOMBRE INCORRECT DE SOMMETS => RETOUR
      NBSOM = MCN( MNXYZS + WNBSOM )
      IF( NBSOM .LE. 0 ) GOTO 9000

C     LA BOUCLE SUR LES NBCOOR COORDONNEES DES NBSOM SOMMETS
C     CALCUL DU MIN ET MAX DE CHAQUE COORDONNEE
      MN = MNXYZS + WYZSOM -1
      DO 1 I=1,NBCOOR
         XYZEFS(I,1) = RMCN( MN + I )
         XYZEFS(I,2) = RMCN( MN + I )
 1    CONTINUE
      DO 3 J=2,NBSOM
         DO 2 I=1,NBCOOR
C           LA COORDONNEE I DU POINT J
            R = RMCN( MN + I )
            IF( R .LT. XYZEFS(I,1) ) THEN
               XYZEFS(I,1) = R
            ENDIF
            IF( R .GT. XYZEFS(I,2) ) THEN
               XYZEFS(I,2) = R
            ENDIF
 2       CONTINUE
         MN = MN + NBCOOR
 3    CONTINUE

      NBTGS = MCN( MNXYZS + WNBTGS )
      IF( NBTGS .LT. 0 ) GOTO 9000

      IF( INTERA .GE. 1 .AND. LCRITR .EQ. 1 ) THEN
C        ON EFFACE LE TRACE DES QUALITES DU MAILLAGE PRECEDENT
         CALL XVCOULEUR( NCOFON )
ccc         CALL XVRECTANGLE( 0, 0, 280, 300 )
         CALL XVRECTANGLE( 0, 0, 450, 500 )
      ENDIF

C     ADRESSE-3 DE LA 1-ERE COORDONNEE DU 1-ER SOMMET DU TMS 'XYZSOMMET'
      MNS = MNXYZS + WYZSOM

C     LES PARAMETRES DES NO SOMMET DU MAILLAGE
      CALL NSEFPA( MCN(MNNSEF),
     %             NUTYMA, NBSOEL, NBSOEF, NBTGEF,
     %             LDAPEF, LDNGEF, LDTGEF, NBEFOB,
     %             NX    , NY    , NZ    ,
     %             IERR   )
      IF( IERR .NE. 0 ) GOTO 9000

C     INITIALISATIONS
      DO I=-3,10
         NBEFQU(I) = 0
      ENDDO
      QUAMOY = 0.0
      QUAMIN = 2.0
      ECATYP = 0.0
      QUAMOB = 0.0
      QUAMIB = 2.0
      ECATYB = 0.0
      NF     = 0
      CALL AZEROI( 10, NBTYEF )

C     LA BOUCLE SUR LES ELEMENTS FINIS DU MAILLAGE
C     ============================================
      SURVOLEF = 0.
      DO 10 N=1,NBEFOB

C        LE NUMERO DES NBSOEF SOMMETS DU SOUSOBJET N
         CALL NSEFNS( N     , NUTYMA, NBSOEF, NBTGEF,
     %                LDAPEF, LDNGEF, LDTGEF,
     %                MNNSEF, NX, NY, NZ,
     %                NCOGEL, NUGEEF, NUEFTG, NOSOEL, IERR )

         IF( NOSOEL(1) .EQ. 0 ) THEN
C           ELEMENT FINI NON ACTIF
            IF( LANGAG .EQ. 0 ) THEN
               WRITE(IMPRIM,*) 'impqua: EF',N,' de SOMMETS:',
     %                        (NOSOEL(I),I=1,NBSOEF),' N''EST PAS ACTIF'
            ELSE
               WRITE(IMPRIM,*) 'impqua: FE',N,' with VERTICES:',
     %                         (NOSOEL(I),I=1,NBSOEF),' IS NOT ACTIVE'
            ENDIF
            GOTO 10
         ENDIF

C        NOMBRE D'EF SELON LE CODE GEOMETRIQUE
         NBTYEF(NCOGEL) = NBTYEF(NCOGEL) + 1

C        LA QUALITE DE L'EF
         CALL QUALEF( NCOGEL,   NOSOEL, NBSOM, RMCN(MNS),
     %                SURFVOLU, QUALIT, IERR )
         SURVOLEF = SURVOLEF + SURFVOLU

         IF( IERR .NE. 0 ) THEN
            IF( LANGAG .EQ. 0 ) THEN
               WRITE(IMPRIM,*) 'impqua: PROBLEME EF',N,' DE SOMMETS:',
     %                         (NOSOEL(I),I=1,NBSOEF),' <=0'
            ELSE
               WRITE(IMPRIM,*) 'impqua: PROBLEM FE',N,' with VERTICES:',
     %                         (NOSOEL(I),I=1,NBSOEF),' <=0'
            ENDIF
         ENDIF

C        COMPTEUR DU NOMBRE D'EF DE QUALITE A AMELIORER
         IF( QUALIT .LT. QUEFMI ) NBEFMQ = NBEFMQ + 1

C        LA QUALITE MOYENNE
         QUAMOY = QUAMOY + MAX( 0.0, QUALIT )
         ECATYP = ECATYP + (1.0-QUALIT)**2

C        LA QUALITE MINIMALE
         IF( QUALIT .LT. QUAMIN ) THEN
            QUAMIN = QUALIT
            NEFMIN = N
         ENDIF

C        CONTRIBUTION A L'HISTOGRAMME
         I = INT( QUALIT * 10.0 )
         IF( I .LE. 0 ) THEN
            IF( QUALIT .GE. 0.01 ) THEN
               I = 0
            ELSE IF( QUALIT .GE. 0.001 ) THEN
               I = -1
            ELSE IF( QUALIT .GE. 0.0 ) THEN
               I = -2
            ELSE
               I = -3
            ENDIF
         ENDIF
         NBEFQU(I) = NBEFQU(I) + 1

 10   CONTINUE

C     LA QUALITE MOYENNE
      QUAMOY = QUAMOY/NBEFOB

C     L'ECART TYPE A 1
      ECATYP = SQRT( ECATYP/NBEFOB )

C     IMPRESSIONS FINALES
C     ===================
      TITRE = NMTYOB( NUTYOB )
      I     = NUDCNB( TITRE  )
      N     = NUDCNB( NMOBJT )
      TITRE = TITRE(1:I) // ' ' // NMOBJT(1:N)

      IF( LANGAG .EQ. 0 ) THEN

      WRITE(IMPRIM,10010) TITRE,
     %                    NBSOM,NBTGS,NBEFOB
10010 FORMAT(/' QUALITE DU MAILLAGE : ',A,/
     %        ' ---------------------'/
     %' NOMBRE DE SOMMETS       =',I9,/,
     %' NOMBRE DE TANGENTES     =',I9,/,
     %' NOMBRE D''ELEMENTS FINIS =',I9)
      IF( NBTYEF(3) .GT. 0 ) WRITE(IMPRIM,10013) NBTYEF(3)
10013 FORMAT(' NOMBRE DE TRIANGLES     =',I9)
      IF( NBTYEF(4) .GT. 0 ) WRITE(IMPRIM,10014) NBTYEF(4)
10014 FORMAT(' NOMBRE DE QUADRANGLES   =',I9)
      IF( NBTYEF(5) .GT. 0 ) WRITE(IMPRIM,10015) NBTYEF(5)
10015 FORMAT(' NOMBRE DE TETRAEDRES    =',I9)
      IF( NBTYEF(6) .GT. 0 ) WRITE(IMPRIM,10016) NBTYEF(6)
10016 FORMAT(' NOMBRE DE PENTAEDRES    =',I9)
      IF( NBTYEF(9) .GT. 0 ) WRITE(IMPRIM,10019) NBTYEF(9)
10019 FORMAT(' NOMBRE DE PYRAMIDES     =',I9)
      IF( NBTYEF(7) .GT. 0 ) WRITE(IMPRIM,10017) NBTYEF(7)
10017 FORMAT(' NOMBRE DE HEXAEDRES     =',I9)

      WRITE(IMPRIM,10030) MCN(MNNSEF+WBEFTG),
     %                    QUAMOY,QUAMIN,NEFMIN,ECATYP
10030 FORMAT(
     %' NOMBRE D''EF AVEC DES TGS=',I9/,
     %' QUALITE MOYENNE  DES EF =',G12.3,/,
     %' QUALITE MINIMALE DES EF =',G12.3,' POUR L''EF ',I9,/,
     %' ECART TYPE / 1   DES EF =',G12.3)

      ELSE

      WRITE(IMPRIM,11010) TITRE,
     %                    NBSOM,NBTGS,NBEFOB
11010 FORMAT(/' QUALITY of the MESH : ',A,/
     %        ' ---------------------'/
     %' NUMBER of VERTICES        =',I9,/,
     %' NUMBER of TANGENTS        =',I9,/,
     %' NUMBER of FINITE ELEMENTS =',I9)

      IF( NBTYEF(3) .GT. 0 ) WRITE(IMPRIM,11013) NBTYEF(3)
11013 FORMAT(' NUMBER of TRIANGLES       =',I9)
      IF( NBTYEF(4) .GT. 0 ) WRITE(IMPRIM,11014) NBTYEF(4)
11014 FORMAT(' NUMBER of QUADRANGLES     =',I9)
      IF( NBTYEF(5) .GT. 0 ) WRITE(IMPRIM,11015) NBTYEF(5)
11015 FORMAT(' NUMBER of TETRAHEDRA      =',I9)
      IF( NBTYEF(6) .GT. 0 ) WRITE(IMPRIM,11016) NBTYEF(6)
11016 FORMAT(' NUMBER of PENTAHEDRA      =',I9)
      IF( NBTYEF(9) .GT. 0 ) WRITE(IMPRIM,11019) NBTYEF(9)
11019 FORMAT(' NUMBER of PYRAMIDS        =',I9)
      IF( NBTYEF(7) .GT. 0 ) WRITE(IMPRIM,11017) NBTYEF(7)
11017 FORMAT(' NUMBER of HEXAHEDRA       =',I9)

      WRITE(IMPRIM,11030) MCN(MNNSEF+WBEFTG),
     %                    QUAMOY,QUAMIN,NEFMIN,ECATYP
11030 FORMAT(
     %' NUMBER of FE WITH TANGENTS=',I9/,
     %' AVERAGE QUALITY of FE     =',G12.3,/,
     %' MINIMUM QUALITY of FE     =',G12.3,' for the FE ',I9,/,
     %' STANDARD DEVIATION/1 of FE=',G12.3)

      ENDIF

      IF( INTERA .LE. 0 ) GOTO 35

      QUALIB = 1.0
      NBEFQU(9) = NBEFQU(9) + NBEFQU(10)
C
C     TRACE DE LA LEGENDE DES QUALITES SUR ECRAN ET POSTSCRIPT
C     --------------------------------------------------------
      NNX =  35
      NNY = 120
      IF( INTERA .GE. 1 .AND. LCRITR .EQ. 1 ) THEN

C        EFFACEMENT DE LA LEGENDE DE LA QUALITE SUR POSTSCRIPT
         IF ( LASOPS.NE.0 ) THEN
           IF ( LASOPS.EQ.1 ) THEN
             LASOPS = -11
           ELSE
             IF ( LASOPS.EQ.2 ) THEN
               LASOPS = -12
             ELSE
               LASOPS = 0
               NBLGRC(NRERR) = 2
               KERR(1) = 'impqua: MAUVAISE VALEUR de LASOPS'
               KERR(2) = '        ARRET du TRACE POSTSCRIPT'
               CALL LEREUR
             ENDIF
           ENDIF
           CALL XVPOSTSCRIPT(LASOPS)
           LASOPS = - LASOPS
           CALL XVPOSTSCRIPT(LASOPS)
         ENDIF

C        PLACEMENT DE LA LEGENDE SUR L'ECRAN
         IF( LANGAG .EQ. 0 ) THEN
            KNOM = ' QUALITE des ELEMENTS FINIS'
         ELSE
            KNOM = ' QUALITY of FINITE ELEMENTS'
         ENDIF
         CALL XVCOULEUR( NCROUG )
         CALL XVTEXTE( KNOM(1:27), 27, NNX, NNY )
         NNY = NNY + 16
      ENDIF

      DO 20 N=9,-3,-1

C        L'INTERVALLE N DES QUALITES
         QUALIT = N/10.0 + 1E-6
         NN = 0
         IF( N .EQ. 0 ) THEN
            QUALIT = 0.01
            NN = 1
         ELSE IF( N .EQ. -1 ) THEN
            QUALIT = 0.001
            NN = 2
         ELSE IF( N .EQ. -2 ) THEN
            QUALIT = 0.0
            NN = 2
         ELSE IF( N .EQ. -3 ) THEN
            QUALIT = -1.0
            NN = 3
         ENDIF

C        LE NOMBRE D'EF POUR CETTE QUALITE
         IF( NBEFQU(N) .LE. 0 ) GOTO 19

C        LE POURCENTAGE ARRONDI D'EF POUR CETTE QUALITE
         I = NINT( ( NBEFQU(N) * 100.0 ) / NBEFOB )

         IF( LANGAG .EQ. 0 ) THEN

         IF((I.GE.2).OR.(NBEFQU(N).EQ.0)) THEN
           WRITE(IMPRIM,10020) QUALIT,QUALIB,NBEFQU(N),I,('*',J=1,I/2)
           WRITE(KNOM,  10020) QUALIT,QUALIB,NBEFQU(N),I,('*',J=1,I/2)
         ELSE IF(I.GT.0) THEN
           WRITE(IMPRIM,10020) QUALIT,QUALIB,NBEFQU(N),I,'+'
           WRITE(KNOM,  10020) QUALIT,QUALIB,NBEFQU(N),I,'+'
         ELSE
           WRITE(IMPRIM,10020) QUALIT,QUALIB,NBEFQU(N),I,'.'
           WRITE(KNOM,  10020) QUALIT,QUALIB,NBEFQU(N),I,'.'
         ENDIF

         ELSE

         IF((I.GE.2).OR.(NBEFQU(N).EQ.0)) THEN
           WRITE(IMPRIM,11020) QUALIT,QUALIB,NBEFQU(N),I,('*',J=1,I/2)
           WRITE(KNOM,  11020) QUALIT,QUALIB,NBEFQU(N),I,('*',J=1,I/2)
         ELSE IF(I.GT.0) THEN
           WRITE(IMPRIM,11020) QUALIT,QUALIB,NBEFQU(N),I,'+'
           WRITE(KNOM,  11020) QUALIT,QUALIB,NBEFQU(N),I,'+'
         ELSE
           WRITE(IMPRIM,11020) QUALIT,QUALIB,NBEFQU(N),I,'.'
           WRITE(KNOM,  11020) QUALIT,QUALIB,NBEFQU(N),I,'.'
         ENDIF

         ENDIF
10020 FORMAT(F6.3,'=< QUALITE <',F5.3,' :',I9,' EF soit ',
     %       I3,' % : ',100(A1))
11020 FORMAT(F6.3,'=< QUALITY <',F5.3,' :',I9,' FE ',
     %       I3,' % : ',100(A1))

         IF( INTERA .GE. 1 .AND. LCRITR .EQ. 1 ) THEN
C           TRACE EN HAUT ET A GAUCHE
            CALL XVCOULEUR( NCOQUA(QUALIT) )

ccc            IF( LANGAG .EQ. 0 ) THEN
ccc               WRITE(KNOM,10021) NBEFQU(N),I,QUALIT,QUALIB
ccc            ELSE
ccc               WRITE(KNOM,11021) NBEFQU(N),I,QUALIT,QUALIB
ccc            ENDIF
ccc            IF( NN .LE. 1 ) KNOM(19:19) = ' '
ccc            IF( NN .LE. 0 ) KNOM(18:18) = ' '

            IF( N .EQ. 9 ) KNOM(17:17) = '='
            NBC = NUDCNB( KNOM )
            CALL XVTEXTE( KNOM(1:NBC), NBC, NNX, NNY )
            NNY = NNY + 16
         ENDIF

ccc10021 FORMAT(I9,' EF',I3,' %',F5.3,'< Q <=',F5.3)
10022 FORMAT('  QUALITE MOYENNE =',G10.3)
10023 FORMAT('  QUALITE MINIMALE=',G10.3)
10024 FORMAT('  ECART TYPE A 1  =',G10.3)
10025 FORMAT('  NOMBRE SOMMETS  =',I9)
10026 FORMAT('  NOMBRE TANGENTES=',I9)
10027 FORMAT('  NOMBRE TOTAL EF =',I9)

ccc11021 FORMAT(I9,' FE',I3,' %',F5.3,'< Q <=',F5.3)
11022 FORMAT('  FE''s AVERAGE QUALITY=',G10.3)
11023 FORMAT('  FE''s MINIMUM QUALITY=',G10.3)
11024 FORMAT('  STANDARD DEVIATION/1=',G10.3)
11025 FORMAT('  NUMBER of VERTICES  =',I9)
11026 FORMAT('  NUMBER of TANGENTS  =',I9)
11027 FORMAT('  TOTAL NUMBER of FE  =',I9)
C
 19      QUALIB = QUALIT

 20   CONTINUE

C     TRACE DU TEXTE EN GRIS FONCE
      CALL XVCOULEUR( NCGRIS )
      IF( INTERA .GE. 1 .AND. LCRITR .EQ. 1 ) THEN
C        TRACE DE LA QUALITE MOYENNE DU MAILLAGE
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(KNOM,10022) QUAMOY
         ELSE
            WRITE(KNOM,11022) QUAMOY
         ENDIF
         NBC = NUDCNB( KNOM )
         CALL XVTEXTE( KNOM(1:NBC), NBC, NNX, NNY )
         NNY = NNY + 16
C        TRACE DE LA QUALITE MINIMALE DU MAILLAGE
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(KNOM,10023) QUAMIN
         ELSE
            WRITE(KNOM,11023) QUAMIN
         ENDIF
         NBC = NUDCNB( KNOM )
         CALL XVTEXTE( KNOM(1:NBC), NBC, NNX, NNY )
         NNY = NNY + 16
C        TRACE DE L'ECART TYPE A 1
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(KNOM,10024) ECATYP
         ELSE
            WRITE(KNOM,11024) ECATYP
         ENDIF
         NBC = NUDCNB( KNOM )
         CALL XVTEXTE( KNOM(1:NBC), NBC, NNX, NNY )
         NNY = NNY + 25

C        TRACE DU NOMBRE TOTAL D'EF
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(KNOM,10027) NBEFOB
         ELSE
            WRITE(KNOM,11027) NBEFOB
         ENDIF
         CALL SANSDBL( KNOM, NBC )
         CALL XVTEXTE( KNOM(1:NBC), NBC, NNX, NNY )
         NNY = NNY + 16
C        TRACE DU NOMBRE DE SOMMETS
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(KNOM,10025) NBSOM
         ELSE
            WRITE(KNOM,11025) NBSOM
         ENDIF
         CALL SANSDBL( KNOM, NBC )
         CALL XVTEXTE( KNOM(1:NBC), NBC, NNX, NNY )
         NNY = NNY + 16
C        TRACE DU NOMBRE DE TANGENTES
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(KNOM,10026) NBTGS
         ELSE
            WRITE(KNOM,11026) NBTGS
         ENDIF
         CALL SANSDBL( KNOM, NBC )
         CALL XVTEXTE( KNOM(1:NBC), NBC, NNX, NNY )

C        RETOUR AU TRACE NORMAL POUR POSTSCRIPT
         IF ( LASOPS.NE.0 ) THEN
           LASOPS = LASOPS -10
           CALL XVPOSTSCRIPT(LASOPS)
         ENDIF
      ENDIF

C     AFFICHAGE DES COORDONNEES DE L'EF DE QUALITE MINIMALE
 35   IF( QUAMIN .LT. 0.3 ) THEN
C
C        LE NUMERO DES NBSOEF SOMMETS DE L'EF NEFMIN
         CALL NSEFNS( NEFMIN, NUTYMA, NBSOEF, NBTGEF,
     %                LDAPEF, LDNGEF, LDTGEF,
     %                MNNSEF, NX, NY, NZ,
     %                NCOGEL, NUGEEF, NUEFTG, NOSOEL, IERR )

         NBS = NBSOEF
 38      IF( NOSOEL(NBS) .LE. 0 ) THEN
            NBS = NBS - 1
            GOTO 38
         ENDIF
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10040) NEFMIN,QUAMIN,NBS
10040       FORMAT(/' EF',I9,' de QUALITE MINIMALE ',G10.3,
     %              ' de',I2,' NOEUDS')
         ELSE
            WRITE(IMPRIM,11040) NEFMIN,QUAMIN,NBS
11040       FORMAT(/' FE',I9,' of MINIMUM QUALITY ',G10.3,
     %              ' of ',I2,' NODES')
         ENDIF

         DO 50 J=1,NBS
C           LE NUMERO DU NOEUD J DE L'EF
            N = NOSOEL(J)
            IF( N .LE. 0 ) GOTO 50
            MN = MNS+3*N-4
            IF( LANGAG .EQ. 0 ) THEN
               WRITE(IMPRIM,10050) N,(RMCN(MN+I),I=1,3)
            ELSE
               WRITE(IMPRIM,11050) N,(RMCN(MN+I),I=1,3)
            ENDIF
   50    CONTINUE
10050 FORMAT(' SOMMET',I9,' X=',G15.7,' Y=',G15.7,' Z=',G15.7)
11050 FORMAT(' VERTEX',I9,' X=',G15.7,' Y=',G15.7,' Z=',G15.7)
      ENDIF

C     AFFICHAGE DU MIN ET MAX DE CHAQUE COORDONNEE
ccc      WRITE(IMPRIM,*)
      DO 4 I=1,NBCOOR
         WRITE(IMPRIM,10004) I, XYZEFS(I,1), XYZEFS(I,2)
 4    CONTINUE
10004 FORMAT(' COORD',I1,' MIN=',G15.7,' MAX=',G15.7)

      NBCSV = NUDCNB( NMOBJT )
      IF( NUTYOB .EQ. 3 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            PRINT*, 'SURFACE des EF=',SURVOLEF
         ELSE
            PRINT*, 'FE SURFACE=',SURVOLEF
         ENDIF
      ENDIF

      IF( NUTYOB .EQ. 4 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            PRINT*, 'VOLUME des EF=', SURVOLEF
         ELSE
            PRINT*, 'FE VOLUME=', SURVOLEF
         ENDIF
      ENDIF

C     LE TYPE DE FERMETURE DU MAILLAGE DE LA LIGNE OU SURFACE
      NUTFMA = MCN( MNNSEF + WUTFMA )

      IF( NUTYOB .EQ. 3 .AND. NUTFMA .EQ. -1 ) THEN
C        RECHERCHER SI LE MAILLAGE DE LA SURFACE EST FERMEE OU NON
C        FERMEE: SI TOUTE ARETE D'UNE FACE APPARTIENT A 2 FACES
         CALL NUOBNM('SURFACE', NMOBJT, NUSURF )
         CALL OBJFER( NUTYOB, NUSURF, 0, NUTFMA )
      ENDIF

      IF( NUTFMA .EQ. 1 ) THEN
         IF( NUTYOB .EQ. 3 ) THEN
C           SURFACE FERMEE (TOUTE ARETE APPARTIENT A 2 FACES)
            IF( LANGAG .EQ. 0 ) THEN
               PRINT*, NMOBJT(1:NBCSV),
     %' est une SURFACE FERMEE de R**3  ...............................'
            ELSE
               PRINT*, NMOBJT(1:NBCSV),
     %' is a CLOSED SURFACE of R**3  ..................................'
            ENDIF
         ENDIF
      ENDIF

C     AFFICHAGE D'UNE ERREUR SI LA QUALITE MINIMALE QUAMIN EST TROP PETITE
C     --------------------------------------------------------------------
      IF( QUAMIN .LE. 0.01 ) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(5)(1:15),'(G15.6)') QUAMIN
         KERR(1) = 'impqua: ' // NMTYOB( NUTYOB ) // ' ' // NMOBJT
         IF( LANGAG .EQ. 0 ) THEN
            KERR(2) ='MAILLAGE avec une TRES MAUVAISE QUALITE MINIMALE='
     %                 // KERR(5)(1:15)
         ELSE
            KERR(2) = 'MESH has a VERY BAD MINIMUM QUALITY='
     %                 // KERR(5)(1:15)
         ENDIF
         CALL SANSDBL( KERR(2), LL )
         CALL SANSDBL( KERR(1), L  )
         KERR(1) = KERR(1)(1:L) // ' : ' // KERR(2)(1:LL)
         CALL SANSDBL( KERR(1), L )
ccc         CALL LAVERT     TROP PENIBLE lors des traces...
         PRINT*
         PRINT*,KERR(1)(9:L)
      ENDIF

C     MISE A JOUR DU CADRE MAXIMAL DES MAILLAGES
C     ------------------------------------------
      CALL MAJXYZEXT( NBCOOR, XYZEFS )

C     SEPARATEUR DE FIN DE PLSV CREE
 9000 WRITE(IMPRIM,19000)
19000 FORMAT(100('+'))

      RETURN
      END
