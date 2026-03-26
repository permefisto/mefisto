      SUBROUTINE OTAJSF( NUVOLU, KNMVOL, HEXAED, NBVOPA, NUVOPA,
     %                   NBSOMM, MXSOMM, PTXYZD, NPSOFR,
     %                   MXARBO, LARBRO, MXARBT, LARBRT, NUOTPT,
     %                   LAVIDE, MXARET, LARETE,
     %                   NBFACO, MXFACO, LEFACO, N1FASC,
     %                   MXAREF, LARETF,
     %                   QTRMIN, QTRMOY, AREMIN, AREMAX, AREMOY,
     %                   COSMIN, COSMAX, IERR  )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    AJOUTER LES SOMMETS DES TRIANGLES DES SURFACES CONTOURS
C -----    METTRE A JOUR NPSOFR NUMERO DE SURFACE DES SOMMETS DES FACES
C          (1 000 000) * (NO SURFACE + 1) + NUMERO DU SOMMET DANS
C          LA NUMEROTATION DES SOMMETS DE LA SURFACE
C
C ENTREES:
C --------
C NUVOLU : NUMERO DU VOLUME DANS LE LEXIQUE DES VOLUMES
C KNMVOL : NOM DU VOLUME
C HEXAED : MIN ET MAX DES 3 COORDONNEES D'UN HEXAEDRE ENGLOBANT
C          LES POINTS DU VOLUME PARTITIONNE A MAILLER
C NBVOPA : NOMBRE DE  VOLUMES V8 DU VOLUME PARTITIONNE V9
C NUVOPA : NUMERO DES VOLUMES V8 DU VOLUME PARTITIONNE V9
C MXSOMM : NOMBRE MAXIMAL DE SOMMETS UTILISABLES DANS PTXYZD
C MXARBO : NOMBRE MAXIMAL D'OCTAEDRES DANS LARBRO
C MXARBT : NOMBRE MAXIMAL DE TETRAEDRES DANS LARBRT
C MXARET : NOMBRE MAXIMAL D'ARETES DES OT
C MXFACO : NOMBRE MAXIMAL DECLARABLE DE FACES DU CONTOUR
C MXAREF : NOMBRE MAXIMAL D'ARETES DES FACES LEFACO
C
C ENTREES ET SORTIES:
C -------------------
C NBSOMM : NOMBRE ACTUEL DE POINTS DECLARES DANS PTXYZD
C PTXYZD : PAR POINT : X  Y  Z  DISTANCE_SOUHAITEE
C NPSOFR :  (1 000 000) * (NO SURFACE + 1) + NUMERO DU SOMMET DANS
C          LA NUMEROTATION DES SOMMETS DE LA SURFACE
C          -4 SI SOMMET DE LA GRILLE DES TETRAEDRES
C          -1 SI SOMMET DE LA GRILLE DES TETRAEDRES ELIMINE
C LARBRO : ARBRE-14 DES OCTAEDRES ( FOND DE LA TETRAEDRISATION )
C      LARBRO(0,0) : NO DU 1-ER OCTAEDRE VIDE DANS LARBRO
C      LARBRO(1,0) : MAXIMUM DU 1-ER INDICE DE LARBRO (ICI -1:20)
C      LARBRO(2,0) : MAXIMUM DECLARE DU 2-EME INDICE DE LARBRO
C                    (ICI = MXARBO)
C
C      LARBRO(-1:20,1) : RACINE DE L'ARBRE (OCTAEDRE SANS PERE)
C
C      LARBRO(-1,J) : NO DU PERE DE L'OCTAEDRE J DANS UN DES 2 ARBRES
C                     >0 => DANS LARBRO
C                     <0 => DANS LARBRT
C      LARBRO(0,J)  : 1 A 14 NO DE FILS DE L'OCTAEDRE J POUR SON PERE
C                     + 100 * NO TYPE DE L'OT J
C                     NO TYPE DE L'OT : 0 SI OCTAEDRE
C                                       1 SI TETRAEDRE T RONDE (T1)
C                                       2 SI TETRAEDRE T       (T2)
C   SI LARBRO(0,J)>0 ALORS J EST UN OCTAEDRE OCCUPE
C      SI LARBRO(1,.)>0 ALORS
C         LARBRO(1:14,J): NO (>0) LARBRO DES 14 SOUS-OCTA-TETRAEDRES
C      SINON
C         LARBRO(1:14,J):-NO PTXYZD DES 0 A 14 POINTS INTERNES DE L'OCTA J
C                         0  SI PAS DE POINT
C                        ( J EST ALORS UNE FEUILLE DE L'ARBRE )
C
C      LARBRO(15:20,J) : NO PTXYZD DES 6 SOMMETS DE L'OCTAEDRE J
C   SINON
C      LARBRO(0,J): -ADRESSE DANS LARBRO DE L'OCTAEDRE VIDE SUIVANT
C
C LARBRT : ARBRE-5 DES TETRAEDRES ( FOND DE LA TETRAEDRISATION )
C      LARBRT(0,0) : NO DU 1-ER TETRAEDRE VIDE DANS LARBRT
C      LARBRT(1,0) : MAXIMUM DU 1-ER INDICE DE LARBRT (ICI -1:9)
C      LARBRT(2,0) : MAXIMUM DECLARE DU 2-EME INDICE DE LARBRT
C                     (ICI = MXARBT)
C
C      LARBRT(-1,J) : NO DU PERE DU TETRAEDRE J DANS UN DES 2 ARBRES
C                     >0 => DANS LARBRO
C                     <0 => DANS LARBRT
C      LARBRT(0,J) : 0 A 4 NO DE FILS DU TETRAEDRE J POUR SON PERE
C                    + 100 * NO TYPE DE L'OT J
C                     NO TYPE DE L'OT : 0 SI OCTAEDRE
C                                       1 SI TETRAEDRE T RONDE (T1)
C                                       2 SI TETRAEDRE T       (T2)
C
C   SI LARBRT(0,J)>0 ALORS J EST UN TETRAEDRE OCCUPE
C      SI LARBRT(1,J)>0 ALORS
C         LARBRT(1:5,J): NO (>0) LARBRT DES 5 SOUS-OCTA-TETRAEDRES
C      SINON
C         LARBRT(1:5,J):-NO PTXYZD DES 0 A 5 POINTS INTERNES AU TETRA J
C                         0  SI PAS DE POINT
C                        ( J EST ALORS UNE FEUILLE DU ARBRE )
C
C      LARBRT(6:9,J) : NO PTXYZD DES 4 SOMMETS DU TETRAEDRE J
C   SINON
C      LARBRT(0,J): ADRESSE DANS LARBRT DU TETRAEDRE VIDE SUIVANT
C
C NUOTPT : NUMERO D'OT (>0 SI LARBRO, <0 SI LARBRT) DE CHAQUE POINT PTXYZD
C          SI SOMMET D'OT, C'EST LE NUMERO D'OT DU PERE
C          SI POINT INTERNE OU FRONTALIER, C'EST LE NUMERO D'OT LE CONTENANT
C
C LAVIDE : PREMIERE ADRESSE LARETE A EXPLORER POUR TROUVER UNE ARETE VIDE
C LARETE : TABLEAU DE HACHAGE DES ARETES DES OT
C          LARETE(1,J) = SOMMET 1 DANS PTXYZD
C          LARETE(2,J) = SOMMET 2 DANS PTXYZD
C          LARETE(3,J) = NUMERO DANS LARETE DE L'ARETE SUIVANTE
C          LARETE(4,J) = MILIEU DE L'ARETE DANS PTXYZD
C          HACHAGE(NS1,NS2) = (NS1 + NS2) MODULO MXARET + 1
C
C NBFACO : NOMBRE DE FACES TRIANGULAIRES RECENSEES SUR LES SURFACES FERMEES
C          ET LES INTERFACES ENTRE VOLUMES
C LEFACO : FACE DU CONTOUR OU INTERFACES ENTRE VOLUMES
C          IL CONTIENT DANS CET ORDRE
C          1:   =0 POUR UNE FACE VIDE
C          123: NO (DANS PTXYZD) DU SOMMET 1 < SOMMET 2 < SOMMET 3
C          45:  NO (DANS NUVOPA 0 SINON) DU VOLUME1, VOLUME2 DE LA FACE
C          678: NO (DANS LEFACO) DE LA FACE ADJACENTE PAR L'ARETE 1 2 3

C          9: CHAINAGE DU HACHAGE DES FACES ISSUES DE LA FONCTION DE HACHAGE
C             LEFACO(9,*) -> FACE SUIVANTE (*=0:VIDE, *<>0:NON VIDE)

C          10: NUMERO DE LA PREMIERE FACE ISSUE DE LA FONCTION DE HACHAGE
C              LA SOMME DES 3 SOMMETS MODULO MXFACO
C              LF = MOD( NOSOTR(1)+NOSOTR(2)+NOSOTR(3), MXFACO ) + 1
C              NF = LEFACO( 10, LF ) LE NUMERO DE LA 1-ERE FACE DANS LEFACO
C              SI LA FACE NE CONVIENT PAS. PASSAGE A LA SUIVANTE
C              NF = LEFACO( 9, NF )  ...

C          11: >0  NO NOTETR D'UN TETRAEDRE AYANT CETTE FACE,
C              =0  SINON
CCCC       LEFACO(12,*) NO DE FACE OpenCascade

C N1FASC : N1FASC(NS)=NUMERO (DANS LEFACO) D'UNE FACE DE SOMMET NS
C
C LARETF : 12: NS1 NS2 NO DES 2 SOMMETS EXTREMITES DE L'ARETE
C           3: ARETE_SUIVANTE
C           4: PREMIERE_FACE LEFACO AYANT CETTE ARETE PAR
C              HACHAGE(NS1,NS2)=(NS1 + NS2) MODULO MXAREF + 1
C
C SORTIES:
C --------
C QTRMIN : QUALITE MINIMALE DES TRIANGLES DES SURFACES FERMEES
C QTRMOY : QUALITE MOYENNE  DES TRIANGLES DES SURFACES FERMEES
C AREMIN : LONGUEUR DE LA PLUS PETITE ARETE DE LA TRIANGULATION INITIALE
C AREMAX : LONGUEUR DE LA PLUS GRANDE ARETE DE LA TRIANGULATION INITIALE
C AREMOY : LONGUEUR DE LA     MOYENNE ARETE DE LA TRIANGULATION INITIALE
C COSMIN : COSINUS MINIMAL DES ANGLES DES TRIANGLES
C COSMAX : COSINUS MAXIMAL DES ANGLES DES TRIANGLES LIMITANT LES ANGLES
C          CORRECTS D'UN TRIANGLE
C IERR   : 0 SI PAS D'ERREUR
C          1 SURFACE DE TYPE TRIANGLE STRUCTURE ou EXISTENCE d'UN QUADRANGLE
C          2 SATURATION DU TABLEAU des SOMMETS
C          3 SI TRIANGULATION AVEC UNE ARETE DE LONGUEUR NULLE
C          4 SI AU MOINS UNE ARETE APPARTIENT A UNE SEULE TRIANGLE
C               ou a AU MOINS 3 TRIANGLES
C          5 SATURATION du TABLEAU LEFACO
C          6 SATURATION du TABLEAU LARETF
C          7 AU MOINS UN TRIANGLE REDUIT A 1 ARETE
C          8 AU MOINS UN COUPLE DE TRIANGLES COLLES
C          MULTIPLE DE 10   SI TRIANGLE DEGENERE EN ARETE
C          MULTIPLE DE 100  SI TRIANGLE APPARTIENT A 3 VOLUMES
C          MULTIPLE DE 1000 SI ARETE APPARTENANT A UNE SEULE TRIANGLE
C                                    SURFACE OUVERTE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 1992
C MODIFS : ALAIN PERRONNET Saint PIERRE DU PERRAY               Mai 2020
C2345X7..............................................................012
      PARAMETER  (MXTR3V=1024,   MXAR1F=1024, MXAR3F=1024, MXAR2SF=1024,
     %            MX2TRCOL=1024, MXTR1AR=1024 )
      include"./incl/langue.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/gsmenu.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/a_surface__definition.inc"
      include"./incl/a_volume__definition.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
      COMMON / TRTETR / STOPTE,TRACTE
      LOGICAL           STOPTE,TRACTE,TRACTE0

      INTRINSIC         DBLE
      CHARACTER*24      KNMVOL, KNMVOLU, KNMSURF
      REAL              HEXAED(6,2)
      DOUBLE PRECISION  PTXYZD(4,MXSOMM)
      INTEGER           NUVOPA(1:NBVOPA),
     &                  LEFACO(1:11,0:MXFACO),
     &                  N1FASC(1:MXSOMM),
     &                  LARBRO(-1:20,0:MXARBO),
     &                  LARBRT(-1:9,0:MXARBT),
     &                  NUOTPT(1:MXSOMM),
     &                  NPSOFR(1:MXSOMM),
     &                  LARETE(1:4,1:MXARET),
     &                  LARETF(1:4,1:MXAREF)
C
      REAL              XYZ(1:3,1:3), DIST(3), VECT(3,2), COSMX, COSMN
      DOUBLE PRECISION  RADEGR, PT(4), SNF, SNFJAD, SURTRD, QD,
     %                  RASF2T, RASFMI, COS2TRD, VECNOR(3)
      INTEGER           NOAR1F(MXAR1F), NOAR3F(MXAR3F),NOAR2SF(MXAR2SF),
     %                  NOTR3V(MXTR3V), NO2TRCOL(2,MX2TRCOL),
     %                  NOTR1AR(MXTR1AR)
      INTEGER           NOSOTR(1:3), NOST(2)
      EQUIVALENCE      (NOST(1),NS1),(NOST(2),NS2)

      IERR    = 0
      TRACTE0 = TRACTE
      NBSOMM0 = NBSOMM
      PRINT*
      IF( LANGAG .EQ. 0 ) THEN
         PRINT*,'otajsf: AJOUT des TRIANGLES du VOLUME',NUVOLU,': ',
     %          KNMVOL,' aux',NBSOMM,' SOMMETS de l''OCTAEDRE ENGLOBANT'
      ELSE
         PRINT*,'otajsf: The TRIANGLES of VOLUME ',KNMVOL,
     %   ' are ADDED to',NBSOMM,' GLOBAL OCTAHEDRON VERTICES'
      ENDIF

C     PRESENCE ou NON de la FONCTION UTILISATEUR TAILLE_IDEALE(x,y,z)
      NOFOTI = NOFOTIEL()

C     DECLARATION INITIALE DU TABLEAU OTFL
      MNOTFL = 0
      MXOTFL = 128
      CALL TNMCDC( 'ENTIER', MXOTFL, MNOTFL )

C     INITIALISATION DE LEFACO
C     ------------------------
C     MISE A ZERO DU PREMIER SOMMET DES FACES DU CONTOUR
C     VALEUR TEST DE NON UTILISATION DE LA FACE
C     LEFACO(9,*) EST ACTUELLEMENT LE CHAINAGE DES FACES VIDES
C     PUIS IL VA DEVENIR CELUI DES FACES DE MEME VALEUR DE HACHAGE
      NBFACO = 0
      NBFAISOL= 0
      DO N=0,MXFACO
         DO K=1,11
            LEFACO( K, N ) = 0
         ENDDO
         LEFACO( 9, N ) = N+1
      ENDDO
C     FIN DU CHAINAGE DES FACES VIDES DANS LEFACO
      LEFACO( 9, MXFACO ) = 0
C     LA PREMIERE FACE VIDE DANS LEFACO
      LEFACO( 10,  0    ) = 1

C     CONVERSION RADIANS DEGRES
      RADEGR = 45D0 / ATAN( 1D0 )

C     LE COSINUS MINIMAL ET MAXIMAL LIMITANT LES ANGLES
C     CORRECTS d'UN TRIANGLE
      COSMX = REAL( COS(   5D0 / RADEGR ) )
      COSMN = REAL( COS( 170D0 / RADEGR ) )

      COSMAX =-2.
      COSMIN = 2.
      NBANMI = 0
      NBANMA = 0
      AREMOY = 0.
      AREMAX = 0.
      AREMIN = RINFO( 'GRAND' )
      QTRMIN = 2
      QTRMOY = 0
      NBAR2SF = 0
C     ====================================================
C     BOUCLE SUR LES VOLUMES 8 DE LA PARTITION DU VOLUME 9
C     ====================================================
      PRINT*
      NBTR3V = 0
      NBSTID = 0
      DO 500 N=1,NBVOPA

C        LE NUMERO DU VOLUME N
         NOVO = NUVOPA( N )
         CALL NMOBNU( 'VOLUME', NOVO, KNMVOLU )
         PRINT*
         IF( LANGAG .EQ. 0 ) THEN
            PRINT*,'otajsf: les TRIANGLES du VOLUME 8: ',KNMVOLU,
     %             ' sont AJOUTES'
         ELSE
            PRINT*,'otajsf: the TRIANGLES of VOLUME 8: ',KNMVOLU,
     %             ' are ADDED'
         ENDIF
C        LE TABLEAU 'DEFINITION' DE CE VOLUME
         CALL LXNLOU( NTVOLU, NOVO, NTLXVO, MNDFVO )
         CALL LXTSOU( NTLXVO, 'DEFINITION', NTDFVO, MNDFVO )
         NBSF1V = MCN( MNDFVO + WBSF1V )
         MNSV   = MNDFVO + WUSF1V - 1

C        BOUCLE SUR LES SURFACES FERMEES DU N-EME VOLUME 8
C        =================================================
         DO 490 IS = 1, NBSF1V

C           LE NUMERO DE LA SURFACE FERMEE DU CONTOUR DU VOLUME N
            NOSU = MCN( MNSV + IS )
            CALL NMOBNU( 'SURFACE', NOSU, KNMSURF )
C           LE LEXIQUE DE LA SURFACE
            CALL LXNLOU( NTSURF, NOSU, NTLXSU, MNFASU )

C           LE TABLEAU 'NSEF' DE CETTE SURFACE NOSU
            CALL LXTSOU( NTLXSU, 'NSEF', NTFASU, MNFASU )
C           LE TYPE STRUCTURE OU NON DE LA SURFACE
            NUTYMA = MCN( MNFASU + WUTYMA )
C           LE NOMBRE DE FACES DE LA SURFACE
            IF( NUTYMA .EQ. 0 ) THEN
C              SURFACE NON STRUCTUREE
C              NOMBRE DE TRIANGLES DE LA SURFACE NOSU
               NBFASU = MCN( MNFASU + WBEFOB )
C              NOMBRE DE SOMMETS STOCKES PAR TRIANGLE (Le 4-EME est 0)
               NBSOEF = MCN( MNFASU + WBSOEF )
            ELSE IF( NUTYMA .EQ. 3 ) THEN
C              TRIANGLE STRUCTURE
               NBLGRC(NRERR) = 1
C              POUR LEVER CETTE RESTRICTION SI BESOIN EST
C              UTILISER ~/util/NSEFPA.f et NSEFNS.f
C              EN FAIT UN TRIANGLE STRUCTURE NE PEUT ETRE UNE SURFACE FERMEE
               IF( LANGAG .EQ. 0 ) THEN
                  KERR(1)='otajsf: TRIANGULATION STRUCTUREE INTERDITE'
               ELSE
               KERR(1)='otajsf: a STRUCTURED TRIANGULATION is FORBIDDEN'
               ENDIF
               CALL LEREUR
               IERR = 1
               GOTO 9999
            ELSE
               IERR = 1
               GOTO 9999
            ENDIF

C           LE TABLEAU 'XYZSOMMET' DE CETTE SURFACE
            CALL LXTSOU( NTLXSU, 'XYZSOMMET', NTSOSU, MNSOSU )
            NBSOSU = MCN( MNSOSU + WNBSOM )

C           DECLARATION DU TABLEAU D'IDENTIFICATION DES SOMMETS DANS PTXYZD
            MXIDST = 1 + NBSOSU
            MNIDST = 0
            CALL TNMCDC( 'ENTIER', MXIDST, MNIDST )
            CALL AZEROI( MXIDST, MCN(MNIDST) )
            IF( LANGAG .EQ. 0 ) THEN
               PRINT*,'otajsf: SURFACE FERMEE',IS,
     %         ' Avant IDENTIFICATION: Nombre SOMMETS INITIAUX=',NBSOSU,
     %         ' Nombre TRIANGLES INITIAUX=',NBFASU
            ELSE
               PRINT*,'otajsf: CLOSED SURFACE',IS,
     %         ' BEFORE IDENTIFICATION: INITIAL VERTEX NUMBER=',NBSOSU,
     %         ' INITIAL TRIANGLE Number=',NBFASU
            ENDIF

C           IDENTIFICATION DES SOMMETS
C           --------------------------
            DO 100 NS1 = 1, NBSOSU

C              LES 3 COORDONNEES DE CE SOMMET DEBUTENT EN
               MNS = MNSOSU + WYZSOM + 3 * NS1 - 3

C              LES 3 COORDONNEES DU SOMMET NS1 INITIAL
C              RECHERCHE DU PLUS PROCHE POINT PTXYZD DE NS1 DANS LES OT
               PT(1) = DBLE( RMCN( MNS ) )
               PT(2) = DBLE( RMCN( MNS + 1 ) )
               PT(3) = DBLE( RMCN( MNS + 2 ) )
               CALL OTPLPRPT( PT, PTXYZD, LARBRO, LARBRT, MXOTFL,MNOTFL,
     %                        NS, DISMIN )

               IF( NS .GT. 0 .AND. DISMIN .EQ. 0 ) THEN
C                 NS1 A ETE IDENTIFIE AU POINT NS DE PTXYZD
                  GOTO 80
               ENDIF

               IF( NS .EQ. 0 ) THEN
                  PRINT*
                  IF( LANGAG .EQ. 0 ) THEN
                     PRINT*,'otajsf: PROBLEME: le SOMMET(',NS1,')=',
     %                     (PT(K),K=1,3),' N''EST PAS DANS l''OT RACINE'
                  ELSE
                     PRINT*,'otajsf: PROBLEM: the VERTEX(',NS1,')=',
     %                     (PT(K),K=1,3),' is NOT in the OT TREE ROOT'
                  ENDIF
                  GOTO 100
               ENDIF

C              NS1 EST UN NOUVEAU SOMMET A AJOUTER A PTXYZD
               IF( NBSOMM .GE. MXSOMM ) THEN
                  PRINT*
                  IF( LANGAG .EQ. 0 ) THEN
                     PRINT*,'otajsf: AUGMENTER LE NOMBRE DE SOMMETS MXSO
     %MM=',MXSOMM
                  ELSE
                     PRINT*,'otajsf: AUGMENT the VERTICES NUMBER MXSOMM=
     %',MXSOMM
                  ENDIF
                  IERR = 2
                  GOTO 9999
               ENDIF
               NBSTID = NBSTID + 1
               NBSOMM = NBSOMM + 1
               NS     = NBSOMM
C              LE SOMMET NS EST AJOUTE CAR NON RETROUVE
               PTXYZD(1,NS) = PT(1)
               PTXYZD(2,NS) = PT(2)
               PTXYZD(3,NS) = PT(3)
C
C              AJOUTER LE POINT NS DANS SON OT FEUILLE
C              EN AJOUTANT EVENTUELLEMENT DES SOMMETS D'OT (NBSOMM CHANGE)
               CALL OTRGPT( HEXAED, NBSOMM, MXSOMM, PTXYZD,
     &                      MXARBO, LARBRO, MXARBT, LARBRT, NUOTPT,
     &                      LAVIDE, MXARET, LARETE,
     &                      NUOTFE, IERR  )
               IF( IERR .NE. 0 ) THEN
                  GOTO 9999
               ENDIF

C              NS LE NUMERO DANS PTXYZD DU SOMMET INITIAL NS1
  80           MCN( MNIDST + NS1 ) = NS

C              LE NUMERO DE LA SURFACE ET DU SOMMET LOCAL
               NPSOFR( NS ) = 1 000 000 * (NOSU+1) + NS1

 100        ENDDO

            PRINT*,'otajsf:',NBSOMM,' SOMMETS des SURFACES et ARBRES OT'
            PRINT*


C           BOUCLE SUR LES TRIANGLES OU FACES DE LA SURFACE FERMEE IS
C           =========================================================
            MNSOEL = MNFASU + WUSOEF - 1 - NBSOEF
            DO 400 NOFA = 1, NBFASU

C              LE NUMERO DES SOMMETS DE LA FACE NOFA
               MNSOEL = MNSOEL + NBSOEF
C
               IF( MCN(MNSOEL+4) .GT. 0 ) THEN
C                 QUADRANGLE INTERDIT!
                  NBLGRC(NRERR) = 2
                  IF( LANGAG .EQ. 0 ) THEN
                     KERR(1)='La PEAU du VOLUME DOIT ETRE TRIANGULEE'
         KERR(2)='L''OPTION 31 SUBDIVISE les QUADRANGLES en 2 TRIANGLES'
                  ELSE
                     KERR(1)='The BOUNDARY MUST BE TRIANGULATED'
         KERR(2)='the OPTION 31 SUBDIVISES QUADRANGLES into 2 TRIANGLES'
                  ENDIF
                  CALL LEREUR
                  IERR = 1
                  GOTO 9999
               ENDIF

C              XYZ DES 3 SOMMETS DE LA FACE NOFA
               DO J=1,3
                  NS1 = MCN( MNSOEL + J )
                  NOSOTR( J ) = NS1
C                 LES 3 COORDONNEES DE CE SOMMET NS1
                  MNS = MNSOSU + WYZSOM + 3 * NS1 - 3
C                 LES COORDONNEES DES 3 SOMMETS
                  XYZ(1,J) = RMCN( MNS )
                  XYZ(2,J) = RMCN( MNS + 1 )
                  XYZ(3,J) = RMCN( MNS + 2 )
               ENDDO

C              QUALITE DU TRIANGLE NOFA
               CALL QUALEF( 3, NOSOTR, MCN(MNSOSU+WNBSOM),
     %                      RMCN(MNSOSU+WYZSOM), SURFVOLU, Q, IER )
               QTRMIN = MIN( QTRMIN, Q )
               QTRMOY = QTRMOY + Q

C              LA LONGUEUR DES 3 COTES DE LA FACE
               DO J=1,3
                  IF( J .EQ. 3 ) THEN
                     K=1
                  ELSE
                     K = J+1
                  ENDIF
                  DIST(J) = DIST2P( XYZ(1,J), XYZ(1,K) )
                  IF( DIST(J) .GT. AREMAX ) AREMAX = DIST(J)
                  IF( DIST(J) .LT. AREMIN ) AREMIN = DIST(J)
               ENDDO
               AREMOY = AREMOY + ( DIST(1)+DIST(2)+DIST(3) ) / 3
C
C              EXISTE-T-IL UN ANGLE > 180-10=170 DEGRES ?
               COSMA = -2.
               COSMI =  2.
               DO J=1,3
                  IF( J .EQ. 1 ) THEN
                     J0 = 3
                  ELSE
                     J0 = J - 1
                  ENDIF
                  IF( J .EQ. 3 ) THEN
                     J1 = 1
                  ELSE
                     J1 = J + 1
                  ENDIF
C                 LE COSINUS DE L'ANGLE AU SOMMET J
                  DO KL=1,3
                     VECT(KL,1) = XYZ(KL,J0) - XYZ(KL,J)
                     VECT(KL,2) = XYZ(KL,J1) - XYZ(KL,J)
                  ENDDO
                  COSM = PROSCR( VECT(1,1), VECT(1,2), 3 )
                  COSM = COSM / ( DIST(J0) * DIST(J) )
                  IF( COSM .LT. COSMI ) THEN
                     COSMI = COSM
                  ENDIF
                  IF( COSM .GT. COSMA ) THEN
                     COSMA = COSM
                  ENDIF
               ENDDO

               IF( COSMI .LT. COSMN ) THEN
C                 => L'ANGLE EST PLUS GRAND QUE 170. DEGRES
                  NBANMA = NBANMA + 1
               ELSE IF( COSMA .GT. COSMX ) THEN
C                 => L'ANGLE EST PLUS PETIT QUE 5 DEGRES
                  NBANMI = NBANMI + 1
               ENDIF

               IF( Q .LT. 0.1 .OR.
     %            (COSMA .GT. COSMX .OR. COSMI .LT. COSMN) ) THEN

C                 TRAITEMENT DES ERREURS D'ARRONDIS
                  IF( COSMA .GT. 1. ) THEN
                     COSMA = 1.
                  ENDIF
                  IF( COSMI .LT. -1. ) THEN
                     COSMI = -1.
                  ENDIF

                  IF( LANGAG .EQ. 0 ) THEN
                     PRINT 10309, KNMVOLU, KNMSURF, NOFA, Q,
     %               ACOS(COSMI)*RADEGR, ACOS(COSMA)*RADEGR,
     %              (NOSOTR(L),(XYZ(K,L),K=1,3),L=1,3)
                  ELSE
                     PRINT 20309, KNMVOLU, KNMSURF, NOFA, Q,
     %               ACOS(COSMI)*RADEGR, ACOS(COSMA)*RADEGR,
     %               (NOSOTR(L),(XYZ(K,L),K=1,3),L=1,3)
                  ENDIF
               ENDIF

10309    FORMAT(' VOLUME ',A,' SURFACE FERMEE ',A,
     %    ' TRIANGLE',I9,' QUALITE=',F8.5,
     %    ' ANGLE MAX= ',F7.2,' ANGLE MIN= ',F7.2,' degres'/
     %    3(' St=',I9,' X=',G15.6,' Y=',G15.6,' Z=',G15.6/))
20309    FORMAT(' VOLUME ',A,' CLOSED SURFACE ',A,
     %    ' TRIANGLE',I9,' QUALITY=',F8.5,
     %    ' ANGLE MAX= ',F7.2,' ANGLE MIN= ',F7.2,' degrees'/
     %    3(' Vx=',I9,' X=',G15.6,' Y=',G15.6,' Z=',G15.6/))

               COSMAX = MAX( COSMA, COSMAX )
               COSMIN = MIN( COSMI, COSMIN )

ccc               DMAX   = MAX( DIST(1), DIST(2), DIST(3) )
ccc               DMOY   = ( DIST(1) + DIST(2) + DIST(3) ) / 3

C              TRAITEMENT DU J-EME SOMMET DU TRIANGLE DE LA SURFACE
C              ====================================================
               DO J=1,3

C                 NS1 NUMERO DU SOMMET J DE LA FACE FRONTIERE NOFA
                  NS1 = NOSOTR(J)

C                 LE NUMERO DANS LA TETRAEDRISATION DU SOMMET NS1
                  NS = MCN( MNIDST + NS1 )
                  NOSOTR( J ) = NS

CCCC                 LA DISTANCE SOUHAITEE EST LA PLUS PETITE DES
CCCC                 PLUS GRANDES ARETES D'UNE FACE DU SOMMET NS
CCC                  IF( PTXYZD(4,NS) .LE. 0D0 ) THEN
CCC                     PTXYZD(4,NS) = DMAX
CCC                  ELSE
CCC                     PTXYZD(4,NS) = MIN( PTXYZD(4,NS), DBLE(DMAX) )
CCC                  ENDIF
cccC                 LA DISTANCE SOUHAITEE EST LA PLUS
cccC                 GRANDE DES ARETES D'UNE FACE DU SOMMET NS
ccc                  PTXYZD(4,NS) = MAX( PTXYZD(4,NS), DBLE(DMAX) )
cccC                 LA DISTANCE SOUHAITEE EST LA PLUS GRANDE
cccC                 ARETE MOYENNE DES FACES DE SOMMET NS
ccc                  PTXYZD(4,NS) = MAX( PTXYZD(4,NS), DBLE(DMOY) )

C                 LA DISTANCE SOUHAITEE EST AU POINT NS
                  IF( NOFOTI .NE. 0 ) THEN

C                    LA TAILLE_IDEALE( PTXYZD(NS) )
                     CALL TAILIDEA( NOFOTI, PTXYZD(1,NS),
     %                              NCODEV, PTXYZD(4,NS) )

                  ELSE

C                    LA DISTANCE SOUHAITEE EST LE MINIMUM DE LA
C                    MOYENNE DES LONGUEURS D'ARETES ISSUES DE CE SOMMET
                     IF( J .EQ. 1 ) THEN
                        J0 = 3
                     ELSE
                        J0 = J-1
                     ENDIF
                     DMOY = ( DIST(J) + DIST(J0) ) / 2

                     IF( PTXYZD(4,NS) .LE. 0D0 ) THEN
                        PTXYZD(4,NS) = DBLE( DMOY )
                     ELSE
                        PTXYZD(4,NS) = MIN( PTXYZD(4,NS), DBLE(DMOY) )
                     ENDIF

                  ENDIF

               ENDDO

C              TRAITEMENT DE LA FACE NOFA
C              ==========================
               IF( NOSOTR(1) .EQ. NOSOTR(2) .OR.
     %             NOSOTR(2) .EQ. NOSOTR(3) .OR.
     %             NOSOTR(3) .EQ. NOSOTR(1) ) THEN
                   IF( LANGAG .EQ. 0 ) THEN
                      PRINT*,'otajsf: Le TRIANGLE',NOFA,
     %                     ' de la SURFACE',KNMSURF,
     %                     'a 2 SOMMETS IDENTIQUES apres IDENTIFICATION'
                   ELSE
                      PRINT*,'otajsf: Le TRIANGLE',NOFA,
     %                      ' of SURFACE',KNMSURF,
     %                   'HAS 2 IDENTICAL VERTICES after IDENTIFICATION'
                   ENDIF
                   PRINT*,'otajsf: TRIANGLE',NOFA,' :',NOSOTR
                   IERR = IERR + 10
                   GOTO 400
               ENDIF

cccC              TRI CROISSANT DES 3 NUMEROS DES SOMMETS DE LA FACE
ccc               CALL TRI3NO( NOSOTR, NOSOTR ) FAIT dans nuletr

C              LA FACE NOSOTR EXISTE-T-ELLE DEJA DANS LEFACO?
C              HACHAGE AVEC LA SOMME DES 3 SOMMETS MODULO MXFACO
               LF = MOD( NOSOTR(1) + NOSOTR(2) + NOSOTR(3), MXFACO ) + 1
               CALL NULETR( NOSOTR, MXFACO, LEFACO,  NF )
C              NF>0 NUMERO LEFACO DE LA FACE NOSOTR SI ELLE EST RETROUVEE
C                =0 SI LA FACE N'EST PAS RETROUVEE
               IF( NF .GT. 0 ) THEN

C                 LA FACE NOSOTR EST RETROUVEE DANS LEFACO
C                 ----------------------------------------
C                 LE SECOND NUMERO DE VOLUME EST MIS A JOUR
                  IF( LEFACO(5,NF) .GT. 0 ) THEN

C                    LE TRIANGLE NF APPARTIENT A 3 VOLUMES
                     NBLGRC(NRERR) = 2
                     WRITE(KERR(MXLGER)(1:10),'(I10)') NF
                     IF( LANGAG .EQ. 0 ) THEN
                        PRINT*,'otajsf: Le TRIANGLE',NF,
     %                         'APPARTIENT A 3 VOLUMES => INTERDIT'
                     ELSE
                        PRINT*,'otajsf: The TRIANGLE',NF,
     %                         'BELONGS to 3 VOLUMES => FORBIDDEN'
                     ENDIF

                     IF( NBTR3V .LT. MXTR3V ) THEN
                        NBTR3V = NBTR3V + 1
                        NOTR3V( NBTR3V ) = NF
                     ENDIF

                     IERR = IERR + 100
                     GOTO 400

                  ENDIF

C                 LE NUMERO DE 1 A NBVOPA DU 2-EME VOLUME DE CETTE FACE NF
                  LEFACO(5,NF) = N

               ELSE

C                 NON: LA FACE NOSOTR N'A PAS ETE RETROUVEE DANS LEFACO
C                      ELLE Y EST AJOUTEE
C                 RECHERCHE D'UNE FACE VIDE DANS LEFACO
                  IF( LEFACO( 10, 0 ) .LE. 0 ) THEN
                     NBLGRC(NRERR) = 1
                     IF( LANGAG .EQ. 0 ) THEN
                        KERR(1) ='otajsf: SATURATION des FACES LEFACO'
                     ELSE
                        KERR(1) ='otajsf: ARRAY LEFACO SATURATED'
                     ENDIF
                     CALL LEREUR
                     IERR = 5
                     GOTO 9999
                  ENDIF

C                 IL EXISTE UNE FACE LEFACO NF VIDE
                  NF = LEFACO( 10, 0 )
C                 MISE A JOUR DE LA PREMIERE FACE VIDE DE LEFACO
                  LEFACO( 10, 0 ) = LEFACO( 9, NF )

C                 LA FACE NF CREEE DANS LEFACO DEVIENT LA PREMIERE
C                 DU HACHAGE DES FACES
                  LEFACO(  9, NF ) = LEFACO( 10, LF )
                  LEFACO( 10, LF ) = NF

C                 LA FACE EST AJOUTEE DANS LE TABLEAU LEFACO
C                 LES 3 SOMMETS TRIES PAR VALEUR CROISSANTE
                  LEFACO(1,NF) = NOSOTR(1)
                  LEFACO(2,NF) = NOSOTR(2)
                  LEFACO(3,NF) = NOSOTR(3)

C                 LE NUMERO DU PREMIER VOLUME
                  LEFACO(4,NF) = N
                  LEFACO(5,NF) = 0

C                 MISE A JOUR DU TABLEAU NO D'UNE FACE AYANT CES SOMMETS
                  N1FASC( NOSOTR(1) ) = NF
                  N1FASC( NOSOTR(2) ) = NF
                  N1FASC( NOSOTR(3) ) = NF

C                 UNE FACE DES SURFACES ET INTERFACES DE PLUS
                  NBFACO = NBFACO + 1

               ENDIF

 400        ENDDO

            IF( MNIDST .GT. 0 ) CALL TNMCDS( 'ENTIER', MXIDST, MNIDST )
            IF( IERR   .NE. 0 ) GOTO 9999

C           FIN SURFACE IS
 490     ENDDO

C        FIN VOLUME V8 N
 500  ENDDO


C     MISE A JOUR DU CHAINAGE DES FACES ADJACENTES PAR LES ARETES
C     ===========================================================
C     NOMBRE DE TRIANGLES REDUITS A UNE ARETE
      NBTR1AR = 0
C     NOMBRE D'ARETES APPARTENANT A 3 FACES
      NBAR3F  = 0
C     INITIALISATION DU HACHAGE DES ARETES DES FACES LEFACO
      CALL AZEROI( 4*MXAREF, LARETF )
C     LA 1-ERE ARETE LIBRE DERRIERE CELLES ADRESSEES PAR LE MINIMUM
      LAVIDF = MXAREF

      DO 590 NF=1,MXFACO

         IF( LEFACO(1,NF) .LE. 0 ) GOTO 590

C        LA FACE EXISTE:  EST ELLE REDUITE A UNE ARETE?
         CALL NORFA3D( PTXYZD(1,LEFACO(1,NF)),
     %                 PTXYZD(1,LEFACO(2,NF)),
     %                 PTXYZD(1,LEFACO(3,NF)),
     %                 VECNOR, IER )

         IF( IER .NE. 0 ) THEN

C           OUI: ELLE EST MISE DANS NOTR1AR POUR ETRE SIGNALEE
            IF( NBTR1AR .GE. MXTR1AR ) GOTO 590
            NBTR1AR = NBTR1AR + 1
            NOTR1AR( NBTR1AR ) = NF
C           MAIS LE TRIANGLE NF DE LEFACO N'EST PAS SUPPRIME
C           POUR CONSERVER LA CONTINUITE DES EF A TRAVERS LES ARETES
            IERR = 7

         ENDIF

C        BOUCLE SUR SES 3 ARETES
         DO 530 J=1,3
            IF( J .EQ. 3 ) THEN
               NS1 = 1
               NS2 = 3
            ELSE
               NS1 = J
               NS2 = J + 1
            ENDIF
            NS1 = LEFACO( NS1, NF )
            NS2 = LEFACO( NS2, NF )
            IF( NS1 .GT. NS2 ) THEN
C              TRI CROISSANT DES 2 NO DE SOMMETS DE L'ARETE
               NOAR = NS1
               NS1  = NS2
               NS2  = NOAR
            ENDIF

C           HACHAGE DE L'ARETE NS1 NS2 POUR LA RETROUVER
            CALL HACHAG( 2, NOST, 4, MXAREF, LARETF, 3,  LAVIDF, NOAR )

            IF( NOAR .EQ. 0 )  THEN
C              SATURATION DES ARETES DE LARETF
               NBLGRC(NRERR) = 1
               IF( LANGAG .EQ. 0 ) THEN
                  KERR(1) ='SP otajsf: TABLEAU LARETF SATURE'
               ELSE
                  KERR(1) ='SP otajsf: ARRAY LARETF SATURATED'
               ENDIF
               CALL LEREUR
               IERR = 6
               GOTO 9999
            ENDIF

            IF( NOAR .LT. 0 ) THEN
C              L'ARETE EST CREEE:  AJOUT DU NUMERO DE FACE LEFACO
               LARETF(4,-NOAR) = NF
            ELSE
C              L'ARETE EXISTE : CHAINAGE CIRCULAIRE DES FACES LEFACO
C              LA PREMIERE FACE CHAINEE EST
               NF0 = LARETF(4,NOAR)
C              RECHERCHE DE L'ARETE I DE NF0  DE SOMMETS NS1 NS2
               DO I=1,3
                  IF( I .EQ. 3 ) THEN
                     II = 1
                  ELSE
                     II = I + 1
                  ENDIF
                  NSS1 = LEFACO(I,NF0)
                  NSS2 = LEFACO(II,NF0)
                  IF( ( NS1 .EQ. NSS1 .AND. NS2 .EQ. NSS2 )   .OR.
     %                ( NS2 .EQ. NSS1 .AND. NS1 .EQ. NSS2 ) ) GOTO 527
               ENDDO
               PRINT*
               PRINT*,'PB otajsf: LEFACO(',NF0,')=',
     %                           (LEFACO(K,NF0),K=1,11)
               PRINT*,'Chainage circulaire INCOMPLET des faces'
               GOTO 530

 527           NFA = LEFACO( 5+I, NF0 )
               IF( NFA .EQ. 0 ) THEN
C                 PREMIERE FACE ADJACENTE DE NF et NF0
                  LEFACO(5+I,NF0) = NF
                  LEFACO(5+J,NF ) = NF0
               ELSE
C                 DEUXIEME FACE ADJACENTE
                  LEFACO(5+I,NF0) = NF
                  LEFACO(5+J,NF ) = NFA
                  IF( NF0 .NE. NFA ) THEN
                     IF( NBAR3F .LT. MXAR3F ) THEN
                        NBAR3F = NBAR3F + 1
                        NOAR3F( NBAR3F ) = NF
                     ENDIF
                     PRINT*,'PB otajsf: 3 FACES ADJACENTES PAR L''ARETE'
     %                     ,NS1,NS2,' NF0=',NF0,' NF=',NF,' NFA=',NFA
                  ENDIF
               ENDIF
            ENDIF

 530     ENDDO
 590  ENDDO

      IF( NBAR3F .GT. 0 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
      PRINT*,'otajsf:',NBAR3F,' ARETES APPARTENANT a AU MOINS 3 TRIANGLE
     %S'
         ELSE
         PRINT*,'otajsf:',NBAR3F,' EDGES BELONGING at LEAST 3 TRIANGLES'
         ENDIF
         IERR = 4
      ENDIF

C     EXISTE T IL DES ARETES LEFACO APPARTENANT A UN SEUL TRIANGLE? et
C     CACUL DU RAPPORT DES SURFACES DES 2 TRIANGLES ADJACENTS PAR CETTE ARETE
C     =======================================================================
      NBAR1F = 0
C     NOMBRE DE COUPLES DE TRIANGLES COLLES (ANGLE DIEDRE TRES FAIBLE)
      NB2TRCOL = 0
C     NOMBRE DE RAPPORT DE SURFACE DE 2 FACES ADJACENTES
      NBAR2SF = 0
      RASFMI = 1D100
      DO NF = 1, MXFACO

         IF( LEFACO(1,NF) .GT. 0 ) THEN

C           LA FACE EXISTE: EST ELLE ISOLEE?
C           I.E. SES 3 ARETES N'ONT PAS D'AUTRE TRIANGLE ADJACENT
            DO J=1,3
C              LA FACE OPPOSEE A L'ARETE J DE NF
               NFJAD = LEFACO(5+J,NF)
               IF( NFJAD .GT. 0 ) GOTO 600
            ENDDO

C           LE TRIANGLE NF N'A PAS DE VOISIN: IL EST SUPPRIME
C           --------------------------------- ---------------
            PRINT*
            PRINT*,'otajsf: Le TRIANGLE',NF,
     %             ' EST ISOLE -> IL est SUPPRIME dans LEFACO'
            PRINT*,'otajsf: TRIANGLE(',NF,'): St=',(LEFACO(K,NF),K=1,3),
     %             ' V1=',LEFACO(4,NF),' V2=',LEFACO(5,NF),
     %             ' FACES Adj=',(LEFACO(K,NF),K=6,8),
     %             ' du tetraedre',LEFACO(11,NF)
            DO K = 1, 3
               NST = LEFACO(K,NF)
               PRINT*,'otajsf: St',NST,': x=',PTXYZD(1,NST),
     %                ' y=',PTXYZD(2,NST),' z=',PTXYZD(3,NST),
     %                ' d=',PTXYZD(4,NST)
            ENDDO
C           LE TRIANGLE NF DE LEFACO EST SUPPRIME
            CALL VDDSFA( NF, MXFACO, LEFACO, NBFACO, N1FASC )
            NBFAISOL = NBFAISOL + 1
            GOTO 620

C           NF N'EST PAS ISOLE. BOUCLE SUR SES 3 ARETES
C           ------------------- -----------------------
 600        DO 610 J=1,3

C              LA FACE OPPOSEE A L'ARETE J DU TRIANGLE NF
               NFJAD = LEFACO(5+J,NF)
               IF( NFJAD .EQ. 0 ) THEN

C                  UNE SEULE FACE DE LEFACO POUR CETTE ARETE J
C                  -> LA SURFACE EST OUVERTE
C                  -------------------------------------------
                   IF( J .EQ. 3 ) THEN
                      NS1 = 1
                      NS2 = 3
                   ELSE
                      NS1 = J
                      NS2 = J + 1
                   ENDIF

C                  LA SURFACE DU TRIANGLE NF
                   SNF = SURTRD( PTXYZD(1,LEFACO(1,NF)),
     %                           PTXYZD(1,LEFACO(2,NF)),
     %                           PTXYZD(1,LEFACO(3,NF)) )

C                  LA QUALITE DU TRIANGLE NF
                   CALL QUATRID( LEFACO(1,NF), PTXYZD, QD )

                   PRINT*
                   IF( LANGAG .EQ. 0 ) THEN
                      PRINT*,'otajsf: PB l''ARETE de St:',
     %                        LEFACO(NS1,NF),LEFACO(NS2,NF),
     %                       'APPARTIENT a UN SEUL TRIANGLE=',NF,
     %                      ' Qualite=',QD,' Surface=',SNF
                   ELSE
                      PRINT*,'otajsf: PB The EDGE of VERTICES',
     %                  LEFACO(NS1,NF),LEFACO(NS2,NF),
     %                 'BELONGS TO ONE and ONLY ONE TRIANGLE=',NF,
     %                      ' Quality=',QD,' Surface=',SNF
                   ENDIF

C                  AFFICHAGE DE LA FACE OUVERTE NF de LEFACO
                   PRINT*,'otajsf: TRIANGLE(',NF,
     %                    '): St=',(LEFACO(K,NF),K=1,3),
     %                    ' V1=',LEFACO(4,NF),' V2=',LEFACO(5,NF),
     %                    ' FACES Adj=',(LEFACO(K,NF),K=6,8),
     %                    ' du tetraedre',LEFACO(11,NF)
                   DO K = 1, 3
                      NST = LEFACO(K,NF)
                      PRINT*,'otajsf: St',NST,': x=',PTXYZD(1,NST),
     %                    ' y=',PTXYZD(2,NST), ' z=',PTXYZD(3,NST),
     %                    ' d=',PTXYZD(4,NST)
                   ENDDO

                   IF( NBAR1F .LT. MXAR1F ) THEN
                      NBAR1F = NBAR1F + 1
                      NOAR1F( NBAR1F ) = NF
                   ENDIF

                   IERR = IERR + 1000
                   GOTO 610

               ENDIF


C              COMPARAISON DES SURFACES des TRIANGLES NF et NFJAD
C              ADJACENTES PAR L'ARETE J DE NF
C              --------------------------------------------------
               SNF    = SURTRD( PTXYZD(1,LEFACO(1,NF)),
     %                          PTXYZD(1,LEFACO(2,NF)),
     %                          PTXYZD(1,LEFACO(3,NF)) )
               SNFJAD = SURTRD( PTXYZD(1,LEFACO(1,NFJAD)),
     %                          PTXYZD(1,LEFACO(2,NFJAD)),
     %                          PTXYZD(1,LEFACO(3,NFJAD)) )
               IF( SNF .LE. SNFJAD ) THEN
                  RASF2T = SNF / SNFJAD
               ELSE
                  RASF2T = SNFJAD / SNF
               ENDIF

C              LE PLUS PETIT RAPPORT DES SURFACES DE 2 FACES ADJACENTES
               RASFMI = MIN( RASFMI, RASF2T )

               IF( RASF2T .LT. 0.0667D0 ) THEN

C                 RAPPORT des 2 SURFACES TROP PETIT: A SIGNALER
                  IF( NBAR2SF .LT. MXAR2SF ) THEN
                     NBAR2SF = NBAR2SF + 1
                     NOAR2SF( NBAR2SF ) = NF
                  ENDIF

                  PRINT*
                  IF( LANGAG .EQ. 0 ) THEN
                     PRINT*,'otajsf: Le RAPPORT',RASF2T,
     %                 ' des SURFACES des TRIANGLES ADJACENTS',NF,NFJAD,
     %                 ' PEUT POSER PROBLEME dans la TETRAEDRISATION'
                  ELSE
                     PRINT*,'otajsf: The SURFACE RATIO',RASF2T,
     %                      ' of ADJACENT TRIANGLES',NF,NFJAD,
     %                      ' MAY BE a PROBLEM TO TETRAHEDRIZE'
                  ENDIF
                  CALL QUATRID( LEFACO(1,NF), PTXYZD, QD )
                  PRINT*,'otajsf: TRIANGLE',NF,' St:',
     %                   (LEFACO(K,NF),K=1,3),
     %                   ' QUALITE=',QD,' SURFACE=',SNF
                  CALL QUATRID( LEFACO(1,NFJAD), PTXYZD, QD )
                  PRINT*,'otajsf: TRIANGLE',NFJAD,' St:',
     %                   (LEFACO(K,NFJAD),K=1,3),
     %                   ' QUALITE=',QD,' SURFACE=',SNFJAD

               ENDIF

C              CALCUL DE L'ANGLE DIEDRE DES 2 TRIANGLES ADJACENTS NF-NFJAD
C              -----------------------------------------------------------
               GOTO( 601, 602, 603 ),J
 601           NS1 = 1
               NS2 = 2
               NS3 = 3
               GOTO 605

 602           NS1 = 2
               NS2 = 3
               NS3 = 1
               GOTO 605

 603           NS1 = 3
               NS2 = 1
               NS3 = 2

 605           NS1 = LEFACO( NS1, NF )
               NS2 = LEFACO( NS2, NF )
               NS3 = LEFACO( NS3, NF )

               DO K = 1, 3
                  NS4 = LEFACO( K, NFJAD )
                  IF( NS4.NE.NS1.AND.NS4.NE.NS2.AND.NS4.NE.NS3 )GOTO 607
               ENDDO
               GOTO 610

 607           CALL COS2TD( PTXYZD(1,NS1), PTXYZD(1,NS2), PTXYZD(1,NS3),
     %                      PTXYZD(1,NS2), PTXYZD(1,NS1), PTXYZD(1,NS4),
     %                      COS2TRD, IERR1, IERR2 )
               IF( IERR1 .NE. 0 .OR. IERR2 .NE. 0 ) GOTO 610

               IF( COS2TRD .LT. -0.99 ) THEN

C              COANPL: SEUIL DU COSINUS DE L'ANGLE FORME PAR LES NORMALES AUX
C              2 FACES ET AU DESSUS DUQUEL LES FACES SONT CONSIDEREES COPLANAIRES
C              ( 0.95     => 18.2  DEGRES )
C              ( 0.96     => 16.3  DEGRES )
C              ( 0.97     => 14.1  DEGRES )
C              ( 0.98     => 11.5  DEGRES )
C              ( 0.99     =>  8.11 DEGRES )
C              ( 0.9962   =>  5    DEGRES )
C              ( 0.99756  =>  4    DEGRES )
C              ( 0.99863  =>  3    DEGRES )
C              ( 0.999    =>  2.56 DEGRES )
C              ( 0.9999   =>  0.8  DEGRES )

C                 LE COUPLE DES TRIANGLES LEFACO NF et NFJAD SONT COLLES

                  IF( LANGAG .EQ. 0 ) THEN
                     PRINT*,'otajsf: les TRIANGLES',NF,NFJAD,
     %               ' COS(ANGLE des NORMALES)=',COS2TRD,' SONT COLLES'
                     PRINT*,'otajsf: TRIANGLE(',NF,
     %                      '): St=',(LEFACO(K,NF),K=1,3),
     %                      ' V1=',LEFACO(4,NF),' V2=',LEFACO(5,NF),
     %                      ' FACES Adj=',(LEFACO(K,NF),K=6,8),
     %                      ' du tetraedre',LEFACO(11,NF)
                     PRINT*,'otajsf: TRIANGLE(',NFJAD,
     %                    '): St=',(LEFACO(K,NFJAD),K=1,3),
     %                    ' V1=',LEFACO(4,NFJAD),' V2=',LEFACO(5,NFJAD),
     %                    ' FACES Adj=',(LEFACO(K,NFJAD),K=6,8),
     %                    ' du tetraedre',LEFACO(11,NFJAD)
                  ELSE
                     PRINT*,'otajsf: the TRIANGLES',NF,NFJAD,
     %               ' COS(NORMAL VECTOR ANGLE)=',COS2TRD,' ARE GLUED'
                     PRINT*,'otajsf: TRIANGLE(',NF,
     %                      '): Vx=',(LEFACO(K,NF),K=1,3),
     %                      ' V1=',LEFACO(4,NF),' V2=',LEFACO(5,NF),
     %                      ' FACES Adj=',(LEFACO(K,NF),K=6,8),
     %                      ' dof tetrahedron',LEFACO(11,NF)
                  ENDIF

                  IF( NB2TRCOL .GE. MX2TRCOL ) GOTO 610
                  NB2TRCOL = NB2TRCOL + 1
                  NO2TRCOL(1,NB2TRCOL) = NF
                  NO2TRCOL(2,NB2TRCOL) = NFJAD

               ENDIF

 610        ENDDO

         ENDIF

 620  ENDDO


      IF( NBAR1F .GT. 0 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
        PRINT*,'otajsf:',NBAR1F,' ARETES APPARTENANT a UN SEUL TRIANGLE'
           PRINT*,'otajsf: la SURFACE est OUVERTE et DOIT ETRE REFERMEE'
         ELSE
            PRINT*,'otajsf:',NBAR1F,' EDGES BELONGING at ONE TRIANGLE'
            PRINT*,'otajsf: the SURFACE is OPENED and HAD to be CLOSED'
         ENDIF
         IERR = 4
      ENDIF


      PRINT*
      IF( LANGAG .EQ. 0 ) THEN
         PRINT*,'otajsf: NOMBRE de RAPPORTS<0.0667 des SURFACES de 2 TRI
     %ANGLES ADJACENTS=',NBAR2SF
         PRINT*,'otajsf: le PLUS PETIT RAPPORT des SURFACES de 2 TRIANGL
     %ES ADJACENTS=',RASFMI
      ELSE
      PRINT*,'otajsf: 2 ADJACENT TRIANGLE SURFACES RATIO<0.0667 NUMBER='
     %      ,NBAR2SF
      PRINT*,'otajsf: 2 ADJACENT TRIANGLE SURFACES RATIO MINIMUM VALUE='
     %      ,RASFMI
      ENDIF


      IF( NB2TRCOL .GT. 0 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            PRINT*,'otajsf:',NB2TRCOL,' COUPLES DE TRIANGLES COLLES'
         ELSE
            PRINT*,'otajsf:',NB2TRCOL,' GLUED TRIANGLE COUPLES'
         ENDIF
         IERR = 8
      ENDIF


      IF( NBTR1AR .GT. 0 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            PRINT*,'otajsf:',NBTR1AR,' TRIANGLES REDUITS a UNE ARETE'
         ELSE
            PRINT*,'otajsf:',NBTR1AR,' REDUCED TRIANGLES to ONE ARETE'
         ENDIF
         IERR = 7
      ENDIF


      IF( NBAR1F.GT.0 .OR. NBAR3F  .GT.0 .OR. NBAR2SF.GT.0 .OR.
     %    NBTR3V.GT.0 .OR. NB2TRCOL.GT.0 .OR. NBTR1AR.GT.0 ) THEN

C        TRACE DES TRIANGLES LEFACO AVEC ANOMALIE
         TRACTE = .TRUE.
         CALL TRASFV9( MXFACO,   LEFACO,   NBSOMM,  PTXYZD, HEXAED,
     %                 NBTR3V,   NOTR3V,   NBAR2SF, NOAR2SF,
     %                 NBAR1F,   NOAR1F,   NBAR3F,  NOAR3F,
     %                 NB2TRCOL, NO2TRCOL, NBTR1AR, NO1TRAR )
         TRACTE = .FALSE.

      ENDIF

C     BILAN SUR LA QUALITE DE LA TRIANGULATION DES SURFACES FERMEES
C     =============================================================
      PRINT*
      AREMOY = AREMOY / NBFACO
      QTRMOY = QTRMOY / NBFACO

      PRINT*
      PRINT*,' otajsf: VOLUME: ', KNMVOL

      IF( LANGAG .EQ. 0 ) THEN
         PRINT 10700, NBSTID,NBSOMM-NBSOMM0-NBSTID,NBSOMM,
     %                NBFAISOL,NBFACO,QTRMIN,QTRMOY,
     %                ACOS(COSMAX)*RADEGR,NBANMI,ACOS(COSMX)*RADEGR,
     %                ACOS(COSMIN)*RADEGR,NBANMA,ACOS(COSMN)*RADEGR,
     %                AREMIN, AREMAX, AREMOY
      ELSE
         PRINT 20700, NBSTID,NBSOMM-NBSOMM0-NBSTID,NBSOMM,
     %                NBFAISOL,NBFACO,QTRMIN,QTRMOY,
     %                ACOS(COSMAX)*RADEGR,NBANMI,ACOS(COSMX)*RADEGR,
     %                ACOS(COSMIN)*RADEGR,NBANMA,ACOS(COSMN)*RADEGR,
     %                AREMIN, AREMAX, AREMOY
      ENDIF

10700 FORMAT(/' NOMBRE DE SOMMETS   IDENTIFIES  =',I7/
     %        ' NOMBRE DE SOMMETS OT   AJOUTES  =',I7/
     %        ' NOMBRE FINAL    DE     SOMMETS  =',I7/
     %        ' NOMBRE de TRIANGLES  SUPPRIMES  =',I7/
     %        ' NOMBRE  TOTAL  de  TRIANGLES    =',I7/
     %        ' QUALITE de la TRIANGULATION de la PEAU DE L''OBJET'/
     %        ' QUALITE MINIMALE des TRIANGLES  =',G15.6/
     %        ' QUALITE MOYENNE  des TRIANGLES  =',G15.6/
     %        ' ANGLE MINIMAL    des TRIANGLES  =',F7.2,' DEGRES  et',
     %        I6,' ANGLES < ',F7.2,' DEGRES' /
     %        ' ANGLE MAXIMAL    des TRIANGLES  =',F7.2,' DEGRES  et',
     %        I6,' ANGLES > ',F7.2,' DEGRES' /
     %        ' LONGUEUR de la PLUS PETITE ARETE=',G15.6/
     %        ' LONGUEUR de la PLUS GRANDE ARETE=',G15.6/
     %        ' LONGUEUR de la     MOYENNE ARETE=',G15.6/)

20700 FORMAT(/' IDENTIFIED VERTICES NUMBER =',I7/
     %        ' ADDED  OT  VERTICES NUMBER =',I7/
     %        ' FINAL      VERTICES NUMBER =',I7/
     %        ' DELETED  TRIANGLES  NUMBER =',I7/
     %        ' BOUNDARY TRIANGLES  NUMBER =',I7/
     %        ' TRIANGLES  MINIMUM  QUALITY=',G15.6/
     %        ' TRIANGLES  AVERAGE  QUALITY=',G15.6/
     %        ' TRIANGLES  MINIMUM  ANGLE  =',F7.2,' DEGREES  and',
     %        I6,' ANGLES < ',F7.2,' DEGREES' /
     %        ' TRIANGLES  MAXIMUM  ANGLE  =',F7.2,' DEGREES  and',
     %        I6,' ANGLES > ',F7.2,' DEGREES' /
     %        ' SMALLEST  EDGE  LENGTH     =',G15.6/
     %        ' GREATEST  EDGE  LENGTH     =',G15.6/
     %        ' AVERAGE   EDGE  LENGTH     =',G15.6/)

      IF( NBANMI + NBANMA .GT. 0 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            PRINT*,'L''OPTION 29 ou 30 des SURFACES est RECOMMANDEE afin
     % de FAVORISER une BONNE TETRAEDRISATION'
         ELSE
            PRINT*,'USE the OPTION 29 or 30 of SURFACES to FAVORIZE a GO
     %OD TETRAHEDRIZATION'
         ENDIF
      ENDIF


      IF( AREMIN .LE. 0 ) IERR = 3

C     DESTRUCTION DU TABLEAU OTFL
 9999 CALL TNMCDS( 'ENTIER', MXOTFL, MNOTFL )


c     debog
      PRINT*
      PRINT*,'otajsf: Volume ',KNMVOL
      IF( LANGAG .EQ. 0 ) THEN
         PRINT*,'HEXAEDRE des Sommets Ajoutes MIN',(hexaed(kk,1),kk=1,3)
         PRINT*,'HEXAEDRE des Sommets Ajoutes MAX',(hexaed(kk,2),kk=1,3)
      ELSE
         PRINT*,'Added Vertices MINIMUM XYZ',(hexaed(kk,1),kk=1,3)
         PRINT*,'Added Vertices MAXIMUM XYZ',(hexaed(kk,2),kk=1,3)
      ENDIF
      PRINT*

      TRACTE = TRACTE0
      RETURN
      END
