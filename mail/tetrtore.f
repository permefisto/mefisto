      SUBROUTINE TETRTORE( QUAMINEX, KTITRE,
     %                     MXSOMM, NBSOMM0,  NBSOMM, PTXYZD,
     %                     NPSOFR, VMOYEN,
     %                     MXTETR, N1TEVI, NOTETR, NUDTETR, N1TETS,
     %                     NBTRCF, NOTRCF, LEFACO, NO0FAR,
     %                     N1FEVI, MXFETO, NFETOI,
     %                     MXSSET, NBSSET0,N1SSET, NBTEETPR,
     %                     MXTECR, NBTECR, NOTECR, VOLET1, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : TETRAEDRISATION DES NBSSET0 SOUS-ETOILES DE L'ETOILE DES FACES
C ----- DEFINIES DANS N1SSET et NFETOI d'une CONFIGURATION EN TORE
C       LE NOMBRE DE SOUS-ETOILES EVOLUE AU COURS DE LA TETRAEDRISATION
C       DANS CHAQUE SOUS-ETOILE LES FACES SONT A NORMALE VERS L'INTERIEUR

C ENTREES:
C --------
C QUAMINEX:QUALITE MINIMALE EXIGEE POUR CREER UN TETRAEDRE
C KTITRE : TITRE DU TRACE
C MXSOMM : NOMBRE MAXIMAL DE SOMMETS DECLARABLES DANS PTXYZD NPSOFR
C MXTETR : NOMBRE DE TETRAEDRES DECLARABLES DANS NOTETR
C MXFETO : NOMBRE DE FACES DECLARABLES DANS NFETOI
C MXTECR : NOMBRE MAXIMAL DE TETRAEDRES CREABLES DANS NOTECR
C VMOYEN : VOLUME MOYEN DES TETRAEDRES INITIAUX DE L'ETOILE
C          UTILE POUR ELIMINER LES TETRAEDRES DE VOLUME PROCHE DE ZERO
C NBSOMM0: NOMBRE DE SOMMETS DANS PTXYZD AVANT TRAITEMENT DE CETTE FACE PERDUE

C MODIFIES:
C ---------
C NBSOMM : NOMBRE DE SOMMETS CREES DANS PTXYZD
C PTXYZD : TABLEAU DES COORDONNEES DES POINTS
C          PAR POINT : X  Y  Z DISTANCE_SOUHAITEE
C NPSOFR : =  0 SI POINT AJOUTE NON SUR LA SURFACE DE L'OBJET
C          =  1 SI LE POINT EST AJOUTE SUR UNE ARETE DE LA SURFACE
C          =  2 SI LE POINT EST AJOUTE COMME BARYCENTRE SUR LA SURFACE
C          =  1 000 000  + NUMERO DU  POINT INTERNE UTILISATEUR
C          = (1 000 000) * (NO SURFACE + 1) + NUMERO DU SOMMET DANS
C              LA NUMEROTATION DES SOMMETS DE LA SURFACE
C          = -4 SI LE POINT EST SOMMET D'OT NON TROP PROCHE PT OU FACE
C          = -1 SI LE POINT EST SOMMET D'OT TROP PROCHE PT OU FACE
C          = -3 SI LE POINT EST SOMMET D'OT REFUSE DANS LA TETRAEDRISATION
C          = -NPSOFR(I) SI POINT DEPLACE SUR LA SURFACE
C            DANS LEUR SURFACE FERMEE INITIALE ou NO DE POINT INTERNE
C N1TEVI : NO DU PREMIER TETRAEDRE VIDE DANS NOTETR
C NOTETR : LISTE DES TETRAEDRES
C          SOMMET1,    SOMMET2,    SOMMET3,    SOMMET4,
C          TETRAEDRE1, TETRAEDRE2, TETRAEDRE3, TETRAEDRE4
C          DE L'AUTRE COTE DE LA FACE
C          1: 123      2: 234      3: 341      4: 412
C NUDTETR: PLUS GRAND NUMERO DANS NOTETR DE TETRAEDRE OCCUPE

C N1TETS : N1TETS(I) NUMERO NOTETR D'UN TETRAEDRE AYANT POUR SOMMET I

C NBTRCF : NOMBRE DE TRIANGLES PERDUS FORMANT LE CF
C NOTRCF : SI NOTRCF(*)>0 NUMERO LEFACO DU TRIANGLE
C                      <0 NUMERO NO0FAR DU TRIANGLE AJOUTE AU CF
C LEFACO : FACE=TRIANGLE DE LA PEAU OU DES INTERFACES ENTRE VOLUMES
C          IL CONTIENT DANS CET ORDRE
C          NUMERO (DANS PTXYZD) DU SOMMET 1, SOMMET 2, SOMMET 3
C          NUMERO (DANS NUVOPA 0 SINON) DU VOLUME1, VOLUME2 DE LA FACE
C          NUMERO (DANS LEFACO) DE LA FACE ADJACENTE PAR L'ARETE 1 2 3
C
C          ATTENTION: UNE ARETE PEUT APPARTENIR A PLUS DE 2 FACES
C          => CHAINAGE CIRCULAIRE DE CES FACES DANS LEFACO
C          LEFACO(9,*)  -> FACE SUIVANTE (*=0:VIDE, *<>0:NON VIDE)
C          HACHAGE AVEC LA SOMME DES 3 SOMMETS MODULO MXFACO
C          LF = MOD( NOSOFA(1)+NOSOFA(2)+NOSOFA(3) , MXFACO ) + 1
C          NF = LEFACO( 10, LF ) LE NUMERO DE LA 1-ERE FACE DANS LEFACO
C          SI LA FACE NE CONVIENT PAS. PASSAGE A LA SUIVANTE
C          NF  = LEFACO( 9, NF )  ...
C          LEFACO(10,*) PREMIERE FACE DANS LE HACHAGE
C          LEFACO(11,.) = NO NOTETR D'UN TETRAEDRE AYANT CETTE FACE, 0 SINON
cccC          LEFACO(12,.) = NO FACEOC DE 1 A NBFACES D'OC
C NO0FAR : NUMERO DES 3 SOMMETS DE LA FACE AJOUTEE AU CF
C          NORMALE VERS L'INTERIEUR DU TETRAEDRE LA CONTENANT

C N1FEVI : POINTEUR SUR LA PREMIERE FACE VIDE DE L'ETOILE
C          CHAINAGE SUIVANT DANS NFETOI(5,*)
C NFETOI : VERSION 2 LES FACES TRIANGULAIRES DE L'ETOILE
C          1: NUMERO DU TETRAEDRE DANS NOTETR OPPOSE A CETTE FACE
C          2: NUMERO PTXYZD DU SOMMET 1 DE LA FACE
C          3: NUMERO PTXYZD DU SOMMET 2 DE LA FACE
C          4: NUMERO PTXYZD DU SOMMET 3 DE LA FACE
C             S1S2xS1S3 EST LA NORMALE DIRIGEE VERS L'INTERIEUR
C             DE L'ETOILE
C          5: CHAINAGE SUIVANT DES FACES OCCUPEES ET VIDES

C MXSSET : NOMBRE MAXIMAL DE SOUS ETOILES DECLARABLES DANS N1SSET
C NBSSET0: NOMBRE INITIAL DE SOUS-ETOILES STOCKEES DANS N1SSET-NFETOI
C N1SSET : NUMERO NFETOI DE LA 1-ERE FACE DES NBSSET0 SOUS ETOILES

C NBTEETPR:NOMBRE DE TETRAEDRES INITIAUX CONSTITUANTS L'ETOILE

C SORTIES:
C --------
C NBTECR : NOMBRE DE TETRAEDRES CREES AVANT ET APRES EXECUTION DU SP
C NOTECR : NUMERO NOTETR DES NBTECR TETRAEDRES CREES
C VOLET1 : VOLUME DES NBTECR TETRAEDRES CREES DANS L'ETOILE
C IERR   : 0 SI PAS D'ERREUR
C          1 ETOILE AVEC MOINS DE 4 FACES
C          2 SATURATION D'UN TABLEAU
C          3 3 SOMMETS NON DANS UN TETRAEDRE
C          4 NOMBRE INCORRECT DE FACES DANS L'ETOILE
C          5 SOMMETS DES 4 FACES DE L'ETOILE INCORRECTS
C          6 UNE ARETE DES FACES DE L'ETOILE APPARTIENT A MOINS DE 2
C            OU PLUS DE 2 FACES DE L'ETOILE
C          7 ou 8 ou 9 ALGORITHME DEFAILLANT A AMELIORER
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Saint Pierre du Perray           Janvier 2018
C2345X7..............................................................012
      PARAMETER         ( MXVPSI=2048, MXASFVP=512, MXCIAS=64 )
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"

      COMMON / TRTETR / STOPTE, TRACTE
      LOGICAL           STOPTE, TRACTE

      DOUBLE PRECISION  PTXYZD(1:4,1:NBSOMM)
      INTEGER           NPSOFR(MXSOMM), NOTETR(8,MXTETR),
     %                  N1TETS(1:MXSOMM), NOTRCF(NBTRCF), LEFACO(11,*),
     %                  NO0FAR(3,*), NFETOI(5,MXFETO), NOTECR(MXTECR),
     %                  N1SSET(MXSSET)

      CHARACTER*(*)     KTITRE
      DOUBLE PRECISION  CBTR(3), VOLET1, V, VOLMAX, VMOYEN, RAVEVM,
     %                  VOLTEC, XYZ(3), VETXYZ,
     %                  ARMIN, ARMAX, SURFTR(4)
      REAL              QUATET, QUAMIN, Q

      INTEGER           NFVPSI0(MXVPSI), NFVPSI(MXVPSI),
     %                  N1CIAS(MXCIAS), NSASFVP(4,MXASFVP)

      INTEGER           NOSOTE(4), NOSOTR2(3), NOSOTR3(3), NSARP2F(1)
      EQUIVALENCE      (NOSOTE(1),NS1), (NOSOTE(2),NS2),
     %                 (NOSOTE(3),NS3), (NOSOTE(4),NS4)
      INTEGER           NTEOPF(4)
      EQUIVALENCE      (NTEOPF(1),NTEOPF1), (NTEOPF(2),NTEOPF2),
     %                 (NTEOPF(3),NTEOPF3), (NTEOPF(4),NTEOPF4)
      INTEGER           NOSOFATE(3,4)
      DATA              NOSOFATE / 1,3,2,  2,3,4,  3,1,4,  4,1,2 /
C     => LA NORMALE A LA FACE EST ORIENTEE VERS L'EXTERIEUR DU TETRAEDRE

ccc      PRINT *
ccc      PRINT *,'tetrtore:',NBSSET0,' SOUS-ETOILES A TETRAEDRISER VMOYEN=',
ccc     %         VMOYEN

C     TRACE DES NBSSET0 SOUS-ETOILES INITIALES
      CALL TRFETO7( PTXYZD, 0, NFETOI, NBSSET0, N1SSET,
     %              0, NSARP2F )

C     NOMBRE DE TRAITEMENTS DE L'ETOILE
      NOAPPEL = 0
      IERR = 0

      MAXTETR = MIN( MAX(6,NBTRCF) * NBTEETPR, MXTECR )
ccc      MAXSOMM = 16 * NBTRCF
      MAXSOMM = 16

      NS4MAX00= -1
      NS4MAX0 = 0

C     NOMBRE DE FACES AVEC UN TETRAEDRE OPPOSE INCONNU
      NBFANR = 0

C     VOLUME DES NBTECR TETRAEDRES CREES
      VOLET1 = 0D0

C     NOMBRE INITIAL DE SOUS-ETOILES DE L'ETOILE
      NBSSET = NBSSET0

C     ======================================================================
C     TETRAEDRISATION DE LA PILE DES NBSSET0 SOUS-ETOILES DE L'ETOILE NFETOI
C     ======================================================================
 5    IF( NBSSET .LE. 0 ) GOTO 9000


C     TETRAEDRISATION DE LA SOUS-ETOILE NBSSET EN HAUT DE PILE N1SSET
C     ===============================================================
C     1-ERE FACE NFETOI VERSION 2 DE LA SOUS-ETOILE NBSSET DE NFETOI
      N1FEOC = N1SSET( NBSSET )
ccc      PRINT*
ccc      PRINT*,'tetrtore: LES FACES DE LA SOUS-ETOILE',NBSSET
ccc      CALL AFETOI( N1FEOC, NFETOI )

C     VERIFIER QUE TOUTE ARETE DE LA SOUS ETOILE APPARTIENT
C     SEULEMENT A 2n FACES DE L'ETOILE avec n=1,2
      CALL VETAFET( NOTETR, N1FEOC, NFETOI, NBFETO, NBARPB )

      IF( NBARPB .GT. 0 ) THEN
C        TRACE DE LA SOUS-ETOILE N1FEOC A TETRAEDRISER
         KTITRE='ATTENTION LA SOUS-ETOILE A TETRAEDRISER avec ARETES A P
     %ROBLEME'
         CALL SANSDBL( KTITRE, L )
         CALL TRFETO4( KTITRE(1:L), PTXYZD, N1FEOC, NFETOI,
     %                 NBTRCF, NOTRCF, LEFACO, NO0FAR,
     %                 NBTECR, NOTECR, NOTETR )
      ENDIF

C     NOMBRE DE TETRAEDRES CREES ACTUELLEMENT DANS L'ETOILE
      NBTECR0 = NBTECR
      NOCHOIXYZ = 0

C     NBFETO NOMBRE DE FACES DE LA SOUS-ETOILE NBSSET
C     -----------------------------------------------
      NBFETO = 0
      NF00   = 0
      NF0    = N1FEOC
 10   IF( NF0 .GT. 0 ) THEN
C        UNE FACE DE PLUS
         NBFETO = NBFETO + 1
C        NO DU 1-ER SOMMET DE LA FACE EST RENDU POSITIF
         NFETOI(2,NF0) = ABS( NFETOI(2,NF0) ) 
C        LA FACE SUIVANTE DE NF0 DANS NFETOI
         NF00 = NF0
         NF0  = NFETOI(5,NF0)
         GOTO 10
      ENDIF
      PRINT *,'tetrtore: SOUS-ETOILE',NBSSET,
     %        ' NB FACES de la SOUS-ETOILE=',NBFETO,' N1FEOC=',N1FEOC

      N = NBSOMM - NBSOMM0
      IF( NBTECR.GT.MAXTETR .OR. N.GT.MAXSOMM )THEN
         print*,'?????????????????????????????????????????????????????'
         print*,'PB tetrtore: TROP DE FACES dans l''ETOILE  ',NBFETO
         print*,'PB tetrtore: TROP DE TETRAEDRES CREES ou  ',NBTECR
         print*,'PB tetrtore: TROP DE SOMMETS    CREES     ',N
         print*,'PB tetrtore: NOMBRE DE TETRAEDRES INITIAUX',NBTEETPR
         print*,'?????????????????????????????????????????????????????'
         IERR = 8
         GOTO 9999
      ENDIF

C     LA SOUS-ETOILE EST DEPILEE CAR EN COURS DE TETRAEDRISATION
      NBSSET = NBSSET - 1

C     --------------------------------------------------------
C     MODE DE CREATION DES TETRAEDRES SELON LE NOMBRE DE FACES
C     --------------------------------------------------------
      IF( NBFETO .LE. 0 ) THEN
C        PLUS AUCUNE FACE DANS LA SOUS-ETOILE NBSSET+1
C        PASSAGE A LA SOUS-ETOILE EN HAUT DE PILE
         GOTO 5
      ENDIF

      IF( NBFETO .LE. 3 ) THEN

C        MOINS DE 4 FACES DANS LA SOUS-ETOILE
C        ------------------------------------
         IF( NBFETO .EQ. 2 ) THEN
C           CES 2 FACES SONT ELLES IDENTIQUES? (CAS 1 SEULE DEMI-ETOILE)
            NF0 = N1FEOC
            NF1 = NFETOI(5,NF0)
            DO K=1,3
               NOSOTR2(K) = NFETOI(1+K,NF0)
               NOSOTR3(K) = NFETOI(1+K,NF1)
            ENDDO
            CALL TRI3NO( NOSOTR2, NOSOTR2 )
            CALL TRI3NO( NOSOTR3, NOSOTR3 )
            IF( NOSOTR2(1) .EQ. NOSOTR3(1) .AND.
     %          NOSOTR2(2) .EQ. NOSOTR3(2) .AND.
     %          NOSOTR2(3) .EQ. NOSOTR3(3) ) THEN
               PRINT *,'tetrtore: SOUS-ETOILE',NBSSET+1,
     %                 ' AVEC 2 FACES IDENTIQUES RETIREES', NOSOTR2
               GOTO 5
            ENDIF
         ENDIF

         PRINT*,'PB tetrtore: NOMBRE DE FACES=',NBFETO,
     %          ' INCORRECT POUR FORMER UN TETRAEDRE'
         KTITRE='PB tetrtore: NBFETO=      '
         WRITE(KTITRE(21:22),'(I2)') NBFETO
         DO K=1,NBSSET+1
            PRINT*
            PRINT*,'SOUS ETOILE',K
            NF = N1SSET( K )
            CALL AFETOI( NF, NFETOI )
         ENDDO

C        TENTATIVE DE JOINDRE CETTE SOUS-ETOILE INCORRECTE A
C        LA SOUS ETOILE QUI LA PRECEDE
         IF( NBSSET .GT. 0 ) THEN
ccc            TRACTE = .TRUE.
            PRINT*,'JONCTION DES SOUS-ETOILES',NBSSET,NBSSET+1
             NFETOI(5,NF00) = N1SSET( NBSSET )
             N1SSET( NBSSET ) = N1FEOC
             GOTO 5
         ENDIF 

C        FIN DE SOUS ETOILE AVEC ERREUR
         IERR = 4
         GOTO 9010

      ENDIF

      IF( NBFETO .EQ. 4 ) THEN

C        -----------------------------------------------------------
C        FIN) LES 4 SEULES FACES DE LA SOUS-ETOILE DOIVENT FORMER UN
C             TETRAEDRE SINON ERREUR
C        -----------------------------------------------------------
C        TETRAEDRE OPPOSE A LA FACE N1FEOC
         NTEOPF(1) = NFETOI(1,N1FEOC)

C        LE NUMERO PTXYZD DES 4 SOMMETS DU TETRAEDRE NOSOTE
         NS1 = NFETOI(2,N1FEOC)
         NS2 = NFETOI(3,N1FEOC)
         NS3 = NFETOI(4,N1FEOC)
C        LE 4-EME SOMMET EST CELUI DES FACES DIFFERENT DES 3 DE LA FACE 1
C        LA FACE 2 DE LA SOUS-ETOILE DE 4 FACES
         NF2 = NFETOI(5,N1FEOC)
         DO K=1,3
            NS4 = NFETOI(1+K,NF2)
            IF( NS4 .NE. NS1 .AND. NS4 .NE. NS2 .AND. NS4 .NE. NS3 )
     %          GOTO 12
         ENDDO

C        RECHERCHE DES TETRAEDRES OPPOSES AUX 3 DERNIERES FACES
C        VERIFICATION: LES 3 DERNIERES FACES DOIVENT AVOIR LEURS 3
C        SOMMETS PARMI LES 4 DU TETRAEDRE
 12      DO 20 K=2,4

C           LES 3 SOMMETS DE LA FACE K DU TETRAEDRE NS1 NS2 NS3 NS4
            DO L=1,3
               NOSOTR2(L) = NOSOTE( NOSOFATE(L,K)  )
            ENDDO
C           TRI CROISSANT DES 3 NUMEROS DE SOMMETS DE LA FACE K
            CALL TRI3NO( NOSOTR2, NOSOTR2 )

            NF = NF2
C           LES 3 SOMMETS DE LA FACE NF DE LA SOUS-ETOILE
 16         DO L=1,3
               NOSOTR3(L) = NFETOI(1+L,NF)
            ENDDO
C           TRI CROISSANT DES 3 NUMEROS DE SOMMETS DE LA FACE NF
            CALL TRI3NO( NOSOTR3, NOSOTR3 )

            IF( NOSOTR2(1) .EQ. NOSOTR3(1) .AND.
     %          NOSOTR2(2) .EQ. NOSOTR3(2) .AND.
     %          NOSOTR2(3) .EQ. NOSOTR3(3) ) THEN
C               LA FACE K DU TETRAEDRE EST LA FACE NF DE LA SOUS-ETOILE
C               LE NO DU TETRAEDRE OPPOSE A LA FACE K DU TETRAEDRE
                NTEOPF(K) = NFETOI(1,NF)
                GOTO 20
            ENDIF

C           FACE SUIVANTE
            NF = NFETOI(5,NF)
            IF( NF .NE. 0 ) GOTO 16

            PRINT*,'PB tetrtore: FACE',K,'DU TETRAEDRE de St:',NOSOTE,
     %' N''EST PAS UNE DES 4 SEULES FACES de l''ETOILE => ETOILE INCORRE
     %CTE'
            IERR  = 5
            NOCAS = 0
            KTITRE='PB tetrtore: FACE NON TROUVEE DANS LE TETRAEDRE    '
            WRITE( KTITRE(49:92),'(4I11)') (NOSOTE(L),L=1,4)
            GOTO 9010

 20      ENDDO

C        AJOUT DU TETRAEDRE FORME DES 4 SEULES FACES DE LA SOUS-ETOILE
C        LES 3 SOMMETS DE LA FACE N1FEOC SONT DANS LE SENS DIRECT
C        DU TETRAEDRE CAR SA NORMALE EST DIRIGEE VERS L'INTERIEUR
         NS1F0 = NFETOI(2,N1FEOC)
         NS2F0 = NFETOI(3,N1FEOC)
         NS3F0 = NFETOI(4,N1FEOC)
         CALL AJTEET( NS1F0,   NS2F0,   NS3F0,   NS4,
     %                NTEOPF1, NTEOPF2, NTEOPF3, NTEOPF4,
     %                NBSOMM,  PTXYZD,  QUAMINEX,
     %                MXTETR,  N1TEVI,  NOTETR,  NUDTETR, N1TETS,
     %                MXTECR,  NBTECR,  NOTECR,  VOLET1,
     %                VOLTEC,  QUATET,  IERR )
         IF( IERR .NE. 0 ) GOTO 9999

C        LE TETRAEDRE NBTECR EST CREE
         NOCAS = 1
         NTE   = NOTECR( NBTECR )
         PRINT*,'tetrtore: NOCAS=',NOCAS,': 4Faces=> TETRAEDRE+',
     %           NTE,'=',(NOTETR(KK,NTE),KK=1,4),' N1TEVI=',N1TEVI,
     %     ' V=',VOLTEC,' Q=',QUATET

         GOTO 1000

      ENDIF


C     ----------------------------------------------------------------------------
C     1) ESSAI DE TETRAEDRISER LES TRIANGLES DE LA SOUS-ETOILE SANS AJOUT DE POINT
C        RECHERCHE DANS LA SOUS-ETOILE DE 3 FACES ADJACENTES FORMANT UN TETRAEDRE
C        DE VOLUME>0 ET D'ARETES SANS INTERSECTION AVEC LES AUTRES FACES
C     ----------------------------------------------------------------------------
      NF0 = N1FEOC
 51   IF( NF0 .GT. 0 ) THEN

C        BOUCLE SUR LES 3 ARETES DE LA FACE NF0
         DO K=1,3

C           L'ARETE K DE LA FACE NF0
C           NUMERO DES 3 SOMMETS AVEC LES 2 DE L'ARETE EN PREMIER
            IF( K .EQ. 3 ) THEN
               KK = 1
            ELSE
               KK = K + 1
            ENDIF
C           LES 2 SOMMETS DE L'ARETE K DE LA FACE NF0
            NS1F0 = NFETOI(1+K ,NF0)
            NS2F0 = NFETOI(1+KK,NF0)
C           LE SOMMET 3 DE LA FACE NF0
            IF( KK .EQ. 3 ) THEN
               KKK = 1
            ELSE
               KKK = KK + 1
            ENDIF
            NS3F0 = NFETOI(1+KKK,NF0)

C           RECHERCHE DE LA 2-EME FACE DE LA SOUS-ETOILE NON NF0 et D'ARETE K
            NF1 = N1FEOC
 52         IF( NF1 .GT. 0 ) THEN
               IF( NF1 .NE. NF0 ) THEN
                  DO L=1,3
C                    L'ARETE L DE LA FACE NF1
                     IF( L .EQ. 3 ) THEN
                        LL = 1
                     ELSE
                        LL = L + 1
                     ENDIF
C                    LES 2 SOMMETS DE L'ARETE L DE LA FACE NF1
                     NS1F1 = NFETOI(1+L ,NF1)
                     NS2F1 = NFETOI(1+LL,NF1)
                     IF( NS1F0 .EQ. NS2F1 .AND. NS2F0 .EQ. NS1F1 ) THEN
ccc     %          .OR. NS1F0 .EQ. NS1F1 .AND. NS2F0 .EQ. NS2F1 ) THEN
C                       L'ARETE K DE NF0 EST L'ARETE L DE NF1
C                       LE SOMMET 3 DE LA FACE NF1
                        IF( LL .EQ. 3 ) THEN
                           LLL = 1
                        ELSE
                           LLL = LL + 1
                        ENDIF
                        NS3F1 = NFETOI(1+LLL,NF1)

                        IF( NS3F1 .EQ. NS3F0 ) THEN
C                          CAS DE 2 FACES DU CF IDENTIQUES DANS L'ENCOCHE
                        print*,'tetrtore: 2FACES de l''etoile avec 3 MEM
     %ES SOMMETS (ENCOCHE POSSIBLE) NF0=',NF0,' ST:',NS1F0,NS2F0,NS3F0,
     %                         ' NF1=',NF1,' ST:',NS1F1,NS2F1,NS3F1
                           CALL AFETOI( N1FEOC, NFETOI )
                           GOTO 69
                        ENDIF

C                       LE VOLUME FACE NF0-SOMMET NS3F1 EST IL POSITIF?
C                       SON VOLUME ET SA QUALITE
                        CALL QUATETD( PTXYZD(1,NFETOI(2,NF0)),
     %                                PTXYZD(1,NFETOI(3,NF0)),
     %                                PTXYZD(1,NFETOI(4,NF0)),
     %                                PTXYZD(1,NS3F1),
     %                                ARMIN, ARMAX, SURFTR, V, Q )

                        IF( V .LE. 0D0 .OR. Q .LE. QUAMINEX ) THEN
C                          TETRAEDRE DE VOLUME<=0 ou QUALITE TROP FAIBLE
C                           => RECHERCHE D'UNE AUTRE FACE NF1
                           GOTO 69
                        ENDIF

C                       RECHERCHE DE LA 3-eme FACE NON NF0 et NON NF1
                        NF2 = N1FEOC
 53                     IF( NF2 .GT. 0 ) THEN
                           IF( NF2 .NE. NF0 .AND. NF2 .NE. NF1 ) THEN

                              DO 65 M=1,3

C                                LES 2 SOMMETS DE L'ARETE M DE LA FACE NF2
                                 NS1F2 = NFETOI(1+M ,NF2)
                                 IF( NS1F2 .NE. NS3F1 ) GOTO 65

C                                SOMMETS NS1F2=NS3F1
                                 IF( M .EQ. 3 ) THEN
                                    MM = 1
                                 ELSE
                                    MM = M + 1
                                 ENDIF
                                 NS2F2 = NFETOI(1+MM,NF2)

                                 IF( MM .EQ. 3 ) THEN
                                    MMM = 1
                                 ELSE
                                    MMM = MM + 1
                                 ENDIF
C                                LE SOMMET 3 DE LA FACE NF2
                                 NS3F2 = NFETOI(1+MMM,NF2)

                                 IF( NS2F2 .EQ. NS1F0 ) THEN
                                    IF( NS3F2 .EQ. NS3F0 ) THEN

C                                      NOCAS 2:  TETRAEDRE NF0-NS3F1 avec
C                                      NS1F2=NS3F1 et NS2F2=NS1F0 et NS3F2=NS3F0
C                                      -----------------------------------------
                                       NOCAS = 2

C                                      UNE DES 3 ARETES DE NF2 INTERSECTE T ELLE
C                                      LA FACE NS3F0 NS2F0 NS3F1?
                                       N0 = NFETOI(4,NF2)
                                       DO N=1,3
                                          N1 = NFETOI(1+N,NF2)
                                          CALL INARTR( PTXYZD(1,N0),
     %                                                 PTXYZD(1,N1),
     %                                                 PTXYZD(1,NS3F0),
     %                                                 PTXYZD(1,NS2F0),
     %                                                 PTXYZD(1,NS3F1),
     %                                                 LINTER,XYZ,CBTR)
                                          IF( LINTER  .EQ. 1 .AND.
     %                                     CBTR(1) .LT. 0.999999D0 .AND.
     %                                     CBTR(2) .LT. 0.999999D0 .AND.
     %                                     CBTR(3) .LT. 0.999999D0 )THEN
C                                          OUI: UN POINT D'INTERSECTION
C                                          ABANDON DE CE CHOIX DE 2 FACES
                                    print *,'bizarre nocas=2 cbtr=',cbtr
                                           GOTO 69
                                          ENDIF
                                          N0 = N1
                                       ENDDO

C                                      RECUPERER LE TETRAEDRE OPPOSE A LA
C                                      FACE NS2F0-NS3F0-NS3F1
C                                      DU TETRAEDRE NS1F0+NS2F0+NS3F0+NS3F1
                                       NTEOPF( 2 )= -1
                                       NOSOTR2(1) = NS2F0
                                       NOSOTR2(2) = NS3F0
                                       NOSOTR2(3) = NS3F1
                                       CALL TRI3NO( NOSOTR2, NOSOTR2 )

                                       NF3 = N1FEOC
 54                                    IF( NF3 .GT. 0 ) THEN
                                      CALL TRI3NO(NFETOI(2,NF3),NOSOTR3)

                                       IF(NOSOTR2(1).EQ.NOSOTR3(1) .AND.
     %                                    NOSOTR2(2).EQ.NOSOTR3(2) .AND.
     %                                    NOSOTR2(3).EQ.NOSOTR3(3) )THEN
C                                         TETRAEDRE OPPOSE A LA FACE 2
C                                         DU TETRAEDRE NS1F0-NS2F0-NS3F0-NS3F1
                                          NTEOPF( 2 ) = NFETOI( 1, NF3 )
                                          GOTO 55
                                       ENDIF
                                       NF3 = NFETOI(5,NF3)
                                       GOTO 54
                                       ENDIF

C                                      TETRAEDRE OPPOSE A LA FACE 3
C                                      DU TETRAEDRE NS1F0-NS2F0-NS3F0-NS3F1
 55                                    NTEOPF( 3 ) = NFETOI( 1, NF2 )
                                       NTEOPF( 4 ) = NFETOI( 1, NF1 )
                                       GOTO 60
                                    ENDIF
                                 ENDIF

C                                AUTRE POSSIBILITE
                                 IF( NS2F2 .EQ. NS3F0 ) THEN
                                    IF( NS3F2 .EQ. NS2F0 ) THEN

C                                      NOCAS 3:  TETRAEDRE NF0-NS3F1 avec
C                                      NS1F2=NS3F1 et NS2F2=NS3F0 et NS3F2=NS2F0
C                                      -----------------------------------------
                                       NOCAS = 3

C                                      UNE DES 3 ARETES DE NF2 INTERSECTE T ELLE
C                                      LA FACE NS3F0 NS3F1 NS1F0?
                                       N0 = NFETOI(4,NF2)
                                       DO N=1,3
                                          N1 = NFETOI(1+N,NF2)
                                          CALL INARTR( PTXYZD(1,N0),
     %                                                 PTXYZD(1,N1),
     %                                                 PTXYZD(1,NS3F0),
     %                                                 PTXYZD(1,NS3F1),
     %                                                 PTXYZD(1,NS1F0),
     %                                                 LINTER,XYZ,CBTR)
                                          IF( LINTER  .EQ. 1 .AND.
     %                                      CBTR(1).LT.0.999999D0 .AND.
     %                                      CBTR(2).LT.0.999999D0 .AND.
     %                                      CBTR(3).LT.0.999999D0)THEN
C                                           OUI: UN POINT D'INTERSECTION
C                                           ABANDON DE CE CHOIX DE 2 FACES
                                    print *,'bizarre nocas=3 cbtr=',cbtr
                                            GOTO 69
                                          ENDIF
                                          N0 = N1
                                       ENDDO

C                                      RECUPERER LE TETRAEDRE OPPOSE A LA
C                                      FACE NS1F0-NS3F0-NS3F1
C                                      DU TETRAEDRE NS1F0+NS2F0+NS3F0+NS3F1
                                       NTEOPF( 3) = -1
                                       NOSOTR2(1) = NS1F0
                                       NOSOTR2(2) = NS3F0
                                       NOSOTR2(3) = NS3F1
                                       CALL TRI3NO( NOSOTR2, NOSOTR2 )

                                       NF3 = N1FEOC
 57                                    IF( NF3 .GT. 0 ) THEN
                                      CALL TRI3NO(NFETOI(2,NF3),NOSOTR3)
                                      IF(NOSOTR2(1) .EQ. NOSOTR3(1).AND.
     %                                   NOSOTR2(2) .EQ. NOSOTR3(2).AND.
     %                                   NOSOTR2(3) .EQ. NOSOTR3(3))THEN
                                         NTEOPF(3) = NFETOI( 1, NF3 )
                                         GOTO 59
                                      ENDIF
                                          NF3 = NFETOI(5,NF3)
                                          GOTO 57
                                       ENDIF

C                                      TETRAEDRE OPPOSE A LA FACE 3
 59                                    NTEOPF( 2 ) = NFETOI( 1, NF2 )
                                       NTEOPF( 4 ) = NFETOI( 1, NF1 )
                                       GOTO 60
                                    ENDIF
                                 ENDIF
                                 GOTO 65

C                                TETRAEDRE OPPOSE PAR LA 4-EME FACE DE NTE
 60                              NTEOPF(4) = NFETOI(1,NF1)

C                                LES 3 FACES NF0 NF1 NF2 FORMENT UN TETRAEDRE
C                                NS1F0-NS2F0-NS3F0-NS3F1 A AJOUTER
C                                CREATION DU TETRAEDRE NS1F0-NS2F0-NS3F0-NS3F1
                                 NFETOI(2,NF0) = NS1F0
                                 NFETOI(3,NF0) = NS2F0
                                 NFETOI(4,NF0) = NS3F0
                                 NTEOPF1 = NFETOI(1,NF0)

                                 CALL AJTEET(NS1F0, NS2F0, NS3F0, NS3F1,
     %                               NTEOPF1, NTEOPF2, NTEOPF3, NTEOPF4,
     %                               NBSOMM,  PTXYZD, QUAMINEX,
     %                               MXTETR,  N1TEVI, NOTETR,
     %                               NUDTETR, N1TETS,
     %                               MXTECR,  NBTECR, NOTECR, VOLET1,
     %                               VOLTEC,  QUATET, IERR )
                                 IF( IERR .NE. 0 ) GOTO 9999

                                 NTE = NOTECR( NBTECR )
                                 print *,'tetrtore: NOCAS=',NOCAS,
     %                           ': 3Faces=> TETRAEDRE+',NTE,'=',
     %                           (NOTETR(KK,NTE),KK=1,4),' N1TEVI=',
     %                            N1TEVI,' V=',VOLTEC,' Q=',QUATET

C                                NOUVELLE ETOILE APRES AJOUT DU TETRAEDRE NTE
                                 GOTO 1000

 65                           CONTINUE
                           ENDIF

C                          LA FACE SUIVANTE DE NF2
                           NF2 = NFETOI(5,NF2)
                           GOTO 53

                        ENDIF
                     ENDIF
                  ENDDO
               ENDIF

C              LA FACE SUIVANTE DE NF1
 69            NF1 = NFETOI(5,NF1)
               GOTO 52
            ENDIF

         ENDDO

C        LA FACE SUIVANTE DE NFETOI
         NF0 = NFETOI(5,NF0)
         GOTO 51

      ENDIF

C     ----------------------------------------------------------------------------
C     2) ESSAI DE TETRAEDRISER LES TRIANGLES DE LA SOUS-ETOILE SANS AJOUT DE POINT
C        RECHERCHE DANS LA SOUS-ETOILE DE 2 FACES ADJACENTES FORMANT UN TETRAEDRE
C        DE VOLUME>0 ET DONT L'ARETE OPPOSEE INTERSECTE AUCUNE FACE
C        DE LA SOUS-ETOILE
C     ----------------------------------------------------------------------------
      NF0 = N1FEOC
 151  IF( NF0 .GT. 0 ) THEN

C        BOUCLE SUR LES 3 ARETES DE LA FACE NF0
         DO K=1,3

C           L'ARETE K DE LA FACE NF0
            IF( K .EQ. 3 ) THEN
               KK = 1
            ELSE
               KK = K + 1
            ENDIF
C           LES 2 SOMMETS DE L'ARETE K DE LA FACE NF0
            NS1F0 = NFETOI(1+K ,NF0)
            NS2F0 = NFETOI(1+KK,NF0)

C           LE SOMMET 3 DE LA FACE NF0
            IF( KK .EQ. 3 ) THEN
               KKK = 1
            ELSE
               KKK = KK + 1
            ENDIF
            NS3F0 = NFETOI(1+KKK,NF0)

C           RECHERCHE DE LA 2-EME FACE NF1 DE LA SOUS-ETOILE NON NF0 et D'ARETE K
            NF1 = N1FEOC
 152        IF( NF1 .GT. 0 ) THEN
               IF( NF1 .NE. NF0 ) THEN
                  DO L=1,3
C                    L'ARETE L DE LA FACE NF1
                     IF( L .EQ. 3 ) THEN
                        LL = 1
                     ELSE
                        LL = L + 1
                     ENDIF
C                    LES 2 SOMMETS DE L'ARETE L DE LA FACE NF1
                     NS1F1 = NFETOI(1+L ,NF1)
                     NS2F1 = NFETOI(1+LL,NF1)
                     IF( NS1F0 .EQ. NS2F1 .AND. NS2F0 .EQ. NS1F1 ) THEN
c     %          .OR. NS1F0 .EQ. NS1F1 .AND. NS2F0 .EQ. NS2F1 ) THEN
C                       L'ARETE K DE NF0 EST L'ARETE L DE NF1
C                       LE SOMMET 3 DE LA FACE NF1
                        IF( LL .EQ. 3 ) THEN
                           LLL = 1
                        ELSE
                           LLL = LL + 1
                        ENDIF
                        NS3F1 = NFETOI(1+LLL,NF1)

                        IF( NS3F1 .EQ. NS3F0 ) THEN
C                          CAS DE 2 FACES DU CF IDENTIQUES DANS L'ENCOCHE
                       print*,'Pb tetrtore: 2FACES avec 3 MEMES SOMMETS'
     %                        ,' NF0=',NF0,' NF1=',NF1,
     %                         ' ST:',NS1F0,NS2F0,NS3F0
                           CALL AFETOI( N1FEOC, NFETOI )
                           GOTO 159
                        ENDIF

C                       NS1F0-NS2F0-NS3F0-NS3F1 TETRAEDRE DE VOLUME>0?
C                       SON VOLUME ET SA QUALITE
                        CALL QUATETD( PTXYZD(1,NS1F0),
     %                                PTXYZD(1,NS2F0),
     %                                PTXYZD(1,NS3F0),
     %                                PTXYZD(1,NS3F1),
     %                                ARMIN, ARMAX, SURFTR, V, Q )

                        IF( V.LE.VMOYEN*1D-3 .OR. Q.LE.QUAMINEX ) THEN
C                          TETRAEDRE DE VOLUME<=0 OU PROCHE 0
C                          ou QUALITE TROP FAIBLE
C                          => RECHERCHE D'UNE AUTRE FACE NF1
                           GOTO 159
                        ENDIF

C                       L'ARETE NS3F0-NS3F1 INTERSECTE T ELLE UNE FACE DE LA SOUS-ETOILE?
                        NF2 = N1FEOC
 153                    IF( NF2 .GT. 0 ) THEN
                           IF( NF2 .NE. NF0 .AND. NF2 .NE. NF1 ) THEN

C                             CETTE FACE NF2 EST ELLE INTERSECTEE PAR 
C                             L'ARETE NS3F0-NS3F1?
                              CALL INARTR( PTXYZD(1,NS3F0),
     %                                     PTXYZD(1,NS3F1),
     %                                     PTXYZD(1,NFETOI(2,NF2)),
     %                                     PTXYZD(1,NFETOI(3,NF2)),
     %                                     PTXYZD(1,NFETOI(4,NF2)),
     %                                     LINTER, XYZ, CBTR )
C                             LINTER : -2 SI S1=S2
C                             -1 SI S1-S2 PARALLELE AU PLAN DU TRIANGLE
C                              0 SI S1-S2 N'INTERSECTE PAS LE TRIANGLE
C                              1 SI S1-S2   INTERSECTE     LE TRIANGLE ET ENTRE S1-S2
C                            XYZ : 3 COORDONNEES DU POINT D'INTERSECTION SI LINTER=1
C                            CBTR: 3 COORDONNEES BARYCENTRIQUES DE XYZ DANS LE TRIANGLE
                              IF(  LINTER  .EQ. 1 .AND.
     %                             CBTR(1) .LT. 0.999999D0 .AND.
     %                             CBTR(2) .LT. 0.999999D0 .AND.
     %                             CBTR(3) .LT. 0.999999D0 ) THEN
C                                OUI: UN POINT D'INTERSECTION
C                                ABANDON DE CE CHOIX DE 2 FACES
                                 GOTO 159
                              ENDIF

C                             UNE DES 3 ARETES DE NF2 INTERSECTE T ELLE
C                             LA FACE NS3F0 NS2F0 NS3F1?
                              N0 = NFETOI(4,NF2)
                              DO M=1,3
                                 N1 = NFETOI(1+M,NF2)
                                 CALL INARTR( PTXYZD(1,N0),
     %                                        PTXYZD(1,N1),
     %                                        PTXYZD(1,NS3F0),
     %                                        PTXYZD(1,NS2F0),
     %                                        PTXYZD(1,NS3F1),
     %                                        LINTER, XYZ, CBTR )
                                 IF(  LINTER  .EQ. 1 .AND.
     %                                CBTR(1) .LT. 0.999999D0 .AND.
     %                                CBTR(2) .LT. 0.999999D0 .AND.
     %                                CBTR(3) .LT. 0.999999D0 ) THEN
C                                     OUI: UN POINT D'INTERSECTION
C                                     ABANDON DE CE CHOIX DE 2 FACES
                                    GOTO 159
                                 ENDIF
                                 N0 = N1
                              ENDDO

C                             UNE DES 3 ARETES DE NF2 INTERSECTE T ELLE
C                             LA FACE NS3F0 NS3F1 NS1F0?
                              N0 = NFETOI(4,NF2)
                              DO M=1,3
                                 N1 = NFETOI(1+M,NF2)
                                 CALL INARTR( PTXYZD(1,N0),
     %                                        PTXYZD(1,N1),
     %                                        PTXYZD(1,NS3F0),
     %                                        PTXYZD(1,NS3F1),
     %                                        PTXYZD(1,NS1F0),
     %                                        LINTER, XYZ, CBTR )
                                 IF(  LINTER  .EQ. 1 .AND.
     %                                CBTR(1) .LT. 0.999999D0 .AND.
     %                                CBTR(2) .LT. 0.999999D0 .AND.
     %                                CBTR(3) .LT. 0.999999D0 ) THEN
C                                     OUI: UN POINT D'INTERSECTION
C                                     ABANDON DE CE CHOIX DE 2 FACES
                                    GOTO 159
                                 ENDIF
                                 N0 = N1
                              ENDDO

                           ENDIF

C                          LA FACE SUIVANTE DE NF2
                           NF2 = NFETOI(5,NF2)
                           GOTO 153

                        ENDIF

C                       CREATION DU TETRAEDRE NS1F0-NS2F0-NS3F0-NS3F1
C                       AVEC LES 2 FACES NF0-NF1 DE LA SOUS-ETOILE
C                       ---------------------------------------------
                        NOCAS = 4

C                       TETRAEDRE INCONNU OPPOSE AUX 2 AUTRES FACES DE NF0-NF1
                        NTEOPF2 = -1
                        NTEOPF3 = -1
                        NTEOPF4 = NFETOI(1,NF1)

C                       CREATION DU TETRAEDRE NS1F0-NS2F0-NS3F0-NS3F1
                        NFETOI(2,NF0) = NS1F0
                        NFETOI(3,NF0) = NS2F0
                        NFETOI(4,NF0) = NS3F0
                        NTEOPF1 = NFETOI(1,NF0)

                        CALL AJTEET( NS1F0,   NS2F0,   NS3F0,   NS3F1,
     %                               NTEOPF1, NTEOPF2, NTEOPF3, NTEOPF4,
     %                               NBSOMM,  PTXYZD,  QUAMINEX,
     %                          MXTETR, N1TEVI, NOTETR, NUDTETR,
     %                          N1TETS, MXTECR, NBTECR, NOTECR,  VOLET1,
     %                          VOLTEC, QUATET, IERR )
                        IF( IERR .NE. 0 ) GOTO 9999

                        NTE = NOTECR( NBTECR )
                        print *,'tetrtore: NOCAS=',NOCAS,
     %                  ': 2Faces=> TETRAEDRE+',NTE,'=',
     %                  (NOTETR(KK,NTE),KK=1,4),' N1TEVI=',N1TEVI,
     %                  ' V=',VOLTEC,' Q=',QUATET

C                       NOUVELLE ETOILE APRES AJOUT DU TETRAEDRE NTE
                        GOTO 1000

                     ENDIF
                  ENDDO
               ENDIF

C              LA FACE SUIVANTE DE NF1
 159           NF1 = NFETOI(5,NF1)
               GOTO 152
            ENDIF

         ENDDO

C        LA FACE SUIVANTE DE NFETOI
         NF0 = NFETOI(5,NF0)
         GOTO 151

      ENDIF

      IF( N1FEOC .LE. 0 ) THEN
C        TOUTE LA SOUS-ETOILE A ETE TETRAEDRISEE
         GOTO 1000
      ENDIF

C     -----------------------------------------------------------------------
C     3) RECHERCHE D'UNE FACE FORMANT AVEC UN SOMMET DE LA SOUS-ETOILE PROCHE
C        UN TETRAEDRE DE VOLUME POSITIF ET SANS INTERSECTION AVEC LES FACES
C     -----------------------------------------------------------------------
      NS4MAX = 0
      NF0MAX = 0
      VOLMAX = -1D123
      QUAMAX = 0.0
      NF0 = N1FEOC
 200  IF( NF0 .GT. 0 ) THEN

C        RECHERCHE D'UN SOMMET NS4 D'UNE FACE NF1 NON ADJACENTE A NF0
C        ET OPPOSEE DANS LA SOUS-ETOILE
C        ET SANS INTERSECTION AVEC LA SOUS-ETOILE
         NF1 = N1FEOC
 210     IF( NF1 .GT. 0 ) THEN
            IF( NF1 .NE. NF0 ) THEN

               DO 250 L=1,3

C                 LE SOMMET L DE LA FACE NF1
                  NS4 = NFETOI(1+L,NF1)
                  IF( NS4 .EQ. NFETOI(2,NF0)  .OR.
     %                NS4 .EQ. NFETOI(3,NF0)  .OR.
     %                NS4 .EQ. NFETOI(4,NF0) ) GOTO 250

C                 NS4 EST SOMMET DE NF1 MAIS PAS SOMMET DE NF0
C                 LE VOLUME DU TETRAEDRE FACE NF0 + NS4 EST-IL POSITIF?
                  CALL QUATETD( PTXYZD(1,NFETOI(2,NF0) ),
     %                          PTXYZD(1,NFETOI(3,NF0) ),
     %                          PTXYZD(1,NFETOI(4,NF0) ),
     %                          PTXYZD(1,NS4),
     %                          ARMIN, ARMAX, SURFTR, V, Q )

                  IF( V .LE. VMOYEN*1D-3 .OR. Q .LE. QUAMINEX ) GOTO 250

C                 LE TETRAEDRE NF0-NS4 EST IL SANS INTERSECTION AVEC LA SOUS-ETOILE?
                  NF2 = N1FEOC
 230              IF( NF2 .GT. 0 ) THEN
                     IF( NF2 .NE. NF0 ) THEN

C                       INTERSECTION TRIANGLE NF2-TETRAEDRE NF0+NS4?
                        CALL INTRITET( PTXYZD(1,NFETOI(2,NF2)),
     %                                 PTXYZD(1,NFETOI(3,NF2)),
     %                                 PTXYZD(1,NFETOI(4,NF2)),
     %                                 PTXYZD(1,NFETOI(2,NF0)),
     %                                 PTXYZD(1,NFETOI(3,NF0)),
     %                                 PTXYZD(1,NFETOI(4,NF0)),
     %                                 PTXYZD(1,NS4),
     %                                 LINTER )
                        IF( LINTER  .NE. 0 ) THEN
C                          OUI: UN POINT D'INTERSECTION=>ABANDON DE CE SOMMET NS4
                           GOTO 250
                        ENDIF

                     ENDIF

C                    LA FACE SUIVANTE DE NF2
                     NF2 = NFETOI(5,NF2)
                     GOTO 230

                  ENDIF

C                 TETRAEDRE FACE NF0 + NS4  de V>0 ET SANS INTERSECTION
ccc                  IF( V .GT. VOLMAX ) THEN  9/9/17
                  IF( Q .GT. QUAMAX ) THEN
C                    TETRAEDRE DE MEILLEURE QUALITE
                     VOLMAX = V
                     QUAMAX = Q
                     NF0MAX = NF0
                     NS4MAX = NS4
                     CALL QUATETD( PTXYZD(1,NFETOI(2,NF0)),
     %                             PTXYZD(1,NFETOI(3,NF0)),
     %                             PTXYZD(1,NFETOI(4,NF0)),
     %                             PTXYZD(1,NS4),
     %                             ARMIN, ARMAX, SURFTR, V, Q )
                     NOCAS = 5
                     print*,'tetrtore: NOCAS=',NOCAS,': 1F+1StPr NS4=',
     %                       NS4,' +NF0=',NF0,' V=',V,' Q=',Q,
     %                      ' EST POSSIBLE'
                  ENDIF
 250           ENDDO
            ENDIF

C           LA FACE SUIVANTE DE NF1
            NF1 = NFETOI(5,NF1)
            GOTO 210

         ENDIF

C        LA FACE SUIVANTE DE NFETOI
         NF0 = NFETOI(5,NF0)
         GOTO 200

      ENDIF

      IF( VOLMAX .EQ. -1D123 ) THEN
C        AUCUN SOMMET NS4 CONVIENT
         GOTO 300
      ENDIF

C     LE TETRAEDRE NF0MAX-NS4MAX EST CREE A TRAVERS NFVPSI
C     S'IL N'A PAS ETE LE POINT CHOISI PRECEDENT
      IF( NS4MAX .EQ. NS4MAX0 .OR. NS4MAX .EQ. NS4MAX00 ) THEN
C        ABANDON POUR EVITER DES CALCULS INUTILES
         PRINT*,'tetrtore: BOUCLE SUR 2 SOMMETS NS4MAX00=',NS4MAX00,
     %          ' NS4MAX0=',NS4MAX0,' et NS4MAX=',NS4MAX
         GOTO 300
      ENDIF

      NBVPSI = 1
      NFVPSI(1) = NF0MAX
      NS4MAX00= NS4MAX0
      NS4MAX0 = NS4MAX
      NS4     = NS4MAX
      NOCAS = 5
      print*,'tetrtore: NOCAS=',NOCAS,': 1F+1StPr avec NS4=',NS4,
     %' +NF0=',NF0,' de St',(NFETOI(M,NF0),M=2,4),
     %  ' V=',VOLMAX,' Q=',QUAMAX,' EST CHOISI'
      GOTO 320


C     ------------------------------------------------------------------
C     4) LA CONFIGURATION DES FACES DE LA SOUS-ETOILE
C        NE PERMET PAS AVEC 3 OU 2 OU 1 FACES DE FORMER UN TETRAEDRE
C      => UN  NOUVEAU POINT EST A LOCALISER PUIS AJOUTER ET TETRAEDRISER
C         AVEC LES FACES DE LA SOUS-ETOILE

C     CHOISIR UN POINT XYZ A AJOUTER ET SIMULER LA TETRAEDRISATION
C     DE TOUT OU UNE PARTIE DE LA SOUS-ETOILE => NBVPSI FACES NFVPSI DE
C     LA SOUS-ETOILE A TETRAEDRISER AVEC LE POINT XYZ
C     ------------------------------------------------------------------
 300  NOAPPEL = NOAPPEL + 1
      IF( NOAPPEL .GE. 8 ) THEN
C        ABANDON POUR EVITER DES CALCULS INUTILES
         IERR = 7
         GOTO 9999
      ENDIF

C     MAXFAC : NOMBRE MINIMUM DE FACES A VOLUME POSITIF AU DELA DUQUEL
C              ON SE CONTENTE DU MAX TROUVE POUR JOINDRE LE POINT XYZ AUX FACES
      MAXFAC = MAX( NBFETO/4, 3 )

      CALL SIPTETO( NOAPPEL,   KTITRE, VMOYEN, PTXYZD, NPSOFR,
     %              N1FEOC,    NFETOI, NBFETO,
     %              MAXFAC,    XYZ,
     %              MXVPSI,    NBVPSI,  NFVPSI, NFVPSI0,
     %              MXCIAS,    NBCIAS,  N1CIAS,
     %              MXASFVP,   NBASFVP, NSASFVP,
     %              QUAMINEX,  QUAMIN,  MIARSICI,
     %              NOCHOIXYZ, VETXYZ,  RAVEVM )
C     MIARSICI: NOMBRE D'ARETES SIMPLES DES CIRCUITS DES FACES NFVPSI

      IF( NBVPSI .LE. 0 ) THEN
C        ALGORITHME DE SIPTETO INSUFFISANT. L'ETOILE EST ECRASEE
C        ABANDON POUR EVITER DES CALCULS INUTILES
         IERR = 7
         GOTO 9999
      ENDIF

C     ------------------------------------------------------------
C     5) TETRAEDRISATION DES NBVPSI FACES NFVPSI DE LA SOUS-ETOILE
C        A PARTIR DU POINT XYZ COMME 4-EME SOMMET
C     ------------------------------------------------------------

C     TRACE EN ORANGE DES FACES+XYZ DONNANT UN VOLUME POSITIF SANS INTERSECTION
C     ET MAXIMISANT LE MINIMUM DES QUALITES DES TETRAEDRES FORMES
      KTITRE='tetrtore:       FACES+XYZ TETRAEDRES de Volume>0 NOCHOIXYZ
     %=        Qualite MIN=         et            FACES de l''ETOILE'
      WRITE(KTITRE(12:16), '(I5)'   ) NBVPSI
      WRITE(KTITRE(61:62), '(I2)'   ) NOCHOIXYZ
      WRITE(KTITRE(80:87), '(F8.5)' ) QUAMIN
      WRITE(KTITRE(92:96), '(I5)'   ) NBFETO
      CALL SANSDBL( KTITRE, L )
      PRINT*, KTITRE(1:L)
      CALL TRFETO3( KTITRE(1:L), XYZ, XYZ, PTXYZD, N1FEOC, NFETOI,
     %              NBVPSI, NFVPSI )

      NOCAS = 6 + NOCHOIXYZ
      IF( NBSOMM .GE. MXSOMM ) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:10),'(I10)') MXSOMM
         KERR(1) = 'tetrtore: MAXIMUM DE SOMMETS ' // KERR(MXLGER)(1:10)
     %             // ' ATTEINT'
         CALL LEREUR
         IERR = 2
         GOTO 9999
      ENDIF

C     CREATION D'UN NOUVEAU SOMMET XYZ DE LA TETRAEDRISATION
C     ------------------------------------------------------
      NBSOMM = NBSOMM + 1

C     XYZ COORDONNEES DU 4-EME SOMMET DES TETRAEDRES A FORMER
      DO L=1,3
         PTXYZD(L,NBSOMM) = XYZ(L)
      ENDDO
      PRINT*,'tetrtore: NOCAS=',NOCAS,': Pt AJOUTE NBSOMM=',NBSOMM,
     %       ' XYZ=',XYZ,' avec les', NBVPSI,' FACES+Pt=Volume>0'

C     POINT AJOUTE NON SUR LA SURFACE DE L'OBJET
      NPSOFR( NBSOMM ) = 0

C     NBSOMM EST LE 4-EME SOMMET DES NBVPSI TETRAEDRES
      NS4 = NBSOMM

C     CONSTRUCTION DES TETRAEDRES NBSOMM+FACES NFVPSI DE LA SOUS-ETOILE
C     -----------------------------------------------------------------
 320  DO 322 K = 1, NBVPSI

C        LA FACE A TETRAEDRISER
         NF0 = NFVPSI( K )
         IF( NF0 .LE. 0 ) GOTO 322

C        LA FACE NF0 + POINT NBSOMM FORMENT UN NOUVEAU TETRAEDRE
         NS1F0 = NFETOI(2,NF0)
         NS2F0 = NFETOI(3,NF0)
         NS3F0 = NFETOI(4,NF0)

C        LE VOLUME NS1 NS2 XYZE XYZ EST IL POSITIF?
         CALL QUATETD( PTXYZD(1,NS1F0), PTXYZD(1,NS2F0),
     %                 PTXYZD(1,NS3F0), PTXYZD(1,NBSOMM),
     %                 ARMIN, ARMAX, SURFTR, V, Q )

         IF( V .LE. 0D0 .OR. Q .LE. QUAMINEX ) THEN
C           TETRAEDRE ABANDONNE DE VOLUME<=0 ou QUALITE TROP FAIBLE
            GOTO 322
         ENDIF

C        NUMERO NOTETR DES TETRAEDRES OPPOSES AUX 4 FACES
         NTEOPF1 = NFETOI(1,NF0)
         NTEOPF2 = -1
         NTEOPF3 = -1
         NTEOPF4 = -1

         CALL AJTEET( NS1F0,   NS2F0,   NS3F0,   NS4,
     %                NTEOPF1, NTEOPF2, NTEOPF3, NTEOPF4,
     %                NBSOMM,  PTXYZD,  QUAMINEX,
     %                MXTETR,  N1TEVI, NOTETR, NUDTETR, N1TETS,
     %                MXTECR,  NBTECR, NOTECR, VOLET1,
     %                VOLTEC,  QUATET, IERR )

         NTE = NOTECR(NBTECR)
         print *,'tetrtore: NOCAS=',NOCAS,
     %           ': 1Face + 1XYZ => TETRAEDRE',NTE,'=',
     %            (NOTETR(KK,NTE),KK=1,4),' V=',VOLTEC,
     %           ' Q=',QUATET

 322  ENDDO


C     MISE A JOUR DES TETRAEDRES OPPOSES DES FACES DES TETRAEDRES
C     CREES AVEC LE SOMMET NS4 ET LES NBVPSI FACES NFVPSI
C     -----------------------------------------------------------
      DO NTC = NBTECR0+1, NBTECR-1

         NTE = NOTECR( NTC )

C        LES 4 FACES DU TETRAEDRE NTE SONT RECHERCHEES
C        DANS LES AUTRES TETRAEDRES DE LA SOUS-ETOILE
         DO 600 L=1,4

C           LA FACE L DU TETRAEDRE NTE
            IF( NOTETR(4+L,NTE) .LE. 0 ) THEN

C              FACE SANS TETRAEDRE OPPOSE
C              LES 3 SOMMETS DE LA FACE L DE NTE
               DO K=1,3
                  NOSOTR3(K) = NOTETR( NOSOFATE(K,L), NTE )
               ENDDO
               CALL TRI3NO( NOSOTR3, NOSOTR3 )

C              PARCOURS DES TETRAEDRES AU DELA DE NTC
               DO NTC1 = NTC+1, NBTECR

                  NTE1 = NOTECR( NTC1 )
C                 LA FACE L DE NTE NOSOTR3 EST ELLE UNE FACE LL DE NTE1
                  CALL NUFATRTE( NOSOTR3, NOTETR(1,NTE1), LL )

                  IF( LL .GT. 0 ) THEN
C                    LA FACE L  DU TETRAEDRE NTE  EST 
C                    LA FACE LL DU TETRAEDRE NTE1
                     NOTETR(4+L ,NTE ) = NTE1
                     NOTETR(4+LL,NTE1) = NTE
                     GOTO 600
                  ENDIF

               ENDDO

            ENDIF

 600     ENDDO

      ENDDO

C     MISE A JOUR DE LA SOUS-ETOILE APRES AJOUT DES TETRAEDRES
C     --------------------------------------------------------
C     NOMBRE DE FACES AJOUTEES A L'ETOILE PAR L'AJOUT DES TETRAEDRES
 1000 NBFAAJ = 0
C     NOMBRE DE FACES RETIREES DE L'ETOILE
      NBFARE = 0

      DO NTC = NBTECR0+1, NBTECR

C        NUMERO NOTETR DU TETRAEDRE NTC
         NTE = NOTECR( NTC )

C        LES 4 FACES DU TETRAEDRE NTE SONT AJOUTEES A L'ETOILE NFETOI
C         => 2 ou 3 ou 4 FACES DE NTE Y SONT EN FAIT RETIREES
C            1 ou 2 FACES PEUVENT LUI ETRE AJOUTEES
         DO K=1,4
C           LA FACE K DU TETRAEDRE NTE
            CALL AJFAET2( NTE, K, NOTETR, N1FEOC, N1FEVI, NFETOI,
     %                    NFA, NBFARE, NBFAAJ )
C           NFA : >0 NUMERO DE CETTE FACE AJOUTEE DANS NFETOI
C                 =0 SUPPRESSION DE LA FACE DE L'ETOILE
         ENDDO

      ENDDO

C     VERIFIER QUE TOUTE ARETE DE L'ETOILE APPARTIENT
C     SEULEMENT A 2n FACES DE L'ETOILE avec n=1,2
      CALL VETAFET( NOTETR, N1FEOC, NFETOI, NBFETO, NBARPB )

      PRINT*,'tetrtore:',NBFARE,' FACES RETIREES',
     %        NBFAAJ,' FACES AJOUTEES et l''ETOILE a maintenant',
     %        NBFETO,' FACES'

C     LE TITRE DU TRACE DES FACES DE LA SOUS-ETOILE ACTUELLE
C     TRACE EN NOIR DES 6 ARETES DU DERNIER TETRAEDRE CREE
      KTITRE='tetrtore:       TETRAEDRES CREES et        FACES DE L''ETO
     %ILE NOCAS=        '
C     AJOUT DANS LE TITRE DU NOMBRE DE FACES DE LA SOUS-ETOILE ET
C     DE TETRAEDRES CREES
      WRITE(KTITRE(10:14),'(I5)') NBTECR
      WRITE(KTITRE(37:41),'(I5)') NBFETO
      WRITE(KTITRE(69:70),'(I2)') NOCAS
      CALL SANSDBL( KTITRE, L )
      CALL TRFETO4( KTITRE(1:L), PTXYZD, N1FEOC, NFETOI,
     %              NBTRCF, NOTRCF, LEFACO, NO0FAR,
     %              NBTECR, NOTECR, NOTETR )

      IF( N1FEOC .GT. 0 ) THEN
C        TOUTE LA SOUS-ETOILE N'A PAS ENCORE ETE TETRAEDRISEE
C        => UN NOUVEAU TRAITEMENT DE LA SOUS-ETOILE
C        LA TETRAEDRISATION A T ELLE CREEE PLUSIEURS SOUS-SOUS-ETOILES?
         CALL SOUSETO( PTXYZD, NOTETR, N1FEOC, NFETOI,
     %                 MXSSET, NBSSET, N1SSET, IERR )
         IF( IERR .NE. 0 ) GOTO 9999
      ENDIF

      GOTO 5


C     MISE A JOUR DES TETRAEDRES OPPOSES DES FACES DES TETRAEDRES
C     SOMMET-FACES NFVPSI
C     -----------------------------------------------------------
 9000 CALL MJOPTE( NBTECR, NOTECR, N1TETS, NOTETR, MXTETR, 
     %             N1TEVI, PTXYZD, NBFANR )

C     TRACE FINAL DES TETRAEDRES CREES FORMANT L'ETOILE
C     -------------------------------------------------
      KTITRE = '      TETRAEDRES CREES en FIN de tetrtore'
      WRITE(KTITRE(1:5),'(I5)') NBTECR

 9010 CALL SANSDBL( KTITRE, L )
      CALL TRFETO4( KTITRE(1:L), PTXYZD, N1FEOC, NFETOI,
     %              NBTRCF, NOTRCF, LEFACO, NO0FAR,
     %              NBTECR, NOTECR, NOTETR )

 9999 PRINT *,'FIN tetrtore: IERR=',IERR,' A CREE',NBTECR,
     %        ' TETRAEDRES de VOLUME',VOLET1,' NBSSET=',NBSSET,
     %        ' NBFANR=',NBFANR,
     %        ' NBSOM0=',NBSOMM0,' NBSOMM=',NBSOMM

      RETURN
      END
