      SUBROUTINE TETRAJPA( ARETGR, NOFOTI, NBVOPA,
     %                     NBSOMM, MXSOMM, PTXYZD, NPSOFR,
     %                     MXTETR, N1TEVI, NOTETR, NUDTETR,
     %                     N1TETS, NOF1VO, NVOLTE,
     %                     MXFACO, LEFACO, MXPILE, LAPILE,
     %                     IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   TETRAEDRISATION DE POINTS SUR L'ARETE MAX DES TETRAEDRES
C -----   SI SA TAILLE EST SUPERIEURE A LA TAILLE SOUHAITEE EN AU
C         MOINS UN DES SOMMETS
C         ESSAI D'AMELIORER LES TETRAEDRES AU DELA DES DECOUPES
C         PAR 2T->3T ou 3T->2T

C ENTREES:
C --------
C ARETGR : LONGUEUR DE L'ARETE SOHAITEE POUR L'ENSEMBLE DES TETRAEDRES
C NOFOTI : NUMERO DE LA FONCTION 'TAILLE_IDEALE' DES ARETES SINON 0
C NBVOPA : NOMBRE DE VOLUMES DE LA PARTITION DE VOLUMES
C MXSOMM : NOMBRE MAXIMAL DE SOMMETS DECLARABLES
C MXTETR : NOMBRE MAXIMAL DE TETRAEDRES DECLARABLES
C MXPILE : MAXIMUM D'ELEMENTS DE LAPILE(2,MXPILE)

C ENTREE ET SORTIE :
C ------------------
C NBSOMM : NOMBRE DES (FUTURS) SOMMETS (DE LA TETRAEDRISATION
C          EN ENTREE ET EN SORTIE APRES AJOUT DES POINTS SUR LES ARETES
C PTXYZD : X Y Z DISTANCE SOUHAITEE DES POINTS
C NPSOFR : NUMERO DE POSITION PAR RAPPORT A LA FRONTIERE DE CHAQUE POINT
C          LE SIGNE DEVIENT NEGATIF SI LE SOMMET EST DEPLACE
C          =  0 SI POINT AJOUTE NON SUR LA SURFACE DE L'OBJET
C          =  1 SI LE POINT EST AJOUTE SUR UNE ARETE DE LA SURFACE
C          =  2 SI LE POINT EST AJOUTE COMME BARYCENTRE SUR LA SURFACE
C          =  1 000 000  + NUMERO DU  POINT INTERNE UTILISATEUR
C          = (1 000 000) * (NO SURFACE + 1) + NUMERO DU SOMMET DANS
C                    LA NUMEROTATION DES SOMMETS DE LA SURFACE
C          = -1 SI LE POINT EST SOMMET RECONNU TROP PROCHE
C          = -4 SI LE POINT EST SOMMET NON TROP PROCHE
C          = -3 SI LE POINT A ETE REFUSE DANS LA TETRAEDRISATION
C          = -NPSOFR(I) SI POINT I DEPLACE SUR LA SURFACE
C             DANS LEUR SURFACE FERMEE INITIALE OU NO DE POINT INTERNE
C N1TEVI : NUMERO DU PREMIER TETRAEDRE VIDE DANS LE CHAINAGE NOTETR(5,*)
C NOTETR : LISTE DES TETRAEDRES
C          SOMMET1,    SOMMET2,    SOMMET3,    SOMMET4,
C          TETRAEDRE1, TETRAEDRE2, TETRAEDRE3, TETRAEDRE4
C          DE L'AUTRE COTE DE LA FACE
C          1: 123      2: 234      3: 341      4: 412
C NUDTETR: PLUS GRAND NUMERO DANS NOTETR DE TETRAEDRE OCCUPE
C N1TETS : NUMERO DANS NOTETR D'UN TETRAEDRE DE CHAQUE SOMMET

C NOF1VO : TABLEAU ENTIER( 1:3 , 1:NBVOPA )
C          NOF1VO(1,NV) = NUMERO DU 1-ER TETRAEDRE DE CHAQUE VOLUME
C          NOF1VO(3,NV) = NUMERO DU VOLUME DE 1 A NBVOPA
C NVOLTE : >0  NUMERO DU VOLUME DE 1 A NBVOPA DE CHAQUE TETRAEDRE
C          =<0 TETRAEDRE INACTIF ou SUPPRIME

C MXFACO : NOMBRE MAXIMAL DECLARABLE DE FACES DU CONTOUR
C LEFACO : FACE DU CONTOUR OU INTERFACES ENTRE VOLUMES
C          IL CONTIENT DANS CET ORDRE
C          NUMERO (DANS PTXYZD) DU SOMMET 1, SOMMET 2, SOMMET 3
C          NUMERO (DANS NUVOPA 0 SINON) DU VOLUME1 , VOLUME2 DE LA FACE
C          NUMERO (DANS LEFACO) DE LA FACE ADJACENTE PAR L'ARETE 1 2 3
C          ATTENTION: UNE ARETE PEUT APPARTENIR A PLUS DE 2 FACES
C          => CHAINAGE CIRCULAIRE DE CES FACES DANS LEFACO
C          LEFACO(9,*)  -> FACE SUIVANTE (*=0:VIDE, *<>0:NON VIDE)
C          HACHAGE AVEC LA SOMME DES 3 SOMMETS MODULO MXFACO
C          LF = MOD( NOSOFA(1)+NOSOFA(2)+NOSOFA(3) , MXFACO ) + 1
C          NF = LEFACO( 10, LF ) LE NUMERO DE LA 1-ERE FACE DANS LEFACO
C          SI LA FACE NE CONVIENT PAS. PASSAGE A LA SUIVANTE
C          NF  = LEFACO( 9, NF )  ...
C          11: >0  NO NOTETR D'UN TETRAEDRE AYANT CETTE FACE,
C              =0  SINON
cccC          12: = NO FACEOC DE 1 A NBFACES D'OC

C LAPILE : PILE AUXILIAIRE

C SORTIE :
C --------
C IERR   : 0 SI PAS D'ERREUR
C          1 SATURATION DES SOMMETS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET Saint PIERRE du PERRAY          Septembre 2018
C23456...............................................................012
      PARAMETER        (MXTEOLD=2048, MXTENEW=2*MXTEOLD)
      PARAMETER        (MXTE1A=1024)
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      COMMON / TRTETR / STOPTE, TRACTE
      LOGICAL           STOPTE, TRACTE
      CHARACTER*80      KTITRE

      DOUBLE PRECISION  PTXYZD(4,MXSOMM),
     %                  AREMAX, V, VOLUMT, D, TI1, TI2,
     %                  LONARE, VOLTET, DBLE, C1, C2, VOLDTA, VOLMOY
      INTEGER           NPSOFR(MXSOMM), LEFACO(1:11,0:MXFACO),
     %                  NOTETR(8,MXTETR), N1TETS(MXSOMM), 
     %                  NOF1VO(3,NBVOPA),
     %                  LAPILE(2,MXPILE), NVOLTE(1:*), NOTE1A(MXTE1A),
     %                  NTEOLD(MXTEOLD), NTENEW(MXTENEW),
     %                  NTEN23(8), NOSOTR(3),
     %                  NTPARE(-1:62), NPTARE(0:63)
      INTRINSIC         MAX

C     NUMERO LOCAL DES 2 SOMMETS DES 6 ARETES DU TETRAEDRE
      INTEGER           NOSOARTE(2,6)
      DATA              NOSOARTE/ 1,2, 2,3, 3,1, 1,4, 2,4, 3,4 /
C     NUMERO LOCAL DE L'ARETE OPPOSEE A UNE ARETE DANS UN TETRAEDRE
      INTEGER           NOAROPTE(6)
      DATA              NOAROPTE/ 6, 4, 5, 2, 3, 1 /

C     LES SOMMETS SONT VUS A PARTIR DU BARYCENTRE DANS LE SENS DIRECT
      INTEGER           NOSOFATE(3,4)
      DATA              NOSOFATE/ 1,2,3, 2,4,3, 3,4,1, 4,2,1 /
C     NO (DE 1 A 9) DES 3 ARETES DES 4 FACES DU TETRAEDRE
      INTEGER           NOARFATE(3,4)
      DATA              NOARFATE / 1,2,3, 2,5,6, 3,6,4, 5,1,4 /

      WRITE(IMPRIM,*)
      WRITE(IMPRIM,*)'tetrajpa: MISE A JOUR de la TAILLE des ARETES par
     %TETRAEDRISATION de POINTS sur l''ARETE MAX TROP LONGUE des TETRAED
     %RES'

      IVOLTE = 1
      NUDTET0= NUDTETR
      VOLDTA = 0D0
      NBT    = 0
      NB2T3T = 0
      NB3T2T = 0
      NBSOM00 = NBSOMM
      ITERAT  = 0

 10   ITERAT = ITERAT + 1

C     BOUCLE SUR LES NBVOPA VOLUMES DE LA PARTITION
C     =============================================
      NBSOM0 = NBSOMM
      VOLUMT = 0D0

      DO 200 NV = 1, NBVOPA

C        LE NO DE VOLUME DE 1 A NBVOPA DU TETRAEDRE NTE ET SUIVANTS
         NOVOLU = NOF1VO( 3, NV )

C        BOUCLE SUR LES TETRAEDRES CHAINES DU VOLUME NOVOLU
C        --------------------------------------------------
         DO 100 NTEV = 1, NUDTET0

            IF( NOTETR(1,NTEV) .EQ. 0 .OR. NVOLTE(NTEV) .NE. NOVOLU )
     %         GOTO 100

C           MISE A JOUR DE LA TAILLE SOUHAITEE DES ARETES AUX 4 SOMMETS DE NTEV
            DO K=1,4

C              NUMERO DU K-EME SOMMET DE NTEV
               NS = NOTETR( K, NTEV )

               IF( NOFOTI .GT. 0 ) THEN
C                 LA FONCTION TAILLE_IDEALE(X,Y,Z) EXISTE
C                 MISE A JOUR DE LA TAILLE IDEALE AUX 4 SOMMETS DE NTEV
C                 CALCUL DE TAILLE_IDEALE( NS )
                  CALL FONVAL( NOFOTI ,3, PTXYZD(1,NS), NCODEV,
     %                         PTXYZD(4,NS) )
                  IF( NCODEV .LE. 0 ) THEN
                     PTXYZD(4,NS) = ARETGR
                  ENDIF
               ENDIF

C              INTERDICTION D'UNE VALEUR NEGATIVE OU NULLE AU SOMMET NS
               IF( PTXYZD(4,NS).LE.0D0 .OR. PTXYZD(4,NS).GT.ARETGR )THEN
                  PTXYZD(4,NS) = ARETGR
               ENDIF

            ENDDO

C           LES COORDONNEES DU CENTRE DE LA BOULE CIRCONSCRITE
C           AU TETRAEDRE NTEV ET CARRE DE SON RAYON
            CALL QV1TET( PTXYZD, NOTETR(1,NTEV), V, Q )
            NBT = NBT + 1
            VOLDTA = VOLDTA + V
C           VOLUME MOYEN DES TETRAEDRES TRAITES
            VOLMOY = VOLDTA / NBT

ccc            IF( Q .LT. 0.001 .OR. V .LT. VOLMOY*1D-3 ) THEN
ccc            IF( Q .LT. 0.001 .OR. V .LT. VOLMOY*1D-4 ) THEN
            IF( Q .LT. 0.005 .OR. V .LT. VOLMOY*1D-4 ) THEN
C              => PAS D'AJOUT D'UN POINT SUR L'ARETE MAX
               GOTO 100
            ENDIF

C           RECHERCHE DE L'ARETE DE TAILLE MAXIMALE DU TETRAEDRE NTEV
            NAREMX = 0
            AREMAX = 0D0
            DO 31 NA = 1, 6

C              NUMERO DES 2 SOMMETS DE L'ARETE NA
               NS1 = NOTETR( NOSOARTE(1,NA), NTEV )
               NS2 = NOTETR( NOSOARTE(2,NA), NTEV )

C              L'ARETE NS1-NS2 EST ELLE SUR LA FRONTIERE?
               IF( NPSOFR(NS1) .GT. 0 .AND. NPSOFR(NS2) .GT. 0 ) THEN
C                 ARETE MAX FRONTALIERE
                  GOTO 31
               ENDIF

               LONARE = ( PTXYZD(1,NS2) - PTXYZD(1,NS1) ) ** 2
     +             + ( PTXYZD(2,NS2) - PTXYZD(2,NS1) ) ** 2
     +             + ( PTXYZD(3,NS2) - PTXYZD(3,NS1) ) ** 2
               IF( LONARE .GT. AREMAX ) THEN
                  NAREMX = NA
                  AREMAX = LONARE
               ENDIF

 31         ENDDO

            IF( NAREMX .EQ. 0 ) THEN
C              LES 4 SOMMETS DU TETRAEDRE SONT SUR LA FRONTIERE
C              => PAS D'AJOUT DU POINT SUR L'ARETE
               GOTO 100
            ENDIF
            AREMAX = SQRT( AREMAX )

C           TAILLE SOUHAITEE AUX 2 SOMMETS DE L'ARETE MAX
            NS1 = NOTETR( NOSOARTE(1,NAREMX), NTEV )
            NS2 = NOTETR( NOSOARTE(2,NAREMX), NTEV )
            TI1 = PTXYZD( 4, NS1 )
            TI2 = PTXYZD( 4, NS2 )

C           NBPARE NOMBRE DE POINTS A AJOUTER ENTRE NS1 et NS2
            IF( AREMAX .LE. (TI1+TI2) * 0.65D0 ) THEN
               NBPARE = 0
C              => PAS D'AJOUT DE POINT SUR L'ARETE NS1-NS2
               GOTO 100
            ENDIF
            IF( AREMAX .LE. (TI1+TI2) * 1.3D0 ) THEN
               NBPARE = 1
            ELSE
               NBPARE = NINT( AREMAX / ( TI1 + TI2 ) )
            ENDIF

C           LISTAGE DES NBTEOLD TETRAEDRES D'ARETE NS1 NS2
            CALL TETR1A( NS1,     NS2,     N1TETS, NOTETR,
     %                   NBTEOLD, MXTEOLD, NTEOLD, IERR )

            IF( IERR .NE. 0 .OR. NBTEOLD .GT. 24 ) THEN

C              TRACE DES TETRAEDRES NTEOLD
               TRACTE = .TRUE.
               KTITRE='tetrajpa:       TETRAEDRES AUTOUR DE L''ARETE
     %           '
               WRITE( KTITRE(11:15),'(I5)') NBTEOLD
               WRITE( KTITRE(47:53),'(I7)') NS1
               WRITE( KTITRE(55:61),'(I7)') NS2
               NTEN23(1) = NTEV
               NTEN23(2) = 1
               NTEN23(3) = NTEV
               NTEN23(4) = 2
               NTEN23(5) = NTEV
               NTEN23(6) = 3
               NTEN23(7) = NTEV
               NTEN23(8) = 4
               CALL TRFETO12( KTITRE,  PTXYZD,
     %                        NBTEOLD, NTEOLD, NOTETR,
     %                        4,  NTEN23 )
               TRACTE = .FALSE.

C              => PAS D'AJOUT D'UN POINT SUR L'ARETE COMMUNE A
C                 TROP DE TETRAEDRES
               GOTO 100

            ENDIF


C           AJOUT DE NBPARE POINTS SUR L'ARETE MAXIMALE DU TETRAEDRE NTEV
C           ---------------------------------------------------------------
C           CALCUL DES XYZ DES NBPARE POINTS SUR L'ARETE MAX DE NTEV
            IF( NBSOMM+NBPARE .GT. MXSOMM ) THEN
               PRINT *,'tetrajpa: TABLEAU PTXYZD SATURE MXSOMM=',MXSOMM
               GOTO 400
            ENDIF

            DO NP = 1, NBPARE

               NBSOMM = NBSOMM + 1
               NPTARE( NP ) = NBSOMM
               IF( MOD(NBSOMM,1000) .EQ. 0 ) THEN
                  IF( LANGAG .EQ. 0 ) THEN
                    WRITE(IMPRIM,*)'tetrajpa: TETRAEDRISATION du POINT',
     %              NBSOMM,' de l''ARETE',NS1,NS2,
     %              ' du TETRAEDRE',NTEV,' VOLMOY=',VOLMOY
                  ELSE
                   WRITE(IMPRIM,*)'tetrajpa: TETRAHEDRIZATION of POINT',
     %             NBSOMM,' of EDGE',NS1,NS2,
     %             ' of TETRAHEDRON',NTEV,' VOLMOY=',VOLMOY
                  ENDIF
               ENDIF

C              NBSOMM NON FRONTALIER EST AJOUTE SUR UNE ARETE INTERNE
               NPSOFR( NBSOMM ) = 0

C              LES 3 COORDONNEES ET DISTANCE SOUHAITEE AU POINT SUR L'ARETE
               IF( NBPARE .EQ. 1 ) THEN
                  D  = TI1 + TI2
                  C1 = TI2 / D
                  C2 = TI1 / D
               ELSE
                  D  = DBLE( NP ) / DBLE( NBPARE+1 )
                  C1 = 1D0 - D
                  C2 = D
               ENDIF

               DO K=1,4
                  PTXYZD(K,NBSOMM) = C1*PTXYZD(K,NS1) + C2*PTXYZD(K,NS2)
               ENDDO

C              CALCUL DE LA TAILLE_IDEALE( NBSOMM )
               IF( NOFOTI .GT. 0 ) THEN
C                 LA FONCTION TAILLE_IDEALE(X,Y,Z) EXISTE
                  CALL FONVAL( NOFOTI ,3, PTXYZD(1,NBSOMM), NCODEV,
     %                         PTXYZD(4,NBSOMM) )
                  IF( NCODEV .LE. 0 ) THEN
                     PTXYZD(4,NBSOMM) = ARETGR
                  ENDIF
               ENDIF
               IF( PTXYZD(4,NS).LE.0D0 .OR. PTXYZD(4,NS).GT.ARETGR )THEN
                  PTXYZD(4,NBSOMM) = ARETGR
               ENDIF

            ENDDO
C           LE NO DE POINT DES 2 EXTREMITES DE L'ARETE
            NPTARE( 0 ) = NS1
            NPTARE( NBPARE+1 ) = NS2

C           DECOUPAGE EN NBPARE+1 SOUS-TETRAEDRES DE CHACUN DES NBTEOLD
C           TETRAEDRES PAR AJOUT DES NBPARE POINTS SUR L'ARETE NS1 NS2
            NBTENEW = 0
            IF( NBPARE+1 .GT. MXTENEW ) THEN
               PRINT*,'tetrajpa: TABLEAU NTENEW SATURE MXTENEW=',MXTENEW
               NBSOMM = NBSOMM - NBPARE
               GOTO 100
            ENDIF

C           INITIALISATION DE LA PILE AVEC LA FACE EXTERNE DES
C           4 TETRAEDRES NTENEW
C           --------------------------------------------------
            NS4    = 0
            LHPILE = 0
            DO N = 1, NBTEOLD

               NTE = NTEOLD( N )

C              QUELS SONT LES NUMEROS DES 2 AUTRES SOMMETS DE NTE?
               DO M = 1, 4
                  NS3 = NOTETR( M, NTE )
                  IF( NS3 .NE. NS1 .AND. NS3 .NE. NS2 ) GOTO 33
               ENDDO

 33            DO MM = M+1, 4
                  NS4 = NOTETR( MM, NTE )
                  IF(NS4 .NE. NS1 .AND. NS4 .NE. NS2 .AND. NS4 .NE. NS3)
     %                 GOTO 34
               ENDDO

C              VOLUME DU TETRAEDRE NTE
 34            V = VOLTET( PTXYZD( 1, NS1 ), PTXYZD( 1, NS2 ),
     %                     PTXYZD( 1, NS3 ), PTXYZD( 1, NS4 ) )
               IF( V .LT. 0D0 ) THEN
C                 VOLUME RENDU POSITIF PAR PERMUTATION NS3-NS4
                  M   = NS3
                  NS3 = NS4
                  NS4 = M
               ENDIF

               DO NP = 0, NBPARE
C                 RESERVATION DU SOUS-TETRAEDRE DE SOMMET NPTARE(NP)
                  IF( N1TEVI .LE. 0 ) THEN
                     PRINT *,'tetrajpa: TABLEAU NOTETR SATURE'
                     GOTO 400
                  ENDIF
C                 LE TETRAEDRE A CREER EST LE PREMIER TETRAEDRE VIDE
                  NBTENEW = NBTENEW + 1
                  NT1     = N1TEVI
                  NTENEW( NBTENEW ) = NT1
C                 MISE A JOUR DU 1-ER TETRAEDRE VIDE
                  N1TEVI = NOTETR( 5, N1TEVI )
C                 NUMERO MAX DES TETRAEDRES UTILISES DANS NOTETR
                  NUDTETR = MAX( NUDTETR, NT1 )
                  NTPARE( NP ) = NT1
               ENDDO

C              LE NUMERO LOCAL NF DE LA FACE DE SOMMETS NS3-NS4-NS1
C              PARMI LES 4 SOMMETS DU TETRAEDRE NTE
               NOSOTR(1) = NS3
               NOSOTR(2) = NS4
               NOSOTR(3) = NS1
               CALL TRI3NO( NOSOTR, NOSOTR )
               CALL NO1F1T( NOSOTR, NOTETR(1,NTE), NF )
C              NO NOTETR DU TETRAEDRE OPPOSE A LA FACE NF DE NTE
               NTOP0 = NOTETR( 4+NF, NTE )
               NTPARE(-1) = NTOP0

               IF( NTOP0 .GT. 0 ) THEN
                  CALL NO1F1T( NOSOTR, NOTETR(1,NTOP0), NFOP )
                  IF( NFOP .GT. 0 ) THEN
C                    NTPARE(0) EST L'OPPOSE DE LA FACE NFOP DE NTOP0
                     NOTETR( 4+NFOP, NTOP0 ) = NTPARE( 0 )
C                    LE TETRAEDRE NTOP0 EST IL DANS LE VOLUME NOVOLU?
                     IF( NVOLTE( NTOP0 ) .EQ. NOVOLU ) THEN
cccC                    LA FACE NFOP DE NTOP0 EST ELLE DANS LEFACO?
ccc                     CALL NULEFT( NFOP, NTOP0, NOTETR, MXFACO, LEFACO,
ccc     %                            NFLE )
ccc                     IF( NFLE .LE. 0 ) THEN
C                       LA FACE N'EST PAS DANS LEFACO. ELLE EST EMPILEE
                        LHPILE = LHPILE + 1
                        LAPILE(1,LHPILE) = NTOP0
                        LAPILE(2,LHPILE) = NFOP
ccc                     ENDIF
                     ENDIF
                  ENDIF
               ENDIF

C              LE NUMERO LOCAL NF DE LA FACE EXTERIEURE DE SOMMETS NS2-NS4-NS3
C              PARMI LES 4 SOMMETS DU TETRAEDRE NTE
               NOSOTR(1) = NS2
               NOSOTR(2) = NS4
               NOSOTR(3) = NS3
               CALL TRI3NO( NOSOTR, NOSOTR )
               CALL NO1F1T( NOSOTR, NOTETR(1,NTE), NF )
C              NO NOTETR DU TETRAEDRE OPPOSE A LA FACE NF DE NTE
               NTOP2 = NOTETR( 4+NF, NTE )
               NTPARE(NBPARE+1) = NTOP2

               IF( NTOP2 .GT. 0 ) THEN
                  CALL NO1F1T( NOSOTR, NOTETR(1,NTOP2), NFOP )
                  IF( NFOP .GT. 0 ) THEN
C                    NTPARE(NBPARE) EST L'OPPOSE DE LA FACE NFOP DE NTOP2
                     NOTETR( 4+NFOP, NTOP2 ) = NTPARE( NBPARE )
C                    LE TETRAEDRE NTOP2 EST ILDANS LE VOLUME NOVOLU?
                     IF( NVOLTE( NTOP2 ) .EQ. NOVOLU ) THEN
cccC                    LA FACE NFOP DE NTOP2 EST ELLE DANS LEFACO?
ccc                     CALL NULEFT( NFOP, NTOP2, NOTETR, MXFACO, LEFACO,
ccc     %                            NFLE )
ccc                     IF( NFLE .LE. 0 ) THEN
C                       LA FACE N'EST PAS DANS LEFACO. ELLE EST EMPILEE
                        LHPILE = LHPILE + 1
                        LAPILE(1,LHPILE) = NTOP2
                        LAPILE(2,LHPILE) = NFOP
ccc                     ENDIF
                     ENDIF
                  ENDIF
               ENDIF

C              LA FACE NS1 NS2 NS3 DE NTE EST ELLE SUR LA FRONTIERE?
               NOSOTR(1) = NS1
               NOSOTR(2) = NS2
               NOSOTR(3) = NS3
               CALL TRI3NO( NOSOTR, NOSOTR )
               CALL NO1F1T( NOSOTR, NOTETR(1,NTE), NF )
C              NO NOTETR DU TETRAEDRE OPPOSE A LA FACE NF DE NTE
               NTOP123 = NOTETR( 4+NF, NTE )

C              LA FACE NS1 NS2 NS4 DE NTE EST ELLE SUR LA FRONTIERE?
               NOSOTR(1) = NS1
               NOSOTR(2) = NS2
               NOSOTR(3) = NS4
               CALL TRI3NO( NOSOTR, NOSOTR )
               CALL NO1F1T( NOSOTR, NOTETR(1,NTE), NF )
C              NO NOTETR DU TETRAEDRE OPPOSE A LA FACE NF DE NTE
               NTOP124 = NOTETR( 4+NF, NTE )

               DO NP = 0, NBPARE

C                 INITIALISATION DU SOUS-TETRAEDRE DE SOMMET NPTARE(NP)
                  NT1 = NTPARE( NP )

C                 SOMMETS DU TETRAEDRE NT1
                  NOTETR( 1, NT1 ) = NPTARE( NP   )
                  NOTETR( 2, NT1 ) = NPTARE( NP+1 )
                  NOTETR( 3, NT1 ) = NS3
                  NOTETR( 4, NT1 ) = NS4

C                 LES TETRAEDRES OPPOSES AUX 4 FACES DE NT1
                  IF( NTOP123 .EQ. 0 ) THEN
C                    FACE FRONTALIERE
                     NOTETR( 5, NT1 ) = 0
                  ELSE
                     NOTETR( 5, NT1 ) = -1
                  ENDIF

                  NOTETR( 6, NT1 ) = NTPARE( NP+1 )
                  NOTETR( 7, NT1 ) = NTPARE( NP-1 )

                  IF( NTOP124 .EQ. 0 ) THEN
C                    FACE FRONTALIERE
                     NOTETR( 8, NT1 ) = 0
                  ELSE
                     NOTETR( 8, NT1 ) = -1
                  ENDIF

C                 LES COORDONNEES DU CENTRE DE LA BOULE CIRCONSCRITE
C                 AU TETRAEDRE NT1 ET CARRE DE SON RAYON
                  CALL QV1TET( PTXYZD, NOTETR(1,NT1), V, Q )

C                 MISE A JOUR D'UN TETRAEDRE CONTENANT CHAQUE SOMMET
                  N1TETS( NPTARE(NP  ) ) = NT1
                  N1TETS( NPTARE(NP+1) ) = NT1
                  N1TETS( NS3 ) = NT1
                  N1TETS( NS4 ) = NT1

C                 AJOUT DE NTE AU VOLUME NOVOLU
                  NVOLTE( NT1 ) = NOVOLU

               ENDDO

            ENDDO


C           COMPLETION DES TETRAEDRES OPPOSES DES TETRAEDRES CREES NTENEW
C           -------------------------------------------------------------
            CALL MJOPTE( NBTENEW, NTENEW, N1TETS, NOTETR, NUDTETR,
     %                   N1TEVI,  PTXYZD, NBFANR )
            IF( NBFANR .GT. 0 ) THEN
               PRINT*,'tetrajpa: Probleme sur les TETRAEDRES OPPOSES aux
     % TETRAEDRES NTENEW'
            ENDIF

C           SUPPRESSION DES TETRAEDRES NTEOLD DU TABLEAU NOTETR ET NVOLTE
C           -------------------------------------------------------------
            DO N = 1, NBTEOLD
               NTE = NTEOLD( N )
C              AJOUT DU TETRAEDRE NTE A LA CHAINE DES TETRAEDRES VIDES
               NOTETR( 1, NTE ) = 0
               NOTETR( 5, NTE ) = N1TEVI
               N1TEVI = NTE
C              SUPPRESSION DE NTE DU VOLUME NOVOLU
               NVOLTE( NTE ) = -1
            ENDDO


C        ESSAI D'AMELIORER LES 4 TETRAEDRES OPPOSES AUX 4 FACES
C        PAR 2T->3T ou 3T->2T et EMPILAGE DES FACES CREEES
C        ======================================================
         NBEC2T3T = 0
         NBEC3T2T = 0

C        TRAITEMENT DU HAUT DE LA PILE DES FACES OPPOSEES
C        ------------------------------------------------
 50      IF( LHPILE .GT. 0 ) THEN

C           RECUPERATION DU TETRAEDRE NTE ET DE SA FACE NFE A TRAITER
            NTE = LAPILE( 1, LHPILE )
            NFE = LAPILE( 2, LHPILE )
            LHPILE = LHPILE - 1

C           LE TETRAEDRE NTE EXISTE-T-IL ENCORE?
            IF( NTE .LE. 0 .OR. NFE .LE. 0 ) GOTO 50
            IF( NOTETR(1,NTE) .LE. 0 ) GOTO 50

C           LE TETRAEDRE NTE EST IL DANS LE VOLUME NOVOLU?
            IF( NVOLTE( NTE ) .NE. NOVOLU ) GOTO 50

cccC           LA FACE NFE DE NTE EST ELLE DANS LEFACO?
ccc            CALL NULEFT( NFE, NTE, NOTETR, MXFACO, LEFACO, NFLE )
ccc            IF( NFLE .GT. 0 ) THEN
cccC              OUI: LA FACE EST DANS LEFACO
cccC              CONSERVATION DE CETTE FACE LEFACO EN NE PASSANT PAS AU DELA
ccc               GOTO 50
ccc            ENDIF

C           TENTATIVE DE CHANGER 2T->3T SI MINIMUM QUALITE EST MAXIMISE
C           -----------------------------------------------------------
C           LE NUMERO NFOP DE LA FACE COMMUNE DANS LE TETRAEDRE OPPOSE NTOP
            CALL NOFAOP( NFE, NTE, NOTETR, NFOP, NTOP )
            IF( NFOP .LE. 0 ) GOTO 50

C           LE TETRAEDRE NTOP EST IL DANS LE VOLUME NOVOLU?
            IF( NVOLTE( NTOP ) .NE. NOVOLU ) GOTO 50

            CALL CH2T3T( 1,      MXFACO, LEFACO, 1, IVOLTE, NVOLTE,
     %                   NTE,    NFE,    NTOP,
     %                   PTXYZD, N1TETS, NOTETR, N1TEVI, NUDTETR,
     %                   NTEN23, IER )

            IF( IER .EQ. 0 ) THEN

C              REUSSITE: CHANGEMENT 2T->3T EFFECTUE
               NBEC2T3T = NBEC2T3T + 1
               NBTEN23 = 3

C              NTEN23=NO NOTETR DES 3 TETRAEDRES QUI REMPLACENT NTE et NTOP
C              LES 2 FACES EXTERIEURES 1 et 4 DE CES 3 TETRAEDRES SONT EMPILEES
C              CF ch2t3t.f
               IF( LHPILE+6 .GT. MXPILE ) GOTO 9990
               DO K = 1, NBTEN23

C                 LE NOUVEAU TETRAEDRE K
                  NT = NTEN23( K )

C                 RECHERCHE DU TETRAEDRE OPPOSE NTOP ET DE SA FACE NFOP
C                 A LA FACE 1 DU TETRAEDRE NT
                  CALL NOFAOP( 1, NT, NOTETR,  NFOP, NTOP )
                  IF( NFOP .GT. 0 ) THEN
C                    LE TETRAEDRE NTOP EST IL DANS LE VOLUME NOVOLU?
                     IF( NVOLTE( NTOP ) .EQ. NOVOLU ) THEN
cccC                       LA FACE NFOP DE NTOP EST ELLE DANS LEFACO?
ccc                        CALL NULEFT( NFOP,NTOP,NOTETR,MXFACO,LEFACO,
ccc     %                               NFLE)
ccc                        IF( NFLE .LE. 0 ) THEN
C                          LA FACE N'EST PAS DANS LEFACO
                        LHPILE = LHPILE + 1
                        LAPILE(1,LHPILE) = NTOP
                        LAPILE(2,LHPILE) = NFOP
ccc                        ENDIF
                     ENDIF
                  ENDIF

C                 RECHERCHE DU TETRAEDRE OPPOSE NTOP ET DE SA FACE NFOP
C                 A LA FACE 4 DU TETRAEDRE NT
                  CALL NOFAOP( 4, NT, NOTETR,  NFOP, NTOP )
                  IF( NFOP .GT. 0 ) THEN
C                    LE TETRAEDRE NTOP EST IL DANS LE VOLUME NOVOLU?
                     IF( NVOLTE( NTOP ) .EQ. NOVOLU ) THEN
cccC                       LA FACE NFOP DE NTOP EST ELLE DANS LEFACO?
ccc                        CALL NULEFT( NFOP,NTOP,NOTETR,MXFACO,LEFACO,
ccc     %                               NFLE )
ccc                        IF( NFLE .LE. 0 ) THEN
C                          LA FACE N'EST PAS DANS LEFACO
                        LHPILE = LHPILE + 1
                        LAPILE(1,LHPILE) = NTOP
                        LAPILE(2,LHPILE) = NFOP
ccc                        ENDIF
                     ENDIF
                  ENDIF

               ENDDO

               GOTO 80

            ENDIF


C           TENTATIVE DE CHANGER 3T->2T SI MIN QUALITE EST MAXIMISE
C           -------------------------------------------------------
C           PARCOURS DES 3 ARETES DE LA FACE NFE DU TETRAEDRE NTE
            DO 70 N = 1, 3

C              NO DE L'ARETE N DE LA FACE NFE
               NAR = NOARFATE( N, NFE )

C              NUMERO PTXYZD DES 2 SOMMETS DE L'ARETE NAR
               NST1 = NOTETR( NOSOARTE( 1, NAR ), NTE )
               NST2 = NOTETR( NOSOARTE( 2, NAR ), NTE )

C              L'ARETE NST1-NST2 PERFORE LA FACE COMMUNE AUX 2 TETRAEDRES?
C              2 TETRAEDRES => 3 TETRAEDRES D'ARETE COMMUNE NST1-NST2
               CALL TETR1A( NST1,   NST2,   N1TETS, NOTETR,
     %                      NBTE1A, MXTE1A, NOTE1A, IERR )

               IF( NBTE1A .NE. 3 ) GOTO 70

               DO M=1,3
                  IF( NVOLTE( NOTE1A(M) ) .NE. NOVOLU ) GOTO 70
               ENDDO

               CALL CH3T2T( 1,      MXFACO, LEFACO, 1, 0, NVOLTE,
     %                      NTE,    NAR,    PTXYZD,
     %                      N1TETS, NOTETR, N1TEVI, NUDTETR,
     %                      NTEN23, IER )

               IF( IER .EQ. 0 ) THEN

C                 REUSSITE: CHANGEMENT 3T->2T EFFECTUE
                  NBEC3T2T = NBEC3T2T + 1
                  NBTEN23 = 2

C                 LES 3 FACES EXTERIEURES (2 3 4) DE CES 2 TETRAEDRES SONT EMPILEES
C                 CF ch3t2t.f
                  IF( LHPILE+6 .GT. MXPILE ) GOTO 9990

                  DO K = 1, NBTEN23

C                    LE NOUVEAU TETRAEDRE K
                     NT = NTEN23( K )

C                    RECHERCHE DU TETRAEDRE OPPOSE NTOP ET DE SA FACE NFOP
C                    A LA FACE 2 DU TETRAEDRE NT
                     CALL NOFAOP( 2, NT, NOTETR,  NFOP, NTOP )
                     IF( NFOP .GT. 0 ) THEN
C                       LE TETRAEDRE NTOP EST IL DANS LE VOLUME NOVOLU?
                        IF( NVOLTE( NTOP ) .EQ. NOVOLU ) THEN
cccC                          LA FACE NFOP DE NTOP EST ELLE DANS LEFACO?
ccc                           CALL NULEFT( NFOP,NTOP,NOTETR,MXFACO,LEFACO,
ccc     %                                  NFLE )
ccc                           IF( NFLE .LE. 0 ) THEN
cccC                             LA FACE N'EST PAS DANS LEFACO
                           LHPILE = LHPILE + 1
                           LAPILE(1,LHPILE) = NTOP
                           LAPILE(2,LHPILE) = NFOP
ccc                           ENDIF
                        ENDIF
                     ENDIF

C                    RECHERCHE DU TETRAEDRE OPPOSE NTOP ET DE SA FACE NFOP
C                    A LA FACE 3 DU TETRAEDRE NT
                     CALL NOFAOP( 3, NT, NOTETR,  NFOP, NTOP )
                     IF( NFOP .GT. 0 ) THEN
C                       LE TETRAEDRE NTOP EST IL DANS LE VOLUME NOVOLU?
                        IF( NVOLTE( NTOP ) .EQ. NOVOLU ) THEN
cccC                          LA FACE NFOP DE NTOP EST ELLE DANS LEFACO?
ccc                           CALL NULEFT( NFOP,NTOP,NOTETR,MXFACO,LEFACO,
ccc     %                                  NFLE)
ccc                           IF( NFLE .LE. 0 ) THEN
cccC                             LA FACE N'EST PAS DANS LEFACO
                           LHPILE = LHPILE + 1
                           LAPILE(1,LHPILE) = NTOP
                           LAPILE(2,LHPILE) = NFOP
ccc                           ENDIF
                        ENDIF
                     ENDIF

C                    RECHERCHE DU TETRAEDRE OPPOSE NTOP ET DE SA FACE NFOP
C                    A LA FACE 4 DU TETRAEDRE NT
                     CALL NOFAOP( 4, NT, NOTETR,  NFOP, NTOP )
                     IF( NFOP .GT. 0 ) THEN
C                       LE TETRAEDRE NTOP EST IL DANS LE VOLUME NOVOLU?
                        IF( NVOLTE( NTOP ) .EQ. NOVOLU ) THEN
cccC                          LA FACE NFOP DE NTOP EST ELLE DANS LEFACO?
ccc                           CALL NULEFT( NFOP,NTOP,NOTETR,MXFACO,LEFACO,
ccc     %                                  NFLE )
ccc                           IF( NFLE .LE. 0 ) THEN
cccC                          LA FACE N'EST PAS DANS LEFACO
                           LHPILE = LHPILE + 1
                           LAPILE(1,LHPILE) = NTOP
                           LAPILE(2,LHPILE) = NFOP
ccc                           ENDIF
                        ENDIF
                     ENDIF

                  ENDDO

               ENDIF

 70         ENDDO

C           RETOUR EN HAUT DE LA PILE
            GOTO 50

         ENDIF

C           SOMME DES ECHANGES DE TETRAEDRES
 80         NB2T3T = NB2T3T + NBEC2T3T
            NB3T2T = NB3T2T + NBEC3T2T

ccc            PRINT*,'tetrajpa: Tetraedre NTEV=',NTEV,
ccc     %      ' Ajout des POINTS',NBSOMM-NBPARE+1,' a',NBSOMM,
ccc     %      ' de l'' ARETE',NS1,NS2,
ccc     %      ' de taille',AREMAX,' mais SOUHAITEE',TI1,TI2,
ccc     %      ' Nb2T3T=',NBEC2T3T,' NB3T2T=',NBEC3T2T
ccc            PRINT*,'NTEOLD=',(NTEOLD(K),K=1,NBTEOLD),
ccc     %      '  =>  NTENEW=',(NTENEW(K),K=1,NBTENEW)

C           PASSAGE AU TETRAEDRE SUIVANT DE NTEV DU VOLUME NOVOLU A TRAITER
C           ---------------------------------------------------------------
 100     ENDDO

C        PASSAGE AU VOLUME SUIVANT DE LA PARTITION
C        -----------------------------------------
 200  ENDDO


C     QUALITE DES TETRAEDRES APRES AJOUT DES POINTS DES ARETES TROP LONGUES
C     ---------------------------------------------------------------------
 400  CALL QUALTETR( PTXYZD, NUDTETR+1, NOTETR,
     %               NBTETR, NUDTETR, QUAMIN, QUAMOY, VOLUMT )

      WRITE(IMPRIM,*)'tetrajpa: ITERATION', ITERAT,
     % NBSOMM-NBSOM0,' POINTS AJOUTES et TETRAEDRISES avec QUAMOY=',
     % QUAMOY,' QUAMIN=',QUAMIN,' NB2T3T=',NB2T3T,' NB3T2T=',NB3T2T,
     %' VOLUME=',VOLUMT,' VOLMOY=',VOLMOY,
     %' NBTETR=',NBTETR,' NUDTETR=',NUDTETR
      WRITE(IMPRIM,*)

      IF( NBSOMM .LT. MXSOMM .AND. NBSOM0 .LT. NBSOMM .AND.
     %    ITERAT .LE. 2 ) THEN

C        UNE ITERATION D'AJOUT DE POINTS SUR L'ARETE MAX EN PLUS
C        -------------------------------------------------------
         NUDTET0 = NUDTETR
         GOTO 10

      ENDIF

C     FIN D'AJOUT DES POINTS DES ARETES MAX TROP LONGUES
      WRITE(IMPRIM,*)'tetrajpa: FIN des ITERATIONS:', NBSOMM-NBSOM00,
     %' POINTS AJOUTES et TETRAEDRISES avec QUAMOY=',QUAMOY,
     %' QUAMIN=',QUAMIN,' NB2T3T=',NB2T3T,' NB3T2T=',NB3T2T,
     %' VOLUME=',VOLUMT,' VOLMOY=',VOLMOY,
     %' NBTETR=',NBTETR,' NUDTETR=',NUDTETR

      IERR = 0
      GOTO 9999

 9990 PRINT *,'tetrajpa: TABLEAU LAPILE SATURE  MXPILE=',MXPILE
      IERR = 1

 9999 RETURN
      END
