      SUBROUTINE TEAJBAPA( ARETGR, NOFOTI, NBVOPA,
     %                     NBSOMM, MXSOMM, PTXYZD, NPSOFR,
     %                     MXTETR, N1TEVI, NOTETR, NUDTETR, N1TETS,
     %                     NOF1VO, NVOLTE,
     %                     MXFACO, LEFACO, MXPILE, LAPILE,
     %                     IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   TETRAEDRISATION DE POINTS SUR L'ARETE MAX DES TETRAEDRES
C -----   SI SA TAILLE EST SUPERIEURE A LA TAILLE SOUHAITEE EN AU
C         MOINS UN DES SOMMETS PUIS
C         ESSAI D'AMELIORER LES TETRAEDRES AU DELA DES DECOUPES
C         PAR 2T->3T ou mT->2m-4T

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
C          1 SATURATION DES SOMMETS    TABLEAU PTXYZD
C          2 SATURATION DES TETRAEDRES TABLEAU NOTETR
C          2 SATURATION DE LA PILE     TABLEAU LAPILE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET Saint PIERRE du PERRAY          Septembre 2018
C23456...............................................................012
      PARAMETER        (MXTOLD=2048, MXTNEW=2*MXTOLD)
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      COMMON / TRTETR / STOPTE, TRACTE
      LOGICAL           STOPTE, TRACTE, TRACTE0
      CHARACTER*80      KTITRE

      DOUBLE PRECISION  PTXYZD(4,MXSOMM), VTE, VOLDES,
     %                  AREMAX, V, VOLUMT, D, D1, D2, D3, D4, TI1, TI2,
     %                  LONARE, VOLTET, DBLE, C1, C2,
     %                  ARMIN, ARMAX, SURFTR(4), ARMOY
      INTEGER           NPSOFR(MXSOMM), LEFACO(1:11,0:MXFACO),
     %                  NOTETR(8,MXTETR), N1TETS(MXSOMM), 
     %                  NOF1VO(3,NBVOPA),
     %                  LAPILE(2,MXPILE), NVOLTE(1:*),
     %                  NOTOLD(MXTOLD), NOTNEW(MXTNEW), NO4TCR(4),
     %                  NOSOTR(3),
     %                  NTPARE(-1:62), NPTARE(0:63)
ccc      INTEGER           NOSOTE(4)
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
      WRITE(IMPRIM,*)'teajbapa: MISE A JOUR de la TAILLE des ARETES par
     % AJOUT de POINTS sur l''ARETE MAX TROP LONGUE des TETRAEDRES et du
     % BARYCENTRE des TETRAEDRES TROP VOLUMINEUX'

      IVOLTE = 1
      NBT    = 0
      NB2T3T = 0
      NBMTNT = 0
      NBSOM00= NBSOMM
      QTEMIT = 0.1
      ITERAT = 0

C     BOUCLE D'ITERATIONS SUR LA RECHERCHE DES ARETES TROP LONGUES
C     ============================================================
 10   ITERAT = ITERAT + 1
      PRINT*
      PRINT*,'teajbapa: Debut ITERATION',ITERAT

C     BOUCLE SUR LES NBVOPA VOLUMES DE LA PARTITION
C     =============================================
      NUDTET0= NUDTETR
      QTEMIT = QTEMIT / ITERAT
      NBSOM0 = NBSOMM
      VOLUMT = 0D0

      DO 200 NV = 1, NBVOPA

C        LE NO DE VOLUME DE 1 A NBVOPA DU TETRAEDRE NTE ET SUIVANTS
         NOVOLU = NOF1VO( 3, NV )

C        BOUCLE SUR LES TETRAEDRES DU VOLUME NOVOLU
C        ------------------------------------------
         DO 100 NTEV = 1, NUDTET0

            IF( NOTETR(1,NTEV) .EQ. 0 .OR. NVOLTE(NTEV) .NE. NOVOLU )
     %         GOTO 100

C           QUALITE ET VOLUME DU TETRAEDRE NTEV
            CALL QUATETD( PTXYZD( 1, NOTETR(1,NTEV) ),
     %                    PTXYZD( 1, NOTETR(2,NTEV) ),
     %                    PTXYZD( 1, NOTETR(3,NTEV) ),
     %                    PTXYZD( 1, NOTETR(4,NTEV) ),
     %                    ARMIN, ARMAX, SURFTR, VTE, QTE )

C           PAS D'AJOUT SUR LES TETRAEDRES DE MAUVAISE QUALITE
            IF( QTE .LE. QTEMIT ) GOTO 100

            NS1AMX = 0
            NS2AMX = 0
            NBARMX = 0

C           MISE A JOUR DE LA TAILLE SOUHAITEE DES ARETES
C           AUX 4 SOMMETS DU TETRAEDRE NTEV
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

C           CALCUL DU VOLUME DESIRE DE CE TETRAEDRE NTEV
            D1 = PTXYZD( 4, NOTETR(1,NTEV) )
            D2 = PTXYZD( 4, NOTETR(2,NTEV) )
            D3 = PTXYZD( 4, NOTETR(3,NTEV) )
            D4 = PTXYZD( 4, NOTETR(4,NTEV) )

ccc         D      = MIN( D1 , D2 , D3 , D4 )  13/7/2010
ccc         VOLDES = ( D1 * D2 * D3 * D4 ) ** 0.75

ccc         D      = MAX( D1 , D2 , D3 , D4 ) 19/12/2017
ccc         VOLDES = D1 * D2 * D3 / D * D4

ccc         D      = ( D1 + D2 + D3 + D4 ) / 4
ccc         VOLDES = D * D * D

C           VOLUME DU TETRAEDRE AVEC L'ARETE MOYENNE DESIREE
            ARMOY = ( D1 + D2 + D3 + D4 ) / 4
            VOLDES = ARMOY**3 / 6

            IF( VTE .LE. VOLDES .AND. ARMAX .LE. MAX(D1,D2,D3,D4) ) THEN
C              VOLUME SUFFISAMMENT FAIBLE et ARETE MAX SUFFISAMMENT COURTE
C              => PAS D'AJOUT DU BARYCENTRE NI DE RECHERCHE D'ARETE MAX
               GOTO 100
            ENDIF

C           ------------------------------------------------------------
C           AJOUT DU BARYCENTRE NBSOMM DU TETRAEDRE NTEV TROP VOLUMINEUX
C           ------------------------------------------------------------
            IF( NBSOMM .GE. MXSOMM ) THEN
               PRINT *,'teajbapa: TABLEAU PTXYZD SATURE'
               IERR = 1
               GOTO 400
            ENDIF

C           LES 3 COORDONNEES ET DISTANCE SOUHAITEE DU BARYCENTRE
            NBSOMM = NBSOMM + 1

C           LE BARYCENTRE EST AVEC LES POIDS 1/TailleIdeale
            D1 = 1D0 / D1
            D2 = 1D0 / D2
            D3 = 1D0 / D3
            D4 = 1D0 / D4
            D  = D1 + D2 + D3 + D4
            DO I=1,4
               PTXYZD(I,NBSOMM) = ( PTXYZD(I,NOTETR(1,NTEV)) * D1
     %                            + PTXYZD(I,NOTETR(2,NTEV)) * D2
     %                            + PTXYZD(I,NOTETR(3,NTEV)) * D3
     %                            + PTXYZD(I,NOTETR(4,NTEV)) * D4 ) / D
            ENDDO

C           CALCUL DE TAILLE_IDEALE( NBSOMM )
            IF( NOFOTI .GT. 0 ) THEN
C              LA FONCTION TAILLE_IDEALE(X,Y,Z) EXISTE
               CALL FONVAL( NOFOTI ,3, PTXYZD(1,NBSOMM), NCODEV,
     %                      PTXYZD(4,NBSOMM) )
               IF( NCODEV .LE. 0 ) THEN
                  PTXYZD(4,NBSOMM) = ARETGR
               ENDIF
            ENDIF
           IF(PTXYZD(4,NBSOMM).LE.0D0.OR.PTXYZD(4,NBSOMM).GT.ARETGR)THEN
              PTXYZD(4,NBSOMM) = ARETGR
            ENDIF

cccC           SIMULATION DES 4 TETRAEDRES DE SOMMET LE BARYCENTRE NBSOMM
cccC           POUR EVITER DE LES CREER SI L'UN D'ENTRE EUX EST
cccC           DE TROP PETIT VOLUME PAR RAPPORT AU VOLUME DE NTEV
ccc            DO NF = 1, 4
cccC              3-ERS SOMMETS DE LA FACE DANS LE SENS DIRECT DU TETRAEDRE NTEV
ccc               NOSOTE( 1 ) = NOTETR( NOSOFATE(1,NF), NTEV )
ccc               NOSOTE( 2 ) = NOTETR( NOSOFATE(2,NF), NTEV )
ccc               NOSOTE( 3 ) = NOTETR( NOSOFATE(3,NF), NTEV )
ccc               NOSOTE( 4 ) = NBSOMM
cccC              LES COORDONNEES DU CENTRE DE LA BOULE CIRCONSCRITE
cccC              AU TETRAEDRE NOSOTE ET CARRE DE SON RAYON
ccc               CALL CEBOQU( 0.0, NOSOTE, PTXYZD, CENTRE,  V, Q )
ccc               IF( V .LT. VOLDES/64 ) THEN
cccC                 TETRAEDRE DE TROP PETIT VOLUME
cccC                 => PAS D'AJOUT DU BARYCENTRE NBSOMM
ccc                  NBSOMM = NBSOMM - 1
ccc                  GOTO 100
ccc               ENDIF
ccc            ENDDO

            IF( MOD(NBSOMM,10000) .EQ. 0 ) THEN
               IF( LANGAG .EQ. 0 ) THEN
               WRITE(IMPRIM,*)'teajbapa: TETRAEDRISATION du BARYCENTRE',
     %                        NBSOMM,' du TETRAEDRE',NTEV,
     %                        ' Iteration',ITERAT
               ELSE
               WRITE(IMPRIM,*)'teajbapa: TETRAHEDRIZATION of BARYCENTRE'
     %                     ,NBSOMM,' of TETRAHEDRON',NTEV,
     %                     ' Iteration',ITERAT
               ENDIF
            ENDIF

C           BARYCENTRE INTERIEUR DE 4 SOMMETS QUELCONQUES NON SUR UNE SURFACE
            NPSOFR( NBSOMM ) = 0

C           GENERATION DES 4 TETRAEDRES DE SOMMET LE BARYCENTRE NBSOMM
C           ----------------------------------------------------------
            NBTOLD = 1
            NOTOLD( 1 ) = NTEV
            NB4TCR = 0
            DO NF = 1, 4

C              LE TETRAEDRE A CREER NT EST LE PREMIER TETRAEDRE VIDE
               NB4TCR = NB4TCR + 1
               NT     = N1TEVI
               NO4TCR( NB4TCR ) = N1TEVI

C              MISE A JOUR DU 1-ER TETRAEDRE VIDE
               N1TEVI = NOTETR( 5, N1TEVI )

C              NUMERO MAX DES TETRAEDRES UTILISES DANS NOTETR
               NUDTETR = MAX( NUDTETR, NT )

C              3-ERS SOMMETS DE LA FACE DANS LE SENS DIRECT DU TETRAEDRE NTEV
               NOTETR( 1, NT ) = NOTETR( NOSOFATE(1,NF), NTEV )
               NOTETR( 2, NT ) = NOTETR( NOSOFATE(2,NF), NTEV )
               NOTETR( 3, NT ) = NOTETR( NOSOFATE(3,NF), NTEV )
               NOTETR( 4, NT ) = NBSOMM

C              LE CHAINAGE DES TETRAEDRES ADJACENTS PAR LES FACES
C              LE TETRAEDRE OPPOSE A LA FACE AU DELA DE LA FACE
C              ATTENTION: IL PEUT NE PAS EXISTER DE TEL TETRAEDRE
C              NTOP=0 POUR UNE FACE FRONTIERE SANS TETRAEDRE OPPOSE
C                  -1 POUR UN TETRAEDRE OPPOSE INCONNU
               NTOP = NOTETR( 4+NF, NTEV )

               NOTETR( 5, NT ) = NTOP
               NOTETR( 6, NT ) = -1
               NOTETR( 7, NT ) = -1
               NOTETR( 8, NT ) = -1

C              MISE A JOUR D'UN TETRAEDRE CONTENANT CHAQUE SOMMET
C              TOUS LES TETRAEDRES ONT PU DISPARAITRE
C              CAS DE L'ETOILE REDUITE AUX 2 TETRAEDRES INITIAUX
               DO N=1,4
                  N1TETS( NOTETR( N, NT ) ) = NT
               ENDDO

C              AJOUT DE NTE AU VOLUME NOVOLU
               NVOLTE( NT ) = NOVOLU

            ENDDO

C           COMPLETION DES 4 TETRAEDRES OPPOSES CREES DANS L'ETOILE
C           -------------------------------------------------------
            CALL MJOPTE( NB4TCR, NO4TCR, N1TETS, NOTETR, NUDTETR,
     %                   N1TEVI, PTXYZD, NBFANR )
            IF( NBFANR .GT. 0 ) THEN
            PRINT*,'teajbapa: Probleme sur les TETRAEDRES OPPOSES aux TE
     %TRAEDRES NO4TCR'
            ENDIF
  
C           SUPPRESSION DU TETRAEDRE NTEV DU TABLEAU NOTETR
C           AJOUT DU TETRAEDRE NTEV A LA CHAINE DES TETRAEDRES VIDES
            NOTETR( 1, NTEV ) = 0
            NOTETR( 5, NTEV ) = N1TEVI
            N1TEVI = NTEV
C           SUPPRESSION DE NTEV DU VOLUME NOVOLU
            NVOLTE( NTEV ) = -1


C           INITIALISATION DE LA PILE DES FACES EXTERNES DES 4 TETRAEDRES NO4TCR
C           --------------------------------------------------------------------
            LHPILE = 0
            DO N = 1, NB4TCR
C              RECHERCHE DU TETRAEDRE OPPOSE NTOP ET DE SA FACE NFOP
C              A LA FACE 1 DU TETRAEDRE NO4TCR(N)
               CALL NOFAOP( 1, NO4TCR(N), NOTETR,  NFOP, NTOP )
               IF( NFOP .GT. 0 ) THEN

C                 LE TETRAEDRE NTOP EST IL DANS LE VOLUME NOVOLU?
                  IF( NVOLTE( NTOP ) .EQ. NOVOLU ) THEN
cccC                    LA FACE NFOP DE NTOP EST ELLE DANS LEFACO?
ccc                     CALL NULEFT( NFOP, NTOP, NOTETR, MXFACO, LEFACO, NFLE)
ccc                     IF( NFLE .LE. 0 ) THEN
cccC                    LA FACE N'EST PAS DANS LEFACO. ELLE EST EMPILEE
                     LHPILE = LHPILE + 1
                     LAPILE(1,LHPILE) = NTOP
                     LAPILE(2,LHPILE) = NFOP
ccc                     ENDIF
                  ENDIF

               ENDIF
            ENDDO

C           ESSAI D'AMELIORER LA QUALITE DES 4 TETRAEDRES NO4TCR
            GOTO 49

C           =====================================================
C           RECHERCHE DE L'ARETE MAX DES NBTNEW TETRAEDRES NOTNEW
C           =====================================================
 19         NBARMX = NBARMX + 1
            NAREMX = 0
            AREMAX = 0D0
            NTEMAX = 0
            NBPARE = 0
            DO 27 NT = 1, NB4TCR

               NTEN = NO4TCR( NT )
               IF( NTEN .LE. 0 .OR. NOTETR(1,NTEN) .EQ. 0 ) GOTO 27

               DO 25 NA = 1, 6

C                 NUMERO DES 2 SOMMETS DE L'ARETE NA
                  NS1 = NOTETR( NOSOARTE(1,NA), NTEN )
                  NS2 = NOTETR( NOSOARTE(2,NA), NTEN )

C                 L'ARETE NS1-NS2 EST ELLE SUR LA FRONTIERE?
                  IF( NPSOFR(NS1) .GT. 0 .AND. NPSOFR(NS2) .GT. 0 ) THEN
C                    ARETE MAX FRONTALIERE
                     GOTO 25
                  ENDIF

                  LONARE = ( PTXYZD(1,NS2) - PTXYZD(1,NS1) ) ** 2
     +                   + ( PTXYZD(2,NS2) - PTXYZD(2,NS1) ) ** 2
     +                   + ( PTXYZD(3,NS2) - PTXYZD(3,NS1) ) ** 2
                  IF( LONARE .GT. AREMAX ) THEN
                     NAREMX = NA
                     AREMAX = LONARE
                     NTEMAX = NTEN
                  ENDIF

 25            ENDDO

 27         ENDDO

            IF( NAREMX .EQ. 0 ) THEN
C              LES 4 SOMMETS DU TETRAEDRE SONT SUR LA FRONTIERE
C              => PAS D'AJOUT DE POINTS SUR L'ARETE
               GOTO 100
            ENDIF
            AREMAX = SQRT( AREMAX )

C           TAILLE SOUHAITEE AUX 2 SOMMETS DE L'ARETE MAX
            NS1AMX = NOTETR( NOSOARTE(1,NAREMX), NTEMAX )
            NS2AMX = NOTETR( NOSOARTE(2,NAREMX), NTEMAX )
            TI1 = PTXYZD( 4, NS1AMX )
            TI2 = PTXYZD( 4, NS2AMX )

C           NBPARE NOMBRE DE POINTS A AJOUTER ENTRE NS1AMX et NS2AMX
            IF( AREMAX .LE. (TI1+TI2) * 0.65D0 ) THEN
               NBPARE = 0
C              => PAS D'AJOUT DE POINT SUR L'ARETE NS1AMX-NS2AMX
               GOTO 100
            ENDIF
            IF( AREMAX .LE. (TI1+TI2) * 1.3D0 ) THEN
               NBPARE = 1
            ELSE
               NBPARE = NINT( AREMAX / ( TI1 + TI2 ) )
            ENDIF

C           RECHERCHE DES NBTOLD TETRAEDRES D'ARETE NS1AMX NS2AMX
            CALL TETR1A( NS1AMX, NS2AMX, N1TETS, NOTETR,
     %                   NBTOLD, MXTOLD, NOTOLD, IERR )

            IF( IERR .NE. 0 .OR. NBTOLD .GT. 28 ) THEN
C              TRACE DES TETRAEDRES NOTOLD
               TRACTE0 = TRACTE
               TRACTE  = .TRUE.
               KTITRE='teajbapa:       TROP de TETRAEDRES AUTOUR DE L''A
     %RETE              '
               WRITE( KTITRE(11:15),'(I5)') NBTOLD
               WRITE( KTITRE(55:61),'(I7)') NS1AMX
               WRITE( KTITRE(63:69),'(I7)') NS2AMX
               CALL TRFETO15( KTITRE, PTXYZD,
     %                        NBTOLD, NOTOLD, NOTETR, NS1AMX, NS2AMX )
               TRACTE = TRACTE0
C               => PAS D'AJOUT D'UN POINT SUR L'ARETE COMMUNE CAR
C                  TROP DE TETRAEDRES
               GOTO 100
            ENDIF

cccC           TRACE DES TETRAEDRES NOTOLD
ccc            TRACTE0 = TRACTE
ccc            TRACTE  = .TRUE.
ccc            KTITRE='teajbapa:       OLD TETRAEDRES AUTOUR DE L''ARETE MA
ccc     %X                 NB Pt/Arete=    '
ccc            WRITE( KTITRE(11:15),'(I5)') NBTOLD
ccc            WRITE( KTITRE(55:61),'(I7)') NS1AMX
ccc            WRITE( KTITRE(63:69),'(I7)') NS2AMX
ccc            WRITE( KTITRE(84:86),'(I3)') NBPARE
ccc            CALL TRFETO15( KTITRE, PTXYZD, NBTOLD, NOTOLD, NOTETR,
ccc     %                     NS1AMX, NS2AMX )
ccc            TRACTE = TRACTE0

C           AJOUT DE NBPARE POINTS SUR L'ARETE MAXIMALE DU TETRAEDRE NTEMAX
C           ---------------------------------------------------------------
C           CALCUL DES XYZ DES NBPARE POINTS SUR L'ARETE MAX DE NTEMAX
            IF( NBSOMM+NBPARE .GT. MXSOMM ) THEN
               PRINT *,'teajbapa: TABLEAU PTXYZD SATURE MXSOMM=',MXSOMM
               IERR = 1
               GOTO 400
            ENDIF

            DO NP = 1, NBPARE

               NBSOMM = NBSOMM + 1
               NPTARE( NP ) = NBSOMM

               IF( MOD(NBSOMM,1000) .EQ. 0 ) THEN
                  IF( LANGAG .EQ. 0 ) THEN
                    WRITE(IMPRIM,*)'teajbapa: TETRAEDRISATION du POINT',
     %              NBSOMM,' de l''ARETE',NS1AMX,NS2AMX,
     %              ' du TETRAEDRE',NTEMAX,' NB Pt/Arete=',NBPARE,
     %              ' Iteration',ITERAT
                  ELSE
                   WRITE(IMPRIM,*)'teajbapa: TETRAHEDRIZATION of POINT',
     %             NBSOMM,' of EDGE',NS1AMX,NS2AMX,
     %                   ' of TETRAHEDRON',NTEMAX,' Pt/Edge Nb=',NBPARE,
     %                   ' Iteration',ITERAT

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
                  PTXYZD(K,NBSOMM) = C1 * PTXYZD(K,NS1AMX)
     %                             + C2 * PTXYZD(K,NS2AMX)
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
               IF( PTXYZD(4,NBSOMM) .LE. 0D0 .OR.
     %             PTXYZD(4,NBSOMM) .GT. ARETGR ) THEN
                   PTXYZD(4,NBSOMM) = ARETGR
               ENDIF

            ENDDO

C           LE NO DE POINT DES 2 EXTREMITES DE L'ARETE
            NPTARE( 0 ) = NS1AMX
            NPTARE( NBPARE+1 ) = NS2AMX

C           DECOUPAGE EN NBPARE+1 SOUS-TETRAEDRES DE CHACUN DES NBTOLD
C           TETRAEDRES PAR AJOUT DES NBPARE POINTS SUR L'ARETE NS1AMX NS2AMX
            IF( NBTOLD*(NBPARE+1) .GT. MXTNEW ) THEN
               PRINT*,'teajbapa: TABLEAU NOTNEW SATURE MXTNEW=',MXTNEW
               NBSOMM = NBSOMM - NBPARE
               GOTO 100
            ENDIF

C           INITIALISATION DE LA PILE AVEC LA FACE EXTERNE DES TETRAEDRES
C           -------------------------------------------------------------
            NS4    = 0
            NBTNEW = 0
            LHPILE = 0
            DO N = 1, NBTOLD

               NTE = NOTOLD( N )

C              QUELS SONT LES NUMEROS DES 2 AUTRES SOMMETS DE NTE?
               DO M = 1, 4
                  NS3 = NOTETR( M, NTE )
                  IF( NS3 .NE. NS1AMX .AND. NS3 .NE. NS2AMX ) GOTO 33
               ENDDO

 33            DO MM = M+1, 4
                  NS4 = NOTETR( MM, NTE )
                  IF( NS4 .NE. NS1AMX .AND. NS4 .NE. NS2AMX .AND.
     %                NS4 .NE. NS3 )  GOTO 34
               ENDDO

C              VOLUME DU TETRAEDRE NTE
 34            V = VOLTET( PTXYZD( 1, NS1AMX ), PTXYZD( 1, NS2AMX ),
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
                     PRINT *,'teajbapa: TABLEAU NOTETR SATURE'
                     IERR = 2
                     GOTO 400
                  ENDIF
C                 LE TETRAEDRE A CREER EST LE PREMIER TETRAEDRE VIDE
                  NBTNEW = NBTNEW + 1
                  NT1    = N1TEVI
                  NOTNEW( NBTNEW ) = NT1
C                 MISE A JOUR DU 1-ER TETRAEDRE VIDE
                  N1TEVI = NOTETR( 5, N1TEVI )
C                 NUMERO MAX DES TETRAEDRES UTILISES DANS NOTETR
                  NUDTETR = MAX( NUDTETR, NT1 )
                  NTPARE( NP ) = NT1
C                 NUMERO DE VOLUME DU TETRAEDRE
                  NVOLTE( NT1 ) = NOVOLU
               ENDDO

C              LE NUMERO LOCAL NF DE LA FACE EXTERIEURE DE SOMMETS
C              NS4-NS3-NS1MAX PARMI LES 4 SOMMETS DU TETRAEDRE NTE
               NOSOTR(1) = NS3
               NOSOTR(2) = NS4
               NOSOTR(3) = NS1AMX
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

C              LE NUMERO LOCAL NF DE LA FACE EXTERIEURE DE SOMMETS
C              NS2AMX-NS4-NS3 PARMI LES 4 SOMMETS DU TETRAEDRE NTE
               NOSOTR(1) = NS2AMX
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

C              LA FACE NS1AMX NS2AMX NS3 DE NTE EST ELLE SUR LA FRONTIERE?
               NOSOTR(1) = NS1AMX
               NOSOTR(2) = NS2AMX
               NOSOTR(3) = NS3
               CALL TRI3NO( NOSOTR, NOSOTR )
               CALL NO1F1T( NOSOTR, NOTETR(1,NTE), NF )
C              NO NOTETR DU TETRAEDRE OPPOSE A LA FACE NF DE NTE
               NTOP123 = NOTETR( 4+NF, NTE )

C              LA FACE NS1AMX NS2AMX NS4 DE NTE EST ELLE SUR LA FRONTIERE?
               NOSOTR(1) = NS1AMX
               NOSOTR(2) = NS2AMX
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
C                    FACE INCONNUE
                     NOTETR( 8, NT1 ) = -1
                  ENDIF

C                 MISE A JOUR D'UN TETRAEDRE CONTENANT CHAQUE SOMMET
                  DO K=1,4
                     N1TETS( NOTETR(K,NT1) ) = NT1
                  ENDDO

C                 AJOUT DE NTE AU VOLUME NOVOLU
                  NVOLTE( NT1 ) = NOVOLU

               ENDDO

            ENDDO

C           COMPLETION DES TETRAEDRES OPPOSES DES TETRAEDRES CREES NOTNEW
C           -------------------------------------------------------------
            CALL MJOPTE( NBTNEW, NOTNEW, N1TETS, NOTETR, NUDTETR,
     %                   N1TEVI, PTXYZD, NBFANR )
            IF( NBFANR .GT. 0 ) THEN
               PRINT*,'teajbapa: Probleme sur les TETRAEDRES OPPOSES aux
     % TETRAEDRES NOTNEW'
            ENDIF

C           SUPPRESSION DES TETRAEDRES NOTOLD DU TABLEAU NOTETR ET NVOLTE
C           -------------------------------------------------------------
            DO N = 1, NBTOLD
               NTE = NOTOLD( N )
C              AJOUT DU TETRAEDRE NTE A LA CHAINE DES TETRAEDRES VIDES
               NOTETR( 1, NTE ) = 0
               NOTETR( 5, NTE ) = N1TEVI
               N1TEVI = NTE
C              SUPPRESSION DE NTE DU VOLUME NOVOLU
               NVOLTE( NTE ) = -1
            ENDDO


C           ===========================================================
C           ESSAI D'AMELIORER LES TETRAEDRES OPPOSES AUX FACES EMPILEES
C           PAR 2T->3T ou mT->2m-4T et EMPILAGE DES FACES CREEES
C           ===========================================================
 49         NBEC2T3T = 0
            NBECMTNT = 0

C           TRAITEMENT DU HAUT DE LA PILE DES FACES OPPOSEES
C           ------------------------------------------------
 50         IF( LHPILE .GT. 0 ) THEN

C              RECUPERATION DU TETRAEDRE NTE ET DE SA FACE NFE A TRAITER
               NTE = LAPILE( 1, LHPILE )
               NFE = LAPILE( 2, LHPILE )
               LHPILE = LHPILE - 1

C              LE TETRAEDRE NTE EXISTE-T-IL ENCORE?
               IF( NTE .LE. 0 .OR. NFE .LE. 0 ) GOTO 50
               IF( NOTETR(1,NTE) .LE. 0 ) GOTO 50

C              LE TETRAEDRE NTE EST IL DANS LE VOLUME NOVOLU?
               IF( NVOLTE( NTE ) .NE. NOVOLU ) GOTO 50

cccC           LA FACE NFE DE NTE EST ELLE DANS LEFACO?
ccc            CALL NULEFT( NFE, NTE, NOTETR, MXFACO, LEFACO, NFLE )
ccc            IF( NFLE .GT. 0 ) THEN
cccC              OUI: LA FACE EST DANS LEFACO
cccC              CONSERVATION DE CETTE FACE LEFACO EN NE PASSANT PAS AU DELA
ccc               GOTO 50
ccc            ENDIF


C              TENTATIVE DE CHANGER 2T->3T SI MINIMUM QUALITE EST MAXIMISE
C              -----------------------------------------------------------
C              LE NUMERO NFOP DE LA FACE COMMUNE DANS LE TETRAEDRE OPPOSE NTOP
               CALL NOFAOP( NFE, NTE, NOTETR, NFOP, NTOP )
               IF( NFOP .LE. 0 ) GOTO 50

C              LE TETRAEDRE NTOP EST IL DANS LE VOLUME NOVOLU?
               IF( NVOLTE( NTOP ) .NE. NOVOLU ) GOTO 50

               CALL CH2T3T( 1,      MXFACO, LEFACO, 1, IVOLTE, NVOLTE,
     %                      NTE,    NFE,    NTOP,
     %                      PTXYZD, N1TETS, NOTETR, N1TEVI, NUDTETR,
     %                      NOTNEW, IER )

               IF( IER .EQ. 0 ) THEN

C                 REUSSITE: CHANGEMENT 2T->3T EFFECTUE
                  NBEC2T3T = NBEC2T3T + 1
                  NBTOLD = 2
                  NOTOLD( 1 ) = NTE
                  NOTOLD( 2 ) = NTOP
                  NVOLTE( NTE  ) = -1
                  NVOLTE( NTOP ) = -1
                  NBTNEW = 3

C                 NOTNEW=NO NOTETR DES 3 TETRAEDRES QUI REMPLACENT NTE et NTOP
C                 LES 2 FACES EXTERIEURES 1 et 4 DE CES 3 TETRAEDRES SONT EMPILEES
C                 CF ch2t3t.f
                  IF( LHPILE+6 .GT. MXPILE ) GOTO 9990
                  DO K = 1, NBTNEW

C                    LE NOUVEAU TETRAEDRE K
                     NT = NOTNEW( K )
                     NVOLTE( NT ) = NOVOLU

C                    RECHERCHE DU TETRAEDRE OPPOSE NTOP ET DE SA FACE NFOP
C                    A LA FACE 1 DU TETRAEDRE NT
                     CALL NOFAOP( 1, NT, NOTETR,  NFOP, NTOP )
                     IF( NFOP .GT. 0 ) THEN
C                       LE TETRAEDRE NTOP EST IL DANS LE VOLUME NOVOLU?
                        IF( NVOLTE( NTOP ) .EQ. NOVOLU ) THEN
cccC                          LA FACE NFOP DE NTOP EST ELLE DANS LEFACO?
ccc                           CALL NULEFT( NFOP,NTOP,NOTETR,MXFACO,LEFACO,
ccc     %                                  NFLE)
ccc                           IF( NFLE .LE. 0 ) THEN
C                             LA FACE N'EST PAS DANS LEFACO
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
C                             LA FACE N'EST PAS DANS LEFACO
                           LHPILE = LHPILE + 1
                           LAPILE(1,LHPILE) = NTOP
                           LAPILE(2,LHPILE) = NFOP
ccc                           ENDIF
                        ENDIF
                     ENDIF

                  ENDDO

                  GOTO 50

               ENDIF


C              TENTATIVE DE CHANGER mT->2m-4T SI MIN QUALITE EST MAXIMISE
C              ----------------------------------------------------------
C              PARCOURS DES 3 ARETES DE LA FACE NFE DU TETRAEDRE NTE
               DO 70 N = 1, 3

C                 NO DE L'ARETE N DE LA FACE NFE
                  NAR = NOARFATE( N, NFE )

C                 NUMERO PTXYZD DES 2 SOMMETS DE L'ARETE NAR
                  NST1 = NOTETR( NOSOARTE( 1, NAR ), NTE )
                  NST2 = NOTETR( NOSOARTE( 2, NAR ), NTE )

C                 L'ARETE NST1-NST2 PERFORE LA FACE COMMUNE AUX 2 TETRAEDRES?
C                 2 TETRAEDRES => 3 TETRAEDRES D'ARETE COMMUNE NST1-NST2
                  CALL TETR1A( NST1,   NST2,   N1TETS, NOTETR,
     %                         NBTOLD, MXTOLD, NOTOLD, IERR )
                  IF( NBTOLD .GT. 8 ) GOTO 70

                  DO M=1,NBTOLD
                     IF( NVOLTE( NOTOLD(M) ) .NE. NOVOLU ) GOTO 70
                  ENDDO

                  CALL CHMTNT( 1,      MXFACO, LEFACO, IVOLTE,NVOLTE,
     %                         PTXYZD, NTE,    NAR, 
     %                         N1TETS, NOTETR, N1TEVI, NUDTETR,
     %                         MXTOLD, NBTOLD, NOTOLD,
     %                         MXTNEW, NBTNEW, NOTNEW, IER )

                  IF( IER .EQ. 0 ) THEN

C                    REUSSITE: CHANGEMENT mT->2m-4T EFFECTUE
                     NBECMTNT = NBECMTNT + 1

ccc                      PRINT*,'teajbapa: mT->2m-4T AMELIORE l''ETOILE'

C                    NOTNEW: NUMERO NOTETR DES 2m-4 TETRAEDRES QUI REMPLACENT NOTOLD(1:NBTOLD)
C                    LES FACES EXTERIEURES DE CES 2m-4 TETRAEDRES SONT EMPILEES
C                    CF chmtnt.f

                     IF( LHPILE+NBTNEW .GT. MXPILE ) GOTO 9990

                     DO K = 1, NBTNEW

C                       LE NOUVEAU TETRAEDRE K
                        NT = NOTNEW( K )
                        NVOLTE( NT ) = NOVOLU

C                       RECHERCHE DU TETRAEDRE OPPOSE NTOP ET DE SA FACE NFOP
C                       A LA FACE 2 DU TETRAEDRE NT
                        CALL NOFAOP( 2, NT, NOTETR,  NFOP, NTOP )
                        IF( NFOP .GT. 0 ) THEN
C                          LE TETRAEDRE NTOP EST IL DANS LE VOLUME NOVOLU?
                           IF( NVOLTE( NTOP ) .EQ. NOVOLU ) THEN
cccC                             LA FACE NFOP DE NTOP EST ELLE DANS LEFACO?
ccc                              CALL NULEFT( NFOP,NTOP,NOTETR,MXFACO,LEFACO,
ccc     %                                     NFLE )
ccc                              IF( NFLE .LE. 0 ) THEN
cccC                                LA FACE N'EST PAS DANS LEFACO
                              LHPILE = LHPILE + 1
                              LAPILE(1,LHPILE) = NTOP
                              LAPILE(2,LHPILE) = NFOP
ccc                              ENDIF
                           ENDIF
                        ENDIF

cccC                       RECHERCHE DU TETRAEDRE OPPOSE NTOP ET DE SA FACE NFOP
cccC                       A LA FACE 3 DU TETRAEDRE NT ( INTERNE AUX 2m-4T )
ccc                        CALL NOFAOP( 3, NT, NOTETR,  NFOP, NTOP )
ccc                        IF( NFOP .GT. 0 ) THEN
cccC                          LE TETRAEDRE NTOP EST IL DANS LE VOLUME NOVOLU?
ccc                           IF( NVOLTE( NTOP ) .EQ. NOVOLU ) THEN
ccccccC                             LA FACE NFOP DE NTOP EST ELLE DANS LEFACO?
cccccc                              CALL NULEFT( NFOP,NTOP,NOTETR,MXFACO,LEFACO,
cccccc     %                                     NFLE)
cccccc                              IF( NFLE .LE. 0 ) THEN
ccccccC                                LA FACE N'EST PAS DANS LEFACO
ccc                              LHPILE = LHPILE + 1
ccc                              LAPILE(1,LHPILE) = NTOP
ccc                              LAPILE(2,LHPILE) = NFOP
cccccc                              ENDIF
ccc                           ENDIF
ccc                        ENDIF

C                       RECHERCHE DU TETRAEDRE OPPOSE NTOP ET DE SA FACE NFOP
C                       A LA FACE 4 DU TETRAEDRE NT
                        CALL NOFAOP( 4, NT, NOTETR,  NFOP, NTOP )
                        IF( NFOP .GT. 0 ) THEN
C                          LE TETRAEDRE NTOP EST IL DANS LE VOLUME NOVOLU?
                           IF( NVOLTE( NTOP ) .EQ. NOVOLU ) THEN
cccC                             LA FACE NFOP DE NTOP EST ELLE DANS LEFACO?
ccc                              CALL NULEFT( NFOP,NTOP,NOTETR,MXFACO,LEFACO,
ccc     %                                     NFLE )
ccc                              IF( NFLE .LE. 0 ) THEN
cccC                             LA FACE N'EST PAS DANS LEFACO
                              LHPILE = LHPILE + 1
                              LAPILE(1,LHPILE) = NTOP
                              LAPILE(2,LHPILE) = NFOP
ccc                              ENDIF
                           ENDIF
                        ENDIF

                     ENDDO

                  ENDIF

 70            ENDDO

C              RETOUR EN HAUT DE LA PILE
               GOTO 50

            ENDIF

C           SOMME DES ECHANGES DE TETRAEDRES
            NB2T3T = NB2T3T + NBEC2T3T
            NBMTNT = NBMTNT + NBECMTNT

            IF( NBARMX .LT. 1 ) THEN
C              ESSAI D'AJOUT DE POINTS SUR L'ARETE MAX DES NBTNEW TETRAEDRES
               GOTO 19
            ENDIF

ccc            PRINT*,'teajbapa: Tetraedre NTEV=',NTEV,
ccc     %      ' Ajout des POINTS',NBSOMM-NBPARE+1,' a',NBSOMM,
ccc     %      ' de l'' ARETE',NS1AMX,NS2AMX,
ccc     %      ' de taille',AREMAX,' mais SOUHAITEE',TI1,TI2,
ccc     %      ' Nb2T3T=',NBEC2T3T,' NBMTNT=',NBECMTNT

cccC           POUR VOIR L'ENVIRONNEMENT DES SOMMETS NS1AMX NS2AMX
cccC           RECHERCHE DES TETRAEDRES DE SOMMETS NS1AMX ou NS2AMX
ccc            CALL TETR1S( NS1AMX, N1TETS, NOTETR,
ccc     %                   NBTNEW, MXTNEW, NOTNEW, IERR )
ccc            CALL TETR1S( NS2AMX, N1TETS, NOTETR,
ccc     %                   NBTNS,  MXTNEW-NBTNEW, NOTNEW(NBTNEW+1), IERR )
ccc            NBTNEW = NBTNEW + NBTNS
cccC  ccc      UNICITE DES NO DE TETRAEDRES DU TABLEAU NOTNEW
ccc            CALL UNITABL( NOTNEW, NBTNEW )


cccC           TRACE DES TETRAEDRES DE SOMMETS NS1AMX ou NS2AMX
ccc            TRACTE0 = TRACTE
ccc            TRACTE  = .TRUE.
ccc            KTITRE='teajbapa:       NEW TETRAEDRES AUTOUR DES SOMMETS
ccc     %              '
ccc            WRITE( KTITRE(11:15),'(I5)') NBTNEW
ccc            WRITE( KTITRE(51:57),'(I7)') NS1AMX
ccc            WRITE( KTITRE(59:65),'(I7)') NS2AMX
ccc            CALL TRFETO15( KTITRE, PTXYZD, NBTNEW, NOTNEW, NOTETR,
ccc     %                     NS1AMX, NS2AMX )
ccc            TRACTE = TRACTE0

C           PASSAGE AU TETRAEDRE SUIVANT DE NTEV DU VOLUME NOVOLU A TRAITER
C           ---------------------------------------------------------------
 100     ENDDO

C        PASSAGE AU VOLUME SUIVANT DE LA PARTITION
C        -----------------------------------------
 200  ENDDO


C     SUPPRESSION DES TETRAEDRES SE TROUVANT DANS AUCUN VOLUME
C     --------------------------------------------------------
      NBELIM = 0
      DO 300 NTEV = 1, NUDTETR
         IF( NOTETR(1,NTEV) .NE. 0 .AND. NVOLTE(NTEV) .EQ. -1 ) THEN
            NBELIM = NBELIM + 1
            NOTETR(1,NTEV) = 0
            NOTETR(5,NTEV) = N1TEVI
            N1TEVI         = NTEV
            NVOLTE( NTEV ) = -1
         ENDIF
 300  ENDDO
      PRINT*,'teajbapa: ITERATION',ITERAT,
     %' NB de TETRAEDRES dans AUCUN VOLUME=',NBELIM


C     QUALITE DES TETRAEDRES APRES AJOUT DES POINTS DES ARETES TROP LONGUES
C     ---------------------------------------------------------------------
 400  CALL QUALTETR( PTXYZD, NUDTETR+1, NOTETR,
     %               NBTETR, NUDTETR, QUAMIN, QUAMOY, VOLUMT )

      WRITE(IMPRIM,*)'teajbapa: ITERATION', ITERAT,
     % NBSOMM-NBSOM0,' POINTS AJOUTES et TETRAEDRISES avec QUAMOY=',
     % QUAMOY,' QUAMIN=',QUAMIN,' NB2T3T=',NB2T3T,' NBMTNT=',NBMTNT,
     %' VOLUME=',VOLUMT,' NBTETR=',NBTETR,' NUDTETR=',NUDTETR
      WRITE(IMPRIM,*)

      IF( IERR .EQ. 0        .AND. ITERAT .LE. 3  .AND.
     %    NBSOMM .LT. MXSOMM .AND. NBSOM0 .LT. NBSOMM ) THEN

C        UNE ITERATION D'AJOUT DE POINTS SUR L'ARETE MAX EN PLUS
C        -------------------------------------------------------
         NUDTET0 = NUDTETR
         GOTO 10

      ENDIF


C     FIN D'AJOUT DES POINTS DES ARETES MAX TROP LONGUES
C     ==================================================
      WRITE(IMPRIM,*)
      WRITE(IMPRIM,*)'teajbapa: FIN des ITERATIONS:', NBSOMM-NBSOM00,
     %' POINTS AJOUTES et TETRAEDRISES avec QUAMOY=',QUAMOY,
     %' QUAMIN=',QUAMIN,' NB2T3T=',NB2T3T,' NBMTNT=',NBMTNT,
     %' VOLUME=',VOLUMT,' NBTETR=',NBTETR,' NUDTETR=',NUDTETR

      IERR = 0
      GOTO 9999

 9990 PRINT *,'teajbapa: TABLEAU LAPILE SATURE  MXPILE=',MXPILE
      IERR = 3

 9999 RETURN
      END
