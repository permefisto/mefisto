      SUBROUTINE DFTOP5( NTLXOB, NDIMEN, NOINTE,
     %                   NBELEM, NTNPEF,
     %                   NBSOM , MNSOMM,
     %                   NTXYZN, MNXYZN, NTXYZP, MNXYZP,
     %                   NDPGST, IERR   )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    AJOUTER LE NUMERO DES NOEUDS NON SOMMETS
C -----    AJOUTER LE NUMERO DES POINTS DES ELEMENTS FINIS
C          CALCULER LES COORDONNEES DES POINTS ET NOEUDS
C          S'ILS DIFFERENT DES SOMMETS
C ENTREES :
C ---------
C NTLXOB : NUMERO DU TMS LEXIQUE DE L'OBJET
C NDIMEN : DIMENSION DE L'ESPACE DE TRAVAIL 6 ou 3 ou 2
C NOINTE : NUMERO DE L'INTERPOLATION CHOISIE  (CF DFTOP0)
C          1  'AXISYMETRIQUE_DEGRE_1',
C          2  'AXISYMETRIQUE_DEGRE_2',
C          3  'LAGRANGE_DEGRE_1',
C          4  'LAGRANGE_DEGRE_2'
C NBELEM : NOMBRE D'ELEMENTS FINIS SELON LEUR TYPE GEOMETRIQUE ( 1 A 9 )
C NTNPEF : NUMERO  DES TMS DES TABLEAUX NPEF"TYPE_EF
C NBSOM  : NOMBRE  TOTAL DE SOMMETS DE L'OBJET
C MNSOMM : ADRESSE MCN DU TABLEAU XYZSOMMET DES SOMMETS DE L'OBJET

C SORTIES :
C ---------
C NTXYZN : NUMERO  TMS DU TABLEAU XYZNOEUD DES NOEUDS DE L'OBJET
C          0 SI TABLEAU NON CREE
C MNXYZN : ADRESSE MCN DU TABLEAU XYZNOEUD DES NOEUDS DE L'OBJET
C          0 SI TABLEAU NON CREE
C NTXYZP : NUMERO  TMS DU TABLEAU XYZPOINT DES POINTS GEOMETRIQUES DE L'OBJET
C          0 SI TABLEAU NON CREE
C MNXYZP : ADRESSE MCN DU TABLEAU XYZPOINT DES POINTS GEOMETRIQUES DE L'OBJET
C          0 SI TABLEAU NON CREE
C NDPGST : CODE GENERAL DE TRAITEMENT DES
C          NOEUDS POINTS GEOMETRIQUES ET SOMMETS
C          0 : NOEUDS=POINTS=SOMMETS
C              LE  TABLEAU  XYZSOMMET EXISTE
C          1 : NOEUDS=POINTS#SOMMETS
C              LES TABLEAUX XYZNOEUD ET XYZSOMMET EXISTENT
C          2 : NOEUDS#POINTS=SOMMETS
C              LES TABLEAUX XYZNOEUD ET XYZSOMMET EXISTENT
C          3 : NOEUDS#POINTS#SOMMETS
C              LES TABLEAUX XYZNOEUD XYZPOINT XYZSOMMET EXISTENT
C IERR   : 0 SI PAS D'ERREUR, >0 SINON
C          3 SI TABLEAU LFACES NON ALLOUABLE FAUTE DE PLACE MEMOIRE
C          4 SI LISTE DE HACHAGE SATUREE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE UPMC PARIS    JANVIER 1997
C MODIF  : ALAIN PERRONNET  TEXAS A & M UNIVERSITY            JULY  2005
C2345X7..............................................................012
      PARAMETER        (MAXNOE=24, MAXPOE=24)
      PARAMETER        (NBTROU=128)
      include"./incl/gsmenu.inc"
      include"./incl/langue.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      include"./incl/donele.inc"
      include"./incl/a___npef.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___xyznoeud.inc"
      include"./incl/a___xyzpoint.inc"
      include"./incl/ponoel.inc"
      include"./incl/inteel.inc"
      INTEGER           NTNPEF(9)
      INTEGER           NBELEM(9)
      INTEGER           NGS(4), NGS1(4), NONOF(MAXNOE),
     %                  NOSOEL(9), NOTGEL(24), NOTGAR(2)
      REAL              RMCN(1), COPOAR(2), XYZ(4,3)
      EQUIVALENCE      (RMCN(1),MCN(1))
      DOUBLE PRECISION  DCOORN(6,MAXNOE), FBASE(MAXNOE), COOR1C(MAXPOE)
      DOUBLE PRECISION  PROSCD
      INTEGER, allocatable, dimension(:,:) :: LARETES, LFACES
      INTRINSIC         INT, REAL
      include"./incl/nomele.inc"

      IERR = 0
      IF( NDIMEN .EQ. 6 ) THEN
C        TRAITEMENT DES 6-CUBES 6Q1C => NOEUDS=POINTS=SOMMETS
         NDPGST = 0
         NTXYZN = 0
         MNXYZN = 0
         NTXYZP = 0
         MNXYZP = 0
         RETURN
      ENDIF

      LFN = 0
      LAN = 0
      MXPOFA = 0
      MXPOAR = 0
      MXNOFA = 0
      MXNOAR = 0
      KVECSO = 0
      IALARETES = 1
      IALFACES  = 1
      NOFACE = 0
      NOARET = 0
      NBP = 0
      NBN = 0

C     LES TYPES D'ELEMENTS FINIS POUR DEFINIR LES NOEUDS ET POINTS
C     ============================================================
      N = 0
      DO 50 NCOGEL = 1, 9
         IF( NCOGEL .NE. 8 .AND. NTNPEF(NCOGEL) .GT. 0 ) THEN
C
C           IL EXISTE DES ELEMENTS DE CODE GEOMETRIQUE NCOGEL
C           LE NUMERO DE TYPE DE L'ELEMENT FINI POUR L'INTERPOLATION CHOISIE
            NUTYEL = NOTYEL( NCOGEL, NOINTE )
C
C           TRAITEMENT PARTICULIER DU TRIANGLE P1 A CAUSE DES ESTIMATEURS
C           D'ERREUR AVEC 2 POINTS SUR L'ARETE A RENDRE COMPATIBLE AVEC
C           LE QUADRANGLE 2Q1C
            IF( NCOGEL .EQ. 3 .AND. NUTYEL .EQ. 13 .AND.
     %          NTNPEF(4) .GT. 0 ) THEN
C               TRIANGLE 2P1D EN PRESENCE DE QUADRANGLES
C               => LE TRIANGLE DEVIENT TRIA 2P1C
C                  AVEC 2 POINTS D'INTEGRATION SUR CHAQUE ARETE
               NUTYEL = 29
            ENDIF
C
            CALL ELTYCA( NUTYEL )
C           RESULTATS DANS ponoel.inc
C           NBPOE  : NOMBRE DE POINTS DE L ELEMENT FINI NOMELE
C           NBNOE  : NOMBRE DE NOEUDS DE L ELEMENT FINI NOMELE
C           NOTRAE : CODE TRAITEMENT  DE L ELEMENT FINI NOMELE
C                    0 : NOEUDS=POINTS=SOMMETS
C                    1 : NOEUDS=POINTS#SOMMETS
C                    2 : NOEUDS#POINTS=SOMMETS
C                    3 : NOEUDS#POINTS#SOMMETS
C           NBNSOM : NOMBRE DE NOEUDS-SOMMETS  DE L ELEMENT FINI
C           NBPOIN : NOMBRE DE POINTS INTERNES DE L ELEMENT FINI
C           NBNOIN : NOMBRE DE NOEUDS INTERNES DE L ELEMENT FINI
C           NARET  : NOMBRE DE SES ARETES
C           NOSOAR : NO DES 2 SOMMETS DE CHACUNE DE SES ARETES
C           NOTYAR : NO DU TYPE       DE CHACUNE DE SES ARETES
C                    ( CF SP TYARCP )
C           NBPOAR : NOMBRE DE POINTS-NON SOMMETS DE CHACUNE DE SES ARETES
C           NOPOAR : NO ELEMENTAIRE DES POINTS NON SOMMETS DE CHACUNE DE SES
C                    ARETES
C           NBNOAR : NOMBRE DE NOEUDS-NON SOMMETS DE CHACUNE DE SES ARETES
C           NONOAR : NO ELEMENTAIRE DES NOEUDS NON SOMMETS DE CHACUNE DE SES
C                    ARETES
C           NFACE  : NOMBRE DE SES FACES
C           NBSOFA : NOMBRE DE SOMMETS DE CHACUNE DE SES FACES
C           NOSOFA : NO ELEMENTAIRE DES SOMMETS DE CHACUNE DE SES FACES
C           NOTYFA : NO DU TYPE       DE CHACUNE DE SES FACES
C           NBPOFA : NOMBRE DE POINTS-NON SUR LES ARETES DE CHACUNE DE SES FACES
C           NBNOFA : NOMBRE DE NOEUDS-NON SUR LES ARETES DE CHACUNE DE SES FACES
C           NONOFA : NO ELEMENTAIRE DES NOEUDS-NON SUR LES ARETES DE CHACUNE DE
C                    SES FACES
C
C           LES COORDONNEES DCOORN DES NOEUDS SUR L'ELEMENT DE REFERENCE
            CALL ELCONO( NUTYEL, NDIM, NBNOE, DCOORN )
C
            IF( NDIM .LE. 2 ) NFACE = 0
            IF( N    .EQ. 0 ) THEN
               NDPGST = NOTRAE
               KVECSO = NBNSOM
               MXPOAR = 0
               MXNOAR = 0
               MXPOFA = 0
               MXNOFA = 0
               N      = 1
            ENDIF
C
            IF(.NOT.( (NBNSOM .EQ. 0 .AND. KVECSO .EQ. 0) .OR.
     &         (NBNSOM .NE. 0 .AND. KVECSO .NE. 0))) THEN
               NBLGRC(NRERR) = 2
               IF( LANGAG .EQ. 0 ) THEN
                  KERR(1) = 'EF '//NOMELE(1,NUTYEL)//NOMELE(2,NUTYEL)
                  KERR(2) = 'INCOMPATIBLE AVEC LES PRECEDENTS EF'
               ELSE
                  KERR(1) = 'FE '//NOMELE(1,NUTYEL)//NOMELE(2,NUTYEL)
                  KERR(2) = 'UNCOMPATIBLE WITH THE PREVIOUS FE'
               ENDIF
               CALL LEREUR
               IERR = 2
            ENDIF
C
            IF( NARET .GT. 0 ) THEN
C              MXPOAR : MAXIMUM DE POINTS NON SOMMETS D UNE ARETE
C              MXNOAR : MAXIMUM DE NOEUDS NON SOMMETS D UNE ARETE
               DO 10 I=1,NARET
                  MXPOAR = MAX0(MXPOAR, NBPOAR(I))
                  MXNOAR = MAX0(MXNOAR, NBNOAR(I))
  10           ENDDO
            ENDIF
            IF( NFACE .GT. 0 ) THEN
C              MXPOFA : MAXIMUM DE POINTS NON SOMMETS D UNE FACE
C              MXNOFA : MAXIMUM DE NOEUDS NON SOMMETS D UNE FACE
               DO 20 I=1,NFACE
                  MXPOFA = MAX0(MXPOFA, NBPOFA(I))
                  MXNOFA = MAX0(MXNOFA, NBNOFA(I))
  20           ENDDO
            ENDIF
C
            IF( NDIMEN .LT. NDIM ) THEN
               NBLGRC(NRERR) = 3
               IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'DFTOP5: DIMENSIONS INCOMPATIBLES DES ESPACES'
               KERR(2) = 'EF '//NOMELE(1,NUTYEL)//NOMELE(2,NUTYEL)
               ELSE
               KERR(1) = 'DFTOP5: UNCOMPATIBLE SPACE DIMENSIONS'
               KERR(2) = 'FE '//NOMELE(1,NUTYEL)//NOMELE(2,NUTYEL)
               ENDIF
               CALL LEREUR
               IERR = 3
            ENDIF
C
            IF(NOTRAE * NDPGST .EQ. 6) NDPGST = 3
C           MISE A JOUR DE NDPGST
            NDPGST = MAX0(NDPGST,NOTRAE)
         ENDIF
  50  ENDDO
C
      KVECSO = MIN0(1,KVECSO)
      IF( IERR .GT. 0 ) RETURN
      IF( NDPGST .EQ. 0 ) THEN
C        NOEUDS=POINTS=SOMMETS
         NTXYZN = 0
         MNXYZN = 0
         NTXYZP = 0
         MNXYZP = 0
         RETURN
      ENDIF

C     ADRESSAGE DES TABLEAUX DES ARETES DE TOUS LES EF
C     ================================================
      LARETE = 0
      IF(MXPOAR + MXNOAR  .GT. 0) THEN
         IF( NDIM .EQ. 1 ) THEN
C           LES ARETES EN 1D
            LARETE = NBELEM(2)
         ELSE IF( NDIM .EQ. 2 ) THEN
C           LES ARETES EN 2D
C           FORMULE EXACTE EN 2-D : LARETE=NS+NTRI+NQUA+NBTROU-1
            DO I=2,4
               LARETE = LARETE + NBELEM(I)
            ENDDO
            LARETE = LARETE + NBSOM + NBTROU
         ELSE IF( NDIM .EQ. 3 ) THEN
C           LES ARETES EN 3D
            LARETE = LARETE + 3*(NBELEM(5) + NBELEM(6))
     %                      + 4*(NBELEM(7) + NBELEM(9))
C           HEXAEDRE AVEC a b c ARETES EN X Y Z
C           NB ARETES = 3 abc + 2 ab + 2 bc + 2 ca + a + b + c
C           NB HEXA   = abc
C           SI PETIT MAILLAGE, LA VALEUR CI DESSUS PEUT ETRE TROP FAIBLE
C           a=b=c=10 => 3630 aretes et la formule ci dessus est une majoration
            IF( LARETE .LE. 3630 ) LARETE=3630
         ENDIF

C        AUGMENTER POUR ASSURER LA NON SATURATION DU HACHAGE DES ARETES
         LARETE = NINT( 1.1 * LARETE )
C
         LAP = MIN( MXPOAR, 1 )
         LAN = MIN( MXNOAR, 1 )
         IF( NDPGST .LE. 1 ) LAN = 0
         LAPN = 3 + 2 * ( LAN + LAP )
         PRINT*,'dftop5: DEMANDE  ALLOCATION LARETES(',LAPN,',',LARETE,
     %          ') entiers'
         ALLOCATE ( LARETES( 1:LAPN, 1:LARETE ), STAT=IALARETES )
         IF( IALARETES .NE. 0 ) THEN
            IERR = 2
            GOTO 9999
         ENDIF
         PRINT*,'dftop5: CORRECTE ALLOCATION LARETES(',LAPN,',',LARETE,
     %          ') entiers'

         LAN = 3 + 2 * LAP
         DO I=1,LARETE
            LARETES( 1, I ) = 0
            LARETES( 3, I ) = 0
         ENDDO
      ENDIF

C     ADRESSAGE DES TABLEAUX DES FACES DE TOUS LES EF
C     ===============================================
      LFACE  = 0
      IF(NDIM .GT. 2) THEN
         IF(MXPOFA + MXNOFA  .GT. 0) THEN
            LFACE  = 1 + NBELEM(3) + NBELEM(4)
     %             + 3 * NBELEM(5) + 4 * NBELEM(6)
     %             + 5 * (NBELEM(7) + NBELEM(9))
C           SI PETIT MAILLAGE, LA VALEUR CI DESSUS PEUT ETRE TROP FAIBLE
            IF( LFACE .LE. 100 ) LFACE=INT( LFACE*1.5 )
            LFP    = MIN0( MXPOFA, 1 )
            LFN    = MIN0( MXNOFA, 1 )
            LFPN   = 5 + (LFP + LFN) * 2
            LFN    = 5 + LFP * 2
            PRINT*,'dftop5: DEMANDE  ALLOCATION LFACES(',LFPN,',',LFACE,
     %             ') entiers'
            ALLOCATE ( LFACES( 1:LFPN, 1:LFACE ), STAT=IALFACES )
            IF( IALFACES .NE. 0 ) THEN
               IERR = 3
               GOTO 9999
            ENDIF
            PRINT*,'dftop5: CORRECTE ALLOCATION LFACES(',LFPN,',',LFACE,
     %             ') entiers'
            DO I=1,LFACE
               LFACES( 1, I ) = 0
               LFACES( 5, I ) = 0
            ENDDO

         ENDIF
      ENDIF

C     DECLARATION DU TABLEAU DES XYZ DES NOEUDS ET/OU POINTS SI NECESSAIRE
C     ====================================================================
C     NOEUDS=POINTS#SOMMETS
C     NOEUDS#POINTS=SOMMETS
      LNOEUD = KVECSO * NBSOM + LARETE * MXNOAR + LFACE * MXNOFA +
     %         NBNOIN * ( NBELEM(2) + NBELEM(3) + NBELEM(4) +
     %                    NBELEM(5) + NBELEM(6) + NBELEM(7) + NBELEM(9))
      MNXYZN = 0
      IF( LNOEUD .GT. 0 ) THEN
         CALL TNMCDC( 'ENTIER', 3*LNOEUD, MNXYZN )
      ENDIF

      IF( NDPGST .EQ. 3 ) THEN
C        NOEUDS#POINTS#SOMMETS
         LPOINT = NBSOM + LARETE * MXPOAR + LFACE * MXPOFA +
     %            NBNOIN * ( NBELEM(2) + NBELEM(3) + NBELEM(4) +
     %                       NBELEM(5) + NBELEM(6) +
     %                       NBELEM(7) + NBELEM(9) )
      ELSE
         LPOINT = 0
      ENDIF
      MNXYZP = 0
      IF( LPOINT .GT. 0 ) THEN
         CALL TNMCDC( 'ENTIER', 3*LPOINT, MNXYZP )
      ENDIF

      IF( NDPGST .EQ. 1 ) THEN
         LPOINT = LNOEUD
         MNXYZP = MNXYZN
      ENDIF

C     COPIE DES NDIM COORDONNEES DES SOMMETS DANS LES NOEUDS
      NBCOOR = MCN(MNSOMM+WBCOOR)
      IF( MNXYZN .GT. 0 ) THEN
         N = MAX( NDIM, 3 )
         CALL TRTATA( MCN(MNSOMM+WYZSOM), MCN(MNXYZN), N*NBSOM )
         IF( NDIM .LE. 2 ) THEN
C           LA 3-EME COMPOSANTE EST MISE A ZERO
            DO I = MNXYZN+2, MNXYZN+N*LPOINT, N
               RMCN(I) = 0
            ENDDO
         ENDIF
         IF( NDIM .EQ. 1 ) THEN
C           LA 2-EME COMPOSANTE EST MISE A ZERO
            DO I = MNXYZN+1, MNXYZN+N*LPOINT, N
               RMCN(I) = 0
            ENDDO
         ENDIF
      ENDIF

C     FORMATION DU TABLEAU DES XYZ DES NOEUDS EVENTUELLEMENT DES POINTS
C     =================================================================
      IF(KVECSO .NE. 0) THEN
C        LES SOMMETS SONT DES NOEUDS
         NONO = NBSOM
      ELSE
C        LES SOMMETS NE SONT PAS DES NOEUDS
         NONO = 0
      ENDIF

      NOPO = NBSOM

C     ADRESSE DANS MCN DES XYZ DES TANGENTES
      MNXYTG = MNSOMM + WYZSOM + 3 * NBSOM
C     LE NOMBRE DE TANGENTES DE L'OBJET
      NBTGS = MCN( MNSOMM + WNBTGS )
C
C     POSITION LIBRE DANS LES TABLEAUX DE HACHAGE
      LIBREA = LARETE
      LIBREF = LFACE
C
C     BOUCLE SUR LES PYRAMIDES, HEXAEDRES, PENTAEDRES, TETRAEDRES, QUADRANGLES
C     ========================================================================
      DO 1000 NCOGEL=9,1,-1
C
C        PAS DE TRAITEMENT DES 6-CUBES
         IF( NCOGEL .EQ. 8 ) GOTO 1000
C
C        LE NUMERO DU TABLEAU TMS ~>OBJET>>NPEF"NCOGEL
         NTELE = NTNPEF(NCOGEL)
         IF( NTELE .LE. 0 ) GOTO 1000
C
C        RECUPERATION DE SON ADRESSE
         CALL TAMSOU( NTELE, MNELE )
         IF( MNELE .LE. 0 ) GOTO 1000
C
C        IL EXISTE DES ELEMENTS DE CODE GEOMETRIQUE NCOGEL
C        LE NOMBRE ACTUEL D'ELEMENTS FINIS AVEC CE CODE GEOMETRIQUE
         NBELE  = MCN( MNELE + WBELEM )
C
C        LE NUMERO DE TYPE DE L'ELEMENT FINI POUR L'INTERPOLATION CHOISIE
         NUTYEL =  NOTYEL( NCOGEL, NOINTE )
C
C        TRAITEMENT PARTICULIER DU TRIANGLE P1 A CAUSE DES ESTIMATEURS
C        D'ERREUR AVEC 2 POINTS SUR L'ARETE A RENDRE COMPATIBLE AVEC
C        LE QUADRANGLE 2Q1C
         IF( NCOGEL .EQ. 3 .AND. NUTYEL .EQ. 13 .AND.
     %       NTNPEF(4) .GT. 0 ) THEN
C            TRIANGLE 2P1D EN PRESENCE DE QUADRANGLES
C            => LE TRIANGLE DEVIENT TRIA 2P1C
C               AVEC 2 POINTS D'INTEGRATION SUR CHAQUE ARETE
            NUTYEL = 29
         ENDIF
C
C        LES CARACTERISTIQUES DE L'ELEMENT
C        ---------------------------------
         CALL ELTYCA( NUTYEL )
C        NBPOE  : NOMBRE DE POINTS DE L ELEMENT FINI NOMELE
C        NBNOE  : NOMBRE DE NOEUDS DE L ELEMENT FINI NOMELE
C        NOTRAE : CODE TRAITEMENT  DE L ELEMENT FINI NOMELE
C        NBNSOM : NOMBRE DE NOEUDS-SOMMETS  DE L ELEMENT FINI
C        NBPOIN : NOMBRE DE POINTS INTERNES DE L ELEMENT FINI
C        NBNOIN : NOMBRE DE NOEUDS INTERNES DE L ELEMENT FINI
C        NARET  : NOMBRE DE SES ARETES
C        NOSOAR : NO DES 2 SOMMETS DE CHACUNE DE SES ARETES
C        NOTYAR : NO DU TYPE       DE CHACUNE DE SES ARETES
C                 ( CF SP TYARCP )
C        NBPOAR : NOMBRE DE POINTS-NON SOMMETS DE CHACUNE DE SES ARETES
C        NOPOAR : NO ELEMENTAIRE DES POINTS NON SOMMETS DE CHACUNE DE SES
C                 ARETES
C        NBNOAR : NOMBRE DE NOEUDS-NON SOMMETS DE CHACUNE DE SES ARETES
C        NONOAR : NO ELEMENTAIRE DES NOEUDS NON SOMMETS DE CHACUNE DE SES
C                 ARETES
C        NFACE  : NOMBRE DE SES FACES
C        NBSOFA : NOMBRE DE SOMMETS DE CHACUNE DE SES FACES
C        NOSOFA : NO ELEMENTAIRE DES SOMMETS DE CHACUNE DE SES FACES
C        NOTYFA : NO DU TYPE       DE CHACUNE DE SES FACES
C        NBPOFA : NOMBRE DE POINTS-NON SUR LES ARETES DE CHACUNE DE SES FACES
C        NBNOFA : NOMBRE DE NOEUDS-NON SUR LES ARETES DE CHACUNE DE SES FACES
C        NONOFA : NO ELEMENTAIRE DES NOEUDS-NON SUR LES ARETES DE CHACUNE DE
C                 SES FACES
C
C        LE NOMBRE DE SOMMETS
         NBSO = NBSOME( NCOGEL )
C
C        LES COORDONNEES DCOORN DES NOEUDS SUR L'ELEMENT DE REFERENCE
         CALL ELCONO( NUTYEL, NDIM, NBNOE, DCOORN )
C
C        LES CARACTERISTIQUES DE L'INTERPOLATION DE F:EF REFERENCE->EF COURANT
         CALL ELINTE( 'THERMIQUE', NUTYEL, NDIMF , NUINTF,
     %                 NBINVA, NUINVA, NUINTI, NBNDIN )
C        NDIMF  : NOMBRE DE COMPOSANTES DE L APPLICATION
C                 ELEMENT REFERENCE --> ELEMENT COURANT
C        NUINTF : NO DE L INTERPOLATION COMMUNE AUX NDIMF COMPOSANTES
C                 ( CF SP INTERP )
C
C        LES VALEURS SUIVANTES INUTILES ICI SONT DONNEES POUR MEMOIRE
C        NBINVA : NOMBRE D INCONNUES VARIATIONNELLES SUR CET ELEMENT
C        NUINVA : NO DE CES NBINVA INCONNUES VARIATIONNELLES
C                 ( CF SP INCVAR )
C        NUINTI : NO D INTERPOLATION DE CES NBINVA INCONNUES VARIATIONNELLES
C                 ( CF SP INTERP )
C        NBNDIN : NOMBRE DE NOEUDS NECESSAIRES A CHAQUE INTERPOLATION
C                 DE CHAQUE INCONNUE VARIATIONNELLE
C                 -1 SIGNIFIE QUE TOUS LES NOEUDS DE L ELEMENT SERVENT
C                    ET QUE NUNOIN DE CETTE INCONNUE N EST PAS DONNE
C        NUNOIN : NO DES NBNDIN NOEUDS DANS LA NUMEROTATION ELEMENTAIRE
C                 NECESSAIRE A CHAQUE INCONNUE VARIATIONNELLE
C                 NUNOIN(I,J) J=1,...,NBINVA . I=1,...,NBNDIN(J) SI NON = -1
C
C        L'ADRESSE DU DEBUT DU TABLEAU DES NUMEROS DES NOEUDS CONTENANT LES SEUL
         MNNE   = MNELE + WUNDEL - 1
         NBNDEL = MCN( MNELE + WBNDEL )
         IF( NDPGST .LE. 1 ) THEN
C           NOEUDS = POINTS
C           ADRESSE DES NUMEROS DES POINTS DANS MCN
            MNPE = MNNE
         ELSE
C           NOEUDS # POINTS
C           ADRESSE DES NUMEROS DES POINTS DANS MCN
            MNPE = MNNE  + NBNDEL * NBELE
         ENDIF
C
         DO 900 NEL=1,NBELE
C           LE PREMIER SOMMET ET POINT DE L'ELEMENT FINI
            MNNE = MNNE + 1
            MNPE = MNPE + 1
C
C           LE NUMERO DES SOMMETS DE L'ELEMENT FINI
            DO 590 I=0,NBSO-1
               NOSOEL(I+1) = MCN( MNNE + I * NBELE )
 590        ENDDO
C
C           LE NOMBRE ET NUMEROS DES EVENTUELLES TANGENTES DE CET EF
            CALL EFTGEF( MNELE,  NEL,
     %                   NBTGEL, NOTGEL )
C
C           TRAITEMENT DES POINTS, NOEUDS NON SOMMETS DES ARETES
C           ----------------------------------------------------
            IF( NARET .EQ. 0 ) GOTO 700
            MNCP = MNXYZN - 1
            DO 690 I=1,NARET
C
C              LE NOMBRE DE POINTS ET NOEUDS NON SOMMETS DE L'ELEMENT FINI
               NBP = NBPOAR(I)
               NBN = NBNOAR(I)
               IF(NBP+NBN .EQ. 0) GOTO 690
C              LE NO GLOBAL DES POINTS-SOMMETS DE L ARETE
               NGS(1) = NOSOEL( NOSOAR(1,I) )
               NGS(2) = NOSOEL( NOSOAR(2,I) )
C
C              REORIENTATION DE L ARETE AFIN QUE LE 1-ER SOMMET AIT LE
C              LE PLUS PETIT NUMERO
               LESENS = 1
               IF( NGS(1) .GT. NGS(2) ) THEN
                  LESENS = -1
                  J      = NGS(1)
                  NGS(1) = NGS(2)
                  NGS(2) = J
               ENDIF
C
C              RECHERCHE DE L ARETE OU ADJONCTION PAR HACHAGE
C              ----------------------------------------------
               CALL HACHAG( 2,NGS,LAPN,LARETE,LARETES,3,LIBREA,
     %                      NOARET )
               IF( NOARET .EQ. 0 ) THEN
C
C                 LISTE DE HACHAGE SATUREE
                  IERR = 4
                  RETURN
C
               ELSE IF( NOARET .LT. 0 ) THEN
C
C                 L'ARETE VIENT D'ETRE CREEE : CALCUL DES
C                 COORDONNEES DES NBP POINTS NON-SOMMETS DE L ARETE
C                 -------------------------------------------------
                  NOARET = -NOARET
                  IF( NBP .GT. 0 ) THEN
                     CALL TYARCP( NOTYAR(I), NBP, COPOAR )
C                    COPOAR COORDONNEE CURVILIGNE NORMALISEE A 1.
C                    DES NBP POINTS DE L ARETE
                     L  = MNCP + 3 * NOPO
                     L1 = MNCP + 3 * NGS(1) - 3
                     L2 = MNCP + 3 * NGS(2) - 3
                     IF( NBTGEL .GT. 0 ) THEN
C
C                       EF AVEC TG: LES 2 TANGENTES DE L'ARETE P3 HERMITE
C                       .................................................
                        CALL TGAREF( NCOGEL, I, NOTGAR )
C                       LE NUMERO GLOBAL DE LA TANGENTE
                        NOTGAR(1) = NOTGEL( NOTGAR(1) )
                        NOTGAR(2) = NOTGEL( NOTGAR(2) )
                        IF(NOTGAR(1).EQ.0 .AND. NOTGAR(2).EQ.0) GOTO 630
C
C                       ARETE P3 HERMITE: POINT SUR L'ARETE COURBE
C                       CONSTRUCTION DES TABLEAUX X Y Z POUR LE PARAMETRAGE
C                       DE L'ARETE P3 HERMITE: XS1 XS2  XTGS1S2 XTGS2S1
                        DO 601 K=1,3
                          XYZ(1,K) = RMCN(L1+K)
                          XYZ(2,K) = RMCN(L2+K)
  601                   ENDDO
C
                        LSIGNE = -1
                        DO 608 M=1,2
C                          CALCUL DES 3 COMPOSANTES DE LA TANGENTE AU SOMMET M D
                           LSIGNE = -LSIGNE
                           IF( NOTGAR(M) .EQ. 0 ) THEN
C                            TANGENTE SELON L'ARETE DROITE
                             DO 602 K=1,3
                               XYZ(2+M,K)=LSIGNE*(RMCN(L2+K)-RMCN(L1+K))
  602                        ENDDO
                           ELSE IF( NOTGAR(M) .LT. 0 ) THEN
C                             TANGENTE OPPOSEE A CELLE STOCKEE
                              DO 604 K=1,3
                                XYZ(2+M,K)=-RMCN(MNXYTG-4-3*NOTGAR(M)+K)
  604                         ENDDO
                           ELSE
C                             TANGENTE STOCKEE
                              DO 606 K=1,3
                                 XYZ(2+M,K)=RMCN(MNXYTG-4+3*NOTGAR(M)+K)
  606                         ENDDO
                           ENDIF
  608                   ENDDO
C
                        DO 620 J=1,NBP
C                          LES 3 COORDONNEES DU POINT SUR L'ARETE P3-HERMITE
                           CALL XYZP3H( COPOAR(J),
     %                                  XYZ(1,1), XYZ(1,2), XYZ(1,3),
     %                                  RMCN(L+1),RMCN(L+2),RMCN(L+3) )
                           L = L + 3
  620                   ENDDO
C
C                       LES 3 COORDONNEES DES POINTS NON SOMMETS DE L'EF SONT CA
                        GOTO 650
                     ENDIF
C
C                    EF SANS TG : POINT SUR L'ARETE DROITE
C                    .....................................
  630                DO 640 J=1,NBP
                        DO 635 K=1,3
                           RMCN(L+K) = RMCN(L1+K) * (1-COPOAR(J))
     %                               + RMCN(L2+K) * COPOAR(J)
  635                   ENDDO
                        L = L + 3
  640                ENDDO
C
C                    ADJONCTION DE L ARETE.MISE A JOUR DES POINTS DE L ARETE
C                    -------------------------------------------------------
  650                LARETES( 4, NOARET ) = NOPO + 1
                     NOPO = NOPO + NBP
                     LARETES( 5, NOARET ) = NOPO

                  ENDIF
C
C                 LES NOEUDS.STOCKAGE DU NO DU 1-ER ET DERNIER DE L ARETE
C                 -------------------------------------------------------
                  IF( NDPGST .GT. 1 ) THEN
                     IF( NBN .GT. 0 ) THEN
                        LARETES( 1+LAN, NOARET ) = NONO + 1
                        NONO = NONO + NBN
                        LARETES( 2+LAN, NOARET ) = NONO
                     ENDIF
                  ENDIF
               ENDIF
C
C              COPIE DU NO DES POINTS NON SOMMETS DE L ARETE DANS NUPGEL
C              ---------------------------------------------------------
               IF(NBP .GT. 0) THEN
                  NO = LARETES( 4, NOARET ) -1
C                 NO = NUMERO DU 1-ER POINT - 1 DE L ARETE
                  DO 670 J=1,NBP
C                    NO ELEMENTAIRE DU J-EME POINT DE L ARETE I
                     I1 = NOPOAR(J,I)
C                    NO DANS NUPGEL
                     I1 = MNPE + (I1-1) * NBELE
                     IF(LESENS .GE. 0) THEN
C                       ARETE DANS LE BON SENS
                        MCN(I1) = NO + J
                     ELSE
C                       ARETE DANS LE SENS INVERSE
                        MCN(I1) = NO + NBP - J + 1
                     ENDIF
  670             ENDDO
               ENDIF
C
C              LE NUMERO DES NOEUDS NON SOMMETS DE L'ARETE
C              -------------------------------------------
               IF(NDPGST .GT. 1) THEN
                  IF(NBN .GT. 0) THEN
                     NO = LARETES( 1+LAN, NOARET ) -1
C                    LE NO DU 1-ER NOEUD - 1
                     DO 685 J=1,NBN
C                       LE NO ELEMENTAIRE DU J-EME NOEUD DE L'ARETE I
                        I1 = NONOAR(J,I)
C                       L'ADRESSE DU NOEUD DANS MCN
                        I1 = MNNE + (I1-1) * NBELE
                        IF( LESENS .GE. 0 ) THEN
                           MCN( I1 ) = NO + J
                        ELSE
                           MCN( I1 ) = NO + NBN - J + 1
                        ENDIF
  685                ENDDO
                  ENDIF
               ENDIF
  690       ENDDO
C
C           TRAITEMENT DES FACES
C           ====================
  700       IF( NFACE .EQ. 0  .OR. NDIM .LE. 2 ) GOTO 800
            DO 790 I=1,NFACE
C
C              LE NOMBRE DE NOEUDS ET POINTS NON SUR LES ARETES DE LA FACE
               NBP = NBPOFA(I)
               NBN = NBNOFA(I)
               IF(NBP+NBN .EQ. 0) GOTO 790
               NS = NBSOFA(I)
C              LE NO GLOBAL DES POINTS SOMMETS DE LA FACE I
               DO 705 J=1,NS
                  NGS(J) = NOSOEL( NOSOFA(J,I) )
  705          ENDDO
C
C              PERMUTATION CIRCULAIRE DES SOMMETS POUR AMENER
C              LE PLUS PETIT NO EN PREMIER
               LESENS = 1
               L      = NGS(1)
               DO 708 J=2,NS
                  IF(L .GT. NGS(J))THEN
                     L      = NGS(J)
                     LESENS = J
                  ENDIF
  708          ENDDO
               CALL TRTATA(NGS(LESENS) , NGS1(1) , NS-LESENS+1)
               IF(LESENS .EQ. 1) GOTO 710
               CALL TRTATA(NGS(1) , NGS1(NS-LESENS+2) , LESENS-1)
C
C              LA NUMEROTATION DES SOMMETS EST MISE A JOUR POUR
C              QUE LA 1-ERE ARETE DE LA FACE SOIT
C              1:NO LE PLUS FAIBLE DES SOMMETS DE LA FACE
C              2:NO MIN(SOMMET2 , SOMMET NS)
C              SI MIN =SOMMET2   FACE DIRECTE
C                      SOMMET NS FACE INDIRECTE
C
  710          IF(NGS1(2) .GE. NGS1(NS)) THEN
C                 FACE INDIRECTE
                  LESENS   = - LESENS
                  L        = NGS1(2)
                  NGS1(2)  = NGS1(NS)
                  NGS1(NS) = L
               ENDIF
C
C              RECHERCHE DE LA FACE DANS LE TABLEAU LFACE PAR HACHAGE
C              ------------------------------------------------------
               CALL HACHAG( NS, NGS1, LFPN, LFACE, LFACES, 5,
     &                      LIBREF, NOFACE )

               IF( NOFACE .LE. 0 ) THEN
C                 NOUVELLE FACE. MISE A JOUR DE LA LISTE DES FACES
                  NOFACE = -NOFACE
C
C                 LES POINTS NON SUR LES ARETES DE LA FACE
C                 ----------------------------------------
                  IF(NBP .GT. 0) THEN
                     LFACES( 6, NOFACE ) = NOPO + 1
                     NOPO = NOPO + NBP
                     LFACES( 7, NOFACE ) = NOPO
                  ENDIF
C
C                 LES NOEUDS NON SUR LES ARETES DE LA FACE
C                 ----------------------------------------
                  IF(NDPGST .GT. 1) THEN
                     IF(NBN .GT. 0) THEN
                        LFACES( 1 + LFN, NOFACE ) = NONO + 1
                        NONO = NONO + NBN
                        LFACES( 2 + LFN, NOFACE ) = NONO
                     ENDIF
                  ENDIF
               ENDIF
C
C              CALCUL DES COORDONNEES DES NBP POINTS DE LA FACE I : A FAIRE
C              POUR L'INSTANT PAS DE CALCUL DES COORDONNEES DES POINTS
C              CAR PAS D'EF AVEC DES POINTS NON SOMMET OU MILIEU D'ARETE
C
C              RECOPIE DU NO DES POINTS DE LA FACE DANS NUPGEL
C              =============================================
C              RENUMEROTATION DES POINTS DE LA FACE DE TYPE NOTYFA(I)
C              VUE DANS LE SENS DIRECT SI NOFACE<0 INDIRECT SI NOFACE>0
C              POUR UNE PERMUTATION INITIALE DE (LESENS-1) POINTS
               CALL RENUFA( LESENS, NOTYFA(I), NBN, NONOFA(1,I), NONOF )
C
C              CE SP DONNE LE NO ELEMENTAIRE DES NBP POINTS DE LA FACE
C              DE NO GLOBAL NOPO+1,NOPO+2,...,NOPO+NBP
C
               IF(NBP .GT. 0) THEN
C
C                 LES POINTS NON SUR LES ARETES DE LA FACE
C                 ----------------------------------------
                  NO = LFACES( 6, NOFACE ) - 1
                  DO J=1,NBP
C                    NO DANS NUPGEL
                     MCN(MNPE+(NONOF(J)-1)*NBELE) = NO + J
                  ENDDO
C
               ENDIF
C
C              LES NOEUDS NON SUR LES ARETES DE LA FACE
C              ----------------------------------------
               IF(NDPGST .GT. 1) THEN
                  IF(NBN .GT. 0) THEN
                     NO = LFACES( 1+LFN, NOFACE ) - 1
                     DO 785 J=1,NBN
                        MCN(MNNE+(NONOF(J)-1)*NBELE) = NO + J
  785                ENDDO
                  ENDIF
               ENDIF
C
  790       ENDDO
C
C           LES NUMEROS DES POINTS ET NOEUDS INTERNES DE L ELEMENT FINI
C           ===========================================================
  800       IF(NBPOIN .GT. 0) THEN
C
C              LES POINTS INTERNES A L'ELEMENT FINI
C              ------------------------------------
C              LE NUMERO DU PREMIER POINT INTERNE DANS L'ELEMENT
               NO = NBPOE - NBPOIN
               DO J=1,NBPOIN
                  MCN(MNPE + (NO+J) * NBELE ) = NOPO + J
               ENDDO
               NOPO = NOPO + NBPOIN
            ENDIF
C
            IF(NDPGST .GT. 1  .AND.  NBPOIN .GT. 0) THEN
C
C              LES NOEUDS INTERNES A L'ELEMENT FINI
C              ------------------------------------
               NO = NBNOE - NBNOIN
               DO J=1,NBNOIN
                  MCN(MNNE + (NO+J) * NBELE) = NONO + J
               ENDDO
               NONO = NONO + NBNOIN
            ENDIF

C           CALCUL DES COORDONNEES DES NOEUDS DE L'ELEMENT FINI
C           ===================================================
            DO 890 J=1,NBNOE

C              LA VALEUR DES NBN FONCTIONS DE BASE DE F EN CE NOEUD J
               L = MCN( MNNE + (J-1) * NBELE )
               CALL INTERP(NUINTF, DCOORN(1,J),DCOORN(2,J),DCOORN(3,J),
     %                     NBPOF, FBASE)

               DO 880 I=0,NDIMF-1
C                 LA COMPOSANTE I DES POINTS DE L'ELEMENT FINI
                  DO K=0,NBPOF-1
                     NO = MCN( MNPE + K * NBELE )
                     COOR1C(K+1) = RMCN( MNXYZN+(NO-1)*3+I )
                  ENDDO

C                 LA COORDONNEE I+1 DU NOEUD L
                  RMCN(MNXYZN+I+3*L-3)= REAL(PROSCD(FBASE,COOR1C,NBPOF))
  880          ENDDO
  890       ENDDO

  900    ENDDO

 1000 ENDDO

C     DESTRUCTION DES TABLEAUX LARETES ET LFACES
C     ------------------------------------------
      IF( IALARETES .EQ. 0 ) THEN
         DEALLOCATE( LARETES )
         IALARETES = 1
      ENDIF

      IF( IALFACES .EQ. 0 ) THEN
         DEALLOCATE( LFACES )
         IALFACES = 1
      ENDIF
C
C     SAUVEGARDE DU TMS XYZNOEUD
C     ==========================
      IF( LNOEUD .GT. 0 ) THEN
 2000    CALL LXTSOU( NTLXOB, 'XYZNOEUD', NT, MN )
         IF( NT .GT. 0 ) THEN
C           LE TABLEAU EXISTANT EST DETRUIT
            CALL LXTSDS( NTLXOB, 'XYZNOEUD' )
            GOTO 2000
         ENDIF
C        DECLARATION ET OUVERTURE DU TABLEAU XYZNOEUD
         IF( NDPGST .EQ. 1 ) THEN
C           NOEUDS=POINTS
            NONO = NOPO
         ENDIF
         CALL LXTNDC( NTLXOB, 'XYZNOEUD', 'MOTS', WYZNOE+3*NONO )
         CALL LXTSOU( NTLXOB, 'XYZNOEUD', NT, MN )
C        COPIE DES 3 COORDONNEES DE CHAQUE NOEUD
         CALL TRTATA( MCN(MNXYZN), MCN(MN+WYZNOE), 3*NONO )
C        DESTRUCTION DU TABLEAU DEVENU INUTILE
         CALL TNMCDS( 'ENTIER', 3*LNOEUD, MNXYZN )
CCCC        COPIE DES 3 COMPOSANTES DE CHAQUE TANGENTE
CCC         CALL TRTATA( MCN(MNXYTG), MCN(MN+WYZNOE+3*NONO), 3*NBTGS )
C        LE NOMBRE DE COORDONNEES D'UN NOEUD
         MCN( MN + WBCOON ) = NBCOOR
C        LE NOMBRE DE NOEUDS
         MCN( MN + WNBNOE ) = NONO
C        PAS DE TG . ELLES SONT STOCKEES DANS XYZSOMMET
         MCN( MN + WNBTGN ) = 0
C        LA DATE
         CALL ECDATE( MCN(MN) )
C        LE NUMERO DU TABLEAU DESCRIPTEUR
         MCN( MN + MOTVAR(6) ) = NONMTD( '~>>>XYZNOEUD' )
C        MISE A JOUR FINALE DU NUMERO TMS ET ADRESSE MCN DE 'XYZNOEUD'
         NTXYZN = NT
         MNXYZN = MN
      ENDIF

C     SAUVEGARDE DU TMS XYZPOINT
C     ==========================
      IF( LPOINT .GT. 0 .AND. NDPGST .GE. 2 ) THEN
 3000    CALL LXTSOU( NTLXOB, 'XYZPOINT', NT, MN )
         IF( NT .GT. 0 ) THEN
C           LE TABLEAU EXISTANT EST DETRUIT
            CALL LXTSDS( NTLXOB, 'XYZPOINT' )
            GOTO 3000
         ENDIF
C        DECLARATION ET OUVERTURE DU TMS XYZPOINT
         CALL LXTNDC( NTLXOB, 'XYZPOINT', 'MOTS', WYZPOI+3*NOPO )
         CALL LXTSOU( NTLXOB, 'XYZPOINT', NT, MN )
C        COPIE DES 3 COORDONNEES DE CHAQUE POINT
         CALL TRTATA( MCN(MNXYZP), MCN(MN+WYZPOI), 3*NOPO )
C        DESTRUCTION DU TABLEAU DEVENU INUTILE
         CALL TNMCDS( 'ENTIER', 3*LPOINT, MNXYZP )
CCCC        COPIE DES 3 COMPOSANTES DE CHAQUE TANGENTE
CCC         CALL TRTATA( MCN(MNXYTG), MCN(MN+WYZPOI+3*NOPO), 3*NBTGS )
C        LE NOMBRE DE COORDONNEES D'UN POINT
         MCN( MN + WBCOOP ) = NBCOOR
C        LE NOMBRE DE POINTS
         MCN( MN + WNBPOI ) = NOPO
C        PAS DE TG . ELLES SONT STOCKEES DANS XYZSOMMET
         MCN( MN + WNBTGP ) = 0
C        LA DATE
         CALL ECDATE( MCN(MN) )
C        LE NUMERO DU TABLEAU DESCRIPTEUR
         MCN( MN + MOTVAR(6) ) = NONMTD( '~>>>XYZPOINT' )
      ENDIF
C
C     MISE A JOUR FINALE DU NUMERO TMS ET ADRESSE MCN DE 'XYZPOINT'
      IF( NDPGST .EQ. 1 ) THEN
C        NOEUDS=POINTS
         NTXYZP = 0
         MNXYZP = 0
      ELSE
C        NOEUDS#POINTS
         NTXYZP = NT
         MNXYZP = MN
      ENDIF
C
 9999 RETURN
      END
