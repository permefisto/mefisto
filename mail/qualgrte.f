      SUBROUTINE QUALGRTE( PTXYZD, MXTETR, NOTETR, NBGRTE, NOGRTE,
     %                     QUAMED, NBTMED,
     %                     QUAMIN, QUAMOY, VOLUMT, VOLMOY )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DES QUALITES MIN MOYENNE et VOLUMES d'un GROUPE de
C -----    TETRAEDRES de la TETRAEDRISATION

C ENTREES:
C --------
C PTXYZD : TABLEAU DES COORDONNEES DES POINTS
C          PAR POINT : X  Y  Z  DISTANCE_SOUHAITEE
C MXTETR : NOMBRE MAXIMAL DE TETRAEDRES DECLARABLES DANS NOTETR
C NOTETR : LISTE DES TETRAEDRES
C          SOMMET1,    SOMMET2,    SOMMET3,    SOMMET4,
C          TETRAEDRE1, TETRAEDRE2, TETRAEDRE3, TETRAEDRE4
C          DE L'AUTRE COTE DE LA FACE
C          1: 123      2: 234      3: 341      4: 412
C NBGRTE : NOMBRE DE TETRAEDRES DU GROUPE
C NOGRTE : NUMERO NOTETR DES TETRAEDRES DU GROUPE
C QUAMED : QUALITE MEDIOCRE POUR RECENSER LES TETRAEDRES DE QUALITE
C          AU DESSOUS DE CETTE VALEUR

C SORTIES:
C --------
C NBTMED : NOMBRE DE TETRAEDRES DE QUALITE INFERIEURE A QUAMED
C QUAMIN : QUALITE MINIMALE DES TETRAEDRES DU GROUPE
C QUAMOY : QUALITE MOYENNE  DES TETRAEDRES DU GROUPE
C VOLUMT : VOLUME  TOTAL    DES TETRAEDRES DU GROUPE
C VOLMOY : VOLUME MOYEN D'UN TETRAEDRE     DU GROUPE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Saint PIERRE du PERRAY          Novembre 2018
C2345X7..............................................................012
ccc      include"./incl/langue.inc"
ccc      COMMON / TRTETR / STOPTE, TRACTE
ccc      LOGICAL           STOPTE, TRACTE, TRACTE0
cccC     STOPTE = FAUX ==> PAS D'ARRET APRES LE TRACE DE CHAQUE TETRAEDRE
cccC              VRAI ==> DEMANDE D'UN CARACTERE POUR REDEMARRER
cccC     TRACTE = FAUX ==> PAS DE TRACE DES TETRAEDRES
cccC              VRAI ==> TRACE DES TETRAEDRES DE L'ETOILE

      DOUBLE PRECISION  PTXYZD(4,*), VOLUMT, VOLMOY, VOLNTE, VTEMIN
      DOUBLE PRECISION  ARMIN, ARMAX, SURFTR(4), ARMINM, ARMAXM
      INTEGER           NOTETR(8,MXTETR), NOGRTE(NBGRTE)
      REAL              QUAMIN, QUAMOY, QUANTE

ccc      DOUBLE PRECISION QL
ccc      CHARACTER*120     KTITRE
ccc      TRACTE0= TRACTE

      VOLUMT = 0D0
      VOLMOY = 0D0
      NBTMED = 0
      QUAMOY = 0.
      QUAMIN = 2.
      NTEMIN = 0
      ARMINM = 1D100
      ARMAXM = 0D0
      NBTE   = 0
      DO NT = 1, NBGRTE

         NTE = NOGRTE( NT )

         IF( NTE .GT. 0 .AND. NOTETR(1,NTE) .GT. 0 ) THEN

            NBTE = NBTE + 1
            CALL QUATETD( PTXYZD(1,NOTETR(1,NTE)),
     %                    PTXYZD(1,NOTETR(2,NTE)),
     %                    PTXYZD(1,NOTETR(3,NTE)),
     %                    PTXYZD(1,NOTETR(4,NTE)),
     %                    ARMIN, ARMAX, SURFTR, VOLNTE, QUANTE )

ccc            IF( VOLNTE.LE.0D0 .OR. QUANTE .LE. 0 ) THEN
ccc               PRINT*,'qualgrte: TETRAEDRE(',NTE,')=',
ccc     %                 (NOTETR(kk,NTE),kk=1,8),
ccc     %                ' V=',VOLNTE,' Q=',QUANTE
ccc               DO N=1,4
ccc                  NST = NOTETR( N, NTE )
ccc                  PRINT*,'PTXYZD(',NST,')=',(PTXYZD(kk,NST),kk=1,4)
ccc               ENDDO
ccc               print*
ccc            ENDIF

            VOLUMT = VOLUMT + VOLNTE
            QUAMOY = QUAMOY + QUANTE

            IF( QUANTE .LT. QUAMIN ) THEN
               QUAMIN = QUANTE
               NTEMIN = NTE
               VTEMIN = VOLNTE
            ENDIF

            IF( ARMIN .LT. ARMINM ) THEN
               ARMINM = ARMIN
               NTEAMI = NTE
            ENDIF

            IF( ARMAX .GT. ARMAXM ) THEN
               ARMAXM = ARMAX
               NTEAMX = NTE
            ENDIF

            IF( QUANTE .LE. QUAMED ) THEN
               NBTMED = NBTMED + 1

ccc               TRACTE  = .TRUE.
ccc               KTITRE='TETRAEDRE               QUALITE=                 
ccc     %VOLUME=             '
ccc               WRITE( KTITRE(11:20),'(I10)'   ) NTE
ccc               WRITE( KTITRE(33:47),'(G15.6)' ) QUANTE
ccc               WRITE( KTITRE(57:71),'(G15.6)' ) VOLNTE
ccc               CALL SANSDBL( KTITRE, NBC )
ccc               PRINT*,'qualgrte: ', KTITRE(1:NBC)
ccc               CALL TRFETO9( KTITRE(1:NBC), PTXYZD, NTE, NOTETR )
ccc               TRACTE = TRACTE0

            ENDIF

         ENDIF

      ENDDO

      IF( NBTE .GT. 0 ) THEN
C        VOLUME MOYEN D'UN TETRAEDRE DU GROUPE DE TETRAEDRES
         VOLMOY = VOLUMT / NBTE

C        QUALITE MOYENNE DU GROUPE DE NBGRTE TETRAEDRES
         QUAMOY = QUAMOY / NBTE
      ENDIF


cccC     QUALITE MIN EN DOUBLE PRECISION POUR AFFICHAGE
ccc      QL = QUAMIN

ccc      IF( LANGAG .EQ. 0 ) THEN
ccc         PRINT*'qualgrte:',NBTE,
ccc     %                  ' TETRAEDRES de VOLUME',VOLUMT,
ccc     %                  ' VOLUME MOYEN 1 TETRAEDRE=',VOLMOY,
ccc     %                  ' QUALITE MOYENNE=',QUAMOY,
ccc     %                ' Nombre TETRAEDRES de QUALITE<',QUAMED,'=',NBTMED
ccc         PRINT*'QUALITE MINIMALE=',QL,
ccc     %    '          du TETRAEDRE',NTEMIN,':',(NOTETR(K,NTEMIN),K=1,8)
ccc         PRINT*'ARETE MINIMALE  =',ARMINM,
ccc     %             ' du TETRAEDRE',NTEAMI,':',(NOTETR(K,NTEAMI),K=1,8)
ccc         PRINT*'ARETE MAXIMALE  =',ARMAXM,
ccc     %             ' du TETRAEDRE',NTEAMX,':',(NOTETR(K,NTEAMX),K=1,8)
ccc      ELSE
ccc         PRINT*'qualgrte:',NBTE,
ccc     %                  ' TETRAHEDRA VOLUME',VOLUMT,
ccc     %                  ' TETRAHEDRON MEAN VOLUME=',VOLMOY,
ccc     %                  ' MEAN QUALITY=',QUAMOY,
ccc     %                  ' QUALITY<',QUAMED,' TETRAHEDRA NUMBER=',NBTMED
ccc         PRINT*'MINIMUM QUALITY=',QL,
ccc     %    '          of TETRAHEDRON',NTEMIN,':',(NOTETR(K,NTEMIN),K=1,8)
ccc         PRINT*'MINIMUM EDGE   =',ARMINM,
ccc     %             ' of TETRAHEDRON',NTEAMI,':',(NOTETR(K,NTEAMI),K=1,8)
ccc         PRINT*'MAXIMUM EDGE   =',ARMAXM,
ccc     %             ' of TETRAHEDRON',NTEAMX,':',(NOTETR(K,NTEAMX),K=1,8)
ccc      ENDIF

cccC     TRACE DU TETRAEDRE DE PLUS MAUVAISE QUALITE ET SES VOISINS
ccc      TRACTE  = .TRUE.
ccc      KTITRE='qualgrte FIN: TETRAEDRE               de QUAMIN=          
ccc     %       VOLUME=                 VOLUME/VOLMOY=                 '
ccc      WRITE( KTITRE(25:34),'(I10)'   ) NTEMIN
ccc      WRITE( KTITRE(50:64),'(G15.6)' ) QUAMIN
ccc      WRITE( KTITRE(74:88),'(G15.6)' ) VTEMIN
ccc      WRITE( KTITRE(106:120),'(G15.6)' ) VTEMIN/VOLMOY
ccc      CALL SANSDBL( KTITRE, NBC )
ccc      PRINT*,KTITRE(1:NBC)
ccc      CALL TRFETO9( KTITRE(1:NBC), PTXYZD, NTEMIN, NOTETR )

cccC     TRACE DES TETRAEDRES FRONTIERE ET DU TETRAEDRE D'ARETE MAX
ccc      KTITRE='FIN qualgrte: TETRAEDRE               ARETE MAX=         '
ccc      WRITE( KTITRE(25:34),'(I10)'   ) NTEAMX
ccc      WRITE( KTITRE(50:64),'(G15.6)' ) ARMAXM
ccc      CALL TRFETO9( KTITRE(1:NBC), PTXYZD, NTEAMX, NOTETR )
ccc      TRACTE = TRACTE0


      RETURN
      END
