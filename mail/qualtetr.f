      SUBROUTINE QUALTETR( PTXYZD, MXTETR,  NOTETR,
     %                     NBTETR, NUDTETR, QUAMIN, QUAMOY, VOLUMT )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  CALCUL DES QUALITES MIN MOYENNE et VOLUME de la TETRAEDRISATION
C -----
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

C SORTIES:
C --------
C NBTETR : NOMBRE DE TETRAEDRES ACTIFS DE NOTETR
C NUDTETR: PLUS GRAND NUMERO DANS NOTETR DE TETRAEDRE OCCUPE
C QUAMIN : QUALITE MINIMALE DES TETRAEDRES DE NOTETR
C QUAMOY : QUALITE MOYENNE  DES TETRAEDRES DE NOTETR
C VOLUMT : VOLUME  TOTAL    DES TETRAEDRES DE NOTETR

C REMARQUE: EN SORTIE LE VOLUME MOYEN D'UN TETRAEDRE VAUT VOLUMT/NBTETR
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC ET ST PIERRE DU PERRAY     Mai 2008
C2345X7..............................................................012
      PARAMETER        (QUAMED=0.005)
      include"./incl/langue.inc"
      COMMON / TRTETR / STOPTE, TRACTE
      LOGICAL           STOPTE, TRACTE, TRACTE0
C     TRACTE = FAUX ==> PAS DE TRACE DES TETRAEDRES
C              VRAI ==> TRACE DES TETRAEDRES DE L'ETOILE

      DOUBLE PRECISION  PTXYZD(1:4,*), VOLUMT, QL
      INTEGER           NOTETR(8,MXTETR)
      REAL              QUAMIN, QUAMOY, QUALTE
      DOUBLE PRECISION  ARMIN,  ARMAX,  SURFTR(4), VOLUTE,
     %                  ARMINM, ARMAXM, VTEMIN, VOLMOY
      CHARACTER*124     KTITRE

      TRACTE0= TRACTE
      VOLUMT = 0D0
      NBTETR = 0
      NBTETM = 0
      QUAMOY = 0
      QUAMIN = 2
      NTEMIN = 0
      ARMINM = 1D100
      ARMAXM = 0D0
      NUDTETR= 0
      NUDSOMM= 0
C     NOMBRE DE TETRAEDRE AYANT AU MOINS UNE FACE DE TETRAEDRE OPPOSE <0
C     C-A-D DE TETRAEDRE OPPOSE INCONNU
      NBTEOPI= 0

      DO NTE = 1, MXTETR

         IF( NOTETR(1,NTE) .GT. 0 ) THEN

            NBTETR  = NBTETR + 1
            NUDTETR = NTE
            NUDSOMM = MAX( NUDSOMM, NOTETR(1,NTE), NOTETR(2,NTE),
     %                              NOTETR(3,NTE), NOTETR(4,NTE) )

            CALL QUATETD( PTXYZD(1,NOTETR(1,NTE)),
     %                    PTXYZD(1,NOTETR(2,NTE)),
     %                    PTXYZD(1,NOTETR(3,NTE)),
     %                    PTXYZD(1,NOTETR(4,NTE)),
     %                    ARMIN, ARMAX, SURFTR, VOLUTE, QUALTE )

            IF( VOLUTE .LE. 0D0 .OR. QUALTE .LE. 0.01 ) THEN

C              TETRAEDRE DE QUALITE MEDIOCRE
               NBTETM  = NBTETM + 1

ccc               PRINT*,'qualtetr: TETRAEDRE(',NTE,')=',
ccc     %                 (NOTETR(kk,NTE),kk=1,8),
ccc     %                ' V=',VOLUTE,' Q=',QUALTE

ccc               DO N=1,4
ccc                  NST = NOTETR( N, NTE )
ccc                  PRINT*,'PTXYZD(',NST,')=',(PTXYZD(kk,NST),kk=1,4)
ccc               ENDDO
ccc               TRACTE  = .TRUE.
ccc               KTITRE='TETRAEDRE               QUALITE=                 
ccc     %VOLUME=             '
ccc               WRITE( KTITRE(11:20),'(I10)'   ) NTE
ccc               WRITE( KTITRE(33:47),'(G15.6)' ) QUALTE
ccc               WRITE( KTITRE(57:71),'(G15.6)' ) VOLUTE
ccc               CALL SANSDBL( KTITRE, NBC )
ccc               CALL TRFETO9( KTITRE(1:NBC), PTXYZD, NTE, NOTETR )

            ENDIF

            QUAMOY = QUAMOY + QUALTE
            IF( QUALTE .LT. QUAMIN ) THEN
               QUAMIN = QUALTE
               NTEMIN = NTE
               VTEMIN = VOLUTE
            ENDIF

            IF( ARMIN .LT. ARMINM ) THEN
               ARMINM = ARMIN
               NTEAMI = NTE
            ENDIF

            IF( ARMAX .GT. ARMAXM ) THEN
               ARMAXM = ARMAX
               NTEAMX = NTE
            ENDIF

C           VOLUME TOTAL DES TETRAEDRES
            VOLUMT = VOLUMT + VOLUTE

C           TETRAEDRE OPPOSE A UNE FACE <0?
            DO N=5,8
               IF( NOTETR(N,NTE) .LT. 0 ) THEN
                  NBTEOPI = NBTEOPI + 1
                  PRINT*,'qualtetr: TETRAEDRE',NTE,
     %            ' St:',(NOTETR(KK,NTE),kk=1,4),
     %            ' Tetra Oppose:',(NOTETR(KK,NTE),kk=5,8)
               ENDIF
            ENDDO
         ENDIF

      ENDDO

C     QUALITE MOYENNE DE LA TETRAEDRISATION
      QUAMOY = QUAMOY / NBTETR

C     QUALITE MIN EN DOUBLE PRECISION POUR AFFICHAGE
      QL = QUAMIN

C     VOLUME MOYEN D'UN TETRAEDRE
      VOLMOY = VOLUMT / NBTETR

      IF( LANGAG .EQ. 0 ) THEN

         PRINT*,'qualtetr:',NBTETR,
     %          ' TETRAEDRES de VOLUME',VOLUMT,
     %          ' No DERNIER SOMMET=',NUDSOMM,
     %          ' No DERNIER TETRAEDRE=',NUDTETR,
     %          ' QUALITE MOYENNE=',QUAMOY,
     %          ' dont',NBTETM,' TETRAEDRES de QUALITE<0.01'
         PRINT*,'QUALITE MINIMALE=',QL,
     %          ' du TETRAEDRE',NTEMIN,':',(NOTETR(K,NTEMIN),K=1,8),
     %          ' de VOLUME=',VTEMIN,' VOLUME MOYEN=',VOLMOY
         PRINT*,'ARETE MINIMALE  =',ARMINM,
     %          ' du TETRAEDRE',NTEAMI,':',(NOTETR(K,NTEAMI),K=1,8)
         PRINT*,'ARETE MAXIMALE  =',ARMAXM,
     %          ' du TETRAEDRE',NTEAMX,':',(NOTETR(K,NTEAMX),K=1,8)
         IF( NBTEOPI .GT. 0 ) THEN
            PRINT*,'Nb de TETRAEDRES OPPOSES INCONNUS=', NBTEOPI
         ENDIF

      ELSE

         PRINT*,'qualtetr:',NBTETR,
     %          ' TETRAHEDRA of VOLUME',VOLUMT,
     %          ' LAST VERTEX NUMBER=',NUDSOMM,
     %          ' LAST TETRAHEDRA NUMBER=',NUDTETR,
     %          ' MEAN QUALITY=',QUAMOY,
     %          ' QUALITY<0.01 TETRAHEDRA NUMBER=',NBTETM
         PRINT*,'MINIMUM QUALITY=',QL,
     %          ' of TETRAHEDRON',NTEMIN,':',(NOTETR(K,NTEMIN),K=1,8),
     %          ' of VOLUME=',VTEMIN,' MEAN VOLUME=',VOLMOY
         PRINT*,'MINIMUM EDGE   =',ARMINM,
     %          ' of TETRAHEDRON',NTEAMI,':',(NOTETR(K,NTEAMI),K=1,8)
         PRINT*,'MAXIMUM EDGE   =',ARMAXM,
     %          ' of TETRAHEDRON',NTEAMX,':',(NOTETR(K,NTEAMX),K=1,8)
         IF( NBTEOPI .GT. 0 ) THEN
            PRINT*,' UNKNOWN OPPOSED TETRAHEDRA Number=', NBTEOPI
         ENDIF
      ENDIF

C     TRACE DU TETRAEDRE DE PLUS MAUVAISE QUALITE ET SES VOISINS
      IF( QUAMIN .LE. QUAMED ) THEN
ccc         TRACTE = .TRUE.
         KTITRE='qualtetr: TETRAEDRE                de QUALITE MIN=     
     %           VOLUME=                 VTEMIN/VOLMOY=                '
         WRITE( KTITRE(25:34),  '(I10)'   ) NTEMIN
         WRITE( KTITRE(51:65),  '(G15.6)' ) QUAMIN
         WRITE( KTITRE(75:89),  '(G15.6)' ) VTEMIN
         WRITE( KTITRE(107:121),'(G15.6)' ) VTEMIN/VOLMOY
         CALL SANSDBL( KTITRE, NBC )
         PRINT*,KTITRE(1:NBC)
         CALL TRFETO9( KTITRE(1:NBC), PTXYZD, NTEMIN, NOTETR )
      ENDIF

cccC     TRACE DES TETRAEDRES FRONTIERE ET DU TETRAEDRE D'ARETE MAX
ccc      IF( ARMAXM .GT. 2.8D0 ) then
ccc      KTITRE='qualtetr: TETRAEDRE               ARETE MAX=         '
ccc      WRITE( KTITRE(25:34),'(I10)'   ) NTEAMX
ccc      WRITE( KTITRE(50:64),'(G15.6)' ) ARMAXM
ccc      CALL TRFETO9( KTITRE(1:NBC), PTXYZD, NTEAMX, NOTETR )
ccc      CALL TRTEFR( KTITRE, 7, NUDSOMM, PTXYZD, NUDTETR, NOTETR, NTEAMX )
ccc      endif

cccC     TRACE DES 3 ARETES DES 4 FACES DE TOUS LES TETRAEDRES DU MAILLAGE
ccc      CALL TRALLTE( KTITRE, PTXYZD, NUDTETR, NOTETR )

      TRACTE = TRACTE0
      RETURN
      END
