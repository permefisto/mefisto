      SUBROUTINE NEF2MEF
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : LIRE LE NOM D'UN OBJET MAILLE PAR NEF
C ----- RELIRE SON MAILLAGE ISSU DE NEF
C       CONSTRUIRE LES TMS XYZSOMMET et NSEF DE SES PLSV
C
C PLUS PRECISEMENT:
C     LECTURE SUR LE FICHIER NOFILE DES XYZ ET POUR CHAQUE EF
C     NO DES SOMMETS: NS1,...,NSm  m=Nb Sommets d'un EF
C     si EF 3d (HEXAEDRE, PENTAEDRE ou TETRAEDRE)
C     NO DES FACES: NF1,...,NFn, n=Nb Faces d'un EF 3d
C     si EF 3d: NO de son VOLUME i.e. MATERIAU
C     si EF 2d: NO de sa SURFACE i.e. MATERIAU
C ATTENTION:
C     LES NUMEROS DES COURBES COMMUNES AUX FACES   EN SONT DEDUITS
C     LES NUMEROS DES POINTS  COMMUNS  AUX COURBES EN SONT DEDUITS
C
C RESULTATS:
C     CREATION DES VOLUMES(si EF 3d), SURFACES, LIGNES et POINTS MEFISTO
C     de L'OBJET NEF lu SUR LE FICHIER NomObjet.XyzNoEF
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE LJLL UPMC PARIS   JUIN 2004
C2345X7..............................................................012
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
c
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/a_objet__definition.inc"
      include"./incl/a_objet__topologie.inc"
      include"./incl/a_volume__definition.inc"
      include"./incl/a_surface__definition.inc"
      include"./incl/a_ligne__definition.inc"
      include"./incl/a_point__definition.inc"
      include"./incl/a___nsef.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/xyzext.inc"
C
      include"./incl/pp.inc"
      COMMON             MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
      CHARACTER*24      KNOMOB, KNOM, KNMFIC
      LOGICAL           LOPEN, LEXIST
      DOUBLE PRECISION  X,Y,Z
      INTEGER           NUSTFA(3,4)
      DATA              NUSTFA / 1,2,3, 2,3,4, 3,4,1, 4,1,2 /
C
      IERR = 0
C
C     LECTURE DU NOM_DE_L'OBJET  LE CARACTERE @ POUR FINIR
 5    CALL INVITE( 48 )
      NCVALS = 0
      CALL LIRCAR( NCVALS , KNOMOB )
      IF( NCVALS .EQ. -1 ) GOTO 9999
C
C     LE NOM DU FICHIER ISSU DE NEF est KNOMOB.XyzNoEF
      N = NUDCNB( KNOMOB )
      KNMFIC = KNOMOB(1:N) // '.XyzNoEF'
      N = NUDCNB( KNMFIC )
      WRITE(IMPRIM,*)
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,*) 'LECTURE DU FICHIER ', KNMFIC(1:N),
     %                   ' MAILLAGE d''un OBJET issu de NEF'
      ELSE
         WRITE(IMPRIM,*) 'READING of FILE ', KNMFIC(1:N),
     %                   ' OBJECT MESH from software NEF'
      ENDIF
      WRITE(IMPRIM,10005)
10005 FORMAT(1X,80('='))
C
C     SI LE FICHIER N'EXISTE PAS, ALORS ERREUR
      INQUIRE( FILE=KNMFIC, EXIST=LEXIST, OPENED=LOPEN )
      IF( .NOT. LEXIST ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(2) = 'FICHIER ' // KNOMOB
            L = NUDCNB( KERR(2) )
            KERR(1) = KERR(2)(1:L) // '.XyzNoEF' // ' INCONNU'
         ELSE
            KERR(2) = 'FILE ' // KNOMOB
            L = NUDCNB( KERR(2) )
            KERR(1) = KERR(2)(1:L) // '.XyzNoEF' // ' UNKNOWN'
         ENDIF
         CALL LEREUR
         GOTO 5
      ENDIF
C
C     LE FICHIER KNMFIC EXISTE
      IF( .NOT. LOPEN ) THEN
C        OUVERTURE DU FICHIER
         CALL TRUNIT( NOFILE )
         OPEN( FILE=KNMFIC, UNIT=NOFILE, STATUS='OLD' )
      ENDIF
C
C     RECHERCHE DU NOM DE L'OBJET DANS LE LEXIQUE DES OBJETS
      IF( NTOBJE .LE. 0 ) RETURN
      CALL LXLXOU( NTOBJE, KNOMOB, NTLXOB, MNLXOB )
C
C     S'IL N'EXISTE PAS IL EST CREE ET OUVERT
 20   IF( NTLXOB .LE. 0 ) THEN
         CALL LXLXDC( NTOBJE, KNOMOB, 24, 8 )
         CALL LXLXOU( NTOBJE, KNOMOB, NTLXOB, MNLXOB )
         ILEXIS = 0
      ELSE
C        S'IL EXISTE CHOIX DE DESTRUCTION OU NON
         NBLGRC(NRERR) = 2
         IF( LANGAG .EQ. 0 ) THEN
            KERR(3) = 'L''OBJET ' // KNOMOB
            L = NUDCNB(KERR(3))
            KERR(1) = KERR(3)(1:L)  // ' EXISTE DEJA'
            KERR(2) = 'VOULEZ VOUS LE DETRUIRE?'
         ELSE
            KERR(3) = 'The OBJECT ' // KNOMOB
            L = NUDCNB(KERR(3))
            KERR(1) = KERR(3)(1:L)  // ' ALREADY EXISTS'
            KERR(2) = 'DO YOU WANT DESTROY IT?'
         ENDIF
         CALL LERESU
         CALL LIMTCL( 'non_oui', N )
         IF( N .EQ. 0 ) GOTO 5
C        DESTRUCTION DE L'OBJET
         CALL LXLXDS( NTOBJE, KNOMOB )
         NTLXOB = 0
         GOTO 20
      ENDIF
C
C     TRACE DU NOM DE L'OBJET
      CALL RECTEF( NRHIST )
      NBLGRC(NRHIST) = 1
      IF( LANGAG .EQ. 0 ) THEN
         KHIST(1) = 'OBJET: ' // KNOMOB
      ELSE
         KHIST(1) = 'OBJECT: ' // KNOMOB
      ENDIF
      CALL LHISTO
      ILEXIS = 1
C
C     LECTURE SUR LE FICHIER NOFILE DES XYZ ET NO DES EF
C     PLUS EXACTEMENT: NS1,...,NSm, NF1, ... , NFn, NV
      READ( NOFILE, * ) KNOM
      READ( NOFILE, * ) NBSOMMET, NBEF, NBNOEF
C     NBNOEF=15 => HEXAEDRES  D'UN COMPSOLID
C     NBNOEF=12 => PENTAEDRES D'UN COMPSOLID
C     NBNOEF= 9 => TETRAEDRES D'UN COMPSOLID
C     NBNOEF= 5 => QUADRANGLES D'UN SHELL
C     NBNOEF= 4 => TRIANGLES   D'UN SHELL
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10020) KNOM, NBSOMMET, NBEF, NBNOEF
      ELSE
         WRITE(IMPRIM,20020) KNOM, NBSOMMET, NBEF, NBNOEF
      ENDIF
10020 FORMAT(' LECTURE de l''OBJET NEF: ', A /
     %' NOMBRE de SOMMETS  =',I12/
     %' NOMBRE d'' EF       =',I12/
     %' NOMBRE de NO par EF=',I12/)
20020 FORMAT(' READING of the NEF OBJECT: ', A /
     %' NUMBER of VERTICES    =',I12/
     %' NUMBER of FE          =',I12/
     %' NUMBER of NUMBER by FE=',I12/)
C
      IF( NBNOEF .EQ. 9 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10021)
         ELSE
            WRITE(IMPRIM,20021)
         ENDIF
10021    FORMAT(' Les EF sont des TETRAEDRES')
20021    FORMAT(' The FE are TETRAHEDRA')
C        NOMBRE DE SOMMETS DU TETRAEDRE
         NBST1EF = 4
C        NOMBRE DE NO DES FACES DU TETRAEDRE
         NBFC1EF = 4
      ELSEIF( NBNOEF .EQ. 4 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10022)
         ELSE
            WRITE(IMPRIM,20022)
          ENDIF
10022    FORMAT(' Les EF sont des TRIANGLES')
20022    FORMAT(' The FE are TRIANGLES')
C        NOMBRE DE SOMMETS DU TRIANGLE
         NBST1EF = 3
C        NOMBRE DE NO DES FACES DU TRIANGLE
         NBFC1EF = 0
      ELSE
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:10),'(I10)') NBNOEF
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,*) 'EF NON TRAITE'
            KERR(1) ='Nb de NUMEROS PAR EF='//KERR(MXLGER)(1:10)
     %             // ' => EF NON TRAITE'
         ELSE
            WRITE(IMPRIM,*) 'UNKNOWN FE'
            KERR(1) ='NUMBER of NUMBERS by FE='//KERR(MXLGER)(1:10)
     %             // ' => UNKNOWN FE'
         ENDIF
         CALL LEREUR
         GOTO 5
      ENDIF
C
C     LE TABLEAU XYZ DES SOMMETS DE L'OBJET NEF
      MNXYZS = 0
      CALL TNMCDC( 'REEL', 3*NBSOMMET, MNXYZS )
      MNS = MNXYZS
      DO 30 N=1, NBSOMMET
         READ( NOFILE, * ) X, Y, Z
         RMCN( MNS     ) = REAL( X )
         RMCN( MNS + 1 ) = REAL( Y )
         RMCN( MNS + 2 ) = REAL( Z )
         MNS = MNS + 3
 30   CONTINUE
C
C     LES XYZ EXTREMES POUR LES TRACES
      COOEXT(1,1) = RMCN(MNXYZS)
      COOEXT(1,2) = RMCN(MNXYZS)
      COOEXT(2,1) = RMCN(MNXYZS+1)
      COOEXT(2,2) = RMCN(MNXYZS+1)
      COOEXT(3,1) = RMCN(MNXYZS+2)
      COOEXT(3,2) = RMCN(MNXYZS+2)
      MNS = MNXYZS + 3
      DO 32 N=2, NBSOMMET
         COOEXT(1,1) = MIN( COOEXT(1,1), RMCN(MNS) )
         COOEXT(1,2) = MAX( COOEXT(1,2), RMCN(MNS) )
         COOEXT(2,1) = MIN( COOEXT(2,1), RMCN(MNS+1) )
         COOEXT(2,2) = MAX( COOEXT(2,2), RMCN(MNS+1) )
         COOEXT(3,1) = MIN( COOEXT(3,1), RMCN(MNS+2) )
         COOEXT(3,2) = MAX( COOEXT(3,2), RMCN(MNS+2) )
         MNS = MNS + 3
 32   CONTINUE
C
C     LE TABLEAU DES NBNOEF NO DES NBEF EF
      MNNOST = 0
      CALL TNMCDC( 'ENTIER', NBNOEF*NBEF, MNNOST )
C
C     RECHERCHE DU NOMBRE DE VOLUMES ET SURFACES DE L'OBJET NEF
C     si EF3d LES NO DES 4 SOMMETS, 4 FACES OC, VOLUME DE CHAQUE TETRAEDRE
      NBVOLUMES = 0
      NBFACES   = 0
C     LE NOMBRE D'ARETES DES TRIANGLES FRONTALIERS
      NBTRFA = 0
      MNNEF  = MNNOST - 1 - NBNOEF
      DO 40 N = 1, NBEF
C        LECTURE SUR NOFILE DES NBNOEF NUMEROS DE L'EF N
         MNNEF = MNNEF + NBNOEF
         READ( NOFILE, * ) (MCN(MNNEF+I),I=1,NBNOEF)
C        LE NO MAX DES VOLUMES DES EF
         NBVOLUMES = MAX( NBVOLUMES, MCN(MNNEF+NBNOEF) )
C        PARCOURS DES NO DES FACES DE L'EF
         DO 35 I=NBST1EF+1, NBST1EF+NBFC1EF
            NF = MCN(MNNEF+I)
            IF( NF .GT. 0 ) THEN
               NBFACES = MAX( NBFACES, NF )
               NBTRFA  = NBTRFA + 1
            ENDIF
 35      CONTINUE
 40   CONTINUE
C
C     FERMETURE DU FICHIER KNOMOB.XyzNoEF
      CLOSE( NOFILE )
C
C     CONSTRUCTION DES VSLP DE L'OBJET SELON L'EXISTENCE DE VOLUMES OU NON
      IF( NBNOEF .GE. 9 ) THEN
C        OBJET AVEC DES VOLUMES
         CALL NEFMEF3D( KNOMOB, NTLXOB, NBVOLUMES, NBFACES,
     %                  NBSOMMET, MNXYZS,
     %                  NBEF, NBNOEF, MNNOST, NBST1EF, NBFC1EF,
     %                  NBTRFA )
      ELSE
C        OBJET SANS VOLUME
         CALL NEFMEF2D( KNOMOB, NTLXOB, NBVOLUMES,
     %                  NBSOMMET, MNXYZS,
     %                  NBEF, NBNOEF, MNNOST, NBST1EF )
      ENDIF
C
C     DESTRUCTIONS DE TABLEAUX MC
      CALL TNMCDS( 'REEL',   3*NBSOMMET,  MNXYZS )
      CALL TNMCDS( 'ENTIER', NBNOEF*NBEF, MNNOST )
CCC
CCCC     CONSTRUCTION DE LA TOPOLOGIE DE L'OBJET
CCC      CALL DFTOP0( KNOMOB, IERR )
C
C     SORTIE
 9999 RETURN
      END
