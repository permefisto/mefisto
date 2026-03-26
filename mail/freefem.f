      SUBROUTINE FREEFEM
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : SORTIR SUR UN FICHIER freefem.NomObjet selon la dimension de l'espace
C ----- si NDIM=2 LES VERTICES, TRIANGLES P1 ET EDGES de la FRONTIERE
C       si NDIM=3 LES VERTICES, TETRAEDRES P1 ET FACES TRIANGULAIRES
C                 de la FRONTIERE
C       D'UN OBJET MAILLE PAR MEFISTO C-A-D A PARTIR DE SES
C       NUMEROS DE NOEUDS DES EF, DES VOLUMES SURFACES LIGNES POINTS
C       ET CELA POUR UNE EXECUTION ULTERIEURE DE FREEFEM
C
C ATTENTION: SOMMETS et NOEUDS et POINTS SONT ICI SUPPOSES IDENTIQUES!
C       i.e. INTERPOLATION LAGRANGE DE DEGRE 1 SUR TRIANGLE ou TETRAEDRE
C
C SORTIE :
C --------
C LE FICHIER freefem.NomObjet DANS LE REPERTOIRE DU PROJET
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET ANALYSE NUMERIQUE LJLL UPMC PARIS OCTOBRE 2003
C23456---------------------------------------------------------------012
      IMPLICIT          INTEGER (W)
      PARAMETER         (MXTYEL=7)
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/a_objet__topologie.inc"
      include"./incl/a_objet__definition.inc"
      include"./incl/a_ligne__definition.inc"
      include"./incl/a_transfo__definition.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___xyznoeud.inc"
      include"./incl/a___xyzpoint.inc"
      include"./incl/a___npef.inc"
      include"./incl/a___texte.inc"
      include"./incl/msvaau.inc"
      include"./incl/ponoel.inc"
      include"./incl/pp.inc"
      COMMON             MCN(MOTMCN)
      COMMON / UNITES /  LECTEU,IMPRIM,INTERA,NUNITE(29)
      REAL               RMCN(1)
      EQUIVALENCE       (MCN(1),RMCN(1))
      CHARACTER*4        NOMELE(2)
      CHARACTER*10       NMTYOB,KNM
      CHARACTER*24       KNOMOB,KNOM
      CHARACTER*32       KNMFIC
      LOGICAL            LEXIST,LOPEN
      INTEGER            NBPLSV(4), NOSOTR(3)
C
C     NOM DE L'OBJET A TRAITER
C     ------------------------
 10   CALL INVITE( 45 )
      IERR   = 0
      NCVALS = 0
      CALL LIRCAR( NCVALS, KNOMOB )
      IF( NCVALS .EQ. -1 ) RETURN
C
      MNAUX  = 0
      MNLABT = 0
      MNLABE = 0
      MNLABV = 0
      MNAUX  = 0
      NBSTEF = 0
      MXTRIA = 0
      MXEDGE = 0
C
C     RECHERCHE DU NOM DE L'OBJET DANS LE LEXIQUE DES OBJETS
      CALL LXLXOU( NTOBJE, KNOMOB, NTLXOB, MNLXOB )
C     S'IL N'EXISTE PAS RETOUR A LA DEMANDE DU NOM DE L'OBJET
      IF( NTLXOB .LE. 0 ) THEN
         NBLGRC(NRERR) = 2
         KERR(1) = KNOMOB
         IF( LANGAG .EQ. 0 ) THEN
            KERR(2) = 'OBJET INCONNU'
         ELSE
            KERR(2) = 'UNKNOWN OBJECT'
         ENDIF
         CALL LEREUR
         CALL LXIM( NTOBJE )
         GOTO 10
      ENDIF
C
C     TRACE DE L'OBJET
      NBLGRC(NRHIST) = 1
      KHIST(1) = 'OBJET: ' // KNOMOB
      CALL LHISTO
      IF( INTERA .GE. 1 ) THEN
C        TRACE DE L'OBJET SI CONSOLE GRAPHIQUE INTERACTIVE
         CALL T1OBJE( KNOMOB )
      ENDIF
C
C     RECHERCHE DU TABLEAU DEFINITION DE L'OBJET
      CALL LXTSOU( NTLXOB, 'DEFINITION', NTDFOB, MNDFOB )
C     S'IL N'EXISTE PAS RETOUR A LA DEMANDE DU NOM DE L'OBJET
      IF( NTDFOB .LE. 0 ) THEN
         NBLGRC(NRERR) = 2
         KERR(1) = KNOMOB
         IF( LANGAG .EQ. 0 ) THEN
            KERR(2) = 'OBJET SANS TMS de DEFINITION'
         ELSE
            KERR(2) = 'OBJECT WITHOUT the TMS DEFINITION'
         ENDIF
         CALL LEREUR
         GOTO 10
      ENDIF
C
C     ADRESSAGE DES ADRESSES DES TABLEAUX NPEF"xxxx DE CET OBJET
      MNELEM = 0
      CALL TNMCDC( 'ENTIER', 2*MXTYEL, MNELEM )
      MNTELE = MNELEM + MXTYEL
C
C     RECHERCHE DES TABLEAUX XYZSOMMET XYZNOEUD XYZPOINT ASSOCIES A L'OBJET
      CALL NDPGEL( NTLXOB, NTTOPO, MNTOPO ,
     %             NTXYZP, MNXYZP, NTXYZN, MNXYZN ,
     %             NBTYEL, MCN(MNTELE), MCN(MNELEM), IERR )
      IF( IERR .NE. 0 ) GOTO 9900
C
C     LE TYPE DE MAILLAGE ET LES OBJETS INTERNES ET AUX LIMITES DE L'OBJET
      NDPGST = MCN( MNTOPO + WDPGST )
      NBOBIN = MCN( MNTOPO + WBOBIN )
      NBOBCL = MCN( MNTOPO + WBOBCL )
C
C     LA VALEUR DE Z DES POINTS DU MAILLAGE DE L'OBJET DEFINIT
C     NDIM LA DIMENSION 2 OU 3 DE L'ESPACE DES COORDONNEES
C     ATTENTION: CETTE PROGRAMMATION SOUS ENTEND NOEUDS=POINTS
      CALL DIMCOO( MCN(MNXYZP+WNBPOI), MCN(MNXYZP+WYZPOI), NDIM )
C     NOMBRE DE SOMMETS 3 POUR UN TRIANGLE, 4 POUR UN TETRAEDRE
      NBSTEF = NDIM + 1
C
C     SANS ou AVEC les TANGENTES de Mefisto dans le FICHIER freefem? => SANS
CCC      CALL LIMTCL( 'tgoupas', NOTG )
CCC      IF( NOTG .LT. 0 ) RETURN
C     PAS DE TANGENTES STOCKEES DANS FREEFEM
      NOTG = 0
C
C     LE NOMBRE TOTAL DE NOEUDS=POINTS DU MAILLAGE DE L'OBJET
      NBNOEU = MCN( MNXYZN + WNBNOE )
C     LE NOMBRE TOTAL DE TANGENTES DU MAILLAGE
      IF( NOTG .GT. 0 ) THEN
         NBTGN = MCN( MNXYZN + WNBTGN )
      ELSE
         NBTGN = 0
      ENDIF
C
      WRITE(IMPRIM,*)
      IF( LANGAG .EQ. 0 ) THEN
      WRITE(IMPRIM,10059) KNOMOB
10059 FORMAT(A24, '  NOM DE L''OBJET pour FreeFem')
      WRITE(IMPRIM,10210) NDIM,NBNOEU,NBTGN,NBOBIN,NBOBCL,NBTYEL
10210 FORMAT(' DIMENSION 2 OU 3 DE L''ESPACE',T32,'=',I9/
     %' NOMBRE DE NOEUDS'                    ,T32,'=',I9/
     %' NOMBRE DE TANGENTES'                 ,T32,'=',I9/
     %' NOMBRE D''OBJETS INTERNES'           ,T32,'=',I9/
     %' NOMBRE D''OBJETS AUX LIMITES'        ,T32,'=',I9/
     %' NOMBRE DE TYPES D''EF'               ,T32,'=',I9)
      ELSE
      WRITE(IMPRIM,11059) KNOMOB
11059 FORMAT(A24, '  NAME of the OBJECT for Freefem')
      WRITE(IMPRIM,11210) NDIM,NBNOEU,NBTGN,NBOBIN,NBOBCL,NBTYEL
11210 FORMAT(' SPACE''s DIMENSION (2 or 3)',T33,'=',I9/
     %' NUMBER of NODES'                   ,T33,'=',I9/
     %' NUMBER of TANGENTS'                ,T33,'=',I9/
     %' NUMBER of MATERIALS'               ,T33,'=',I9/
     %' NUMBER of PL(S) at the BOUNDARY'   ,T33,'=',I9/
     %' NUMBER of TYPES of FE'             ,T33,'=',I9)
      ENDIF
C
C     BOUCLE SUR LES OBJETS "INTERNES": calcul du nombre de TRIANGLES Freefem
C     =======================================================================
      WRITE(IMPRIM,*)
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10100) NBOBIN
      ELSE
         WRITE(IMPRIM,11100) NBOBIN
      ENDIF
10100 FORMAT(/,I10, '  {NBOBIN NOMBRE DE MATERIAUX}' )
11100 FORMAT(/,I10, '  {NBOBIN NUMBER of MATERIALS}' )
C
C     L'ADRESSE MCN DU DEBUT DU TABLEAU NUOBIN
      MNOBIN = MNTOPO + WMTYEL + NBTYEL
      MN     = MNOBIN - 2
      DO 100 I=1,NBOBIN
C        LE TYPE DE L'OBJET
         MN   = MN + 2
         NYOB = MCN( MN )
         NUOB = MCN( MN + 1 )
C        OUVERTURE DE L'OBJET
         CALL LXNLOU( NTMN(NYOB), NUOB, NTOB, MNOB )
         KNM = NMTYOB( NYOB )
         CALL NMOBNU( KNM, NUOB, KNOM )
         IF( NTOB .LE. 0 ) THEN
            NBLGRC(NRERR) = 2
            KERR(1) = KNOM
            IF( LANGAG .EQ. 0 ) THEN
               KERR(2) = 'MATERIAU INCONNU'
            ELSE
               KERR(2) = 'UNKNOWN MATERIAL'
            ENDIF
            CALL LEREUR
            IERR = IERR + 1
            GOTO 100
         ENDIF
C        ECRITURE DU NOM DE L'OBJET INTERNE ou MATERIAU I
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10102) NUOB, KNM, KNOM
         ELSE
            WRITE(IMPRIM,11102) NUOB, KNM, KNOM
         ENDIF
 100  CONTINUE
10102 FORMAT(I10,'   {est le NUMERO ou LABEL dans ',A10,' de } ', A24 )
11102 FORMAT(I10,'   {is the NUMBER or LABEL in ',A10,' of } ', A24 )
      IF( IERR .NE. 0 ) GOTO 9900
      NBPLSV(NDIM+1) = NBOBIN
C
C     BOUCLE SUR LES OBJETS AUX LIMITES DE L'OBJET
C     RECHERCHE DU NOMBRE DE POINTS LIGNES SURFACES DE L'OBJET
C     ========================================================
C     RECENSEMENT DU NOMBRE DE POINTS, LIGNES, ... DE L'OBJET
C     BOUCLE SUR LES SURFACES(SI NDIM=3) LIGNES ET POINTS
C     L'ADRESSE MCN DU DEBUT DU TABLEAU NUOBCL
      MNOBCL = MNOBIN + MOTVAR(13) * NBOBIN
      DO 130 L=NDIM,1,-1
         NBPLSV(L) = 0
         KNM = NMTYOB( L )
         MN  = MNOBCL - 2
         DO 120 I=1,NBOBCL
C           LE TYPE DE L'OBJET
            MN   = MN + 2
            NYOB = MCN( MN )
            IF( NYOB .NE. L ) GOTO 120
            NUOB = MCN( MN + 1 )
C           OUVERTURE DE L'OBJET
            CALL LXNLOU( NTMN(NYOB), NUOB, NTOB, MNOB )
            CALL NMOBNU( KNM, NUOB, KNOM )
            IF( NTOB .LE. 0 ) THEN
               NBLGRC(NRERR) = 1
               IF( LANGAG .EQ. 0 ) THEN
                  KERR(1) = KNM // ':' // KNOM // ' INCONNU'
               ELSE
                  KERR(1) = KNM // ': UNKNOWN ' // KNOM
               ENDIF
               CALL LEREUR
               IERR = IERR + 1
               GOTO 120
            ENDIF
C           UN PLS DE TYPE L DE PLUS
            NBPLSV(L) = NBPLSV(L) + 1
 120     CONTINUE
 130  CONTINUE
C
C     BOUCLE SUR LES SURFACES(SI NDIM=3) LIGNES ET POINTS AUX LIMITES
      WRITE(IMPRIM,*)
      DO 150 L=NDIM,1,-1
C        ECRITURE DU NOMBRE D'OBJETS AUX LIMITES DE TYPE L
         KNM = NMTYOB( L )
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10140) NBPLSV(L),KNM
         ELSE
            WRITE(IMPRIM,11140) NBPLSV(L),KNM
         ENDIF
10140    FORMAT(/,I10, '  {NOMBRE DE ',A10,'}' )
11140    FORMAT(/,I10, '  {NUMBER of ',A10,'}' )
         MN  = MNOBCL - 2
         DO 140 I=1,NBOBCL
C           LE TYPE DE L'OBJET
            MN   = MN + 2
            NYOB = MCN( MN )
            IF( NYOB .NE. L ) GOTO 140
            NUOB = MCN( MN + 1 )
C           OUVERTURE DE L'OBJET
            CALL LXNLOU( NTMN(NYOB), NUOB, NTOB, MNOB )
            CALL NMOBNU( KNM, NUOB, KNOM )
            IF( NTOB .LE. 0 ) THEN
               NBLGRC(NRERR) = 1
               IF( LANGAG .EQ. 0 ) THEN
                  KERR(1) = KNM // ':' // KNOM // ' INCONNU'
               ELSE
                  KERR(1) = KNM // ': UNKNOWN ' // KNOM
               ENDIF
               CALL LEREUR
               GOTO 140
            ENDIF
C           ECRITURE DU NOM DE L'OBJET AUX LIMITES
            IF( LANGAG .EQ. 0 ) THEN
               WRITE(IMPRIM,10102) NUOB, KNM, KNOM
            ELSE
               WRITE(IMPRIM,11102) NUOB, KNM, KNOM
            ENDIF
 140     CONTINUE
 150  CONTINUE
      IF( IERR .GT. 0 ) GOTO 9900
C
C     LA BOUCLE SUR LES DIFFERENTS TYPES D'EF (TABLEAUX NPEF"xxxx)
C     CALCUL DU NOMBRE DE TRIANGLES POUR FREEFEM
C     ============================================================
C     ECRITURE DU NOMBRE DE TYPE D'EF
      WRITE(IMPRIM,*)
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10390) NBTYEL
      ELSE
         WRITE(IMPRIM,11390) NBTYEL
      ENDIF
10390 FORMAT(/,I10,'  {NBTYEL NOMBRE DE TYPES D''EF de l''OBJET}' )
11390 FORMAT(/,I10,'  {NBTYEL NUMBER of TYPES of FE of the OBJECT}' )
C
      NBEF = 0
      DO 403 NOTYEL = 1, NBTYEL
C
C        L'ADRESSE DU TABLEAU NPEF"xxxx
         MNELE = MCN( MNELEM - 1 + NOTYEL )
C
C        LE NUMERO DU TYPE DE L'ELEMENT FINI
         NUTYEL = MCN( MNELE + WUTYEL )
C
C        NCOGEL LE CODE GEOMETRIQUE DE L'ELEMENT FINI EST CALCULE
         CALL ELNUCG( NUTYEL, NCOGEL )
C
C        LE NOMBRE D'ELEMENTS FINIS DE CE TYPE
         NBELEM = MCN( MNELE + WBELEM )
C
         IF( NUTYEL .EQ.  1 .OR.
     %       NUTYEL .EQ. 13 .OR.
     %       NUTYEL .EQ. 29 ) THEN
C
C           TRIANGLE AVEC 3 SOMMETS
            NBEF = NBEF + NBELEM
C
         ELSE IF( NUTYEL .EQ.  3 .OR.
     %            NUTYEL .EQ. 16 ) THEN
C
C           QUADRANGLE AVEC 4 SOMMETS => 2 TRIANGLES
            NBEF = NBEF + NBELEM * 2
C
C
         ELSE IF( NUTYEL .EQ. 19 ) THEN
C
C           TETRAEDRE AVEC 4 SOMMETS
            NBEF = NBEF + NBELEM
C
         ENDIF
C
 403  CONTINUE
C
C     VERIFICATIONS
      IF( NDIM .EQ. 2 .AND. NBEF .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'AUCUN TRIANGLE ou QUADRANGLE POUR CET OBJET'
         ELSE
            KERR(1) = 'NO TRIANGLE or QUADRANGLE FOR THIS OBJECT'
         ENDIF
         GOTO 9900
      ELSE IF( NDIM .EQ. 3 .AND. NBEF .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'AUCUN ELEMENT FINI 3D POUR CET OBJET'
         ELSE
            KERR(1) = 'NO 3D FINITE ELEMENT FOR THIS OBJECT'
         ENDIF
         GOTO 9900
      ENDIF
C
      IF( LANGAG .EQ. 0 ) THEN
         IF( NDIM .EQ. 2 ) THEN
          WRITE(IMPRIM,*)NBEF,' TRIANGLES pour le maillage vers FreeFem'
         ELSE
         WRITE(IMPRIM,*)NBEF,' TETRAEDRES pour le maillage vers FreeFem'
         ENDIF
      ELSE
         IF( NDIM .EQ. 2 ) THEN
            WRITE(IMPRIM,*) NBEF, ' TRIANGLES for the FreeFem file'
         ELSE
            WRITE(IMPRIM,*) NBEF, ' TETRAHEDRA for the FreeFem file'
         ENDIF
      ENDIF
C
C     TABLEAU DES TRIANGLES ou TETRAEDRES DE FREEFEM NS1 NS2 NS3 (NS4) LABEL
      IF( NDIM .EQ. 2 ) THEN
C        1 QUADRANGLE => 2 TRIANGLES
         MXTRIA = NBEF * 2
      ELSE
         MXTRIA = NBEF
      ENDIF
      CALL TNMCDC( 'ENTIER', (NBSTEF+1)*MXTRIA, MNLABT )
C
C     TABLEAU DES EDGES ON BOUNDARY DE FREEFEM NS1 NS2 LABEL
      IF( NBEF .LE. 100 ) THEN
         MXEDGE = NBSTEF * NBEF
      ELSE
         MXEDGE = NBEF
      ENDIF
      CALL TNMCDC( 'ENTIER', NBSTEF*MXEDGE, MNLABE )
C
C     TABLEAU DU LABEL DES VERTICES DE FREEFEM
      CALL TNMCDC( 'ENTIER', NBNOEU, MNLABV )
      CALL AZEROI( NBNOEU, MCN(MNLABV) )
C
C     RESERVATION D'UN TABLEAU AUXILIAIRE POUR LES NOEUDS,..., TGS D'UN EF
      CALL TNMCDC( 'ENTIER', 27+8+12+6+1+24, MNAUX )
      MNAUX1 = MNAUX - 1
      MNPS   = MNAUX + 27
      MNPS1  = MNPS  - 1
      MNLA   = MNPS  + 8
      MNLA1  = MNLA  - 1
      MNSF   = MNLA  + 12
      MNSF1  = MNSF  - 1
      MNVC   = MNSF  + 6
      MNVC1  = MNVC  - 1
      MNTGEL = MNVC  + 1
      MNTGE1 = MNTGEL -1
C
CCCC     NO DES 3 SOMMETS DES 4 FACES DU TETRAEDRE deja fait dans call eltyca
CCC      IF( NDIM .EQ. 3 ) CALL SOFACU( 5, NBSOFA, NOSOFA )
CCCC     NOSOFA : NOSOFA(I,J) NO DU I-EME SOMMET DE LA FACE J DU "CUBE"
CCCC     F1:132  F2:143  F3:124  F4:234
C
C     LA BOUCLE SUR LES DIFFERENTS TYPES D'EF (TABLEAUX NPEF"xxxx)
C     CREATION DES TRIANGLES ET EDGES POUR FREEFEM DANS MNLAB...
C     ============================================================
      NBEF   = 0
      NBEDGE = 0
      MNLV   = MNLABV - 1
      MNLT   = MNLABT - 1
      MNLE   = MNLABE - 1
C
      DO 500 NOTYEL = 1, NBTYEL
C
C        L'ADRESSE DU TABLEAU NPEF"xxxx
         MNELE = MCN( MNELEM - 1 + NOTYEL )
C
C        LE NUMERO DU TYPE DE L'ELEMENT FINI
         NUTYEL = MCN( MNELE + WUTYEL )
C
C        NCOGEL LE CODE GEOMETRIQUE DE L'ELEMENT FINI EST CALCULE
         CALL ELNUCG( NUTYEL, NCOGEL )
C
C        LE NOMBRE D'ELEMENTS FINIS DE CE TYPE
         NBELEM = MCN( MNELE + WBELEM )
C
C        L'ADRESSE MCN DES NOEUDS ET POINTS GEOMETRIQUES DES ELEMENTS FINIS
         MNNDEL = MNELE + WUNDEL
         MNPGEL = MNNDEL
         IF( NDPGST .GE. 2 ) THEN
            MNPGEL = MNPGEL + NBELEM * MCN(MNELE+WBNDEL)
         ENDIF
C
C        LES CARACTERISTIQUES DE L'ELEMENT FINI
         CALL ELNUNM( NUTYEL, NOMELE )
         CALL ELTYCA( NUTYEL )
C
C        LE NOMBRE DE VOLUMES DE CE TYPE D'EF
         IF( NDIM .EQ. 3 ) THEN
            NBVOL = 1
            NBFAC = NFACE
         ELSE
C           EN DIMENSION 2 PAS DE VOLUME ET 1 FACE
            NBVOL = 0
            NBFAC = 1
         ENDIF
C        LE NOMBRE DE SOMMETS DE CE TYPE D'EF
         NBSOE = NBSOME(NCOGEL)
C
         IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10400) NOTYEL
10400    FORMAT(/,I10, '  {NUMERO DU TYPE D''EF}' )
         WRITE(IMPRIM,10401) NBELEM
10401    FORMAT(I10, '  {NBELEM NOMBRE D''EF DE CE TYPE D''EF}' )
         WRITE(IMPRIM,10402) NBNOE
10402    FORMAT(I10, '  {NBNOE  NOMBRE DE NOEUDS DE CE TYPE D''EF}')
         WRITE(IMPRIM,10403) NBSOE
10403    FORMAT(I10, '  {NBSOE  NOMBRE DE SOMMETS DE CE TYPE D''EF}')
         WRITE(IMPRIM,10404) NARET
10404    FORMAT(I10, '  {NBARET NOMBRE D''ARETES DE CE TYPE D''EF}')
         WRITE(IMPRIM,10405) NBFAC
10405    FORMAT(I10, '  {NBFACE NOMBRE DE FACES DE CE TYPE D''EF}')
         WRITE(IMPRIM,10406) NBVOL
10406    FORMAT(I10, '  {NBVOL  NOMBRE DE VOLUME DE CE TYPE D''EF}')
         ELSE
         WRITE(IMPRIM,11400) NOTYEL
11400    FORMAT(/,I10, '  {NUMBER OF THIS TYPE of FE}' )
         WRITE(IMPRIM,11401) NBELEM
11401    FORMAT(I10, '   {NBELEM NUMBER of FE of THIS TYPE}' )
         WRITE(IMPRIM,11402) NBNOE
11402    FORMAT(I10, '   {NBNOE  NUMBER of NODES of a FE of THIS TYPE}')
         WRITE(IMPRIM,11403) NBSOE
11403    FORMAT(I10,
     %'   {NBSOE  NUMBER of VERTICES of a FE of THIS TYPE}')
         WRITE(IMPRIM,11404) NARET
11404    FORMAT(I10, '   {NBARET NUMBER of EDGES of a FE of THIS TYPE}')
         WRITE(IMPRIM,11405) NBFAC
11405    FORMAT(I10, '   {NBFACE NUMBER of FACES of a FE of THIS TYPE}')
         WRITE(IMPRIM,11406) NBVOL
11406    FORMAT(I10,
     %   '   {NBVOL  NUMBER of VOLUMES of a FE of THIS TYPE}')
         ENDIF
C
         CALL AZEROI(  8, MCN(MNPS) )
         CALL AZEROI( 12, MCN(MNLA) )
         CALL AZEROI(  6, MCN(MNSF) )
         MCN( MNVC ) = 0
C
C        LES ADRESSAGES POUR RECUPERER LES INFORMATIONS
C        SELON LE TYPE DE L'ELEMENT FINI lIMITE ICI A TRIA et QUAD DEGRE 1
C        OU TETRAEDRE P1
         GOTO( 401, 400, 401, 400, 400, 400, 400, 400, 400, 400,
     &         400, 400, 401, 400, 400, 401, 400, 400, 401, 400,
     &         400, 400, 400, 400, 400, 400, 400, 400, 401, 400),NUTYEL
C
C        ERREUR
 400     NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'ELEMENT FINI '// NOMELE(1)
     &             // NOMELE(2) //' NON PROGRAMME pour FREEFEM'
         ELSE
            KERR(1) = 'FINITE ELEMENT '// NOMELE(1)
     &             // NOMELE(2) //' NOT PROGRAMMED for FREEFEM'
         ENDIF
         CALL LEREUR
         GOTO 9900
C
C        --------------------------------------
C        LA BOUCLE SUR LES EF DE CE TYPE NUTYEL
C        --------------------------------------
 401     DO 450 NUELEM = 1, NBELEM
C
C           LES NOEUDS DE L'ELEMENT FINI
            CALL EFNOEU( MNELE, NUELEM, NBNDEL, MCN(MNAUX) )
C
CCCC           LES POINTS GEOMETRIQUES DE L'ELEMENT FINI
CCC            CALL EFPOGE( MNELE, NUELEM, NBPGEF, MCN(MNPOEF) )
C
C           LE NUMERO DE VOLUME  DE L'EF
C           LE NUMERO DE SURFACE DES FACES   DE L'EF
C           LE NUMERO DE LIGNE   DES ARETES  DE L'EF
C           LE NUMERO DE POINT   DES SOMMETS DE L'EF
C           ----------------------------------------
            CALL EFPLSV( MNELE , NUELEM ,
     %                   NOVCEL, NOSFEL, NOLAEL, NOPSEL,
     %                   MCN(MNVC), MCN(MNSF), MCN(MNLA), MCN(MNPS) )
CCCC
CCCC           LE NOMBRE ET LES NO EVENTUELS DES TANGENTES AUX ARETES DE L'EF
CCC            CALL EFTGEF( MNELE,  NUELEM,
CCC     %                   NBTGEL, MCN(MNTGEL) )
CCCC           NBTGEL : =0 SI L'EF N'A PAS DE TANGENTE
CCCC                    >0 SI L'EF A NBTGEL TANGENTES
CCCC           ( 0 : '0 tg par sommet' ,      2 : '2 tg par arete' ,
CCCC             6 : '6 tg par triangle' ,    8 : '8 tg par quadrangle' ,
CCCC            12 : '12 tg par tetraedre' , 18 : '18 tg par pentaedre' ,
CCCC            24 : '24 tg par hexaedre' )
CCCC           NOTGEL : +-NUMERO DES NBTGEL TANGENTES DE L'EF NUELEM
CCCC
            IF( NUTYEL .EQ.  1 .OR.
     %          NUTYEL .EQ. 13 .OR.
     %          NUTYEL .EQ. 29 ) THEN
C
C              UN EF FREEFEM TRIANGLE DE PLUS
               NBEF = NBEF + 1
C
C              LES 3 SOMMETS DU TRIANGLE
               MCN( MNLT + 1 ) = MCN(MNAUX1 + 1 )
               MCN( MNLT + 2 ) = MCN(MNAUX1 + 2 )
               MCN( MNLT + 3 ) = MCN(MNAUX1 + 3 )
C
C              LE LABEL OU NUMERO DE SURFACE DU TRIANGLE
               MCN( MNLT + 4 ) = MCN(MNSF1  + 1 )
               MNLT = MNLT + 4
C
C              LES EDGES on BOUNDARY DU TRIANGLE
               DO 410 L=1,3
                  IF( MCN(MNLA1+L) .GT. 0 ) THEN
C
C                    UN EDGE DE PLUS
C                   (ON SUPPOSE VOIR UNE SEULE FOIS UN EDGE BOUNDARY)
                     NBEDGE = NBEDGE + 1
C
C                    LE NO DES 2 SOMMETS DE L'ARETE FRONTALIERE
                     LL = MOD(L,3) + 1
                     MCN( MNLE + 1 ) = MCN(MNAUX1 + L  )
                     MCN( MNLE + 2 ) = MCN(MNAUX1 + LL )
C
C                    LE LABEL OU NO DE LIGNE DE L'ARETE
                     MCN( MNLE + 3 ) = MCN(MNLA1 + L )
                     MNLE = MNLE + 3
C
C                    LE NO DE FRONTIERE AUX 2 SOMMETS EXTREMITES
C                    VALEUR ECRASEE ENSUITE SI LE SOMMET EST UN POINT
                     MCN( MNLV + MCN(MNAUX1+L ) ) = MCN(MNLA1 + L)
                     MCN( MNLV + MCN(MNAUX1+LL) ) = MCN(MNLA1 + L)
C
                  ENDIF
 410           CONTINUE
C
C              LE NO DE POINT DES 3 VERTICES DU TRIANGLE
               DO 414 L=1,3
                  IF( MCN(MNPS1+L) .GT. 0 ) THEN
C                    LE LABEL OU -NO DE POINT DU SOMMET
                     MCN( MNLV + MCN(MNAUX1+L) ) = -MCN(MNPS1 + L )
                  ENDIF
 414           CONTINUE
C
            ELSE IF( NUTYEL .EQ.  3 .OR.
     %               NUTYEL .EQ. 16  ) THEN
C
C              QUADRANGLE AVEC 4 SOMMETS => 2 TRIANGLES
C              CHOIX DU DECOUPAGE EN 2 TRIANGLES POUR MAXIMISER
C              LE MINIMUM DES QUALITES
               MNXY = MNXYZN + WYZNOE
C              TRIANGLE 123
               CALL QUATRI( MCN(MNAUX), RMCN(MNXY), QUA123 )
C              TRIANGLE 134
               NOSOTR(1) = MCN(MNAUX1+1)
               NOSOTR(2) = MCN(MNAUX1+3)
               NOSOTR(3) = MCN(MNAUX1+4)
               CALL QUATRI( NOSOTR, RMCN(MNXY), QUA134 )
C
C              TRIANGLE 124
               NOSOTR(2) = MCN(MNAUX1+2)
               CALL QUATRI( NOSOTR, RMCN(MNXY), QUA124 )
C              TRIANGLE 234
               NOSOTR(1) = MCN(MNAUX1+2)
               NOSOTR(2) = MCN(MNAUX1+3)
               CALL QUATRI( NOSOTR, RMCN(MNXY), QUA234 )
C
               IF(MIN(QUA123,QUA134)+0.001 .GE. MIN(QUA124,QUA234)) THEN
C                 +0.001 POUR FAVORISER UN CHOIX STABLE SI QUALITES EGALES
C
C                 UN EF FREEFEM TRIANGLE DE PLUS
                  NBEF = NBEF + 1
C                 LES 3 NO DES SOMMETS DU TRIANGLE 1 2 3
                  MCN( MNLT + 1 ) = MCN(MNAUX1 + 1 )
                  MCN( MNLT + 2 ) = MCN(MNAUX1 + 2 )
                  MCN( MNLT + 3 ) = MCN(MNAUX1 + 3 )
C                 LE LABEL OU NUMERO DE SURFACE DU TRIANGLE
                  MCN( MNLT + 4 ) = MCN(MNSF1  + 1 )
                  MNLT = MNLT + 4
C
C                 UN EF FREEFEM TRIANGLE DE PLUS
                  NBEF = NBEF + 1
C                 LES 3 NO DES SOMMETS DU TRIANGLE 1 3 4
                  MCN( MNLT + 1 ) = MCN(MNAUX1 + 1 )
                  MCN( MNLT + 2 ) = MCN(MNAUX1 + 3 )
                  MCN( MNLT + 3 ) = MCN(MNAUX1 + 4 )
C                 LE LABEL OU NUMERO DE SURFACE DU TRIANGLE
                  MCN( MNLT + 4 ) = MCN(MNSF1  + 1 )
                  MNLT = MNLT + 4
C
               ELSE
C
C                 UN EF FREEFEM TRIANGLE DE PLUS
                  NBEF = NBEF + 1
C                 LES 3 NO DES SOMMETS DU TRIANGLE 1 2 4
                  MCN( MNLT + 1 ) = MCN(MNAUX1 + 1 )
                  MCN( MNLT + 2 ) = MCN(MNAUX1 + 2 )
                  MCN( MNLT + 3 ) = MCN(MNAUX1 + 4 )
C                 LE LABEL OU NUMERO DE SURFACE DU TRIANGLE
                  MCN( MNLT + 4 ) = MCN(MNSF1  + 1 )
                  MNLT = MNLT + 4
C
C                 UN EF FREEFEM TRIANGLE DE PLUS
                  NBEF = NBEF + 1
C                 LES 3 NO DES SOMMETS DU TRIANGLE 2 3 4
                  MCN( MNLT + 1 ) = MCN(MNAUX1 + 2 )
                  MCN( MNLT + 2 ) = MCN(MNAUX1 + 3 )
                  MCN( MNLT + 3 ) = MCN(MNAUX1 + 4 )
C                 LE LABEL OU NUMERO DE SURFACE DU TRIANGLE
                  MCN( MNLT + 4 ) = MCN(MNSF1  + 1 )
                  MNLT = MNLT + 4
C
               ENDIF
C
C              LES BOUNDARY EDGES DU QUADRANGLE
               DO 420 L=1,4
                  IF( MCN(MNLA1+L) .GT. 0 ) THEN
C
C                    UNE ARETE FRONTALIERE DE PLUS
C                    (ON SUPPOSE VOIR UNE SEULE FOIS UN BOUNDARY EDGE)
                     NBEDGE = NBEDGE + 1
C
C                    LE NO DES 2 SOMMETS DE L'ARETE L
                     LL = MOD(L,4) + 1
                     MCN( MNLE + 1 ) = MCN(MNAUX1 + L  )
                     MCN( MNLE + 2 ) = MCN(MNAUX1 + LL )
C                    LE LABEL OU NO DE LIGNE DE L'ARETE
                     MCN( MNLE + 3 ) = MCN(MNLA1 + L )
                     MNLE = MNLE + 3
C
C                    LE NO DE LIGNE FRONTIERE AUX 2 SOMMETS EXTREMITES
C                    VALEUR ECRASEE SI LE SOMMET EST UN POINT
                     MCN( MNLV + MCN(MNAUX1+L ) ) = MCN(MNLA1 + L)
                     MCN( MNLV + MCN(MNAUX1+LL) ) = MCN(MNLA1 + L)
C
                  ENDIF
 420           CONTINUE
C
C              LE NO DE POINT DES 4 VERTICES DU QUADRANGLE
               DO 424 L=1,4
                  IF( MCN(MNPS1+L) .GT. 0 ) THEN
C                    LE LABEL OU -NO DE POINT DU SOMMET
                     MCN( MNLV + MCN(MNAUX1+L) ) = -MCN(MNPS1 + L )
                  ENDIF
 424           CONTINUE
C
            ELSE IF( NUTYEL .EQ.  19 ) THEN
C
C              UN EF FREEFEM TETRAEDRE DE PLUS
               NBEF = NBEF + 1
C
C              LES 4 SOMMETS DU TETRAEDRE
               MCN( MNLT + 1 ) = MCN(MNAUX1 + 1 )
               MCN( MNLT + 2 ) = MCN(MNAUX1 + 2 )
               MCN( MNLT + 3 ) = MCN(MNAUX1 + 3 )
               MCN( MNLT + 4 ) = MCN(MNAUX1 + 4 )
C
C              LE LABEL OU NUMERO DE VOLUME DU TETRAEDRE
               MCN( MNLT + 5 ) = MCN(MNVC1  + 1 )
               MNLT = MNLT + 5
C
C              LES FACES on BOUNDARY DU TETRAEDRE
               DO L=1,4
                  IF( MCN(MNSF1+L) .GT. 0 ) THEN
C
C                    UNE FACE FRONTALIERE DE PLUS
C                   (ON SUPPOSE VOIR UNE SEULE FOIS UNE FACE BOUNDARY)
                     NBEDGE = NBEDGE + 1
C
C                    LE NO DES 3 SOMMETS
                     MCN( MNLE + 1 ) = MCN( MNAUX1 + NOSOFA(1,L) )
                     MCN( MNLE + 2 ) = MCN( MNAUX1 + NOSOFA(2,L) )
                     MCN( MNLE + 3 ) = MCN( MNAUX1 + NOSOFA(3,L) )
C
C                    LE LABEL OU NO DE SURFACE DE LA FACE
                     MCN( MNLE + 4 ) = MCN(MNSF1 + L )
                     MNLE = MNLE + 4
C
C                    LE NO DE SURFACE AUX 3 SOMMETS DE LA FACE
C                    VALEUR ECRASEE SI LE SOMMET EST UN POINT
                     DO LL=1,3
                        MCN(MNLV+MCN(MNAUX1+NOSOFA(LL,L)))=MCN(MNSF1+L)
                     ENDDO
C
                  ENDIF
               ENDDO
C
C              LE NO DE POINT DES 4 VERTICES DU TETRAEDRE
               DO L=1,4
                  IF( MCN(MNPS1+L) .GT. 0 ) THEN
C                    LE LABEL OU -NO DE POINT DU SOMMET
                     MCN( MNLV + MCN(MNAUX1+L) ) = -MCN( MNPS1+L )
                  ENDIF
               ENDDO
C
            ENDIF
 450     CONTINUE
 500  CONTINUE
C
C     DECLARATION OUVERTURE DU FICHIER NOM_OBJET.freefem.msh
C     ======================================================
C     LE NOM DU FICHIER
      I = NUDCNB( KNOMOB )
ccc      KNMFIC = 'freefem' // '.' // KNOMOB(1:I)
      KNMFIC = KNOMOB(1:I) // '.freefem.msh'
      WRITE(IMPRIM,*) 'Creation du fichier: ', KNMFIC
C
C     SI LE FICHIER KNMFIC EXISTE ALORS IL EST DETRUIT PUIS RECONSTRUIT
      INQUIRE( FILE=KNMFIC, EXIST=LEXIST, OPENED=LOPEN )
      IF( LEXIST ) THEN
C        LE FICHIER KNMFIC EXISTE
         IF( .NOT. LOPEN ) THEN
C           OUVERTURE DU FICHIER KNMFIC
            CALL TRUNIT( NF )
            OPEN( FILE=KNMFIC, UNIT=NF, STATUS='OLD' )
         ENDIF
C        DESTRUCTION DU FICHIER
         CLOSE( NF, STATUS='DELETE' )
      ENDIF
C
C     CREATION DU FICHIER KNMFIC
      CALL TRUNIT( NF )
      OPEN( UNIT=NF, ERR=9900, STATUS='NEW',
     %      FILE=KNMFIC, ACCESS='SEQUENTIAL', FORM='FORMATTED' )
C
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10001) KNMFIC
10001 FORMAT(/' CREATION du FICHIER ',A)
      ELSE
         WRITE(IMPRIM,10002) KNMFIC
10002 FORMAT(/' CREATION of FILE ',A)
      ENDIF
C
C     ECRITURE DES DONNEES DU MAILLAGE SUR LE FICHIER freefem.msh
C     ===========================================================
C     NOMBRE DE VERTICES, TRIANGLES, EDGES BOUNDARY
CCC      WRITE(IMPRIM,12010) NDIM
CCC      WRITE(NF,    12010) NDIM
CCC      WRITE(IMPRIM,12010) NBNOEU, NBEF, NBEDGE
      WRITE(NF,12010) NBNOEU, NBEF, NBEDGE
12010 FORMAT( 10I10 )
C
C     LES COORDONNEES DE CHAQUE VERTICE  X Y (Z) LABEL
      MN = MNXYZN + WYZNOE
      IF( NDIM .EQ. 2 ) THEN
         DO I=1,NBNOEU
CCC            WRITE(IMPRIM,12020) RMCN(MN), RMCN(MN+1), MCN(MNLABV-1+I)
            WRITE(NF,    12020) RMCN(MN), RMCN(MN+1), MCN(MNLABV-1+I)
            MN = MN + 3
         ENDDO
      ELSE
         DO I=1,NBNOEU
CCC            WRITE(IMPRIM,12023) (RMCN(MN+L),L=0,2), MCN(MNLABV-1+I)
            WRITE(NF,    12023) (RMCN(MN+L),L=0,2), MCN(MNLABV-1+I)
            MN = MN + 3
         ENDDO
      ENDIF
12020 FORMAT( 2E15.7, I10 )
12023 FORMAT( 3E15.7, I10 )
C
C     LES TRIANGLES ou TETRAEDRES  NS1 ... NSnbstef  LABEL
      MN = MNLABT - 1
      DO 2030 I=1,NBEF
CCC         WRITE(IMPRIM,12010) (MCN(MN+L),L=1,NBSTEF+1)
         WRITE(NF,    12010) (MCN(MN+L),L=1,NBSTEF+1)
         MN = MN + NBSTEF+1
 2030 CONTINUE
C
C     LES EDGES ou FACES NS1 ... NSnbstef-1 LABEL
      MN = MNLABE - 1
      DO 2040 I=1,NBEDGE
CCC         WRITE(IMPRIM,12010) (MCN(MN+L),L=1,NBSTEF)
         WRITE(NF,    12010) (MCN(MN+L),L=1,NBSTEF)
         MN = MN + NBSTEF
 2040 CONTINUE
C
C     ECRITURE DU NOM DE L'OBJET SUR LE FICHIER freefem.NOM_PLSV
      IF( LANGAG .EQ. 0 ) THEN
CCC         WRITE(IMPRIM,14050)
         WRITE(NF,14050)
         WRITE(NF,14059) KNOMOB
         WRITE(NF,14060) NDIM
      ELSE
CCC         WRITE(IMPRIM,13050)
         WRITE(NF,13050)
         WRITE(NF,13059) KNOMOB
         WRITE(NF,13060) NDIM
      ENDIF
14050 FORMAT(/,'FIN DU FICHIER FREEFEM: SUITE NON LUE PAR FREEFEM',/)
13050 FORMAT(/,'END of FILE FREEFEM: UNDER NOT READ BY FREEFEM',/)
14059 FORMAT(A24, '   {NOM DE L''OBJET}')
13059 FORMAT(A24, '   {NAME of the OBJECT}')
14060 FORMAT(I10, '   {DIMENSION 2 ou 3 DE L''ESPACE DE L''OBJET}' )
13060 FORMAT(I10, '   {SPACE''s DIMENSION (2 or 3) of the OBJECT}' )
C
C     ECRITURE DU NOMBRE DE SOMMETS DE L'OBJET
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(NF,14061) NBNOEU
         WRITE(NF,14062) NBTGN
      ELSE
         WRITE(NF,13061) NBNOEU
         WRITE(NF,13062) NBTGN
      ENDIF
14061 FORMAT( I10, '   {NBNOEUDS NOMBRE DE NOEUDS DE L''OBJET}' )
14062 FORMAT( I10, '   {NBTG     NOMBRE DE TANGENTES DE L''OBJET}' )
13061 FORMAT( I10, '   {NBNOEUDS NUMBER of NODES of the OBJECT}' )
13062 FORMAT( I10, '   {NBTG     NUMBER of TANGENTS of the OBJECT}' )
C
C     ECRITURE DU NOMBRE D'OBJETS INTERNES OU MATERIAUX
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(NF,14100) NBOBIN
      ELSE
         WRITE(NF,13100) NBOBIN
      ENDIF
14100 FORMAT( I10, '   {NBOBIN NOMBRE DE MATERIAUX}' )
13100 FORMAT( I10, '   {NBOBIN NUMBER of MATERIALS}' )
C
C     CORRESPONDANCE LABEL <=> NOM de MATERIAU
C     L'ADRESSE MCN DU DEBUT DU TABLEAU NUOBIN
      MNOBIN = MNTOPO + WMTYEL + NBTYEL
      MN     = MNOBIN - 2
      DO 600 I=1,NBOBIN
C        LE TYPE DE L'OBJET
         MN   = MN + 2
         NYOB = MCN( MN )
         NUOB = MCN( MN + 1 )
C        OUVERTURE DE L'OBJET
         CALL LXNLOU( NTMN(NYOB), NUOB, NTOB, MNOB )
         KNM = NMTYOB( NYOB )
         CALL NMOBNU( KNM, NUOB, KNOM )
C        ECRITURE DU NOM DE L'OBJET INTERNE I
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(NF,14602) NUOB, KNM, KNOM
         ELSE
            WRITE(NF,13602) NUOB, KNM, KNOM
         ENDIF
 600  CONTINUE
14602 FORMAT( I10, '   {est le NUMERO dans ',A10,' de } ',A )
13602 FORMAT( I10, '   {is the NUMBER in ',A10,' of } ',A )
C
C     ECRITURE DU NOMBRE D'OBJETS AUX LIMITES DE L'OBJET
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(NF,14620) NBOBCL
      ELSE
         WRITE(NF,13620) NBOBCL
      ENDIF
14620 FORMAT( I10, '   {NBOBCL NOMBRE D''OBJETS AUX LIMITES}' )
13620 FORMAT( I10, '   {NBOBCL NUMBER of PL(S) at the BOUNDARY}' )
C
C     BOUCLE SUR LES SURFACES(SI NDIM=3) LIGNES ET POINTS
      DO 650 L=NDIM,1,-1
C        ECRITURE DU NOMBRE D'OBJETS AUX LIMITES DE TYPE L
         KNM = NMTYOB( L )
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(NF,14640) NBPLSV(L),KNM
         ELSE
            WRITE(NF,13640) NBPLSV(L),KNM
         ENDIF
14640    FORMAT( I10, '   {NOMBRE DE ',A10,'}' )
13640    FORMAT( I10, '   {NUMBER of ',A10,'}' )
         MN  = MNOBCL - 2
         DO 640 I=1,NBOBCL
C           LE TYPE DE L'OBJET
            MN   = MN + 2
            NYOB = MCN( MN )
            IF( NYOB .NE. L ) GOTO 640
            NUOB = MCN( MN + 1 )
C           OUVERTURE DE L'OBJET
            CALL LXNLOU( NTMN(NYOB), NUOB, NTOB, MNOB )
            CALL NMOBNU( KNM, NUOB, KNOM )
C           ECRITURE DU NOM DE L'OBJET AUX LIMITES
            IF( LANGAG .EQ. 0 ) THEN
               WRITE(NF,10102) NUOB, KNM, KNOM
            ELSE
               WRITE(NF,11102) NUOB, KNM, KNOM
            ENDIF
 640     CONTINUE
 650  CONTINUE
C
C     ECRITURE DU NOMBRE D'EF
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(NF,14801) NBEF, KNOMOB
      ELSE
         WRITE(NF,13801) NBEF, KNOMOB
      ENDIF
14801 FORMAT(I10,'   {NOMBRE TOTAL D''EF de l''OBJET ',A,'}')
13801 FORMAT(I10,'   {TOTAL NUMBER of FE of the OBJECT ',A,'}')
C
C     FERMETURE DU FICHIER freefem
      CLOSE( NF )
C
C     DESTRUCTION DES TABLEAUX AUXILIAIRES
 9900 IF(MNAUX .GT.0) CALL TNMCDS( 'ENTIER', 27+8+12+6+1+24,    MNAUX )
      IF(MNLABT.GT.0) CALL TNMCDS( 'ENTIER', (NBSTEF+1)*MXTRIA, MNLABT)
      IF(MNLABT.GT.0) CALL TNMCDS( 'ENTIER', NBSTEF*MXEDGE,     MNLABE)
      IF(MNLABT.GT.0) CALL TNMCDS( 'ENTIER', NBNOEU,            MNLABV)
      IF(MNELEM.GT.0) CALL TNMCDS( 'ENTIER', 2*MXTYEL,          MNELEM)
      RETURN
C
C     PB A L'OUVERTURE DU FICHIER
      NBLGRC(NRERR) = 1
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = 'IMPOSSIBLE OUVERTURE DU FICHIER freefem'
      ELSE
         KERR(1) = 'FILE freefem CAN NOT BE OPENED'
      ENDIF
      CALL LEREUR
      GOTO 10
C
      END
