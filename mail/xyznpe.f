      SUBROUTINE XYZNPE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : SORTIR SUR UN FICHIER LES 3 COORDONNEES DES SOMMETS ET TANGENTES
C ----- ET LES NUMEROS DES NOEUDS DES VOLUMES SURFACES LIGNES POINTS
C       DES EF D'UN OBJET DANS UN FICHIER FORME DE CARACTERES ASCII
C
C ATTENTION: NOEUDS et POINTS SONT SUPPOSES IDENTIQUES!
C
C SORTIE :
C --------
C LE FICHIER xyznpef DANS LE REPERTOIRE DU PROJET
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     OCTOBRE 1997
C23456---------------------------------------------------------------012
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
      COMMON / UNITES /  LECTEU, IMPRIM, INTERA, NUNITE(29)
      REAL               RMCN(1)
      EQUIVALENCE       (MCN(1),RMCN(1))
      CHARACTER*4        NOMELE(2)
      CHARACTER*10       NMTYOB,KNM
      CHARACTER*24       KNOMOB,KNOM
      CHARACTER*32       KNMFIC
      LOGICAL            LEXIST,LOPEN
      INTEGER            NBPLSV(4)
C
C     NOM DE L'OBJET A TRAITER
 10   CALL INVITE( 45 )
      IERR   = 0
      NCVALS = 0
      CALL LIRCAR( NCVALS, KNOMOB )
      IF( NCVALS .EQ. -1 ) RETURN
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
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10001) KNOMOB
10001 FORMAT(/' CREATION du FICHIER xyznpef.',A)
      ELSE
         WRITE(IMPRIM,10002) KNOMOB
10002 FORMAT(/' CREATION of the FILE xyznpef.',A)
      ENDIF
C
C     RECHERCHE DU TABLEAU DEFINITION
      CALL LXTSOU( NTLXOB, 'DEFINITION', NTDFOB, MNDFOB )
C     S'IL N'EXISTE PAS RETOUR A LA DEMANDE DU NOM DE L'OBJET
      IF( NTDFOB .LE. 0 ) THEN
         NBLGRC(NRERR) = 2
         KERR(1) = KNOMOB
         IF( LANGAG .EQ. 0 ) THEN
            KERR(2) = 'OBJET SANS TMS DEFINITION'
         ELSE
            KERR(2) = 'OBJECT WITHOUT TMS DEFINITION'
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
C     =====================================================================
      CALL NDPGEL( NTLXOB, NTTOPO, MNTOPO ,
     %             NTXYZP, MNXYZP, NTXYZN, MNXYZN ,
     %             NBTYEL, MCN(MNTELE), MCN(MNELEM), IERR )
      IF( IERR .NE. 0 ) GOTO 9999
      NBTGN = MCN( MNXYZN + WNBTGN )
      NBTGP = MCN( MNXYZP + WNBTGP )
      print *,'xyznpe: Nombre de Tangentes aux NOEUDS NBTGN=',NBTGN
      print *,'xyznpe: Nombre de Tangentes aux POINTS NBTGP=',NBTGP
C
C     LE TYPE DE MAILLAGE ET LES OBJETS INTERNES ET AUX LIMITES
      NDPGST = MCN( MNTOPO + WDPGST )
      NBOBIN = MCN( MNTOPO + WBOBIN )
      NBOBCL = MCN( MNTOPO + WBOBCL )
C
C     LA VALEUR DE Z DES POINTS DU MAILLAGE DE L'OBJET DEFINIT
C     NDIM LA DIMENSION 1 OU 2 OU 3 OU 6 DE L'ESPACE DES COORDONNEES
C     ATTENTION: CETTE PROGRAMMATION SOUS ENTEND NOEUDS=POINTS
      NBCOOR = MCN( MNXYZP + WBCOOP )
      IF( NBCOOR .EQ. 3 ) THEN
         CALL DIMCOO( MCN(MNXYZP+WNBPOI), MCN(MNXYZP+WYZPOI), NDIM )
         NDIMV = NDIM
      ELSE IF( NBCOOR .EQ. 6 ) THEN
         NDIM  = 6
         NDIMV = 3
      ELSE
         WRITE(IMPRIM,*) 'xyznpe: PROBLEME avec NBCOOR=',NBCOOR
         RETURN
      ENDIF
C
C     SANS ou AVEC les TANGENTES dans le FICHIER?
C     ===========================================
      CALL LIMTCL( 'tgoupas', NOTG )
      IF( NOTG .LT. 0 ) RETURN
C
C     LE NOMBRE TOTAL DE NOEUDS DU MAILLAGE DE L'OBJET
      NBNOEU = MCN( MNXYZN + WNBNOE )
C     LE NOMBRE TOTAL DE TANGENTES DU MAILLAGE
      IF( NOTG .GT. 0 ) THEN
         NBTGN = MCN( MNXYZN + WNBTGN )
      ELSE
         NBTGN = 0
      ENDIF
C
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10210) NDIM,NBNOEU,NBTGN,NBOBIN,NBOBCL,NBTYEL
10210 FORMAT(' DIMENSION de L''ESPACE'       ,T32,'=',I6/
     %' NOMBRE de NOEUDS'                    ,T32,'=',I6/
     %' NOMBRE de TANGENTES'                 ,T32,'=',I6/
     %' NOMBRE d''OBJETS INTERNES'           ,T32,'=',I6/
     %' NOMBRE d''OBJETS AUX LIMITES'        ,T32,'=',I6/
     %' NOMBRE de TYPES d''EF'               ,T32,'=',I6)
      ELSE
         WRITE(IMPRIM,11210) NDIM,NBNOEU,NBTGN,NBOBIN,NBOBCL,NBTYEL
11210 FORMAT(' SPACE''s DIMENSION'         ,T33,'=',I6/
     %' NUMBER of NODES'                   ,T33,'=',I6/
     %' NUMBER of TANGENTS'                ,T33,'=',I6/
     %' NUMBER of MATERIALS'               ,T33,'=',I6/
     %' NUMBER of PL(S) at the BOUNDARY'   ,T33,'=',I6/
     %' NUMBER of TYPES of FE'             ,T33,'=',I6)
      ENDIF
C
C     DECLARATION OUVERTURE DU FICHIER xyznpef.NOM_PLSV
C     =================================================
C     LE NOM DU FICHIER
      I = NUDCNB( KNOMOB )
      KNMFIC = 'xyznpef' // '.' // KNOMOB(1:I)
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
      OPEN( UNIT=NF, ERR=9998, STATUS='NEW',
     %      FILE=KNMFIC, ACCESS='SEQUENTIAL', FORM='FORMATTED' )
C
C     ADRESSE-NBCOOR DE LA 1-ERE COORDONNEE DU 1-ER SOMMET DU TMS 'XYZNOEUD'
      MNS = MNXYZN + WYZSOM - NBCOOR
C
C     ECRITURE DU NOM DE L'OBJET SUR LE FICHIER xyznpef.NOM_PLSV
C     ==========================================================
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(NF,10059) KNOMOB
         WRITE(NF,13060) NBCOOR
         WRITE(NF,10060) NDIM
      ELSE
         WRITE(NF,11059) KNOMOB
         WRITE(NF,14060) NBCOOR
         WRITE(NF,11060) NDIM
      ENDIF
10059 FORMAT(A24,'   {NOM de L''OBJET}')
11059 FORMAT(A24,'   {NAME of the OBJECT}')
13060 FORMAT( I8,'   {NBCOOR NOMBRE de COORDONNEES D UN NOEUD (3 ou 6)}
     %')
14060 FORMAT( I8,'   {NBCOOR COORDINATE NUMBER of a NODE (3 or 6) }')
10060 FORMAT( I8,'   {DIMENSION de l''ESPACE de l''OBJET (1 2 3 ou 6)}')
11060 FORMAT( I8,'   {SPACE''s DIMENSION of the OBJECT (1 2 3 or 6)}')

C     ECRITURE DES COORDONNEES DES NOEUDS ET TANGENTES DE L'OBJET
C     ===========================================================
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(NF,10061) NBNOEU
         WRITE(NF,10062) NBTGN
      ELSE
         WRITE(NF,11061) NBNOEU
         WRITE(NF,11062) NBTGN
      ENDIF
10061 FORMAT( I8, '   {NBNOEUDS NOMBRE de NOEUDS de l''OBJET}' )
10062 FORMAT( I8, '   {NBTG     NOMBRE de TANGENTES de l''OBJET}' )
11061 FORMAT( I8, '   {NBNOEUDS NUMBER of NODES of the OBJECT}' )
11062 FORMAT( I8, '   {NBTG     NUMBER of TANGENTS of the OBJECT}' )

C     ECRITURE DE XYZ DES NBNOEU NOEUDS
10065 FORMAT( 6E15.7 )
C     ADRESSE DE LA 1-ERE COORDONNEE DU 1-ER SOMMET DU TMS 'XYZNOEUD'
      MNS = MNXYZN + WYZNOE
      DO I=1,NBNOEU
         DO K=0,NBCOOR-1
C           MISE A ZERO DE LA VALEUR ABSOLUE de XYZ < 1E-20
            IF( ABS( RMCN(MNS+K) ) .LE. 1E-20 ) RMCN(MNS+K)=0.
         ENDDO
         WRITE(NF,10065) (RMCN(MNS+K),K=0,NBCOOR-1)
         MNS = MNS + NBCOOR
      ENDDO

C     ECRITURE DE XYZ DES NBTGN TANGENTES AUX ARETES DU MAILLAGE
C     MNS ADRESSE DE LA 1-ERE COORDONNEE DE LA 1-ERE TG DU TMS 'XYZNOEUD'
      IF( NBTGN .GT. 0 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(NF,10066) (RMCN(MNS+K),K=0,NBCOOR-1)
         ELSE
            WRITE(NF,11066) (RMCN(MNS+K),K=0,NBCOOR-1)
         ENDIF
10066    FORMAT( 3E15.7, ' {XYZ de la PREMIERE TANGENTE}' )
11066    FORMAT( 3E15.7, ' {XYZ of the FIRST TANGENT VECTOR}' )
         MNS = MNS + NBCOOR
         DO 66 I=2,NBTGN
            WRITE(NF,10065) (RMCN(MNS+K),K=0,NBCOOR-1)
            MNS = MNS + NBCOOR
 66      CONTINUE
      ENDIF

C     BOUCLE SUR LES OBJETS "INTERNES"
C     ================================
C     ECRITURE DU NOMBRE D'OBJETS INTERNES OU MATERIAUX
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(NF,10100) NBOBIN
      ELSE
         WRITE(NF,11100) NBOBIN
      ENDIF
10100 FORMAT( I8, '   {NBOBIN NOMBRE de MATERIAUX}' )
11100 FORMAT( I8, '   {NBOBIN NUMBER of MATERIALS}' )
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
C        ECRITURE DU NOM DE L'OBJET INTERNE I
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(NF,10102) NUOB, KNM
         ELSE
            WRITE(NF,11102) NUOB, KNM
         ENDIF
         WRITE(NF,10103) KNOM
 100  CONTINUE
10102 FORMAT( I8, '   {est le NUMERO dans ',A10,' de }' )
11102 FORMAT( I8, '   {is the NUMBER in ',A10,' of }' )
10103 FORMAT( A24 )
      IF( IERR .NE. 0 ) GOTO 9999
      NBPLSV(NDIMV+1) = NBOBIN
C
C     BOUCLE SUR LES OBJETS AUX LIMITES DE L'OBJET
C     RECHERCHE DU NOMBRE DE POINTS LIGNES SURFACES DE L'OBJET
C     ========================================================
C     ECRITURE DU NOMBRE D'OBJETS AUX LIMITES DE L'OBJET
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(NF,10120) NBOBCL
      ELSE
         WRITE(NF,11120) NBOBCL
      ENDIF
10120 FORMAT( I8, '   {NBOBCL NOMBRE d''OBJETS AUX LIMITES}' )
11120 FORMAT( I8, '   {NBOBCL NUMBER of PL(S) at the BOUNDARY}' )
C
C     RECENSEMENT DU NOMBRE DE POINTS, LIGNES, ... DE L'OBJET
C     BOUCLE SUR LES SURFACES(SI NDIMV=3) et LIGNES et POINTS
C     L'ADRESSE MCN DU DEBUT DU TABLEAU NUOBCL
      MNOBCL = MNOBIN + MOTVAR(13) * NBOBIN
      DO 130 L=NDIMV,1,-1
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
C     BOUCLE SUR LES SURFACES(SI NDIMV=3) et LIGNES et POINTS
      DO 150 L=NDIMV,1,-1
C        ECRITURE DU NOMBRE D'OBJETS AUX LIMITES DE TYPE L
         KNM = NMTYOB( L )
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(NF,10140) NBPLSV(L),KNM
         ELSE
            WRITE(NF,11140) NBPLSV(L),KNM
         ENDIF
10140    FORMAT( I8, '   {NOMBRE de ',A10,'}' )
11140    FORMAT( I8, '   {NUMBER of ',A10,'}' )
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
               WRITE(NF,10102) NUOB, KNM
            ELSE
               WRITE(NF,11102) NUOB, KNM
            ENDIF
            WRITE(NF,10143) KNOM
 140     CONTINUE
 150  CONTINUE
10143 FORMAT( A24 )
      IF( IERR .GT. 0 ) GOTO 9999
C
C     RESERVATION D'UN TABLEAU AUXILIAIRE POUR LES NOEUDS,..., TGS D'UN EF
      MNAUX = 0
      CALL TNMCDC( 'ENTIER', 64+64+12+6+1+24, MNAUX )
      MNAUX1 = MNAUX - 1
      MNPS   = MNAUX + 64
      MNPS1  = MNPS  - 1
      MNLA   = MNPS  + 64
      MNLA1  = MNLA  - 1
      MNSF   = MNLA  + 12
      MNSF1  = MNSF  - 1
      MNVC   = MNSF  + 6
      MNVC1  = MNVC  - 1
      MNTGEL = MNVC  + 1
      MNTGE1 = MNTGEL -1
C
C     LA BOUCLE SUR LES DIFFERENTS TYPES D'EF (TABLEAUX NPEF"xxxx)
C     ============================================================
C     ECRITURE DU NOMBRE DE TYPE D'EF
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(NF,10400) NBTYEL
      ELSE
         WRITE(NF,11400) NBTYEL
      ENDIF
10400 FORMAT( I8, '   {NBTYEL NOMBRE de TYPES d''EF de l''OBJET}' )
11400 FORMAT( I8, '   {NBTYEL NUMBER of TYPES of FE of the OBJECT}' )
      NBEF   = 0
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
C        LE NOMBRE DE SOMMETS DE CE TYPE D'EF
         NBSOE = NBSOME(NCOGEL)
C
C        LE NOMBRE DE VOLUMES DE CE TYPE D'EF
         NBVCEL = MCN( MNELE + WBVCEL )
         IF( NDIM .EQ. 6 ) THEN
            NBSOE = 0
            NBNSOM= 0
            NARET = 0
            NBFAC = 0
            NFACE = 0
            IF( NBVCEL .GT. 0 ) THEN
                NBVOL = 1
            ELSE
                NBVOL = 0
            ENDIF
         ELSEIF( NDIM .EQ. 3 ) THEN
            IF( NBVCEL .GT. 0 ) THEN
                NBVOL = 1
            ELSE
                NBVOL = 0
            ENDIF
            NBFAC = NFACE
         ELSEIF( NDIM .EQ. 2 ) THEN
C           EN DIMENSION 2 PAS DE VOLUME ET 1 FACE
            NBVOL = 0
            NBFAC = 1
         ELSE
C           EN DIMENSION 1 PAS DE VOLUME ET PAS DE FACE
            NBVOL = 0
            NBFAC = 0
         ENDIF
C
         IF( LANGAG .EQ. 0 ) THEN
         WRITE(NF,10401) NBELEM
10401    FORMAT( I8, '   {NBELEM NOMBRE d''EF de CE TYPE d''EF}' )
         WRITE(NF,10402) NBNOE
10402    FORMAT( I8, '   {NBNOE  NOMBRE de NOEUDS de CE TYPE d''EF}')
         WRITE(NF,10403) NBSOE
10403    FORMAT( I8, '   {NBSOE  NOMBRE de SOMMETS de CE TYPE d''EF}')
         WRITE(NF,10404) NARET
10404    FORMAT( I8, '   {NBARET NOMBRE d''ARETES de CE TYPE d''EF}')
         WRITE(NF,10405) NBFAC
10405    FORMAT( I8, '   {NBFACE NOMBRE de FACES de CE TYPE d''EF}')
         WRITE(NF,10406) NBVOL
10406    FORMAT( I8, '   {NBVOL  NOMBRE de VOLUME de CE TYPE d''EF}')
         ELSE
         WRITE(NF,11401) NBELEM
11401    FORMAT( I8, '   {NBELEM NUMBER of FE of THIS TYPE}' )
         WRITE(NF,11402) NBNOE
11402    FORMAT( I8, '   {NBNOE  NUMBER of NODES of a FE of THIS TYPE}')
         WRITE(NF,11403) NBSOE
11403    FORMAT( I8,
     %          '   {NBSOE  NUMBER of VERTICES of a FE of THIS TYPE}')
         WRITE(NF,11404) NARET
11404    FORMAT( I8, '   {NBARET NUMBER of EDGES of a FE of THIS TYPE}')
         WRITE(NF,11405) NBFAC
11405    FORMAT( I8, '   {NBFACE NUMBER of FACES of a FE of THIS TYPE}')
         WRITE(NF,11406) NBVOL
11406    FORMAT( I8,
     %   '   {NBVOL  NUMBER of VOLUMES of a FE of THIS TYPE}')
         ENDIF
C
         CALL AZEROI( 64, MCN(MNPS) )
         CALL AZEROI( 12, MCN(MNLA) )
         CALL AZEROI(  6, MCN(MNSF) )
         MCN( MNVC ) = 0
C
C        LES ADRESSAGES POUR RECUPERER LES INFORMATIONS
C        SELON LE TYPE DE L'ELEMENT FINI
         GOTO( 401, 401, 401, 401, 400, 400, 400, 400, 400, 400,
     &         400, 400, 401, 400, 401, 401, 400, 401, 401, 401,
     &         401, 401, 401, 401, 400, 400, 400, 401, 401, 401,
     %         401, 401, 401, 401                         ),NUTYEL
C
C        ERREUR
 400     NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'ELEMENT FINI '// NOMELE(1)
     &             // NOMELE(2) //' NON PROGRAMME'
         ELSE
            KERR(1) = 'FINITE ELEMENT '// NOMELE(1)
     &             // NOMELE(2) //' NOT PROGRAMMED'
         ENDIF
         CALL LEREUR
         GOTO 9999
C
C        --------------------------------------
C        LA BOUCLE SUR LES EF DE CE TYPE NUTYEL
C        --------------------------------------
 401     DO 450 NUELEM = 1, NBELEM
C
C           LES NOEUDS DE L'ELEMENT FINI
C           ----------------------------
            CALL EFNOEU( MNELE, NUELEM, NBNDEL, MCN(MNAUX) )
C
C           LES POINTS GEOMETRIQUES DE L'ELEMENT FINI
C           -----------------------------------------
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
C
C           LE NOMBRE ET LES NUMEROS EVENTUELS DES TANGENTES AUX ARETES DE L'EF
C           -------------------------------------------------------------------
            CALL EFTGEF( MNELE,  NUELEM,
     %                   NBTGEL, MCN(MNTGEL) )
C           NBTGEL : =0 SI L'EF N'A PAS DE TANGENTE
C                    >0 SI L'EF A NBTGEL TANGENTES
C           ( 0 : '0 tg par sommet' ,      2 : '2 tg par arete' ,
C             6 : '6 tg par triangle' ,    8 : '8 tg par quadrangle' ,
C            12 : '12 tg par tetraedre' , 18 : '18 tg par pentaedre' ,
C            24 : '24 tg par hexaedre' )
C           NOTGEL : +-NUMERO DES NBTGEL TANGENTES DE L'EF NUELEM
C
C           LE NOMBRE D'EF
            NBEF = NBEF + 1
C           ECRITURE SUR LE FICHIER xyznpef :
C           DES NUMEROS DES NOEUDS
C           DES POINTS DES SOMMETS
C           DES LIGNES DES ARETES
C           DES SURFACES DES FACES
C           DU VOLUME
            WRITE(NF,10500) (MCN(MNAUX1+L),L=1,NBNOE),
     %                      (MCN(MNPS1+L) ,L=1,NBSOE),
     %                      (MCN(MNLA1+L) ,L=1,NARET),
     %                      (MCN(MNSF1+L) ,L=1,NBFAC),
     %                      (MCN(MNVC1+L) ,L=1,NBVOL)
C
C           NOMBRE DE TANGENTES DE L'EF ET NUMEROS DES TANGENTES
            IF( NOTG .GT. 0 .AND. (NBTGN+NBTGP .GT. 0) ) THEN
               WRITE(NF,10500) NBTGEL,
     %                        (MCN(MNTGE1+L),L=1,NBTGEL)
            ENDIF
 450     CONTINUE
 500  CONTINUE
10500 FORMAT( 10I8 )
C
C     DESTRUCTION DU TABLEAU AUXILIAIRE
      CALL TNMCDS( 'ENTIER', 64+64+12+6+1+24, MNAUX )
C
C     ECRITURE DU NOMBRE D'EF
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(NF,10501) NBEF, KNOMOB
      ELSE
         WRITE(NF,11501) NBEF, KNOMOB
      ENDIF
10501 FORMAT( I8,' {NOMBRE TOTAL d''EF de l''OBJET ',A,'}' )
11501 FORMAT( I8,' {TOTAL NUMBER of FE of the OBJECT ',A,'}' )
C
C     FERMETURE DU FICHIER xyznpef
      CLOSE( NF )
      RETURN
C
C     PB A L'OUVERTURE DU FICHIER
 9998 NBLGRC(NRERR) = 1
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = 'IMPOSSIBLE OUVERTURE DU FICHIER xyznpef'
      ELSE
         KERR(1) = 'FILE xyznpef CAN NOT BE OPENED'
      ENDIF
      CALL LEREUR
      GOTO 10
C
C     DESTRUCTION DU FICHIER
 9999 CLOSE( NF, STATUS='DELETE' )
      GOTO 10
      END
