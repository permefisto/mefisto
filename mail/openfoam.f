      SUBROUTINE OpenFOAM
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : SORTIR LES FICHIERS NECESSAIRES POUR CONSTRUIRE UN MAILLAGE
C ----- POUR LE LOGICIEL OPENFOAM
C
C ATTENTION: L'OBJET DOIT ETRE TRIDIMENSIONNEL ET
C            INTERPOLE PAR DES POLYNOMES DE DEGRE 1
C
C SORTIE :
C --------
C LES FICHIERS openfoam.points  DANS LE REPERTOIRE DU PROJET
C              openfoam.faces
C              openfoam.owner
C              openfoam.neighbour
C              openfoam.boundary
C
C ATTENTION:
C ----------
C LE TMS a_objet__face est modifie en position 9
C NO TYPE EF  est devenu NO TYPE EF * 1 000 000 + NO de la SURFACE
C POUR LES FACES FRONTALIERES SUR LES SURFACES DE L'OBJET
C
C OpenFOAM demande que les faces du maillage 3D (uniquement)
C SOIENT RENUMEROTEES SELON L'ORDRE:
C 1:          FAEX(NF)>0 LES FACES EXISTANTES INTERNES APPARTENANT A 2 EF
C 2:          FAEX(NF)=-NOSF1 LES FACES DE LA SURFACE 1 DE L'OBJET
C 3:          FAEX(NF)=-NOSF2 LES FACES DE LA SURFACE 2 DE L'OBJET
C    ...
C 1+NBSURF:   FAEX(NF)=-NOSF NBSURF LES FACES DE LA SURFACE NBSURF DE L'OBJET
C 1+NBSURF+1: FAEX(NF)=- 999 000 LES FACES DANS 1 SEUL EF ET
C                       NON SUR UNE SURFACE DE L'OBJET
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET TEXAS A & M UNIVERSITY at QATAR     Mars 2011
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
      include"./incl/a___face.inc"
      include"./incl/a___texte.inc"
      include"./incl/msvaau.inc"
      include"./incl/ponoel.inc"
C
      include"./incl/pp.inc"
      COMMON             MCN(MOTMCN)
      COMMON / UNITES /  LECTEU,IMPRIM,INTERA,NUNITE(29)
      REAL               RMCN(1)
      EQUIVALENCE       (MCN(1),RMCN(1))
      CHARACTER*4        NOMELE(2)
      CHARACTER*10       NMTYOB,KNM
      CHARACTER*24       KNOMOB
      CHARACTER*24       KNOMSF
      CHARACTER*32       KNOM
      CHARACTER*32       KNMFIC
      LOGICAL            LEXIST,LOPEN
      INTEGER            NBPLSV(4)
      INTEGER            NGS(4), NGS1(4)
C
      WRITE(IMPRIM,*)
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,*)'Creation des fichiers openfoam.* du maillage'
      ELSE
         WRITE(IMPRIM,*)'Creation of openfoam.* files'
      ENDIF
C
C     INITIALISATIONS
      MNELEM = 0
      MNFAEX = 0
      MNFARE = 0
      MNNUSF = 0
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
C
C     LE TYPE DE MAILLAGE ET LES OBJETS INTERNES ET AUX LIMITES
      NDPGST = MCN( MNTOPO + WDPGST )
C
C     OpenFOAM utilise seuelement les SOMMETS
      IF( NDPGST .NE. 0 ) THEN
         NBLGRC(NRERR) = 2
         KERR(1) = KNOMOB
         IF( LANGAG .EQ. 0 ) THEN
            KERR(2) = 'INTERPOLATION DEGRE 1 UNIQUEMENT PERMISE'
         ELSE
            KERR(2) = 'ONLY 1 DEGREE POLYNOMIAL INTERPOLATION ACCEPTED'
         ENDIF
         CALL LEREUR
         RETURN
      ENDIF
C
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
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,*) 'openfoam: PROBLEME avec NBCOOR=',NBCOOR
         ELSE
            WRITE(IMPRIM,*) 'openfoam: PROBLEM with NBCOOR=',NBCOOR
         ENDIF
         RETURN
      ENDIF
C
C     OPENFOAM EST EN 3D UNIQUEMENT
      IF( NDIM .LT. 3 .OR. NDIM .GT. 3 ) THEN
         NBLGRC(NRERR) = 2
         KERR(1) = KNOM
         IF( LANGAG .EQ. 0 ) THEN
            KERR(2) = 'OBJET NON TRIDIMENSIONNEL'
         ELSE
            KERR(2) = 'NOT a 3D OBJECT'
         ENDIF
         CALL LEREUR
         RETURN
      ENDIF
C
C     SANS les TANGENTES dans le FICHIER
C     ==================================
      NOTG = 0
C
C     LE NOMBRE TOTAL DE NOEUDS DU MAILLAGE DE L'OBJET
      NBNOEU = MCN( MNXYZN + WNBNOE )
C
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10210) NDIM,NBNOEU,NBOBIN,NBOBCL,NBTYEL
10210 FORMAT(' DIMENSION DE L''ESPACE'  ,T32,'=',I10/
     %' NOMBRE DE NOEUDS'               ,T32,'=',I10/
     %' NOMBRE D''OBJETS INTERNES'      ,T32,'=',I10/
     %' NOMBRE D''OBJETS AUX LIMITES'   ,T32,'=',I10/
     %' NOMBRE DE TYPES D''EF'          ,T32,'=',I10)
      ELSE
         WRITE(IMPRIM,11210) NDIM,NBNOEU,NBOBIN,NBOBCL,NBTYEL
11210 FORMAT(' SPACE''s DIMENSION'      ,T33,'=',I10/
     %' NUMBER of NODES'                ,T33,'=',I10/
     %' NUMBER of MATERIALS'            ,T33,'=',I10/
     %' NUMBER of PL(S) at the BOUNDARY',T33,'=',I10/
     %' NUMBER of TYPES of FE'          ,T33,'=',I10)
      ENDIF
C
C     DECLARATION OUVERTURE DU FICHIER openfoam.points
C     ================================================
C     LE NOM DU FICHIER
      LNOMOB = NUDCNB( KNOMOB )
      KNMFIC = 'openfoam.points.' //  KNOMOB(1:LNOMOB)
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
C     ECRITURE DES 3 COORDONNEES DES SOMMETS SUR LE FICHIER openfoam.points
      CALL MEFIFOAM( NF )
C
      WRITE(NF,10020)
10020 FORMAT('FoamFile'/'{'/
     % T5,'version',  T16,'2.0;'/
     % T5,'format',   T16,'ascii;'/
     % T5,'class',    T16,'vectorField;'/
     % T5,'location', T16,'"constant\polyMesh";'/
     % T5,'object',   T16,'points;'/
     %'}'/)
C
C     ECRITURE DE XYZ DES NBNOEU SOMMETS
10011 FORMAT( '//',80('*'),'//' )
10012 FORMAT( '(' )
10013 FORMAT( ')' )
C
      WRITE(NF,10011)
      WRITE(NF,10021) NBNOEU
10021 FORMAT( I10 )
      WRITE(NF,10012)
10022 FORMAT( '(',3E16.7,')' )
C
C     ADRESSE DE LA 1-ERE COORDONNEE DU 1-ER SOMMET DU TMS 'XYZNOEUD'
      MNS = MNXYZN + WYZNOE
      DO I=1,NBNOEU
         WRITE(NF,10022) (RMCN(MNS+K),K=0,2)
         MNS = MNS + NBCOOR
      ENDDO
C
C     FIN DU FICHIER openfoam.points
      WRITE(NF,10013)
      WRITE(NF,10011)
C
C     FERMETURE DU FICHIER openfoam.points
      CLOSE( NF )
C
C
C     BOUCLE SUR LES OBJETS "INTERNES"
C     ================================
      WRITE(IMPRIM,*)
C     ECRITURE DU NOMBRE D'OBJETS INTERNES OU MATERIAUX
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10100) NBOBIN
      ELSE
         WRITE(IMPRIM,11100) NBOBIN
      ENDIF
10100 FORMAT( I8, ' NBOBIN NOMBRE DE MATERIAUX' )
11100 FORMAT( I8, ' NBOBIN NUMBER of MATERIALS' )
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
            WRITE(IMPRIM,10102) NUOB, KNM, KNOM
         ELSE
            WRITE(IMPRIM,11102) NUOB, KNM, KNOM
         ENDIF
 100  CONTINUE
10102 FORMAT( I8, ' est le NUMERO dans ',A10,' de ',A)
11102 FORMAT( I8, ' is the NUMBER in ',A10,' of ',A )
C
      IF( IERR .NE. 0 ) GOTO 9999
      NBPLSV(NDIMV+1) = NBOBIN
C
C     CREATION OU REDECOUVERTE DU TMS OBJET>>>FACE DES FACES DE L'OBJET
C     =================================================================
      CALL HACHOB( KNOMOB, 4, NTFAOB, MNFAOB, IERR )
      IF( IERR .NE. 0 ) GOTO 9999
C
C     EN SORTIE LE TABLEAU a___face  CONTIENT :
C     variable MOFACE 'nombre d''entiers stockes par face(>=8)' entier ;
C     variable MXFACE 'nombre maximal de faces'                 entier ;
C     variable L1FAFR 'numero de la 1-ere face frontaliere'     entier ;
C     variable NBFAFR 'nombre de faces frontalieres'            entier ;
C     variable L1FA2M 'numero de la 1-ere face interface'       entier ;
C     variable NBFA2M 'nombre de faces interfaces'              entier ;
C     variable NBFATG 'nombre de faces avec des tangentes'      entier ;
C     tableau  LFACES(1..MOFACE,1..MXFACE)
C     LFACES(1,I)= NO DU 1-ER  SOMMET DE LA FACE
C     LFACES(2,I)= NO DU 2-EME SOMMET > 1-ER  SOMMET
C     LFACES(3,I)= NO DU 3-EME SOMMET DE LA FACE
C     LFACES(4,I)= NO DU 4-EME SOMMET DE LA FACE OU 0 SI TRIANGLE
C                  0 SI TRIANGLE
C
C     LFACES(5,I)= CHAINAGE HACHAGE SUR FACE SUIVANTE
C
C     ATTENTION: POUR UNE FACE
C        SI SOMMET 2 < DERNIER SOMMET  => FACE   DIRECTE DANS LE CUBE
C                    >                 => FACE INDIRECTE
C        UNE FACE DIRECTE EST VUE DE L EXTERIEUR AU CUBE
C            SOUS LA FORME DIRECTE
C           (CUBE SIGNIFIE ICI UN EF 3D: tetra, penta, hexaedre)
C     LFACES(6,I)= NUMERO DU 1-ER  CUBE CONTENANT CETTE FACE
C                  >0 SI FACE   DIRECTE DANS CE CUBE
C                  <0 SI FACE INDIRECTE DANS CE CUBE
C     SI LA FACE APPARTIENT A 2 CUBES ALORS
C        LFACES(7,I)= NUMERO DU 2-EME CUBE CONTENANT CETTE FACE
C                     >0 SI FACE   DIRECTE DANS CE CUBE
C                     <0 SI FACE INDIRECTE DANS CE CUBE
C     SINON
C        LFACES(7,I)= NUMERO DE LA FACE FRONTALIERE SUIVANTE
C                     0 SI C'EST LA DERNIERE
C     L1FAFR = NUMERO DE LA PREMIERE FACE FRONTALIERE DANS LFACES
C     L1FA2M = NUMERO DE LA PREMIERE FACE INTERFACE   DANS LFACES
C
C     LFACES(8,I)= NUMERO DE LA FACE A TANGENTE DANS LE TABLEAU NUTGFA
C
C     LFACES(9,I)= NUMERO DU TYPE DU DERNIER EF DE CETTE FACE
C                  CETTE INFORMATION EXISTE SEULEMENT POUR UN OBJET
C                  PAS POUR UNE SURFACE OU UN VOLUME!
C
C     TABLEAU NOFACE( MXFACE ) NO DE LA FACE EXISTANTE OU ZERO SINON
      MOFACE = MCN( MNFAOB + WOFACE )
      MXFACE = MCN( MNFAOB + WXFACE )
      LIBREF = MXFACE
      CALL TNMCDC( 'ENTIER', MXFACE, MNFAEX )
C
C     CALCUL DU NOMBRE D'EF
C     ---------------------
      NBEF   = 0
      DO NOTYEL = 1, NBTYEL
C
C        L'ADRESSE DU TABLEAU NPEF"xxxx
         MNELE = MCN( MNELEM - 1 + NOTYEL )
C
C        LE NOMBRE D'ELEMENTS FINIS DE CE TYPE
         NBELEM = MCN( MNELE + WBELEM )
C
C        LE NOMBRE TOTAL D'EF
         NBEF = NBEF + NBELEM
C
      ENDDO
C
C     BOUCLE SUR LES NBFACES FACES DE L'OBJET POUR CONSTRUIRE FAEX
C     ------------------------------------------------------------
C     FAEX(NF) =  0 SI FACE NON UTILISEE DANS LFACES
C              = +NUMERO DE FACE UTILISEE SI APPARTIENT A 2EF
C              = -NUMERO DE FACE UTILISEE SI APPARTIENT A 1EF DONC FRONTIERE
C     ENSUITE LE NO NEGATIF DES FACES FRONTIERES EST REMPLACE PAR
C                -NUMERO DE SURFACE SI ELLE   APPARTIENT     A UNE SURFACE
C                - 999000           SI ELLE N'APPARTIENT PAS A UNE SURFACE
      MAFAEX = MNFAEX - 1
      NBFAEX = 0
      MNFACE = MNFAOB + WFACES
      MNF    = MNFACE - MOFACE - 1
C
C     COMPTAGE DU NOMBRE DE FACES EXISTANTES ET CONSTRUCTION DE
C     FAEX( K ) = NO DE FACE EXISTANTE
C                +NO SI APPARTIENT A 2 EF
C                -NO SI APPARTIENT A 1 EF
      DO K = 1, MXFACE
         MNF = MNF + MOFACE
         IF( MCN(MNF+1) .EQ. 0 ) THEN
C           FACE VIDE
            MCN( MAFAEX + K ) = 0
         ELSE
C           FACE EXISTANTE
            NBFAEX = NBFAEX + 1
            MCN( MAFAEX + K ) = NBFAEX
         ENDIF
      ENDDO
C
C     LES FACES FRONTALIERES SONT MARQUEES DANS FAEX
C     PAR - LE NUMERO DE FACE EXISTANTE DANS LFACES
      MNF = MNFACE - MOFACE - 1
      L1FAFR = MCN( MNFAOB + W1FAFR )
      NBFAFR = 0
 23   IF( L1FAFR .GT. 0 ) THEN
         NBFAFR = NBFAFR + 1
         MN = MNF + MOFACE * L1FAFR
C        LE NUMERO FAEX EST INVERSE
         MCN( MAFAEX + L1FAFR ) = -MCN( MAFAEX + L1FAFR )
C        FACE FRONTALIERE SUIVANTE
         L1FAFR = MCN( MN + 7 )
         GOTO 23
      ENDIF
C
      NBFA2F = NBFAEX - NBFAFR
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,*) 'NOMBRE DE FACES          =',NBFAEX
         WRITE(IMPRIM,*) 'NOMBRE DE FACES DANS 1 EF=',NBFAFR
         WRITE(IMPRIM,*) 'NOMBRE DE FACES DANS 2 EF=',NBFA2F
      ELSE
         WRITE(IMPRIM,*) 'NUMBER of FACES        =',NBFAEX
         WRITE(IMPRIM,*) 'NUMBER of FACES in 1 FE=',NBFAFR
         WRITE(IMPRIM,*) 'NUMBER of FACES in 2 FE=',NBFA2F
      ENDIF
C
C     TABLEAU DU NUMERO DES SURFACES DE L'OBJET
C     NUSF(I,1) = NO DE LA SURFACE DE LA I-EME SURFACE DE L'OBJET
C     NUSF(I,2) = NO DE LA PREMIERE FACE DE LA SURFACE APRES RENUMEROTATION
C          I = 1 A NBSURF
      CALL TNMCDC( 'ENTIER', 2*NBOBCL+2, MNNUSF )
C
C     RECENSEMENT DU NOMBRE DE POINTS, LIGNES, ... DE L'OBJET
C     =======================================================
      NBSURF = 0
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
            IF( NYOB .EQ. 3 ) THEN
C              STOCKAGE DE SON NUMERO DE SURFACE DANS LE LEXIQUE SURFACES
               MCN( MNNUSF + NBSURF ) = NUOB
C              UNE SURFACE DE PLUS
               NBSURF = NBSURF + 1
            ENDIF
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
C     BOUCLE SUR LES SURFACES(SI NDIMV=3) puis LIGNES puis POINTS
      DO 150 L=NDIMV,1,-1
C        ECRITURE DU NOMBRE D'OBJETS AUX LIMITES DE TYPE L
         KNM = NMTYOB( L )
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10140) NBPLSV(L), KNM
         ELSE
            WRITE(IMPRIM,11140) NBPLSV(L), KNM
         ENDIF
10140    FORMAT( I8, ' NOMBRE DE ',A10 )
11140    FORMAT( I8, ' NUMBER of ',A10 )
         MN = MNOBCL - 2
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
C
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
      NBFASF = 0
      NBFANF = 0
C
C     LA BOUCLE SUR LES DIFFERENTS TYPES D'EF (TABLEAUX NPEF"xxxx)
C     ============================================================
C     ECRITURE DU NOMBRE DE TYPE D'EF
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10400) NBTYEL
      ELSE
         WRITE(IMPRIM,11400) NBTYEL
      ENDIF
10400 FORMAT( I8, ' NBTYEL NOMBRE DE TYPES D''EF de l''OBJET' )
11400 FORMAT( I8, ' NBTYEL NUMBER of TYPES of FE of the OBJECT' )
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
         WRITE(IMPRIM,10401) NBELEM
10401    FORMAT( I8, ' NBELEM NOMBRE D''EF DE CE TYPE D''EF' )
         WRITE(IMPRIM,10402) NBNOE
10402    FORMAT( I8, ' NBNOE  NOMBRE DE NOEUDS DE CE TYPE D''EF')
         WRITE(IMPRIM,10403) NBSOE
10403    FORMAT( I8, ' NBSOE  NOMBRE DE SOMMETS DE CE TYPE D''EF')
         WRITE(IMPRIM,10404) NARET
10404    FORMAT( I8, ' NBARET NOMBRE D''ARETES DE CE TYPE D''EF')
         WRITE(IMPRIM,10405) NBFAC
10405    FORMAT( I8, ' NBFACE NOMBRE DE FACES DE CE TYPE D''EF')
         WRITE(IMPRIM,10406) NBVOL
10406    FORMAT( I8, ' NBVOL  NOMBRE DE VOLUME DE CE TYPE D''EF')
         ELSE
         WRITE(IMPRIM,11401) NBELEM
11401    FORMAT( I8, ' NBELEM NUMBER of FE of THIS TYPE' )
         WRITE(IMPRIM,11402) NBNOE
11402    FORMAT( I8, ' NBNOE  NUMBER of NODES of a FE of THIS TYPE')
         WRITE(IMPRIM,11403) NBSOE
11403    FORMAT( I8,
     %          ' NBSOE  NUMBER of VERTICES of a FE of THIS TYPE')
         WRITE(IMPRIM,11404) NARET
11404    FORMAT( I8, ' NBARET NUMBER of EDGES of a FE of THIS TYPE')
         WRITE(IMPRIM,11405) NBFAC
11405    FORMAT( I8, ' NBFACE NUMBER of FACES of a FE of THIS TYPE')
         WRITE(IMPRIM,11406) NBVOL
11406    FORMAT( I8,
     %   ' NBVOL  NUMBER of VOLUMES of a FE of THIS TYPE')
         ENDIF
C
         CALL AZEROI( 64, MCN(MNPS) )
         CALL AZEROI( 12, MCN(MNLA) )
         CALL AZEROI(  6, MCN(MNSF) )
         MCN( MNVC ) = 0
C
C        LES ADRESSAGES POUR RECUPERER LES INFORMATIONS
C        SELON LE TYPE NUTYEL DE L'ELEMENT FINI
         GOTO( 401, 401, 401, 401, 400, 400, 400, 400, 400, 400,
     &         400, 400, 401, 400, 401, 401, 400, 401, 401, 401,
     &         401, 401, 401, 401, 400, 400, 400, 401, 401, 401,
     %         401, 401, 401, 401                         ),NUTYEL
C
C        ERREUR
 400     NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'ELEMENT FINI '// NOMELE(1)
     &             // NOMELE(2) //' NON PROGRAMME dans openfoam.f'
         ELSE
            KERR(1) = 'FINITE ELEMENT '// NOMELE(1)
     &             // NOMELE(2) //' NOT PROGRAMMED in openfoam.f'
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
C           LE NUMERO DE VOLUME  DE L'EF
C           LE NUMERO DE SURFACE DES FACES   DE L'EF
C           LE NUMERO DE LIGNE   DES ARETES  DE L'EF
C           LE NUMERO DE POINT   DES SOMMETS DE L'EF
C           ----------------------------------------
            CALL EFPLSV( MNELE , NUELEM ,
     %                   NOVCEL, NOSFEL, NOLAEL, NOPSEL,
     %                   MCN(MNVC), MCN(MNSF), MCN(MNLA), MCN(MNPS) )
C
cccC           ECRITURE:
cccC           DES NUMEROS DES NOEUDS
cccC           DES POINTS DES SOMMETS
cccC           DES LIGNES DES ARETES
cccC           DES SURFACES DES FACES
cccC           DU VOLUME
ccc            WRITE(IMPRIM,10500) (MCN(MNAUX1+L),L=1,NBNOE),
ccc     %                          (MCN(MNPS1+L) ,L=1,NBSOE),
ccc     %                          (MCN(MNLA1+L) ,L=1,NARET),
ccc     %                          (MCN(MNSF1+L) ,L=1,NBFAC),
ccc     %                          (MCN(MNVC1+L) ,L=1,NBVOL)
C
C           PARCOURS DES FACES DE L'EF POUR STOCKER DANS LFACES
C           LE NO DE SURFACE DES FACES FRONTALIERES
C           ===================================================
            DO NFA=1,NBFAC
C
C              NOSF = NO DE LA SURFACE DE LA FACE NFA
C                     0 SI LA FACE N'APPARTIENT PAS A UNE SURFACE
               NOSF = MCN(MNSF1+NFA)
C
C              LA FACE APPARTIENT A LA SURFACE NOSF
C              NBSTF LE NOMBRE DE SES SOMMETS
               NBSTF = NBSOFA(NFA)
C
C              RECHERCHE DE CETTE FACE DANS LE TABLEAU DES FACES
               ISENS = 0
               NOMIN = NBNOEU
               DO J=1,NBSTF
C                 LE NUMERO J DE SOMMET DE LA FACE NFA DE L'EF
                  NGS(J) = MCN( MNAUX1 + NOSOFA(J,NFA) )
                  IF( NGS(J) .LT. NOMIN ) THEN
                     NOMIN = NGS(J)
                     ISENS = J
                  ENDIF
               ENDDO
C
C              PERMUTATION CIRCULAIRE DES SOMMETS POUR AMENER NGS(ISENS)
C              LE PLUS PETIT NO EN PREMIERE POSITION DE NGS1
               DO J=1,NBSTF
                  NGS1(J) = NGS(ISENS)
                  ISENS   = ISENS + 1
                  IF( ISENS .GT. NBSTF ) ISENS = 1
               ENDDO
C
C              LA NUMEROTATION DES SOMMETS EST MISE A JOUR POUR
C              QUE LA 1-ERE ARETE DE LA FACE SOIT
C              1:NO LE PLUS FAIBLE DES SOMMETS DE LA FACE
C              2:NO MIN(SOMMET2, SOMMET NBSTF)
C              SI MIN = SOMMET 2  FACE DIRECTE
C                       SOMMET NBSTF FACE INDIRECTE
C              UNE FACE DIRECTE EST VUE DE L EXTERIEUR AU SOLIDE
C              SOUS LA FORME DIRECTE
               J = NGS1( 2 )
               IF( J .GT. NGS1(NBSTF) ) THEN
C                 FACE INDIRECTE
                  ISENS   = -1
                  NGS1(2) = NGS1(NBSTF)
                  NGS1(NBSTF) = J
               ELSE
                  ISENS = 1
               ENDIF
C
C              RECHERCHE DE LA FACE
C              --------------------
               CALL HACHAG( NBSTF, NGS1, MOFACE, MXFACE, MCN(MNFACE),
     %                      5, LIBREF,   NOFAC )
               IF( NOFAC .LE. 0 ) THEN
                  print *,'ERREUR dans openfoam: FACE NON RETROUVEE'
                  print *,'NO TYPE ELEMENT=',NOTYEL,' NO EF=',NUELEM,
     %                    '  SOMMETS=',(NGS1(J),J=1,NBSTF)
C                 ADRESSE-NBCOOR-1 DE LA 1-ERE COORDONNEE DU
C                 1-ER SOMMET DU TMS 'XYZNOEUD'
                  MNS = MNXYZN + WYZSOM - NBCOOR -1
                  DO J=1,NBSTF
                     MN = MNS + NGS1(J) * NBCOOR
                     print *,'SOMMET',NGS1(J),(RMCN(MN+K),K=1,NBCOOR)
                  ENDDO
                  GOTO 9999
               ENDIF
C              NOFAC>0 NUMERO DE LA FACE DANS LE TABLEAU LFACES
               MN = MNFACE + MOFACE * ( NOFAC - 1 ) -1
C
               IF( MCN(MAFAEX+NOFAC) .LT. 0 ) THEN
C
C                 LA FACE EST FRONTALIERE
C                 ATTENTION:
C                 ====================================================
C                 LE NO (1 A NBTYEL) DU TYPE DE L'EF DANS LE MAILLAGE
                  NTYEL = MCN( MN + 9 )
                  IF( NTYEL .GT. 1 000 000 ) NTYEL = NTYEL / 1 000 000
C                 EST TRANSFORME EN 1000000 * NTYEL + NOSF !
C                 ====================================================
                  MCN( MN + 9 ) = 1 000 000 * NTYEL + NOSF
C
C                 LE NUMERO FAEX( FACE ) = - NO de SURFACE SI NO>0
C                                          - 999 000 SI FACE FRONTALIERE
C                                                NON SUR UNE SURFACE
                  IF( NOSF .GT. 0 ) THEN
C
C                    FACE DE LA SURFACE NOSF
                     MCN( MAFAEX + NOFAC ) = -NOSF
C                    UNE FACE DES SURFACES DE L'OBJET DE PLUS
                     NBFASF = NBFASF + 1
C
                  ELSE IF( NOSF .EQ. 0 ) THEN
C
C                    FACE FRONTALIERE SUR AUCUNE SURFACE
                     MCN( MAFAEX + NOFAC ) = - 999 000
                     NBFANF = NBFANF + 1
C
                  ENDIF
C
               ENDIF
C
            ENDDO
 450     CONTINUE
 500  CONTINUE
CCC10500 FORMAT( 10I8 )
C
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,*)NBFASF,' FACES des SURFACES        dans LFACES'
         WRITE(IMPRIM,*)NBFANF,' FACES NON sur une SURFACE dans LFACES'
      ELSE
         WRITE(IMPRIM,*)NBFASF,' FACES of SURFACES      in LFACES'
         WRITE(IMPRIM,*)NBFANF,' FACES NOT on a SURFACE in LFACES'
      ENDIF
C
C     DESTRUCTION DU TABLEAU AUXILIAIRE
      CALL TNMCDS( 'ENTIER', 64+64+12+6+1+24, MNAUX )
C
C
C     DECLARATION OUVERTURE DU FICHIER openfoam.faces
C     ===============================================
C     LE NOM DU FICHIER
      KNMFIC = 'openfoam.faces.' //  KNOMOB(1:LNOMOB)
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
C     ECRITURE DE L'ENTETE DU FICHIER openfoam.faces
      CALL MEFIFOAM( NF )
      WRITE(NF,10031)
10031 FORMAT('FoamFile'/'{'/
     % T5,'version',  T16,'2.0;'/
     % T5,'format',   T16,'ascii;'/
     % T5,'class',    T16,'faceList;'/
     % T5,'location', T16,'"constant\polyMesh";'/
     % T5,'object',   T16,'faces;'/
     %'}'/)
C
C     LE TABLEAU FARE(N)=M N DE 1 A NBFAEX M DE 1 A MXFACE
C                        NO DANS LFACES DE LA N-EME FACE RENUMEROTEE
      CALL TNMCDC( 'ENTIER', NBFAEX, MNFARE )
      MAFARE = MNFARE - 1
C
C     BOUCLE SUR LES NBFACES FACES DE L'OBJET
C     ---------------------------------------
10030 FORMAT(I1,'(',3I10,')')
10040 FORMAT(I1,'(',4I10,')')
C
C     ICI FAEX(NF) CONTIENT POUR CHAQUE FACE NF DU TABLEAU LFACES:
C     FAEX(NF) =  0 SI LA FACE NF EST NON UTILISEE DANS LFACES
C              = +NUMERO DE FACE UTILISEE SI APPARTIENT A 2EF
C                -NUMERO DE SURFACE SI ELLE   APPARTIENT     A UNE SURFACE
C                - 999000           SI ELLE N'APPARTIENT PAS A UNE SURFACE
C
C     RENUMEROTATION DES FACES SELON L'ORDRE:
C     1. FAEX(NF)>0 LES FACES EXISTANTES INTERNES APPARTENANT A 2 EF
C     2. FAEX(NF)=-( NOSF1 )
C     1+NBSURF:   FAEX(NF)=-( NOSF 2...NBSURF )
C     1+NBSURF+1: FAEX(NF)=- 999 000
C
C     ECRITURE  //***//    NBFACES    (
      WRITE(NF,10011)
      WRITE(NF,10021) NBFAEX
      WRITE(NF,10012)
C
C     RENUMEROTATION DES FACES EXISTANTES APPARTENANT A 2 EF
C     ET ECRITURE DU NO DES SOMMETS DES FACES EXISTANTES
C     ------------------------------------------------------
      MNFAC1  = MNFACE - MOFACE - 1
      NBFACER = 0
      DO NFA = 1, MXFACE
C
C        NO DE FACE EXISTANTE
         NOFAE = MCN( MAFAEX + NFA )
C
         IF( NOFAE .GT. 0 ) THEN
C
C           UNE FACE EXISTANTE APPARTENANT A 2 EF DE PLUS
            NBFACER = NBFACER + 1
C
C           ANCIEN  NO => NOUVEAU NO
            MCN( MAFAEX + NFA ) = NBFACER
C           NOUVEAU NO => ANCIEN  NO
            MCN( MAFARE + NBFACER ) = NFA
C
C           ADRESSE DE CETTE FACE DANS LFACES
            MNF = MNFAC1 + MOFACE * NFA
C
C           ECRITURE SUR LE FICHIER SELON NOMBRE DE SOMMETS DE LA FACE
            IF( MCN(MNF+4) .EQ. 0 ) THEN
C
C              TRIANGLE
               IF( MCN(MNF+6) .GT. 0 ) THEN
C                 TRIANGLE DIRECT
                  WRITE(NF,10030) 3, (MCN(MNF+KK)-1,KK=1,3)
               ELSE
C                 TRIANGLE INDIRECT => ECRIT DANS LE SENS DIRECT
                  WRITE(NF,10030) 3, MCN(MNF+1)-1,MCN(MNF+3)-1,
     %                               MCN(MNF+2)-1
               ENDIF
C
            ELSE
C
C              QUADRANGLE
               IF( MCN(MNF+6) .GT. 0 ) THEN
C                 QUADRANGLE DIRECT
                  WRITE(NF,10040) 4, (MCN(MNF+KK)-1,KK=1,4)
               ELSE
C                 TRIANGLE INDIRECT => ECRIT DANS LE SENS DIRECT
                  WRITE(NF,10040) 4, MCN(MNF+1)-1,MCN(MNF+4)-1,
     %                               MCN(MNF+3)-1,MCN(MNF+2)-1
               ENDIF
C
            ENDIF
C
         ENDIF
C
      ENDDO
      NBFA2F = NBFACER
C
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,*)'NOMBRE DE FACES DANS 2 EF            =',NBFA2F
      ELSE
         WRITE(IMPRIM,*)'NUMBER of FACES in 2 FE            =',NBFA2F
      ENDIF
C
C     RENUMEROTATION DES FACES EXISTANTES APPARTENANT A UNE SURFACE
C     -------------------------------------------------------------
C     NBTSFASF NOMBRE TOTAL DE FACES SUR LES SURFACES DE L'OBJET
      NBFASF = NBFACER
      DO NSURF = 1, NBSURF
C
C       -NO DE LA SURFACE NSURF
         MN = MNNUSF -1 + NSURF
         NOSURF = -MCN( MN )
C
C        NO DE LA 1-ERE FACE DE LA SURFACE NSURF
         MCN( MN + NBSURF + 1 ) = NBFACER+1
C
C        NO DE LA FACE QUI PRECEDE LA 1-ERE DE LA SURFACE
         NOFSM1 = NBFACER
C
         DO NFA = 1, MXFACE
C           NO DE FACE EXISTANTE
            NOFAE = MCN( MAFAEX + NFA )
            IF( NOFAE .EQ. NOSURF ) THEN
C
C              UNE FACE EXISTANTE APPARTENANT A 1 SURFACE DE L'OBJET
               NBFACER = NBFACER + 1
C
C              ANCIEN  NO => NOUVEAU NO
               MCN( MAFAEX + NFA ) = NBFACER
C              NOUVEAU NO => ANCIEN  NO
               MCN( MAFARE + NBFACER ) = NFA
C
C              ADRESSE DE CETTE FACE DANS LFACES
               MNF = MNFAC1 + MOFACE * NFA
C
C              ECRITURE SUR LE FICHIER SELON NOMBRE DE SOMMETS DE LA FACE
               IF( MCN(MNF+4) .EQ. 0 ) THEN
C
C                 TRIANGLE
                  IF( MCN(MNF+6) .GT. 0 ) THEN
C                    TRIANGLE DIRECT
                     WRITE(NF,10030) 3, (MCN(MNF+KK)-1,KK=1,3)
                  ELSE
C                    TRIANGLE INDIRECT => ECRIT DANS LE SENS DIRECT
                     WRITE(NF,10030) 3, MCN(MNF+1)-1,MCN(MNF+3)-1,
     %                                  MCN(MNF+2)-1
                  ENDIF
C
               ELSE
C
C                 QUADRANGLE
                  IF( MCN(MNF+6) .GT. 0 ) THEN
C                    QUADRANGLE DIRECT
                     WRITE(NF,10040) 4, (MCN(MNF+KK)-1,KK=1,4)
                  ELSE
C                    TRIANGLE INDIRECT => ECRIT DANS LE SENS DIRECT
                     WRITE(NF,10040) 4, MCN(MNF+1)-1,MCN(MNF+4)-1,
     %                                  MCN(MNF+3)-1,MCN(MNF+2)-1
                  ENDIF
C
               ENDIF
C
            ENDIF
C
         ENDDO
C
         CALL NMOBNU( 'SURFACE', ABS(NOSURF), KNOMSF )
         NB = NBFACER-NOFSM1
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,*) NB, ' FACES dans la SURFACE ', KNOMSF
         ELSE
            WRITE(IMPRIM,*) NB, ' FACES in the SURFACE ', KNOMSF
         ENDIF
C
      ENDDO
      NBFASF = NBFACER - NBFASF
C
C     RENUMEROTATION DES FACES EXISTANTES APPARTENANT A 1 EF
C     ET N'APPARTENANT PAS AUX FACES DES SURFACES DE L'OBJET
C     ELLES SONT DITES APPARTENIR A LA SURFACE NBSURF+1
C     ------------------------------------------------------
      NBFANF = NBFACER
      NSURF  = NBSURF + 1
C
C     NO DE LA 1-ERE FACE DE LA SURFACE NBSURF+1
      MCN( MNNUSF -1 + 2*NBSURF + 2 ) = NBFACER+1
C
      DO NFA = 1, MXFACE
C        NO DE FACE EXISTANTE
         NOFAE = MCN( MAFAEX + NFA )
         IF( NOFAE .EQ. -999 000 ) THEN
C
C           UNE FACE EXISTANTE APPARTENANT A 1 EF ET PAS DE SURFACE
            NBFACER = NBFACER + 1
C
C           ANCIEN  NO => NOUVEAU NO
            MCN( MAFAEX + NFA ) = NBFACER
C           NOUVEAU NO => ANCIEN  NO
            MCN( MAFARE + NBFACER ) = NFA
C
C           ADRESSE DE CETTE FACE DANS LFACES
            MNF = MNFAC1 + MOFACE * NFA
C
C           ECRITURE SUR LE FICHIER SELON NOMBRE DE SOMMETS DE LA FACE
            IF( MCN(MNF+4) .EQ. 0 ) THEN
C
C              TRIANGLE
               IF( MCN(MNF+6) .GT. 0 ) THEN
C                 TRIANGLE DIRECT
                  WRITE(NF,10030) 3, (MCN(MNF+KK)-1,KK=1,3)
               ELSE
C                 TRIANGLE INDIRECT => ECRIT DANS LE SENS DIRECT
                  WRITE(NF,10030) 3, MCN(MNF+1)-1,MCN(MNF+3)-1,
     %                               MCN(MNF+2)-1
               ENDIF
C
            ELSE
C
C              QUADRANGLE
               IF( MCN(MNF+6) .GT. 0 ) THEN
C                 QUADRANGLE DIRECT
                  WRITE(NF,10040) 4, (MCN(MNF+KK)-1,KK=1,4)
               ELSE
C                 TRIANGLE INDIRECT => ECRIT DANS LE SENS DIRECT
                  WRITE(NF,10040) 4, MCN(MNF+1)-1,MCN(MNF+4)-1,
     %                               MCN(MNF+3)-1,MCN(MNF+2)-1
               ENDIF
C
            ENDIF
         ENDIF
      ENDDO
      NBFANF = NBFACER - NBFANF
C
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,*)'NOMBRE DE FACES DANS 1 EF & 0 SURFACE=',NBFANF
         WRITE(IMPRIM,*)'NOMBRE TOTAL DE FACES RENUMEROTEES   =',NBFACER
      ELSE
         WRITE(IMPRIM,*)'NUMBER of FACES in 1 FE & 0 SURFACE=',NBFANF
         WRITE(IMPRIM,*)'TOTAL NUMBER of RENUMBERED FACES   =',NBFACER
      ENDIF
C
C     FIN DU FICHIER openfoam.faces
      WRITE(NF,10013)
      WRITE(NF,10011)
C
C     FERMETURE DU FICHIER openfoam.faces
      CLOSE( NF )
      WRITE(IMPRIM,*) 'openfoam.faces is closed'
C
C
C     DECLARATION OUVERTURE DU FICHIER openfoam.owner
C     NO D'UN EF CONTENANT LA FACE
C     ===============================================
C     LE NOM DU FICHIER
      KNMFIC = 'openfoam.owner.' //  KNOMOB(1:LNOMOB)
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
C     ECRITURE DU PREMIER NUMERO D'EF DE CHAQUE FACE
C     SUR LE FICHIER openfoam.owner
      CALL MEFIFOAM( NF )
      WRITE(NF,10032) NBNOEU, NBEF, NBFAFR, NBFA2F
10032 FORMAT('FoamFile'/'{'/
     % T5,'version',  T16,'2.0;'/
     % T5,'format',   T16,'ascii;'/
     % T5,'class',    T16,'labelList;'/
     % T5,'note',     T16,'"nPoints:',I10,' nCells:',I10,' nFaces:',
     %I10,' nInternalFaces:',I10,'";'/
     % T5,'location', T16,'"constant\polyMesh";'/
     % T5,'object',   T16,'owner;'/
     %'}'/)
C
C     ECRITURE DE XYZ DES NBFAEX NUMERO DU 1-ER EF CONTENANT LA FACE
      WRITE(NF,10011)
      WRITE(NF,10021) NBFAEX
      WRITE(NF,10012)
C
C     ECRITURE DU NO DU PREMIER EF DE CHAQUE FACE
      MNFAC1 = MNFACE - MOFACE - 1
      DO K = 1, NBFAEX
C
C        NO DE LA FACE DANS LFACES
         NLFAC = MCN( MAFARE + K )
C
C        ADRESSE DANS LFACES
         MNF = MNFAC1 + MOFACE * NLFAC
C
C        LA NUMEROTATION DEBUTE A ZERO => -1
         WRITE(NF,'(I10)') ABS(MCN(MNF+6))-1
C
      ENDDO
C
C     FIN DU FICHIER openfoam.owner
      WRITE(NF,10013)
      WRITE(NF,10011)
C
C     FERMETURE DU FICHIER openfoam.owner
      CLOSE( NF )
      WRITE(IMPRIM,*) 'openfoam.owner is closed'
C
C
C     DECLARATION OUVERTURE DU FICHIER openfoam.neighbour
C     NO DU SECOND EF DE LA FACE
C     ===================================================
C     LE NOM DU FICHIER
      KNMFIC = 'openfoam.neighbour.' //  KNOMOB(1:LNOMOB)
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
C     ECRITURE DU PREMIER NUMERO D'EF DE CHAQUE FACE
C     SUR LE FICHIER openfoam.neighbour
      CALL MEFIFOAM( NF )
      WRITE(NF,10033) NBNOEU, NBEF, NBFAFR, NBFA2F
10033 FORMAT('FoamFile'/'{'/
     % T5,'version',  T16,'2.0;'/
     % T5,'format',   T16,'ascii;'/
     % T5,'class',    T16,'labelList;'/
     % T5,'note',     T16,'"nPoints:',I10,' nCells:',I10,' nFaces:',
     %I10,' nInternalFaces:',I10,'";'/
     % T5,'location', T16,'"constant\polyMesh";'/
     % T5,'object',   T16,'neighbour;'/
     %'}'/)
C
C     ECRITURE DE XYZ DES NBFAEX NUMERO DU 1-ER EF CONTENANT LA FACE
      WRITE(NF,10011)
      WRITE(NF,10021) NBFAEX
      WRITE(NF,10012)
C
C     ECRITURE DU NO DU SECOND EF DE CHAQUE FACE APPARTENANT A 2 EF
      DO K = 1, NBFA2F
C
C        NO DE LA FACE DANS LFACES
         NLFAC = MCN( MAFARE + K )
C
C        ADRESSE DANS LFACES
         MNF = MNFAC1 + MOFACE * NLFAC
C
C        FACE APPARTENANT A 2 EF
C        LA NUMEROTATION DEBUTE A ZERO => -1
         WRITE(NF,'(I10)') ABS(MCN(MNF+7))-1
      ENDDO
C
C     ECRITURE DU NO (-1) DU SECOND EF DE CHAQUE FACE APPARTENANT A 1 EF
      DO K = NBFA2F+1, NBFAEX
C
C        NO DE LA FACE DANS LFACES
         NLFAC = MCN( MAFARE + K )
C
C        ADRESSE DANS LFACES
         MNF = MNFAC1 + MOFACE * NLFAC
C
C        FACE APPARTENANT 1 SEUL EF => PAR CONVENTION -1
         WRITE(NF,'(I10)') -1
C
      ENDDO
C
C     FIN DU FICHIER openfoam.neighbour
      WRITE(NF,10013)
      WRITE(NF,10011)
C
C     FERMETURE DU FICHIER openfoam.neighbour
      CLOSE( NF )
      WRITE(IMPRIM,*) 'openfoam.neighbour is closed'
C
C     DECLARATION OUVERTURE DU FICHIER openfoam.boundary
C     POUR OPENFOAM SEULES LES SURFACES SUPPORTENT UNE CONDITION AUX LIMITES
C     1 SURFACE MEFISTO => 1 BOUNDARY POUR OpenFOAM
C     ======================================================================
C     LE NOM DU FICHIER
      KNMFIC = 'openfoam.boundary.' //  KNOMOB(1:LNOMOB)
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
C     ECRITURE DU PREMIER NUMERO D'EF DE CHAQUE FACE
C     SUR LE FICHIER openfoam.boundary
      CALL MEFIFOAM( NF )
      WRITE(NF,10035)
10035 FORMAT('FoamFile'/'{'/
     % T5,'version',  T16,'2.0;'/
     % T5,'format',   T16,'ascii;'/
     % T5,'class',    T16,'polyBoundaryMesh;'/
     % T5,'location', T16,'"constant\polyMesh";'/
     % T5,'object',   T16,'boundary;'/
     %'}'/)
C
C     ECRITURE DE XYZ DES NBFAEX NUMERO DU 1-ER EF CONTENANT LA FACE
C     1 SURFACE MEFISTO => 1 BOUNDARY POUR OpenFOAM
C     NBOUNDARY NOMBRE DE SURFACES DE L'OBJET
      NBOUNDARY = NBSURF+1
      WRITE(NF,10011)
      WRITE(NF,10021) NBOUNDARY
      WRITE(NF,10012)
C
10051 FORMAT(T5,A/T5,'{'/
     %   T10,'type',     T22,'patch;')
C
10052 FORMAT(T10,'nFaces',   T22,I10,';'/
     %       T10,'startFace',T22,I10,'; // ',A,'//')
C
      DO NSURF = 1, NBSURF
C
C        NOM DE SURFACE = NOM DE BOUNDARY
C        ---------------------------------------------
C        NOSU NUMERO DE LA SURFACE DANS LE LEXIQUE DES SURFACES
         MNF  = MNNUSF - 1 + NSURF
         NOSU = MCN( MNF )
C
C        OUVERTURE DE LA SURFACE
         CALL LXNLOU( NTMN(3), NOSU, NTOB, MNOB )
C
C        NOM DE LA SURFACE
         CALL NMOBNU( 'SURFACE', NOSU, KNOM )
C
C        ECRITURE DU NOM DE LA SURFACE DANS openfoam.boundary
         L = NUDCNB( KNOM )
         WRITE(NF,10051) KNOM(1:L)
         WRITE(IMPRIM,*)
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,*)
     %     'SURFACE ',KNOM(1:L),' est ajoute a openfoam.boundary'
         ELSE
            WRITE(IMPRIM,*) 'SURFACE ',KNOM(1:L),
     %      ' put in openfoam.boundary'
         ENDIF
C
C        NUMERO DE LA PREMIERE FACE DE CETTE SURFACE DANS LE FICHIER faces
         N1F = MCN( MNF + NBSURF + 1 )
C        NUMERO+1 DE LA DERNIERE FACE DE CETTE SURFACE DANS LE FICHIER faces
         NDF = MCN( MNF + NBSURF + 2 )
C
C        ECRITURE DU NOMBRE DE FACES DE LA SURFACE ET
c        NUMERO DE LA PREMIERE FACE
         WRITE(NF,10052) NDF-N1F, N1F-1,KNOM
C
C        FIN DE LA SURFACE=BOUNDARY
         WRITE(NF,'(T5,''}'')' )
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,*) NDF-N1F,' FACES DE LA SURFACE ',KNOM
         ELSE
            WRITE(IMPRIM,*) NDF-N1F,' FACES of SURFACE ',KNOM
         ENDIF
C
      ENDDO
      IF( LANGAG .EQ. 0 ) THEN
        WRITE(IMPRIM,*) NBFASF,' FACES DES SURFACES DE L''OBJET ',KNOMOB
      ELSE
       WRITE(IMPRIM,*) NBFASF,' FACES of SURFACES of the OBJECT ',KNOMOB
      ENDIF
C
C
C     LE RESTANT DES FACES FRONTALIERES NON SUR LES SURFACES DE L'OBJET
C     -----------------------------------------------------------------
10061 FORMAT(T5,A/T5,'{'/
     %   T10,'type',     T22,'empty;')
      WRITE(NF,10061) 'BoundaryEmpty'
      WRITE(IMPRIM,*)
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,*) 'BoundaryEmpty est ajoute a openfoam.boundary'
      ELSE
         WRITE(IMPRIM,*) 'BoundaryEmpty put in openfoam.boundary'
      ENDIF
C
C     BOUCLE SUR LES FACES FRONTALIERES NON SUR UNE SURFACE DE L'OBJET
      KNOM = 'BOUNDARY FACES NOT on a SURFACE'
C
C     NUMERO DE LA PREMIERE FACE DE CETTE SURFACE DANS LE FICHIER faces
      N1F = MCN( MNNUSF -1 + 2*NBSURF + 2 )
C
C     ECRITURE DU NOMBRE DE FACES D'UN SEUL EF ET NON SUR UNE SURFACE
c     ET NUMERO DE LA PREMIERE FACE
      WRITE(NF,10052) NBFANF, N1F-1, KNOM
C
C     FIN DE LA SURFACE=BOUNDARY
      WRITE(NF,'(T5,''}'')' )
      IF( LANGAG .EQ. 0 ) THEN
       WRITE(IMPRIM,*) NBFANF,' FACES NON sur une SURFACE de l''OBJET ',
     %                 KNOMOB
      ELSE
         WRITE(IMPRIM,*) NBFANF,' FACES NOT on a SURFACE of OBJECT ',
     %                   KNOMOB
      ENDIF
C
C     FIN DU FICHIER openfoam.boundary
      WRITE(NF,10013)
      WRITE(NF,10011)
C
C     FERMETURE DU FICHIER openfoam.boundary
      CLOSE( NF )
      WRITE(IMPRIM,*) 'openfoam.boundary is closed'
      GOTO 9999
C
C     PB A L'OUVERTURE D'UN FICHIER openfoam.*
 9998 NBLGRC(NRERR) = 1
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = 'IMPOSSIBLE OUVERTURE d''un FICHIER openfoam'
      ELSE
         KERR(1) = 'A FILE openfoam CAN NOT BE OPENED'
      ENDIF
      CALL LEREUR
C
C     SUPPRESSION DES TABLEAUX INUTILES
 9999 IF( MNFAEX .GT. 0 ) CALL TNMCDS( 'ENTIER', MXFACE,     MNFAEX )
      IF( MNNUSF .GT. 0 ) CALL TNMCDS( 'ENTIER', 2*NBOBCL+2, MNNUSF )
      IF( MNFARE .GT. 0 ) CALL TNMCDS( 'ENTIER', NBFAEX,     MNFARE )
      IF( MNELEM .GT. 0 ) CALL TNMCDS( 'ENTIER', 2*MXTYEL, MNELEM )
C
      RETURN
      END
