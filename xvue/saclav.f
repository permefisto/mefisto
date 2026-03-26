      SUBROUTINE SACLAV
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : SAISIR INTERACTIVEMENT UN TEXTE POUR KLG a l'AIDE de
C ----- la SOURIS ou du CLAVIER
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS         MAI 1994
C23456---------------------------------------------------------------012
      PARAMETER     (DIST0=1E6)
      include"./incl/gsmenu.inc"
      include"./incl/lu.inc"
      include"./incl/pilect.inc"
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      include"./incl/ntmnlt.inc"
CCC      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      include"./incl/pp.inc"
      COMMON         MCN(MOTMCN)

      COMMON / LECLXQ / LECLEX
C     SERT AU DIALOGUE AVEC LE SP LIRLEX
C     LECLEX = 0     PAS DE NOM DANS UN LEXIQUE DEMANDE
C              1 a 5 NOM DANS le LEXIQUE POINT LIGNE SURFACE VOLUME OBJET
C              6 a 9 NUMERO de SOMMET ou ARETE ou FACE ou EF3D
C              1000  * NO DU TYPE DE L'OBJET SINON

      CHARACTER*1   CHAR,CARLU
      CHARACTER*10  NMTOBJ,NMTYOB
      CHARACTER*24  NMOBJ
      INTRINSIC     INT

      NSORTI = 0
      NBLANC = 0

 10   IF( LECLEX .GT. 0 .AND. LECLEX .LE. 5 ) THEN

C        CE QUI SUIT EST UNE VARIANTE DE CALL SAIPTC(NOCODE,NX,NY,NASCII)
C        AVEC LA GESTION DES ACCROCHAGES SUR LES ITEMS DU LEXIQUE LECLEX
C        EN ATTENTE d'UN NOM dans un des LEXIQUES
C        1:POINT,  2:LIGNE, 3:SURFACE, 4:VOLUME, 5:OBJET
C     ou d'UN NUMERO de
C        6:SOMMET, 7:ARETE, 8:FACE,    9:EF3D  dans UN MAILLAGE TRACE
C          attention: SANS LEXIQUE PLSVO!
C
C        SAISIE D'UN POINT PAR CLIC DE LA SOURIS DANS LE MENU
C        OU DANS LA FENETRE GRAPHIQUE OU PAR FRAPPE SUR LE CLAVIER PHYSIQUE
C        ATTENDRE LA FRAPPE D'UNE TOUCHE DE LA SOURIS OU DU CLAVIER
C        ==================================================================
C        ADRESSE DE DEBUT DES ITEMS DE CES OBJETS PRESENTS SUR L'ECRAN
         MNIT = MNITEM( LECLEX )
         IF( MNIT .LE. 0 ) THEN
            LECLEX = 0
            GOTO 10
         ENDIF
C        PAS D'ITEM ACCROCHE
         NMIN0  = -2
         NMIN00 = NMIN0
C
C        ATTENTE D'UN CLIC EN MONTRANT UN ACCROCHAGE SUR LES POIGNEES
C        DES ITEMS DU LEXIQUE
C        ------------------------------------------------------------
 11      CALL XVSOURIS2( MCN( MNIT ), NMIN0,
     %                   NOCODE, NASCII, NX, NY )
C
C        RETOURNE  NOCODE=0 SI PAS D'EVENEMENT
C                        =1 SI CLIC D'UN BOUTON DE LA SOURIS
C                           => NASCII=NUMERO DU BOUTON
C                        =2 SI FRAPPE D'UN CARACTERE AU CLAVIER
C                           => NASCII=NUMERO DU CARACTERE DANS LA TABLE ASCII
C                        =5 SI UN BOUTON EST ENFONCE, NON RELACHE
C                           L'ITEM LE PLUS PROCHE DU CLIC A ETE CALCULE
C                           NMIN0 CONTIENT LE DECALAGE DANS MCN(MNIT)
C                           DU (X,Y) PIXEL DE LA POIGNEE DE L'ITEM

         IF( NOCODE .EQ. 0 ) THEN

            PRINT *,'SACLAV: PROBLEME AVEC LA SOURIS'
            CALL XVPAUSE
            GOTO 11

         ELSE IF( NOCODE .EQ. 1 ) THEN

C           UN DES 3 BOUTONS DE LA SOURIS
            IF( NASCII .EQ. 2 ) THEN

C              BOUTON 2 DE LA SOURIS <=> FRAPPE AU CLAVIER DU CARACTERE '@'
               NOCODE = -1
               NASCII = 64

            ELSE

C              BOUTON 1 OU 3 => PAS DE MODIFICATION
               NOCODE = NASCII

            ENDIF

CCC         PRINT *,'SACLAV: EVENEMENT SOURIS : BOUTON ',NOCODE,
CCC     %           ' POSITION DU CLIC X=',NX,' Y=',NY

         ELSE IF( NOCODE .EQ. 2 ) THEN

C           FRAPPE D'UN CARACTERE AU CLAVIER
            NOCODE = -1

CCC         PRINT *,'SACLAV: EVENEMENT CLAVIER : NO ASCII CARACTERE ',
CCC     %            NASCII,' => CARACTERE=',CHAR(NASCII)

            IF( NASCII .EQ. 27 ) THEN

C              LE CARACTERE 'ECHAPPEMENT' DEVIENT LE CARACTERE '@'
               NASCII = 64

            ENDIF

         ELSE IF( NOCODE .EQ. 5 ) THEN

C           UN ITEM PLSVO A ETE ACCROCHE
            IF( NMIN00 .NE. NMIN0 .AND. NMIN0 .GE. 0 ) THEN
C              NUMERO DU PLSVO ou NUMERO de SOMMET ou ARETE ou ...
               NUOBJ = MCN(MNIT+NMIN0+2)

               IF( LECLEX .GT. 0 .AND. LECLEX .LE. 5 ) THEN
C                 NOM DU TYPE DE PLSVO
                  NMTOBJ = NMTYOB( LECLEX )
C                 NOM DU PLSVO
                  CALL NMOBNU( NMTOBJ, NUOBJ, NMOBJ )
C                 LEXIQUE DU PLSVO
                  CALL LXLXOU( NTMN(LECLEX), NMOBJ, NTPLSV, MNPLSV )
                  IF( NTPLSV .LE. 0 ) GOTO 11
                  IF( LECLEX .LT. 5 ) THEN
                     IF( LECLEX .GT. 1 ) THEN
C                       LE TABLEAU 'NSEF' DE CE LSV
                        CALL LXTSOU( NTPLSV, 'NSEF', NTNSEF, MNNSEF )
                        IF( NTNSEF .LE. 0 ) GOTO 11
                     ELSE
C                       LE TABLEAU 'NSEF' DE CE POINT N'EXISTE PAS
                        NTNSEF = 0
                        MNNSEF = 0
                     ENDIF
C                    LE TABLEAU 'XYZSOMMET' DE CE PLSV
                     CALL LXTSOU( NTPLSV, 'XYZSOMMET', NTXYZS,MNXYZS)
                     IF( NTXYZS .LE. 0 ) GOTO 11
C
C                    TRACE DU PLSV
                     IAVF   = IAVFAC
                     IAVFAC = 0
                     NCL    = NCOUAL
                     NCOUAL = NCBLAN
                     NCF    = NCOUAF
                     NCOUAF = NCBLAN
CCC                  CALL T1SOBJ( LECLEX, NMOBJ,  NUOBJ,
CCC     %                         NTNSEF, MNNSEF, NTXYZS, MNXYZS )
                     CALL XVVOIR
                     NCOUAF = NCF
                     NCOUAL = NCL
                     IAVFAC = IAVF
CCC                  ELSE
CCCC                    TRACE DE L'OBJET
CCC                     CALL T1OBJE( NMOBJ )
                  ENDIF
               ENDIF
C
C              SAUVEGARDE POUR FAIRE UN SEUL TRACE
               NMIN00 = NMIN0
C
            ENDIF
            GOTO 11
C
         ENDIF
C
      ELSE
C
C        ACTIVATION DE LA SAISIE A LA SOURIS OU AU CLAVIER
C        SAISIE D'UN POINT PAR CLIC DE LA SOURIS OU PAR CLAVIER PHYSIQUE
         CALL SAIPTCSU( NOCODE, NX, NY, NASCII )
C
      ENDIF
C
      IF( NOCODE .EQ. 0 ) THEN
C        ABANDON
         KLG(LHKLG) = '@;'
         GOTO 9900
      ENDIF

      IF( NOCODE .GT. 0 ) THEN
C
C        CLIC RELACHE DANS L'ECRAN GRAPHIQUE EN (NX,NY) PIXELS
         IF( LECLEX .GT. 0 .AND. LECLEX .LE. 5 ) THEN

C           EN ATTENTE D'UN NOM DANS UN LEXIQUE PLSVO
C           -----------------------------------------
C           RECHERCHE DU PLUS PROCHE OBJET DU CLIC
C           OBJET POSSIBLE DE TYPE D'OBJET
            DISTMI = DIST0
            NUOBJ  = 0
C
C           ADRESSE DE DEBUT DES ITEMS DE CES OBJETS PRESENTS SUR L'ECRAN
            MNIT = MNITEM( LECLEX )
            IF( MNIT .LE. 0 ) GOTO 10
            MOTS = MCN( MNIT )
C           MNIT+2 POINTEUR SUR LE NOMBRE ACTUEL D'ITEMS
            DO 90 I=1,MCN(MNIT+2)
C              LES COORDONNEES ECRAN DE L'OBJET
               MNIT = MNIT + MOTS
               DIST = (NX-MCN(MNIT)) ** 2  + (NY-MCN(MNIT+1)) ** 2
               IF( DIST .GT. 40000 ) GOTO 90
C              CLIC DANS LE CERCLE DE RAYON 200 PIXELS
               IF( DIST .LT. DISTMI ) THEN
C                 OBJET PLUS PROCHE
                  NUOBJ  = MCN(MNIT+2)
                  DISTMI = DIST
               ENDIF
 90         CONTINUE
C
            IF( NUOBJ .GT. 0 ) THEN

               IF( LECLEX .GT. 0 .AND. LECLEX .LE. 5 ) THEN
C                 PLSVO RETROUVE DE NOM
                  CALL TAMSOU( NTMN(LECLEX) , MNLX )
C                 L'ADRESSE DU 1-ER MOT DE CET OBJET DANS LE LEXIQUE
                  MN = MNLX + MCN( MNLX ) * NUOBJ
C                 LE NOM CONVERTI D'ENTIERS EN CARACTERES
                  CALL ENTNOM( MCN(MNLX+2), MCN(MN), KLG(LHKLG) )
C                 AJOUT DU ;
                  MN = NUDCNB( KLG(LHKLG) )
                  KLG(LHKLG)(MN+1:MN+1) = ';'
                  GOTO 9900
               ENDIF

ccc               IF( LECLEX .GE. 6 .AND. LECLEX .LE. 9 ) THEN
cccC                 NUOBJ est un NUMERO de SOMMET ou ARETE ou FACE ou EF3D
ccc                  WRITE( KLG(LHKLG)(1:10), '(I10)' ), NUOBJ
cccC                 AJOUT DU ;
ccc                  MN = NUDCNB( KLG(LHKLG) )
ccc                  KLG(LHKLG)(MN+1:MN+1) = ';'
ccc                  GOTO 9900
ccc               ENDIF

            ENDIF
C           PAS DE DESIGNATION SUFFISAMMENT PRECISE
C
         ELSE IF( NBLGRC(NRMENU) .GT. 0 ) THEN
C
C           LECTURE D'UN ENTIER PARMI PLUSIEURS POSSIBILITES D'UN MENU
C           ----------------------------------------------------------
            IF( (NX.GE.XRECT(NRMENU)) .AND.
     &          (NX.LE.XRECT(NRMENU)+DXRECT(NRMENU)) .AND.
     &          (NY.GE.YRECT(NRMENU)) .AND.
     &          (NY.LE.YRECT(NRMENU)+DYRECT(NRMENU)) ) THEN
C              CLIC DANS LE MENU : NUMERO DE LA CASE ?
               NUMERO = INT( (NY-YRECT(NRMENU)) / DYLGRC(NRMENU) + 1.0 )
               NUMERO = NAMENU( NUMERO )
               IF( NUMERO .EQ. INMENU+1 ) THEN
C                 DEMANDE DE DOCUMENTATION
                  KLG(LHKLG) = '?;'
               ELSE IF( NUMERO .EQ. INMENU-1 ) THEN
C                 ABANDON
                  KLG(LHKLG) = '@;'
               ELSE IF( NUMERO .NE. INMENU ) THEN
C                 L'ENTIER ASSOCIE A L'OPTION
                  WRITE( KLG(LHKLG)(1:12) , '(I12)' ) NUMERO
                  KLG(LHKLG)(13:13) = ';'
               ENDIF
               GOTO 9900
            ENDIF
C           CLIC EN DEHORS DES LIMITES DU RECTANGLE
         ENDIF
         GOTO 10
      ENDIF
C
C     LE CARACTERE CARLU A ETE LU AU CLAVIER PHYSIQUE
      CARLU = CHAR( NASCII )
      IF( NASCII .EQ. 13  ) THEN
C        CR RETOUR A LA LIGNE
         GOTO 9900
      ENDIF
C
CCC      WRITE(IMPRIM,*)
CCC      WRITE(IMPRIM,*) 'SACLAV: LECTURE DU CARACTERE:',CARLU
CCC      WRITE(IMPRIM,*) 'SACLAV: NUMERO ASCII        :',NASCII
C
C
C     LE NUMERO DU DERNIER CARACTERE NON BLANC DE LA LIGNE DE SAISIE
      NDC = NUDCNB( KLG(LHKLG) ) + NBLANC
C
C     EST CE LA TOUCHE BACKSPACE ou SUPPRESSION
      IF( NASCII .EQ. 8 .OR. NASCII .EQ. 127 ) THEN
C        OUI : PERTE D'UN CARACTERE  (8=BACKSPACE)
         IF( NDC .GT. 0 ) THEN
            KLG(LHKLG)(NDC:NDC) = ' '
         ENDIF
         NBLANC = 0
         GOTO 400
      ELSE IF( NASCII .EQ. 172 .OR. CARLU .EQ. ' ' ) THEN
C        OUI : AJOUT D'UN BLANC
         IF( NDC .GE. NBCALI ) GOTO 9000
         KLG(LHKLG)(NDC+1:NDC+1) = ' '
         NBLANC = NBLANC + 1
         GOTO 400
      ELSE IF( CARLU .EQ. '^' .OR. NASCII .EQ. 168 ) THEN
C        PASSAGE A LA LIGNE DU DESSUS DANS KLG
         LHKLG  = MAX( 1, LHKLG - 1 )
         NBLANC = 0
         GOTO 400
      ELSE IF( CARLU .EQ. '@' .OR. CARLU .EQ. '?' ) THEN
C        ABANDON DE LA LECTURE  OU  DEMANDE DE DOCUMENTATION
         KLG(LHKLG)(NDC+1:NDC+1) = CARLU
         KLG(LHKLG)(NDC+2:NDC+2) = ';'
         NBLANC = 0
         NSORTI = 1
      ELSE IF( NASCII .GE. 32 .AND. NASCII .LE. 126 ) THEN
C        UN CARACTERE STANDARD
C        RECHERCHE DU DERNIER CARACTERE NON BLANC
         KLG(LHKLG)(NDC+1:NDC+1) = CARLU
         NBLANC = 0
      ELSE
C        UN CARACTERE INUTILISABLE
         GOTO 10
      ENDIF
C
C     TRACE DE LA LIGNE DE SAISIE
C     ---------------------------
 400  NBLGRC(NRLGSA) = 1
      IF( NUDCNB( KLG(LHKLG) ) .GE. NBCALI ) THEN
C        LA LIGNE DE SAISIE EST EFFACEE
         CALL RECTEF( NRLGSA )
C        AFFICHAGE DE L'ERREUR
         KLG(LHKLG) = 'LIGNE TROP LONGUE. A RETAPER '
         CALL LIGSAI
C        REINITIALISATION DE LA LIGNE SAISIE
         KLG(LHKLG) = ' '
         GOTO 10
      ENDIF
      IF( NSORTI .NE. 0 ) GOTO 9900
      CALL LIGSAI
      GOTO 10
C
C     ERREUR
C     ------
 9000 NBLGRC(NRERR) = 1
      KERR(1) = 'ERREUR: TEXTE TROP LONG! A REDUIRE'
      CALL RECTTX( NRERR , KERR , 0 , NA )
      GOTO 400
C
C     SORTIE
C     ------
 9900 RETURN
      END
