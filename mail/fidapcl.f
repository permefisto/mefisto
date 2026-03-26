      SUBROUTINE FIDAPCL
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  DONNEES ET FICHIER DE SORTIE fidapcl DES CONDITIONS AUX LIMITES
C -----  D'UN OBJET 2D POUR LE LOGICIEL FIDAP
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEURS: PERRONNET RICOU ANALYSE NUMERIQUE UPMC PARIS         MAI 1998
C23456---------------------------------------------------------------012
      PARAMETER         (MXTYEL=7)
      PARAMETER         (MOPAGE=512)
      include"./incl/gsmenu.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/donthe.inc"
      include"./incl/donele.inc"
      include"./incl/a_objet__topologie.inc"
      include"./incl/a_objet__definition.inc"
      include"./incl/a___arete.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___xyznoeud.inc"
      include"./incl/a___xyzpoint.inc"
      include"./incl/a___npef.inc"
      include"./incl/msvaau.inc"
      include"./incl/ponoel.inc"
C
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE      (MCN(1), RMCN(1), DMCN(1))
C
      INTEGER           NOOBVC, NOOBSF(6), NOOBLA(12), NOOBPS(8)
      INTEGER           MNDOEL(4), MXDOEL(4)
      CHARACTER*4       NOMELE(2)
      CHARACTER*24      KNOMOB
      CHARACTER*24      KMALIGNE
C    O.R
      DOUBLE PRECISION  DPARAF(4),VAL
      INTEGER           NOFOCL
C
C     NOM_DE_L'OBJET  LE CARACTERE @ POUR FINIR
 100  CALL INVITE( 45 )
      NCVALS = 0
      CALL LIRCAR( NCVALS , KNOMOB )
      IF( NCVALS .EQ. -1 ) RETURN
      ILEXIS = 1
C
C     RECHERCHE DU NOM DE L'OBJET DANS LE LEXIQUE DES OBJETS
      IF( NTOBJE .LE. 0 ) RETURN
      CALL LXLXOU( NTOBJE , KNOMOB , NTLXOB , MNLXOB )
C
      IF( NTLXOB .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'FIDAPCL: OBJET INCONNU ' // KNOMOB
         CALL LEREUR
         CALL LXIM( NTOBJE )
         GOTO 100
      ENDIF
C
C     LECTURE DE LA LIGNE CONDITION AUX LIMITES
 200  CALL INVITE( 40 )
      NCVALS = 0
      CALL LIRCAR( NCVALS , KMALIGNE )
      IF( NCVALS .EQ. -1 ) RETURN
C
C     RECHERCHE DU NOM DE L'OBJET DANS LE LEXIQUE DES OBJETS
      IF( NTLIGN .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'FIDAPCL: AUCUNE LIGNE CONNUE'
         CALL LEREUR
         RETURN
      ENDIF
      CALL LXLXOU( NTLIGN , KMALIGNE , NTLXLG , MNLXLG )
C
      IF( NTLXLG .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'FIDAPCL: LIGNE INCONNUE ' // KMALIGNE
         CALL LEREUR
         GOTO 200
      ENDIF
C
C     LECTURE DE L'ENTIER NOMBRE DE DONNEES EN CHAQUE DE LA LIGNE
      CALL INVITE( 72 )
      NCVALS = 0
      CALL LIRENT( NCVALS , NBCL )
      IF( NCVALS .EQ. -1 ) RETURN
      IF( NBCL .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'FIDAPCL: LE NOMBRE DE CL DOIT ETRE >0'
         CALL LEREUR
         GOTO 200
      ENDIF
C
C     VERIFICATION QUE LES FONCTIONS EXISTENT (ET AFFECTATION DE LEUR NUMERO)
C
      CALL LXNMNO(NTFONC,'CL',NOFOCL,MNLX)
      IF (NOFOCL.LE.0) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'FIDAPCL: FONCTION CL INCONNUE'
         CALL LEREUR
      ENDIF
C
C     RECHERCHE D'UNE UNITE LIBRE
      CALL TRUNIT( NF )
C     SAUVEGARDE SUR FICHIER NF
      OPEN( UNIT=NF, FILE='fidapcl', STATUS='UNKNOWN',
     %      FORM='FORMATTED', ACCESS='SEQUENTIAL',
CCC     %      FORM='FORMATTED', ACCESS='APPEND',
     %      ERR=9999)
 3    READ(NF,*,END=4)
      GOTO 3
C
C     QUELQUES INITIALISATIONS
 4    ITEMP  = 0
      IERR   = 0
C     LE NOMBRE DE MOTS D'UNE VARIABLE REEL2
      MOREE2 = MOTVAR(6)
C
C     PROTECTION DES ADRESSES POUR EVITER DES PROBLEMES LORS
C     DE LA DESTRUCTION DES TABLEAUX
      MNNODL = 0
      DO 5 I=1,4
         MNDOEL(I) = 0
         MXDOEL(I) = 0
 5    CONTINUE
C
C     RECHERCHE DU TABLEAU DEFINITION DE L'OBJET
      CALL LXTSOU( NTLXOB , 'DEFINITION' , NTDFOB , MNDFOB )
C
C     S'IL N'EXISTE PAS RETOUR A LA DEMANDE DU NOM DE L'OBJET
      IF( NTDFOB .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'ERREUR: DEFINITION INCONNUE OBJET ' // KNOMOB
         CALL LEREUR
         GOTO 100
      ENDIF
C
C     L'OBJET EST-IL STRUCTURE EN SOUS-DOMAINES ?
      NDOUNO = MCN(MNDFOB+WDOUNO)
C
      IF ( NDOUNO .EQ. 1 ) THEN
C
C        RESOLUTION PAR LA METHOOE DES SOUS-DOMAINES
C        ===========================================
         NBLGRC(NRERR) = 1
         KERR(1) = 'ERREUR : METHODE DES SOUS-DOMAINES NON SUPPORTEE'
         CALL LEREUR
         GOTO 100
C
      ELSE IF ( NDOUNO .EQ. 2 ) THEN
C
C        RESOLUTION PAR LA METHOOE DES JOINTS
C        ====================================
         NBLGRC(NRERR) = 1
         KERR(1) = 'ERREUR : METHODE DES JOINTS NON SUPPORTEE'
         CALL LEREUR
         GOTO 100
C
      ENDIF
C
C     SINON : RESOLUTION CLASSIQUE
C     ============================
C     ADRESSAGE DES ADRESSES DES TABLEAUX NPEF"... DE CET OBJET
      MNELEM = 0
      CALL TNMCDC( 'ENTIER' , 2*MXTYEL , MNELEM )
      MNTELE = MNELEM + MXTYEL
C
C     RECHERCHE DES TABLEAUX SOMMETS NOEUDS POINTS ASSOCIES A L'OBJET
C     ===============================================================
      CALL NDPGEL( NTLXOB , NTTOPO , MNTOPO ,
     %             NTXYZP , MNXYZP , NTXYZN , MNXYZN ,
     %             NBTYEL , MCN(MNTELE) , MCN(MNELEM) , IERR )
      IF( IERR .NE. 0 ) RETURN
C
      NDPGST = MCN( MNTOPO + WDPGST )
      NBOBIN = MCN( MNTOPO + WBOBIN )
      NBOBCL = MCN( MNTOPO + WBOBCL )
C
C     INITIALISATIONS DE VARIABLES ET AFFICHAGES
C     ==========================================
C     LE NOMBRE DE POINTS DU MAILLAGE DE L'OBJET
      NBPOI = MCN( MNXYZP + WNBPOI )
C     NDIM LA DIMENSION 2 OU 3 DE L'ESPACE DES COORDONNEES
      CALL DIMCOO( NBPOI , MCN(MNXYZP+WYZPOI) , NDIM )
C     LE NOMBRE DE NOEUDS DU MAILLAGE DE L'OBJET
      NBNOEU = MCN(MNXYZN+WNBNOE)
C     LE NOMBRE TOTAL DE DEGRES DE LIBERTE THERMIQUES
      NTDL   = NBNOEU
C
      WRITE(IMPRIM,10210) NDIM,NBNOEU
10210 FORMAT(' DIMENSION 2 OU 3 DE L''ESPACE',T32,'=',I6/
     %' NOMBRE DE NOEUDS'                    ,T32,'=',I6)
C
C    CREATION D'UN TABLEAU POUR EVITER D'ECRIRE 2 FOIS LE MEME NEOUD
C
      CALL TNMCDC('ENTIER',NBNOEU,NODOUB)
      IF (NODOUB.LE.0) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'ERREUR FIDAPCL: PB D''ALLOCATION DE MEMOIRE'
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF
C     NBMIS EST L'ADRESSE DU DERNIER NOEUD DEJA ECRIT
      NBMIS=NODOUB-1
C
C     LA BOUCLE SUR LES TYPES D'ELEMENTS FINIS
C     ----------------------------------------
      CALL LXNMNO( NTLIGN , KMALIGNE , MALIGNE , MN )
      IF( MALIGNE .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'ERREUR : MALIGNE INCONNU'
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF
C
C     ECRITURE DU NOM DE L'OBJET ET DE LA LIGNE
      WRITE(NF,12000) KNOMOB, KMALIGNE
12000 FORMAT(1X,A24,1X,A24)
      NUTTEF = 0
      DO 2000 NOTYEL = 1 , NBTYEL
C
C        L'ADRESSE DU TABLEAU NPEF"TYPE_EF
         MNELE = MCN( MNELEM - 1 + NOTYEL )
C        LE NUMERO DU TYPE DE L'ELEMENT FINI
         NUTYEL = MCN( MNELE + WUTYEL )
C
C        LE NOMBRE D'ELEMENTS FINIS DE CE TYPE
         NBELEM = MCN( MNELE + WBELEM )
C
C        L'ADRESSE MCN DES NOEUDS ET POINTS GEOMETRIQUES DES ELEMENTS
         MNNDEL = MNELE + WUNDEL
         MNPGEL = MNNDEL
         IF( NDPGST .GE. 2 ) THEN
            MNPGEL = MNPGEL + NBELEM * MCN(MNELE+WBNDEL)
         ENDIF
C
C        LES CARACTERISTIQUES DE L'ELEMENT FINI
         CALL ELNUNM( NUTYEL , NOMELE )
         CALL ELTYCA( NUTYEL )
C
C        LA BOUCLE SUR LES ELEMENTS DE CE TYPE NUTYEL
C        ============================================
         DO 1900 NUELEM = 1 , NBELEM
C           LE NUMERO GLOBAL DE CET EF
            NUTTEF = NUTTEF + 1
C
C           LES NOEUDS DEGRE 1 DE L'ELEMENT
C           -------------------------------
            CALL EFNOEU( MNELE , NUELEM , NBNDEL , MCN(MNNODL) )
C
C           LES POINTS GEOMETRIQUES DE L'ELEMENT
C           ------------------------------------
CCC            CALL EFPOGE( MNELE , NUELEM , NBPGEF , MCN(IAPOEF) )
C
C           LE NUMERO DE VOLUME  DE L'EF
C           LE NUMERO DE SURFACE DES FACES   DE L'EF
C           LE NUMERO DE LIGNE   DES ARETES  DE L'EF
C           LE NUMERO DE POINT   DES SOMMETS DE L'EF
C           ----------------------------------------
            CALL EFPLSV( MNELE  , NUELEM ,
     %                   NOVCEL , NOSFEL , NOLAEL , NOPSEL ,
     %                   NOOBVC , NOOBSF , NOOBLA , NOOBPS )
C
            DO 9000 I=1,NARET
C                NUMERO DE LA LIGNE DE L'ARETE I
                 NL = NOOBLA(I)
                 IF( NL .NE. MALIGNE ) GOTO 9000
C                NUMERO DES 2 SOMMETS DE L'ARETE
                 IF( I .NE. NARET ) THEN
                    NS2 = I + 1
                 ELSE
                    NS2 = 1
                 ENDIF
                 NS1 = MCN(MNNODL-1+I)
                 NS2 = MCN(MNNODL-1+NS2)
                 DO 9088 L=1,2
                    IF( L .EQ. 1 ) THEN
                       NS = NS1
                    ELSE
                       NS = NS2
                    ENDIF
C
C                   TESTER SI LE NOEUD 1 A DEJA ETE ECRIT
                    DO 9010 K=NBMIS,NODOUB,-1
                      IF (NS.EQ.MCN(K)) GOTO 9088
 9010               CONTINUE
                    NBMIS = NBMIS + 1
                    MCN(NBMIS)=NS
C
C                   PUIS L'ECRIRE AVEC LES CONDITIONS A LA FRONTIERE
                    IACE = MNXYZP + WYZPOI + (NS-1) * 3
C                   XYZNS1= RMCN(IACE+0,+1,+2))
                    DPARAF(1)=RMCN(IACE)
                    DPARAF(2)=RMCN(IACE+1)
                    DPARAF(3)=RMCN(IACE+2)
                    DPARAF(4)=0D0
                     CALL FONVAL(NOFOCL,4,DPARAF,NCODEV,VAL)
                  IF (NCODEV.EQ.0) THEN
                       NBLGRC(NRERR) = 1
                        KERR(1) =
     %                'PB LORS DE L''EXECUTION DE LA FONCTION CL'
                        CALL LEREUR
                     ENDIF
                    WRITE(NF,19100) NS,VAL
C                   SANS PASSAGE A LA LIGNE (RC) SVP
      	    DO 9015 NOCL=1,NBCL-1
                       DPARAF(4)=NOCL
                        CALL FONVAL(NOFOCL,4,DPARAF,NCODEV,VAL)
                     IF (NCODEV.EQ.0) THEN
                          NBLGRC(NRERR) = 1
                           KERR(1) =
     %                   'PB LORS DE L''EXECUTION DE LA FONCTION CL'
                           CALL LEREUR
                        ENDIF
                       WRITE(NF,'(E15.7,$)') VAL
C                      $ => SANS PASSAGE A LA LIGNE
 9015               CONTINUE
                    WRITE (NF,'(1X)')
C                   ICI AVEC PASSAGE A LA LIGNE
 9088            CONTINUE
19100            FORMAT(I8,E15.7,$)
C
 9000       CONTINUE
 1900    CONTINUE
 2000 CONTINUE
      CALL TNMCDS('ENTIER',NBNOEU,NODOUB)
C
      CLOSE( NF )
      RETURN
C
C     PB A L'OUVERTURE DU FICHIER
 9999 NBLGRC(NRERR) = 2
      KERR(1) = 'PB A L''OUVERTURE DU FICHIER fidapcl'
      CALL LEREUR
      GOTO 100
      END
