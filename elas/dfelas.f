      SUBROUTINE DFELAS( KNOMOB , NERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    DEFINIR LES CARACTERISTIQUES PHYSIQUES DE L'ELASTICITE
C -----    D'UN OBJET LES FORCES EXERCEES ET LES CONDITIONS AUX LIMITES
C
C ENTREE :
C --------
C KNOMOB : NOM DE L'OBJET A ELASTICITER
C
C SORTIE :
C --------
C NERR   : CODE D'ERREUR =0 SI PAS DE PROBLEME
C                        =1 SINON
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE UPMC PARIS    JUILLET 1989
C2345X7..............................................................012
      PARAMETER        (LLIST1=9, LLIST2=9)
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NFDOCU,NFFRAP,NUNITE(27)
      include"./incl/gsmenu.inc"
      include"./incl/a_objet__topologie.inc"
      include"./incl/a_objet__definition.inc"
      include"./incl/a___fixation.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/xyzext.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
      CHARACTER*80      NOMTS
      CHARACTER*10      NMTYOB,KNM
      CHARACTER*24      KNOMOB,KNOM
C
      CHARACTER*15      LISTE1(LLIST1)
      CHARACTER*15      LISTE2(LLIST2)
      DATA LISTE1 / 'MASSE','YOUNG','DILATATION',
     %              'FORCE',' ','COEFDEPLACEMENT','CONTRINIT',
     %              'DEPLACTINIT','VITESSEINIT' /
      DATA LISTE2 / 'FORCE','FIXATION',' ',' ',' ',' ',' ',
     %              'DEPLACTINIT','VITESSEINIT' /
C
C     VERIFICATIONS
C     =============
C     RECHERCHE DU NOM DE L'OBJET DANS LE LEXIQUE DES OBJETS
      CALL LXLXOU( NTOBJE , KNOMOB , NTLXOB , MNLXOB )
C
C     S'IL N'EXISTE PAS RETOUR A LA DEMANDE DU NOM DE L'OBJET
      IF( NTLXOB .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = ' ERREUR : OBJET INCONNU ' // KNOMOB
         CALL LEREUR
         GOTO 9900
      ENDIF
C     LE NUMERO DE L'OBJET KNOMOB DANS LE LEXIQUE
      CALL LXNMNO( NTOBJE , KNOMOB , NUOBJE , MNLXOB )
C
C     RECHERCHE DU TABLEAU DEFINITION
      CALL LXTSOU( NTLXOB , 'DEFINITION' , NTDFOB , MNDFOB )
C     S'IL N'EXISTE PAS RETOUR A LA DEMANDE DU NOM DE L'OBJET
      IF( NTDFOB .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'ERREUR: DEFINITION INCONNUE OBJET ' // KNOMOB
         CALL LEREUR
         GOTO 9900
      ENDIF
C
C     RECHERCHE DU TABLEAU XYZSOMMET
      CALL LXTSOU( NTLXOB , 'XYZSOMMET' , NTXYZ , MNXYZ )
C     S'IL N'EXISTE PAS RETOUR A LA DEMANDE DU NOM DE L'OBJET
      IF( NTDFOB .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'ERREUR : XYZSOMMET INCONNU OBJET ' // KNOMOB
         CALL LEREUR
         GOTO 9900
      ENDIF
      NBSOM  = MCN( MNXYZ + WNBSOM )
      NBCOOR = MCN( MNXYZ + WBCOOR )
C
C     L'OBJET DOIT-IL ETRE STRUCTURE EN SOUS-DOMAINES ?
      IF ( MCN(MNDFOB+WDOUNO) .EQ. 1 ) THEN
        GO TO 101
      END IF
C     LE TABLEAU TOPOLOGIE EST RECHERCHE
      CALL LXTSOU( NTLXOB , 'TOPOLOGIE' , NTTOPO , MNTOPO )
      IF( NTTOPO .LE. 0 ) THEN
         NBLGRC(NRERR) = 2
         KERR(1) = 'ERREUR: OBJET ' // KNOMOB
         KERR(2) = 'SANS INTERPOLATION DEFINIE'
         CALL LEREUR
         GOTO 9900
      ENDIF
C
      IF( INTERA .GE. 1 ) THEN
C        TRACE DE L'OBJET SI CONSOLE GRAPHIQUE INTERACTIVE
         CALL EFFACE
         CALL MIMXPT( NBCOOR, NBSOM, MCN(MNXYZ+WYZSOM),  COOEXT )
         CALL REDUIRE
         CALL VISEE0
         CALL T1OBJE( KNOMOB )
         CALL AUGMENTER
         CALL CLICSO
      ENDIF
C
CCCC     LECTURE DES DONNEES INTERNES AU NIVEAU DE L'OBJET
CCC      CALL SDDEF1( KNOMOB , LISTE1 , 'ELASTICITE' , 'coefelas' )
C
C     BOUCLE SUR LES OBJETS INTERNES
C     ==============================
      IERR   = 0
      NERR   = 0
      NBTYEL = MCN( MNTOPO + WBTYEL )
      NBOBIN = MCN( MNTOPO + WBOBIN )
      MN     = MNTOPO + WMTYEL + NBTYEL - 2
      NDIM   = 2
      DO 200 I=0,NBOBIN-1
C
C        RECHERCHE DE L'OBJET DANS LE LEXIQUE
         MN   = MN + 2
         NYOB = MCN( MN )
         IF( NYOB .EQ. 4 ) NDIM=3
         NUOB = MCN( MN + 1 )
         CALL LXNLOU( NTMN(NYOB) , NUOB , NTOB , MNOB )
         KNM = NMTYOB( NYOB )
C        LE NOM DE L'OBJET
         CALL NMOBNU( KNM , NUOB , KNOM )
         IF( NTOB .LE. 0 ) THEN
            NBLGRC(NRERR) = 1
            KERR(1) = ' ERREUR : OBJET INCONNU ' // KNOM
            CALL LEREUR
            GOTO 9900
         ENDIF
C
         IF( INTERA .GE. 1 ) THEN
C           SI CONSOLE GRAPHIQUE INTERACTIVE
            CALL LXTSOU( NTOB, 'NSEF', NTNSEF, MNNSEF )
            CALL LXTSOU( NTOB, 'XYZSOMMET', NTXYZS, MNXYZS )
            IF( MNNSEF.GT.0 .AND. MNXYZS.GT.0 ) THEN
               CALL EFFACE
               IF( NYOB .EQ. 3 ) THEN
C                 TRACE DE LA SURFACE
                  CALL TRAFAC( KNOM, NUOB, MNNSEF, MNXYZS )
               ELSE IF( NYOB .EQ. 4 ) THEN
C                 TRACE DU VOLUME
                  CALL TRACUB( KNOM, NUOB, MNNSEF, MNXYZS )
               ENDIF
            ENDIF
         ENDIF
C
C        LECTURE EFFECTIVE DES DONNEES "INTERNES"
         CALL SDDEF2( NUOBJE,NYOB,NUOB,LISTE1,LLIST1,
     %               'ELASTICITE','coefelas', 1 )
 200  CONTINUE
C
C     BOUCLE SUR LES OBJETS "AUX LIMITES" OU FRONTALIERS
C     ==================================================
      NERR   = 0
      NBOBCL = MCN( MNTOPO + WBOBCL )
      MN     = MNTOPO + WMTYEL + NBTYEL - 2 + MOTVAR(13)*NBOBIN
      DO 300 I=0,NBOBCL-1
C        LE TYPE DE L'OBJET
         MN   = MN + 2
         NYOB = MCN( MN )
         NUOB = MCN( MN + 1 )
         CALL LXNLOU( NTMN(NYOB) , NUOB , NTOB , MNOB )
         KNM = NMTYOB( NYOB )
         CALL NMOBNU( KNM , NUOB , KNOM )
         IF ( NTOB .LE. 0 ) THEN
            NBLGRC(NRERR) = 1
            KERR(1) = ' ERREUR : OBJET INCONNU ' // KNOM
            CALL LEREUR
            GOTO 9900
         ENDIF
C
         IF( INTERA .GE. 1 ) THEN
            CALL EFFACE
            IF( NDIM .EQ. 3 ) THEN
C              TRACE DE L'OBJET SI CONSOLE GRAPHIQUE INTERACTIVE
               CALL TOPLSV( KNOMOB, NYOB, KNOM, IERR )
               IERR = 0
            ELSE
               CALL LXTSOU( NTOB, 'XYZSOMMET', NTXYZS, MNXYZS )
               IF( MNXYZS.GT.0 ) THEN
                  IF( NYOB .EQ. 1 ) THEN
                     CALL T21SOM( KNOM, NUOB, MNXYZS )
                  ELSE IF( NYOB .EQ. 2 ) THEN
                     CALL LXTSOU( NTOB, 'NSEF', NTNSEF, MNNSEF )
                     IF( MNNSEF .GT. 0 ) THEN
                        CALL T21ARE( KNOM, NUOB, MNNSEF, MNXYZS )
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
C
C        RECHERCHE DU 1-ER BLANC DANS LES NOMS
         L1 = INDEX( KNM , ' ' )
         IF( L1 .GT. 0 ) THEN
            L1 = L1 - 1
         ELSE
            L1 = LEN( KNM )
         ENDIF
         L2 = INDEX( KNOM , ' ' )
         IF( L2 .GT. 0 ) THEN
            L2 = L2 - 1
         ELSE
            L2 = LEN( KNOM )
         ENDIF
C
C        IMPRESSION D'UN MESSAGE
         NBLGRC(NRHIST) = 1
         KHIST(1) = KNM(1:L1) // ': ' // KNOM(1:L2)
         CALL LHISTO
         L3 = NUDCNB( KHIST(1) )
         WRITE( NFFRAP, * ) '{ ' , KHIST(1)(1:L3), ' }'
         WRITE( IMPRIM, * )
         WRITE( IMPRIM, * ) KHIST(1)(1:L3)
C
C        LE NOM DU TABLEAU A ENTRER
 210     CALL LIMTCL( 'cl_elas' , NMTCL )
         IF( NMTCL .LE. 0 ) GO TO 300
C
C        LE NOM DU TABLEAU TMS A REMPLIR
         NOMTS = '~>' // KNM(1:L1) // '>' // KNOM(1:L2) //
     %            '>' // LISTE2(NMTCL)
         L3 = INDEX( NOMTS , ' ' ) - 1
         CALL MOTSTD( '~>>>' // LISTE2(NMTCL) , NOMTS(1:L3) , IERR )
         IF( IERR .NE. 0 ) THEN
            NERR = NERR + 1
         ENDIF
C
C        LE TABLEAU DOIT IL ETRE SUPPRIME ?
         CALL LXTSOU( NTOB , LISTE2(NMTCL) , NTTS , MNTS )
C        LE TYPE DE CONDITIONS AUX LIMITES
         IF( NTTS .GT. 0 ) THEN
C
C           ATTENTION : TOUS LES TYPES DOIVENT ETRE A LA MEME ADRESSE
C                       DANS LE TABLEAU DESCRIPTEUR
C
            IF( IERR .NE. 0 .OR. MCN(MNTS+WTFIXA) .EQ. 0 ) THEN
C               LE TYPE EST NUL => DESTRUCTION DU TABLEAU
                CALL LXTSDS( NTOB , LISTE2(NMTCL) )
                NBLGRC(NRERR) = 2
                KERR(1) = LISTE2(NMTCL)
                KERR(2) = 'CONDITION AUX LIMITES SUPPRIMEE'
                CALL LERESU
                GOTO 210
            ENDIF
         ENDIF
C
C        AFFICHAGE DES VALEURS DU TABLEAU LU
         IF( IERR .EQ. 0 ) CALL AFTSTD( NOMTS(1:L3) )
         GOTO 210
 300  CONTINUE
      RETURN
C
C     ERREUR SI BATCH ou LECTURE DE FICHIER
 9900 IF( INTERA .GE. 3 ) RETURN
      CALL ARRET(100)
C
C     RESOLUTION PAR SOUS-DOMAINES
C     ============================
 101  CALL SDDFEL(KNOMOB,LISTE1,LLIST1,LISTE2,LLIST2,
     %            'ELASTICITE','coefelas','cl_elas',IERR)
      RETURN
      END
