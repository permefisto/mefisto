      SUBROUTINE DFNLSE( KNOMOB, NERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    DEFINIR LES DONNEES DE NLSEQ
C -----    D'UN OBJET ET LES FORCES EXERCEES
C
C          INPUT DATA of the NLSE CHARACTERISTICS of EVERY MATERIAL
C          INPUT DATA of the NLSE BOUNDARY CONDITIONS
C
C ENTREE :
C --------
C KNOMOB : NOM DE L'OBJET A TRAITER   (WORK OBJECT NAME)
C
C SORTIE :
C --------
C NERR   : CODE D'ERREUR =0 SI PAS DE PROBLEME ,=1 SINON (ERROR CODE)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET TEXAS A & M at DOHA QATAR        Fevrier 2011
C23456---------------------------------------------------------------012
      PARAMETER        (LLIST1=9, LLIST2=8)
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NFDOCU,NFFRAP,NUNITE(27)
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/mecoit.inc"
      include"./incl/a_objet__topologie.inc"
      include"./incl/a_objet__definition.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___fixation.inc"
      include"./incl/xyzext.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/pp.inc"
      COMMON             MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (RMCN(1),MCN(1))
      CHARACTER*80      NOMTS
      CHARACTER*10      NMTYOB,KNM
      CHARACTER*24      KNOMOB,KNOM
      CHARACTER*1       KNOJEU
      CHARACTER*2       KAJOUT
      CHARACTER*24      KLISTE
      CHARACTER*16      LISTE1(LLIST1)
      CHARACTER*11      LISTE2(LLIST2)
      DATA LISTE1/'MASSE', 'CONDUCTIVITE',  '    ',  'FORCE', '    ',
     %            'COEFNLSE', '    ', 'ONDEINIT', 'VITESSEANGULAIRE' /
      DATA LISTE2/'    ', 'FORCE', 'FIXATION', '    ',
     %            '    ',   '    ',   '    ', 'ONDEINIT' /
C
C     REMARQUE: Alfa = Conductivite
C               Beta = COEFNLSE = Coef devant la temperature
C               Force    a 2 composantes
C               Fixation a 2 composantes
C               OndeInit a 2 composantes
C
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10000) KNOMOB
      ELSE
         WRITE(IMPRIM,20000) KNOMOB
      ENDIF
10000 FORMAT(/1X,80('*')/
     %' LECTURE DES DONNEES NLSE DE L''OBJET ',A/1X,80('*'))
20000 FORMAT(/1X,80('*')/
     %' READING of NLSE INPUT DATA of the OBJECT ',A/1X,80('*'))
C
C     VERIFICATIONS
C     =============
C     RECHERCHE DU NOM DE L'OBJET DANS LE LEXIQUE DES OBJETS
C     SEARCH THE OBJECT NAME?
      CALL LXLXOU( NTOBJE , KNOMOB , NTLXOB , MNLXOB )
C
C     S'IL N'EXISTE PAS RETOUR A LA DEMANDE DU NOM DE L'OBJET
      IF( NTLXOB .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'ERREUR: OBJET INCONNU ' // KNOMOB
         ELSE
            KERR(1) = 'ERROR: UNKNOWN OBJECT ' // KNOMOB
         ENDIF
         CALL LEREUR
         GOTO 9900
      ENDIF
C     LE NUMERO DE L'OBJET KNOMST DANS LE LEXIQUE
      CALL LXNMNO( NTOBJE , KNOMOB , NUOBJE , MNLXOB )
C     RECHERCHE DU TABLEAU DEFINITION
      CALL LXTSOU( NTLXOB , 'DEFINITION' , NTDFOB , MNDFOB )
C     S'IL N'EXISTE PAS RETOUR A LA DEMANDE DU NOM DE L'OBJET
      IF( NTDFOB .LE. 0 ) THEN
         NBLGRC(NRERR) = 2
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'ERREUR: OBJET ' // KNOMOB
            KERR(2) = 'SANS DEFINITION => A REDEFINIR'
         ELSE
            KERR(1) = 'ERROR: OBJECT ' // KNOMOB
            KERR(2) = 'WITHOUT DEFINITION => DEFINE AGAIN'
         ENDIF
         CALL LEREUR
         GOTO 9900
      ENDIF
C
C     RECHERCHE DU TMS (SECONDARY MEMORY ARRAY) XYZSOMMET
C     SEARCH THE TMS (SECONDARY MEMORY ARRAY) XYZSOMMET
      CALL LXTSOU( NTLXOB , 'XYZSOMMET' , NTXYZ , MNXYZ )
C     S'IL N'EXISTE PAS RETOUR A LA DEMANDE DU NOM DE L'OBJET
      IF( NTXYZ .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
          KERR(1)='ERREUR: tms XYZSOMMET INCONNU pour l''OBJET '//KNOMOB
         ELSE
            KERR(1) = 'ERROR: UNKNOWN XYZSOMMET tms for OBJECT '//KNOMOB
         ENDIF
         CALL LEREUR
         GOTO 9900
      ENDIF
      NBSOM  = MCN( MNXYZ + WNBSOM )
      NBCOOR = MCN( MNXYZ + WBCOOR )
C
C     RESOLUTION CLASSIQUE
C     CLASSIC FORMULATION
C     ====================
C     LE TMS TOPOLOGIE EST RECHERCHE
C     SEARCH THE TMS (SECONDARY MEMORY ARRAY) TOPOLOGIE
      CALL LXTSOU( NTLXOB , 'TOPOLOGIE' , NTTOPO , MNTOPO )
      IF( NTTOPO .LE. 0 ) THEN
         NBLGRC(NRERR) = 3
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'ERREUR: OBJET ' // KNOMOB
            KERR(2) = 'SANS INTERPOLATION DEFINIE'
            KERR(3) = '=> INTERPOLATION a DEFINIR dans MAILLER'
         ELSE
            KERR(1) = 'ERROR: OBJECT ' // KNOMOB
            KERR(2) = 'WITHOUT DEFINED INTERPOLATION'
            KERR(3) = '=> DEFINE INTERPOLATION with MAILLER'
         ENDIF
         CALL LEREUR
         GOTO 9900
      ENDIF
C
      IF( INTERA .GE. 1 ) THEN
C        TRACE DE L'OBJET SI CONSOLE GRAPHIQUE INTERACTIVE
C        DRAWING OF THE OBJECT IF SUFFICENT INTERACTIVITY
         CALL EFFACEMEMPX
         CALL MIMXPT( NBCOOR, NBSOM, MCN(MNXYZ+WYZSOM), COOEXT )
         CALL REDUIRE
         CALL VISEE0
         CALL T1OBJE( KNOMOB )
         CALL AUGMENTER
         CALL CLICSO
      ENDIF
C
C     NOMBRE DE JEUX DE DONNEES = 1
      NBJEUX = 1
C
CCCC     LECTURE DES DONNEES INTERNES AU NIVEAU DE L'OBJET GLOBAL
CCCC     READING OF INTERNAL DATA FOR THE WHOLE OBJECT
CCCC     ========================================================
CCC      CALL SDDEF1( KNOMOB , LISTE1 , 'THERMIQUE', 'coefther' )
C
C     BOUCLE SUR LES OBJETS INTERNES
C     LOOP ON THE INTERIOR MATERIALS
C     ------------------------------
      IERR   = 0
      NERR   = 0
      NBTYEL = MCN( MNTOPO + WBTYEL )
      NBOBIN = MCN( MNTOPO + WBOBIN )
      MN     = MNTOPO + WMTYEL + NBTYEL - 2
      NDIM   = 2
      DO 100 I=0,NBOBIN-1
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
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'ERREUR: OBJET INCONNU ' // KNOM
            ELSE
               KERR(1) = 'ERROR: UNKNOWN OBJECT ' // KNOM
            ENDIF
            CALL LEREUR
            GOTO 9900
         ENDIF
C
         IF( INTERA .GE. 1 ) THEN
C           SI CONSOLE GRAPHIQUE INTERACTIVE
            CALL LXTSOU( NTOB, 'NSEF', NTNSEF, MNNSEF )
            CALL LXTSOU( NTOB, 'XYZSOMMET', NTXYZS, MNXYZS )
            IF( MNNSEF.GT.0 .AND. MNXYZS.GT.0 ) THEN
               CALL EFFACEMEMPX
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
C        LECTURE DES DONNEES AU NIVEAU DE L'OBJET "INTERNE"
C        READING OF INTERNAL DATA FOR THE I-th MATERIAL
C         variable COEFNLSE 'SCHRODINGER Equation INPUT' entier
C       ( 1 : 'Rho MASS Density'
C       , 2 : 'Alfa Laplacian Coefficient'
C       , 4 : 'FOmega Complex Force Second Member 2 Components'
C       , 6 : 'Beta Coefficient of (ur**2+ui**2) ur and ui'
C       , 8 : 'INITIAL Complex Wave 2 Components '
C       , 9 : 'Omega  Rotation Angular Velocity' ) ;
C        ======================================================
         CALL SDDEF2( NUOBJE, NYOB, NUOB, LISTE1, LLIST1,
     %               'NLSE','coefnlse', NBJEUX )
 100  CONTINUE
C
C     BOUCLE SUR LES OBJETS "AUX LIMITES" OU FRONTALIERS
C     LOOP ON THE BOUNDARY OBJECTS as SLP in 3D
C     --------------------------------------------------
      NERR   = 0
      NBOBCL = MCN( MNTOPO + WBOBCL )
      MN     = MNTOPO + WMTYEL + NBTYEL - 2 + MOTVAR(13)*NBOBIN
      DO 400 I=0,NBOBCL-1
C        LE TYPE DE L'OBJET
         MN   = MN + 2
         NYOB = MCN( MN )
         NUOB = MCN( MN + 1 )
         CALL LXNLOU( NTMN(NYOB) , NUOB , NTOB , MNOB )
         KNM = NMTYOB( NYOB )
         CALL NMOBNU( KNM , NUOB , KNOM )
         IF ( NTOB .LE. 0 ) THEN
            NBLGRC(NRERR) = 1
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'ERREUR : Point ou Ligne ou Surface INCONNU '
     %                 // KNOM
            ELSE
               KERR(1) = 'ERROR: UNKNOWN Point or Line or Surface '
     %                 // KNOM
            ENDIF
            CALL LEREUR
            GOTO 9900
         ENDIF
C
         IF( INTERA .GE. 1 ) THEN
            CALL EFFACEMEMPX
            IF( NDIM .EQ. 3 ) THEN
C              TRACE DE L'OBJET FRONTIERE SI CONSOLE GRAPHIQUE INTERACTIVE
C              DRAWING OF THE BOUNDARY OBJECT IF SUFFICENT INTERACTIVITY
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
C        SEARCH THE FIRST SPACE ' ' IN THE NAME
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
C        DISPLAY A MESSAGE
         NBLGRC(NRHIST) = 1
C
ccc         IF( LANGAG .GT. 0 ) THEN
cccC           EN ANGLAIS LINE REMPLACE LIGNE
ccc            M = INDEX( KNM(1:L1), 'LIGNE' )
ccc            IF( M .GT. 0 ) THEN
ccc               L1  = 4
ccc               KNM = 'LINE'
ccc            ENDIF
ccc         ENDIF
C
         KHIST(1) = KNM(1:L1) // ': ' // KNOM(1:L2)
         CALL LHISTO
         L3 = NUDCNB( KHIST(1) )
         WRITE( NFFRAP, * ) '{ ' , KHIST(1)(1:L3), ' }'
         WRITE( IMPRIM, * )
         WRITE( IMPRIM, * ) KHIST(1)(1:L3)
C
C        BOUCLE SUR LES JEUX DE DONNEES DE CHACUN DES OBJETS AUX LIMITES
C        LOOP ON THE NBJEUX DATA GAMES OF EVERY BOUNDARY OBJECT
C        ===============================================================
         DO 300 K = 1, NBJEUX
            NOJEU = K
            WRITE( KNOJEU, '(I1)' ) NOJEU
            WRITE(IMPRIM,*)
            NBLGRC(NRERR) = 1
            IF( LANGAG .EQ. 0 ) THEN
           KERR(1)='ATTENTION: Entree du JEU de DONNEES Numero='//KNOJEU
            ELSE
            KERR(1)='ATTENTION: INPUT DATA GAME Number='//KNOJEU
            ENDIF
            CALL LERESU
C
C           LE NOM DU TABLEAU A ENTRER EST UNE OPTION de cl_nlse
C           THE NAME of the ARRAY to INPUT Is an option of cl_nlse
C            variable CL_NLSE 'SCHRODINGER Equation Boundary Condition' entier
C          ( 2 : 'FGamma IMPOSED Complex WAVE FLUX  at Boundary'
C          , 3 : 'WAVE_D IMPOSED Complex WAVE VALUE at Boundary'
C          , 8 : 'WAVE_0 INITIAL Complex WAVE VALUE at Boundary') ;
 210        CALL LIMTCL( 'cl_nlse' , NMTCL )
            IF( NMTCL .LE. 0 ) GO TO 300
C
C           LE NOM DU TABLEAU TMS A REMPLIR
            NOMTS = '~>' // KNM(1:L1) // '>' // KNOM(1:L2) //
     %               '>' // LISTE2(NMTCL)
            L3 = INDEX( NOMTS, ' ' ) - 1
C
            CALL MOTSTD( '~>>>' // LISTE2(NMTCL), NOMTS(1:L3), IERR )
            IF( IERR .NE. 0 ) THEN
               NERR = NERR + 1
            ENDIF
C
C           LE TABLEAU EST IL CREE ?
            CALL LXTSOU( NTOB, LISTE2(NMTCL), NTTS, MNTS )
C           LE TYPE DE CONDITIONS AUX LIMITES
C           THE TYPE OF BOUNDARY CONDITIONS
            IF( NTTS .GT. 0 ) THEN
C
               IF( NBJEUX .GT. 1 ) THEN
C                 AJOUT DE "NOJEU au NOM du TMS si NBJEUX>1
                  KAJOUT = '"' // KNOJEU
                  NOMTS(L3+1:L3+2) = KAJOUT
                  L3 = L3+2
C                 LE TABLEAU EST RENOMME AVEC "NOJEU
                  L4 = NUDCNB( LISTE2(NMTCL) )
                  KLISTE = LISTE2(NMTCL)(1:L4) // KAJOUT
                  CALL LXNMNM( NTOB, LISTE2(NMTCL),  KLISTE )
               ELSE
                  KAJOUT = '  '
               ENDIF
C
               L4 = NUDCNB( LISTE2(NMTCL) )
               KLISTE = LISTE2(NMTCL)(1:L4) // KAJOUT
               L4 = NUDCNB( KLISTE )
               CALL LXTSOU( NTOB, KLISTE(1:L4), NTTS, MNTS )
C
C              ATTENTION : TOUS LES TYPES DOIVENT ETRE A LA MEME ADRESSE
C                          DANS LE TABLEAU DESCRIPTEUR
C              ATTENTION : ALL TYPES MUST BE AT THE SAME ADDRESS IN THE
C                          DESCRIPTOR ARRAY
               IF( IERR .NE. 0 .OR. MCN(MNTS+WTFIXA) .EQ. 0 ) THEN
C                 LE TYPE EST NUL => DESTRUCTION DU TABLEAU
                  CALL LXTSDS( NTOB, KLISTE(1:L4) )
                  NBLGRC(NRERR) = 2
                  KERR(1)(1:L3) = NOMTS(1:L3)
                  IF( LANGAG .EQ. 0 ) THEN
                     KERR(2)='Cette CONDITION aux LIMITES est SUPPRIMEE'
                  ELSE
                     KERR(2)='This BOUNDARY CONDITION is DELETED'
                  ENDIF
                  CALL LERESU
                  GOTO 210
               ENDIF
            ENDIF
C
C           AFFICHAGE DES VALEURS DU TABLEAU LU
C           DISPLAY THE READ VALUES
            IF( IERR .EQ. 0 ) CALL AFTSTD( NOMTS(1:L3) )
            GOTO 210
 300     CONTINUE
 400  CONTINUE
C
C     SUPPRESSION DE L'AFFICHAGE DU NOM D'OBJET PREMIER
C     SUPPRESS THE DISPLAY OF THE PRIME OBJECT
      CALL RECTEF( NRHIST )
      NBLGRC(NRHIST) = 1
      CALL LHISTO
      GOTO 9999
C
C     ERREUR ERROR
 9900 NERR = 1
C
 9999 RETURN
      END
