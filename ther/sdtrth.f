      SUBROUTINE SDTRTH( KNOMSD, NTLXSD, MNDFSD, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACER LES TEMPERATURES ET/OU LES FLUX
C -----    D'UN OBJET 2D OU 3D DE NOM KNOMSD APRES CALCUL
C          LA METHODE DES SOUS-DOMAINES
C
C ENTREE :
C --------
C KNOMSD : NOM DE L'OBJET
C NTLXSD : NUMERO DE SON LEXIQUE DANS LE LEXIQUE DES OBJETS
C MNDFSD : ADRESSE MCN DU TABLEAU 'DEFINITION' DE CET OBJET
C
C SORTIE :
C --------
C IERR   : 0 SI PAS D'ERREUR , NON NUL SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY   ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 1994
C23456---------------------------------------------------------------012
      PARAMETER   (MXTYEL=7, LIGCON=0)
      include"./incl/gsmenu.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/msvaau.inc"
      include"./incl/trvari.inc"
      include"./incl/xyzext.inc"
      include"./incl/a_objet__definition.inc"
      include"./incl/a_objet__topologie.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___xyznoeud.inc"
      include"./incl/a___xyzpoint.inc"
      include"./incl/a___npef.inc"
      include"./incl/a___vecteur.inc"
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              GRADMX,COIN(6,2),COIN0(6,2),COIM(6,2),ACOIN(6,2)
      DOUBLE PRECISION  DTMIN0, DTMAX0
      CHARACTER*10      NMTYOB,KNOMTY
      CHARACTER*80      KNOMOB,KNOM
      CHARACTER*(*)     KNOMSD
      LOGICAL           DEJAVU
C
C     NOM DE L'OBJET "SOUS-DOMAINES"
C     ------------------------------
      L = INDEX( KNOMSD , ' ' )
      IF( L .GT. 0 ) THEN
         L = L - 1
      ELSE
         L = LEN( KNOMSD )
      ENDIF
      KNOMSD = KNOMSD(1:L) // '_SD'
C     RECHERCHE DE L'OBJET KNOMSD DANS LE LEXIQUE DES OBJETS
      CALL LXLXOU( NTOBJE , KNOMSD , NTLXSD , MNLXSD )
C     S'IL N'EXISTE PAS RETOUR A LA DEMANDE DU NOM DE L'OBJET
      IF( NTLXSD .LE. 0 ) THEN
         NBLGRC(NRERR) = 2
         KERR(1) = 'ERREUR: OBJET INCONNU'
         KERR(2) = KNOMSD
         CALL LEREUR
         CALL LXIM( NTOBJE )
         IERR = 1
         RETURN
      ENDIF
C     LE TABLEAU DEFINITION DE L'OBJET KNOMSD
      CALL LXTSOU( NTLXSD , 'DEFINITION' , NTDFSD , MNDFSD )
C     IL N'EXISTE PAS
      IF( NTDFSD .LE. 0 ) THEN
         NBLGRC(NRERR) = 2
         KERR(1) = 'ERREUR: TMS DEFINITION INCONNU de l''OBJET'
         KERR(2) = KNOMSD
         CALL LEREUR
         CALL ARRET(100)
      ENDIF
C     LE NOMBRE DE SOUS-DOMAINES
      NBDOBJ = MCN(MNDFSD+WBDOBJ)
C
C     INITIALISATIONS POUR LA REMANENCE DES VALEURS
C     ---------------------------------------------
      NCAS   = 1
      NOPT   = 1
      CMFLEC = 2.5
      CMPGRA = 2.5
      CMPFLU = 2.5
      IERR   = 0
      NBISO  = 11
      DEJAVU = .FALSE.
C     LES VALEURS EXTREMES
      RMAX   = RINFO( 'GRAND' )
      GRADMX = 0.
      FLUXMX = 0.
      TMIN = RMAX
      TMAX = -TMIN
      XMI  = RMAX
      YMI  = XMI
      XMA  = -XMI
      YMA  = -YMI
C
C     RECHERCHE DU MIN ET MAX DES COORDONNEES DES NOEUDS
C     --------------------------------------------------
      MNDF = MNDFSD + WTYOBJ
      NBSD = 0
      DO 149 I=1,3
C        LE MINIMUM
         COIM(I,1) =  RMAX
C        LE MAXIMUM
         COIM(I,2) = -RMAX
149   CONTINUE
      DO 150 NO = 0 , NBDOBJ - 1
         NUTYSD = MCN(MNDF)
         NUOBSD = MCN(MNDF+1)
         KNOMTY = NMTYOB(NUTYSD)
         MNDF   = MNDF + 2
         IF (NUTYSD .EQ. 5 ) THEN
C
C           LE LEXIQUE DE L'OBJET
C
            NBSD = NBSD + 1
            CALL NMOBNU( KNOMTY , NUOBSD , KNOMOB )
            CALL LXLXOU( NTOBJE , KNOMOB , NTLXOB , MNLXOB )
            IF( NTLXOB .LE. 0 ) THEN
               NBLGRC(NRERR) = 2
               KERR(1) = ' ERREUR : OBJET INCONNU'
               KERR(2) = KNOMOB
               CALL LEREUR
               GOTO 9000
            ENDIF
C           RECHERCHE DES TEMPERATURES DE L'OBJET
            CALL  LXLXOU( NTLXOB , 'VECTEUR"TEMPERATURE' ,
     %                    NTTEMP , MNTEMP )
            IF( NTTEMP .LE. 0 ) THEN
               NBLGRC(NRERR) = 2
               KERR(1) = 'ERREUR : OBJET ' // KNOMOB
               KERR(2) = 'SANS TEMPERATURES CALCULEES'
               CALL LEREUR
               GOTO 9000
            ENDIF
C           NOMBRE DE DEGRES DE LIBERTE
            NTDL = MCN( MNTEMP + WBCOVE )
C           NOMBRE DE CAS
            NDSM = MCN( MNTEMP + WBVECT )
C
C           RECHERCHE DES TABLEAUX SOMMETS NOEUDS POINTS ASSOCIES A L'OBJET
C
            CALL TNMCDC( 'ENTIER' , 2*MXTYEL , MNELEM )
            MNTELE = MNELEM + MXTYEL
            CALL NDPGEL( NTLXOB , NTTOPO , MNTOPO ,
     %                   NTPOGE , MNPOGE , NTNOEU , MNNOEU ,
     %                   NBTYEL , MCN(MNTELE) , MCN(MNELEM) , IERR )
C           RECHERCHE DU MIN ET MAX DES COORDONNEES DES POINTS
            IF( COIN0(3,1) .NE. 0. .OR. COIN0(3,2) .NE. 0. ) THEN
                NDIMLI = 3
            ELSE
                NDIMLI = 2
            ENDIF
            DO 170  I=1,3
               COIM(I,1) = MIN( COIM(I,1) , COIN0(I,1) )
               COIM(I,2) = MAX( COIM(I,2) , COIN0(I,2) )
 170        CONTINUE
C           RECHERCHE DES TEMPERATURES MIN ET MAX
            CALL MXVECT( NTDL, NDSM, MCN(MNTEMP+WECTEU) ,
     %                   DTMIN0, NOEMIN, NCAMIN, DTMAX0, NOEMAX, NCAMAX)
            TMIN0 = REAL( DTMIN0 )
            TMAX0 = REAL( DTMAX0 )
            TMIN = MIN ( TMIN , TMIN0 )
            TMAX = MAX ( TMAX , TMAX0 )
C           RECHERCHE DU MAXIMUM EN VALEUR ABSOLUE DES GRADIENTS
            CALL MXGRSD( NTLXOB, NBTYEL, MNELEM,
     %                   GRAMX0, COIN0, IERR )
            IF( GRAMX0 .LE. 0 ) GRAMX0 = MAX( GRAMX0 , 1.0 )
            GRADMX = MAX ( GRADMX , GRAMX0 )
C           RECHERCHE DU MAXIMUM EN VALEUR ABSOLUE DES FLUX
            CALL MXFLSD( NTLXOB, NBTYEL, MNELEM,  FLUMX0, COIN0, IERR )
            IF( FLUMX0 .LE. 0 ) FLUMX0 = MAX( FLUMX0 , 1.0 )
            FLUXMX = MAX ( FLUXMX , FLUMX0 )
            CALL TNMCDS( 'ENTIER' , 2*MXTYEL , MNELEM )
         END IF
 150  CONTINUE
      X    = ( TMAX - TMIN ) / 100.
      IF( X .EQ. 0 ) X = 1.0
      TMIN = TMIN + X
      TMAX = TMAX - X
      ECAMAX = 0.
      DO 171  I=1,NDIMLI
         COIN(I,1) = COIM(I,1)
         COIN(I,2) = COIM(I,2)
         ECAMAX = MAX( ECAMAX , COIN(I,2) - COIN(I,1) )
 171  CONTINUE
      IF( ECAMAX .LE. 0. ) ECAMAX = 1.
C     UNE MARGE DE 6%
      ECAMAX = ECAMAX * 0.53
      DO 160 I=1,NDIMLI
         COORMI    = ( COIN(I,1) + COIN(I,2) ) * 0.5
         COIN(I,1) = COORMI - ECAMAX
         COIN(I,2) = COORMI + ECAMAX
 160  CONTINUE
C
C     MEMORISATION DES DIMENSIONS DE L'OBJET
C     --------------------------------------
      DO 165 I=1,NDIMLI
         ACOIN(I,1) = COIN(I,1)
         ACOIN(I,2) = COIN(I,2)
 165  CONTINUE
C
C     LA FENETRE EST EFFACEE
      CALL EFFACE
C
C     TRACE DES TEMPERATURES OU DES FLUX
C     ----------------------------------
 10   CALL XVEPAISSEUR( 1 )
      CALL XVTYPETRAIT( LIGCON )
      NBLGRC(NRHIST) = 1
      L = NUDCNB( KNOMOB )
      KHIST(1) = 'OBJET: ' // KNOMOB(1:L)
      CALL LHISTO
C
      CALL LIMTCL( 'tempgrad' , NMTC1 )
      IF( NMTC1 .LE. 0 ) GOTO 9000
      NMTC2 = 1
      GOTO ( 100, 200, 300, 300, 10, 10, 10, 10, 10, 10 ) , NMTC1
C
C     NUMERO DU CAS A VISUALISER
C     ==========================
 100  NCVALS = 4
      CALL INVITE( 84 )
      CALL LIRENT( NCVALS, NCAS )
C     VALEUR FORCEE POUR L'INSTANT. 1 SEUL SECOND MEMBRE EST CALCULE
      NCAS = 1
      GOTO 10
C
C     TRACE DE LA SOLUTION NCAS
C     =========================
C     RECHERCHE DES TEMPERATURES DE L'OBJET
C     TRACE PAR SOUS-DOMAINE
C
C     OPTIONS DU TRACE DE LA TEMPERATURE NCAS
 200  CALL LIMTCL( 'tractemp' , NMTC2 )
      IF( NMTC2 .LE. 0 ) GOTO 100
      IF( NMTC2 .EQ. 3 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'ERREUR : OPTION IMPOSSIBLE '
         CALL LEREUR
         GOTO 200
      ENDIF
C
 300  NBTRTE = 0
      NBTRFL = 0
      TMI = RINFO( 'GRAND' )
      TMA = - TMI
      IF (DEJAVU) THEN
         DO 166 I=1,NDIMLI
            COIN(I,1) = ACOIN(I,1)
            COIN(I,2) = ACOIN(I,2)
 166     CONTINUE
      ELSE
         DEJAVU = .TRUE.
      END IF
C
      MNDF = MNDFSD + WTYOBJ
      DO 210 NO = 0 , NBDOBJ - 1
         NUTYSD = MCN(MNDF)
         NUOBSD = MCN(MNDF+1)
         KNOMTY = NMTYOB(NUTYSD)
         MNDF   = MNDF + 2
         IF (NUTYSD .EQ. 5 ) THEN
C
C           LE LEXIQUE DE L'OBJET
            CALL NMOBNU( KNOMTY , NUOBSD , KNOMOB )
            CALL LXLXOU( NTOBJE , KNOMOB , NTLXOB , MNLXOB )
C
C           RECHERCHE DES TABLEAUX SOMMETS NOEUDS POINTS ASSOCIES A L'OBJET
C           ADRESSES DES TABLEAUX ELEMENTS DE CET OBJET
            CALL TNMCDC( 'ENTIER' , 2*MXTYEL , MNELEM )
            MNTELE = MNELEM + MXTYEL
            CALL NDPGEL( NTLXOB , NTTOPO , MNTOPO ,
     %                   NTPOGE , MNPOGE , NTNOEU , MNNOEU ,
     %                   NBTYEL , MCN(MNTELE) , MCN(MNELEM) , IERR )
            IF( IERR .NE. 0 ) GOTO 10
            NBTRTE = NBTRTE + 1
            NBTRFL = NBTRFL + 1
C
            GOTO  ( 210, 220, 230, 240 ) , NMTC1
C
C           TRACE DE LA TEMPERATURE
C           =======================
 220           GOTO ( 221, 222 ), NMTC2
C              TRACE DES ISOTHERMES (LIGNES EN 2D SEULEMENT)
 221           CALL SDTRIT( KNOMOB, NTLXOB, NBTRTE, NBSD,
     %                      NBTYEL, MNTOPO, MNELEM, MNPOGE,
     %                      NCAS  , NBISO , TMIN  , TMAX,
     %                      TMI, XMI, YMI, TMA, XMA, YMA )
               CALL TNMCDC( 'ENTIER', 2*MXTYEL, MNELEM )
               GOTO 210
C
C              TRACE DES ZONES DE COULEURS ISOTHERMES
 222           IF( NOIBLA .EQ. 1 ) GO TO 220
               CALL SDTRZT( KNOMOB, NTLXOB, NBTRTE,  NBSD,
     %                      NBTYEL, MNTOPO, MNELEM, MNPOGE,
     %                      NCAS  , NBISO, TMIN, TMAX,
     %                      TMI, XMI, YMI, TMA, XMA, YMA )
               CALL TNMCDC( 'ENTIER', 2*MXTYEL, MNELEM )
               GOTO 210
C
C           TRACE DU GRADIENT DE LA TEMPERATURE
C            ===================================
 230            CALL SDTRGT( NTLXOB, NBTRFL,  NBSD, KNOMSD,
     %                       NBTYEL, MNTOPO, MNELEM, MNPOGE,
     %                       NCAS  , NOPT  , CMFLEC, CMPGRA, GRADMX )
                CALL TNMCDC( 'ENTIER', 2*MXTYEL, MNELEM )
                GOTO 210
C
C           TRACE DU GRADIENT DE LA TEMPERATURE (FLUX)
C            =========================================
 240            CALL SDTRFL( -1,1,KNOMOB, NTLXOB, NBTRFL,  NBSD, KNOMSD,
     %                      NBTYEL, MNTOPO, MNELEM, MNPOGE, NDIMLI,
     %                      NCAS  , NOPT  , CMFLEC, CMPFLU, FLUXMX )
                CALL TNMCDC( 'ENTIER', 2*MXTYEL, MNELEM )
C
         ENDIF
 210  CONTINUE
C
C      AFFICHAGE DES TEMPERATURES EXTREMES
C     -----------------------------------
      IF (NMTC2 .EQ. 2) THEN
         CALL SYMBOLE2D( NCCYAN, XMI, YMI, '+m' )
         CALL SYMBOLE2D( NCCYAN, XMA, YMA, '+M' )
      ENDIF
C
C     LES LIGNES SONT REMISES A LEUR EPAISSEUR NORMALE
c     ------------------------------------------------
      CALL XVEPAISSEUR( 1 )
      CALL XVTYPETRAIT( LIGCON )
C
C     FIN DU TRACE
C     ------------
      IF ( NMTC1.EQ.1 ) THEN
         KNOM = 'CARTE       DES TEMPERATURES'
         WRITE( KNOM(7:10) , '(I4)' ) NCAS
         CALL TRFINS( KNOM )
      ELSE IF ( NMTC1.EQ.2 ) THEN
         CALL TRFINS( 'GRADIENTS DE TEMPERATURE' )
      END IF
      GO TO 10
C
C     FIN DE L'EXECUTION
 9000 CALL XVEPAISSEUR( 1 )
      CALL XVTYPETRAIT( LIGCON )
      RETURN
      END
