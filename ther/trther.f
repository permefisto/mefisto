      SUBROUTINE TRTHER( KNOMOB, NTYSOL, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACER LES VECTEURS SOLUTIONS + GRADIENTS + FLUX
C -----    D'UN OBJET 1D ou 2D ou 3D ou 6D de NOM KNOMOB
C          DRAW the SOLUTION VECTORS + GRADIENTS + FLUXES
C          of an 1D or 2D or 3D or 6D OBJECT of NAME KNOMOB
C
C ENTREE :
C --------
C KNOMOB : NOM DE L'OBJET                      (OBJECT NAME)
C NTYSOL : NO DU TYPE DE SOLUTION A TRACER     (TYPE of SOLUTIONS)
C          1 VECTEURS TEMPERATURES             (TEMPERATURES)
C          2 VECTEURS PROPRES                  (EIGEN-SOLUTIONS)
C          3 VECTEURS DEPLACEMENTS D'UNE ONDE  (WAVE DISPLACEMENTS)
C          4 VECTEURS DU MODULE de l'ONDENLSE          (Complex WAVE}
C          5 VECTEURS DE LA PARTIE REELLE ONDENLSE     (Complex WAVE}
C          6 VECTEURS DE LA PARTIE IMAGINAIRE ONDENLSE (Complex WAVE}
C SORTIE :
C --------
C IERR   : 0 SI PAS D'ERREUR, NON NUL SINON    (ERROR CODE)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     OCTOBRE 1990
C MODIFS : ALAIN PERRONNET Laboratoire JLL   UPMC PARIS     OCTOBRE 2005
C MODIFS : ALAIN PERRONNET Laboratoire JLL   UPMC PARIS    NOVEMBRE 2006
C MODIFS : ALAIN PERRONNET LJLL UPMC PARIS St Pierre du Perray JUIN 2009
C MODIFS : ALAIN PERRONNET TIMS NTU TAIPEI TAIWAN           OCTOBRE 2009
C MODIFS : ALAIN PERRONNET LJLL UPMC & ST PIERRE DU PERRAY  OCTOBRE 2010
C MODIFS : ALAIN PERRONNET TEXAS A&M University at DOHA QATAR  MARS 2011
C MODIFS : ALAIN PERRONNET Saint Pierre du Perray              Mars 2021
C23456---------------------------------------------------------------012
      PARAMETER        (MXTYEL=7, LIGCON=0)
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/msvaau.inc"
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      include"./incl/xyzext.inc"
      include"./incl/donthe.inc"
      include"./incl/a_objet__definition.inc"
      include"./incl/a_objet__topologie.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___xyznoeud.inc"
      include"./incl/a___xyzpoint.inc"
      include"./incl/a___npef.inc"
      include"./incl/a___vecteur.inc"
      include"./incl/ctemps.inc"
C
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      DOUBLE PRECISION DMCN(1)
      REAL             RMCN(1)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      EQUIVALENCE      (MCN(1),RMCN(1),DMCN(1))
      CHARACTER*(*)     KNOMOB
      DOUBLE PRECISION  DTMIN,  DTMAX
      DOUBLE PRECISION  TECMOY, TECMIN, TECMAX, TEXMIN, TEXMAX
      INTEGER           NUMIOB(4), NUMAOB(4), MNDOEL(4), MXDOEL(4)
      type typ_dptab
         DOUBLE PRECISION, dimension(:), pointer :: dptab
      end type typ_dptab
      type( typ_dptab ),dimension(0:1) :: dptemp

      IF( INTERA .LE. 0 ) THEN
C        DEMANDE DE TRACE EN MODE BATCH => ARRET DE MEFISTO
         CALL ARRET( 100 )
      ELSE IF( INTERA .GE. 3 ) THEN
         LORBITE = 1
         NORBITE = 0
      ENDIF

      IF( NTYSOL .EQ. 3 ) THEN
         CALL TRONDE( KNOMOB, IERR )
         RETURN
      ENDIF
C
C     INITIALISATIONS POUR LA REMANENCE DES VALEURS
      MNERRE = 0
      MNTIMP = 0
      MNTEER = 0
      MNNPEF = 0
      MNVALP = 0
      MNMOWA = 0
      MOMOWA = 0
      NCA0   = 0
      NCAS   = 1
      NCAS0  = 1
      NCAS1  = 1
      TEMPS  = 0.0
      NOPT   = 1
      CMFLEC = 2.5
      CMPGRA = 0.0
      NOPROJ0= 0
      NOPROJ =-1
      IERR   = 0
C     LE NOMBRE DE MOTS D'UNE VARIABLE REEL2
      MOREE2 = MOTVAR(6)
C
      MNXYZ3Q1C = 0
      MOXYZ3Q1C = 0
      MNEF3Q1C  = 0
      MNTE3Q1C  = 0
      NDIMES    = 0
      NBCPIN    = 0
      MNXYZ6Q1C = 0
      MNTE6Q1C  = 0
      MNTIME    = 0
C
C     RETRIEVE THE OBJECT
      CALL LXLXOU( NTOBJE, KNOMOB, NTLXOB, MNLXOB )
      IF( NTLXOB .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'ERREUR: OBJET INCONNU ' // KNOMOB
         ELSE
            KERR(1) = 'ERROR: UNKNOWN OBJECT ' // KNOMOB
         ENDIF
         CALL LEREUR
         IERR = 1
         GOTO 9000
      ENDIF
C
C     RECHERCHE DU TABLEAU TMS DEFINITION DE L'OBJET
      CALL LXTSOU( NTLXOB, 'DEFINITION', NTDFOB, MNDFOB )
C
C     S'IL N'EXISTE PAS RETOUR A LA DEMANDE DU NOM DE L'OBJET
      IF( NTDFOB .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'ERREUR: DEFINITION INCONNUE OBJET ' // KNOMOB
         ELSE
            KERR(1)='ERROR: UNKNOWM DEFINITION of the OBJECT ' // KNOMOB
         ENDIF
         CALL LEREUR
         IERR = 2
         GOTO 9000
      ENDIF
C
C     L'OBJET EST-IL STRUCTURE EN SOUS-DOMAINES 2D?
C     STRUCTURE BY 2D SUB-DOMAINS?
      NDOUNO = MCN(MNDFOB+WDOUNO)
      IF (NDOUNO .GT. 0) THEN
C        TRACE PAR SOUS-DOMAINES
         CALL SDTRTH( KNOMOB, NTLXOB, MNDFOB, IERR )
         RETURN
      ENDIF
C
C     RESOLUTION CLASSIQUE
C     RECHERCHE DES TABLEAUX XYZSOMMET XYZNOEUD XYZPOINT NPEF" DE l'OBJET
C     FIND TMS XYZSOMMET XYZNOEUD XYZPOINT NPEF"xxxx  of the OBJECT
C     ===================================================================
      CALL MIMAOB( 1,      NTLXOB, MXDOTH, NTTOPO, MNTOPO,
     %             NTXYZP, MNXYZP, NTXYZN, MNXYZN,
     %             MXTYEL, NBTYEL, MTNPEF, MNNPEF,
     %             NUMIOB, NUMAOB,
     %             NDPGST, NBOBIN, MNOBIN, NBOBCL, MNOBCL,
     %             MXDOEL, MNDOEL, IERR )
C     ICI NUMIOB CONTIENT LES 4 NUMEROS MINIMA DES OBJETS
C         NUMAOB          LES 4 NUMEROS MAXIMA DES OBJETS
C     MNDOEL LES 4 ADRESSES MCN DES TABLEAUX DES ADRESSES DES
C     TABLEAUX DECRIVANT LA THERMIQUE DE L'OBJET COMPLET
C     NTTOPO : NUMERO      DU TMS 'TOPOLOGIE' DE L'OBJET
C     MNTOPO : ADRESSE MCN DU TMS 'TOPOLOGIE' DE L'OBJET
C     NTXYZP : NUMERO      DU TMS 'XYZPOINT'  DE L'OBJET
C     MNXYZP : ADRESSE MCN DU TMS 'XYZPOINT'  DE L'OBJET
C     NTXYZN : NUMERO      DU TMS 'XYZNOEUD'  DE L'OBJET
C     MNXYZN : ADRESSE MCN DU TMS 'XYZNOEUD'  DE L'OBJET
C     NBTYEL : NOMBRE DE TYPES D'ELEMENTS FINIS DU MAILLAGE
C     MTNPEF : ADRESSE MCN DU TMC CONTENANT LE NUMERO     DU TMS 'NPEF'
C              DES NBTYEL TYPES D'ELEMENTS FINIS
C     MNNPEF : ADRESSE MCN DU TMC CONTENANT L ADRESSE MCN DU TMS 'NPEF'
C              DES NBTYEL TYPES D'ELEMENTS FINIS
      IF( IERR .NE. 0 ) GOTO 10
C
C     NDPGST : CODE TRAITEMENT DES XYZ DES SOMMETS POINTS NOEUDS DU MAILLAGE
C              0 : NOEUDS=POINTS=SOMMETS
C              1 : NOEUDS=POINTS#SOMMETS
C              2 : NOEUDS#POINTS=SOMMETS
C              3 : NOEUDS#POINTS#SOMMETS
      NDPGST = MCN( MNTOPO + WDPGST )
C
C     NBCOOR : NOMBRE DE COORDONNEES DES NOEUDS=POINTS (3 ou 6)
      NBCOOR = MCN( MNXYZP + WBCOOP )
      NBSOM6 = MCN( MNXYZP + WNBPOI )
C
C     NDIM   DIMENSION EFFECTIVE DE L'ESPACE DES COORDONNEES DES NOEUDS 2 3 ou 6
C     NDIMLI DIMENSION EFFECTIVE DE L'ESPACE DES TRACES 2 ou 3
      IF( NBCOOR .LE. 3 ) THEN
C        CAS DES DIMENSIONS 1 2 et 3
         CALL DIMCOO( MCN(MNXYZP+WNBPOI), MCN(MNXYZP+WYZPOI), NDIMLI )
      ELSE
C        CAS DE LA DIMENSION 6 => TRACES EN 3D
         NDIMLI = 3
      ENDIF
      NDIMES = NDIMLI
C
C     NOMBRE DE NOEUDS DU MAILLAGE
      NBNOMA = MCN(MNXYZN+WNBNOE)
C
C     SAUVEGARDE DE L'ADRESSE DES POINTS
      MNXYZ6Q1C = MNXYZP
C
      IF( NTYSOL .EQ. 1 ) THEN
C
C        RECHERCHE DES VECTEURS TEMPERATURES DE L'OBJET
         CALL  LXTSOU( NTLXOB, 'VECTEUR"TEMPERATURE', NTTEMP, MNTEMP )
         IF( NTTEMP .LE. 0 ) THEN
C           ERREUR PAS DE VECTEUR A VISUALISER
            NBLGRC(NRERR) = 2
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'ERREUR: OBJET '// KNOMOB
               KERR(2) = 'VECTEURS SOLUTION NON CALCULEES'
            ELSE
               KERR(1) = 'ERROR: OBJECT '// KNOMOB
               KERR(2) = 'SOLUTION VECTORS ARE NOT COMPUTED'
            ENDIF
            CALL LEREUR
            IERR = 3
            GOTO 9000
         ELSE
C           MODE DE TRACE DES TEMPERATURES SOLUTIONS
            MODECO = 1
         ENDIF
C
      ELSE IF( NTYSOL .EQ. 2 ) THEN
C
C        RECHERCHE DES VALEURS ET VECTEURS PROPRES DE L'OBJET
         CALL LXTSOU( NTLXOB, 'VECTEUR"VALEURPROPRE', NTTEMP, MNTEMP )
         IF( NTTEMP .LE. 0 ) THEN
C           ERREUR PAS DE VECTEUR A VISUALISER
            NBLGRC(NRERR) = 2
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'ERREUR: OBJET '// KNOMOB
               KERR(2) = 'VECTEURS PROPRES NON CALCULES'
            ELSE
               KERR(1) = 'ERROR: OBJECT '// KNOMOB
               KERR(2) = 'EIGENVECTORS NOT COMPUTED'
            ENDIF
            CALL LEREUR
            IERR = 3
            GOTO 9000
         ENDIF
C        MODE DE TRACE DES VALEURS ET VECTEURS PROPRES
         MODECO = 2
         MNVALP = MNTEMP
C
      ELSE IF( NTYSOL .GE. 4 ) THEN
C
C        RECHERCHE DES VECTEURS DE L'ONDENLSE
         CALL  LXTSOU( NTLXOB, 'VECTEUR"ONDENLSE', NTTEMP, MNTEMP )
         IF( NTTEMP .LE. 0 ) THEN
C           ERREUR PAS DE VECTEUR A VISUALISER
            NBLGRC(NRERR) = 2
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'ERREUR: OBJET '// KNOMOB
               KERR(2) = 'VECTEURS ONDENLSE NON CALCULES'
            ELSE
               KERR(1) = 'ERROR: OBJECT '// KNOMOB
               KERR(2) = 'VECTORS"ONDENLSE NOT COMPUTED'
            ENDIF
            CALL LEREUR
            IERR = 3
            GOTO 9000
         ENDIF
C        MODE DE TRACE DES VECTEURS du MODULE de l'ONDENLSE COMPLEXE
         MODECO = 4 + NTYSOL
      ENDIF
C
      IF( NDIMLI .EQ. 1 .AND. MODECO .EQ. 1 ) THEN
C        CAS PARTICULIER DU 1D: MODE DE TRACE DES VECTEURS TEMPERATURES
         MODECO = 5
      ENDIF
C
C     SAUVEGARDE DE MNTEMP
      MNTE6Q1C = MNTEMP
C
C     NOMBRE TOTAL DE CAS
      NDSM = MCN( MNTEMP + WBVECT )
C
C     NOMBRE DE DEGRES DE LIBERTE
      NTDL = MCN( MNTEMP + WBCOVE )
      IF( NTYSOL .LT. 4 .AND. NTDL .NE. NBNOMA ) THEN
         NBLGRC(NRERR) = 4
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'ERREUR: OBJET ' // KNOMOB
            KERR(2) = 'NOMBRE DES NOEUDS DU MAILLAGE NON EGAL AU'
            KERR(3) = 'NOMBRE DES NOEUDS DE CALCUL DE LA TEMPERATURE'
            KERR(4) = 'RESOUDRE A NOUVEAU LE PROBLEME THERMIQUE'
         ELSE
            KERR(1) = 'ERROR: OBJECT ' // KNOMOB
            KERR(2) = 'MESH NODES NUMBER NOT EQUAL TO'
            KERR(3) = 'COMPUTED TEMPERATURE NODES NUMBER'
            KERR(4) = 'SOLVE AGAIN THE HEAT PROBLEM'
         ENDIF
         CALL LEREUR
         IERR = 3
         GOTO 9000
      ENDIF
C
      NBCPIN = MCN( MNTEMP + WBCPIN )
      IF( NBCPIN .GT. 0 .AND. MODECO .NE. 2 ) THEN
C        SI PB INSTATIONNAIRE DU PREMIER AU DERNIER VECTEUR
         NCAS0 = 1
         NCAS1 = NDSM
      ELSE
C        SINON LE DERNIER VECTEUR SOLUTION EST CHOISI PAR DEFAUT
         NCAS0 = NDSM
         NCAS1 = NDSM
      ENDIF
      NCAS = NCAS1
C
C     LES TEMPS DU CALCUL DES NBVECT VECTEURS OU LEUR NUMERO
      MNTIME = MNTEMP + WECTEU + NTDL * NDSM * MOREE2
      CALL LESTEMPS( KNOMOB, MNTEMP, NDSM, MNTIMP, IERR )
      IF( IERR .NE. 0 ) GOTO 9000
      MNTIME = MNTIMP - 1
C
cccC     AFFICHAGE DES VALEURS COMPLEMENTAIRES DES VECTEURS
ccc      WRITE(IMPRIM,*)
ccc      IF( LANGAG .EQ. 0 ) THEN
ccc         WRITE(IMPRIM,*) 'Les',NBCPIN,' Valeurs Complementaires des SOLU
ccc     %TIONS'
ccc      ELSE
ccc       WRITE(IMPRIM,*) 'The',NBCPIN,' Complementary Values of SOLUTIONS'
ccc      ENDIF
ccc      DO 5 N=1,NBCPIN,5
ccc         IF( N+4 .LE. NBCPIN ) THEN
ccc             NK = 4
ccc         ELSE
ccc             NK = NBCPIN - N
ccc         ENDIF
ccc         IF( LANGAG .EQ. 0 ) THEN
ccc            WRITE(IMPRIM,10005) (N+K, RMCN(MNTIME+N+K), K=0,NK)
ccc         ELSE
ccc            WRITE(IMPRIM,20005) (N+K, RMCN(MNTIME+N+K), K=0,NK)
ccc         ENDIF
ccc 5    CONTINUE
ccc10005 FORMAT(5(' Valeur',I5,' : ',G15.7,3X))
ccc20005 FORMAT(5(' Value' ,I5,' : ',G15.7,3X))
C
      IF( MODECO .GE. 8 ) THEN
C
C        CAS DE L'ONDE NLSE COMPLEXE: CALCUL DE SON MODULE
C        U = V + i W => |U|  et 2 DL => 1 DL
         IF( MNMOWA .GT. 0 ) CALL TNMCDS( 'REEL2', MOMOWA, MNMOWA )
         NTDLMO = NTDL/2
         MOMOWA = NTDLMO * NDSM
         CALL TNMCDC( 'REEL2', MOMOWA, MNMOWA )
C
         IF( MODECO .EQ. 8 ) THEN
C           CALCUL DU MODULE(NBNOEU,NDSM) DE L'ONDE(NBNOEU,2,NDSM)
            MNVE = ( MNTEMP + WECTEU - 1 ) / MOREE2
            MNMO = ( MNMOWA          - 1 ) / MOREE2
            DO NCAS=1,NDSM
               DO I=1,NTDLMO
                  MNVE = MNVE + 1
                  MNMO = MNMO + 1
                  DMCN(MNMO)=SQRT( DMCN(MNVE)**2+DMCN(MNVE+NTDLMO)**2 )
               ENDDO
               MNVE = MNVE + NTDLMO
            ENDDO
C
         ELSE IF( MODECO .EQ. 9 ) THEN
C           COPIE DE LA PARTIE REELLE DE L'ONDE(NBNOEU,2,NDSM)
            MNVE = ( MNTEMP + WECTEU - 1 ) / MOREE2
            MNMO = ( MNMOWA          - 1 ) / MOREE2
            DO NCAS=1,NDSM
               DO I=1,NTDLMO
                  MNVE = MNVE + 1
                  MNMO = MNMO + 1
                  DMCN(MNMO)= DMCN(MNVE)
               ENDDO
               MNVE = MNVE + NTDLMO
            ENDDO
C
         ELSE IF( MODECO .EQ. 10 ) THEN
C           COPIE DE LA PARTIE IMAGINAIRE DE L'ONDE(NBNOEU,2,NDSM)
            MNVE = ( MNTEMP + WECTEU - 1 ) / MOREE2 + NTDLMO
            MNMO = ( MNMOWA          - 1 ) / MOREE2
            DO NCAS=1,NDSM
               DO I=1,NTDLMO
                  MNVE = MNVE + 1
                  MNMO = MNMO + 1
                  DMCN(MNMO)= DMCN(MNVE)
               ENDDO
               MNVE = MNVE + NTDLMO
            ENDDO
C
         ENDIF
C
C        ATTENTION ON ECRASE MNTEMP ET NTDL!...
         NTDL   = NTDLMO
         MNTEMP = MNMOWA - WECTEU
C
      ENDIF
C
cccC     LA FENETRE EST EFFACEE
ccc      CALL EFFACEMEMPX
C
C     LE TITRE EST TRACE
      IAVTIT = 1
C
C     ****************************************************************
C     NO DU CAS OU TRACE DE LA TEMPERATURE OU DE SON GRADIENT OU FLUX?
C     ****************************************************************
 10   NDIMLI = NDIMES
      CALL RECTEF( NRHIST )
      CALL XVEPAISSEUR( 1 )
      CALL XVTYPETRAIT( LIGCON )
      NBLGRC(NRHIST) = 1
      IF( LANGAG .EQ. 0 ) THEN
         KHIST(1) = 'OBJET: ' // KNOMOB
      ELSE
         KHIST(1) = 'OBJECT: ' // KNOMOB
      ENDIF
      CALL LHISTO
C
      CALL LIMTCL( 'tempgrad', NMTCL )
      IF( NMTCL .LE. 0 ) GOTO 9000
      IF( NBCOOR .GT. 3 .AND. NMTCL .GT. 2 ) GOTO 10
      GOTO ( 100, 200, 300, 400, 500, 10, 200, 800, 900, 1000 ), NMTCL
C
C     NUMERO DU CAS A VISUALISER
C     ==========================
 100  CALL LIRNOCAS( NDSM, NCAS0, NCAS1, N )
      IF( N .LT. 0 ) GOTO 10
      IF( NCAS0 .EQ. NCAS1 ) THEN
         NCAS = NCAS1
      ELSE
         NCAS = NDSM
      ENDIF
C     LES TEMPS ONT ILS ETE STOCKES?
      IF( NBCPIN .GT. 0 ) THEN
C        OUI: LE TEMPS INITIAL EST CELUI DU VECTEUR"TEMPERATURE
         TEMPS = RMCN( MNTIME + NCAS )
         IF( MODECO .NE. 2 ) THEN
            IF( LANGAG .EQ. 0 ) THEN
               WRITE(IMPRIM,*) 'TRACE AU TEMPS ',TEMPS
            ELSE
               WRITE(IMPRIM,*) 'DRAWING AT TIME ',TEMPS
            ENDIF
         ELSE
            IF( LANGAG .EQ. 0 ) THEN
               WRITE(IMPRIM,*) 'VALEUR PROPRE ',TEMPS
            ELSE
               WRITE(IMPRIM,*) 'EIGENVALUE ',TEMPS
            ENDIF
         ENDIF
      ELSE
C        TEMPS INITIAL SUPPOSE NUL
         TEMPS = 0
      ENDIF
      GOTO 10
C
C     TRACE DE LA SOLUTION NCAS DE LA TEMPERATURE
C     ===========================================
 200  IF( NBCOOR .EQ. 6 ) THEN
C
C        MAILLAGE de 6CUBES
         CALL LIMTCL( 'proj6cub', NOPROJ )
         IF( NOPROJ .LT. 0 ) GOTO 10
         IF( NOPROJ .EQ. NOPROJ0 .AND. NCAS .EQ. NCA0 ) GOTO 201
         NOPROJ0 = NOPROJ
C
C        RECHERCHE DES MIN-MAX DES 6 COORDONNEES DE L'OBJET
         CALL MAJEXT( MNXYZP )
C
C        CONSTRUCTION DU SOUS-MAILLAGE EN 3-CUBES DU MAILLAGE DE 6-CUBES
C        EN PROJECTION SELON LA DIRECTION DEFINIE PAR NOPROJ
         MNTEMP = MNTE6Q1C
         MNXYZN = MNXYZ6Q1C
         MNXYZP = MNXYZ6Q1C
         CALL C6X123000( NOPROJ, KNOMOB, NTLXOB,
     %                   MCN(MTNPEF), MCN(MNNPEF), MNTEMP,
     %                   NTEF3Q1C, MNEF3Q1C, NTTE3Q1C, MNTE3Q1C,
     %                   MOXYZ3Q1C, MNXYZ3Q1C, NBTYEL, IERR )
         IF( IERR .NE. 0 ) GOTO 10
C
C        REMISE A JOUR
         CALL NDPGEL( NTLXOB, NTTOPO, MNTOPO,
     %                NTXYZP, MNXYZP, NTXYZN, MNXYZN,
     %                NBTYEL, MCN(MTNPEF), MCN(MNNPEF), IERR )
C
C        LES TABLEAUX DU MAILLAGE 3Q1C EXTRAIT DES 6CUBES
C        POUR ETRE HOMOGENE AVEC LA SUITE
         MNXYZN = MNXYZ3Q1C
         MNXYZP = MNXYZ3Q1C
         MNTEMP = MNTE3Q1C
         NTDL   = MCN( MNTE3Q1C + WBCOVE )
C        LE TEMPS SI PB INSTATIONNAIRE
         IF( NBCPIN .GT. 0 ) TEMPS = RMCN( MNTIME + NCAS )
      ENDIF
C
C     CALCUL DES TEMPERATURES MIN ET MAX, DE LEURS NOEUDS ET CAS
      CALL MXVECT( NTDL,  NDSM,   RMCN(MNTEMP+WECTEU),
     %             DTMIN, NOEMIN, NCAMIN, DTMAX, NOEMAX, NCAMAX )
      TMIN = REAL( DTMIN )
      TMAX = REAL( DTMAX )
      NCA0 = NCAS
      MNS  = MNXYZP + WYZPOI - 4
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10020) NCAMIN, DTMIN, NOEMIN,
     %                      (RMCN(MNS+3*NOEMIN+K),K=1,3)
         WRITE(IMPRIM,10021) NCAMAX, DTMAX, NOEMAX,
     %                      (RMCN(MNS+3*NOEMAX+K),K=1,3)
      ELSE
         WRITE(IMPRIM,20020) NCAMIN, DTMIN, NOEMIN,
     %                      (RMCN(MNS+3*NOEMIN+K),K=1,3)
         WRITE(IMPRIM,20021) NCAMAX, DTMAX, NOEMAX,
     %                      (RMCN(MNS+3*NOEMAX+K),K=1,3)
      ENDIF
C
10020 FORMAT(/' VECTEUR',I9,
     %        ' VALEUR MINIMALE=',G15.7,' au NOEUD',I12,' COOR ',3G15.7)
10021 FORMAT( ' VECTEUR',I9,
     %        ' VALEUR MAXIMALE=',G15.7,' au NOEUD',I12,' COOR ',3G15.7)
20020 FORMAT(/' VECTOR',I9,
     %        ' MINIMUM VALUE=',G15.7,' at NODE',I12,' COOR ',3G15.7)
20021 FORMAT( ' VECTOR',I9,
     %        ' MAXIMUM VALUE=',G15.7,' at NODE',I12,' COOR ',3G15.7)
C
C     ----------------------------------------------------
C     OPTIONS DU TRACE DE LA TEMPERATURE OU AUTRE SOLUTION
C     ----------------------------------------------------
 201  IF( NMTCL .EQ. 7 ) GOTO 290
C     AFFICHAGE DES TEMPERATURES DU CAS NON DEMANDE
C     UN TRACE  DES TEMPERATURES DU CAS EST DEMANDE
      NDIMLI = NDIMES
      IF( NDIMLI .EQ. 1 ) THEN
         CALL LIMTCL( 'tractem1', NMTCL0 )
         IF( NMTCL0 .LE. 0 ) GOTO 10
         GOTO ( 110, 120, 201, 201, 201, 201, 201, 280, 290 ), NMTCL0
      ELSE IF( NDIMLI .EQ. 2 ) THEN
         CALL LIMTCL( 'tractem2', NMTCL0 )
      ELSE
         CALL LIMTCL( 'tractem3', NMTCL0 )
      ENDIF
      IF( NMTCL0 .LE. 0 ) THEN
         IF( NBCOOR .EQ. 6 ) THEN
            GOTO 200
         ELSE
            GOTO 10
         ENDIF
      ENDIF
      GOTO ( 210, 220, 230, 240, 250, 201, 270, 280, 290 ), NMTCL0
C
C     TRACE DE Z=TEMPERATURE(t,x) AVEC X=X Y=TEMPS Z=TEMPERATURE
C     ================================================================
 110  MODECI = MODECO
      MODECO = 5
      CALL TR1DTEZ( KNOMOB, MODECO, NDIMLI,
     %              NBTYEL, MNNPEF, MNXYZP, NDPGST,
     %              NDSM,   NTDL,   RMCN(MNTEMP+WECTEU),
     %              RMCN(MNTIMP) )
      MODECO = MODECI
      GOTO 201
C
C     TRACE DE Y=TEMPERATURE(t,x) AVEC ZONES DE COULEUR
C     ================================================================
 120  MODECI = MODECO
      MODECO = 5
      CALL TR1DTEY( KNOMOB, MODECO, NDIMLI,
     %              NBTYEL, MNNPEF, MNXYZP, NDPGST,
     %              NCAS0,  NCAS1,  NTDL,   RMCN(MNTEMP+WECTEU),
     %              TMIN,   NOEMIN, NCAMIN, TMAX,  NOEMAX, NCAMAX,
     %              RMCN(MNTIMP) )
      MODECO = MODECI
      GOTO 201
C
C     TRACE DES ISOTHERMES (LIGNES EN 2D ET SURFACES EN 3D ET 3-CUBES)
C     ================================================================
 210  CALL TRISOT( NDIMLI, KNOMOB, MODECO,
     %             NBTYEL, MNNPEF, MNXYZP, NDPGST,
     %             NCAS0,  NCAS1,  NTDL,
     %             0,      RMCN(MNTEMP+WECTEU),  dptemp,
     %             TMIN,   NOEMIN, NCAMIN, TMAX, NOEMAX, NCAMAX,
     %             RMCN(MNTIMP) )
      GOTO 201
C
C     TRACE DES ZONES DE COULEURS ISOTHERMES
C     ======================================
 220  CALL TRZONT( NOPROJ, NDIMLI, KNOMOB, MODECO,
     %             NBTYEL, MNNPEF, MNXYZP, MNXYZN, NDPGST,
     %             NCAS0,  NCAS1,  NTDL,
     %             0,      RMCN(MNTEMP+WECTEU),   dptemp,
     %             TMIN,   NOEMIN, NCAMIN, TMAX,  NOEMAX, NCAMAX,
     %             RMCN(MNTIMP) )
      GOTO 201
C
C     TRACE DES ZONES DE COULEURS ISOTHERMES PAR SECTIONS X ou Y ou Z=CTE
C     ===================================================================
 230  IF( NDIMLI .LE. 2 ) GOTO 220
      CALL TRPLSE( 0,      KNOMOB, MODECO,
     %             NBTYEL, MNNPEF, MNXYZP, NDPGST,
     %             NCAS0,  NCAS1,  NTDL,
     %             0,      RMCN(MNTEMP+WECTEU),   dptemp,
     %             TMIN,   NOEMIN, NCAMIN, TMAX,  NOEMAX, NCAMAX,
     %             RMCN(MNTIMP) )
      GOTO 201
C
C     TRACE DES PROFILS DE COULEURS PAR SECTIONS X ou Y ou Z=CTE
C     ==========================================================
 240  IF( NDIMLI .LE. 2 ) GOTO 220
      CALL TRPLSE( 1,      KNOMOB, MODECO,
     %             NBTYEL, MNNPEF, MNXYZP, NDPGST,
     %             NCAS0,  NCAS1,  NTDL,
     %             0,      RMCN(MNTEMP+WECTEU),   dptemp,
     %             TMIN,   NOEMIN, NCAMIN, TMAX,  NOEMAX, NCAMAX,
     %             RMCN(MNTIMP) )
      GOTO 201
C
C     TRACE DE LA TEMPERATURE LE LONG D'UNE DROITE DEFINIE PAR 2 POINTS
C     =================================================================
 250  IF( NDIMLI .LE. 2 ) GOTO 201
      CALL TRLLDR( NOPROJ, NDIMLI, KNOMOB, MODECO,
     %             NBTYEL, MNTOPO, MNNPEF, MNXYZP, MNXYZN, NDPGST,
     %             NCAS0,  NCAS1,  NTDL,
     %             0,      RMCN(MNTEMP+WECTEU),   dptemp,
     %             TMIN,   NOEMIN, NCAMIN, TMAX,  NOEMAX, NCAMAX,
     %             RMCN(MNTIMP) )
      GOTO 201
C
C     TRACE EN 2D SURFACE(X,Y,TEMPERATURE(X,Y))
C     =========================================
 270  IF( NDIMLI .NE. 2 ) GOTO 10
C     ICI PAS DE TEMPERATURE_EXACTE A PRENDRE EN COMPTE
      IF( MNVALP .LE. 0 ) THEN
         MODE = 1
      ELSE
         MODE = 2
      ENDIF
      CALL TRZTXY( 0,      NDIMLI, KNOMOB, MODE,
     %             NBTYEL, MNNPEF, MNXYZP,
     %             NCAS0,  NCAS1,  NTDL,
     %             0,      RMCN(MNTEMP+WECTEU),   dptemp,
     %             TMIN,   NOEMIN, NCAMIN, TMAX,  NOEMAX, NCAMAX,
     %             RMCN(MNTIMP) )
C     RETOUR A SA VALEUR INITIALE
      NDIMLI = 2
      GOTO 201
C
C     TRACE DE L'ERREUR ABSOLUE EN CHAQUE NOEUD DU MAILLAGE
C     =====================================================
 280  NOFOTI = 0
      IF( NTYSOL .LT. 4 ) THEN
C        EXISTENCE OU NON DE LA FONCTION 'TEMPERATURE_EXACTE'
         NOFOTI = NOFOTEEX()
C        NOFOTI>0 SI CETTE FONCTION EXISTE
         IF( NOFOTI .LE. 0 ) THEN
            NBLGRC(NRERR) = 1
            IF( LANGAG .EQ. 0 ) THEN
             KERR(1) = 'FONCTION TEMPERATURE_EXACTE(t,x,y,z) NON DONNEE'
            ELSE
               KERR(1) = 'FUNCTION EXACT_TEMPERATURE(t,x,y,z) NOT GIVEN'
            ENDIF
            CALL LEREUR
            GOTO 201
         ENDIF
      ELSE IF( NTYSOL .EQ. 4 .OR. NTYSOL .EQ. 5 ) THEN
C        EXISTENCE OU NON DE LA FONCTION ONDE NLSE 'PARTIE_REELLE_EXACTE'
         NOFOTI = NOFOPREX()
C        NOFOTI>0 SI CETTE FONCTION EXISTE
         IF( NOFOTI .LE. 0 ) THEN
            NBLGRC(NRERR) = 1
            IF( LANGAG .EQ. 0 ) THEN
           KERR(1) = 'FONCTION PARTIE_REELLE_EXACTE(t,x,y,z) NON DONNEE'
            ELSE
               KERR(1) = 'FUNCTION EXACT_REAL_PART(t,x,y,z) NOT GIVEN'
            ENDIF
            CALL LEREUR
            GOTO 201
         ENDIF
      ELSE IF( NTYSOL .EQ. 4 .OR. NTYSOL .EQ. 6 ) THEN
C        EXISTENCE OU NON DE LA FONCTION ONDE NLSE 'PARTIE_IMAGINAIRE_EXACTE'
         NOFOTI = NOFOPIEX()
C        NOFOTI>0 SI CETTE FONCTION EXISTE
         IF( NOFOTI .LE. 0 ) THEN
            NBLGRC(NRERR) = 1
            IF( LANGAG .EQ. 0 ) THEN
         KERR(1)='FONCTION PARTIE_IMAGINAIRE_EXACTE(t,x,y,z) NON DONNEE'
            ELSE
            KERR(1)='FUNCTION EXACT_IMAGINARY_PART(t,x,y,z) NOT GIVEN'
            ENDIF
            CALL LEREUR
            GOTO 201
         ENDIF
      ENDIF
      IF( NOFOTI .EQ. 0 ) GOTO 201
C
      MODECI = MODECO
      IF( MODECO .NE. 8 ) THEN
         MODECR = 4
      ELSE
         MODECR = 8
      ENDIF
C
C     CALCUL DES ERREURS MIN ET MAX, DE LEURS NOEUDS ET CAS
C     -----------------------------------------------------
      CALL VERREUR( NBCOOR, NOFOTI, MODECR, MNXYZP,
     %              NDSM,   NTDL,   RMCN(MNTEMP+WECTEU), RMCN(MNTIMP),
     %              ERRMIN, NERMIN, NCRMIN, ERRMAX, NERMAX, NCRMAX,
     %              MOERRE, MNERRE )
      IF( MNERRE .LE. 0 ) GOTO 201
C
C     TRACE DE L'ERREUR SUR LE MAILLAGE SELON LA DIMENSION DE L'ESPACE
C     ----------------------------------------------------------------
      MODECR = 4
      IF( NDIMLI .EQ. 1 ) THEN
C
C        TRACE EN 1D DE L'ERREUR AUX NOEUDS
         CALL TR1DTER( KNOMOB, MODECR, NDIMLI,
     %                 NBTYEL, MNNPEF, MNXYZP, NDPGST,
     %                 NDSM,   NTDL, RMCN(MNTIMP), RMCN(MNERRE+WECTEU),
     %                 ERRMIN, ERRMAX )
C
      ELSE IF( NDIMLI .EQ. 2 ) THEN
C
C        TRACE EN 2D DE L'ERREUR SUR LE MAILLAGE
 281     CALL LIMTCL( 'tracerr2', NMTCLE )
         IF( NMTCLE .LE. 0 ) GOTO 201
C
         GOTO ( 282, 284, 281, 285, 281, 286 ), NMTCLE
C
C        TRACE DES LIGNES ISOERREURS
 282     CALL TRISOT( NDIMLI, KNOMOB, MODECR,
     %                NBTYEL, MNNPEF, MNXYZP, NDPGST,
     %                NCAS0,  NCAS1,  NTDL,
     %                0,      RMCN(MNERRE+WECTEU), dptemp,
     %                ERRMIN, NERMIN, NCRMIN, ERRMAX, NERMAX, NCRMAX,
     %                RMCN(MNTIMP) )
         GOTO 281
C
C        TRACE DES ZONES DE COULEURS ISOERREURS
 284     CALL TRZONT( 0,      NDIMLI, KNOMOB, MODECR,
     %                NBTYEL, MNNPEF, MNXYZP, MNXYZN, NDPGST,
     %                NCAS0,  NCAS1,  NTDL,
     %                0,      RMCN(MNERRE+WECTEU), dptemp,
     %                ERRMIN, NERMIN, NCRMIN, ERRMAX, NERMAX, NCRMAX,
     %                RMCN(MNTIMP) )
         GOTO 281
C
C        TRACE DU PROFIL Z=ERREUR(X,Y) ou SURFACE(X,Y,ERREUR(X,Y))
 285     CALL TRZTXY( NOFOTI, NDIMLI, KNOMOB, MODECR,
     %                NBTYEL, MNNPEF, MNXYZP,
     %                NCAS0,  NCAS1,  NTDL,
     %                0,      RMCN(MNERRE+WECTEU), dptemp,
     %                ERRMIN, NERMIN, NCRMIN, ERRMAX, NERMAX, NCRMAX,
     %                RMCN(MNTIMP) )
         GOTO 281
C
C        TRACE de la COURBE ERREUR(Temps)
 286     CALL TRNLSERR( NTLXOB )
         GOTO 281
C
      ELSE
C
C        TRACE EN 3D DE L'ERREUR
         CALL TRER3D( NDIMLI, KNOMOB, MODECR,
     %                NBTYEL, MNTOPO, MNNPEF, MNXYZP, MNXYZN, NDPGST,
     %                NCAS0,  NCAS1,  NTDL, RMCN(MNERRE+WECTEU), dptemp,
     %                ERRMIN, NERMIN, NCRMIN, ERRMAX, NERMAX, NCRMAX,
     %                RMCN(MNTIMP) )
C
      ENDIF
C
      MODECO = MODECI
      GOTO 201
C
C     AFFICHAGE DES TEMPERATURES
C     ==========================
 290  CALL AFTEMP( NTDL,   NCAS,   MNXYZN,
     %             NTDL,   NDSM,   RMCN(MNTEMP+WECTEU),
     %             TECMOY, TECMIN, NOTMIN, TECMAX, NOTMAX,
     %             NOFOTI, TEXMIN, TEXMAX )
      IF( NMTCL .EQ. 7 ) GOTO 10
      GOTO 201
C
C     TRACE DU GRADIENT DE LA TEMPERATURE
C     ===================================
 300  CALL TRGRAD( NOPROJ, NCAS0, NCAS1, NDIMLI, KNOMOB, NTLXOB, MODECO,
     %             NBTYEL, MNNPEF, MNXYZP, NDPGST, RMCN(MNTIMP))
      GOTO 10
C
C     TRACE DU FLUX NORMAL DE LA TEMPERATURE
C     ======================================
 400  CALL TRFLUX( NOPROJ, NCAS0, NCAS1, NDIMLI, KNOMOB, NTLXOB, MODECO,
     %             NBTYEL, MNNPEF, MNXYZP, NDPGST, RMCN(MNTIMP))
      GOTO 10
C
C     TRACE DES ESTIMATEURS D'ERREUR EN 2D ET EN TEMPERATURE SEULEMENT
C     ================================================================
 500  IF( NDIMLI .EQ. 2 .AND. MODECO .EQ. 1 ) THEN
         CALL TRERTH( NCAS1,  NDIMLI, NTLXOB,
     %                NBTYEL, MNNPEF, MNXYZP, NDPGST )
      ENDIF
      GOTO 10
C
C     AFFICHAGE DU GRADIENT DE LA TEMPERATURE DE TOUS LES EF
C     ======================================================
 800  CALL AFGREF( KNOMOB, NTLXOB, NCAS0, NCAS1, MODECO )
      GOTO 10
C
C     AFFICHAGE DU FLUX NORMAL DE LA TEMPERATURE
C     AUX POINTS DES INTERFACES ENTRE EF
C     ==========================================
 900  CALL AFFLEF( KNOMOB, NTLXOB, NCAS0, NCAS1, MODECO )
      GOTO 10
C
C     AFFICHAGE DU FLUX NORMAL DE LA TEMPERATURE
C     SUR LES PLS FRONTIERE DE L'OBJET
C     ==========================================
 1000 CALL AFFLUXFR( KNOMOB )
      GOTO 10
C
C     FIN DE L'EXECUTION
C     ==================
C     RETOUR AUX PARAMETRES INITIAUX
 9000 CALL XVEPAISSEUR( 1 )
      CALL XVTYPETRAIT( LIGCON )
      IF( MNNPEF .GT. 0 ) CALL TNMCDS( 'ENTIER', 2*MXTYEL, MNNPEF )
      IF( NBCOOR .EQ. 6 .AND. MNXYZ3Q1C .GT. 0 ) THEN
         CALL TNMCDS( 'MOTS', MOXYZ3Q1C, MNXYZ3Q1C )
      ENDIF
      IF( MNTIMP .GT. 0 ) CALL TNMCDS( 'REEL',   NDSM,   MNTIMP )
      IF( MNMOWA .GT. 0 ) CALL TNMCDS( 'REEL2',  MOMOWA, MNMOWA )
      IF( MNERRE .GT. 0 ) CALL TNMCDS( 'ENTIER', MOERRE, MNERRE )
CC
      RETURN
      END
