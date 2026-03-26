      SUBROUTINE STOKESTA( KNOMOB, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LES VITESSES ET PRESSIONS DANS UN FLUIDE
C -----    SOLUTIONS DE L'EQUATION de STOKES STATIONNAIRE

C ENTREE :
C --------
C KNOMOB : NOM DE L'OBJET FLUIDE STATIONNAIRE A TRAITER

C SORTIE :
C --------
C IERR   : =0 SI PAS D'ERREUR D'EXECUTION
C          >0 SI ERREUR DECTECTEE
C          20 SI PAS ASSEZ DE PLACE MEMOIRE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : BENHAMADOUCHE BOYER DEA ANALYSE NUMERIQUE PARIS  JANVIER 2000
C MODIFS : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 2000
C MODIFS : FREDERIC LEGOLL DEA ANALYSE NUMERIQUE PARIS      FEVRIER 2001
C MODIFS : YANN  REVALOR   DEA ANALYSE NUMERIQUE PARIS      FEVRIER 2003
C MODIFS : ALAIN PERRONNET Laboratoire J-L.LIONS UPMC Paris    JUIN 2007
C MODIFS : ALAIN PERRONNET Laboratoire J-L.LIONS UPMC Paris OCTOBRE 2007
C MODIFS : ALAIN PERRONNET LJLL UPMC & St Pierre du Perray Decembre 2008
C MODIFS : ALAIN PERRONNET LJLL UPMC & St Pierre du Perray     Mars 2009
C MODIFS : ALAIN PERRONNET  Saint Pierre du Perray             Mars 2021
C MODIFS : ALAIN PERRONNET  Saint Pierre du Perray          Fevrier 2022
C23456---------------------------------------------------------------012
      PARAMETER      (MXTYEL=7)
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/donflu.inc"
      include"./incl/donele.inc"
      include"./incl/a_objet__topologie.inc"
      include"./incl/a_objet__definition.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___xyznoeud.inc"
      include"./incl/a___xyzpoint.inc"
      include"./incl/a___npef.inc"
      include"./incl/a___force.inc"
      include"./incl/a___vecteur.inc"
      include"./incl/a___profil.inc"
      include"./incl/a___morse.inc"
      include"./incl/msvaau.inc"
      include"./incl/ponoel.inc"
      include"./incl/ctemps.inc"
      include"./incl/cnonlin.inc"

      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1),DMCN(1))
C
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      EXTERNAL          ETTAEL
C
      CHARACTER*(*)     KNOMOB
      CHARACTER*4       NOMELE(2)
      CHARACTER*80      KNOMFIC
      LOGICAL           COMP
C
      REAL              COORDP(30)
      INTEGER           NUMIOB(4), NUMAOB(4), MNDOEL(4),  MXDOEL(4)
      INTEGER           NOOBVC,    NOOBSF(6), NOOBLA(12), NOOBPS(8)
      INTEGER           NONOEF(10), NOGLDL(34)
      DOUBLE PRECISION  DINFO, DCPU, DFABG, DFACTO, DVITPR
ccc      DOUBLE PRECISION  DIAMIN, DIAMAX
      DOUBLE PRECISION  VITMIN, VITMAX, VITMOY, PREMIN, PREMAX, PREMOY,
     %                  TEMPER(1)
      DOUBLE PRECISION  RELMIN
      DATA              RELMIN / -1D28 /

C     QUELQUES INITIALISATIONS
      DCPU   = DINFO( 'DELTA CPU' )
      TEMPS  = 0.0
      IDEPL  = 0
      IERR   = 0
C     LE NOMBRE DE MOTS D'UNE VARIABLE REEL2
      MOREE2 = MOTVAR(6)
C
C     PROTECTION DES ADRESSES POUR EVITER DES PROBLEMES LORS
C     DE LA DESTRUCTION DES TABLEAUX
      MNAUX  = 0
      MNNPEF = 0
      MNNPEF = 0
      MNTPOB = 0
      MNTAUX = 0
      MNTAEL = 0
      MNNODL = 0
      MNIP   = 0
      MNFLUI = 0
      MNFORC = 0
      MNX    = 0
      MNPTDG = 0
      MNNFNX = 0
      MNVFNX = 0
      DO 5 I=1,4
         MNDOEL(I) = 0
         MXDOEL(I) = 0
 5    CONTINUE
      NBTYEL = 0
      MOAUX  = 0
      MOTAEL = 0
      NBDLMX = 0
      NTDL   = 0
      NTDLF  = 0
      NBRFNX = 0
      NBFNFX = 0
      NBDLFX = 0
      MOFLTO = 0
      MOFLPT = 0
      MONFNX = 0
      MOVFNX = 0
      MNLPCO = 0
      MNVG   = 0
      MNNDDL = 0
      MNSMGL = 0
      DCPU   = 0D0
      DFABG  = 0D0
      DFACTO = 0D0
      DVITPR = 0D0
      MNNDLFX= 0
      MNVDLFX= 0
      MNDPDP = 0
      NBLGAU = 0
      MNLGAU = 0
      MNLPLC = 0
      MNLPCC = 0
      MNLPLU = 0
      MNVGC  = 0
      NIVMAX = 10
      MNBDIR = 0
      MNNDFX = 0
      MNPRPR = 0
      MNAUX1 = 0
      MNAUX2 = 0
      MNAUX3 = 0
      MNAUX4 = 0
      MNAUX5 = 0
      MNAUX6 = 0
      MNAUX7 = 0
      MNAUX8 = 0
      MNDIR  = 0
      MNDAD  = 0
      MNBET  = 0
      MNADIR = 0
      MOTSGC = 0
      NBNOTE = 0
      NTDLTE = 0
C
      DO I=1,10
        NONOEF(I) = 0
      ENDDO
C
C     VERIFICATION DU NOM_DE_L'OBJET
C     ==============================
C     RECHERCHE DU NOM DE L'OBJET DANS LE LEXIQUE DES OBJETS
 100  CALL LXLXOU( NTOBJE, KNOMOB, NTLXOB, MNLXOB )
C     S'IL N'EXISTE PAS RETOUR A LA DEMANDE DU NOM DE L'OBJET
      IF( NTLXOB .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'ERREUR: OBJET INCONNU ' // KNOMOB
         ELSE
            KERR(1) = 'ERROR: UNKNOWN OBJECT ' // KNOMOB
         ENDIF
         CALL LEREUR
         CALL LXIM( NTOBJE )
         IERR = 1
         RETURN
      ENDIF
C
C     RECHERCHE DU TABLEAU DEFINITION
      CALL LXTSOU( NTLXOB, 'DEFINITION', NTDFOB, MNDFOB )
C     S'IL N'EXISTE PAS RETOUR A LA DEMANDE DU NOM DE L'OBJET
      IF( NTDFOB .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) ='ERREUR: DEFINITION INCONNUE de l''OBJET ' //KNOMOB
         ELSE
            KERR(1) ='ERROR: UNKNOWN DEFINITION for the OBJECT '//KNOMOB
         ENDIF
         CALL LEREUR
         GOTO 100
      ENDIF
C
C     L'ANCIEN HISTORIQUE EST EFFACE
      CALL RECTEF( NRHIST )
      NBLGRC(NRHIST) = 1
      IF( LANGAG .EQ. 0 ) THEN
         KHIST(1) = 'OBJET: ' // KNOMOB
      ELSE
         KHIST(1) = 'OBJECT: ' // KNOMOB
      ENDIF
      CALL LHISTO
C
 10   IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10009) KNOMOB
      ELSE
         WRITE(IMPRIM,20009) KNOMOB
      ENDIF
10009 FORMAT(/,1X,100('=')/,
     %' stokesta: ECOULEMENT de STOKES d''un FLUIDE STATIONNAIRE dans ',
     % A,/1X,100('='))
20009 FORMAT(/,1X,100('=')/,
     %' stokesta: STUDY of a STOKES STEADY FLUID FLOW in ',A,/,
     %  1X,100('='))
C
C     RECHERCHE DES TABLEAUX XYZNOEUD XYZPOINT ASSOCIES A L'OBJET
C     -----------------------------------------------------------
      CALL MIMAOB( 1,      NTLXOB, MXDOFL, NTTOPO, MNTOPO,
     %             NTXYZP, MNXYZP, NTXYZN, MNXYZN,
     %             MXTYEL, NBTYEL, MTNPEF, MNNPEF,
     %             NUMIOB, NUMAOB,
     %             NDPGST, NBOBIN, MNOBIN, NBOBCL, MNOBCL,
     %             MXDOEL, MNDOEL, IERR )
      IF( IERR .NE. 0 ) RETURN
C
C     UN SEUL TYPE D'EF PERMIS
      IF( NBTYEL .LE. 0 .OR. NBTYEL .GT. 1 ) THEN
         NBLGRC(NRERR) = 2
         WRITE(KERR(MXLGER)(1:4),'(I4)') NBTYEL
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'STOKESTA:' // KERR(MXLGER)(1:4)
     %             //' TYPES d''ELEMENTS FINIS'
            KERR(2) = '1 SEUL TYPE D''ELEMENT FINI PERMIS'
         ELSE
            KERR(1) = 'STOKESTA:' // KERR(MXLGER)(1:4)
     %             //' TYPES of FINITE ELEMENTS'
            KERR(2) = 'ONLY ONE TYPE OF FINITE ELEMENT ALLOWED'
         ENDIF
         CALL LEREUR
         IERR = 3
         RETURN
      ENDIF
C
C     CHOIX DE LA METHODE DE RESOLUTION DU SYSTEME LINEAIRE
C     -----------------------------------------------------
C     ENTREE du CHOIX de la METHODE de RESOLUTION des SYSTEMES Ax=b
C     ( => le MODE de STOCKAGE de la MATRICE GLOBALE )
      CALL CHMERESO( NAVSTO, NORESO, NIVOGC, IERR )
C     NORESO : CODE DE RESOLUTION DES SYSTEMES LINEAIRES
C              1 FACTORISATION CROUT avec stockage PROFIL
C              2 GRADIENT CONJUGUE avec CROUT INCOMPLET
C              3 CGS DOUBLE GRADIENT CONJUGUE ACCELERE ou
C                CONJUGATE GRADIENT SQUARED METHOD sur MATRICES MORSES
C              4 BICG Stab ou Bi-CONJUGATE GRADIENT STABILISED 
C                sur MATRICES MORSES CONDENSEES
C     NIVOGC : SI NORESO=2 alors
C              0 ou 1 ou 2 NIVEAU DU STOCKAGE DE LA MATRICE DE
C                          PRECONDITIONNEMENT DU GC

C*********************************************************************
C     LE GRADIENT CONJUGUE NE CONVERGE PAS
C     CAR LA MATRICE N'EST PAS DEFINIE POSITIVE
C     PAR CONTRE LA FACTORISATION DE CROUT L D tL
C     CONVIENT ET EST IMPOSEE
      IF( NORESO .NE. 1 ) THEN
         NORESO = 1
         IF( LANGAG .EQ. 0 ) THEN
         PRINT*,'stokesta: Le GRADIENT CONJUGUE NE CONVERGE PAS'
         PRINT*,'stokesta: La FACTORISATION L D tL de CROUT est IMPOSEE'
         ELSE
         PRINT*,'stokesta: The CONJUGATE GRADIENT DOES NOT CONVERGE'
         PRINT*,'stokesta: The CROUT FACTORISATION L D tL IS IMP0SED'
         ENDIF
      ENDIF
C*********************************************************************

      IF( NORESO .EQ. 1 ) THEN
         NIVOGC = 0
         IF( LANGAG .EQ. 0 ) THEN
           PRINT*,' RESOLUTION par FACTORISATION COMPLETE de CROUT LDtL'
         ELSE
            PRINT*,' SOLUTION by COMPLETE CROUT FACTORIZATION L D tL'
         ENDIF
      ELSE IF( NORESO .GE. 2 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            PRINT*,'RESOLUTION par GRADIENT CONJUGUE PRECPNDITIONNE'
            PRINT*,'AVEC UN NIVEAU',NIVOGC,' DE FACTORISATION'
         ELSE
            PRINT*,'SOLUTION by PRECONDITIONED CONJUGATE GRADIENT'
            PRINT*,'with a LEVEL',NIVOGC,' of FACTORIZATION'
         ENDIF
      ELSE
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'METHODE de RESOLUTION Ax=b NON PROGRAMMEE'
         ELSE
            KERR(1) = 'Type of SOLVER Ax=b NOT PROGRAMMED'
         ENDIF
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF

      IF( NBOBIN .GT. 1 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'RESOLUTION POUR 1 FLUIDE SEULEMENT'
         ELSE
            KERR(1) = 'ONLY ONE FLUID IS PROGRAMMED'
         ENDIF
         CALL LEREUR
         IF( INTERA .GE. 3 ) GOTO 10
         IERR = 1
         RETURN
      ENDIF

C     RETROUVER LES ADRESSES MCN DES DONNEES FLUIDE
C     DES   SV "OBJETS INTERNES"    DE L'OBJET
C     DES PLS  "OBJETS AUX LIMITES" DE L'OBJET
C     =============================================
      CALL FLUDON( NUMIOB, NBOBIN, MNOBIN, NBOBCL, MNOBCL, MNDOEL,
     %             IEMASS, IEVISC, IECOPR, IEFOIN, IEVTIN, IEPRIN,
     %             IEFOCL, IEFOPO, IEBLVI, IEBLPR, IEVIAN, IECOBO,
     %             IERR )
C
      IF( IEVISC .NE. NBOBIN ) THEN
          NBLGRC(NRERR) = 2
          IF( LANGAG .EQ. 0 ) THEN
             KERR(1) = 'Une SURFACE ou VOLUME SANS VISCOSITE'
             KERR(2) = '=> PROBLEME SANS SOLUTION'
          ELSE
             KERR(1) = '1 SURFACE or VOLUME WITHOUT VISCOSITY'
             KERR(2) = '=> PROBLEM WITHOUT SOLUTION'
          ENDIF
          CALL LEREUR
          IERR = 2
      ENDIF
C
C     ON N'AUTORISE PAS UN PB SANS CL POUR LA VITESSE
      IF ( IEBLVI .LE. 0 ) THEN
          NBLGRC(NRERR) = 2
          IF( LANGAG .EQ. 0 ) THEN
             KERR(1) = 'OBJET SANS CONDITION sur la VITESSE'
             KERR(2) = '=> PROBLEME SANS SOLUTION'
          ELSE
             KERR(1) = 'OBJECT WITHOUT BOUNDARY CONDITION on VELOCITY'
             KERR(2) = '=> PROBLEM WITHOUT SOLUTION'
          ENDIF
          CALL LEREUR
          IERR = 2
      ENDIF

cccC     UN PB SANS CL SUR LA PRESSION EST AVERTI
ccc      IF ( IEBLPR .LE. 0 ) THEN
ccc         WRITE(IMPRIM,*)
ccc         NBLGRC(NRERR) = 2
ccc         IF( LANGAG .EQ. 0 ) THEN
ccc             KERR(1) = 'PAS DE CONDITION AUX LIMITES sur la PRESSION'
ccc             KERR(2) = '=> PRESSION DEFINIE A UNE CONSTANTE PRES'
ccc         ELSE
ccc             KERR(1) = 'NO BOUNDARY CONDITION on PRESSURE'
ccc             KERR(2) = '=> PRESSURE DEFINED MODULO A CONSTANT'
ccc         ENDIF
ccc         CALL LERESU
ccc      ENDIF

      IF( IERR .NE. 0 ) RETURN

C     INITIALISATIONS DE VARIABLES ET AFFICHAGES
C     ==========================================
C     LE NOMBRE DE POINTS DU MAILLAGE DE L'OBJET
      NBPOI = MCN( MNXYZP + WNBPOI )

C     NDIM LA DIMENSION 2 OU 3 DE L'ESPACE DES COORDONNEES
      CALL DIMCOO( NBPOI, MCN(MNXYZP+WYZPOI), NDIM )

C     MNELE : ADRESSE DU TABLEAU NPEF"TYPE EF
      MNELE = MCN( MNNPEF )

C     NOMBRE D'ELEMENTS FINIS DE CE TYPE
      NBELEM = MCN( MNELE + WBELEM )

C     LE NUMERO DU TYPE DE L'ELEMENT FINI
      NUTYEL = MCN( MNELE + WUTYEL )

C     LE TYPE D'EF DE FLUIDE
      IF( NUTYEL .EQ. 13 ) THEN

C        TRIANGLE 2P1D  =>  EF BREZZI-FORTIN 2D
C        --------------------------------------
C        NOMBRE DE SOMMETS DU MAILLAGE ET NOMBRE DE NOEUDS VITESSE
         NBSOM = NBPOI
 3       IF( NBSOM .GT. NBELEM ) THEN
            NBSOM = NBSOM - NBELEM
            GOTO 3
         ENDIF
C        NOMBRE DE NOEUDS DU MAILLAGE
         NBNOEU = NBSOM
         NBPOI  = NBSOM
C        NOMBRE DE NOEUDS SUPPORT DE LA VITESSE
         NBNOVI = NBNOEU + NBELEM
C        NOMBRE TOTAL DE DEGRES DE LIBERTE EF TYPE BREZZI-FORTIN
C        CHAQUE SOMMET SUPPORTE V1 V2 P SANS LES VITESSE
C        AU BARYCENTRE DU TRIANGLE SUPPRIMEES PAR GAUSS
         NTDL  = 3 * NBSOM
C        ET AVEC LA VITESSE AU BARYCENTRE DU TRIANGLE
         NTDLF = NTDL + 2 * NBELEM
C        NOMBRE DE DEGRES DE LIBERTE ELEMENTAIRES DE L'EF (P1*(V+P)+BL*V)
         NBDLMX = 11

      ELSE IF( NUTYEL .EQ. 19 ) THEN

C        TETRAEDRE 3P1D =>  EF BREZZI-FORTIN 3D
C        --------------------------------------
C        NOMBRE DE SOMMETS DU MAILLAGE ET NOMBRE DE NOEUDS VITESSE
         NBSOM = NBPOI
 4       IF( NBSOM .GT. NBELEM ) THEN
            NBSOM = NBSOM - NBELEM
            GOTO 4
         ENDIF
C        NOMBRE DE NOEUDS DU MAILLAGE
         NBNOEU = NBSOM
C        NOMBRE DE NOEUDS SUPPORT DE LA VITESSE
         NBNOVI = NBNOEU + NBELEM
C        CHAQUE NOEUD=SOMMET SUPPORTE V1 V2 V3 P
C        SANS LES VITESSES MOYENNES DU TETRAEDRE SUPPRIMEES PAR GAUSS
         NTDL  = 4 * NBSOM
C        ET AVEC LES VITESSES MOYENNES DU TETRAEDRE
         NTDLF = NTDL + 3 * NBELEM
C        NOMBRE DE DEGRES DE LIBERTE ELEMENTAIRES DE L'EF (P1*(V+P)+BL*V)
         NBDLMX = 19
C
      ELSE IF( (NUTYEL.EQ.15) .OR. (NUTYEL.EQ.20) ) THEN
C
C        TRIANGLE  TAYLOR-HOOD avec TRIA 2P2C ou
C        TETRAEDRE TAYLOR-HOOD avec TETR 3P2C
C        ---------------------------------------
C        NOMBRE DE NOEUDS DU MAILLAGE
         NBNOEU = MCN( MNXYZN + WNBNOE )
C        NOMBRE DE NOEUDS SUPPORT DE LA VITESSE
         NBNOVI = NBNOEU
C
C        DETERMINATION DU NOMBRE NBSOM DE SOMMETS DE LA TRIANGULATION
         MNTAUZ=0
         CALL TNMCDC('ENTIER', NBNOEU, MNTAUZ)
         CALL AZEROI( NBNOEU,  MCN(MNTAUZ) )
         DO NUELEM=1,NBELEM
            CALL EFNOEU( MNELE, NUELEM, NBNOEF, NONOEF )
            DO I=1,3
               MCN(MNTAUZ+NONOEF(I)-1) = 1
            ENDDO
         ENDDO
         NBSOM = 0
         DO I=1,NBNOEU
            IF( MCN(MNTAUZ+I-1) .GT. 0 ) THEN
               NBSOM = NBSOM + 1
            ENDIF
         ENDDO
         CALL TNMCDS( 'ENTIER', NBNOEU, MNTAUZ )
C
C        LE NOMBRE TOTAL DE DEGRES DE LIBERTE EF TYPE TAYLOR-HOOD
C        CHAQUE SOMMET SUPPORTE V1 V2 P et le MILIEU DES ARETES SUPPORTE V1 V2
         NTDL  = NBNOEU * NDIM + NBSOM
         NTDLF = NTDL
C        NOMBRE DE DEGRES DE LIBERTE ELEMENTAIRES DE L'EF (EN 2D)
         NBDLMX = 15
C
C        RECUPERATION DES TABLEAUX POBA NECESSAIRES A LA
C        CONSTRUCTION DES TABLEAUX ELEMENTAIRES DU TETRAEDRE P2
         IF( NUTYEL .EQ. 20 ) THEN
            CALL TAPOBA( NBTYEL, MNNPEF, ETTAEL,
     %                   MNTPOB, NBDLMX, MOAUX, NBTTEF, NOAXIS, I, IERR)
            IF( IERR .NE. 0 ) GOTO 9900
C           NOMBRE DE DL DE L'EF
            NBDLMX = 34
         ENDIF
C
      ELSE
C
C        EF NI BREZZI-FORTIN, NI TAYLOR-HOOD => ERREUR
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'EF ni BREZZI-FORTIN ni TAYLOR-HOOD'
         ELSE
            KERR(1) = 'FE neither BREZZI-FORTIN nor TAYLOR-HOOD'
         ENDIF
         CALL LEREUR
         IERR = 1
         GOTO 9000
C
      ENDIF
C
C     QQ AFFICHAGES
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,8300) NDIM,NUTYEL,NBSOM,NBPOI,NBNOEU,NBELEM,NTDLF
      ELSE
         WRITE(IMPRIM,28300) NDIM,NUTYEL,NBSOM,NBPOI,NBNOEU,NBELEM,NTDLF
      ENDIF
 8300 FORMAT(/' ESPACE du PROBLEME R',I1,/,
     %        ' TYPE DE L''EF' ,T25,I7,/,
     %        ' NB DE SOMMETS ',T25,I7,/,
     %        ' NB DE POINTS ', T25,I7,/,
     %        ' NB DE NOEUDS ', T25,I7,/,
     %        ' NB D''EF ',     T25,I7,/,
     %        ' NB DL VITESSE-PRESSION',T25,I7)
28300 FORMAT(/' SPACE DIMENSION of PB ',T30,I7,/
     %        ' TYPE of FE',     T30,I7,/,
     %        ' NB of VERTICES ',T30,I7,/,
     %        ' NB of POINTS ',  T30,I7,/,
     %        ' NB of NODES ',   T30,I7,/,
     %        ' NB of FE ',      T30,I7,/,
     %        ' NB of DoF VELOCITY-PRESSURE',T30,I7)
C
C
C     =====================================================================
C     TRAITEMENT DU PROBLEME DE STOCKES PAR EF BREZZI-FORTIN ou TAYLOR-HOOD
C     =====================================================================
C
C     CONSTRUCTION DU TABLEAU NDDLNO(0:NBNOVI)
C     POINTEUR SUR LE DERNIER DL DE CHAQUE NOEUD VITESSE DU FLUIDE
C     ============================================================
      MNNDDL = 0
      CALL TNMCDC( 'ENTIER', 1+NBNOVI, MNNDDL )
      CALL AZEROI( 1+NBNOVI, MCN(MNNDDL) )
      CALL PTDLFL( MNELE,  NDIM, NBELEM, NBNOEU, NUTYEL,
     %             NBSOMT, NBNOVI, MCN(MNNDDL), IERR )
      IF( IERR .NE. 0 ) RETURN
C
C     CONSTRUCTION DU POINTEUR SUR LA DIAGONALE DE LA MATRICE GLOBALE
C     ===============================================================
      IF( NORESO .EQ. 1 ) THEN
C
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10015)
         ELSE
            WRITE(IMPRIM,20015)
         ENDIF
10015    FORMAT(' RESOLUTION Ax=b PAR CROUT PROFIL L D tL')
20015    FORMAT(' SOLUTION Ax=b by CROUT FACTORIZATION L D Lt' )
C
C        METHODE DE CROUT avec STOCKAGE PROFIL:
C        ======================================
C        CALCUL DU PROFIL DE LA MATRICE
C        (LA MATRICE PROFIL EST ICI SYMETRIQUE)
         NCODSV = 1
C
C        CREATION DU TABLEAU POINTEUR SUR DIAGONALE POUR LE FLUIDE
         CALL TNMCDC( 'ENTIER', 1+NTDL, MNPTDG )
         CALL AZEROI( 1+NTDL, MCN(MNPTDG) )
C
         IF( NUTYEL .EQ. 13 .OR. NUTYEL .EQ. 19 ) THEN
C
C           BREZZI-FORTIN TRIANGLE 2D ou TETRAEDRE 3D
C           CALCUL DU POINTEUR SUR CHAQUE COEFFICIENT DIAGONAL DE LA
C           MATRICE PROFIL GLOBALE (SANS LES BARYCENTRES DES EF)
C           MNELE : ADRESSE DU TABLEAU NPEF"TYPE EF
            MNELE   = MCN( MNNPEF )
C           MNUNDEL : ADRESSE DU TABLEAU NUNDEL DES EF DE BREZZI-FORTIN
            MNUNDEL = MNELE + WUNDEL
C           NBSTEF  = NOMBRE DE SOMMETS DE L'EF = NB DE DL DE CHAQUE SOMMET
            NBSTEF  = NDIM + 1
            print *,'STOKESTA BF: NBSOM=',NBSOM,' NBSTEF=',NBSTEF
            CALL DIAGBF( NBSTEF, NBSOM,  NBELEM, NBSTEF, MCN(MNUNDEL),
     %                   NCODSV, MCN(MNPTDG), IERR )
ccc            CALL PRPRMC( MNTOPO, MCN(MNNPEF), MNXYZP, NDIM+1, NCODSV,
ccc     %                   MCN(MNPTDG), IERR )
            IF( IERR .GT. 0 ) GOTO 9000
C
C           NBCOVG DONNE LE NOMBRE DE TERMES NON NULS
            NBCOVG = MCN( MNPTDG + NBSOM*NBSTEF )
C
         ELSE
C
C           TAYLOR-HOOD TRIANGLE 2D ou TETRAEDRE 3D
C           CREATION DU TABLEAU AUXILIAIRE POUR CALCULER LE POINTEUR
C           SUR DIAGONALE "CLASSIQUE " A UTILISER POUR NOTRE POINTEUR
C           SUR LA DIAGONALE (POUR L'ELEMENT FINI TAYLOR-HOOD P2 P1)
            CALL TNMCDC( 'ENTIER', NBNOEU+1, MNTAUX )
            CALL AZEROI( NBNOEU+1, MCN(MNTAUX) )

C           CALCUL DU POINTEUR SUR CHAQUE COEFFICIENT DIAGONAL DE LA
C           MATRICE PROFIL GLOBALE
            CALL PRPRMC( MNTOPO, MCN(MNNPEF), MNXYZN, 1, NCODSV,
     %                   MCN(MNTAUX), IERR )
            IF( IERR .GT. 0 ) GOTO 9001

C           A PARTIR DU POINTEUR SUR LA DIAGONALE DU PROFIL AUX NOEUDS
C           GENERER LE POINTEUR SUR LA DIAGONALE DES DEGRES DE LIBERTE
C           DE CHAQUE NOEUD DEFINI PAR NDDLNO POINTEUR SUR LE NO DU
C           DERNIER DL DE CHAQUE NOEUD
            CALL LPLINODL( NBNOEU, NTDL, MCN(MNNDDL),NCODSV,MCN(MNTAUX),
     %                     MCN(MNPTDG), NBCOVG, IERR )
cccC           NBCOVG DONNE LE NOMBRE DE TERMES NON NULS
ccc            NBCOVG = MCN( MNPTDG + NTDL )

C           DESTRUCTION DU TABLEAU AUXILIAIRE
 9001       CALL TNMCDS('ENTIER',NBNOEU+1,MNTAUX)
            IF( IERR .GT. 0 ) GOTO 9000

         ENDIF
C
C        DECLARATION INITIALISATION DE LA MATRICE PROFIL SYMETRIQUE
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10290) NBCOVG
         ELSE
            WRITE(IMPRIM,20290) NBCOVG
         ENDIF
10290    FORMAT(' MATRICE PROFIL DE',I15,' REELS DOUBLE PRECISION')
20290    FORMAT(' SKYLINE MATRIX of',I15,' DOUBLE PRECISION REALS')
C
C        TEST DE LA PLACE MEMOIRE DISPONIBLE
         CALL TNMCMX( 'REEL2', MAXVAR )
         LO = NBCOVG
         IF (MAXVAR.LT.LO) THEN
C            PLACE MEMOIRE INSUFFISANTE
             NBLGRC(NRERR) = 2
             IF( LANGAG .EQ. 0 ) THEN
                KERR(1) = 'ERREUR PLACE MEMOIRE INSUFFISANTE'
                KERR(2) = 'UTILISEZ LA METHODE DU GRADIENT CONJUGUE'
             ELSE
                KERR(1) = 'ERROR NOT ENOUGH MEMORY'
                KERR(2) = 'USE the CONJUGATE GRADIENT METHOD'
             ENDIF
             CALL LEREUR
             IERR = 20
             GOTO 9000
         ENDIF
C
C        DECLARATION DE LA MATRICE GLOBALE DE VISCOSITE
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10002)
         ELSE
            WRITE(IMPRIM,20002)
         ENDIF
10002    FORMAT(' CONSTRUCTION de la MATRICE de VISCOSITE [V]')
20002    FORMAT(' CONSTRUCTION of the VISCOSITY MATRIX [V]')

      ELSE IF( NORESO .GE. 2 ) THEN

C        METHODES DU GRADIENT CONJUGUE avec STOCKAGE MORSE
C        =================================================
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10016)
         ELSE
            WRITE(IMPRIM,20016)
         ENDIF
10016    FORMAT(/' RESOLUTION PAR GRADIENT CONJUGUE MATRICE CONDENSEE')
20016    FORMAT(/' RESOLUTION by SPARSE MATRIX CONJUGATE GRADIENT')

C        CALCUL DES TABLEAUX POINTEURS DE LA MATRICE MORSE PAR NDDLNO
C        (LA MATRICE EST ICI SYMETRIQUE)
         NCODSV = 1
         MNPTDG = 0
         MNLPCO = 0
         CALL PRGCFL0( MNTOPO, MCN(MNNPEF), MNXYZN,
     %                 NBNOEU, MCN(MNNDDL), NCODSV,
     %                 MNPTDG, MNLPCO, IERR)
         IF( IERR .NE. 0 ) THEN
            IERR = 20
            GOTO 9000
         ENDIF

C        DECLARATION INITIALISATION DE LA MATRICE MORSE SYMETRIQUE
         NBCOVG = MCN(MNPTDG+NTDL)
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10291) NBCOVG
         ELSE
            WRITE(IMPRIM,20291) NBCOVG
         ENDIF
10291    FORMAT(' MATRICE MORSE ',I15,' VARIABLES DOUBLE PRECISION')
20291    FORMAT(' SPARSE MATRIX of',I15,' DOUBLE PRECISION REALS')

      ENDIF

C     DECLARATION MATRICE GLOBALE DE VISCOSITE POUR LES 2 METHODES
      MNVG = 0
      CALL TNMCDC( 'REEL2', NBCOVG, MNVG )
      CALL AZEROD( NBCOVG, MCN(MNVG) )

C     ADRESSAGE INITIALISATION DU VECTEUR VITESSEPRESSION
C     ===================================================
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10003)
      ELSE
         WRITE(IMPRIM,20003)
      ENDIF
10003 FORMAT(' CONSTRUCTION DU VECTEUR {b} SECOND MEMBRE')
20003 FORMAT(' CONSTRUCTION of the SECOND MEMBER VECTOR {b}')

      CALL LXTSOU( NTLXOB, 'VECTEUR"VITESSEPRESSION', NTVIPE, MNVIPE )
      IF( NTVIPE .GT. 0 ) THEN
C        LE VECTEUR EST DETRUIT POUR ETRE REDECLARE
         CALL LXTSDS( NTLXOB, 'VECTEUR"VITESSEPRESSION' )
      ENDIF

      L = WECTEU + NTDLF * MOREE2
      CALL LXTNDC( NTLXOB, 'VECTEUR"VITESSEPRESSION', 'MOTS', L )
      CALL LXTSOU( NTLXOB, 'VECTEUR"VITESSEPRESSION', NTVIPE, MNVIPE )
C     ADRESSE DANS MCN DU TMS VITESSEPRESSION
      MNU = MNVIPE + WECTEU
      CALL AZEROD( NTDLF, MCN(MNU) )

C     ADRESSE DU SECOND MEMBRE GLOBAL
      MNSMGL = 0
      CALL TNMCDC( 'REEL2', NTDLF, MNSMGL )
      CALL AZEROD(  NTDLF, MCN(MNSMGL) )

C     DIMENSION MAXIMALE MATRICE ELEMENTAIRE SYMETRIQUE + SECOND MEMBRE
C     ELEMENTAIRE POUR LES METTRE L'UN DERRIERE L'AUTRE
      MOTAEL = NBDLMX * (NBDLMX+1) / 2 + NBDLMX
      MNTAEL = 0
      CALL TNMCDC( 'REEL2', MOTAEL, MNTAEL )

C     TABLEAU QUI DONNERA LE NUMERO GLOBAL D'UN DL DANS UN EF
      MNNODL = 0
      CALL TNMCDC( 'ENTIER', NBDLMX, MNNODL )
C
C     MNNDEL ADRESSE MCN DES NUMEROS NOEUDS DES EF
      MNNDEL = MNELE + WUNDEL
C     MNPGEL ADRESSE MCN DES NUMEROS POINTS GEOMETRIQUES DES EF
      MNPGEL = MNNDEL
ccc      Dessous mis en commentaire car pas de stockage du no des points
ccc      des EF de BREZZI-FORTIN
ccc      IF( NDPGST .GE. 2 .AND. NUTYEL .NE. 13 .AND. NUTYEL .NE. 19 ) THEN
ccc         MNPGEL = MNPGEL + NBELEM * MCN(MNELE+WBNDEL)
ccc      ENDIF

C     LES CARACTERISTIQUES DE L'ELEMENT FINI
      CALL ELNUNM( NUTYEL, NOMELE )
      CALL ELTYCA( NUTYEL )
C
C     LE NOMBRE DE DEGRES DE LIBERTE FIXES
      NBDLFX = 0
C
      IF( NUTYEL .EQ. 13 .OR. NUTYEL .EQ. 19 ) THEN
C
C        TRIANGLE ou TETRAEDRE BREZZI-FORTIN
C        ===================================
C        CONSTRUCTION DU TABLEAU PERMETTANT APRES RESOLUTION DE CALCULER
C        LES DL VITESSES MOYENNES SUR LES EF BREZZI-FORTIN
C        STOCKAGE DES LIGNES DE AE et BE DES DL ELIMINES PAR GAUSS
         NBLGAU = NDIM * NBELEM * ( NBDLMX + 1 )
         CALL TNMCDC( 'REEL2', NBLGAU, MNLGAU )
C
C        CONSTRUCTION DU TMS XYZNOEUD des XYZ des SOMMETS+BARYCENTRES DES EF
         CALL XYZNOEBF( NTLXOB, MNNPEF, NTXYZN, MNXYZN )
C
      ELSE IF( NUTYEL .EQ. 20 ) THEN
C
C        TAYLOR-HOOD 3D  TETRAEDRE 3P2C CALCUL des INTEGRALES DP2 DP2 dX
C        ===============================================================
C        LES ADRESSAGES POUR RECUPERER LES INFORMATIONS
C        DU TABLEAU DES POLYNOMES(POINTS INTEGRATION), ...
C        RECHERCHE DU TABLEAU DE POBA ET PARTAGE EN P ET DP
         L = MNTPOB
C
C        SAUT DU TABLEAU POUR LES EF 2D
C        LES VALEURS DES POLYNOMES DE L'EF DE DIMENSION TOTALE
C        EN 3D sur le TETRAEDRE DE REFERENCE
         L  = L + 1
         IA = MCN( L )
C        DIMENSION DE L ESPACE
C        NDIM   = MCN( IA )
C        NOMBRE DE POLYNOMES D INTERPOLATION
         NBPOLY = MCN( IA + 1 )
C        NOMBRE DE POINTS D INTEGRATION NUMERIQUE
         NPI    = MCN( IA + 2 )
C        ADRESSE DES POIDS DES POINTS DE LA FORMULE D INTEGRATION
         MNPOID = IA + 8
C        ADRESSE DES VALEURS DES POLYNOMES AUX POINTS D INTEGRATION
         MNPOL  = MNPOID + MCN( IA + 3 ) * MOREE2 * NPI
C        ADRESSE DES VALEURS DES DERIVEES DES POLYNOMES AUX POINTS D'INTEGRATION
         MNDPOL = MNPOL + MCN( IA + 4 ) * MOREE2 * NBPOLY * NPI
C
C        CALCUL DE L'INTEGRALE DP2 DP2 dx dy dz sur LE TETRAEDRE
C        POUR L'EF de TAYLOR-HOOD en FLUIDE
         CALL TNMCDC( 'REEL2', 3*10*3*10, MNDPDP )
         CALL INTDPDP( NPI, MCN(MNPOID), NDIM, NBPOLY, MCN(MNDPOL),
     %                 MCN(MNDPDP))
C
      ENDIF
C
C     LE TABLEAU ELEMENTAIRE BE
      MNSE = MNTAEL + NBDLMX * (NBDLMX+1) / 2  * MOREE2
C
C     LE STOCKAGE DES LIGNES ELIMINEES PAR GAUSS POUR LES EF BREZZI-FORTIN
      MNLGA = MNLGAU
C
C     CONSTRUCTION DE LA MATRICE [VG] ET DU SECOND MEMBRE {BG}
C     ========================================================
C     LA BOUCLE SUR LES ELEMENTS FINIS
      DO 22 NUELEM = 1, NBELEM
C
C        NO DES NOEUDS DE L'ELEMENT FINI NUELEM
         CALL EFNOEU( MNELE, NUELEM, NBNOEF, NONOEF )
C
C        NO DES POINTS LIGNES SURFACES VOLUMES DES SOMMETS ARETES FACES VOLUME
         CALL EFPLSV( MNELE , NUELEM,
     %                NOVCEL, NOSFEL, NOLAEL, NOPSEL,
     %                NOOBVC, NOOBSF, NOOBLA, NOOBPS )
C
C        INTERPOLATION DU FLUIDE avec F: echapeau->e  P1 ndim
C        COORDONNEES DES POINTS=NOEUDS DE L'ELEMENT FINI
         CALL EFXYZP( NDIM, MNXYZP, NBELEM, NUELEM, MNPGEL, NBNOEF,
     %                COORDP )
C
         IF( NUTYEL .EQ. 13 ) THEN
C
C           TRIANGLE BREZZI-FORTIN P1BULLEP3 EN VITESSE et P1 EN PRESSION
C           -------------------------------------------------------------
C           MATRICE ELEMENTAIRE DE VISCOSITE
            CALL F2RP1BP1( COORDP,
     %                     NOOBSF,NUMIOB(3),NUMAOB(3),MCN(MNDOEL(3)),
     %                     MCN(MNTAEL) )
C
C           SECOND MEMBRE ELEMENTAIRE
            CALL F2SP1BP1( COORDP,
     %                     NOOBLA,NUMIOB(2),NUMAOB(2),MCN(MNDOEL(2)),
     %                     NOOBSF,NUMIOB(3),NUMAOB(3),MCN(MNDOEL(3)),
     %                     MCN(MNSE) )
C
C           GAUSS SUR LE DEGRE DE LIBERTE 4 DES 2 VITESSES
C           (ELIMINATION DES 2 DL INTERNES au TRIANGLE)
            CALL GA2P1BP1( MCN(MNTAEL), MCN(MNSE), MCN(MNLGA) )
C           AVANT      4 + 4 + 3 = 11 DEGRES DE LIBERTE
C           MAINTENANT 3 + 3 + 3 =  9 DEGRES DE LIBERTE
            MNLGA = MNLGA + MOREE2 * (11+1) * 2
C
         ELSE IF( NUTYEL .EQ. 19 ) THEN
C
C           TETRAEDRE BREZZI-FORTIN P1BULLEP4 EN VITESSE et P1 EN PRESSION
C           --------------------------------------------------------------
C           MATRICE ELEMENTAIRE DE VISCOSITE
            CALL F3RP1BP1( COORDP,
     %                     NOOBVC,NUMIOB(4),NUMAOB(4),MCN(MNDOEL(4)),
     %                     MCN(MNTAEL) )
C
C           SECOND MEMBRE ELEMENTAIRE
            CALL F3SP1BP1( COORDP,
     %                     NOOBSF,NUMIOB(3),NUMAOB(3),MCN(MNDOEL(3)),
     %                     NOOBVC,NUMIOB(4),NUMAOB(4),MCN(MNDOEL(4)),
     %                     MCN(MNSE) )
C
C           GAUSS SUR LE DEGRE DE LIBERTE 5 DES 3 COMPOSANTES DE LA VITESSE
C           (ELIMINATION DES 3 DL INTERNES au TETRAEDRE)
            CALL GA3P1BP1( MCN(MNTAEL), MCN(MNSE), MCN(MNLGA) )
C           AVANT      5 + 5 + 5 + 4 = 19 DEGRES DE LIBERTE
C           MAINTENANT 4 + 4 + 4 + 4 = 16 DEGRES DE LIBERTE
            MNLGA = MNLGA + MOREE2 * (19+1) * 3
C
         ELSE IF( NUTYEL .EQ. 15 ) THEN
C
C           TRIANGLE TAYLOR-HOOD P2 EN VITESSE et P1 EN PRESSION
C           ----------------------------------------------------
C           MATRICE ELEMENTAIRE DE VISCOSITE
            CALL F2RP2P1( COORDP,
     %                    NOOBSF,NUMIOB(3),NUMAOB(3),MCN(MNDOEL(3)),
     %                    MCN(MNTAEL) )
C
C           SECOND MEMBRE ELEMENTAIRE
            CALL F2SP2P1( COORDP,
     %                    NOOBLA,NUMIOB(2),NUMAOB(2),MCN(MNDOEL(2)),
     %                    NOOBSF,NUMIOB(3),NUMAOB(3),MCN(MNDOEL(3)),
     %                    MCN(MNSE) )
C
         ELSE IF( NUTYEL .EQ. 20 ) THEN
C
C           TETRAEDRE TAYLOR-HOOD P2 EN VITESSE et P1 EN PRESSION
C           -----------------------------------------------------
C           MATRICE ELEMENTAIRE DE VISCOSITE
            CALL F3RP2P1( COORDP, MCN(MNDPDP),
     %                    NOOBVC,NUMIOB(4),NUMAOB(4),MCN(MNDOEL(4)),
     %                    MCN(MNTAEL) )
C
C           SECOND MEMBRE ELEMENTAIRE
            CALL F3SP2P1( COORDP,
     %                    NOOBSF,NUMIOB(3), NUMAOB(3), MCN(MNDOEL(3)),
     %                    NOOBVC,NUMIOB(4), NUMAOB(4), MCN(MNDOEL(4)),
     %                    MCN(MNSE) )
C
         ENDIF
C
C        ASSEMBLAGE DES TABLEAUX ELEMENTAIRES DANS LES TABLEAUX GLOBAUX
C        ==============================================================
C        NOGLDL: NUMERO GLOBAL DES DEGRES DE LIBERTE DE L'ELEMENT FINI
         CALL NUGDLE( NDIM,   NUTYEL, NBDLMX, NBNOEU,
     %                NONOEF, MCN(MNNDDL),
     %                NBDLEF, NOGLDL, IERR )
C
         IF( NORESO .EQ. 1 ) THEN
C
C           ASSEMBLAGE EN STOCKAGE PROFIL DE LA MATRICE DE VISCOSITE
C           --------------------------------------------------------
            CALL ASMEPC( NBDLEF, NOGLDL,
     %                   NCODSV, MCN(MNTAEL), MCN(MNTAEL),
     %                   NCODSV, MCN(MNPTDG), MCN(MNVG)  )
C
         ELSE IF( NORESO .GE. 2 ) THEN
C
C          ASSEMBLAGE EN STOCKAGE MORSE DE LA MATRICE DE VISCOSITE
C          -------------------------------------------------------
            CALL ASMEGC( NBDLEF, NOGLDL,
     %                   NCODSV, MCN(MNTAEL), MCN(MNTAEL),
     %                   NCODSV, MCN(MNPTDG), MCN(MNLPCO),  MCN(MNVG) )
         ENDIF
C
C        ASSEMBLAGE DU SECOND MEMBRE ELEMENTAIRE
C        ---------------------------------------
         CALL AS1BEBG( NBDLEF, NOGLDL, MCN(MNSE),  MCN(MNSMGL) )
C
 22   CONTINUE
C
C     FIN DE LA BOUCLE SUR LES EF
      IF( MNDPDP .GT. 0 ) CALL TNMCDS( 'REEL2', 3*10*3*10, MNDPDP )
C
C     CONSTRUCTION DES TABLEAUX DES NO ET VALEURS VITESSES-PRESSIONS FIXES
C     ====================================================================
      CALL VIPRFX( RELMIN, NTDL,   NDIM,   MNXYZN, MCN(MNNDDL),
     %             NBTYEL, MNNPEF, NUMIOB, MNDOEL,
     %             NBVCFX, NBDLFX, MNNDLFX, MNVDLFX )

      IF( NBDLFX .LE. 0 ) THEN
         NBLGRC(NRERR) = 2
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'ERREUR: PAS DE CONDITION AUX LIMITES'
            KERR(2) = 'SUR LA VITESSE ET PRESSION'
         ELSE
            KERR(1) = 'ERROR: NO BOUNDARY CONDITION'
            KERR(2) = 'ON VELOCITIES and PRESSURES'
         ENDIF
         CALL LEREUR
         IERR = 1
         GOTO 9000
      ENDIF

      IF( NBVCFX .GT. 0 ) THEN
         NBLGRC(NRERR) = 2
         IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = 'ERREUR: CONDITION AUX LIMITES SUR VITESSE CONVECTEE'
         KERR(2) = '        INTERDITE SUR UN FLUIDE STATIONNAIRE'
         ELSE
         KERR(1) = 'ERROR: BOUNDARY CONDITION ON CONVECTED VELOCITY'
         KERR(2) = '       FORBIDDEN ON A STEADY FLUID'
         ENDIF
         CALL LEREUR
         IERR = 2
         GOTO 9000
      ENDIF
C
C     TEMPS CPU DE FABRICATION DES TABLEAUX GLOBAUX
      DFABG  = DINFO( 'DELTA CPU' )
C
      IF( NORESO .EQ. 1 ) THEN
C
C        STOCKAGE PROFIL et FACTORISATION DE CROUT V = L D tL
C        ====================================================
C        PRISE EN COMPTE DES CONDITIONS AUX LIMITES. DL FIXES
         CALL BLDLPC( NTDL, 1, NBDLFX, MCN(MNNDLFX), MCN(MNVDLFX),
     %                NCODSV, MCN(MNPTDG), MCN(MNVG), MCN(MNSMGL) )
ccc
ccc      print *,'sortie bldlpc'
ccc      print *,(dmcn((mnsmgl-1)/2+i),i=1,ntdl)
ccc
cccc        nombre de coeff<0 de AG
ccc         call diagne( ntdl,  MCN(MNPTDG), DMCN((MNVG+1)/2), NBDINE, DIAMIN, D
C
C        AFFICHAGE DU CHOIX ici FACTORISATION DE CROUT
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'FACTORISATION DE CROUT. PATIENCE...'
         ELSE
            KERR(1) = 'CROUT FACTORIZATION. WAIT PATIENTLY...'
         ENDIF
         CALL LERESU
C
C        FACTORISATION L * D * TL
         EPSFLU=0.0
         CALL CRMC1D( MCN(MNPTDG), MCN(MNVG), NTDL, EPSFLU, 1,
     %                MCN(MNVG), NRETOU )
C
cccc        nombre de coeff<0 de AG
ccc         call diagne( ntdlr,  MCN(MNPTDG), DMCN((MNVG+1)/2), NBDINE, DIAMIN,
C
C        TEMPS CALCUL DE FORMATION DU SYSTEME LINEAIRE
         DFACTO = DINFO( 'DELTA CPU' )
         IF( NRETOU .NE. 0 ) THEN
C            MATRICE NON INVERSIBLE
             NBLGRC(NRERR) = 2
             IF( LANGAG .EQ. 0 ) THEN
                KERR(1) = 'ERREUR MATRICE NON INVERSIBLE'
                KERR(2) = 'REVOYEZ LES CONDITIONS AUX LIMITES'
             ELSE
                KERR(1) = 'ERROR NOT INVERSIBLE MATRIX'
                KERR(2) = 'SEE AGAIN the BOUNDARY CONDITIONS'
             ENDIF
             CALL LEREUR
             IF( INTERA .LE. 1 ) CALL ARRET( 100 )
             IERR = 7
             IF( MNVG .GT. 0 ) CALL TNMCDS( 'REEL2', NBCOVG, MNVG )
             GOTO 9000
         ENDIF
C
C        RESOLUTION DU SYSTEME FACTORISE
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'RESOLUTION de L D tL x = b'
         ELSE
            KERR(1) = 'SOLUTION of L D Lt x = b'
         ENDIF
         CALL LERESU
         CALL DRCRPR( NTDL,1,MCN(MNPTDG),MCN(MNVG),MCN(MNSMGL),3,
     %                MCN(MNU),IERR)

      ELSE

C        METHODE DU GRADIENT CONJUGUE AVEC STOCKAGE MORSE
C        ================================================

C        PRISE EN COMPTE SUR VG et SMGL DES CONDITIONS AUX LIMITES. DL FIXES
C        -------------------------------------------------------------------
         CALL BLDLGC( NTDL,   1, NBDLFX, MCN(MNNDLFX), MCN(MNVDLFX),
     %                NCODSV, MCN(MNPTDG), MCN(MNLPCO), MCN(MNVG),
     %                MCN(MNSMGL) )

cccc        nombre de coeff<0 de VG
ccc         call diagne( ntdl, MCN(MNPTDG), MCN(MNVG), NBDINE, DIAMIN, DIAMAX )
C
C        CALCUL DE LA MATRICE DE PRECONDITIONNEMENT SELON NIVOGC
C        -------------------------------------------------------
C        CALCUL DU SQUELETTE (POINTEURS DIAGONALE ET COLONNES)
C        DE LA MATRICE DE PRECONDITIONNEMENT DU GRADIENT CONJUGUE
         IERR = 0
         CALL CALPNT( NTDL,   NIVOGC, NCODSV, MNPTDG, MNLPCO,
     %                MNLPLC, MNLPCC, MNLPLU, COMP,   IERR )

C        VERIFICATION DE L'ETAT DE LA FACTORISATION
         LOLPCC = MCN( MNLPLC + NTDL )
         IF( IERR .NE. 0 ) THEN
            WRITE(KERR(MXLGER)(1:8),'(I8)') NIVOGC
            NBLGRC(NRERR) = 2
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'PROBLEME lors de la FACTORISATION de NIVEAU '
     %              // KERR(MXLGER)(1:8)
               KERR(2) = 'ABANDON de la METHODE du GRADIENT CONJUGUE'
             ELSE
               KERR(1)='PROBLEM of INCOMPLETE FACTORIZATION of LEVEL '
     %               // KERR(MXLGER)(1:8)
               KERR(2) = 'EXIT of the CONJUGATE GRADIENT METHOD'
            ENDIF
            CALL LEREUR
            GOTO 9000
         ENDIF

C        DECLARATION INITIALISATION DE LA MATRICE MORSE DE CONDITIONNEMENT
C        TEST DE LA PLACE MEMOIRE DISPONIBLE
         CALL TNMCMX( 'REEL2', MAXVAR )
         LO = LOLPCC
         IF( MAXVAR .LT. LO ) THEN
C            PLACE MEMOIRE INSUFFISANTE
             NBLGRC(NRERR) = 2
             WRITE(KERR(MXLGER)(1:2),'(I2)') NIVOGC
             IF( LANGAG .EQ. 0 ) THEN
                KERR(1) = 'STOKESTA: PLACE MEMOIRE INSUFFISANTE'
                KERR(2) = 'DIMINUEZ la VALEUR du NIVEAU '
     %                 // KERR(MXLGER)(1:2)
             ELSE
                KERR(1) = 'STOKESTA: NOT ENOUGH MEMORY'
                KERR(2) = 'DIMINISH the LEVEL VALUE'
     %                 // KERR(MXLGER)(1:2)
             ENDIF
             GOTO 9000
         ENDIF
C
C        DECLARATION DE LA MATRICE DE PRECONDITIONNEMENT
         CALL TNMCDC( 'REEL2', LOLPCC, MNVGC )
C
C        MATRICE MORSE
C        CONSTRUCTION DU TABLEAU DL FIXES OU NON?
         CALL TNMCDC( 'ENTIER', NTDL, MNNDFX )
C        A PRIORI TOUS LES DL SONT LIBRES => 0
         CALL AZEROI( NTDL, MCN(MNNDFX) )
C        MISE A 1 DES DL FIXES
         DO I=1, NBDLFX
C           LE NO GLOBAL DU DL FIXE
            NDL = MCN( MNNDLFX - 1 + I )
C           MISE A 1 DU DL FIXE NDL
            MCN( MNNDFX-1 + NDL ) = 1
         ENDDO
C
C        AFFICHAGE
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'CROUT FACTORISATION INCOMPLETE PAR NIVEAUX'
         ELSE
            KERR(1) = 'CROUT INCOMPLETE FACTORIZATION by LEVELS'
         ENDIF
         CALL LERESU
C
C        CONSTRUCTION DE LA MATRICE DE PRECONDITIONNEMENT
         CALL INCRCO( NTDL,        MCN(MNNDFX),
     &                MCN(MNPTDG), MCN(MNLPCO), MCN(MNVG),
     &                MCN(MNLPLC), MCN(MNLPCC), MCN(MNVGC),
     &                IERR )
C
C        FACTORISATION INCOMPLETE DE CROUT L D tL
         CALL INCRFA( NTDL, MCN(MNLPLC), MCN(MNLPCC), MCN(MNVGC), IERR )
C
C        TEMPS CALCUL DE FORMATION DU SYSTEME LINEAIRE
         DFACTO = DINFO( 'DELTA CPU' )
C
C        VERIFICATION DE LA STABILITE
         IF( IERR .GT. 0 ) THEN
C           LA FACTORISATION EST INSTABLE : ABANDON DES CALCULS
            WRITE(KERR(MXLGER)(1:8),'(I8)') NIVOGC
            NBLGRC(NRERR) = 2
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'FACTORISATION INSTABLE AU NIVEAU '
     &              // KERR(MXLGER)(1:8)
               KERR(2) = 'ABANDON DE LA METHODE DU GRADIENT CONJUGUE'
            ELSE
               KERR(1)='UNSTABLE INCOMPLETE FACTORIZATION at LEVEL '
     %              // KERR(MXLGER)(1:8)
               KERR(2) = 'EXIT of the CONJUGATE GRADIENT METHOD'
            ENDIF
            CALL LEREUR
            IERR = 8
            GO TO 9000
         ENDIF
C
C        RESOLUTION DU SYSTEME FACTORISE
C        ===============================
         IF( NORESO .EQ. 2 ) THEN
C
C           LES TABLEAUX AUXILIAIRES DE GCPRCR
            LOAUX = NTDL * 4
            CALL TNMCDC( 'REEL2', LOAUX, MNAUX )
C           MISE A ZERO DU TABLEAU AUXILIAIRE
            CALL AZEROD( LOAUX, MCN(MNAUX) )
C           DECOUPAGE EN 4 TABLEAUX DOUBLE PRECISION
            MNAUX1 = MNAUX
            MNAUX2 = MNAUX1 + NTDL * MOREE2
            MNAUX3 = MNAUX2 + NTDL * MOREE2
            MNAUX4 = MNAUX3 + NTDL * MOREE2
C           DECLARATION DE TABLEAUX AUXILIAIRES
            NBDIR = 3
C           NOMBRE DE DIRECTIONS >=1 POUR STABILISER LE GC
            LODIR = (NTDL+1) * NBDIR * 2
            CALL TNMCDC( 'REEL2', LODIR, MNBDIR )
C           MISE A ZERO
            CALL AZEROD( LODIR, MCN(MNBDIR) )
C           DECOUPAGE EN 3 TABLEAUX DOUBLE PRECISION
            MNDIR  = MNBDIR
            MNADIR = MNDIR  + NTDL * NBDIR * MOREE2
            MNDAD  = MNADIR + NTDL * NBDIR * MOREE2
            MNBET  = MNDAD  + NBDIR * MOREE2
C
C           AFFICHAGE DU DEBUT DES ITERATIONS DE GC
            NBLGRC(NRERR) = 1
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'DEBUT DES ITERATIONS DE GRADIENT CONJUGUE'
            ELSE
               KERR(1) = 'START of CONJUGATE GRADIENT ITERATIONS'
            ENDIF
            CALL LERESU
C
         ELSE IF( NORESO .EQ. 3 ) THEN
C
C           LES TABLEAUX AUXILIAIRES DE CGSQSY
            LOAUX = NTDL * 8
            CALL TNMCDC( 'REEL2', LOAUX, MNAUX )
C           MISE A ZERO DU TABLEAU AUXILIAIRE
            CALL AZEROD( LOAUX, MCN(MNAUX) )
C           DECOUPAGE EN 2 TABLEAUX DOUBLE PRECISION
            MNAUX1 = MNAUX
            MNAUX2 = MNAUX1 + NTDL * MOREE2
C
C           AFFICHAGE DU DEBUT DES ITERATIONS DE GC
            NBLGRC(NRERR) = 1
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'DEBUT DES ITERATIONS DE CGS'
            ELSE
               KERR(1) = 'START of CGS ITERATIONS'
            ENDIF
            CALL LERESU
C
         ELSE IF( NORESO .EQ. 4 ) THEN
C
C           LES TABLEAUX AUXILIAIRES DE BICGST
            LOAUX = NTDL * 8
            CALL TNMCDC( 'REEL2', LOAUX, MNAUX )
C           MISE A ZERO DU TABLEAU AUXILIAIRE
            CALL AZEROD( LOAUX, MCN(MNAUX) )
C           DECOUPAGE EN 8 TABLEAUX DOUBLE PRECISION
            MNAUX1 = MNAUX
            MNAUX2 = MNAUX1 + NTDL * MOREE2
            MNAUX3 = MNAUX2 + NTDL * MOREE2
            MNAUX4 = MNAUX3 + NTDL * MOREE2
            MNAUX5 = MNAUX4 + NTDL * MOREE2
            MNAUX6 = MNAUX5 + NTDL * MOREE2
            MNAUX7 = MNAUX6 + NTDL * MOREE2
            MNAUX8 = MNAUX7 + NTDL * MOREE2
C
C           AFFICHAGE DU DEBUT DES ITERATIONS DE GC
            NBLGRC(NRERR) = 1
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'DEBUT DES ITERATIONS DE BICG STABILISE'
            ELSE
               KERR(1) = 'START of BICG STABILISED ITERATIONS'
            ENDIF
            CALL LERESU
         ENDIF
C
C        PRISE EN COMPTE DES CONDITIONS AUX LIMITES SUR
C        INITIALISATION DE LA SOLUTION POUR LES GC
         IF( NBDLFX .NE. 0 ) THEN
            CALL BLDLX0( NTDL, 1, NBDLFX, MCN(MNNDLFX), MCN(MNVDLFX),
     %                   MCN(MNAUX1) )
         ENDIF
C
         IF( NORESO .EQ. 2 ) THEN
C
C           RESOLUTION PAR GRADIENT CONJUGUE PRECONDITIONNE PAR CROUT INCOMPLET
            CALL GCPRCR( NTDL, 1,  NBDIR,    MCN(MNNDFX),
     &                MCN(MNPTDG), MCN(MNLPCO), MCN(MNVG),  MCN(MNSMGL),
     &                MCN(MNLPLC), MCN(MNLPCC), MCN(MNVGC),
     &                MCN(MNAUX1), MCN(MNAUX2), MCN(MNAUX3),MCN(MNAUX4),
     &                MCN(MNDIR),  MCN(MNADIR),
     &                MCN(MNDAD),  MCN(MNBET),
     &                MCN(MNU),    IERR )
C
C           DESTRUCTION DES TABLEAUX AUXILIAIRES DE GC CROUT INCOMPLET
            CALL TNMCDS( 'REEL2', LODIR, MNBDIR )
C
         ELSE IF( NORESO .EQ. 3 ) THEN
C
C           RESOLUTION PAR 'CONJUGATE GRADIENTS SQUARED' SUR MATRICE SYMETRIQUE
            CALL CGSQSY( NTDL,     MCN(MNNDFX),
     &                 MCN(MNPTDG), MCN(MNLPCO), MCN(MNVG), MCN(MNSMGL),
     &                 MCN(MNLPLC), MCN(MNLPCC), MCN(MNVGC),
     &                 MCN(MNAUX1), MCN(MNAUX2),
     &                 MCN(MNU),    IERR )
C
         ELSE IF( NORESO .EQ. 4 ) THEN
C
C           RESOLUTION 'BI-CONJUGATE GRADIENT STABILISED' SUR MATRICE SYMETRIQE
            CALL  BICGST( NTDL,     MCN(MNNDFX),
     &               MCN(MNPTDG), MCN(MNLPCO), MCN(MNVG),   MCN(MNSMGL),
     &               MCN(MNLPLC), MCN(MNLPCC), MCN(MNVGC),  MCN(MNAUX1),
     &               MCN(MNAUX2), MCN(MNAUX3), MCN(MNAUX4), MCN(MNAUX5),
     &               MCN(MNAUX6), MCN(MNAUX7), MCN(MNAUX8),
     &               MCN(MNU),    IERR )
C
         ENDIF
C
C        DESTRUCTION DES TABLEAUX AUXILIAIRES
         CALL TNMCDS( 'REEL2', LOAUX, MNAUX )
C
C        VERIFICATION DE LA CONVERGENCE
         IF( IERR .NE. 0 ) THEN
C
C           LE GC NE CONVERGE PAS : ABANDON DES CALCULS
            WRITE(KERR(MXLGER)(1:8),'(I8)') NIVOGC
            NBLGRC(NRERR) = 2
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'NON CONVERGENCE AU NIVEAU '
     &              // KERR(MXLGER)(1:8)
               KERR(2) = 'ABANDON DE LA METHODE ITERATIVE'
            ELSE
               KERR(1) = 'NO CONVERGENCE at LEVEL '
     &              // KERR(MXLGER)(1:8)
               KERR(2) = 'EXIT of the ITERATIVE METHOD'
            ENDIF
            CALL LEREUR
            IERR = 9
            GOTO 9000
         ENDIF
C
      ENDIF
C
      IF( NUTYEL .EQ. 13 .OR. NUTYEL .EQ. 19 ) THEN
C
C        TRIANGLE ou TETRAEDRE BREZZI-FORTIN
C        ===================================
C        CALCUL DES DEGRES DE LIBERTE ELIMINES DANS LES EF
C        PAR LA METHODE DE GAUSS
         CALL DLVIBF( NDIM, MNNPEF, NBDLMX+1, NBELEM, MCN(MNLGAU),
     %                NTDL, NTDLF, MCN(MNU) )
C
      ENDIF
C
      IF ( IEBLPR .LE. 0 ) THEN
C
C        PAS DE CONDITION AUX LIMITES SUR LA PRESSION
C        ============================================
C        => ELLE EST DEFINIE A UNE CONSTANTE PRES
C        RECHERCHE DE LA PLUS BASSE PRESSION POUR SERVIR DE TRANSLATION
C        ELLE EST SOUSTRAITE DES AUTRES PRESSIONS ET
C        LA PLUS BASSE PRESSION DEVIENT ZERO
C
         MN = (MNU-1)/2
         PREMIN = 1D100
         NDL1 = MCN(MNNDDL)
         DO N=1,NBNOEU
            NDL2 = MCN(MNNDDL+N)
            IF( NDL2 - NDL1 .EQ. NDIM+1 ) THEN
               IF( DMCN(MN+NDL2) .LT. PREMIN ) THEN
                  PREMIN = DMCN(MN+NDL2)
               ENDIF
            ENDIF
            NDL1 = NDL2
         ENDDO
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,*)'PRESSION MINIMALE ',PREMIN,' RAMENEE a ZERO'
         ELSE
         WRITE(IMPRIM,*)'MINIMUM PRESSURE ',PREMIN,' TRANSLATED to ZERO'
         ENDIF
C
         NDL1 = MCN(MNNDDL)
         DO N=1,NBNOEU
            NDL2 = MCN(MNNDDL+N)
            IF( NDL2 - NDL1 .EQ. NDIM+1 ) THEN
               DMCN(MN+NDL2) = DMCN(MN+NDL2) - PREMIN
            ENDIF
            NDL1 = NDL2
         ENDDO
C
      ENDIF
C
C     MISE A JOUR DU TMS 'VECTEUR"VITESSEPRESSION'
C     ============================================
      MCN( MNVIPE + WBCOVE ) = NTDLF
      MCN( MNVIPE + WBVECT ) = 1
      MCN( MNVIPE + WBCPIN ) = 0
C     LA DATE
      CALL ECDATE( MCN(MNVIPE) )
C     LE NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNVIPE + MOREE2 ) = NONMTD( '~>>>VECTEUR' )

C     COUT CALCUL DES VITESSES PRESSIONS
C     ==================================
      DVITPR = DINFO( 'DELTA CPU' )

C     AFFICHAGE DU DERNIER VECTEUR"VITESSEPRESSION CALCULE
C     ====================================================
      CALL AFVIPR( NUTYEL, NDIM,   NBNOVI, MNXYZN, MCN(MNNDDL),
     %             TEMPS,  NTDLF,  MIN(4,NBNOVI),  MCN(MNU),
     %             VITMIN, VITMAX, VITMOY, PREMIN, PREMAX, PREMOY )

C     MISE SUR FICHIER DU VECTEUR NOVVIPR VITESSEPRESSION UG0
C     AU TEMPS INITIAL TPSINI
      NOVVIPR = 1
      NAVSTO  = -1
      NBPASDT = 0
C     LE NOM DU FICHIER du VECTEUR VITESSE-PRESSION
      CALL NMFIVIPRTE( KNOMOB, TEMPS,  NOVVIPR, KNOMFIC, NBK )
      CALL ECFIVIPRTE( KNOMOB, TEMPS,  NOVVIPR, NAVSTO, NBPASDT,
     %                 NDIM,   NBNOVI, NBSOM, NBNOTE,
     %                 NTDLF,  MCN(MNU), NTDLTE, TEMPER,
     %                 KNOMFIC, IERR )

C .....................................................................

C     FERMETURE DE LA MATRICE POUR REDONNER DE LA PLACE EN MC
C     =======================================================
 9000 IF( MNNDFX .GT. 0 ) CALL TNMCDS( 'ENTIER', NTDL,   MNNDFX )
      IF( MNPTDG .GT. 0 ) CALL TNMCDS( 'ENTIER', NTDL+1, MNPTDG )
      IF( MNLPCO .GT. 0 ) CALL TNMCDS( 'ENTIER', NBCOVG, MNLPCO )
      IF( MNVG   .GT. 0 ) CALL TNMCDS( 'REEL2',  NBCOVG, MNVG   )
      IF( MNLPLC .GT. 0 ) CALL TNMCDS( 'ENTIER', NTDL+1, MNLPLC )
      IF( MNLPCC .GT. 0 ) CALL TNMCDS( 'ENTIER', LOLPCC, MNLPCC )
      IF( MNVGC  .GT. 0 ) CALL TNMCDS( 'REEL2',  LOLPCC, MNVGC  )

C     DESTRUCTION DES TABLEAUX TEMPORAIRES
C     ====================================
      IF( MNNPEF .GT. 0 ) CALL TNMCDS( 'ENTIER',2*MXTYEL, MNNPEF )
      DO I=1,4
         IF( MNDOEL(I) .GT. 0 .AND. MXDOEL(I) .GT. 0 ) THEN
            CALL TNMCDS( 'ENTIER', MXDOEL(I), MNDOEL(I) )
         ENDIF
      ENDDO
      IF( MNLGAU .GT. 0 ) CALL TNMCDS( 'REEL2',  NBLGAU,   MNLGAU )
      IF( MNDPDP .GT. 0 ) CALL TNMCDS( 'REEL2', 3*10*3*10, MNDPDP )
      IF( MNNDDL .GT. 0 ) CALL TNMCDS( 'ENTIER', 1+NBNOVI, MNNDDL )
      IF( MNSMGL .GT. 0 ) CALL TNMCDS( 'REEL2',  NTDL,     MNSMGL )
      IF( MNTAEL .GT. 0 ) CALL TNMCDS( 'REEL2',  MOTAEL,   MNTAEL )
      IF( MNNODL .GT. 0 ) CALL TNMCDS( 'ENTIER', NBDLMX,   MNNODL )
      IF( MNNDLFX.GT. 0 ) CALL TNMCDS( 'ENTIER', NBDLFX,   MNNDLFX)
      IF( MNVDLFX.GT. 0 ) CALL TNMCDS( 'REEL2',  NBDLFX,   MNVDLFX)
C
C     GESTION DES ERREURS
C     ===================
      IF( IERR .EQ. 7 ) THEN
C        RETOUR SI MATRICE NON INVERSIBLE
         IERR = 0
         RETURN
      ELSE IF( IERR .EQ. 8 ) THEN
C        RETOUR SI FACTORISAION INCOMPLETE INSTABLE ET TRAVAIL INTERACTIF
         IERR = 0
         GOTO 10
      ELSE IF( IERR .EQ. 9 ) THEN
C        RETOUR SI NON CONVERGENCE DU G.C. ET TRAVAIL INTERACTIF
         IERR = 0
         GOTO 10
      ELSE IF( IERR .NE. 0 ) THEN
         RETURN
      ENDIF
C
C     AFFICHAGE DES TEMPS CALCUL
      DCPU = DINFO( 'DELTA CPU' )
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,12001)DFABG,DFACTO,DVITPR,DFABG+DFACTO+DVITPR+DCPU
      ELSE
         WRITE(IMPRIM,22001)DFABG,DFACTO,DVITPR,DFABG+DFACTO+DVITPR+DCPU
      ENDIF
12001 FORMAT(/
     % 'TEMPS FORMATION     DE LA MATRICE  =',F12.2,' SECONDES CPU'/,
     % 'TEMPS FACTORISATION DE LA MATRICE  =',F12.2,' SECONDES CPU'/,
     % 'TEMPS CALCUL DES VITESSES PRESSIONS=',F12.2,' SECONDES CPU'/
     % 'TEMPS CALCUL TOTAL                 =',F12.2,' SECONDES CPU'/)
22001 FORMAT(/
     % 'MATRIX FORMATION              TIME=',F12.2,' CPU SECONDS'/,
     % 'MATRIX FACTORIZATION          TIME=',F12.2,' CPU SECONDS'/,
     % 'VELOCITY-PRESSURE COMPUTATION TIME=',F12.2,' CPU SECONDS'/
     % 'TOTAL COMPUTATION             TIME=',F12.2,' CPU SECONDS'/)
      RETURN


C     ERREUR A L'OUVERTURE DU FICHIER POBA
 9900 NBLGRC(NRERR) = 1
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = 'IMPOSSIBLE D''OUVRIR LE FICHIER POBA'
      ELSE
         KERR(1) = 'IMPOSSIBLE to OPEN the FILE POBA'
      ENDIF
      CALL LEREUR

      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10999) KNOMOB
      ELSE
         WRITE(IMPRIM,20999) KNOMOB
      ENDIF
10999 FORMAT(/,1X,100('=')/,
     %' stokesta: FIN RESOLUTION de l''EQUATION de STOKES STATIONNAIRE d
     %ans',A,/1X,100('='))
20999 FORMAT(/,1X,100('=')/,
     %' stokesta: END of STUDY of the STOKES STEADY FLUID FLOW in ',A,/,
     %  1X,100('='))

      RETURN
      END
