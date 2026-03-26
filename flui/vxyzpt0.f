      SUBROUTINE VXYZPT0( NAVSTO, NTCALCUL, KNOMOB, RELMIN, KNOMFIC,
     %                    NTLXOB, NDIM,   MNXYZN, MNNPEF, NUMIOB,
     %                    NBTYEL, NUTYEL, MNDOEL, MNDTEL, NDPGST,
     %                    NBELEM, NBSOM,  NBNOHB, NBNOVI, NBNOTE,
     %                    NTDLHB, NTDLVP, NTDLTE,
     %                    IEVTIN, IEPRIN, IETEIN,
     %                    NBVCFX, NBVPFX,  MNNVPFX, MNVVPFX,
     %                    NBTEFX, MONTEFX, MNNTEFX, MNVTEFX,
     %                    MNTIMES,TEMPS0,  NOVVIPR, NDDLNO, NBPAST0,
     %                    VXYZP0, TEMPER0, IERR )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : PROBLEME DE NAVIER et/ou STOKES INSTATIONNAIRE ou BOUSSINESQ
C ----- CALCUL DU VECTEUR VITESSE+PRESSION VXYZP0 AU TEMPS0 INITIAL
C       VXYZP0 = VX1,VY1,VZ1,PR1(AU NOEUD1) ,  ... ,
C                VXn,VYn,VZn,PRn(AU NOEUD nbnoeu)
C       et a la suite pour EF BREZZI-FORTIN (PAS pour les EF TAYLOR-HOOD)
C                VXn+1,VYn+1,VZn+1(AU BARYCENTRE DE L'EF 1), ...
C                Vm,   VYm,  VZm  (AU BARYCENTRE du DERNIER EF)

C       3 ETAPES DE L'INITIALISATION de VXYZP0:
C       1) RECHERCHE du FICHIER du VECTEUR VITESSE+PRESSION
C          LE PLUS PROCHE du TEMPS INITIAL
C          Si UN FICHIER EXISTE, SON CONTENU EST LU POUR VXYZP0
C                                et PASSAGE AUX CONDITIONS AUX LIMITES
C          sinon, TOUTES LES COMPOSANTES de VXYZP0 SONT MISES A ZERO.

C       2) SI un ou les TMS VITFLUIN et PREFLUIN des OVSLP de l'OBJET EXISTENT
C          ALORS LEURS VALEURS sont APPLIQUEES SUR VXYZP0

C       3) LES DEGRES DE LIBERTE VITESSE-PRESSION FIXES du MAILLAGE
C          SONT IMPOSES

C       SI NAVSTO=3 ALORS
C          CALCUL DU VECTEUR TEMPERATURE TEMPER0 AU TEMPS0 INITIAL

C ENTREES:
C --------
C NAVSTO :(-1POUR RESOUDRE LE PROBLEME DE STOKES   STATIONNAIRE Cf stokesta.f)
C          0 POUR RESOUDRE LE PROBLEME DE STOKES INSTATIONNAIRE Cf stokesins.f
C          1 POUR RESOUDRE LE PROBLEME DE NAVIER-STOKES SOUS FORME IMPLICITE
C          2 POUR RESOUDRE LE PROBLEME DE NAVIER-STOKES PAR PAS FRACTIONNAIRES
C          3 POUR RESOUDRE LE PROBLEME DE NAVIER-STOKES + THERMIQUE BOUSSINESQ

C          4 POUR RESOUDRE LE PROBLEME DE NAVIER-STOKES PAR PISO+CARACTERISTIQUE
C          5 POUR RESOUDRE LE PROBLEME DE NAVIER-STOKES PAR PISO+(V.D) V)

C NTCALCUL:1 DEPART du CALCUL avec LECTURE du FICHIER VITESSE+PRESSION+TEMPERATURE
C            de TEMPS le PLUS PROCHE du TEMPS INITIAL TEMPS0 et ZERO SINON
C            puis, MODIFICATION a partir des DONNEES a TEMPS0 du FLUIDE des VSLP
C            puis, PRISE EN COMPTE DES CONDITIONS AUX LIMITES a TEMPS0 des SLP
C          2 DEPART avec la DONNEE aux NOEUDS du VECTEUR VITESSE PRESSION a TEMPS0
C          3 REPRISE avec le VECTEUR VITESSE PRESSION de TEMPS le + PROCHE TEMPS0

C KNOMOB : NOM DE L'OBJET A TRAITER
C RELMIN : REEL DOUBLE PRECISION TEMOIN DE NON-INITIALISATION

C KNOMFIC: NOM DU FICHIER SUPPORT DU VECTEUR NOVVIPR des
C          VITESSES+PRESSIONS(+TEMPERATURES)
C NTLXOB : NUMERO DU TMS DU LEXIQUE DE L'OBJET KNOMOB
C MNTOPO : NUMERO DU TMS TOPOLOGIE  DE L'OBJET KNOMOB
C NDIM   : DIMENSION DES COORDONNEES DES POINTS ( 2 OU 3 )
C MNXYZP : ADRESSE MCN DU TABLEAU XYZPOINT DE L'OBJET KNOMOB
C MNXYZN : ADRESSE MCN DU TABLEAU XYZNOEUD DE L'OBJET KNOMOB

C NBTYEL : NOMBRE DE TYPES D'EF DU MAILLAGE DE CET OBJET
C MNNPEF : ADRESSE MCN DU TABLEAU DES ADRESSES MCN DES TMS NPEF"TYPE EF
C NBDLMX : NOMBRE MAXIMAL DE DEGRES DE LIBERTE D'UN EF
C
C NUMIOB : NUMERO MINIMAL DU VSLP DANS LA DEFINITION DE L'OBJET
C MXTYEL : NOMBRE MAXIMAL DE TYPES D'EF (7)
C MNDOEL : LES 4 ADRESSES MCN DES TABLEAUX DES ADRESSES DES
C          TABLEAUX DECRIVANT LES DONNEES DU FLUIDE DE L'OBJET COMPLET
C MNDTEL : LES 4 ADRESSES MCN DES TABLEAUX DES ADRESSES DES
C          TABLEAUX DECRIVANT LES DONNEES THERMIQUES DU FLUIDE DE L'OBJET

C NBELEM : NOMBRE D'ELEMENTS FINIS DE CET UNIQUE TYPE D'EF
C NBSOM  : NOMBRE DE SOMMETS DU MAILLAGE
C NBPOI  : NOMBRE DE POINTS  DU MAILLAGE
C NBNOHB : NOMBRE DE NOEUDS DU MAILLAGE HORS BARYCENTRES POUR EF BREZZI-FORTIN
C          SOMMETS                      des EF BREZZI-FORTIN
C          SOMMETS + MILIEUX DES ARETES des EF TAYLOR-HOOD
C NBNOVI : NOMBRE DE NOEUDS DU MAILLAGE AVEC BARYCENTRES POUR BREZZI-FORTIN
C          SOMMETS + BARYCENTRES        des TETRAEDRES POUR BREZZI-FORTIN
C          SOMMETS + MILIEUX DES ARETES des TETRAEDRES POUR TAYLOR-HOOD
C NBNOTE : NOMBRE DE NOEUDS SUPPORT DE LA TEMPERATURE
C          = 0 SI NAVSTO NON EGAL 3  BOUSSINESQ
C          = NBNOVI SINON (=NOMBRE DE NOEUDS P2 POUR EF TAYLOR-HOOD)

C NBVCFX : NOMBRE DE DEGRES DE LIBERTE DE VITESSE CONVECTEE ou PRESSION IMPOSE
C NBVPFX : NOMBRE DE DL VITESSE OU PRESSION FIXES
C MNNVPFX: ADRESSE MCN DU TABLEAU DU NO GLOBAL DE CHAQUE DL FIXE
C MNVVPFX: ADRESSE MCN DU TABLEAU DE LA VALEUR DE CHAQUE DL FIXE

C NBTEFX : NOMBRE DE DL SUPPORT D'UNE TEMPERATURE FIXEE
C MONTEFX: NOMBRE DE MOTS DECLARES DU TABLEAU MC NO DES DL FIXES
C MNNTEFX: ADRESSE MCN DU TABLEAU MC DES NUMEROS DES DL FIXES, 0 SINON
C MNVTEFX: ADRESSE MCN DU TABLEAU MC DES VALEURS DES DL FIXES, 0 SINON

C DT     : PAS CONSTANT DU TEMPS
C THETA  : COEFFICIENT DU SCHEMA= THETA(1), THETA(0)=1D0-THETA(0)
C DTSTOC : PAS CONSTANT DU TEMPS ENTRE 2 STOCKAGES DU VECTEUR"VITESSEPRESSION
C TEMPS0 : TEMPS INITIAL DU CALCUL DU VECTEUR VITESSE+PRESSION

C NOVVIPR: >0 NUMERO du VECTEUR VITESSE PRESSION de TEMPS le PLUS PROCHE
C             de TEMPS0 QUI EXISTE SUR FICHIER DANS LE REPERTOIRE PROJET
C             DEJA RETROUVE dans fluideNS.f
C          =0 SI PAS DE TEL FICHIER

C NTDLHB : NOMBRE TOTAL DE DEGRES DE LIBERTE ASSOCIES AUX NOEUDS
C          SANS LES DL VITESSE AUX BARYCENTRES SI EF BF
C NTDLVP : NOMBRE TOTAL DE DEGRES DE LIBERTE VITESSE+PRESSION DE L'OBJET
C          AVEC LES DL VITESSES AUX BARYCENTRES SI EF BREZZI-FORTIN
C          POUR LES EF TAYLOR-HOOD NTDLHB=NTDLVP
C NTDLTE : >0 NOMBRE TOTAL DE DEGRES DE LIBERTE TEMPERATURE AUX NOEUDS
C          =0 SI NAVSTO NON EGAL 3 (CAR TEMPERATURE NON CALCULEE)

C IEVTIN : NOMBRE DE TMS VITESSE INITIALE DES SV DE L'OBJET RETROUVES
C IEPRIN : NOMBRE DE TMS PRESSION INITIALEDES SV DE L'OBJET RETROUVES

C NBVPFILE: NOMBRE DE FICHIERS VECTEURS VITESSE+PRESSION RETROUVES
C           NOMBRE DE TEMPS DU TABLEAU MCN D'ADRESSE MNTIMES(1:NBVPFILE)
C NDDLNO : TABLEAU POINTEURS SUR LE DERNIER DL DE CHAQUE NOEUD

C SORTIES:
C --------
C NBPAST0: NOMBRE DE PAS DE TEMPS DEJA CALCULES
C VXYZP0 : TABLEAU VITESSEPRESSION PAR NOEUDS A L'INSTANT INITIAL
C TEMPER0: TABLEAU TEMPERATURE     PAR NOEUDS A L'INSTANT INITIAL
C IERR   :=0  SI PAS D'ERREUR D'EXECUTION, NON NUL SINON
C         >0  SI ABANDON DANS L'ENTREE DES DONNEES
C         =1  SI OPTION DE DEMARRAGE 3 SANS FICHIER VITESSE+PRESSION
C         =2  SI INCOHERENCE NDIM et NDIM LU
C         =3  SI INCOHERENCE NBNOVI et NBNOVI LU
C         =4  SI INCOHERENCE NBNOPR=NBSOM et NBNOPR LU
C         =5  SI INCOHERENCE NBNOTE et NBNOTE LU
C         =20 SI PAS ASSEZ DE MEMOIRE POUR STOCKER UNE MATRICE GLOBALE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET Veulettes & St Pierre du Perray      Aout 2020
C MODIFS: ALAIN PERRONNET             St Pierre du Perray   Fevrier 2021
C MODIFS: ALAIN PERRONNET             St Pierre du Perray   Fevrier 2022
C23456---------------------------------------------------------------012
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
      include"./incl/a___vecteur.inc"
      include"./incl/a___profil.inc"
      include"./incl/a___force.inc"
      include"./incl/a___vitfluin.inc"
      include"./incl/a___prefluin.inc"
      include"./incl/a___morse.inc"
      include"./incl/msvaau.inc"
      include"./incl/ponoel.inc"
      include"./incl/ctemps.inc"
      include"./incl/cthet.inc"
      include"./incl/cnonlin.inc"

      include"./incl/pp.inc"
      COMMON             MCN(MOTMCN)
      REAL              RMCN(1)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE      (MCN(1), RMCN(1), DMCN(1))
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)

      INTEGER           NUMIOB(4), MNDOEL(4), MNDTEL(4)
      INTEGER           NDDLNO(0:NBNOVI)
      INTEGER           NOTYOB(3,MXNOEL)
      INTEGER           NOOBSF(6), NOOBLA(12), NOOBPS(8)

      CHARACTER*(*)     KNOMOB, KNOMFIC
      CHARACTER*80      NOMTS
      CHARACTER*80      KNOMTD

      REAL              TEMPS0

      DOUBLE PRECISION  VXYZP0(NTDLVP), VITES, VITESSE(3), PRESSION
      DOUBLE PRECISION  TEMPER0(NTDLTE)
      DOUBLE PRECISION  PARFON(7)
      DOUBLE PRECISION  XPOIN,  YPOIN, ZPOIN, DVAL, RELMIN


C     QUELQUES INITIALISATIONS
C     ========================
      PRINT*
      IF(LANGAG .EQ. 0 ) THEN
       PRINT*,'vxyzpt0: Recherche du VECTEUR Vitesse+Pression(+Temperatu
     %re) au TEMPS',TEMPS0,' pour l''option NAVSTO=',NAVSTO,
     %' NO FICHIER LU=',NOVVIPR
      ELSE
       PRINT*,'vxyzpt0: Search of the Velocity+Pressure(+Temperature) VE
     %CTOR at TIME',TEMPS0,' with the option NAVSTO=',NAVSTO,
     %' READ FILE NUMBER=',NOVVIPR
      ENDIF

      IERR = 0

C     PROTECTION DES ADRESSES POUR EVITER DES PROBLEMES LORS
C     DE LA DESTRUCTION DES TABLEAUX
      MNNFNX = 0
      MNVFNX = 0
      MNNDLX = 0
      MNVDLX = 0
      MNFLPT = 0
      MNFLTO = 0
      MNERTH = 0
      NBFNFX = 0
      MONFNX = 0
      MNNFNX = 0
      MOVFNX = 0
      MNVFNX = 0
      NTVECT = 0
      MNVECT = 0
      NBVECT = 0

C     ==============================================================================
C     REPRISE ou INITIALISATION DU VECTEUR VITESSES PRESSIONS TEMPERATURES INITIALES
C     ==============================================================================
      IF( NOVVIPR .GT. 0 ) THEN

C        LECTURE du VECTEUR VITESSE PRESSION NOVVIPR au TEMPS0
C        QUI EXISTE DEJA SUR FICHIER DANS LE REPERTOIRE du PROJET
C        ========================================================
C        LE TEMPS INITIAL EST CELUI DU FICHIER LU
         TEMPS   = RMCN( MNTIMES-1 + NOVVIPR )
         TEMPS0  = TEMPS
         TEMPSINI= TEMPS

C        LECTURE du VECTEUR VITESSE+PRESSION+TEMPERATURE PAR NOEUDS AU TEMPS TEMPS
C        FICHIER DU REPERTOIRE PROJET DANS LE TABLEAU VXYZP0
         CALL LIFIVIPRTE( KNOMOB,  TEMPS,   NOVVIPR, NAVSTO0, NBPAST0,
     %                    NDIM0,   NBNOVI0, NBNOPR0, NBNOTE0,
     %                    NTDLVP,  VXYZP0,  NTDLTE0, TEMPER0,
     %                    KNOMFIC, IERR )

         IF( IERR .GT. 0 ) THEN
C           FICHIER NOVVIPR VITESSE+PRESSION+TEMPERATURE N'EXISTE PAS ou EST INCORRECT
            IF(LANGAG .EQ. 0 ) THEN
               PRINT*,'INCORRECT FICHIER VITESSE+PRESSION+TEMPERATURE au
     % TEMPS',TEMPS0
               PRINT*,'REVOIR LE CALCUL de VITESSE+PRESSION+TEMPERATURE 
     %INITIAL'
            ELSE
               PRINT*,'INCORRECT VELOCITY+PRESSURE+TEMPERATURE FILE at T
     %IME',TEMPS0
               PRINT*,'COMPUTE AGAIN THE INITIAL VELOCITY+PRESSURE+TEMPE
     %RATURE'
            ENDIF
            IERR = 1
            GOTO 9999
         ENDIF

         IF( NDIM0 .NE. NDIM ) THEN
            IF(LANGAG .EQ. 0 ) THEN
               PRINT*,'INCOHERENCE NDIM=',NDIM,' NDIM LU=',NDIM0
            ELSE
               PRINT*,'INCORRECT NDIM=',NDIM,' READ NDIM=',NDIM0
            ENDIF
            IERR = 2
            GOTO 9999
         ENDIF

         IF( NBNOVI0 .NE. NBNOVI ) THEN
            IF(LANGAG .EQ. 0 ) THEN
               PRINT*,'INCOHERENCE NBNOVI=',NBNOVI,' NBNOVI LU=',NBNOVI0
            ELSE
               PRINT*,'INCORRECT NBNOVI=',NBNOVI,' READ NBNOVI=',NBNOVI0
            ENDIF
            IERR = 3
            GOTO 9999
         ENDIF

         IF( NBSOM .NE. NBNOPR0 ) THEN
            IF(LANGAG .EQ. 0 ) THEN
               PRINT*,'INCOHERENCE NBSOM=',NBSOM,' NBNOPR LU=',NBNOPR0
            ELSE
               PRINT*,'INCORRECT NBSOM=',NBSOM,' READ NBNOPR=',NBNOPR0
            ENDIF
            IERR = 4
            GOTO 9999
         ENDIF

         IF( NTDLTE .GT. 0 .AND. NBNOTE0 .NE. NBNOTE ) THEN
            IF(LANGAG .EQ. 0 ) THEN
               PRINT*,'INCOHERENCE NBNOTE=',NBNOTE,' NBNOTE LU=',NBNOTE0
            ELSE
               PRINT*,'INCORRECT NBNOTE=',NBNOTE,' READ NBNOTE=',NBNOTE0
            ENDIF
            IERR = 5
            GOTO 9999
         ENDIF

C        PASSAGE AUX CONDITIONS AUX LIMITES AU TEMPS DU FICHIER LU NOVVIPR
         GOTO 300

      ENDIF


C     INITIALISATION A ZERO DU VECTEUR VP AU TEMPS INITIAL TEMPS0
C     -----------------------------------------------------------
      TEMPS = TEMPS0
      CALL AZEROD( NTDLVP, VXYZP0 )

      IF( NTCALCUL .EQ. 2 ) THEN

C        SUR l'OBJET FLUIDE:

C        CONSTRUCTION DE LA VITESSE INITIALE DU FLUIDE-OBJET
C        A PARTIR DE LA DONNEE DES VALEURS DU TMS VITFLUIN DE L'OBJET
C        ------------------------------------------------------------
         L1 = NUDCNB(KNOMOB)
         NOMTS  = '~>OBJET>' // KNOMOB(1:L1) // '>VITFLUIN'
         KNOMTD = '~>>>VITFLUIN'
         L1 = NUDCNB(NOMTS)
         CALL MOTSTD( KNOMTD, NOMTS(1:L1), IERR )
         IF( IERR .NE. 0 ) THEN
            NBLGRC(NRERR) = 1
            IF( LANGAG .EQ. 0 ) THEN
            KERR(1)='ERREUR: IMPOSSIBLE RECUPERER VITESSE INITIALE'
            ELSE
            KERR(1)='ERROR: IMPOSSIBLE to RETRIEVE the INITIAL VELOCITY'
            ENDIF
            CALL LEREUR
            IERR = 15
            GOTO 9999
         ENDIF

C        CONSTRUCTION DE LA PRESSION INITIALE DU FLUIDE-OBJET
C        A PARTIR DE LA DONNEE DES VALEURS DU TMS PREFLUIN DE L'OBJET
C        ------------------------------------------------------------
         L1 = NUDCNB(KNOMOB)
         NOMTS  = '~>OBJET>' // KNOMOB(1:L1) // '>PREFLUIN'
         KNOMTD = '~>>>PREFLUIN'
         L1 = NUDCNB(NOMTS)
         CALL MOTSTD( KNOMTD, NOMTS(1:L1), IERR )
         IF( IERR .NE. 0 ) THEN
            NBLGRC(NRERR) = 1
            IF( LANGAG .EQ. 0 ) THEN
            KERR(1)='ERREUR: IMPOSSIBLE RECUPERER PRESSION INITIALE'
            ELSE
            KERR(1)='ERROR: IMPOSSIBLE to RETRIEVE the INITIAL PRESSURE'
            ENDIF
            CALL LEREUR
            IERR = 16
            GOTO 9999
         ENDIF

C        INITIALISATION DANS VXYZP0 DES VITESSES INITIALES DE L'OBJET
C        A PARTIR DU TMS VITFLUIN DONNE JUSTE AVANT
C        ------------------------------------------------------------
         CALL LXTSOU( NTLXOB, 'VITFLUIN', NTVTIN, MNVTIN )
         IF( MNVTIN .LE. 0 ) THEN
            IF( LANGAG .EQ. 0 ) THEN
               PRINT *, 'vxyzpt0: PAS de VITESSE INITIALE sur l''OBJET'
            ELSE
               PRINT *, 'NO INITIAL SPEED for the OBJECT'
            ENDIF
            GOTO 50
         ENDIF

C        LE CODE DE CALCUL DE LA VITESSE INITIALE
         LTVIF0 = MCN( MNVTIN + WTVIF0 )
C        LES DIFFERENTES POSSIBILITES D'INITIALISATION DE LA VITESSE INITIALE
         IF( LTVIF0 .EQ. 1 ) THEN

C           VITESSE INITIALE CONSTANTE SUR L'OBJET
C           ......................................
C           TEST SUR LES COMPOSANTES DE LA VITESSE
            NBCOVI = MCN(MNVTIN+WBCOVI)
            DO K=1,NBCOVI

C              TRAITEMENT DE LA COMPOSANTE NOC DES VITESSES INITIALES
               NOC = MCN( MNVTIN + WUCOVI - 1 + K )

C              VITESSE CONSTANTE SUIVANT LA DIRECTION K
               VITES = RMCN( MNVTIN + WUCOVI + NBCOVI - 1 + K )

C              ON ASSIGNE LA VALEUR VITES A TOUS LES NOEUDS SUPPORT
C              DE VITESSE A PARTIR DU TABLEAU NDDLNO
C              ADRESSE DE VX EST L'ADRESSE DU DERNIER DDL
C              DU POINT PRECEDENT+NOC
               DO I=1,NBNOVI
                  INDX = NDDLNO(I-1)+NOC
                  IF( (INDX.LE.0) .OR. (INDX.GT.NTDLVP) ) THEN
                     IERR = 17
                     GOTO 9999
                  ENDIF
                  VXYZP0( INDX ) = VITES
               ENDDO
C
            ENDDO
            GOTO 50

         ELSE IF( LTVIF0 .EQ. -1 ) THEN

C           VITESSE INITIALE DONNEE PAR UNE FONCTION SUR L'OBJET
C           ....................................................
            CALL LXNMNO( NTOBJE, KNOMOB, NOOB, I )
            NBCOVI = MCN(MNVTIN+WBCOVI)
C           LE NUMERO DE LA FONCTION
            NOFONC = MCN(MNVTIN + WUCOVI + NBCOVI)

            DO K=1,NBCOVI

C              TRAITEMENT DES COMPOSANTES DES VITESSES INITIALES
               NOC = MNVTIN + WUCOVI - 1 + K
C              ON SUPPOSE QUE LA FONCTION A EXACTEMENT 7 ARGUMENTS:
C              Temps,X,Y,Z, NoTypeObjet, NoObjet, No Composante
C              OU NOC INDIQUE LE NUMERO DE LA COMPOSANTE DE LA VITESSE TRAITEE
               NBPAFO=7
ccc               write(*,*) 'Numero de la fonction Vitesse0:',NOFONC
               MM = MNXYZN + WYZNOE - 3
               DO I=1,NBNOVI
                  MM   = MM + 3
                  INDX = NDDLNO(I-1)+NOC
                  XPOIN = RMCN(MM  )
                  YPOIN = RMCN(MM+1)
                  ZPOIN = RMCN(MM+2)
                  PARFON(1) = TEMPS0
                  PARFON(2) = XPOIN
                  PARFON(3) = YPOIN
                  PARFON(4) = ZPOIN
                  PARFON(5) = 5
                  PARFON(6) = NOOB
                  PARFON(7) = NOC
                  NCODEV = 1
                  CALL FONVAL( NOFONC, NBPAFO, PARFON, NCODEV,
     %                         VXYZP0(INDX) )
               ENDDO

            ENDDO
            GOTO 50
         ENDIF

C        INITIALISATION DANS VXYZP0 DES PRESSIONS INITIALES
C        DU FLUIDE-OBJET A PARTIR DU TMS PREFLUIN DONNE JUSTE AVANT
C        ----------------------------------------------------------
 50      CALL LXTSOU( NTLXOB, 'PREFLUIN', NTPRIN, MNPRIN )

         IF( MNPRIN .LE. 0 ) THEN
            IF( LANGAG .EQ. 0 ) THEN
               WRITE(IMPRIM,*) 'PAS de PRESSION INITIALE sur l''OBJET'
            ELSE
               WRITE(IMPRIM,*) 'NO INITIAL PRESSURE on the OBJECT'
            ENDIF
            GOTO 200
         ENDIF

C        LE CODE DE CALCUL DE LA PRESSION INITIALE SUR L'OBJET
         LTPRE0 = MCN( MNPRIN + WTPRE0 )
C        LES DIFFERENTES POSSIBILITES D'INITIALISATION DE LA PRESSION SUR L'OBJE
         IF( LTPRE0 .EQ. 1 ) THEN

C           PRESSION INITIALE CONSTANTE SUR L'OBJET
C           .......................................
            PRESSION = RMCN(MNPRIN+WAPRE0)

C           ON ASSIGNE LA VALEUR PRESSION AUX NOEUDS QUI ONT TROIS DDL
C           ON UTILISE LE TABLEAU NDDLNO (TABLEAU PAR NOEUDS)
            DO I=1,NBNOHB
               INDY = NDDLNO(I-1)
               INDX = NDDLNO(I) - INDY
               IF( INDX .EQ. NDIM+1 ) THEN
                  VXYZP0(INDY+NDIM+1) = PRESSION
               ENDIF
            ENDDO
            GOTO 200

         ELSE IF( LTPRE0 .EQ. -1 ) THEN

C           PRESSION INITIALE DONNEE PAR FONCTION SUR L'OBJET
C           .................................................
            CALL LXNMNO( NTOBJE, KNOMOB, NOOB, I )
C           LE NUMERO DE LA FONCTION
            NOFONC = MCN(MNPRIN +WFPRE0)
C           ON SUPPOSE QUE LA FONCTION A EXACTEMENT SIX ARGUMENTS:
C           Temps, X,Y,Z, NoTypeObjet, NoObjet
            NBPAFO=6
            MM = MNXYZN + WYZNOE - 3
            DO  I=1,NBNOHB
               MM = MM + 3
C              S'AGIT-IL d'UN NOEUD avec un DL PRESSION ?
               INDY = NDDLNO(I-1)
               INDX = NDDLNO(I) - INDY
               IF( INDX .EQ. NDIM+1 ) THEN
C                 OUI
                  XPOIN = RMCN(MM  )
                  YPOIN = RMCN(MM+1)
                  ZPOIN = RMCN(MM+2)
                  PARFON(1) = TEMPS0
                  PARFON(2) = XPOIN
                  PARFON(3) = YPOIN
                  PARFON(4) = ZPOIN
                  PARFON(5) = 5
                  PARFON(6) = NOOB
                  CALL FONVAL( NOFONC, NBPAFO, PARFON, NCODEV,
     %                         PRESSION )
                  VXYZP0(INDY+NDIM+1) = PRESSION
               ENDIF
            ENDDO
            GOTO 200
         ENDIF

      ENDIF


C     SUR les VSLP du FLUIDE:
C     MODIFICATION DU VECTEUR VITESSE ET PRESSION VXYZP0
C     A PARTIR DES TMS VITFLUIN ET PREFLUIN DES VSLP EXISTANTS
C     ========================================================
 200  IF( IEVTIN .LE. 0 ) GOTO 250

C     AU MOINS UN VSLP SUPPORTE UN TMS 'VITFLUIN'
C     -------------------------------------------

      WRITE(IMPRIM,*)
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,*) 'VITESSES INITIALES sur VSLP PRISES EN COMPTE'
      ELSE
        WRITE(IMPRIM,*)'INITIAL VELOCITIES on VSLP are TAKEN in ACCOUNT'
      ENDIF

C     LA BOUCLE SUR LES TYPES D'ELEMENTS FINIS
      DO NOTYEL = 1, NBTYEL

C        L'ADRESSE DU TABLEAU NPEF"No TYPE EF
         MNELE = MCN( MNNPEF - 1 + NOTYEL )
C        LE NUMERO DU TYPE DE L'ELEMENT FINI
         NUTYEL = MCN( MNELE + WUTYEL )
C        LE NOMBRE D'ELEMENTS DE CE TYPE
         NBELEM = MCN( MNELE + WBELEM )

C        L'ADRESSE MCN DES NOEUDS ET POINTS GEOMETRIQUES DES EF
         MNNDEL = MNELE + WUNDEL
         MNPGEL = MNNDEL

C        LES CARACTERISTIQUES DE L'ELEMENT FINI
C        ON RETROUVE: NBPOE, NBNOE, NARET
         CALL ELTYCA( NUTYEL )

C        LA BOUCLE SUR LES ELEMENTS FINIS DE CE TYPE NUTYEL
         DO NUELEM = 1, NBELEM

C           LES NOEUDS DE L'ELEMENT FINI
C           LE NUMERO DE VOLUME  DE L'EF
C           LE NUMERO DE SURFACE DES FACES   DE L'EF
C           LE NUMERO DE LIGNE   DES ARETES  DE L'EF
C           LE NUMERO DE POINT   DES SOMMETS DE L'EF
            CALL EFPLSV( MNELE , NUELEM,
     %                   NOVCEL, NOSFEL, NOLAEL, NOPSEL,
     %                   NOOBVC, NOOBSF, NOOBLA, NOOBPS )

C           LE CALCUL DU TYPE OBJET DE CHAQUE NOEUD DE L'ELEMENT FINI
C           EN COMMENCANT PAR LE VOLUME, PUIS LES FACES, PUIS LES ARETES
C           PUIS LES SOMMETS CE QUI ASSURE LA PRIORITE DES POINTS SUR
C           LES LIGNES, SURFACES ET VOLUMES
            CALL EFTNND( NOOBVC, NOOBSF, NOOBLA, NOOBPS,
     %                   NUMIOB, MNDOEL,
     %                  'VITFLUIN', MXDOFL, LPVITF,
     %                   NOTYOB )

C           LE RECENSEMENT DE LA VITESSE INITIALE
C           AUX NBNOE NOEUDS DE CET ELEMENT FINI
            DO J=1,NBNOE

C              LE NUMERO DU NOEUD DANS LE MAILLAGE DE L'OBJET
               NONOE = MCN( MNNDEL-1 + NUELEM + NBELEM*(J-1) )

C              L'ADRESSE MCN DU TABLEAU VITFLUIN
               MNVIIN = NOTYOB(3,J)

C              EXISTE-T-IL UNE "VITFLUIN" EN CE NOEUD ?
               IF( MNVIIN .GT. 0 ) THEN
C                 CALCUL DE LA VITESSE INITIALE EN CE NOEUD
C                 LE TYPE OBJET DU NOEUD J DE L'EF
                  NTYOB = NOTYOB(1,J)
                  NOOB  = NOTYOB(2,J)
                  N = MNXYZN + WYZNOE + 3 * NONOE - 3
                  XPOIN = RMCN(N)
                  YPOIN = RMCN(N+1)
                  ZPOIN = RMCN(N+2)
                  CALL REVITIN( NTYOB,NOOB,XPOIN,YPOIN,ZPOIN,MNVIIN,
     %                          VITESSE )
                  DO K=1,MCN( MNVIIN + WBCOVI )
C                    NUMERO DE LA COMPOSANTE DE LA VITESSE
                     N = MCN( MNVIIN + WUCOVI - 1 + K )
C                    NO DU DL
                     INDX = NDDLNO(NONOE-1)+N
C                    LA VALEUR DE LA VITESSE INITIALISEE
                     VXYZP0(INDX) = VITESSE(N)
                  ENDDO
               ENDIF

            ENDDO

            IF( NUTYEL .EQ. 13 .OR. NUTYEL .EQ. 19 ) THEN

C              TRIANGLE ou TETRAEDRE BREZZI-FORTIN:
C              LA VITESSE INITIALE AU BARYCENTRE EST ICI PRISE EN COMPTE
C              .........................................................
C              LE NUMERO DU NOEUD BARYCENTRE DANS LES NOEUDS VITESSE DU FLUIDE
               NONOE = NBSOM + NUELEM

C              L'ADRESSE MCN DU TABLEAU VITFLUIN EST CELLE
               IF( NUTYEL .EQ. 13 ) THEN
C                 DE LA SURFACE POUR LE TRIANGLE
                  NTYOB = 3
                  NOOB  = NOOBSF(1)
                  MNVIIN = MNDOEL(3)
     %                   + MXDOFL * ( NOOB - NUMIOB(3) ) - 1
     %                   + LPVITF
                  MNVIIN = MCN( MNVIIN )
               ELSE
C                 DU VOLUME POUR LE TETRAEDRE
                  NTYOB = 4
                  NOOB  = NOOBVC
                  MNVIIN = MNDOEL(4)
     %                   + MXDOFL * ( NOOB - NUMIOB(4) ) - 1
     %                   + LPVITF
                  MNVIIN = MCN( MNVIIN )
               ENDIF

C              EXISTE-T-IL UNE "VITFLUIN" EN CE NOEUD BARYCENTRE?
               IF( MNVIIN .GT. 0 ) THEN
C                 OUI: CALCUL DE LA VITESSE INITIALE EN CE NOEUD BARYCENTRE
                  N = MNXYZN + WYZNOE + 3 * NONOE - 3
                  XPOIN = RMCN(N)
                  YPOIN = RMCN(N+1)
                  ZPOIN = RMCN(N+2)
                  CALL REVITIN( NTYOB,NOOB,XPOIN,YPOIN,ZPOIN,MNVIIN,
     %                          VITESSE )
                  DO K=1,MCN( MNVIIN + WBCOVI )
C                    NUMERO DE LA COMPOSANTE DE LA VITESSE
                     N = MCN( MNVIIN + WUCOVI - 1 + K )
C                    NO DU DL
                     INDX = NDDLNO(NONOE-1)+N
C                    LA VALEUR DE LA VITESSE INITIALISEE
                     VXYZP0(INDX) = VITESSE(N)
                  ENDDO
               ENDIF
            ENDIF

         ENDDO
      ENDDO


C     TRAITEMENT DE LA PRESSION INITIALE SUR LES VSLP
C     ===============================================
 250  IF( IEPRIN .GT. 0 ) THEN

C        AU MOINS UN VSLP SUPPORTE UN TMS 'PREFLUIN'
C        -------------------------------------------

         WRITE(IMPRIM,*)
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,*) 'PRESSION INITIALE sur VSLP PRISE EN COMPTE'
         ELSE
          WRITE(IMPRIM,*)'INITIAL PRESSURE on VSLP are TAKEN in ACCOUNT'
         ENDIF

C        LA BOUCLE SUR LES TYPES D'ELEMENTS FINIS
         DO NOTYEL = 1, NBTYEL

C           L'ADRESSE DU TABLEAU NPEF"No TYPE EF
            MNELE = MCN( MNNPEF - 1 + NOTYEL )
C           LE NUMERO DU TYPE DE L'ELEMENT FINI
            NUTYEL = MCN( MNELE + WUTYEL )
C           LE NOMBRE D'ELEMENTS DE CE TYPE
            NBELEM = MCN( MNELE + WBELEM )

C           L'ADRESSE MCN DES NOEUDS ET POINTS GEOMETRIQUES DES EF
            MNNDEL = MNELE + WUNDEL
            MNPGEL = MNNDEL

C           LES CARACTERISTIQUES DE L'ELEMENT FINI
C           ON TROUVE: NBPOE, NBNOE, NARET
            CALL ELTYCA( NUTYEL )

C           LA BOUCLE SUR LES ELEMENTS FINIS DE CE TYPE NUTYEL
C           ==================================================
            DO NUELEM = 1, NBELEM

C              LES NOEUDS DE L'ELEMENT FINI
C              LE NUMERO DE VOLUME  DE L'EF
C              LE NUMERO DE SURFACE DES FACES   DE L'EF
C              LE NUMERO DE LIGNE   DES ARETES  DE L'EF
C              LE NUMERO DE POINT   DES SOMMETS DE L'EF
               CALL EFPLSV( MNELE , NUELEM,
     %                      NOVCEL, NOSFEL, NOLAEL, NOPSEL,
     %                      NOOBVC, NOOBSF, NOOBLA, NOOBPS )

C              LE CALCUL DU TYPE OBJET DE CHAQUE NOEUD DE L'ELEMENT FINI
C              EN COMMENCANT PAR LE VOLUME, PUIS LES FACES, PUIS LES ARETES
C              PUIS LES SOMMETS CE QUI ASSURE LA PRIORITE DES POINTS SUR
C              LES LIGNES, SURFACES ET VOLUMES
               CALL EFTNND( NOOBVC, NOOBSF, NOOBLA, NOOBPS,
     %                      NUMIOB, MNDOEL,
     %                     'PREFLUIN', MXDOFL, LPPREF,
     %                      NOTYOB )

C              LE RECENSEMENT DE LA PRESSION INITIALE
C              AUX NBNSOM SOMMETS DE CET ELEMENT FINI TAYLOR-HOOD ou BREZZI-FORT
C              CAR LA PRESSION Y EST TOUJOURS UNE P1-INTERPOLATION
C              ET LES SOMMETS SONT TOUJOURS LES PREMIERS NOEUDS DE L'EF
               DO J=1,NBNSOM

C                 LE NUMERO DU NOEUD DANS LE MAILLAGE DE L'OBJET
                  NONOE = MCN( MNNDEL-1 + NUELEM + NBELEM*(J-1) )
C
C                 L'ADRESSE MCN DU TABLEAU PREFLUIN
                  MNPRIN = NOTYOB(3,J)
C
C                 EXISTE-T-IL UNE "PREFLUIN" EN CE NOEUD ?
                  IF( MNPRIN .GT. 0 ) THEN
C                    OUI: MAIS TRAITER SEULEMENT LES SOMMETS
                     INDX = NDDLNO(NONOE-1)
                     INDY = NDDLNO(NONOE)
C                    POUR LES NOEUDS PORTANT UN DDL EN PRESSION SEULEMENT
                     IF( INDY-INDX .EQ. NDIM+1 ) THEN
C
C                       CALCUL DE LA PRESSION INITIALE EN CE NOEUD
C                       LE TYPE OBJET DU NOEUD J DE L'EF
                        NTYOB = NOTYOB(1,J)
                        NOOB  = NOTYOB(2,J)
                        N = MNXYZN + WYZNOE + 3 * NONOE - 3
                        XPOIN = RMCN(N)
                        YPOIN = RMCN(N+1)
                        ZPOIN = RMCN(N+2)
                        CALL REPREIN( NTYOB,NOOB,XPOIN,YPOIN,ZPOIN,
     %                                MNPRIN, PRESSION )
C                       LA VALEUR DE LA PRESSION FIXEE DERNIER DL DU SOMMET
                        VXYZP0(INDY) = PRESSION

                     ENDIF
                  ENDIF

               ENDDO
            ENDDO
         ENDDO
      ENDIF


C     CONSTRUCTION DES TABLEAUX DES NO ET DES VALEURS VITESSES-PRESSIONS FIXEES
C     A L'INSTANT INITIAL TEMPS0 POUR IMPOSER LES CONDITIONS AUX LIMITES
C     =========================================================================
 300  CALL VIPRFX( RELMIN, NTDLHB, NDIM,    MNXYZN, NDDLNO,
     %             NBTYEL, MNNPEF, NUMIOB,  MNDOEL,
     %             NBVCFX, NBVPFX, MNNVPFX, MNVVPFX )

      IF( NBVPFX .LE. 0 ) THEN
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
         GOTO 9999
      ENDIF


C     PRISE EN COMPTE DES CONDITIONS AUX LIMITES VITESSES-PRESSIONS FIXEES
C     SUR LE VECTEUR INITIAL VXYZP0
C     ====================================================================
      IF( NOVVIPR .EQ. 0 ) THEN
C        PAS DE FICHIER RECUPERE AU TEMPS INITIAL
         MOREE2 = MOTVAR(6)
         DO I=1, NBVPFX
C           LE NO DU DL FIXE
            NDL = MCN( MNNVPFX - 1 + I )
C           LA VALEUR A FIXER
            DVAL = DMCN( (MNVVPFX-1)/MOREE2 + I )
C           MODIFICATION DU VECTEUR VXYZP0
            IF( DVAL .NE. 1D222 ) THEN
C              VITESSE NON CONVECTEE AU TEMPS TEMPS0
               VXYZP0( NDL ) = DVAL
            ENDIF
         ENDDO
      ENDIF
  
      IF( NAVSTO .EQ. 3 .AND. NOVVIPR .EQ. 0 ) THEN

C        PAS DE FICHIER VITESSE+PRESSION+TEMPERATURE RECUPERE AU TEMPS INITIAL
C        INITIALISATION DU VECTEUR TEMPER0 DES TEMPERATURES INITIALES
C        AUX NOEUDS DU MAILLAGE ET A L'INSTANT TEMPS INITIAL TEMPS0
C        =====================================================================
C        TEMPER0 EST LE VECTEUR GLOBAL DES TEMPERATURES INITIALES
         CALL TEMPERT0( KNOMOB, NTLXOB, MOREE2, RELMIN, NDIM,
     %                  NTDLTE, TEMPS0, IETEIN,
     %                  NBTYEL, MNNPEF, NDPGST,
     %                  MNXYZN, NUMIOB, MNDTEL,
     %                  NTVECT, MNVECT, NBVECT,
     %                  NBTEFX, MONTEFX, MNNTEFX, MNVTEFX,
     %                  TEMPER0, IERR )

      ENDIF

 9999 IF( LANGAG .EQ. 0 ) THEN
         PRINT*,'vxyzpt0: Fin calcul du VECTEUR Vitesse+Pression+Tempera
     %ture du TEMPS',TEMPS0,' option NAVSTO=',NAVSTO,
     %' Nb DL VP Fixes=',NBVPFX
      ELSE
         PRINT*,'vxyzpt0: End computation of Velocity+Pressure+Temperatu
     %re VECTOR at TIME',TEMPS0,' option NAVSTO=',NAVSTO,
     %' Nb FIXED VP DoF=',NBVPFX
      ENDIF

      RETURN
      END
