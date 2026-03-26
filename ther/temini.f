      SUBROUTINE TEMINI( KNOMOB, NTLXOB, MOREE2, NTDL,   TPSINI, IETEIN,
     %                   NBTYEL, MNNPEF, NDPGST,
     %                   MNXYZN, NUMIOB, MNDOEL,
     %                   NTVECT, MNVECT, NBTEMP, MNUG0,  IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CONSTRUIRE LE VECTEUR DE LA TEMPERATURE INITIALE
C -----    AUX NOEUDS DU MAILLAGE DE L'OBJET THERMIQUE
C          LA TEMPERATURE INITIALE A ETE DONNEE EVENTUELLEMENT SUR
C          L'OBJET soit deja CALCULEE, soit CONSTANTE, soit FONCTION
C          LES VOLUMES  soit Constante, soit Fonction
C          LES SURFACES soit Constante, soit Fonction
C          LES LIGNES   soit Constante, soit Fonction
C          LES POINTS   soit Constante, soit Fonction
C
C ENTREES:
C --------
C KNOMOB : NOM DE L'OBJET A TRAITER
C NTLXOB : NUMERO DU TMS DU LEXIQUE DE L'OBJET KNOMOB
C MOREE2 : NOMBRE DE MOTS   D'UNE VARIABLE REELLE DOUBLE PRECISION
C NTDL   : NOMBRE TOTAL DE DEGRES DE LIBERTE SUR LE MAILLAGE
C          EGAL AU NOMBRE DE NOEUDS DU MAILLAGE
C
C TPSINI : TEMPS INITIAL DU CALCUL (PEUT ETRE MODIFIE!)
C IETEIN :=0 SI AUCUN PLSV DE L'OBJET SUPPORTE UN TMS TEMPERINIT
C         >0 SINON
C
C NBTYEL : NOMBRE DE TYPES D'EF DE L'OBJET
C MNNPEF : ADRESSE MCN DES TABLEAUX NPEF"TYPE EF
C NDPGST : CODE DE STOCKAGE DES SOMMETS POINTS NOEUDS
C
C MNXYZN : ADRESSE MCN DU TMS XYZNOEUD DE L'OBJET
C NUMIOB : NUMERO MINIMAL DES PLSV
C MNDOEL : ADRESSE MCN DES DONNEES THERMIQUES DE L'OBJET
C
C SORTIES:
C --------
C NTVECT : NUMERO      DU TMS VECTEUR"TEMPERATURE DE L'OBJET
C MNVECT : ADRESSE MCN DU TMS VECTEUR"TEMPERATURE DE L'OBJET
C NBTEMP : NUMERO DU DERNIER VECTEUR TEMPERATURE STOCKE = VECTEUR A TPSINI
C MNUG0  : ADRESSE MCN DE UG0 TEMPERATURE A L'INSTANT INITIAL
C IERR   : 0 SI PAS D'ERREUR RENCONTREE, 1 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS       AVRIL 1999
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/gsmenu.inc"
      include"./incl/a___temperinit.inc"
      include"./incl/a___vecteur.inc"
      include"./incl/a___npef.inc"
      include"./incl/a___xyznoeud.inc"
      include"./incl/ponoel.inc"
      include"./incl/donele.inc"
      include"./incl/donthe.inc"
      include"./incl/ctemps.inc"
      include"./incl/cthet.inc"
      include"./incl/cnonlin.inc"
C
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE      (MCN(1), RMCN(1), DMCN(1))
C
      CHARACTER*(*)     KNOMOB
C
      DOUBLE PRECISION  TEMPE0
      DOUBLE PRECISION  XYZP(3)
      INTEGER           NOTYOB(3,MXNOEL)
      INTEGER           NOOBSF(6),NOOBLA(12),NOOBPS(8)
      INTEGER           NUMIOB(4),MNDOEL(4)
C
C     INITIALISATION DU VECTEUR UG0 DES TEMPERATURES INITIALES
C     AUX NOEUDS DU MAILLAGE DE L'OBJET ET A L'INSTANT TEMPS INITIAL
C     ==============================================================
C     LE VECTEUR GLOBAL DES TEMPERATURES INITIALES EST DECLARE
      IF( MNUG0 .LE. 0 ) THEN
         CALL TNMCDC( 'REEL2', NTDL, MNUG0 )
         CALL AZEROD( NTDL, MCN(MNUG0) )
      ENDIF
C
C     ---------------------------------------------
C     TABLEAU DES TEMPERATURES INITIALES DE L'OBJET
C     ---------------------------------------------
      NBTEMP = 0
      CALL LXTSOU( NTLXOB, 'TEMPERINIT', NTTEIN, MNTEIN )
      IF( MNTEIN .LE. 0 ) GOTO 80
C
C     ICI LE TABLEAU DES TEMPERATURES INITIALES DE L'OBJET EXISTE
C     LE CODE DE CALCUL DE LA TEMPERATURE INITIALE
      LTTEM0 = MCN( MNTEIN + WTTEM0 )
      IF( LTTEM0 .EQ. 1 ) THEN
C
C        TEMPERATURE INITIALE CONSTANTE
C        ------------------------------
         TEMPE0 = RMCN( MNTEIN + WATEM0 )
         MN     = ( MNUG0 - 1 ) / MOREE2
         DO 10 I=1,NTDL
            DMCN( MN + I ) = TEMPE0
 10      CONTINUE
         GOTO 80
C
      ELSE IF( LTTEM0 .EQ. -1 ) THEN
C
C        LA TEMPERATURE EST OBTENUE PAR FONCTION
C        ---------------------------------------
         CALL LXNMNO( NTOBJE, KNOMOB, NOOB, I )
         MN = (MNUG0-1) / MOREE2
         MM = MNXYZN + WYZNOE - 3
         DO 20 I=1,NTDL
            MM = MM + 3
            XYZP(1) = RMCN(MM  )
            XYZP(2) = RMCN(MM+1)
            XYZP(3) = RMCN(MM+2)
            CALL RETEIN( 5,NOOB, 3,XYZP, MNTEIN, TEMPE0 )
            DMCN( MN + I ) = TEMPE0
 20      CONTINUE
ccc
ccc      print *,'TABLEAU temperature initiale='
ccc      MN = (MNUG0-1) / MOREE2
ccc      mm=MNXYZN+wyznoe
ccc      do 222 k=1,ntdl
ccc         print *,'XYZ=',(RMCN(mm+kk),kk=0,2),
ccc     %           ' UG00(',k,')=',DMCN(MN+k)
ccc         mm = mm + 3
ccc 222  continue
ccc
         GOTO 80
C
      ELSE
C
C        LA TEMPERATURE INITIALE EST STOCKEE DANS L'ACTUEL (LTTEM0=2)
C        ~>OBJET>NomObjet>VECTEUR"TEMPERATURE
C        ------------------------------------------------------------
         CALL LXTSOU( NTLXOB, 'VECTEUR"TEMPERATURE', NTVECT, MNVECT )
         IF( MNVECT .LE. 0 ) THEN
            NBLGRC(NRERR) = 1
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'temini: PAS DE VECTEUR"TEMPERATURE INITIALE'
            ELSE
               KERR(1) = 'temini: NO VECTOR INITIAL TEMPERATURE'
            ENDIF
            CALL LEREUR
            IERR = 1
            GOTO 9999
         ENDIF
         IF( NTDL .NE. MCN(MNVECT+WBCOVE) ) THEN
            NBLGRC(NRERR) = 2
            WRITE(KERR(5)(1:8),'(I8)') NTDL
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'NOMBRE DE COMPOSANTES DU VECTEUR"TEMPERATURE'
               KERR(2) = 'DIFFERENT DE ' // KERR(5)(1:8) // ' NOEUDS'
            ELSE
               KERR(1) = 'NUMBER of COMPONENTS of VECTOR"TEMPERATURE'
               KERR(2) = 'NOT EQUAL TO ' // KERR(5)(1:8) // ' NODES'
            ENDIF
            CALL LEREUR
            IERR = 1
            GOTO 9999
         ENDIF
C
C        RECHERCHE DU TEMPS INITIAL TPSINI COMME PLUS PROCHE TEMPS
C        PARMI LES TEMPS STOCKES DU CALCUL PRECEDENT
C        NOMBRE DE TEMPS STOCKES
         N      = MCN( MNVECT + WBCPIN )
         MNTIME = MNVECT + WECTEU + MOREE2*NTDL*MCN(MNVECT+WBVECT) - 1
C        NBTEMP LE NUMERO DU VECTEUR STOCKE AU TEMPS LE PLUS PROCHE DE TPSINI
         NBTEMP = 1
         R      = RINFO( 'GRAND' )
         DO 30 I=1,N
            RR = ABS( RMCN(MNTIME+I) - TPSINI )
            IF( RR .LT. R ) THEN
               NBTEMP = I
               R      = RR
            ENDIF
 30      CONTINUE
C
         MNTEMP = MNVECT + WECTEU + MOREE2 * NTDL * (NBTEMP-1)
C        TRANSFERT DE LA TEMPERATURE INITIALE NBTEMP DANS LE TABLEAU MNUG0
         CALL TRTATD( MCN(MNTEMP), MCN(MNUG0), NTDL )
C
C        AFFICHAGE DU NOUVEAU TEMPS INITIAL RETENU
         TPSINI   = RMCN(MNTIME+NBTEMP)
         TEMPS    = TPSINI
         TEMPSINI = TPSINI
C
         IF( LANGAG .EQ. 0 ) THEN
           WRITE(IMPRIM,*)'TEMPERATURE INITIALE au TEMPS INITIAL=',TEMPS
         ELSE
            WRITE(IMPRIM,*) 'INITIAL TEMPERATURE at INITIAL TIME=',TEMPS
         ENDIF
C        INITIALISATION TERMINEE
         GOTO 9999
C
      ENDIF
C
C     -----------------------------------------------------------------
C     TRAITEMENT DE LA TEMPERATURE INITIALE PARTICULARISEE SUR DES PLSV
C     -----------------------------------------------------------------
 80   IF( IETEIN .LE. 0 ) GOTO 9999
C     AU MOINS UN PLSV SUPPORTE UN TMS 'TEMPERINIT'
C
C     LA BOUCLE SUR LES TYPES D'ELEMENTS FINIS
C     ========================================
      MM0 = MNXYZN + WYZNOE - 3
      MN0 = (MNUG0-1) / MOREE2
      DO 100 NOTYEL = 1, NBTYEL
C
C        L'ADRESSE DU TABLEAU NPEF"No TYPE EF
         MNELE = MCN( MNNPEF - 1 + NOTYEL )
C        LE NUMERO DU TYPE DE L'ELEMENT FINI
         NUTYEL = MCN( MNELE + WUTYEL )
C        LE NOMBRE D'ELEMENTS DE CE TYPE
         NBELEM = MCN( MNELE + WBELEM )
C
C        L'ADRESSE MCN DES NOEUDS ET POINTS GEOMETRIQUES DES EF
         MNNDEL = MNELE + WUNDEL
         MNPGEL = MNNDEL
         IF( NDPGST .GE. 2 ) THEN
            MNPGEL = MNPGEL + NBELEM * MCN(MNELE+WBNDEL)
         ENDIF
C
C        LES CARACTERISTIQUES DE L'ELEMENT FINI
C        ON TROUVE: NBPOE, NBNOE, NARET
         CALL ELTYCA( NUTYEL )
C
C        LA BOUCLE SUR LES ELEMENTS FINIS DE CE TYPE NUTYEL
C        ==================================================
         DO 95 NUELEM = 1, NBELEM
C
C           LES NOEUDS DE L'ELEMENT FINI
C           LE NUMERO DE VOLUME  DE L'EF
C           LE NUMERO DE SURFACE DES FACES   DE L'EF
C           LE NUMERO DE LIGNE   DES ARETES  DE L'EF
C           LE NUMERO DE POINT   DES SOMMETS DE L'EF
            CALL EFPLSV( MNELE , NUELEM,
     %                   NOVCEL, NOSFEL, NOLAEL, NOPSEL,
     %                   NOOBVC, NOOBSF, NOOBLA, NOOBPS )
C
C           LE CALCUL DU TYPE OBJET DE CHAQUE NOEUD DE L'ELEMENT FINI
C           EN COMMENCANT PAR LE VOLUME, PUIS LES FACES, PUIS LES ARETES
C           PUIS LES SOMMETS CE QUI ASSURE LA PRIORITE DES POINTS SUR
C           LES LIGNES, SURFACES ET VOLUMES
            CALL EFTNND( NOOBVC, NOOBSF, NOOBLA, NOOBPS,
     %                   NUMIOB, MNDOEL,
     %                  'TEMPERINIT', MXDOTH, LPTEIN,
     %                   NOTYOB )
C
C           LE RECENSEMENT DE LA TEMPERATURE INITIALE
C           AUX NBNOE NOEUDS DE CET ELEMENT FINI
            DO 90 J=1,NBNOE
C
C              L'ADRESSE MCN DU TABLEAU TEMPERINIT
               MN1 = NOTYOB(3,J)
C
C              EXISTE-T-IL UNE "TEMPERINIT" EN CE NOEUD ?
               IF( MN1 .GT. 0 ) THEN
C
C                 CALCUL DE LA TEMPERATURE INITIALE EN CE NOEUD
C                 LE TYPE OBJET DU NOEUD J DE L'EF
                  NTYOB = NOTYOB(1,J)
                  NOOB  = NOTYOB(2,J)
C                 LE NUMERO DU NOEUD DANS LE MAILLAGE DE L'OBJET
                  NONOE = MCN( MNNDEL-1 + NUELEM + NBELEM*(J-1) )
                  MM    = MM0 + 3 * NONOE
                  XYZP(1) = RMCN(MM)
                  XYZP(2) = RMCN(MM+1)
                  XYZP(3) = RMCN(MM+2)
                  CALL RETEIN( NTYOB,NOOB, 3, XYZP,
     %                         MN1, DMCN(MN0+NONOE) )
C
               ENDIF
C
 90         CONTINUE
 95      CONTINUE
100   CONTINUE
C
9999  RETURN
ccc
ccc      print *,'TABLEAU de la TEMPERATURE INITIALE='
ccc      mnn=MNXYZN+wyznoe
ccc      MN = (MNUG0-1) / MOREE2
ccc      do 200 k=1,ntdl
ccc         print *,'XYZ=',(RMCN(mnn+kk),kk=0,2),
ccc     %           ' UG01(',k,')=',DMCN(MN+k)
ccc         mnn = mnn + 3
ccc 200  continue
ccc
      END
