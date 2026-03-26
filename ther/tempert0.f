      SUBROUTINE TEMPERT0( KNOMOB, NTLXOB, MOREE2, RELMIN, NDIM,
     %                     NTDLTE, TPSINI, IETEIN,
     %                     NBTYEL, MNNPEF, NDPGST,
     %                     MNXYZN, NUMIOB, MNDTEL,
     %                     NTVECT, MNVECT, NBTEMP,
     %                     NBTEFX, MONTEFX, MNNTEFX, MNVTEFX,
     %                     TEMPER0, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CONSTRUIRE LE VECTEUR DE LA TEMPERATURE INITIALE
C -----    AUX NOEUDS DU MAILLAGE DE L'OBJET FLUIDE THERMIQUE
C          LA TEMPERATURE INITIALE A ETE DONNEE EVENTUELLEMENT SUR
C          L'OBJET soit deja CALCULEE, soit CONSTANTE, soit FONCTION
C          LES VOLUMES  soit Constante, soit Fonction
C          LES SURFACES soit Constante, soit Fonction
C          LES LIGNES   soit Constante, soit Fonction
C          LES POINTS   soit Constante, soit Fonction

C ENTREES:
C --------
C KNOMOB : NOM DE L'OBJET A TRAITER
C NTLXOB : NUMERO DU TMS DU LEXIQUE DE L'OBJET KNOMOB
C MOREE2 : NOMBRE DE MOTS   D'UNE VARIABLE REELLE DOUBLE PRECISION
C RELMIN : REEL DOUBLE PRECISION TEMOIN DE NON-INITIALISATION
C NDIM   : DIMENSION DE L'ESPACE DE L'OBJET 2 ou 3
C NTDLTE : NOMBRE TOTAL DE DEGRES DE LIBERTE SUR LE MAILLAGE
C          EGAL AU NOMBRE DE NOEUDS DU MAILLAGE

C TPSINI : TEMPS INITIAL DU CALCUL (PEUT ETRE MODIFIE!)
C IETEIN :=0 SI AUCUN PLSV DE L'OBJET SUPPORTE UN TMS TEMPERINIT
C         >0 SINON
C
C NBTYEL : NOMBRE DE TYPES D'EF DE L'OBJET
C MNNPEF : ADRESSE MCN DES TABLEAUX NPEF"TYPE EF
C NDPGST : CODE DE STOCKAGE DES SOMMETS POINTS NOEUDS
C MNXYZN : ADRESSE MCN DU TMS XYZNOEUD DE L'OBJET
C NUMIOB : NUMERO MINIMAL DES PLSV
C MNDTEL : ADRESSE MCN DES DONNEES THERMIQUES DE L'OBJET FLUIDE

C SORTIES:
C --------
C NTVECT : NUMERO      DU TMS VECTEUR"TEMPERATURE DE L'OBJET
C MNVECT : ADRESSE MCN DU TMS VECTEUR"TEMPERATURE DE L'OBJET
C NBTEMP : NUMERO DU DERNIER VECTEUR TEMPERATURE STOCKE = VECTEUR A TPSINI

C NBTEFX : NOMBRE DE NOEUDS DE TEMPERATURE FIXEE
C MONTEFX: NOMBRE DE MOTS DECLARES DU TABLEAU MC NO DES TEMPERATURE FIXES
C MNNTEFX: ADRESSE MCN DU TABLEAU MC DES NUMEROS DES TEMPERATURE FIXES
C          =0 SINON
C MNVTEFX: ADRESSE MCN DU TABLEAU MC DES VALEURS DES TEMPERATURE FIXES
C          =0 SINON

C TEMPER0: VECTEUR TEMPERATURE A L'INSTANT INITIAL
C IERR   : 0 SI PAS D'ERREUR RENCONTREE, 1 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE UPMC PARIS      Avril 1999
C MODIFS : ALAIN PERRONNET  Saint Pierre du Perray            Avril 2022
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

      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE      (MCN(1), RMCN(1), DMCN(1))

      CHARACTER*(*)     KNOMOB

      DOUBLE PRECISION  TEMPE0, TEMPER0(NTDLTE), XYZP(3), RELMIN,
     %                  TECMOY, TECMIN, TECMAX, TEXMIN, TEXMAX
      INTEGER           NOTYOB(3,MXNOEL)
      INTEGER           NOOBSF(6), NOOBLA(12), NOOBPS(8)
      INTEGER           NUMIOB(4), MNDTEL(4)

      PRINT*
      PRINT*,'tempert0: INITIALISATION DU VECTEUR TEMPER0 DES TEMPERATUR
     %ES INITIALES DES NOEUDS'

C     INITIALISATION DU VECTEUR TEMPER0 DES TEMPERATURES INITIALES
C     AUX NOEUDS DU MAILLAGE DE L'OBJET ET A L'INSTANT TEMPS INITIAL
C     ==============================================================
C     LE VECTEUR GLOBAL DES TEMPERATURES INITIALES EST ANNULE
      CALL AZEROD( NTDLTE, TEMPER0 )

C     --------------------------------------------------------
C     TABLEAU DES TEMPERATURES INITIALES AUX NOEUDS DE L'OBJET
C     --------------------------------------------------------
      NBTEMP = 0
      CALL LXTSOU( NTLXOB, 'TEMPERINIT', NTTEIN, MNTEIN )
      IF( MNTEIN .LE. 0 ) GOTO 80

C     LE CODE DE CALCUL DE LA TEMPERATURE INITIALE
      LTTEM0 = MCN( MNTEIN + WTTEM0 )
      IF( LTTEM0 .EQ. 1 ) THEN
C
C        TEMPERATURE INITIALE CONSTANTE
C        ------------------------------
         TEMPE0 = RMCN( MNTEIN + WATEM0 )
         DO I=1,NTDLTE
            TEMPER0( I ) = TEMPE0
         ENDDO
         GOTO 80
C
      ELSE IF( LTTEM0 .EQ. -1 ) THEN
C
C        LA TEMPERATURE EST OBTENUE PAR FONCTION
C        ---------------------------------------
         CALL LXNMNO( NTOBJE, KNOMOB, NOOB, I )
         MM = MNXYZN + WYZNOE - 3
         DO I=1,NTDLTE
            MM = MM + 3
            XYZP(1) = RMCN(MM  )
            XYZP(2) = RMCN(MM+1)
            XYZP(3) = RMCN(MM+2)
            CALL RETEIN( 5,NOOB, 3,XYZP, MNTEIN, TEMPE0 )
            TEMPER0( I ) = TEMPE0
         ENDDO

ccc      print *,'TABLEAU temperature initiale='
ccc      mm=MNXYZN+wyznoe
ccc      do k=1,ntdlte
ccc         print *,'XYZ=',(RMCN(mm+kk),kk=0,2),
ccc     %           ' Temper0(',k,')=',TEMPER0(k)
ccc         mm = mm + 3
ccc      enddo

         GOTO 80

      ELSE

C        LA TEMPERATURE INITIALE EST STOCKEE DANS L'ACTUEL (LTTEM0=2)
C        ~>OBJET>NomObjet>VECTEUR"TEMPERATURE
C        ------------------------------------------------------------
         CALL LXTSOU( NTLXOB, 'VECTEUR"TEMPERATURE', NTVECT, MNVECT )
         IF( MNVECT .LE. 0 ) THEN
            NBLGRC(NRERR) = 1
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'tempert0: PAS DE VECTEUR"TEMPERATURE INITIALE'
            ELSE
               KERR(1) = 'tempert0: NO VECTOR INITIAL TEMPERATURE'
            ENDIF
            CALL LEREUR
            IERR = 1
            GOTO 9999
         ENDIF
         IF( NTDLTE .NE. MCN(MNVECT+WBCOVE) ) THEN
            NBLGRC(NRERR) = 2
            WRITE(KERR(5)(1:8),'(I8)') NTDLTE
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
         MNTIME = MNVECT + WECTEU + MOREE2*NTDLTE*MCN(MNVECT+WBVECT) - 1
C        NBTEMP LE NUMERO DU VECTEUR STOCKE AU TEMPS LE PLUS PROCHE DE TPSINI
         NBTEMP = 1
         R      = RINFO( 'GRAND' )
         DO I=1,N
            RR = ABS( RMCN(MNTIME+I) - TPSINI )
            IF( RR .LT. R ) THEN
               NBTEMP = I
               R      = RR
            ENDIF
         ENDDO
C
         MNTEMP = MNVECT + WECTEU + MOREE2 * NTDLTE * (NBTEMP-1)
C        TRANSFERT DE LA TEMPERATURE INITIALE NBTEMP DANS LE TABLEAU TEMPER0
         CALL TRTATD( MCN(MNTEMP), TEMPER0, NTDLTE )
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


C     -----------------------------------------------------------------
C     TRAITEMENT DE LA TEMPERATURE INITIALE PARTICULARISEE SUR DES PLSV
C     -----------------------------------------------------------------
 80   IF( IETEIN .LE. 0 ) GOTO 9999
C     AU MOINS UN PLSV SUPPORTE UN TMS 'TEMPERINIT'
C
C     LA BOUCLE SUR LES TYPES D'ELEMENTS FINIS
C     ========================================
      MM0 = MNXYZN + WYZNOE - 3
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
     %                   NUMIOB, MNDTEL,
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
     %                         MN1, TEMPER0(NONOE) )

               ENDIF

 90         ENDDO
 95      ENDDO
100   ENDDO


C     PRISE EN COMPTE DES CONDITIONS AUX LIMITES AU TEMPS INITIAL
C     ===========================================================
C     LISTE DES NUMEROS DES TEMPERATURES FIXEES
C     -----------------------------------------
      CALL THDLFX(      1, NTDLTE, NDIM,
     %             NBTYEL, MNNPEF, NDPGST,
     %             MNXYZN, NUMIOB, MNDTEL, RELMIN,
     %             NBTEFX, MONTEFX, MNNTEFX, MNVTEFX,  IERR )

C     PRISE EN COMPTE DES CONDITIONS AUX LIMITES de TEMPERATURE
C     ---------------------------------------------------------
      CALL BLDLFX( NTDLTE, 1, NBTEFX, MCN(MNNTEFX), MCN(MNVTEFX),
     %             1D0, TEMPER0 )

      print *,'tempert0: TABLEAU de la TEMPERATURE INITIALE TEMPER0'

C     AFFICHAGE DU TABLEAU TEMPER0
C     ----------------------------
      CALL AFTEMP( 3,   1, MNXYZN, NTDLTE, 1,      TEMPER0,
     %             TECMOY, TECMIN, IDLMIN, TECMAX, IDLMAX,
     %             NOFOTI, TEXMIN, TEXMAX )

9999  RETURN
      END
