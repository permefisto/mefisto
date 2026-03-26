      SUBROUTINE PRESFX( RELMIN, NTDLPR, NDIM,   MNXYZN, MNNDDL,
     %                   NBTYEL, MNNPEF, NUMIOB, MNDOEL,
     %                   NBRDLX, MNNDLX, MNVDLX )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    RECENSER LES NUMEROS ET VALEURS DES DEGRES DE LIBERTE FIXES
C -----    EN PRESSION du FLUIDE
C
C ENTREES:
C --------
C RELMIN : REEL DOUBLE PRECISION TEMOIN DE NON-INITIALISATION
C NTDLPR : NOMBRE TOTAL DE DEGRES DE LIBERTE DES PRESSIONS
C NDIM   : NOMBRE DE COMPOSANTES DE LA VITESSE=DIMENSION DE L'ESPACE
C MNXYZN : ADRESSE MCN DU TMS XYZNOEUD
C MNNDDL : ADRESSE DU TABLEAU POINTEURS SUR LE DERNIER DL DE CHAQUE NOEUD
C NBTYEL : NOMBRE DE TYPES D'EF DE L'OBJET
C MNNPEF : ADRESSE MCN DES TABLEAUX NPEF"TYPE EF
C NUMIOB : NUMERO MINIMAL DES OBJETS
C NUMAOB : NUMERO MAXIMAL DES OBJETS
C MNDOEL : ADRESSE MCN DES DONNEES DE L'OBJET
C
C SORTIES:
C --------
C NBRDLX : NOMBRE DE DL PRESSION FIXEES
C          ERREUR SI NBRDLX=0 CAR PAS DE CONDITION AUX LIMITES!
C MNNDLX : ADRESSE MCN DU TABLEAU MC DES NUMEROS DES DL FIXES
C MNVDLX : ADRESSE MCN DU TABLEAU MC DES VALEURS DES PRESSIONS FIXEES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY Novembre 2010
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/donflu.inc"
      include"./incl/donele.inc"
      include"./incl/a___npef.inc"
      include"./incl/a___xyznoeud.inc"
      include"./incl/ponoel.inc"
      include"./incl/gsmenu.inc"
C
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE      (MCN(1), RMCN(1), DMCN(1))
C
      INTEGER           NOTYOB(3,MXNOEL)
      INTEGER           NOOBSF(6), NOOBLA(12), NOOBPS(8)
      INTEGER           NUMIOB(4), MNDOEL(4)
C
      DOUBLE PRECISION  RELMIN, XPOIN, YPOIN, ZPOIN, PRESSION
C
C     POUR EVITER DE DOUBLER LES TABLEAUX
      IF( NBRDLX .GT. 0 .AND. MNNDLX .GT. 0 ) THEN
         CALL TNMCDS( 'ENTIER', NBRDLX, MNNDLX )
         CALL TNMCDS( 'REEL2',  NBRDLX, MNVDLX )
      ENDIF
C
      NBRDLX = 0
      MNNDLX = 0
      MNVDLX = 0
C
C     LES NUMEROS DES DEGRES DE LIBERTE FIXES
      CALL TNMCDC( 'ENTIER', NTDLPR, MNNDLX )
      CALL AZEROI( NTDLPR, MCN(MNNDLX) )
C
C     LES TABLEAUX REELS DES VALEURS DES DL FIXES
C     VALEUR DE NON INITIALISATION = RELMIN
      CALL TNMCDC( 'REEL2', NTDLPR, MNVDLX )
      MNVDX  = ( MNVDLX - 1 ) / 2
      DO 20 I=1,NTDLPR
         DMCN(MNVDX + I) = RELMIN
 20   CONTINUE
C
C     ========================================
C     LA BOUCLE SUR LES TYPES D'ELEMENTS FINIS
C     ========================================
      DO 100 NOTYEL = 1, NBTYEL
C
C        L'ADRESSE DU TABLEAU ELEMENTS (NPEF)
         MNELE = MCN( MNNPEF - 1 + NOTYEL )
C        LE NUMERO DU TYPE DE L'ELEMENT FINI
         NUTYEL = MCN( MNELE + WUTYEL )
C        LE NOMBRE D'ELEMENTS DE CE TYPE
         NBELEM = MCN( MNELE + WBELEM )
C        L'ADRESSE MCN DES NOEUDS ET POINTS GEOMETRIQUES DES ELEMENTS
         MNNDEL = MNELE + WUNDEL
C
C        LES CARACTERISTIQUES DE L'ELEMENT FINI
C        ON TROUVE: NBPOE, NBNOE, NARET
         CALL ELTYCA( NUTYEL )
C
C        LA BOUCLE SUR LES ELEMENTS FINIS DE CE TYPE NUTYEL
C        ==================================================
         DO 90 NUELEM = 1, NBELEM
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
C           TRAITEMENT DES PRESSIONS FIXEES
C           -------------------------------
C           LE CALCUL DU TYPE OBJET DE CHAQUE NOEUD DE L'ELEMENT FINI
C           EN COMMENCANT PAR LE VOLUME, PUIS LES FACES, PUIS LES ARETES
C           PUIS LES SOMMETS CE QUI ASSURE LA PRIORITE DES POINTS SUR
C           LES LIGNES, SURFACES ET VOLUMES
            CALL EFTNND( NOOBVC, NOOBSF, NOOBLA, NOOBPS,
     %                   NUMIOB, MNDOEL,
     %                  'BLPRESSION', MXDOFL, LPBLPR,
     %                   NOTYOB )
C
C           LE RECENSEMENT DES PRESSIONS FIXEES
C           AUX NBNOE NOEUDS DE CET ELEMENT FINI
            DO 80 J=1,NBNOE
C
C              LE NUMERO DU NOEUD DANS LE MAILLAGE DE L'OBJET
               NONOE = MCN( MNNDEL-1 + NUELEM + NBELEM*(J-1) )
C
C              NOMBRE DE DL DU NOEUD
               ND = MCN( MNNDDL + NONOE ) - MCN( MNNDDL + NONOE - 1 )
               IF( ND .NE. NDIM+1 ) GOTO 80
C
C              ICI LE NOEUD SUPPORTE UN DL DE PRESSION
C
C              LE TYPE OBJET DU NOEUD J DE L'EF
               NTYOB = NOTYOB(1,J)
               NOOB  = NOTYOB(2,J)
               MN1   = NOTYOB(3,J)
C
               IF( NTYOB .LE. NDIM .AND. NOOB .GT. 0 ) THEN
C
C                 EXISTE-T-IL UNE PRESSION IMPOSEE EN CE NOEUD ?
C                 ----------------------------------------------
                  IF( MN1 .GT. 0 ) THEN
C                    CALCUL DE LA PRESSION FIXEE EN CE NOEUD
                     N = MNXYZN + WYZNOE + 3 * NONOE - 3
                     XPOIN = RMCN(N)
                     YPOIN = RMCN(N+1)
                     ZPOIN = RMCN(N+2)
                     CALL REBLPR( NTYOB, NOOB, XPOIN, YPOIN, ZPOIN, MN1,
     %                            PRESSION )
C                    LE TEMOIN DE FIXATION DU D.L.
                     MCN(MNNDLX-1+NONOE) = 1
C                    LA VALEUR DE LA PRESSION  FIXEE
                     DMCN(MNVDX+NONOE) = PRESSION
                  ENDIF
               ENDIF
C
 80        CONTINUE
 90      CONTINUE
100   CONTINUE
C
C     NBRDLX NOMBRE DE DL PRESSION FIXES
C     ==================================
      MN = MNNDLX - 1
      DO I=1,NTDLPR
         IF( MCN( MN + I ) .NE. 0 ) THEN
            NBRDLX = NBRDLX + 1
            MCN( MN + NBRDLX ) = I
            DMCN(MNVDX+NBRDLX) = DMCN(MNVDX+I)
            NONOE  = I
            PRESSION = DMCN(MNVDX+I)
         ENDIF
      ENDDO
C
C     REDUCTION EVENTUELLE DES TABLEAUX DE DONNEES
C     ============================================
      IF( NBRDLX .EQ. 0 ) THEN
C        DESTRUCTION DES TABLEAUX
         CALL TNMCDS( 'ENTIER', NTDLPR, MNNDLX )
         CALL TNMCDS( 'REEL2',  NTDLPR, MNVDLX )
         IERR = 1
         GOTO 9999
      ELSE
C        REDUCTION DES 2 TABLEAUX
         CALL TNMCRA('ENTIER', NTDLPR, NBRDLX, MNNDLX )
         CALL TNMCRA('REEL2' , NTDLPR, NBRDLX, MNVDLX )
      ENDIF
C
C     AFFICHAGE FINAL
C     ===============
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10170) NBRDLX, PRESSION, NONOE
      ELSE
         WRITE(IMPRIM,20170) NBRDLX, PRESSION, NONOE
      ENDIF
C
10170 FORMAT(' NOMBRE DE PRESSIONS FIXEES=',I8,
     %       '  DERNIERE PRESSION',G15.6,' FIXEE AU SOMMET',I9)
20170 FORMAT(' NUMBER of FIXED PRESSURES =',I8,
     %       '  LAST FIXED PRESSURE',G15.6,' at VERTEX',I9)
C
 9999 RETURN
      END
