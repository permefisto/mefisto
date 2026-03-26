      SUBROUTINE VIPRFX( RELMIN, NTDL,   NDIM,   MNXYZN, NDDLNO,
     %                   NBTYEL, MNNPEF, NUMIOB, MNDOEL,
     %                   NBVCFX, NBDLFX, MNNDLFX, MNVDLFX )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    RECENSER LES NUMEROS ET VALEURS DES DEGRES DE LIBERTE FIXES
C -----    EN VITESSE et/ou PRESSION du FLUIDE
C          LA VITESSE AU BARYCENTRE DES EF N'EST PAS PRISE EN COMPTE
C          POUR LES EF DE BREZZI-FORTIN CAR LA VALEUR AU BARYCENTRE
C          N'EST PAS SUR LA FRONTIERE

C ENTREES:
C --------
C RELMIN : REEL DOUBLE PRECISION TEMOIN DE NON-INITIALISATION
C NTDL   : NOMBRE TOTAL DE DEGRES DE LIBERTE DES VITESSES-PRESSION
C NDIM   : NOMBRE DE COMPOSANTES DE LA VITESSE=DIMENSION DE L'ESPACE
C MNXYZN : ADRESSE MCN DU TMS XYZNOEUD
C NDDLNO : TABLEAU DES POINTEURS SUR LE DERNIER DL DE CHAQUE NOEUD
C NBTYEL : NOMBRE DE TYPES D'EF DE L'OBJET
C MNNPEF : ADRESSE MCN DES TABLEAUX NPEF"TYPE EF
C NUMIOB : NUMERO MINIMAL DES OBJETS
C MNDOEL : ADRESSE MCN DES DONNEES DE L'OBJET

C SORTIES:
C --------
C NBVCFX  : NOMBRE DE DEGRES DE LIBERTE DE VITESSE CONVECTEE IMPOSEE
C NBDLFX  : NOMBRE DE DL VITESSE-PRESSION FIXES
C           ERREUR SI NBDLFX=0 CAR PAS DE CONDITION AUX LIMITES!
C MNNDLFX : ADRESSE MCN DU TABLEAU MC DES NUMEROS DES DL FIXES
C MNVDLFX : ADRESSE MCN DU TABLEAU MC DES VALEURS DES DL FIXES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Laboratoire J-L.LIONS UPMC Paris     Mai 2007
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/donflu.inc"
      include"./incl/donele.inc"
      include"./incl/a___npef.inc"
      include"./incl/a___xyznoeud.inc"
      include"./incl/a___blvitesse.inc"
      include"./incl/a___force.inc"
      include"./incl/ponoel.inc"
      include"./incl/gsmenu.inc"

      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE      (MCN(1), RMCN(1), DMCN(1))

      INTEGER           NDDLNO(0:*), NOTYOB(3,MXNOEL)
      INTEGER           NOOBSF(6), NOOBLA(12), NOOBPS(8)
      INTEGER           NUMIOB(4), MNDOEL(4)

      DOUBLE PRECISION  D, RELMIN, XPOIN,YPOIN,ZPOIN, BLVIT(3), PRESSION

C     POUR EVITER DE DOUBLER LES TABLEAUX
      IF( NBDLFX .GT. 0 .AND. MNNDLFX .GT. 0 ) THEN
         CALL TNMCDS('ENTIER', NBDLFX, MNNDLFX)
         CALL TNMCDS('REEL2',  NBDLFX, MNVDLFX)
      ENDIF

      NBVCFX = 0
      NBDLFX = 0
      MNNDLFX = 0
      MNVDLFX = 0

C     LES NUMEROS DES DEGRES DE LIBERTE FIXES
      CALL TNMCDC( 'ENTIER', NTDL, MNNDLFX )
      CALL AZEROI( NTDL, MCN(MNNDLFX) )

C     LES TABLEAUX REELS DES VALEURS DES DL FIXES
C     VALEUR DE NON INITIALISATION = RELMIN
      CALL TNMCDC( 'REEL2', NTDL, MNVDLFX )
      MOREE2 = MOTVAR(6)
      MNVDX  = ( MNVDLFX - 1 ) / MOREE2
      DO I=1,NTDL
         DMCN(MNVDX + I) = RELMIN
      ENDDO

C     ========================================
C     LA BOUCLE SUR LES TYPES D'ELEMENTS FINIS
C     ========================================
      DO 100 NOTYEL = 1, NBTYEL

C        L'ADRESSE DU TABLEAU ELEMENTS (NPEF)
         MNELE = MCN( MNNPEF - 1 + NOTYEL )
C        LE NUMERO DU TYPE DE L'ELEMENT FINI
         NUTYEL = MCN( MNELE + WUTYEL )
C        LE NOMBRE D'ELEMENTS DE CE TYPE
         NBELEM = MCN( MNELE + WBELEM )

C        L'ADRESSE MCN DES NOEUDS ET POINTS GEOMETRIQUES DES ELEMENTS
         MNNDEL = MNELE + WUNDEL

C        LES CARACTERISTIQUES DE L'ELEMENT FINI
C        ON TROUVE: NBPOE, NBNOE, NARET
         CALL ELTYCA( NUTYEL )

C        LA BOUCLE SUR LES ELEMENTS FINIS DE CE TYPE NUTYEL
C        ==================================================
         DO 90 NUELEM = 1, NBELEM

C           LES NOEUDS DE L'ELEMENT FINI NUELEM
C           LE NUMERO DE VOLUME  DE L'EF
C           LE NUMERO DE SURFACE DES FACES   DE L'EF
C           LE NUMERO DE LIGNE   DES ARETES  DE L'EF
C           LE NUMERO DE POINT   DES SOMMETS DE L'EF
            CALL EFPLSV( MNELE , NUELEM,
     %                   NOVCEL, NOSFEL, NOLAEL, NOPSEL,
     %                   NOOBVC, NOOBSF, NOOBLA, NOOBPS )

C           TRAITEMENT DES COMPOSANTES FIXEES DE LA VITESSE
C           -----------------------------------------------
C           LE CALCUL DU TYPE OBJET DE CHAQUE NOEUD DE L'ELEMENT FINI
C           EN COMMENCANT PAR LE VOLUME, PUIS LES FACES, PUIS LES ARETES
C           PUIS LES SOMMETS CE QUI ASSURE PAR ECRASEMENT LA PRIORITE DES
C           SURFACES SUR LE VOLUME
C           DES LIGNES SUR LES SURFACES ET LES VOLUMES
C           DES POINTS SUR LES LIGNES, LES SURFACES ET LES VOLUMES
            CALL EFTNND( NOOBVC,    NOOBSF, NOOBLA, NOOBPS,
     %                   NUMIOB,    MNDOEL,
     %                  'BLVITESSE', MXDOFL, LPBLVI,
     %                   NOTYOB )

C           LE RECENSEMENT DES VITESSES BLOQUEES AUX NOEUDS DE CET EF
            DO 40 J=1,NBNOE

C              LE NUMERO DU NOEUD
               NONOE = MCN( MNNDEL-1 + NUELEM + NBELEM*(J-1) )
C              LE TYPE OBJET DU NOEUD
               NTYOB = NOTYOB(1,J)
               NOOB  = NOTYOB(2,J)

               IF( NTYOB .LE. NDIM .AND. NOOB .GT. 0 ) THEN

C                 EXISTE-T-IL UNE FIXATION VITESSE EN CE NOEUD ?
C                 ----------------------------------------------
                  MN1 = NOTYOB(3,J)
                  IF( MN1 .GT. 0 ) THEN
C                    CALCUL DES VITESSES FIXES EN CE NOEUD
                     N = MNXYZN + WYZNOE + 3 * NONOE - 3
                     XPOIN = RMCN(N)
                     YPOIN = RMCN(N+1)
                     ZPOIN = RMCN(N+2)
                     CALL REBLVI( NTYOB,  NOOB, XPOIN,YPOIN,ZPOIN, MN1,
     %                            NBCOBV, BLVIT )
C                    NBCOBV EST LE NOMBRE DE COMPOSANTES FIXEES
                     DO 30 K=1,NBCOBV
C                       LE NUMERO DE LA COMPOSANTE K
                        NOCOMP = MCN( MN1 + WUCOBV - 1 + K )
C                       LE NO DU DL EN VITESSE-PRESSION
                        NODLVP = NDDLNO( NONOE - 1 ) + NOCOMP
C                       LE TEMOIN DE FIXATION DU D.L.
                        MCN(MNNDLFX-1+NODLVP) = NOCOMP
C                       LA VALEUR DE LA VITESSE FIXEE
                        DMCN(MNVDX+NODLVP) = BLVIT(K)
 30                  CONTINUE
                  ENDIF
               ENDIF

 40         CONTINUE

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

C           LE RECENSEMENT DES CONTACTS NODAUX
C           AUX NBNOE NOEUDS DE CET ELEMENT FINI
            DO 80 J=1,NBNOE

C              LE NUMERO DU NOEUD DANS LE MAILLAGE DE L'OBJET
               NONOE = MCN( MNNDEL-1 + NUELEM + NBELEM*(J-1) )
C
C              NOMBRE DE DL DU NOEUD
               ND = NDDLNO( NONOE ) - NDDLNO( NONOE - 1 )
               IF( ND .NE. NDIM+1 ) GOTO 80

C              ICI LE NOEUD SUPPORTE UN DL DE PRESSION

C              LE TYPE OBJET DU NOEUD J DE L'EF
               NTYOB = NOTYOB(1,J)
               NOOB  = NOTYOB(2,J)
               MN1   = NOTYOB(3,J)

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
C                    LA VALEUR DE LA PRESSION FIXEE
C                    LE NO DU DL EN VITESSE-PRESSION
                     NODLVP = NDDLNO( NONOE - 1 ) + NDIM+1
C                    LE TEMOIN DE FIXATION DU D.L.
                     MCN(MNNDLFX-1+NODLVP) = NDIM+1
C                    LA VALEUR DE LA PRESSION  FIXEE
                     DMCN(MNVDX+NODLVP) = PRESSION
                  ENDIF
               ENDIF

 80        CONTINUE
 90      CONTINUE
100   CONTINUE

C     COMPRESSION DU TABLEAU DES DL VITESSE-PRESSION FIXES
C     ====================================================
      NBDLFX = 0
      NBVCFX = 0
      MN = MNNDLFX - 1
      DO I=1,NTDL
         IF( MCN( MN + I ) .NE. 0 ) THEN
            NBDLFX = NBDLFX + 1
            MCN( MN + NBDLFX ) = I
            D = DMCN(MNVDX+I)
            IF( D .EQ. 1D222 ) NBVCFX = NBVCFX + 1
            DMCN(MNVDX+NBDLFX) = D
         ENDIF
      ENDDO

C     REDUCTION EVENTUELLE DES TABLEAUX DE DONNEES
C     ============================================
      IF( NBDLFX .EQ. 0 ) THEN
C        DESTRUCTION DES TABLEAUX SI AUCUN DL FIXE
         CALL TNMCDS( 'ENTIER', NTDL, MNNDLFX )
         CALL TNMCDS( 'REEL2',  NTDL, MNVDLFX )
         GOTO 9999
      ELSE
C        REDUCTION DES 2 TABLEAUX de NTDL a NBDLFX
         CALL TNMCRA('ENTIER', NTDL, NBDLFX, MNNDLFX )
         CALL TNMCRA('REEL2' , NTDL, NBDLFX, MNVDLFX )
      ENDIF

C     AFFICHAGE FINAL
C     ===============
      NBV   = MIN(10,NBDLFX)
      MN    = MNNDLFX - 1
      MNVDX = ( MNVDLFX - 1 ) / MOREE2
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10170) NBDLFX, NBVCFX, NBV
         WRITE(IMPRIM,10180) (MCN(MN+I),DMCN(MNVDX+I),
     %                        I=NBDLFX-NBV+1,NBDLFX)
      ELSE
         WRITE(IMPRIM,20170) NBDLFX, NBVCFX, NBV
         WRITE(IMPRIM,20180) (MCN(MN+I),DMCN(MNVDX+I),
     %                        I=NBDLFX-NBV+1,NBDLFX)
      ENDIF

10170 FORMAT(' NOMBRE DE COMPOSANTES VITESSES-PRESSIONS  FIXEES=',I9/
     %       ' NOMBRE DE COMPOSANTES VITESSES CONVECTEES FIXEES=',I9/
     %       ' NOMBRE DERNIERES COMPOSANTES AFFICHEES          =',I9)
20170 FORMAT(' NUMBER of FIXED VELOCITY-PRESSURE COMPONENTS =',I9/
     %       ' NUMBER of FIXED CONVECTED VELOCITY COMPONENTS=',I9/
     %       ' NUMBER of LAST PRINTED COMPONENTS            =',I9)

10180 FORMAT(2(' VITESSE ou PRESSION',I9,  ' FIXEE A', G15.7,' | '))
20180 FORMAT(2(' VELOCITY or PRESSURE DoF',I9,' FIXED to',G15.7,' | '))

 9999 RETURN
      END
