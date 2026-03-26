      SUBROUTINE PRESFXST( RELMIN, NTDLPR, NDIM,   MNXYZN, NONOSO,
     %                     NBTYEL, MNNPEF, NUMIOB, MNDOEL,
     %                     NBPRFX, MNNPRFX, MNVPRFX )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    RECENSER LES NUMEROS ET VALEURS DES DEGRES DE LIBERTE FIXES
C -----    EN PRESSION du FLUIDE

C ENTREES:
C --------
C RELMIN : REEL DOUBLE PRECISION TEMOIN DE NON-INITIALISATION
C NTDLPR : NOMBRE TOTAL DE DEGRES DE LIBERTE DES PRESSIONS
C NDIM   : NOMBRE DE COMPOSANTES DE LA VITESSE=DIMENSION DE L'ESPACE
C MNXYZN : ADRESSE MCN DU TMS XYZNOEUD
C NONOSO : TABLEAU NONOSO NO NOEUD => NO SOMMET
C NBTYEL : NOMBRE DE TYPES D'EF DE L'OBJET
C MNNPEF : ADRESSE MCN DES TABLEAUX NPEF"TYPE EF
C NUMIOB : NUMERO MINIMAL DES OBJETS
C NUMAOB : NUMERO MAXIMAL DES OBJETS
C MNDOEL : ADRESSE MCN DES DONNEES DE L'OBJET

C SORTIES:
C --------
C NBPRFX : NOMBRE DE DL PRESSION FIXEES
C          ERREUR SI NBPRFX=0 CAR PAS DE CONDITION AUX LIMITES!
C MNNPRFX: ADRESSE MCN DU TABLEAU MC DES NUMEROS DES DL FIXES
C MNVPRFX: ADRESSE MCN DU TABLEAU MC DES VALEURS DES PRESSIONS FIXEES
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

      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE      (MCN(1), RMCN(1), DMCN(1))

      INTEGER           NONOSO(*), NOTYOB(3,MXNOEL)
      INTEGER           NOOBSF(6), NOOBLA(12), NOOBPS(8)
      INTEGER           NUMIOB(4), MNDOEL(4)

      DOUBLE PRECISION  RELMIN, XPOIN, YPOIN, ZPOIN, PRESSION

C     POUR EVITER DE DOUBLER LES TABLEAUX
      IF( NBPRFX .GT. 0 .AND. MNNPRFX .GT. 0 ) THEN
         CALL TNMCDS( 'ENTIER', NBPRFX, MNNPRFX )
         CALL TNMCDS( 'REEL2',  NBPRFX, MNVPRFX )
      ENDIF

      NBPRFX  = 0
      MNNPRFX = 0
      MNVPRFX = 0

C     LES NUMEROS DES DEGRES DE LIBERTE FIXES
      CALL TNMCDC( 'ENTIER', NTDLPR, MNNPRFX )
      CALL AZEROI( NTDLPR, MCN(MNNPRFX) )

C     LES TABLEAUX REELS DES VALEURS DES DL FIXES
C     VALEUR DE NON INITIALISATION = RELMIN
      CALL TNMCDC( 'REEL2', NTDLPR, MNVPRFX )
      MNVDX  = ( MNVPRFX - 1 ) / 2
      DO I=1,NTDLPR
         DMCN(MNVDX + I) = RELMIN
      ENDDO

C     ========================================
C     LA BOUCLE SUR LES TYPES D'ELEMENTS FINIS
C     ========================================
      DO NOTYEL = 1, NBTYEL

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
         DO NUELEM = 1, NBELEM

C           LES NOEUDS DE L'ELEMENT FINI
C           LE NUMERO DE VOLUME  DE L'EF
C           LE NUMERO DE SURFACE DES FACES   DE L'EF
C           LE NUMERO DE LIGNE   DES ARETES  DE L'EF
C           LE NUMERO DE POINT   DES SOMMETS DE L'EF
            CALL EFPLSV( MNELE , NUELEM,
     %                   NOVCEL, NOSFEL, NOLAEL, NOPSEL,
     %                   NOOBVC, NOOBSF, NOOBLA, NOOBPS )

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

C           LE RECENSEMENT DES PRESSIONS FIXEES
C           AUX NBNSOM SOMMETS DE CET ELEMENT FINI
            DO J=1,NBNSOM

C              LE NUMERO DU NOEUD DANS LE MAILLAGE DE L'OBJET
               NONOE = MCN( MNNDEL-1 + NUELEM + NBELEM*(J-1) )

C              LE NUMERO DU SOMMET
               NOSOE = NONOSO( NONOE )
C
C              ICI LE SOMMET SUPPORTE UN DL DE PRESSION
C              LE TYPE OBJET DU NOEUD J DE L'EF
               NTYOB = NOTYOB(1,J)
               NOOB  = NOTYOB(2,J)
               MN1   = NOTYOB(3,J)

               IF( NTYOB .LE. NDIM .AND. NOOB .GT. 0 ) THEN

C                 EXISTE-T-IL UNE PRESSION IMPOSEE EN CE SOMMET ?
C                 -----------------------------------------------
                  IF( MN1 .GT. 0 ) THEN
C                    CALCUL DE LA PRESSION FIXEE EN CE SOMMET
                     N = MNXYZN + WYZNOE + 3 * NONOE - 3
                     XPOIN = RMCN(N)
                     YPOIN = RMCN(N+1)
                     ZPOIN = RMCN(N+2)
                     CALL REBLPR( NTYOB, NOOB, XPOIN, YPOIN, ZPOIN, MN1,
     %                            PRESSION )
C                    LE TEMOIN DE FIXATION DU D.L. AU SOMMET NOSOE
                     MCN(MNNPRFX-1+NOSOE) = 1
C                    LA VALEUR DE LA PRESSION FIXEE AU SOMMET NOSOE
                     DMCN(MNVDX+NOSOE) = PRESSION
                  ENDIF
               ENDIF

            ENDDO
         ENDDO
      ENDDO

C     NBPRFX NOMBRE DE DL PRESSION FIXES
C     ==================================
      MN = MNNPRFX - 1
      DO I=1,NTDLPR
         IF( MCN( MN + I ) .NE. 0 ) THEN
            NBPRFX = NBPRFX + 1
C           NUMERO DU SOMMET
            MCN( MN + NBPRFX ) = I
C           VALEUR FIXEE AU SOMMET
            DMCN(MNVDX+NBPRFX) = DMCN(MNVDX+I)
C           POUR AFFICHER LE DERNIER
            NOSOE  = I
            PRESSION = DMCN(MNVDX+I)
         ENDIF
      ENDDO

C     REDUCTION EVENTUELLE DES TABLEAUX DE DONNEES
C     ============================================
      IF( NBPRFX .EQ. 0 ) THEN
C        DESTRUCTION DES TABLEAUX
         CALL TNMCDS( 'ENTIER', NTDLPR, MNNPRFX )
         CALL TNMCDS( 'REEL2',  NTDLPR, MNVPRFX )
         IERR = 1
         GOTO 9999
      ELSE
C        REDUCTION DES 2 TABLEAUX
         CALL TNMCRA('ENTIER', NTDLPR, NBPRFX, MNNPRFX )
         CALL TNMCRA('REEL2' , NTDLPR, NBPRFX, MNVPRFX )
      ENDIF

C     AFFICHAGE FINAL
C     ===============
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10170) NBPRFX, PRESSION, NOSOE
      ELSE
         WRITE(IMPRIM,20170) NBPRFX, PRESSION, NOSOE
      ENDIF

10170 FORMAT(' NOMBRE DE PRESSIONS FIXEES=',I8,
     %       '  DERNIERE PRESSION',G15.6,' FIXEE AU SOMMET',I9)
20170 FORMAT(' NUMBER of FIXED PRESSURES =',I8,
     %       '  LAST FIXED PRESSURE',G15.6,' at VERTEX',I9)

 9999 RETURN
      END
