      SUBROUTINE THSNFX( NTDL,   NBTYEL, MNNPEF, NDPGST,
     &                   MNXYZN, NUMIOB, MNDOEL, RELMIN,
     &                   NBFNFX, MONFNX, MNNFNX, MOVFNX, MNVFNX )
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DES SOURCES NODALES DEFINIES AUX POINTS UTILISATEUR
C -----
C         (EN FAIT, CES POINTS NE PEUVENT ETRE QUE DES SOMMETS DES EF
C          CAR LES POINTS UTILISATEUR SONT IDENTIFIES AUX SEULS SOMMETS
C          CF LE SOUS-PROGRAMME DFTOPO )
C
C ENTREES:
C --------
C NTDL   : NOMBRE TOTAL DE DEGRES DE LIBERTE DES TEMPERATURES
C NBTYEL : NOMBRE DE TYPES D'EF DE L'OBJET
C MNNPEF : ADRESSE MCN DES TABLEAUX NPEF"TYPE EF
C NDPGST : CODE DE STOCKAGE DES SOMMETS POINTS NOEUDS
C MNXYZN : ADRESSE MCN DU TMS XYZNOEUD DE L'OBJET
C NUMIOB : NUMERO MINIMAL DES PLSV
C NUMAOB : NUMERO MAXIMAL DES PLSV
C MNDOEL : ADRESSE MCN DES DONNEES THERMIQUES DE L'OBJET
C RELMIN : REEL DOUBLE PRECISION TEMOIN DE NON-INITIALISATION
C
C MODIFIES:
C ---------
C NBFNFX : NOMBRE DE SOURCES NODALES
C MONFNX : NOMBRE DE MOTS DECLARES DU TABLEAU MC NO DL  DES SOURCES NODALES
C MNNFNX : ADRESSE MCN             DU TABLEAU MC NO DL  DES SOURCES NODALES
C MOVFNX : NOMBRE DE MOTS DECLARES DU TABLEAU MC VALEUR DES SOURCES NODALES
C MNVFNX : ADRESSE MCN             DU TABLEAU MC VALEUR DES SOURCES NODALES
C          ATTENTION SI NBFNFX=0 CES 2 ADRESSES SONT NULLES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET  ANALYSE NUMERIQUE UPMC PARIS     FEVRIER 1999
C MODIFS: ALAIN PERRONNET  Texas A & M University at QATAR  Fevrier 2011
C23456---------------------------------------------------------------012
      include"./incl/gsmenu.inc"
      include"./incl/donele.inc"
      include"./incl/a___npef.inc"
      include"./incl/a___xyznoeud.inc"
      include"./incl/ponoel.inc"
      include"./incl/donthe.inc"
      include"./incl/ctemps.inc"
      include"./incl/cthet.inc"
      include"./incl/cnonlin.inc"
C
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1),DMCN(1))
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
C
      INTEGER           NOTYOB(3,MXNOEL)
      INTEGER           NOOBSF(6),NOOBLA(12),NOOBPS(8)
      INTEGER           NUMIOB(4),MNDOEL(4)
C
      DOUBLE PRECISION  RELMIN, XYZP(3), FORCE(3)
C
C     REMISE A JOUR POUR PERMETTRE DES APPELS SUCCESSIFS
C     --------------------------------------------------
      IF( MONFNX .GT. 0 .AND. MNNFNX .GT. 0 ) THEN
          CALL TNMCDS( 'ENTIER', MONFNX, MNNFNX )
      ENDIF
      IF( MOVFNX .GT. 0 .AND. MNVFNX .GT. 0 ) THEN
          CALL TNMCDS( 'REEL2', MOVFNX, MNVFNX )
      ENDIF
C
C     DECLARATIONS
C     ------------
      NBFNFX = 0
      MOREE2 = MOTVAR(6)
      IF( TESTNL .LE. 5 ) THEN
C        THERMIQUE
         NBNODE = NTDL
      ELSE
C        NLSE AVEC UNE PARTIE REELLE ET UNE IMAGINAIRE
         NBNODE = NTDL / 2
      ENDIF
C
C     LES NUMEROS DES DEGRES DE LIBERTE DES SOURCES NODALES
      MONFNX = NTDL
      CALL TNMCDC( 'ENTIER', MONFNX, MNNFNX )
      CALL AZEROI( NTDL, MCN(MNNFNX) )
C
C     LES TABLEAUX REELS DES VALEURS DES SOURCES NODALES
      MOVFNX = NTDL
      CALL TNMCDC( 'REEL2', MOVFNX, MNVFNX )
      MNVFX  = ( MNVFNX - 1 ) / MOREE2
      DO I=1,MOVFNX
         DMCN(MNVFX + I) = RELMIN
      ENDDO
C
C     ========================================
C     LA BOUCLE SUR LES TYPES D'ELEMENTS FINIS
C     ========================================
      DO 100 NOTYEL = 1, NBTYEL
C
C        L'ADRESSE DU TABLEAU NPEF"No TYPE EF
         MNELE = MCN( MNNPEF - 1 + NOTYEL )
C        LE NUMERO DU TYPE DE L'ELEMENT FINI
         NUTYEL = MCN( MNELE + WUTYEL )
C        LE NOMBRE D'EF DE CE TYPE
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
CCC            CALL EFNTND( NOOBVC, NOOBSF, NOOBLA, NOOBPS, NOTYOB )
C           LES LIGNES, SURFACES ET VOLUMES
            CALL EFTNND( NOOBVC, NOOBSF, NOOBLA, NOOBPS,
     %                   NUMIOB, MNDOEL,
     %                  'SOURCE',MXDOTH, LPSOUR,
     %                   NOTYOB )
C
C           LE RECENSEMENT DES SOURCES NODALES AUX POINTS
C           AUX NBNOE NOEUDS DE CET ELEMENT FINI
            DO 90 J=1,NBNOE
C
C              LE NUMERO DU NOEUD DANS LE MAILLAGE DE L'OBJET
               NONOE = MCN( MNNDEL-1 + NUELEM + NBELEM*(J-1) )
C
C              LE TYPE OBJET DU NOEUD J DE L'EF
               NTYOB = NOTYOB(1,J)
               NOOB  = NOTYOB(2,J)
C
               IF( NTYOB .EQ. 1 ) THEN
C
C                 EXISTE-T-IL UNE "SOURCE" NODALE EN CE NOEUD ?
C                 ---------------------------------------------
                  MN1 = NOTYOB(3,J)
                  IF( MN1 .GT. 0 ) THEN
C
C                    OUI: CALCUL DE LA SOURCE NODALE EN CE NOEUD
                     N = MNXYZN + WYZNOE + 3 * NONOE - 3
                     XYZP(1) = RMCN(N)
                     XYZP(2) = RMCN(N+1)
                     XYZP(3) = RMCN(N+2)

                     IF( TESTNL .LE. 5 ) THEN
C                       LA VALEUR DE LA SOURCE NODALE
                        CALL RESOUR( NTYOB, NOOB, 3, XYZP,
     %                               MN1,   DMCN(MNVFX+NONOE) )
C                       LE TEMOIN DE SOURCE NODALE EXISTANTE
                        MCN(MNNFNX-1+NONOE) = 1
                     ELSE
                        FORCE(1) = 0D0
                        FORCE(2) = 0D0
                        FORCE(3) = 0D0
                        CALL REFORC( NTYOB,   NOOB,    3,
     %                               XYZP(1), XYZP(2), XYZP(3),
     %                               0D0,     0D0,     0D0,
     %                               MN1,     FORCE )
                        DO K = 1, 3
                           IF( FORCE(K) .NE. 0D0 ) THEN
C                             NUMERO DU DL SELON K PAR COMPOSANTES
                              NUDL = NONOE + (K-1) * NBNODE
                              DMCN(MNVFX+NUDL) = FORCE(K)
C                             LE TEMOIN DE SOURCE NODALE EXISTANTE
                              MCN(MNNFNX-1+NUDL) = 1
                           ENDIF
                        ENDDO
                     ENDIF
C
                  ENDIF
               ENDIF
 90         CONTINUE
 95      CONTINUE
 100  CONTINUE
C
C     COMPRESSION DES TABLEAUX DES SOURCES NODALES
C     ============================================
      MNN = MNNFNX - 1
      DO I = 1, NTDL
C
         IF( MCN(MNN+I) .NE. 0 ) THEN
C
            IF( DMCN( MNVFX + I ) .NE. 0D0 .AND.
     %          DMCN( MNVFX + I ) .NE. RELMIN ) THEN
C              UNE SOURCE NODALE DE PLUS
               NBFNFX = NBFNFX + 1
               MCN( MNN + NBFNFX ) = I
               DMCN( MNVFX + NBFNFX ) = DMCN( MNVFX + I )
            ENDIF
C
         ENDIF
C
      ENDDO
C
C     REDUCTION EVENTUELLE DES TABLEAUX DES SOURCES NODALES
C     =====================================================
      IF( NBFNFX .EQ. 0 ) THEN
C        DESTRUCTION DES 2 TABLEAUX
         CALL TNMCDS( 'ENTIER', MONFNX, MNNFNX )
         CALL TNMCDS( 'REEL2',  MOVFNX, MNVFNX )
         MONFNX = 0
         MOVFNX = 0
      ELSE
C        REDUCTION DES 2 TABLEAUX
         MONFNX = NBFNFX
         CALL TNMCRA( 'ENTIER', NTDL, MONFNX, MNNFNX )
         MOVFNX = NBFNFX
         CALL TNMCRA( 'REEL2' , NTDL, MOVFNX, MNVFNX )
      ENDIF
C
C     AFFICHAGE FINAL
C     ===============
      WRITE(IMPRIM,10150) NBFNFX
10150 FORMAT(' NOMBRE DE SOURCES NODALES =',I5)
      DO 160 I=1,NBFNFX
         WRITE(IMPRIM,10160) MCN(MNN+I),DMCN(MNVFX+I)
10160    FORMAT(' SOURCE NODALE AU NOEUD =',I6,(T32,4G15.7))
 160  CONTINUE
      RETURN
      END
