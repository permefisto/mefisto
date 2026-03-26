      SUBROUTINE ELFNFX( NTDL,   NBCODE, NBTYEL, MNNPEF, NDPGST,
     &                   MNXYZN, NUMIOB, MNDOEL, RELMIN,
     &                   NBRFNX, MONFNX, MNNFNX, MNVFNX,
     &                   IERR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    RECENSEMENT DES FORCES NODALES A FIXER ENSUITE
C -----
C
C ENTREES:
C --------
C NTDL   : NOMBRE TOTAL DE DEGRES DE LIBERTE DES DEPLACEMENTS
C NBCODE : NOMBRE DE COMPOSANTES DU DEPLACEMENT EN CHAQUE NOEUD
C NBTYEL : NOMBRE DE TYPES D'EF DE L'OBJET
C MNNPEF : ADRESSE MCN DES TABLEAUX NPEF"TYPE EF
C NDPGST : CODE DE STOCKAGE DES SOMMETS POINTS NOEUDS
C MNXYZN : ADRESSE MCN DU TMS XYZNOEUD DE L'OBJET
C NUMIOB : NUMERO MINIMAL DES OBJETS
C MNDOEL : ADRESSE MCN DES DONNEES DE L'OBJET
C RELMIN : REEL DOUBLE PRECISION TEMOIN DE NON-INITIALISATION
C
C SORTIES:
C --------
C NBRFNX : NOMBRE DE COMPOSANTES DES FORCES NODALES
C MONFNX : NOMBRE DE MOTS DECLARES DU TABLEAU NO DL  DES FORCES NODALES
C          0 SI NBRFNX=0
C MNNFNX : ADRESSE MCN DU TABLEAU MC DES NUMEROS DES FORCES NODALES
C          0 SI NBRFNX=0
C MNVFNX : ADRESSE MCN DU TABLEAU MC DES VALEURS DES FORCES NODALES
C          0 SI NBRFNX=0
C IERR   : =0 SI PAS D'ERREUR
C          >1 SINON ET TOUTES LES ADRESSES MN DES FNX SONT NULLES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET   ANALYSE NUMERIQUE UPMC PARIS       MARS 1999
C23456---------------------------------------------------------------012
      include"./incl/donela.inc"
      include"./incl/donele.inc"
      include"./incl/a___npef.inc"
      include"./incl/a___xyznoeud.inc"
      include"./incl/a___fixation.inc"
      include"./incl/a___force.inc"
      include"./incl/ponoel.inc"
      include"./incl/gsmenu.inc"
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
      DOUBLE PRECISION  RELMIN, D, XPOIN, YPOIN, ZPOIN, FORCE(3)
C
      IERR   = 0
      NBRFNX = 0
      MNNFNX = 0
      MNVFNX = 0
C
C     LES NUMEROS DES DEGRES DE LIBERTE DES FORCES NODALES
      MONFNX = NTDL
      CALL TNMCDC( 'ENTIER', MONFNX, MNNFNX )
      CALL AZEROI( NTDL, MCN(MNNFNX) )
C
C     LES TABLEAUX REELS DES VALEURS DES FORCES NODALES
      CALL TNMCDC( 'REEL2', MONFNX, MNVFNX )
      MNVFX  = ( MNVFNX - 1 ) / 2
      DO 20 I=1,MONFNX
         DMCN(MNVFX + I) = RELMIN
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
C
C        L'ADRESSE MCN DES NOEUDS ET POINTS GEOMETRIQUES DES ELEMENTS
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
            CALL EFTNND( NOOBVC, NOOBSF, NOOBLA, NOOBPS,
     %                   NUMIOB, MNDOEL,
     %                  'FORCE', MXDOEL, LPFORC,
     %                   NOTYOB )
C
C           LE RECENSEMENT DES FORCES NODALES AUX NOEUDS DE CET ELEMENT FINI
C           ================================================================
            DO 90 J=1,NBNOE
C              LE NUMERO DU NOEUD
               NONOE = MCN( MNNDEL-1 + NUELEM + NBELEM*(J-1) )
C              LE TYPE OBJET DU NOEUD
               NTYOB = NOTYOB(1,J)
               NOOB  = NOTYOB(2,J)
C
               IF( NTYOB .EQ. 1 ) THEN
C
C                 UNE FORCE PONCTUELLE EST ELLE FIXEE ?
C                 -------------------------------------
                  MN1 = NOTYOB(3,J)
                  IF( MN1 .GT. 0 ) THEN
C
C                    CALCUL DES FORCES PONCTUELLES EN CE NOEUD
                     N     = MNXYZN + WYZNOE + 3 * NONOE - 3
                     XPOIN = RMCN(N)
                     YPOIN = RMCN(N+1)
                     ZPOIN = RMCN(N+2)
C                    PAS DE NORMALE EN UN POINT
                     CALL REFORC( NTYOB, NOOB,  NBCODE,
     %                            XPOIN, YPOIN, ZPOIN,
     %                            0D0,   0D0,   0D0,
     %                            MN1,   FORCE )
C                    ATTENTION: FORCES RAMENE NBCODE COMPOSANTES!
C
C                    LA BOUCLE SUR LES COMPOSANTES DE LA FORCE NODALE
                     NODL = (NONOE - 1) * NBCODE
                     DO 60 K=1,NBCODE
C                       PARCOURS DES NBCODE COMPOSANTES
C                       LA VALEUR DE LA FORCE PONCTUELLE
                        D = FORCE(K)
                        IF( D .NE. 0D0 ) THEN
C                          LA VALEUR DE LA FORCE PONCTUELLE
                           DMCN(MNVFX+NODL+K) = D
C                          LE TEMOIN SUR LE DL ASSOCIE A LA FORCE
                           MCN(MNNFNX-1+NODL+K) = 1
                        ENDIF
 60                  CONTINUE
C
                  ENDIF
               ENDIF
 90        CONTINUE
 95      CONTINUE
100   CONTINUE
      IF( IERR .NE. 0 ) GOTO 9999
C
C     COMPRESSION DES TABLEAUX DES FORCES NODALES
C     ===========================================
      MNN = MNNFNX - 1
      DO 150 I=1,NTDL
         IF( MCN(MNN+I) .NE. 0 ) THEN
            NBRFNX = NBRFNX + 1
            MCN( MNN + NBRFNX ) = I
            KK  = MNVFX + I - 1
            LL  = MNVFX + NBRFNX-1
            DMCN( LL + 1 ) = DMCN( KK + 1 )
         ENDIF
 150  CONTINUE
      IF( IERR .GT. 0 ) GOTO 9999
C
C     REDUCTION EVENTUELLE DES TABLEAUX DE DONNEES
C     ============================================
      IF( NBRFNX .EQ. 0 ) THEN
C        DESTRUCTION DES TABLEAUX
         GOTO 9999
      ELSE IF( NBRFNX .LT. NTDL ) THEN
C        REDUCTION DES 2 TABLEAUX
         MONFNX = NBRFNX
         CALL TNMCRA('ENTIER', NTDL, MONFNX, MNNFNX )
         CALL TNMCRA('REEL2',  NTDL, MONFNX, MNVFNX )
      ENDIF
C
C     AFFICHAGE FINAL
C     ===============
      WRITE(IMPRIM,10150) NBRFNX
10150 FORMAT(' NOMBRE DE COMPOSANTES DE FORCES NODALES =',I5)
      DO 160 I=1,NBRFNX
         WRITE(IMPRIM,10160) MCN(MNN+I),DMCN(MNVFX+I-1+1)
10160    FORMAT(' FORCE NODALE ',I6,' =',G15.7)
 160  CONTINUE
      RETURN
C
C     ERREUR RENCONTREE
 9999 CALL TNMCDS( 'ENTIER', MONFNX, MNNFNX )
      CALL TNMCDS( 'REEL2',  MONFNX, MNVFNX )
      MONFNX = 0
      RETURN
      END
