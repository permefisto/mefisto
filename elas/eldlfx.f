      SUBROUTINE ELDLFX( NTDL,   NBCODE, NBTYEL, MNNPEF, NDPGST,
     &                   MNXYZN, NUMIOB, MNDOEL, RELMIN,
     &                   NBRDLX, MONDLX, MNNDLX, MNVDLX,
     &                   IERR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    RECENSER LES NUMEROS ET VALEURS DES DEGRES DE LIBERTE FIXES
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
C NBRDLX : NOMBRE DE DL SUPPORT D'UN DEGRE DE LIBERTE FIXE
C MONDLX : NOMBRE DE MOTS DECLARES DU TABLEAU NO DES DL FIXES
C          =0 SI NBRDLX=0
C MNNDLX : ADRESSE MCN DU TABLEAU MC DES NUMEROS DES DL FIXES
C          =0 SI NBRDLX=0
C MNVDLX : ADRESSE MCN DU TABLEAU MC DES VALEURS DES DL FIXES
C          =0 SI NBRDLX=0
C IERR   : =0 SI PAS D'ERREUR
C          >1 SINON ET TOUTES LES ADRESSES MN DES TMS DLX SONT NULLES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET   ANALYSE NUMERIQUE UPMC PARIS       MARS 1999
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
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
      DOUBLE PRECISION  RELMIN, XPOIN, YPOIN, ZPOIN, FIXAT(3)
C
      IERR   = 0
      NBRDLX = 0
      MNNDLX = 0
      MNVDLX = 0
C
C     LES NUMEROS DES DEGRES DE LIBERTE FIXES
      MONDLX = NTDL
      CALL TNMCDC( 'ENTIER', MONDLX, MNNDLX )
      CALL AZEROI( MONDLX, MCN(MNNDLX) )
C
C     LES TABLEAUX REELS DES VALEURS DES DL FIXES
      CALL TNMCDC( 'REEL2', MONDLX, MNVDLX )
      MNVDX  = ( MNVDLX - 1 ) / 2
      DO 20 I=1,MONDLX
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
            CALL EFTNND( NOOBVC,    NOOBSF, NOOBLA, NOOBPS,
     %                   NUMIOB,    MNDOEL,
     %                  'FIXATION', MXDOEL, LPFIXA,
     %                   NOTYOB )
C
C           LE RECENSEMENT DES FIXATIONS AUX NOEUDS DE CET ELEMENT FINI
            DO 90 J=1,NBNOE
C
C              LE NUMERO DU NOEUD
               NONOE = MCN( MNNDEL-1 + NUELEM + NBELEM*(J-1) )
C              LE TYPE OBJET DU NOEUD
               NTYOB = NOTYOB(1,J)
               NOOB  = NOTYOB(2,J)
C
               IF( NTYOB .LE. NBCODE .AND. NOOB .GT. 0 ) THEN
C
C                 EXISTE-T-IL UNE FIXATION EN CE NOEUD ?
C                 --------------------------------------
                  MN1 = NOTYOB(3,J)
                  IF( MN1 .GT. 0 ) THEN
C                    CALCUL DES DEPLACEMENTS FIXES EN CE NOEUD
                     N = MNXYZN + WYZNOE + 3 * NONOE - 3
                     XPOIN = RMCN(N)
                     YPOIN = RMCN(N+1)
                     ZPOIN = RMCN(N+2)
                     CALL REFIXA( NTYOB,  NOOB, XPOIN,YPOIN,ZPOIN, MN1,
     %                            NBCOFI, FIXAT )
C                    NBCOFI EST LE NOMBRE DE COMPOSANTES FIXEES
C                    LA BOUCLE SUR LES D.L. FIXES
                     NODL = (NONOE - 1) * NBCODE
                     DO 40 K=1,NBCOFI
C                       LE NUMERO D.L. DE LA COMPOSANTE K
                        NODLF = NODL + MCN( MN1 + WUCOFI - 1 + K )
C                       LA VALEUR DU DEPLACEMENT FIXE
                        DMCN(MNVDX+NODLF) = FIXAT(K)
C                       LE TEMOIN DE FIXATION DU D.L.
                        MCN(MNNDLX-1+NODLF) = 1
 40                  CONTINUE
                  ENDIF
               ENDIF
C
 90        CONTINUE
 95      CONTINUE
100   CONTINUE
      IF( IERR .NE. 0 ) GOTO 9999
C
C     COMPRESSION DES TABLEAUX DES DL FIXES ET DES FORCES NODALES
C     ===========================================================
      MN = MNNDLX - 1
      DO 150 I=1,NTDL
         IF( MCN( MN + I ) .NE. 0 ) THEN
            NBRDLX = NBRDLX + 1
            MCN( MN + NBRDLX ) = I
            DMCN(MNVDX+NBRDLX) = DMCN(MNVDX+I)
         ENDIF
 150  CONTINUE
      IF( IERR .GT. 0 ) GOTO 9999
C
C     REDUCTION EVENTUELLE DES TABLEAUX DE DONNEES
C     ============================================
      IF( NBRDLX .EQ. 0 ) THEN
C        DESTRUCTION DES TABLEAUX
         GOTO 9999
      ELSE IF( NBRDLX .LT. NTDL ) THEN
C        REDUCTION DES 2 TABLEAUX
         MONDLX = NBRDLX
         CALL TNMCRA('ENTIER', NTDL, MONDLX, MNNDLX )
         CALL TNMCRA('REEL2' , NTDL, MONDLX, MNVDLX )
      ENDIF
C
C     AFFICHAGE FINAL
C     ===============
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10170) NBRDLX
         WRITE(IMPRIM,10180) (MCN(MN+I),DMCN(MNVDX+I),I=1,NBRDLX)
      ELSE
         WRITE(IMPRIM,20170) NBRDLX
         WRITE(IMPRIM,20180) (MCN(MN+I),DMCN(MNVDX+I),I=1,NBRDLX)
      ENDIF
10170 FORMAT(' NOMBRE de COMPOSANTES de DEPLACEMENTS FIXES =',I5)
20170 FORMAT(' NUMBER of COMPONENTS of FIXED DISPLACEMENTS =',I5)
10180 FORMAT(' DEPLACEMENT',I6,' FIXE a',G15.7,2X,
     %       ' DEPLACEMENT',I6,' FIXE a',G15.7,2X,
     %       ' DEPLACEMENT',I6,' FIXE a',G15.7)
20180 FORMAT(' DISPLACEMENT',I6,' FIXED to',G15.7,2X,
     %       ' DISPLACEMENT',I6,' FIXED to',G15.7,2X,
     %       ' DISPLACEMENT',I6,' FIXED to',G15.7)
      RETURN
C
C     ERREUR RENCONTREE
 9999 CALL TNMCDS( 'ENTIER', MONDLX, MNNDLX )
      CALL TNMCDS( 'REEL2',  MONDLX, MNVDLX )
      MONDLX = 0
      RETURN
      END
