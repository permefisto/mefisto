      SUBROUTINE DLFXLI( NCDLF0, NBNOMA, NDIM,   NBTYEL, MNNPEF, NDPGST,
     %                   MNXYZN, NUMIOB, MNDOEL, RELMIN,
     %                   NBRDLX, MONDLX, MNNDLX, MNVDLX, IERR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    RECENSER TOUS LES DEGRES DE LIBERTE FIXES (CONDITION DIRICHLET)
C -----    A PARTIR DU TMS a___fixation ET
C          STOCKER LEUR NUMERO DE DL (RANGES PAR NOEUDS) ET
C          VALEUR POUR UNE ONDE NLSE avec PARTIE REELLE et IMAGINAIRE
C
C ENTREES:
C --------
C NCDLF0 : CODE 0 SI TOUT DL FIXE NON NUL DOIT DECLENCHER IERR NON NUL
C               1 SINON
C NBNOMA : NOMBRE DE NOEUDS DU MAILLAGE DE L'OBJET
C NDIM   : DIMENSION DE L'ESPACE DES COORDONNEES 2 ou 3
C         (SI AXISYMETRIE NDIM=2 X=>R>=0 Y=>Z Z=0)
C NBTYEL : NOMBRE DE TYPES D'EF DE L'OBJET
C MNNPEF : ADRESSE MCN DES TABLEAUX NPEF"TYPE EF
C NDPGST : CODE DE STOCKAGE DES SOMMETS POINTS NOEUDS
C MNXYZN : ADRESSE MCN DU TMS XYZNOEUD DE L'OBJET
C NUMIOB : NUMERO MINIMAL DES PLSV
C MNDOEL : ADRESSE MCN DES DONNEES DE L'OBJET
C RELMIN : REEL DOUBLE PRECISION TEMOIN DE NON-INITIALISATION
C
C SORTIES:
C --------
C NBRDLX : NOMBRE DE DL SUPPORT D'UN DEGRE DE LIBERTE FIXE
C MONDLX : NOMBRE DE MOTS DECLARES DU TABLEAU MC NO DES DL FIXES
C MNNDLX : ADRESSE MCN DU TABLEAU MC DES NUMEROS DES DL FIXES, 0 SINON
C MNVDLX : ADRESSE MCN DU TABLEAU MC DES VALEURS DES DL FIXES, 0 SINON
C IERR   : =0 SI PAS D'ERREUR
C          >1 SINON ET LES ADRESSES MNNDLX ET MNVDLX SONT NULLES
C
C ATTENTION: NUMERO DU DL PAR NOEUD DANS LE VECTEUR GLOBAL(2,NBNOMA)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET TEXAS A & M University at QATAR     MARS 2011
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/gsmenu.inc"
      include"./incl/donele.inc"
      include"./incl/a___npef.inc"
      include"./incl/a___xyznoeud.inc"
      include"./incl/ponoel.inc"
      include"./incl/donthe.inc"
      include"./incl/xyzext.inc"
      include"./incl/cnonlin.inc"
      include"./incl/a___fixation.inc"
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
      INTEGER           NONOEF(20)
C
      DOUBLE PRECISION  RELMIN, XYZP(3), FIXA(3)
C
      IERR   = 0
      MOREE2 = MOTVAR(6)
      NBRDLX = 0
      MONDLX = 0
      MNNDLX = 0
      MNVDLX = 0
C
      IF( TESTNL .EQ. 5 ) THEN
C        NTDL   : NOMBRE TOTAL DE DEGRES DE LIBERTE DE L'ONDE
C                 1 FOIS LE NOMBRE DE NOEUDS DU MAILLAGE
         NTDL = NBNOMA
      ELSE
C        NTDL   : NOMBRE TOTAL DE DEGRES DE LIBERTE DE L'ONDE
C                 2 FOIS LE NOMBRE DE NOEUDS DU MAILLAGE
         NTDL = 2 * NBNOMA
      ENDIF
C
C     CONDITION DE DIRICHLET HOMOGENE SUR LA FRONTIERE
C     LES NUMEROS DES DEGRES DE LIBERTE FIXES
      MONDLX = NTDL
      CALL TNMCDC( 'ENTIER', MONDLX, MNNDLX )
      CALL AZEROI( MONDLX, MCN(MNNDLX) )
      MNDLX = MNNDLX - 1
C
C     LES TABLEAUX REELS DES VALEURS DES DL FIXES
      CALL TNMCDC( 'REEL2', MONDLX, MNVDLX )
      MNVDX = ( MNVDLX - 1 ) / MOREE2
      DO I=1,MONDLX
         DMCN(MNVDX + I) = RELMIN
      ENDDO
C
C     RECHERCHE DES MIN MAX DES COORDONNEES DE L'OBJET
      NBNOEU = MCN(MNXYZN+WNBNOE)
      NBCOOR = MCN(MNXYZN+WBCOON)
      CALL MIMXPT( NBCOOR, NBNOEU, RMCN(MNXYZN+WYZNOE), COOEXT )
C
C     =================================================
C     LA BOUCLE SUR LES TYPES D'ELEMENTS FINIS 2D ou 3D
C     =================================================
      DO 100 NOTYEL = 1, NBTYEL
C
C        L'ADRESSE DU TABLEAU NPEF"No TYPE EF
         MNELE = MCN( MNNPEF - 1 + NOTYEL )
C        LE NUMERO DU TYPE DE L'ELEMENT FINI
         NUTYEL = MCN( MNELE + WUTYEL )
C        LE NOMBRE D'ELEMENTS FINIS DE CE TYPE
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
C           LES NBNOE NOEUDS DE L'ELEMENT FINI
            CALL EFNOEU( MNELE, NUELEM, NBNOE, NONOEF )
C
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
            CALL EFTNND( NOOBVC,    NOOBSF, NOOBLA, NOOBPS,
     %                   NUMIOB,    MNDOEL,
     %                  'FIXATION', MXDOTH, LPCONT,
     %                   NOTYOB )
C
C           LE RECENSEMENT DES FIXATIONS NODAUX
C           AUX NBNOE NOEUDS DE CET ELEMENT FINI
            DO 90 J = 1, NBNOE
C
C              LE NUMERO DU NOEUD J DE L'EF NUELEM DANS LE MAILLAGE DE L'OBJET
               NONOE = NONOEF(J)
C
C              LE TYPE OBJET DU NOEUD J DE L'EF
               NTYOB = NOTYOB(1,J)
               NOOB  = NOTYOB(2,J)
               MN1   = NOTYOB(3,J)
C
               IF( NTYOB .LE. NDIM .AND. NOOB .GT. 0 ) THEN
C
C                 EXISTE-T-IL UNE "FIXATION" EN CE NOEUD ?
C                 ----------------------------------------
                  IF( MN1 .GT. 0 ) THEN
C
C                    CALCUL DE L'ONDE FIXEE EN CE NOEUD
                     N = MNXYZN + WYZNOE + 3 * NONOE - 3
                     XYZP(1) = RMCN(N)
                     XYZP(2) = RMCN(N+1)
                     XYZP(3) = RMCN(N+2)
                     FIXA(1) = 0D0
                     FIXA(2) = 0D0
                     CALL REFIXA( NTYOB,   NOOB,
     %                            XYZP(1), XYZP(2), XYZP(3), MN1,
     %                            NBCOFI,  FIXA )
C
C                    LA VALEUR DE LA COMPOSANTE FIXEE
                     IF( NCDLF0 .EQ. 0 .AND.
     %                 ( ABS(FIXA(1)) .GT. 1D-6 .OR.
     %                   ABS(FIXA(2)) .GT. 1D-6 ) ) THEN
                        NBLGRC(NRERR) = 1
                        IF( LANGAG .EQ. 0 ) THEN
                          KERR(1)='FIXATION NON NULLE EST INTERDITE ICI'
                        ELSE
                           KERR(1)='NOT NULL FIXATION IS FORBIDDEN HERE'
                        ENDIF
                        CALL LEREUR
                        IERR = 1
                        GOTO 9999
                     ENDIF
C
C                    NBCOFI LE NOMBRE DE COMPOSANTES FIXEES
                     DO 85 I = 1, NBCOFI
C
C                       LE NUMERO DE LA COMPOSANTE FIXEE
                        NU = MCN( MN1 + WUCOFI - 1 + I )
C
                        IF( TESTNL .EQ. 5 ) THEN
C
C                          CALCUL DES VALEURS ET VECTEURS PROPRES DE
C                          - Alfa LAPLACIEN u + Beta u**2 u = E u
C                          1 SEULE COMPOSANTE
C                          SI 2 COMPOSANTES DONNEES, LA SECONDE SERA UTILISEE
                           NUDL = NONOE
C
                        ELSE
C
C                          NLSE: PARTIE REELLE ET IMAGINAIRE EXISTENT
C                                ET SON NUMEROTEES PAR NOEUDS (2n-1,2n)
                           IF( NU .EQ. 1 ) THEN
C                             LA PARTIE REELLE EST FIXEE A FIXA(I)
C                             NUMERO DU DL DANS LE VECTEUR GLOBAL
                              NUDL = 2 * NONOE - 1
                           ELSE IF( NU .EQ. 2 ) THEN
C                             LA PARTIE IMAGINAIRE EST FIXEE A FIXA(I)
C                             NUMERO DU DL DANS LE VECTEUR GLOBAL
                              NUDL = 2 * NONOE
                           ELSE
                              GOTO 85
                           ENDIF
C
                        ENDIF
C
C                       LA VALEUR FIXEE DU DL NUDL
                        DMCN(MNVDX+NUDL) = FIXA(I)
C
C                       LE TEMOIN DE VALEUR FIXEE DU DL NUDL
                        MCN(MNDLX+NUDL) = 1
C
 85                  CONTINUE
                  ENDIF
               ENDIF
C
 90        CONTINUE
 95      CONTINUE
100   CONTINUE
C
C     COMPRESSION DES TABLEAUX DES DL FIXES
C     =====================================
      DO 150 I=1,NTDL
C
         IF( MCN( MNDLX + I ) .NE. 0 ) THEN
C           UN DEGRE DE LIBERTE FIXE DE PLUS
            NBRDLX = NBRDLX + 1
C           LE NUMERO GLOBAL DU DL
            MCN( MNDLX + NBRDLX ) = I
C           LA VALEUR FIXEE
            DMCN(MNVDX + NBRDLX ) = DMCN(MNVDX+I)
         ENDIF
C
 150  CONTINUE
C
C     REDUCTION EVENTUELLE DES TABLEAUX DES DEGRES DE LIBERTE FIXES
C     =============================================================
      IF( NBRDLX .EQ. 0 ) THEN
C        DESTRUCTION DES 2 TABLEAUX
         CALL TNMCDS( 'ENTIER', MONDLX, MNNDLX )
         CALL TNMCDS( 'REEL2',  MONDLX, MNVDLX )
         MONDLX = 0
      ELSE
C        REDUCTION DE LA TAILLE DES 2 TABLEAUX
         MONDLX = NBRDLX
         CALL TNMCRA('ENTIER', NTDL, MONDLX, MNNDLX )
         CALL TNMCRA('REEL2' , NTDL, MONDLX, MNVDLX )
      ENDIF
C
C     AFFICHAGE FINAL
C     ===============
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10170) NTDL, NBRDLX, NTDL-NBRDLX
      ELSE
         WRITE(IMPRIM,20170) NTDL, NBRDLX, NTDL-NBRDLX
      ENDIF
10170 FORMAT('DLFXLI: NOMBRE DEGRES de LIBERTE=',I8,'  FIXES =',I8,
     %        '  NON FIXES =',I8)
20170 FORMAT('DLFXLI: NUMBER of DEGREES of FREEDOM=',I8,'  FIXED =',I8,
     %        '  NOT FIXED =',I8)
C
CCC      WRITE(IMPRIM,10180) (MCN(MNDLX+I),DMCN(MNVDX+I),I=1,NBRDLX)
CCC10180 FORMAT(3(' DL',I8,' FIXE A',G15.7,2X))
C
      RETURN
C
C     ERREUR RENCONTREE
C     =================
 9999 CALL TNMCDS( 'ENTIER', MONDLX, MNNDLX )
      CALL TNMCDS( 'REEL2',  MONDLX, MNVDLX )
      RETURN
      END
