      SUBROUTINE SUEXQ3( NBSOCT , NONSTR , LINSTR , LISTR  ,
     %                   MNSOCT , NUCOTE , NSENS  ,
     %                   NBSOBO , COSOBO , IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    GENERER LE MAILLAGE D'UN QUADRANGLE DEFINI PAR DES
C -----    LIGNES STRUCTUREES ET UNE OU DEUX LIGNES FERMEES
C
C ENTREES:
C --------
C NBSOCT : NOMBRE DE NOEUDS DE CHAQUE LIGNE
C NONSTR : NOMBRE DES LIGNES NON STRUCTUREES
C LINSTR : NUMERO DES LIGNES NON STRUCTUREES
C LISTR  : NUMERO DES LIGNES STRUCTUREES
C MNSOCT : ADRESSE MCN DES COORDONNES DE CHAQUE LIGNE
C NBSOBO : DIMENSION DU TABLEAU COSOBO
C
C SORTIES:
C --------
C NUCOTE : ADRESSE DES 4 LIGNES ORDONNEES
C NSENS  : SENS DE LECTURE DES CES QUATRE LIGNES
C IERR   : 0 SI PAS D'ERREUR
C        > 0 SINON
C
C TRAVAIL:
C --------
C COSOBO : TABLEAU DES COORDONNES DU BORD (SERT AUSSI DE TABLEAU COPIE)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : CHRISTOPHE DOURSAT  ANALYSE NUMERIQUE PARIS   JANVIER   1991
C.......................................................................
      IMPLICIT INTEGER (W)
      include"./incl/gsmenu.inc"
      COMMON / EPSSSS / EPZERO,EPSXYZ
C
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE       (MCN(1),RMCN(1))
C
      INTEGER           NBSOCT(4),MNSOCT(4),NUCOTE(4),NSENS(4)
      INTEGER           LISTR(4),LINSTR(4),LVOIS(2,2,2)
      DIMENSION         COSOBO(3*NBSOBO)
      REAL              XYZI(3,4),XYZF(3,4)
      INTEGER           SENSFR
      DOUBLE PRECISION  D2D3(3,3)
      REAL              POINTO(3),POINTI(3),POINTJ(3),XYZ(3),XY3(3)
C
C     TENTATIVE DE RANGEMENT DES COTES SUIVANT L'ORDRE SUIVANT :
C     ==========================================================
C
C        333333333333333333         3333333333333333333|
C       3                  3       3                   |
C      3      11111         3     3      11111         2
C      3     1     1---2----3     3     1     111111111|
C      3     1     1+++4++++3     3     1     111111111+
C      3      11111         3     3      11111         4
C       3                  3       3                   +
C        333333333333333333         3333333333333333333+
C
C
C        2 LIGNES FERMEES           1 LIGNE FERMEE
C          MAILLAGE EN O              MAILLAGE EN C
C
      IERR = 0
C
C     MISE EN MEMOIRE DES COORDONNES DES EXTREMITES DES LIGNES STRUCT.
C     ----------------------------------------------------------------
      DO 10 LL=1,4-NONSTR
        LS = LISTR(LL)
C
C       LES COORDONNEES DES POINTS INITIAUX
        IAI = MNSOCT(LS)
        XYZI(1,LS) = RMCN(IAI  )
        XYZI(2,LS) = RMCN(IAI+1)
        XYZI(3,LS) = RMCN(IAI+2)
C
C       LES COORDONNEES DES POINTS FINAUX
        IAF = MNSOCT(LS) + 3 * (NBSOCT(LS)-1)
        XYZF(1,LS) = RMCN(IAF  )
        XYZF(2,LS) = RMCN(IAF+1)
        XYZF(3,LS) = RMCN(IAF+2)
   10 CONTINUE
C
C     RECHERCHE DE CES POINTS (XYZI ET XYZF) DANS LES LIGNES FERMEES
C     --------------------------------------------------------------
      DO 140 LL=1,NONSTR
C       REVUE DES LIGNES FERMEES
C
        LNS = LINSTR(LL)
        INDIC = 0
        DO 70 N=1,NBSOCT(LNS)-1
C         REVUE DES POINTS DE CHAQUE LIGNE
C
          NIDENT = 0
C         COMPARAISON AVEC LES EXTREMITES DEBUT (XYZI)
          DO 60 I=1,4-NONSTR
            LS = LISTR(I)
            DO 20 J=1,3
              R1 = RMCN(MNSOCT(LNS)+(N-1)*3+(J-1))
              R2 = XYZI  (J,LS)
              IF     ( ABS(R1) .LE. EPZERO ) THEN
                IF  ( ABS(R2) .GT. EPZERO ) GOTO 30
              ELSE IF( ABS(R2) .LE. EPZERO ) THEN
                GOTO 30
              ELSE IF( ABS(R1-R2) .GT. ABS(R1) * EPSXYZ ) THEN
                 GOTO 30
              ENDIF
   20       CONTINUE
C           IDENTIFIACTION
            NIDENT = NIDENT+1
            IF (NIDENT.LE.2) THEN
              LVOIS(LL,NIDENT,1) = LS
              LVOIS(LL,NIDENT,2) = 1
            ELSE
              IERR = 1
            ENDIF
            IF (INDIC.EQ.0) THEN
              INDIC = N
            ELSE
              IF (INDIC.NE.N) THEN
                IERR = 1
              ENDIF
            ENDIF
   30       CONTINUE
C           COMPARAISON AVEC LES EXTREMITES FIN (XYZF)
            DO 40 J=1,3
              R1 = RMCN(MNSOCT(LNS)+(N-1)*3+(J-1))
              R2 = XYZF (J,LS)
              IF     ( ABS(R1) .LE. EPZERO ) THEN
                IF  ( ABS(R2) .GT. EPZERO ) GOTO 50
              ELSE IF( ABS(R2) .LE. EPZERO ) THEN
                GOTO 50
              ELSE IF( ABS(R1-R2) .GT. ABS(R1) * EPSXYZ ) THEN
                 GOTO 50
              ENDIF
   40       CONTINUE
C           IDENTIFIACTION
            NIDENT = NIDENT+1
            IF (NIDENT.LE.2) THEN
              LVOIS(LL,NIDENT,1) = LS
              LVOIS(LL,NIDENT,2) = 2
            ELSE
              IERR = 1
            ENDIF
            IF (INDIC.EQ.0) THEN
              INDIC = N
            ELSE
              IF (INDIC.NE.N) THEN
                IERR = 1
              ENDIF
            ENDIF
   50       CONTINUE
   60     CONTINUE
   70   CONTINUE
C
C       SORTIE SI IL N'Y A PAS EU IDENTIFICATION
C
        IF (INDIC.EQ.0) THEN
          IERR = 1
        ENDIF
        IF (IERR.NE.0) THEN
          RETURN
        ENDIF
C
C       RECONSTRUCTION DE LA SUITE DES COORDONNEES DES LIGNES NON
C       STRUCTUREES EN PARTANT DES POINTS IDENTIFIES DEBUT OU FIN
C
        DO 90 I=1,NBSOCT(LNS)
          DO 80 J=1,3
            COSOBO(3*(I-1)+J) = RMCN(MNSOCT(LNS)+(I-1)*3+J-1)
   80     CONTINUE
   90   CONTINUE
        DO 110 I=INDIC,NBSOCT(LNS)
          DO 100 J=1,3
            RMCN(MNSOCT(LNS)+(I-INDIC)*3+J-1) = COSOBO(3*(I-1)+J)
  100     CONTINUE
  110   CONTINUE
        DO 130 I=2,INDIC
          DO 120 J=1,3
            RMCN(MNSOCT(LNS)+(NBSOCT(LNS)-INDIC+I-1)*3+J-1) =
     %              COSOBO(3*(I-1)+J)
  120     CONTINUE
  130   CONTINUE
  140 CONTINUE
C
C     PREPARATION DU RECOLLEMENT DES LIGNES
C
      NUCOTE(1) = LVOIS(1,1,1)
      IF (LVOIS(1,1,2).EQ.1) THEN
        NSENS(NUCOTE(1)) = -1
      ELSE
        NSENS(NUCOTE(1)) = 1
      ENDIF
      NUCOTE(2) = LINSTR(1)
      NSENS(NUCOTE(2))  = 1
      NUCOTE(3) = LVOIS(1,2,1)
      IF (LVOIS(1,2,2).EQ.1) THEN
        NSENS(NUCOTE(3)) =  1
      ELSE
        NSENS(NUCOTE(3)) = -1
      ENDIF
      IF (NONSTR.EQ.2) THEN
        NUCOTE(4) = LINSTR(2)
        NSENS (NUCOTE(4)) =  1
      ELSE
        LS = 0
        DO 150 I=1,4-NONSTR
          IF ((LISTR(I).NE.NUCOTE(1)).AND.(LISTR(I).NE.NUCOTE(3))) THEN
            LS = LISTR(I)
            NUCOTE(4) = LISTR(I)
          ENDIF
  150   CONTINUE
C
C       RECHERCHE DU RECOLLEMENT DE LA DERNIERE LIGNE
C       AUX DEUX AUTRES LIGNES STRUCTURES 1 ET 3
C
        L1        = NUCOTE(1)
        NSENS(LS) = 0
        DO 160 J=1,3
          R1 =  XYZI  (J,LS)
          IF (NSENS(L1).EQ.1) THEN
            R2 =  XYZI  (J,L1)
          ELSE
            R2 =  XYZF  (J,L1)
          ENDIF
          IF     ( ABS(R1) .LE. EPZERO ) THEN
            IF  ( ABS(R2) .GT. EPZERO ) GOTO 170
          ELSE IF( ABS(R2) .LE. EPZERO ) THEN
            GOTO 170
          ELSE IF( ABS(R1-R2) .GT. ABS(R1) * EPSXYZ ) THEN
            GOTO 170
          ENDIF
  160   CONTINUE
        NSENS(LS) = -1
        GOTO 190
  170   CONTINUE
        DO 180 J=1,3
          R1 =  XYZF  (J,LS)
          IF (NSENS(L1).EQ.1) THEN
            R2 =  XYZI  (J,L1)
          ELSE
            R2 =  XYZF  (J,L1)
          ENDIF
          IF     ( ABS(R1) .LE. EPZERO ) THEN
            IF  ( ABS(R2) .GT. EPZERO ) GOTO 190
          ELSE IF( ABS(R2) .LE. EPZERO ) THEN
            GOTO 190
          ELSE IF( ABS(R1-R2) .GT. ABS(R1) * EPSXYZ ) THEN
            GOTO 190
          ENDIF
  180   CONTINUE
        NSENS(LS) =  1
  190   CONTINUE
        IF (NSENS(LS).EQ.0) THEN
          NBLGRC(NRERR) = 1
          KERR(1) =  " LE BORD N'EST PAS FERME"
          CALL LEREUR
          IERR = 2
          RETURN
        ELSE
C
C         VERIFICATION A L'AUTRE EXTREMITE
          L1 = NUCOTE(3)
          DO 200 J=1,3
            IF (NSENS(LS).EQ.1) THEN
              R1 =  XYZI  (J,LS)
            ELSE
              R1 =  XYZF  (J,LS)
            ENDIF
            IF (NSENS(L1).EQ.1) THEN
              R2 =  XYZF  (J,L1)
            ELSE
              R2 =  XYZI  (J,L1)
            ENDIF
            IF     ( ABS(R1) .LE. EPZERO ) THEN
              IF  ( ABS(R2) .GT. EPZERO ) IERR = 2
            ELSE IF( ABS(R2) .LE. EPZERO ) THEN
              IERR = 2
            ELSE IF( ABS(R1-R2) .GT. ABS(R1) * EPSXYZ ) THEN
              IERR = 2
            ENDIF
  200     CONTINUE
          IF (IERR.NE.0) THEN
            NBLGRC(NRERR) = 1
            KERR(1) =  " LE BORD N'EST PAS FERME"
            CALL LEREUR
            RETURN
          ENDIF
        ENDIF
      ENDIF
C
C     REMPLISSAGE DU TABLEAU COSOBO DES COORDONNES ET DE LA FRONTIERE
C
      NTESTI = 0
  210 CONTINUE
C
C    * LE COTE 1
      N=NUCOTE(1)
      N1 = NBSOCT(N)
      DO 230 NP=1,N1
         IF(NSENS(N).GE.0) THEN
            IA  = MNSOCT(N) - 1 + 3 * (NP-1)
         ELSE
            IA  = MNSOCT(N) - 1 + 3 * (NBSOCT(N)-NP)
         ENDIF
         IN  = 3 * (NP-1)
         DO 220 I=1,3
            COSOBO(IN+I) = RMCN(IA+I)
  220   CONTINUE
  230 CONTINUE
C
C    * LE COTE 2
      N=NUCOTE(2)
      N2 = NBSOCT(N)-2
      DO 250 NP=2,NBSOCT(N)-1
         IF(NSENS(N).GE.0) THEN
            IA  = MNSOCT(N) - 1 + 3 * (NP-1)
         ELSE
            IA  = MNSOCT(N) - 1 + 3 * (NBSOCT(N)-NP)
         ENDIF
         IN  = 3*N1+3*(NP-2)
         DO 240 I=1,3
            COSOBO(IN+I) = RMCN(IA+I)
  240   CONTINUE
  250 CONTINUE
C
C    * LE COTE 3
      N=NUCOTE(3)
      N3 = NBSOCT(N)
      DO 270 NP=1,NBSOCT(N)
         IF(NSENS(N).GE.0) THEN
            IA  = MNSOCT(N) - 1 + 3 * (NP-1)
         ELSE
            IA  = MNSOCT(N) - 1 + 3 * (NBSOCT(N)-NP)
         ENDIF
         IN  = 3*(N1+N2) + 3 * (NP-1)
         DO 260 I=1,3
            COSOBO(IN+I) = RMCN(IA+I)
  260   CONTINUE
  270 CONTINUE
C
C    * LE COTE 4
      N=NUCOTE(4)
      N4 = NBSOCT(N)-2
      DO 290 NP=2,NBSOCT(N)
         IF(NSENS(N).GE.0) THEN
            IA  = MNSOCT(N) - 1 + 3 * (NP-1)
         ELSE
            IA  = MNSOCT(N) - 1 + 3 * (NBSOCT(N)-NP)
         ENDIF
         IN  = 3*(N1+N2+N3) + 3 * (NP-2)
         DO 280 I=1,3
            COSOBO(IN+I) = RMCN(IA+I)
  280   CONTINUE
  290 CONTINUE
C
C     PASSAGE DANS LE PLAN
      NEUO = 3*(1-1)
      NEUI = 3*(2-1)
      NEUJ = 3*(N1+N2+N3+N4-1)
      DO 300 K=1,3
        POINTO(K) = COSOBO(NEUO+K)
        POINTI(K) = COSOBO(NEUI+K)
        POINTJ(K) = COSOBO(NEUJ+K)
  300 CONTINUE
      CALL DF3D2D(POINTO,POINTI,POINTJ,D2D3,IERR)
      IF (IERR.EQ.1) THEN
        NBLGRC(NRERR) = 2
        KERR(1) =  " LE PASSAGE DANS LE PLAN"
        KERR(2) =  " (O,x,y)  EST IMPOSSIBLE"
        CALL LEREUR
        IERR = 2
        RETURN
      ENDIF
      NBS1 = N1+N2+N3+N4+1
      CALL VERPLA(POINTO,D2D3,NBS1,1,COSOBO,IERR)
      IF (IERR.EQ.1) THEN
        NBLGRC(NRERR) = 1
        KERR(1) =  "LA SURFACE DEFINIE EST GAUCHE"
        CALL LEREUR
        IERR = 2
        RETURN
      ELSE
C
C       PASSAGE DANS LE PLAN
         DO 310 I=1,NBSOBO
            XYZ(1) = COSOBO(3*(I-1)+1)
            XYZ(2) = COSOBO(3*(I-1)+2)
            XYZ(3) = COSOBO(3*(I-1)+3)
            CALL CH3D3D( POINTO , D2D3 ,XYZ ,XY3 )
            COSOBO(3*(I-1)+1) = XY3(1)
            COSOBO(3*(I-1)+2) = XY3(2)
            COSOBO(3*(I-1)+3) = XY3(3)
310      CONTINUE
      ENDIF
C
C     VERIFICATION DE LA NON DEGENERESCENCE DU BORD
C
      LALP   = NBS1
      NADALP = 0
      CALL TNMCDC( 'REEL' , LALP , NADALP )
      CALL CALANG(NBS1,COSOBO,0,SENSFR,RMCN(NADALP))
      CALL TNMCDS( 'REEL' , LALP , NADALP )
      IF (SENSFR.EQ.0) THEN
        IF (NTESTI.EQ.0) THEN
          NSENS(LINSTR(NONSTR)) = -NSENS(LINSTR(NONSTR))
          NTESTI = 1
          GOTO 210
        ELSE
          NBLGRC(NRERR) = 1
          KERR(1) =  "LE BORD DEFINI EST DEGENERE"
          CALL LEREUR
          IERR = 3
          RETURN
        ENDIF
      ENDIF
      RETURN
      END
