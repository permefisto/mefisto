      SUBROUTINE VERINJ(NBS1,NBS2,COSO,COSOBO,NBPEXT,NBPFRT,NEUEXT)
C***********************************************************************
C BUT :  VERIFICATION QUE TOUS LES POINTS SONT INTERIEURS AU MAILLAGE
C***********************************************************************
C
C ENTREES:
C           NBS1   : NOMBRE DE POINTS SUR LIGNES "HORIZONTALES"
C           NBS2   : NOMBRE DE POINTS SUR LIGNES "VERTICALES"
C           COSO   : COORDONNEES DE TOUS LES POINTS DU MAILLAGE
C
C TRAVAIL:
C           COSOBO : MATRICE DES COORDONNEES DE LA FRONTIERE
C
C SORTIES:
C           NBPEXT : NOMBRE DE POINTS EXTERIEURS AU MAILLAGE
C           NBPFRT : NOMBRE DE POINTS SUR LA FRONTIERE
C           NEUEXT : LISTE DES POINTS EXTERIEURS AU MAILLAGE
C                    (1 A NBPEXT) SUIVIE DE LA LISTE DES POINTS
C                    SUR LA FRONTIERE (NBPEXT+1 A NBPEXT+NBPFRT)
C
C***********************************************************************
C AUTEUR: CHRISTOPHE DOURSAT NOVEMBRE 1988
C****+7**************************************************************012
C***********************************************************************
C             DECLARATION DES VARIABLES ET DES TABLEAUX
C                  ET INITIALISATION DES VARIABLES
C***********************************************************************
      COMMON / UNITES / LECTEU , IMPRIM , INTERA , NUNITE(29)
      DIMENSION COSO(3,NBS1*NBS2)
      DIMENSION COSOBO(3,2*(NBS1+NBS2)-3)
      DIMENSION NEUEXT(NBS1*NBS2)
      DIMENSION A(2,2),B(2)
      LOGICAL TESTFR
      DATA EPS/0.0001/
C***********************************************************************
C                RECHERCHE D UN NOEUD EXTERIEUR AU MAILLAGE
C***********************************************************************
      X0 = COSOBO(1,1)
      Y0 = COSOBO(2,1)
      Y1 = COSOBO(2,1)
      DO 10 I=2,2*(NBS1+NBS2-2)
        X0 = AMIN1(X0,COSOBO(1,I))
        Y0 = AMIN1(Y0,COSOBO(2,I))
        Y1 = AMAX1(Y1,COSOBO(2,I))
   10 CONTINUE
C***********************************************************************
C      LE NOEUD DE COORDONNEES X0 Y0 EST EXTERIEUR AU MAILLAGE
C***********************************************************************
      X0 = X0-1.
      Y0 = (Y0+Y1)/2
      NUMMEM = 0
C***********************************************************************
C    RECHERCHE DU NOMBRE DE POINTS D INTERSECTION DU SEGMENT (M0,M)
C    AVEC LA FRONTIERE; ON PASSE EN REVUE LES SEGMENT (M(K),M(K+1))
C
C                            M(K) * BET=0.
C                                  \
C                                   \
C                                    \
C       ALP=0.               ALP=1.   \
C         *--------------------*       \
C         M0                   M        \
C                                 M(K+1) *----------*
C                                      BET=1.     M(K+2)
C
C LE POINT D INTERSECTION AURA POUR COORDONNEES BARYCENTRIQUES ALP BET
C***********************************************************************
      NBAR   = 2*(NBS1+NBS2)-4
      NBPEXT = 0
      NBPFRT = 0
C-----------------------------------------------------------------------
C           BOUCLE SUR LES NOEUDS INTERNES AU MAILLAGE
C-----------------------------------------------------------------------
      DO 100 J=2,NBS2-1
        DO 90 I=2,NBS1-1
          NE = (J-1)*NBS1+I
          NBPINT = 0
          ITEST  = 0
          TESTFR = .FALSE.
C-----------------------------------------------------------------------
C             BOUCLE SUR LES SEGMENTS DE LA FRONTIERE
C    ET CALCUL DES COORDONNEES BARYCENTRIQUES DU POINT D INTERSECTION
C-----------------------------------------------------------------------
   50     DO 60 K=1,NBAR
            A(1,1) = COSO(1,NE)-X0
            A(1,2) = COSOBO(1,K)-COSOBO(1,K+1)
            A(2,1) = COSO(2,NE)-Y0
            A(2,2) = COSOBO(2,K)-COSOBO(2,K+1)
            B(1)   = COSOBO(1,K)-X0
            B(2)   = COSOBO(2,K)-Y0
            DET = A(1,1)*A(2,2)-A(1,2)*A(2,1)
            IF (DET.NE.0.) THEN
              ALP = (B(1)*A(2,2)-B(2)*A(1,2))/DET
              BET = (A(1,1)*B(2)-A(2,1)*B(1))/DET
            ELSE
C-----------------------------------------------------------------------
C        ON EST DANS LE CAS OU LES DEUX SEGMENTS SONT PARALELLES
C      ON CALCULE LES DEUX AUTRES DETERMINANTS. SI ILS SONT TOUS LES
C   DEUX NULS LES SEGMENTS SONT PORTES PAR LA MEME DROITE ET ON DEPLACE
C ALORS M0. SINON LES DEUX SEGMENTS SONT STRICT. // ET PAS D INTERSEC.
C-----------------------------------------------------------------------
              DET1 = B(1)*A(2,2)-B(2)*A(1,2)
              DET2 = A(1,1)*B(2)-A(2,1)*B(1)
              IF ((DET1.EQ.0.).AND.(DET2.EQ.0.)) THEN
                DIST = SQRT(A(1,1)**2+A(2,1)**2)
                A(1,1) = A(1,1)/DIST
                A(2,1) = A(2,1)/DIST
                IF (A(2,1).LE.0.) THEN
                  X0 = X0+A(2,1)
                  Y0 = Y0-A(1,1)
                ELSE
                  X0 = X0-A(2,1)
                  Y0 = Y0+A(1,1)
                ENDIF
                NBPINT = 0
                ITEST  = 0
                GOTO 50
              ELSE
                ALP = 2.
                BET = 2.
              ENDIF
            ENDIF
C-----------------------------------------------------------------------
C       CAS OU LE POINT D INTERSECTION EST SUR LE SEGMENT (M0,M)
C-----------------------------------------------------------------------
            IF ((ALP.GT.0.).AND.(ALP.LT.(1.+EPS))) THEN
              NBPIM = NBPINT
              IF ((BET.GE.EPS).AND.(BET.LE.(1.-EPS))) THEN
                NBPINT = NBPINT+1
              ELSE
C-----------------------------------------------------------------------
C                  ON CHERCHE DANS QUEL CAS ON EST
C
C                         *
C                  *     /                       *
C                   \   /                         \
C                    \ /                           \
C          *----------*-----------*   *-------------*---------*
C         M0                     NE  M0            /         NE
C                                                 /
C                                                *
C
C       PAS DE POINT D'INTERSECTION      POINT D'INTERSECTION
C-----------------------------------------------------------------------
C               CAS OU LE POINT D INTERSECTION EST M(K+1)
C      ON GARDE EN MEMOIRE LE NUMERO DU NOEUD POUR NE PAS COMPTER
C    L INTERSECTION 2 FOIS SI ON LA RETROUVE POUR LE SEGMENT SUIVANT
C-----------------------------------------------------------------------
                IF (ABS(BET-1.).LT.EPS) THEN
                  IF (K.NE.NBAR) THEN
                    DIST = SQRT(A(1,1)**2+A(2,1)**2)
                    B(1) = A(1,1)/DIST
                    B(2) = A(2,1)/DIST
                    A(1,1) = COSOBO(1,K+1)-COSOBO(1,K)
                    A(2,1) = COSOBO(2,K+1)-COSOBO(2,K)
                    A(1,2) = COSOBO(1,K+2)-COSOBO(1,K+1)
                    A(2,2) = COSOBO(2,K+2)-COSOBO(2,K+1)
                    PSCAL  = A(1,1)*B(1)+A(2,1)*B(2)
                    A(1,1) = A(1,1)-PSCAL*B(1)
                    A(2,1) = A(2,1)-PSCAL*B(2)
                    PSCAL  = A(1,2)*B(1)+A(2,2)*B(2)
                    A(1,2) = A(1,2)-PSCAL*B(1)
                    A(2,2) = A(2,2)-PSCAL*B(2)
                    PSCAL  = A(1,1)*A(1,2)+A(2,1)*A(2,2)
                    IF (PSCAL.GT.0) THEN
                      NBPINT = NBPINT+1
                    ENDIF
                    NUMMEM = K
                  ELSE
                    IF (ITEST.NE.1) THEN
                      DIST = SQRT(A(1,1)**2+A(2,1)**2)
                      B(1) = A(1,1)/DIST
                      B(2) = A(2,1)/DIST
                      A(1,1) = COSOBO(1,K+1)-COSOBO(1,K)
                      A(2,1) = COSOBO(2,K+1)-COSOBO(2,K)
                      A(1,2) = COSOBO(1,2)-COSOBO(1,1)
                      A(2,2) = COSOBO(2,2)-COSOBO(2,1)
                      PSCAL  = A(1,1)*B(1)+A(2,1)*B(2)
                      A(1,1) = A(1,1)-PSCAL*B(1)
                      A(2,1) = A(2,1)-PSCAL*B(2)
                      PSCAL  = A(1,2)*B(1)+A(2,2)*B(2)
                      A(1,2) = A(1,2)-PSCAL*B(1)
                      A(2,2) = A(2,2)-PSCAL*B(2)
                      PSCAL  = A(1,1)*A(1,2)+A(2,1)*A(2,2)
                      IF (PSCAL.GT.0) THEN
                        NBPINT = NBPINT+1
                      ENDIF
                    ENDIF
                  ENDIF
                ENDIF
C-----------------------------------------------------------------------
C               CAS OU LE POINT D INTERSECTION EST M(K)
C-----------------------------------------------------------------------
                IF (ABS(BET).LT.EPS) THEN
                  IF (K.EQ.1) THEN
                    DIST = SQRT(A(1,1)**2+A(2,1)**2)
                    B(1) = A(1,1)/DIST
                    B(2) = A(2,1)/DIST
                    A(1,1) = COSOBO(1,NBAR+1)-COSOBO(1,NBAR)
                    A(2,1) = COSOBO(2,NBAR+1)-COSOBO(2,NBAR)
                    A(1,2) = COSOBO(1,2)-COSOBO(1,1)
                    A(2,2) = COSOBO(2,2)-COSOBO(2,1)
                    PSCAL  = A(1,1)*B(1)+A(2,1)*B(2)
                    A(1,1) = A(1,1)-PSCAL*B(1)
                    A(2,1) = A(2,1)-PSCAL*B(2)
                    PSCAL  = A(1,2)*B(1)+A(2,2)*B(2)
                    A(1,2) = A(1,2)-PSCAL*B(1)
                    A(2,2) = A(2,2)-PSCAL*B(2)
                    PSCAL  = A(1,1)*A(1,2)+A(2,1)*A(2,2)
                    IF (PSCAL.GT.0) THEN
                      NBPINT = NBPINT+1
                      ITEST  = 1
                    ENDIF
                  ELSE
                    IF (K.NE.(NUMMEM+1)) THEN
                      DIST = SQRT(A(1,1)**2+A(2,1)**2)
                      B(1) = A(1,1)/DIST
                      B(2) = A(2,1)/DIST
                      A(1,1) = COSOBO(1,K)-COSOBO(1,K-1)
                      A(2,1) = COSOBO(2,K)-COSOBO(2,K-1)
                      A(1,2) = COSOBO(1,K+1)-COSOBO(1,K)
                      A(2,2) = COSOBO(2,K+1)-COSOBO(2,K)
                      PSCAL  = A(1,1)*B(1)+A(2,1)*B(2)
                      A(1,1) = A(1,1)-PSCAL*B(1)
                      A(2,1) = A(2,1)-PSCAL*B(2)
                      PSCAL  = A(1,2)*B(1)+A(2,2)*B(2)
                      A(1,2) = A(1,2)-PSCAL*B(1)
                      A(2,2) = A(2,2)-PSCAL*B(2)
                      PSCAL  = A(1,1)*A(1,2)+A(2,1)*A(2,2)
                      IF (PSCAL.GT.0) THEN
                        NBPINT = NBPINT+1
                        ITEST  = 1
                      ENDIF
                    ENDIF
                  ENDIF
                ENDIF
              ENDIF
C----------------------------------------------------------------------
C            CAS OU LE POINT M EST SUR LE SEGMENT FRONTIERE
C-----------------------------------------------------------------------
              IF ((NBPINT.NE.NBPIM).AND.(ABS(ALP-1).LT.EPS)) THEN
                NBPINT = NBPIM
                TESTFR = .TRUE.
                GOTO 70
              ENDIF
            ENDIF
   60     CONTINUE
   70     CONTINUE
          IF (TESTFR) THEN
            NBPFRT = NBPFRT+1
            NEUEXT(NBPEXT+NBPFRT) = NE
          ELSE
            IF (MOD(NBPINT,2).NE.1) THEN
              DO 80 K=NBPFRT,1,-1
                NEUEXT(NBPEXT+K+1) = NEUEXT(NBPEXT+K)
   80         CONTINUE
              NBPEXT = NBPEXT+1
              NEUEXT(NBPEXT) = NE
            ENDIF
          ENDIF
   90   CONTINUE
  100 CONTINUE
C***********************************************************************
C           IMPRESSION DES NOEUDS EXTERIEURS AU MAILLAGE
C***********************************************************************
      IF ((NBPEXT.NE.0).OR.(NBPFRT.NE.0)) THEN
        WRITE (IMPRIM,1000)
        WRITE (IMPRIM,1001)
        WRITE (IMPRIM,1002) NBPEXT
        WRITE (IMPRIM,1003) NBPFRT
        WRITE (IMPRIM,1000)
      ENDIF
1000  FORMAT('**************************************')
1001  FORMAT('*        TOPOLOGIE DES NOEUDS        *')
1002  FORMAT('*   ',I4,' POINT(S) A L EXTERIEUR      *')
1003  FORMAT('*  ',I4,' POINT(S) SUR LA FRONTIERE    *')
      RETURN
      END
