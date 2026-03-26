      SUBROUTINE T33LAG( NBPOLY, NBPOF,  NOPOF,
     &                   POLYF,  DPOLYF, X,
     &                   NPIF,   POIDSF,
     &                   NOOBSU, NUMISU, NUMASU, NBJEUX, JEU, LTDESU,
     &                   IECHAN, PENALI, CONFAC)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   CONTRIBUTION D UNE FACE A LA MATRICE DE CONDUCTIVITE
C -----   D UN ELEMENT TRIDIMENSIONNEL LAGRANGE DEGRE 2 ISOPARAMETRIQUE
C
C ENTREES:
C --------
C NBPOLY : NBRE DE POINTS DE L ELEMENT TRIDIMENSIONNEL COMPLET
C NBPOF  : NBRE DE POINTS DEFINISSANT L APPLICATION G DE LA FACE
C NOPOF  : NO DANS L ELEMENT DES NBPOF POINTS DE LA FACE
C POLYF  : POLYF(I,J)=VALEUR DU I-EME POLYNOME AU J-EME POINT D INTEGRAT
C DPOLYF : DPOLYF(I,J,L)=VALEUR DE LA I-EME DERIVEE DU J-EME
C                        POLYNOME AU L-EME POINT D INTEGRATION DE LA FACE
C X      : COORDONNEES DES POINTS DE L ELEMENT  (NBPOLY,3)
C NPIF   : NOMBRE DE POINTS D INTEGRATION SUR CETTE FACE
C POIDSF : LES NPIF POIDS DE LA FORMULE D INTEGRATION DE LA FACE
C NOOBSU : NUMERO DE L'OBJET SURFACE DE LA FACE DE L'EF
C NUMISU : NUMERO MINIMAL DES OBJETS SURFACES
C NUMASU : NUMERO MAXIMAL DES OBJETS SURFACES
C NBJEUX : NOMBRE DE JEUX DE DONNEES
C JEU    : NUMERO DU JEU  DE DONNEES POUR CE CALCUL DE LA MATRICE ELEMENTAIRE
C LTDESU : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES THERMIQUES
C          DES OBJETS SURFACES
C IECHAN : 1 SI EXISTE UN ECHANGE SUR CETTE FACE
C          2 SI EXISTE UN CONTACT SUR CETTE FACE
C PENALI : COEFFICIENT DE PENALISATION DE LA COONDITION DE DIRICHLET
C          NON 0D0 SI LE TMS "CONTACT" EST A PENALISER ET PENALI DOIT ETRE GRAND
C          0D0 PAS DE PENALISATION DU TMS "CONTACT"
C
C SORTIES:
C --------
C CONFAC : CONTRIBUTION DE LA FACE A LA MATRICE DE CONDUCTIVITE
C          MATRICE (NPIF,NPIF) REELLE DOUBLE PRECISION SYMETRIQUE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     OCTOBRE 1994
C23456---------------------------------------------------------------012
      include"./incl/donthe.inc"
      include"./incl/cthet.inc"
      include"./incl/cnonlin.inc"
C
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE      (MCN(1), RMCN(1), DMCN(1))
C
      REAL              X(NBPOLY,3)
      INTEGER           NOPOF(*)
      INTEGER           LTDESU( 1:MXDOTH, 1:NBJEUX, NUMISU:NUMASU )
      DOUBLE PRECISION  POLYF(NBPOF,NPIF),
     &                  DPOLYF(2,NBPOF,NPIF),
     &                  CONFAC(*),
     &                  POIDSF(NPIF),
     &                  PENALI,
     &                  PROSCD
C
      DOUBLE PRECISION  ECHANG,
     &                  DELTA,
     &                  GL(3),
     &                  DGL(2,3)
C
C     MISE A ZERO DU TABLEAU CONFAC
C     -----------------------------
      CALL AZEROD( NBPOF*(NBPOF+1)/2 , CONFAC )
C
      IF( TESTNL .GT. 0 ) THEN
C        RECUPERATION DE LA TEMPERATURE AUX NBPOF DL DE LA FACE
         MN = (MNTHET-1)/2
         DO 10 I=1,NBPOF
            DMCN((MNTHDL-1)/2+I)=
     &      DMCN( MN+MCN(MNNODL+(NOPOF(I)-1)) )
10       CONTINUE
      ENDIF
C
C     BOUCLE SUR LES NPIF POINTS D INTEGRATION DE LA FACE
C     ---------------------------------------------------
      DO 100 L=1,NPIF
C
C        RECHERCHE DES COORDONNEES DU POINT ET DU JACOBIEN DE G
C        ------------------------------------------------------
         CALL E23LAG( NBPOLY, NBPOF, NOPOF,
     &               POLYF(1,L), DPOLYF(1,1,L), X,
     &               GL, DGL, DELTA )
C
C        RECHERCHE DU COEFFICIENT D ECHANGE EN CE POINT DE LA FACE
C        ---------------------------------------------------------
         IF( TESTNL .GT. 0 ) THEN
C           CALCUL DE LA TEMPERATURE AU POINT D'INTEGRATION L
            TEMPEL=PROSCD( POLYF(1,L), MCN(MNTHDL), NBPOF )
         ENDIF
C
         IF( IECHAN .EQ. 1 ) THEN
C           IL EXISTE UN ECHANGE SUR CETTE FACE
            CALL REECHA( 3, NOOBSU, 3, GL,
     &                   LTDESU(LPECHA,JEU,NOOBSU), ECHANG )
         ELSE
C           IL EXISTE UN CONTACT SUR CETTE FACE
            ECHANG = PENALI
         ENDIF
C
C        CALCUL DE TRANSPOSEE(P(GL)) * (P(GL))
C        -------------------------------------
         DELTA = DELTA * POIDSF(L) * ECHANG
         K = 0
         DO 40 I=1,NBPOF
            DO 30 J=1,I
               K = K + 1
               CONFAC(K) = CONFAC(K) + DELTA * POLYF(I,L)
     &                                       * POLYF(J,L)
   30       CONTINUE
   40    CONTINUE
C
  100 CONTINUE
      END
