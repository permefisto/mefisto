      SUBROUTINE NSF3LAG( NBPOLY, X,      NBPOF,  NOPOF,
     %                    POLYF,  DPOLYF, NPIF,   POIDSF,
     %                    NOOBSF, NUMISU, NUMASU, NBJEUX, JEU, LTDESU,
     %                    IECHAN, PENALI,  BE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   NLSE: CONTRIBUTION D UNE FACE AU SECOND MEMBRE ELEMENTAIRE
C -----   D UN ELEMENT TRIDIMENSIONNEL LAGRANGE DEGRE 2 ISOPARAMETRIQUE
C
C ENTREES:
C --------
C NBPOLY : NBRE DE POINTS DE L ELEMENT TRIDIMENSIONNEL COMPLET
C X      : COORDONNEES DES POINTS DE L ELEMENT  (NBPOLY,3)
C NBPOF  : NBRE DE POINTS DEFINISSANT L APPLICATION G DE LA FACE
C NOPOF  : NO DANS L ELEMENT FINI DES NBPOF POINTS DE LA FACE
C POLYF  : POLYF(I,J)=VALEUR DU I-EME POLYNOME AU J-EME POINT D INTEGRAT
C DPOLYF : DPOLYF(I,J,L)=VALEUR DE LA I-EME DERIVEE DU J-EME
C                        POLYNOME AU L-EME POINT D INTEGRATION DE LA FACE
C NPIF   : NOMBRE DE POINTS D INTEGRATION SUR CETTE FACE
C POIDSF : LES NPIF POIDS DE LA FORMULE D INTEGRATION DE LA FACE
C NOOBSF : NUMERO DE SURFACE DE LA FACE DE L'EF
C NUMISU : NUMERO MINIMAL DES OBJETS SURFACES
C NUMASU : NUMERO MAXIMAL DES OBJETS SURFACES
C NBJEUX : NOMBRE DE JEUX DE DONNEES
C JEU    : NUMERO DU JEU  DE DONNEES POUR CE CALCUL DE LA MATRICE ELEMENTAIRE
C LTDESU : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES THERMIQUES
C          DES OBJETS SURFACES
C IECHAN : 1 SI EXISTE UNE FORCE    SUR CETTE FACE
C          2 SI EXISTE UNE FIXATION SUR CETTE FACE
C PENALI : COEFFICIENT DE PENALISATION DE LA COONDITION DE DIRICHLET
C          NON 0D0 SI LE TMS "CONTACT" EST A PENALISER ET PENALI DOIT ETRE GRAND
C          0D0 PAS DE PENALISATION DU TMS "CONTACT"
C
C SORTIES:
C --------
C BE     : LE SECOND MEMBRE ELEMENTAIRE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     OCTOBRE 1994
C23456---------------------------------------------------------------012
      include"./incl/donthe.inc"
      include"./incl/cthet.inc"
      include"./incl/cnonlin.inc"
      include"./incl/a___fixation.inc"
C
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE      (MCN(1), RMCN(1), DMCN(1))
C
      REAL              X(NBPOLY,3)
      INTEGER           NOPOF(NBPOF)
      INTEGER           LTDESU( 1:MXDOTH, 1:NBJEUX, NUMISU:NUMASU )
      DOUBLE PRECISION  POLYF(NBPOF,NPIF),
     %                  DPOLYF(2,NBPOF,NPIF),
     %                  POIDSF(NPIF),
     %                  PENALI,
     %                  FORCE(2), FIXA(2),
     %                  BE(NBPOLY,2),
     %                  PROSCD
C
      DOUBLE PRECISION  DELTA, D,
     %                  GL(3),
     %                  DGL(2,3),
     %                  DGLN,
     %                  VN(3),
     %                  U(9,2)
C
C     BOUCLE SUR LES NPIF POINTS D INTEGRATION DE LA FACE
C     ---------------------------------------------------
C     RECUPERATION DE L'ONDE AUX NBPOF DL DE LA FACE
      MN = (MNTHET-1)/2
      DO I=1,NBPOF
C        NUMERO DU DERNIER DL DU NOEUD I DE L'EF
         NU = MCN( MNNODL+NOPOF(I)-1 ) * 2
C        UG0(NUNOEUD,2) EST RANGE PAR NOEUDS
         U(I,1) = DMCN( MN + NU -1 )
         U(I,2) = DMCN( MN + NU    )
      ENDDO
C
      DO 100 L=1,NPIF
C
C        CALCUL DE LA TEMPERATURE AU POINT D'INTEGRATION L
         TEMPEL = PROSCD( POLYF(1,L), U(1,1), NBPOF )
         ONDEPI = PROSCD( POLYF(1,L), U(1,2), NBPOF )
C
C        RECHERCHE DES COORDONNEES DU POINT ET DU JACOBIEN DE G
         CALL E23LAG( NBPOLY, NBPOF, NOPOF,
     %                POLYF(1,L), DPOLYF(1,1,L), X,
     %                GL, DGL, DELTA )
C        EN SORTIE GL=LES 3 COORDONNEES DU POINT D'INTEGRATION L
C
C        LE VECTEUR NORMAL UNITAIRE A LA FACE
         CALL VECNOR( DGL, DGLN, VN )
C
         IF( IECHAN .EQ. 1 ) THEN
C
C           FORCE(2) REQUISE
C           LE VECTEUR NORMAL UNITAIRE EST REMPLACE PAR L'ONDE
            CALL REFORC( 3, NOOB, 2, GL(1), GL(2), GL(3),
     %                               TEMPEL, ONDEPI, 0D0,
ccc     %                               VN(1), VN(2), VN(3),
     %                   LTDESU(LPSOUR,JEU,NOOB), FORCE )
C
         ELSE
C
C           FIXATION(2) PENALISEE
            MN =  LTDESU(LPCONT,JEU,NOOBSF)
            CALL REFIXA( 3, NOOB, GL(1), GL(2), GL(3), MN,
     %                   NBCOFI, FIXA )
C           NBCOFI LE NOMBRE DE COMPOSANTES FIXEES
            DO I = 1, NBCOFI
C              LE NUMERO DE LA COMPOSANTE FIXEE
               NU = MCN( MN + WUCOFI - 1 + I )
               FORCE(NU) = FIXA(I) * PENALI
            ENDDO
C
         ENDIF
C
C        PRODUIT AVEC LE POIDS
         D = DELTA * POIDSF(L)
C
C        ASSEMBLAGE DANS BE
         DO I=1,NBPOF
            II = NOPOF(I)
            BE(II,1) = BE(II,1) + D * POLYF(I,L) * FORCE(1)
            BE(II,2) = BE(II,2) + D * POLYF(I,L) * FORCE(2)
         ENDDO
C
100   CONTINUE
C
      RETURN
      END
