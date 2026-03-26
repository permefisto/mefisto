      SUBROUTINE T43LAG( NBPOLY, X,      NBPOF,  NOPOF,
     &                   POLYF,  DPOLYF, NPIF,   POIDSF,
     &                   NOOBSF, NUMISU, NUMASU, NBJEUX, JEU, LTDESU,
     &                   IECHAN, PENALI,
     &                   NOPART, BE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   CONTRIBUTION D UNE FACE AU SECONDS MEMBRE ELEMENTAIRE
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
C IECHAN : 1 SI EXISTE UN ECHANGE SUR CETTE FACE
C          2 SI EXISTE UN CONTACT SUR CETTE FACE
C PENALI : COEFFICIENT DE PENALISATION DE LA COONDITION DE DIRICHLET
C          NON 0D0 SI LE TMS "CONTACT" EST A PENALISER ET PENALI DOIT ETRE GRAND
C          0D0 PAS DE PENALISATION DU TMS "CONTACT"
C NOPART : POUR NLSE SEULEMENT AU NIVEAU DE SOURCE=FORCE et CONTACT=FIXATION
C          1 SI PARTIE REELLE TRAITEE ou 2 SI PARTIE IMAGINAIRE TRAITEE
C          0 SI INACTIF (CAS THERMIQUE STANDARD D'UNE SOURCE)
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
     &                  DPOLYF(2,NBPOF,NPIF),
     &                  POIDSF(NPIF),
     &                  PENALI,
     &                  FLUXN(3),
     &                  BE(NBPOLY),
     &                  PROSCD
C
      DOUBLE PRECISION  DELTA,
     &                  GL(3),
     &                  DGL(2,3),
     &                  DGLN,
     %                  VN(3)
C
C     BOUCLE SUR LES NPIF POINTS D INTEGRATION DE LA FACE
C     ---------------------------------------------------
      IF( TESTNL .GT. 0 ) THEN
C        RECUPERATION DE LA TEMPERATURE AUX NBPOF DL DE LA FACE
         MN = (MNTHET-1)/2
         DO 10 I=1,NBPOF
            DMCN((MNTHDL-1)/2+I)=
     &      DMCN( MN+MCN(MNNODL+NOPOF(I)-1) )
 10      CONTINUE
      ENDIF
C
      DO 100 L=1,NPIF
C
         IF( TESTNL .GT. 0 ) THEN
C           CALCUL DE LA TEMPERATURE AU POINT D'INTEGRATION L
            TEMPEL=PROSCD( POLYF(1,L), MCN(MNTHDL), NBPOF )
         ENDIF
C
C        RECHERCHE DES COORDONNEES DU POINT ET DU JACOBIEN DE G
         CALL E23LAG( NBPOLY, NBPOF, NOPOF,
     &                POLYF(1,L), DPOLYF(1,1,L), X,
     &                GL, DGL, DELTA )
C        EN SORTIE GL=LES 3 COORDONNEES DU POINT D'INTEGRATION L
C
C        LE VECTEUR NORMAL UNITAIRE A LA FACE
         CALL VECNOR( DGL, DGLN, VN )
C
         IF( IECHAN .EQ. 1 ) THEN
            IF( TESTNL .LE. 5 ) THEN
C              IL EXISTE UNE SOURCE SUR CETTE FACE
C              CALCUL DU FLUX NORMAL
               CALL RESOUR( 3, NOOBSF, 3, GL,
     &                      LTDESU(LPSOUR,JEU,NOOBSF), FLUXN )
            ELSE
C              FORCE(2) REQUISE
C              LE VECTEUR NORMAL UNITAIRE EST UTILISE
               CALL REFORC( 3, NOOB, 3, GL(1), GL(2), GL(3),
     %                                  VN(1), VN(2), VN(3),
     %                      LTDESU(LPSOUR,JEU,NOOB), FLUXN )
               FLUXN(1) = FLUXN(NOPART)
            ENDIF
         ELSE
            IF( TESTNL .LE. 5 ) THEN
C              IL EXISTE UN CONTACT SUR CETTE FACE
               CALL RECONT( 3, NOOBSF, 3, GL,
     %                      LTDESU(LPCONT,JEU,NOOBSF), FLUXN )
               FLUXN(1) = FLUXN(1) * PENALI
            ELSE
C              FIXATION(2) PENALISEE
               MN =  LTDESU(LPCONT,JEU,NOOBSF)
               CALL REFIXA( 3, NOOB, GL(1), GL(2), GL(3), MN,
     %                      NBCOFI, FLUXN )
C              NBCOFI LE NOMBRE DE COMPOSANTES FIXEES
               DO I = 1, NBCOFI
C                 LE NUMERO DE LA COMPOSANTE FIXEE
                  NU = MCN( MN + WUCOFI - 1 + I )
                  IF( NU .EQ. NOPART ) THEN
                     FLUXN(1) = FLUXN(I) * PENALI
                  ENDIF
               ENDDO
            ENDIF
         ENDIF
C
C        PRODUIT AVEC LE POIDS
         DELTA = DELTA * POIDSF(L)
C
C        ASSEMBLAGE DANS BE
         DO 30 I=1,NBPOF
            II = NOPOF(I)
            BE(II) = BE(II) + DELTA * POLYF(I,L) * FLUXN(1)
 30      CONTINUE
C
100   CONTINUE
C
      RETURN
      END
