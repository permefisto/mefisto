      SUBROUTINE INTDPDP( NBPOID, POIDS, NDIM, NBPOL, DPOLYP, DPDP )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  CALCUL DE L'INTEGRALE DP DP dx dy dz sur l'EF de REFERENCE
C -----  DES GRADIENTS DES POLYNOMES LAGRANGE
C      x DES GRADIENTS DES POLYNOMES LAGRANGE
C      ( NECESSAIRE POUR LE TETRAEDRE DE TAYLOR HOOD )
C
C ENTREES:
C --------
C NBPOID : NOMBRE DE  POIDS DE LA FORMULE D'INTEGRATION NUMERIQUE
C POIDS  : VALEUR DES POIDS DE LA FORMULE D'INTEGRATION NUMERIQUE
C NDIM   : DIMENSION DE L'ESPACE OU NOMBRE DE COORDONNEES
C NBPOL  : NOMBRE DE POLYNOMES
C DPOLYP : DPOLYP(NDIM,NBPOLY,NPI)
C          DPOLYP( I  , J    , L ) = DPJ/DXI (XL,YL)
C
C SORTIE :
C --------
C DPDP   : INTEGRALE DPi/DXk DPj/DXl dX SUR LE TETRAEDRE UNITE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Laboratoire J-L.LIONS UPMC Paris    Juin 2007
C23456---------------------------------------------------------------012
      DOUBLE PRECISION POIDS(NBPOID), S
      DOUBLE PRECISION DPOLYP(NDIM,NBPOL,NBPOID)
      DOUBLE PRECISION DPDP(NDIM,NBPOL,NDIM,NBPOL)
C
      DO 50 J=1,NBPOL
         DO 40 L=1,NDIM
            DO 30 I=1,NBPOL
               DO 20 K=1,NDIM
                  S = 0D0
                  DO 10 M=1,NBPOID
                     S = S + POIDS(M) * DPOLYP(K,I,M) * DPOLYP(L,J,M)
 10               CONTINUE
                  IF( ABS(S) .LT. 1D-14 ) S=0D0
                  DPDP(K,I,L,J) = S
 20            CONTINUE
 30         CONTINUE
 40      CONTINUE
 50   CONTINUE
C
cccC     AFFICHAGE DES INTEGRALES DP2 DP2 dX
ccc      PRINT 10000
ccc10000 FORMAT(//'INTEGRALE EXACTE des DP2 DP2 dX pour le TETRAEDRE de TAYL
ccc     %OR-HOOD')
ccc      DO 90 J=1,NBPOL
ccc         DO 80 L=1,NDIM
ccc            DO 70 I=1,NBPOL
ccc               PRINT 10001, (K,I,L,J,DPDP(K,I,L,J),K=1,NDIM)
ccc 70         CONTINUE
ccc 80      CONTINUE
ccc 90   CONTINUE
ccc10001 FORMAT(2(4I3,1X,D25.17,4X))
C
      RETURN
      END
