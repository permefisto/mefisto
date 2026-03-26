      SUBROUTINE PLONAD(NSTOC,I1,I2,NBLOCI,NBLOCJ,IJT,GS,GN,AE)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : PLONGER DANS LA MATRICE AE SYMETRIQUE STOCKEE TRIANGULAIREMENT
C ---   UNE MATRICE GN(I1,I2) SI NSTOC = -1
C       UNE MATRICE GS(I1,I1) SI NSTOC =  1   STOCKEE TRIANGULAIREMENT
C       DE HAUT EN BAS ET DE LA GAUCHE VERS LA DROITE  1  2  4  ...
C                                                         3  5  ...
C                                                            6  ...
C ENTREES:
C --------
C NSTOC  : 1 SI GS(1 A I1*(I1+1)/2) , -1 SI GN(I1,I2)
C I1     : NOMBRE DE LIGNES DE GS(N)
C I2     : NOMBRE DE COLONNES DE GN SI NSTOC = -1
C NBLOCI : NUMERO DU BLOC DE LIGNES I DANS LA MATRICE AE CORRESPONDANT
C          AUX I1 LIGNES DE GS(N)
C NBLOCJ : NUMERO DU BLOC DE COLONNES DANS LA MATRICE AE CORRESPONDANT
C          AUX I2 COLONNES DE GN
C IJT    : PERMUTATION DES D.L. PAR COMPOSANTE => PAR NOEUD
C GS     : MATRICE POUR NSTOC = 1
C GN     : MATRICE POUR NSTOC =-1
C
C SORTIES:
C --------
C AE     : MATRICE TRIANGULAIRE SUPERIEURE SOMMEE A GS(N)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     OCTOBRE 1980
C23456---------------------------------------------------------------012
      DOUBLE PRECISION  GS(1:*),GN(I1,I2),AE(1:*)
      INTEGER           IJT(1:*)
C
      NI = I1 * (NBLOCI - 1)
      NJ = I2 * (NBLOCJ - 1)
      IF( NSTOC.LE. 0 ) GOTO 5
C
C     GS SYMETRIQUE STOCKEE TRIANGULAIREMENT
C     ======================================
      IG = 0
      DO 4 J=1,I1
         JJJ = IJT( J + NJ )
         DO 3 I=1,J
            II = IJT( I + NI )
            JJ = JJJ
            IF ( II .LE. JJ ) GO TO 2
            K = II
            II= JJ
            JJ = K
    2       K = JJ * ( JJ - 1 ) / 2 + II
            IG = IG + 1
            AE( K ) = AE( K ) + GS( IG )
    3    CONTINUE
    4 CONTINUE
      RETURN
C
C     MATRICE NON SYMETRIQUE  GN(II,I2)
C     =================================
    5 IG = I1
      DO 8 J=1,I2
         IF ( NBLOCI .EQ. NBLOCJ ) IG= J
         JJJ = IJT( J + NJ )
         DO 7 I= 1,IG
            II = IJT( I + NI )
            JJ = JJJ
            IF ( II .LE. JJ ) GO TO 6
            K = II
            II= JJ
            JJ = K
    6       K = JJ * ( JJ - 1 ) / 2 + II
            AE( K ) = AE( K ) + GN( I , J )
    7    CONTINUE
    8 CONTINUE
C
      RETURN
      END
