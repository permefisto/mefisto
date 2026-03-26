      SUBROUTINE A2DPD ( NTDL, ALFA1D,    A1,
     +                         ALFA2D,MU2,A2,
     +                                    A3  )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : FORMER LA MATRICE A3 COMBINAISON LINEAIRE DES MATRICES A1 ET A2
C ----- A3 = ALFA1D * A1 + ALFA2D * A2  MATRICES STOCKEES PROFIL
C       A1 EST DIAGONALE ET A2 N EST PAS DIAGONALE
C       => A3 A LE PROFIL DE A2
C
C ENTREES:
C --------
C NTDL   : ORDRE DES MATRICES A
C ALFA1D : COEFFICIENT REEL DOUBLE PRECISION
C A1     : MATRICE DIAGONALE STOCKEE DANS UN VECTEUR
C ALFA2D : COEFFICIENT REEL DOUBLE PRECISION
C MU2    : MU2(1)=0 MU2(I+1)=ADRESSE DANS A2 DU I-EME COEFFICIENT DIAGONAL
C A2     : MATRICE STOCKEE SOUS FORME PROFIL
C
C SORTIE:
C -------
C A3    : MATRICE STOCKEE SOUS FORME PROFIL
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 1979
C2345X7..............................................................012
      DOUBLE PRECISION  ALFA1D,ALFA2D,A1(*),A2(*),A3(*)
      INTEGER           MU2(*)
C
      IA2   = 0
      IA3   = 0
      MU2I  = 0
C
      DO 1 I=1,NTDL
         MU2I1 = MU2(I + 1)
C
C        LES COEFFICIENTS NON DIAGONAUX
         DO 2 J = 1, MU2I1-MU2I-1
            IA2 = IA2 + 1
            IA3 = IA3 + 1
            A3(IA3) = ALFA2D * A2(IA2)
    2    CONTINUE
C
C        LE COEFFICIENT DIAGONAL I
         IA2 = IA2 + 1
         IA3 = IA3 + 1
         A3(IA3) = ALFA1D * A1(I) + ALFA2D * A2(IA2)
         MU2I = MU2I1
C
    1 CONTINUE
C
      END
