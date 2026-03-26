      SUBROUTINE A2NSSD( NTDL,ALFA1,MU1,A1,
     +                        ALFA2,MU2,A2,  A3)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : FORMER LA MATRICE A3 COMBINAISON LINEAIRE DES MATRICES A1 ET A2
C ----- A3 = ALFA1 * A1 + ALFA2 * A2  MATRICES STOCKEES PROFIL
C       A1 NON SYMETRIQUE A2 SYMETRIQUE A3 NON SYMETRIQUE
C       TOUTES STOCKEES SOUS FORME PROFIL
C
C ENTREES:
C --------
C NTDL  : ORDRE DES MATRICES A
C ALFA1 : COEFFICIENT REEL DOUBLE PRECISION
C MU1   : MU1(1)=0 MU1(I+1)=ADRESSE DANS A1 DU I-EME COEFFICIENT DIAGONA
C A1    : MATRICE NON SYMETRIQUE STOCKEE SOUS FORME PROFIL
C ALFA2 : COEFFICIENT REEL DOUBLE PRECISION
C MU2   : MU2(1)=0 MU2(I+1)=ADRESSE DANS A2 DU I-EME COEFFICIENT DIAGONA
C A2    : MATRICE SYMETRIQUE STOCKEE SOUS FORME PROFIL
C
C SORTIE:
C -------
C A3    : MATRICE NON SYMETRIQUE STOCKEE SOUS FORME PROFIL
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 1979
C2345X7..............................................................012
      DOUBLE PRECISION ALFA1,ALFA2,A1,A2,A3,S
      DIMENSION MU1(*),A1(*),MU2(*),A2(*),A3(*)
C
      IA1   = 0
      IA2   = 0
      IA3   = 0
      MU1I  = 0
      MU2I  = 0
      DO 1 I=1,NTDL
            I1 = I + 1
            MU1I1 = MU1(I1)
            MU2I1 = MU2(I1)
            N1 = (MU1I1 - MU1I) / 2
            N2 =  MU2I1 - MU2I  - 1
            N3 = MAX0(N1,N2)
            L1 = I - N1
            L2 = I - N2
            IF(L1 .LE. L2) GOTO 2
C
C           LA LIGNE DE A1 EST PLUS COURTE QUE CELLE DE A2
C
            L = L1 - 1
            DO 3 J=L2,L
                 IA2 = IA2 + 1
                 IA3 = IA3 + 1
                 A3(IA3     ) = ALFA2 * A2(IA2)
                 A3(IA3 + N3) = A3(IA3)
    3       CONTINUE
            GOTO 24
C
C           LA LIGNE DE A2 EST PLUS COURTE QUE CELLE DE A1
C
    2       IF(L1 .EQ. I) GOTO 7
            L = L2 - 1
            IF(L .LT. L1) GOTO 24
            DO 4 J=L1,L
                 IA1 = IA1 + 1
                 IA3 = IA3 + 1
                 A3(IA3   ) = ALFA1 * A1(IA1   )
                 A3(IA3+N3) = ALFA1 * A1(IA1+N1)
    4       CONTINUE
C
C           LA PARTIE COMMUNE A LA LIGNE I DE A1 ET A2
C
   24       L = L + 1
            IF(L.EQ.I) GOTO 7
            I1 = I - 1
            DO 5 J=L,I1
                 IA1 = IA1 + 1
                 IA2 = IA2 + 1
                 IA3 = IA3 + 1
                 S = ALFA2 * A2(IA2)
                 A3(IA3   ) = ALFA1 * A1(IA1   ) + S
                 A3(IA3+N3) = ALFA1 * A1(IA1+N1) + S
    5       CONTINUE
C
C     LE COEFFICIENT DIAGONAL
C
    7       IA1 = MU1I1
            IA2 = MU2I1
            IA3 = IA3 + N3 + 1
            A3(IA3) = ALFA1 * A1(IA1) + ALFA2 * A2(IA2)
            MU1I = MU1I1
            MU2I = MU2I1
    1 CONTINUE
C
      END
