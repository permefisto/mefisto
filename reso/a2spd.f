      SUBROUTINE A2SPD ( NTDL, ALFA1,MU1,A1,
     +                         ALFA2,MU2,A2,  A3)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : FORMER LA MATRICE A3 COMBINAISON LINEAIRE DES MATRICES A1 ET A2
C ----- A3 = ALFA1 * A1 + ALFA2 * A2  MATRICES STOCKEES PROFIL
C       A1 ET A2 SONT SYMETRIQUES DE PROFILS  DIFFERENTS
C
C ENTREES:
C --------
C NTDL  : ORDRE DES MATRICES A
C ALFA1 : COEFFICIENT REEL DOUBLE PRECISION
C MU1   : MU1(1)=0 MU1(I+1)=ADRESSE DANS A1 DU I-EME COEFFICIENT DIAGONA
C A1    : MATRICE SYMETRIQUE STOCKEE SOUS FORME PROFIL
C ALFA2 : COEFFICIENT REEL DOUBLE PRECISION
C MU2   : MU2(1)=0 MU2(I+1)=ADRESSE DANS A2 DU I-EME COEFFICIENT DIAGONA
C A2    : MATRICE SYMETRIQUE STOCKEE SOUS FORME PROFIL
C
C SORTIE:
C -------
C A3    : MATRICE SYMETRIQUE STOCKEE SOUS FORME PROFIL
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 1979
C2345X7..............................................................012
      DOUBLE PRECISION  ALFA1,ALFA2,A1,A2,A3
      DIMENSION         MU1(*),A1(*),MU2(*),A2(*), A3(*)
C
      IA1   = 0
      IA2   = 0
      IA3   = 0
      MU1I  = 0
      MU2I  = 0
      DO 11 I=1,NTDL
            I1 = I + 1
            MU1I1 = MU1(I1)
            MU2I1 = MU2(I1)
            L1 = I1 - MU1I1 + MU1I
            L2 = I1 - MU2I1 + MU2I
            IF(L1 .LE. L2) GOTO 12
C
C           LA LIGNE DE A1 EST PLUS COURTE QUE CELLE DE A2
            L = L1 - 1
            DO 13 J=L2,L
                 IA2 = IA2 + 1
                 IA3 = IA3 + 1
                 A3(IA3) = ALFA2 * A2(IA2)
   13       CONTINUE
            GOTO 14
C
C           LA LIGNE DE A2 EST PLUS COURTE QUE CELLE DE A1
   12       L = L2 - 1
            IF(L1 .EQ. I .OR. L.LT. L1) GOTO 14
            DO 15 J=L1,L
                 IA1 = IA1 + 1
                 IA3 = IA3 + 1
                 A3(IA3) = ALFA1 * A1(IA1)
   15       CONTINUE
C
C           LA PARTIE COMMUNE A LA LIGNE I DE A1 ET A2
   14       L = L + 1
            DO 16 J=L,I
                 IA1 = IA1 + 1
                 IA2 = IA2 + 1
                 IA3 = IA3 + 1
                 A3(IA3) = ALFA1 * A1(IA1) + ALFA2 * A2(IA2)
   16       CONTINUE
            MU1I = MU1I1
            MU2I = MU2I1
   11 CONTINUE
C
      END
