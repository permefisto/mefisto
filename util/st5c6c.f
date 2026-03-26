      SUBROUTINE ST5C6C( NBST5C, NB5C6C, NOST5C6C )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    FOURNIR LE NUMERO DES 32 SOMMETS DES FACES=5-CUBES D'UN 6-CUBE
C -----
C
C SORTIES:
C --------
C NBST5C : 32=NOMBRE DE SOMMETS D'UNE FACE=5-CUBE D'UN 6-CUBE
C NB5C6C : 12=NOMBRE DE FACES=5-CUBES             D'UN 6-CUBE
C NOST5C6C : NOST5C6C(I,J) NO DU I-EME SOMMET DE LA FACE J D'UN 6-CUBE
C
C RESULTAT AFFICHE:
C       1,  3,  5,  7,  9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31,
C      33, 35, 37, 39, 41, 43, 45, 47, 49, 51, 53, 55, 57, 59, 61, 63,
C       2,  4,  6,  8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32,
C      34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62, 64,
C       1,  2,  5,  6,  9, 10, 13, 14, 17, 18, 21, 22, 25, 26, 29, 30,
C      33, 34, 37, 38, 41, 42, 45, 46, 49, 50, 53, 54, 57, 58, 61, 62,
C       3,  4,  7,  8, 11, 12, 15, 16, 19, 20, 23, 24, 27, 28, 31, 32,
C      35, 36, 39, 40, 43, 44, 47, 48, 51, 52, 55, 56, 59, 60, 63, 64,
C       1,  2,  3,  4,  9, 10, 11, 12, 17, 18, 19, 20, 25, 26, 27, 28,
C      33, 34, 35, 36, 41, 42, 43, 44, 49, 50, 51, 52, 57, 58, 59, 60,
C       5,  6,  7,  8, 13, 14, 15, 16, 21, 22, 23, 24, 29, 30, 31, 32,
C      37, 38, 39, 40, 45, 46, 47, 48, 53, 54, 55, 56, 61, 62, 63, 64,
C       1,  2,  3,  4,  5,  6,  7,  8, 17, 18, 19, 20, 21, 22, 23, 24,
C      33, 34, 35, 36, 37, 38, 39, 40, 49, 50, 51, 52, 53, 54, 55, 56,
C       9, 10, 11, 12, 13, 14, 15, 16, 25, 26, 27, 28, 29, 30, 31, 32,
C      41, 42, 43, 44, 45, 46, 47, 48, 57, 58, 59, 60, 61, 62, 63, 64,
C       1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16,
C      33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48,
C      17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32,
C      49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64,
C       1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16,
C      17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32,
C      33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48,
C      49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64,
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: : ALAIN PERRONNET  TEXAS A & M UNIVERSITY         15 JULY 2005
C2345X7..............................................................012
      INTEGER  NOST5C6C(32,12)
C
      NBST5C = 32
      NB5C6C = 12
C
C     FACE 1 => I=1 EST A SUPPRIMER (<=> X=1 N'APPARTIENT PAS A LA FACE)
      NUFACE = 1
      NSC = 0
      NSF = 0
      DO 16 N=0,1
         DO 15 M=0,1
            DO 14 L=0,1
               DO 13 K=0,1
                  DO 12 J=0,1
                     DO 11 I=0,1
C                       LE NO DE SOMMET DU 6-CUBE
                        NSC = NSC + 1
                        IF( I .EQ. 0 ) THEN
C                          LE SOMMET EST SUR LA FACE 1
                           NSF = NSF + 1
                           NOST5C6C( NSF, NUFACE ) = NSC
                        ENDIF
 11                  CONTINUE
 12               CONTINUE
 13            CONTINUE
 14         CONTINUE
 15      CONTINUE
 16   CONTINUE
C
C     FACE 2 => I=0 EST A SUPPRIMER (<=> X=0 N'APPARTIENT PAS A LA FACE)
      NUFACE = 2
      NSC = 0
      NSF = 0
      DO 26 N=0,1
         DO 25 M=0,1
            DO 24 L=0,1
               DO 23 K=0,1
                  DO 22 J=0,1
                     DO 21 I=0,1
C                       LE NO DE SOMMET DU 6-CUBE
                        NSC = NSC + 1
                        IF( I .EQ. 1 ) THEN
C                          LE SOMMET EST SUR LA FACE 2
                           NSF = NSF + 1
                           NOST5C6C( NSF, NUFACE ) = NSC
                        ENDIF
 21                  CONTINUE
 22               CONTINUE
 23            CONTINUE
 24         CONTINUE
 25      CONTINUE
 26   CONTINUE
C
C     FACE 3 => J=1 EST A SUPPRIMER (<=> Y=1 N'APPARTIENT PAS A LA FACE)
      NUFACE = 3
      NSC = 0
      NSF = 0
      DO 36 N=0,1
         DO 35 M=0,1
            DO 34 L=0,1
               DO 33 K=0,1
                  DO 32 J=0,1
                     DO 31 I=0,1
C                       LE NO DE SOMMET DU 6-CUBE
                        NSC = NSC + 1
                        IF( J .EQ. 0 ) THEN
C                          LE SOMMET EST SUR LA FACE 3
                           NSF = NSF + 1
                           NOST5C6C( NSF, NUFACE ) = NSC
                        ENDIF
 31                  CONTINUE
 32               CONTINUE
 33            CONTINUE
 34         CONTINUE
 35      CONTINUE
 36   CONTINUE
C
C     FACE 4 => J=0 EST A SUPPRIMER (<=> Y=0 N'APPARTIENT PAS A LA FACE)
      NUFACE = 4
      NSC = 0
      NSF = 0
      DO 46 N=0,1
         DO 45 M=0,1
            DO 44 L=0,1
               DO 43 K=0,1
                  DO 42 J=0,1
                     DO 41 I=0,1
C                       LE NO DE SOMMET DU 6-CUBE
                        NSC = NSC + 1
                        IF( J .EQ. 1 ) THEN
C                          LE SOMMET EST SUR LA FACE 4
                           NSF = NSF + 1
                           NOST5C6C( NSF, NUFACE ) = NSC
                        ENDIF
 41                  CONTINUE
 42               CONTINUE
 43            CONTINUE
 44         CONTINUE
 45      CONTINUE
 46   CONTINUE
C
C     FACE 5 => K=1 EST A SUPPRIMER (<=> Z=1 N'APPARTIENT PAS A LA FACE)
      NUFACE = 5
      NSC = 0
      NSF = 0
      DO 56 N=0,1
         DO 55 M=0,1
            DO 54 L=0,1
               DO 53 K=0,1
                  DO 52 J=0,1
                     DO 51 I=0,1
C                       LE NO DE SOMMET DU 6-CUBE
                        NSC = NSC + 1
                        IF( K .EQ. 0 ) THEN
C                          LE SOMMET EST SUR LA FACE 5
                           NSF = NSF + 1
                           NOST5C6C( NSF, NUFACE ) = NSC
                        ENDIF
 51                  CONTINUE
 52               CONTINUE
 53            CONTINUE
 54         CONTINUE
 55      CONTINUE
 56   CONTINUE
C
C     FACE 6 => K=0 EST A SUPPRIMER (<=> Z=0 N'APPARTIENT PAS A LA FACE)
      NUFACE = 6
      NSC = 0
      NSF = 0
      DO 66 N=0,1
         DO 65 M=0,1
            DO 64 L=0,1
               DO 63 K=0,1
                  DO 62 J=0,1
                     DO 61 I=0,1
C                       LE NO DE SOMMET DU 6-CUBE
                        NSC = NSC + 1
                        IF( K .EQ. 1 ) THEN
C                          LE SOMMET EST SUR LA FACE 6
                           NSF = NSF + 1
                           NOST5C6C( NSF, NUFACE ) = NSC
                        ENDIF
 61                  CONTINUE
 62               CONTINUE
 63            CONTINUE
 64         CONTINUE
 65      CONTINUE
 66   CONTINUE
C
C     FACE 7 => L=1 EST A SUPPRIMER (<=> U=1 N'APPARTIENT PAS A LA FACE)
      NUFACE = 7
      NSC = 0
      NSF = 0
      DO 76 N=0,1
         DO 75 M=0,1
            DO 74 L=0,1
               DO 73 K=0,1
                  DO 72 J=0,1
                     DO 71 I=0,1
C                       LE NO DE SOMMET DU 6-CUBE
                        NSC = NSC + 1
                        IF( L .EQ. 0 ) THEN
C                          LE SOMMET EST SUR LA FACE 7
                           NSF = NSF + 1
                           NOST5C6C( NSF, NUFACE ) = NSC
                        ENDIF
 71                  CONTINUE
 72               CONTINUE
 73            CONTINUE
 74         CONTINUE
 75      CONTINUE
 76   CONTINUE
C
C     FACE 8 => L=0 EST A SUPPRIMER (<=> U=0 N'APPARTIENT PAS A LA FACE)
      NUFACE = 8
      NSC = 0
      NSF = 0
      DO 86 N=0,1
         DO 85 M=0,1
            DO 84 L=0,1
               DO 83 K=0,1
                  DO 82 J=0,1
                     DO 81 I=0,1
C                       LE NO DE SOMMET DU 6-CUBE
                        NSC = NSC + 1
                        IF( L .EQ. 1 ) THEN
C                          LE SOMMET EST SUR LA FACE 8
                           NSF = NSF + 1
                           NOST5C6C( NSF, NUFACE ) = NSC
                        ENDIF
 81                  CONTINUE
 82               CONTINUE
 83            CONTINUE
 84         CONTINUE
 85      CONTINUE
 86   CONTINUE
C
C     FACE 9 => M=1 EST A SUPPRIMER (<=> V=1 N'APPARTIENT PAS A LA FACE)
      NUFACE = 9
      NSC = 0
      NSF = 0
      DO 96 N=0,1
         DO 95 M=0,1
            DO 94 L=0,1
               DO 93 K=0,1
                  DO 92 J=0,1
                     DO 91 I=0,1
C                       LE NO DE SOMMET DU 6-CUBE
                        NSC = NSC + 1
                        IF( M .EQ. 0 ) THEN
C                          LE SOMMET EST SUR LA FACE 9
                           NSF = NSF + 1
                           NOST5C6C( NSF, NUFACE ) = NSC
                        ENDIF
 91                  CONTINUE
 92               CONTINUE
 93            CONTINUE
 94         CONTINUE
 95      CONTINUE
 96   CONTINUE
C
C     FACE 10 => M=0 EST A SUPPRIMER (<=> V=0 N'APPARTIENT PAS A LA FACE)
      NUFACE = 10
      NSC = 0
      NSF = 0
      DO 106 N=0,1
         DO 105 M=0,1
            DO 104 L=0,1
               DO 103 K=0,1
                  DO 102 J=0,1
                     DO 101 I=0,1
C                       LE NO DE SOMMET DU 6-CUBE
                        NSC = NSC + 1
                        IF( M .EQ. 1 ) THEN
C                          LE SOMMET EST SUR LA FACE 10
                           NSF = NSF + 1
                           NOST5C6C( NSF, NUFACE ) = NSC
                        ENDIF
 101                 CONTINUE
 102              CONTINUE
 103           CONTINUE
 104        CONTINUE
 105     CONTINUE
 106  CONTINUE
C
C     FACE 11 => N=1 EST A SUPPRIMER (<=> W=1 N'APPARTIENT PAS A LA FACE)
      NUFACE = 11
      NSC = 0
      NSF = 0
      DO 116 N=0,1
         DO 115 M=0,1
            DO 114 L=0,1
               DO 113 K=0,1
                  DO 112 J=0,1
                     DO 111 I=0,1
C                       LE NO DE SOMMET DU 6-CUBE
                        NSC = NSC + 1
                        IF( N .EQ. 0 ) THEN
C                          LE SOMMET EST SUR LA FACE 11
                           NSF = NSF + 1
                           NOST5C6C( NSF, NUFACE ) = NSC
                        ENDIF
 111                 CONTINUE
 112              CONTINUE
 113           CONTINUE
 114        CONTINUE
 115     CONTINUE
 116  CONTINUE
C
C     FACE 12 => N=0 EST A SUPPRIMER (<=> W=0 N'APPARTIENT PAS A LA FACE)
      NUFACE = 12
      NSC = 0
      NSF = 0
      DO 126 N=0,1
         DO 125 M=0,1
            DO 124 L=0,1
               DO 123 K=0,1
                  DO 122 J=0,1
                     DO 121 I=0,1
C                       LE NO DE SOMMET DU 6-CUBE
                        NSC = NSC + 1
                        IF( N .EQ. 1 ) THEN
C                          LE SOMMET EST SUR LA FACE 12
                           NSF = NSF + 1
                           NOST5C6C( NSF, NUFACE ) = NSC
                        ENDIF
 121                 CONTINUE
 122              CONTINUE
 123           CONTINUE
 124        CONTINUE
 125     CONTINUE
 126  CONTINUE
C
      RETURN
      END
