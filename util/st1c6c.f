      SUBROUTINE ST1C6C( NBST1C, NB1C6C, NOST1C6C )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :FOURNIR LE NUMERO DES 2 SOMMETS DES 192 ARETES=1-CUBE D'UN 6-CUBE
C -----
C
C SORTIES:
C --------
C NBST1C   :   2=NOMBRE DE SOMMETS D'UNE ARETE=1-CUBE D'UN 6-CUBE
C NB1C6C   : 192=2**5 C16 = NOMBRE D'ARETES=1-CUBES D'UN 6-CUBE
C            soit 32 ARETES POUR CHACUNE DES DIRECTIONS Xi i=1,...,6
C NOST1C6C : NOST1C6C(I,J) NO DU I-EME SOMMET DE L'ARETE J D'UN 6-CUBE
C
C RESULTAT AFFICHE:
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: : ALAIN PERRONNET  Laboratoire JLL UPMC PARIS     OCTOBRE 2005
C2345X7..............................................................012
      INTEGER  NOST1C6C(2,192)
C
      NBST1C = 2
      NB1C6C = 192
C     NUMERO DE L'ARETE CONSTRUITE
      NA     = 0
C
C     ARETE 1 => I=0 et 1 POUR N=M=L=K=0 J=0
C     ARETE 2 => I=0 et 1 POUR N=M=L=K=0 J=1
C     ARETE = 2 SOMMETS DANS UNE DIRECTION FIXEE ( C16 possibles ) ET
C             LES 5 AUTRES DIRECTIONS SONT PARCOURUES (=> 2**5 possibles)
C
C     ARETES DANS LA DIRECTION X1
      DO 16 N=0,1
         DO 15 M=0,1
            DO 14 L=0,1
               DO 13 K=0,1
                  DO 12 J=0,1
C                          UNE ARETE DE PLUS
                           NA = NA + 1
C                          LE NO DU PREMIER SOMMET DE L'ARETE SELON X1
                           NOST1C6C( 1, NA )= N1ST6C( 0, J, K, L, M, N )
C                          LE NO DU SECOND SOMMET DE L'ARETE SELON X1
                           NOST1C6C( 2, NA )= N1ST6C( 1, J, K, L, M, N )
 12               CONTINUE
 13            CONTINUE
 14         CONTINUE
 15      CONTINUE
 16   CONTINUE
C
C     ARETES DANS LA DIRECTION X2
      DO 26 N=0,1
         DO 25 M=0,1
            DO 24 L=0,1
               DO 23 K=0,1
                     DO 21 I=0,1
C                          UNE ARETE DE PLUS
                           NA = NA + 1
C                          LE NO DU PREMIER SOMMET DE L'ARETE SELON X2
                           NOST1C6C( 1, NA )= N1ST6C( I, 0, K, L, M, N )
C                          LE NO DU SECOND SOMMET DE L'ARETE SELON X2
                           NOST1C6C( 2, NA )= N1ST6C( I, 1, K, L, M, N )
 21                  CONTINUE
 23            CONTINUE
 24         CONTINUE
 25      CONTINUE
 26   CONTINUE
C
C     ARETES DANS LA DIRECTION X3
      DO 36 N=0,1
         DO 35 M=0,1
            DO 34 L=0,1
                  DO 32 J=0,1
                     DO 31 I=0,1
C                          UNE ARETE DE PLUS
                           NA = NA + 1
C                          LE NO DU PREMIER SOMMET DE L'ARETE SELON X3
                           NOST1C6C( 1, NA )= N1ST6C( I, J, 0, L, M, N )
C                          LE NO DU SECOND SOMMET DE L'ARETE SELON X3
                           NOST1C6C( 2, NA )= N1ST6C( I, J, 1, L, M, N )
 31                  CONTINUE
 32               CONTINUE
 34         CONTINUE
 35      CONTINUE
 36   CONTINUE
C
C     ARETES DANS LA DIRECTION X4
      DO 46 N=0,1
         DO 45 M=0,1
               DO 43 K=0,1
                  DO 42 J=0,1
                     DO 41 I=0,1
C                          UNE ARETE DE PLUS
                           NA = NA + 1
C                          LE NO DU PREMIER SOMMET DE L'ARETE SELON X4
                           NOST1C6C( 1, NA )= N1ST6C( I, J, K, 0, M, N )
C                          LE NO DU SECOND SOMMET DE L'ARETE SELON X4
                           NOST1C6C( 2, NA )= N1ST6C( I, J, K, 1, M, N )
 41                  CONTINUE
 42               CONTINUE
 43            CONTINUE
 45      CONTINUE
 46   CONTINUE
C
C     ARETES DANS LA DIRECTION X5
      DO 56 N=0,1
            DO 54 L=0,1
               DO 53 K=0,1
                  DO 52 J=0,1
                     DO 51 I=0,1
C                          UNE ARETE DE PLUS
                           NA = NA + 1
C                          LE NO DU PREMIER SOMMET DE L'ARETE SELON X5
                           NOST1C6C( 1, NA )= N1ST6C( I, J, K, L, 0, N )
C                          LE NO DU SECOND SOMMET DE L'ARETE SELON X5
                           NOST1C6C( 2, NA )= N1ST6C( I, J, K, L, 1, N )
 51                  CONTINUE
 52               CONTINUE
 53            CONTINUE
 54         CONTINUE
 56   CONTINUE
C
C     ARETES DANS LA DIRECTION X6
         DO 65 M=0,1
            DO 64 L=0,1
               DO 63 K=0,1
                  DO 62 J=0,1
                     DO 61 I=0,1
C                          UNE ARETE DE PLUS
                           NA = NA + 1
C                          LE NO DU PREMIER SOMMET DE L'ARETE SELON X6
                           NOST1C6C( 1, NA )= N1ST6C( I, J, K, L, M, 0 )
C                          LE NO DU SECOND SOMMET DE L'ARETE SELON X6
                           NOST1C6C( 2, NA )= N1ST6C( I, J, K, L, M, 1 )
 61                  CONTINUE
 62               CONTINUE
 63            CONTINUE
 64         CONTINUE
 65      CONTINUE
C
      print *,'     NB ARETES=',NA
      RETURN
      END
