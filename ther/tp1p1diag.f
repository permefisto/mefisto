      SUBROUTINE TP1P1DIAG( NBSOM, XYZSOM, NBTRIA, NOSOTR, TPPDIAG )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LA MATRICE ELEMENTAIRE ET L'ASSEMBLER DANS LA
C -----    MATRICE DIAGONALE DE MASSE DU MAILLAGE DE TRIANGLES 2P1D
C
C ENTREES:
C --------
C NBSOM  : NOMBRE DE SOMMETS DU MAILLAGE
C XYZSOM : COORDONNEES DES NBSOM SOMMETS
C NBTRIA : NOMBRE DE TRIANGLES
C NOSOTR : NUMERO DES 3 SOMMETS DES NBTRIA TRIANGLES
C
C SORTIES:
C --------
C TPPDIAG: MATRICE DIAGONALE DE MASSE DU MAILLAGE DE TRIANGLES 2P1D
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Texas A & M University at QATAR  Fevrier 2011
C23456---------------------------------------------------------------012
      REAL              XYZSOM(3,NBSOM)
      DOUBLE PRECISION  TPPDIAG(NBSOM)
      INTEGER           NOSOTR(NBTRIA,3)
      DOUBLE PRECISION  X21, Y21, X31, Y31, X32, Y32, DELTA, DELTAS6
      INTEGER           NS1, NS2, NS3, NOSTR(3)
      EQUIVALENCE      (NOSTR(1),NS1), (NOSTR(2),NS2), (NOSTR(3),NS3)
C
C     MISE A ZERO DE LA MATRICE DIAGONALE
C     ===================================
      DO NS = 1, NBSOM
         TPPDIAG( NS ) = 0D0
      ENDDO
C
C     LA BOUCLE SUR LES TRIANGLES 2P1D DU MAILLAGE
C     ============================================
      DO NUTR = 1, NBTRIA
C
C        NUMERO DES 3 SOMMETS DU TRIANGLE
         NS1 = NOSOTR( NUTR, 1 )
         NS2 = NOSOTR( NUTR, 2 )
         NS3 = NOSOTR( NUTR, 3 )
C
C        DETERMINANT DE LA TRANSFORMATION F: e ref --> e
C        -----------------------------------------------
         X21 = XYZSOM(1,NS2) - XYZSOM(1,NS1)
         X31 = XYZSOM(1,NS3) - XYZSOM(1,NS1)
         X32 = XYZSOM(1,NS3) - XYZSOM(1,NS2)
C
         Y21 = XYZSOM(2,NS2) - XYZSOM(2,NS1)
         Y31 = XYZSOM(2,NS3) - XYZSOM(2,NS1)
         Y32 = XYZSOM(2,NS3) - XYZSOM(2,NS2)
C
         DELTA = ABS( X21 * Y31 - X31 * Y21 )
C
C        Integrale t[P][P] dX = TPP MATRICE DIAGONALE
C        --------------------------------------------
         DELTAS6 = DELTA / 6D0
ccc      tPP(1) = DELTAS6
ccc      tPP(2) = DELTAS6
ccc      tPP(3) = DELTAS6
C
C        ASSEMBLAGE DANS LE VECTEUR DE LA DIAGONALE DE LA MATRICE
         DO K=1,3
C
C           NUMERO DU SOMMET K DE L'EF
            NS = NOSTR( K )
C
C           ASSEMBLAGE DE LA MATRICE DIAGONALE DE MASSE
            TPPDIAG( NS ) = TPPDIAG( NS ) + DELTAS6
C
         ENDDO
C
C     FIN DE BOUCLE SUR LES EF DE TYPE NUTYEL
      ENDDO
C
      RETURN
      END
