      SUBROUTINE GV2P1BP1( AE, AGAUSS )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    ELIMINATION PAR GAUSS DES DEGRES DE LIBERTE VITESSE 4 et 8
C -----    (VITESSE1 et 2 en moyenne sur e)
C          INTERNES ATTACHES AU TRIANGLE BREZZI-FORTIN
C
C MODIFIES:
C --------
C AE     : MATRICE DE VISCOSITE (11,11) STOCKEE TRIANGULAIRE SUPERIEURE
C
C SORTIES:
C --------
C AE     : MATRICE DE VISCOSITE (9,9) STOCKEE TRIANGULAIRE SUPERIEURE
C AGAUSS : LES 2 LIGNES DE AE et BE DES DEGRES DE LIBERTE ELIMINES
C          POUR LES CALCULER APRES RESOLUTION
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC Saint Pierre du Perray Decembre 2008
C23456---------------------------------------------------------------012
      DOUBLE PRECISION AE(66),   AGAUSS(11,2)
      DOUBLE PRECISION A(11,11), PIVOT, S
C
C     MISE SOUS FORME CARREE DE LA MATRICE ELEMENTAIRE
      M = 0
      DO 20 J=1,11
         DO 10 I=1,J
            M = M + 1
            A(I,J) = AE(M)
            A(J,I) = AE(M)
 10      CONTINUE
 20   CONTINUE
C
C     LA LIGNE DE AE DU DL 4 ELIMINE PAR GAUSS EST STOCKEE DANS AGAUSS
C     POUR PERMETTRE LE CALCUL DE LA VITESSE AU BARYCENTRE APRES RESOLUTION
      DO 40 J=1,11
         AGAUSS(J,1) = A(4,J)
 40   CONTINUE
C
C     TRIANGULATION DE GAUSS SUR LA LIGNE 4
      PIVOT = A(4,4)
      DO 45 I=1,11
         IF( I .NE. 4 ) THEN
            S = A(I,4) / PIVOT
            DO 42 J=1,11
               IF( J .NE. 4 ) THEN
                  A(I,J) = A(I,J) - S * A(4,J)
               ENDIF
 42         CONTINUE
         ENDIF
 45   CONTINUE
C
C     LA LIGNE DE AE DU DL 8 ELIMINE PAR GAUSS EST STOCKEE DANS AGAUSS
C     POUR PERMETTRE LE CALCUL DES VITESSES MOYENNES APRES RESOLUTION
      DO 80 J=1,11
         AGAUSS(J,2) = A(8,J)
 80   CONTINUE
C     EN FAIT A(8,4)=0 SUITE A GAUSS MAIS CALCUL DANS A NON FAIT
      AGAUSS(4,2) = 0D0
C
C     TRIANGULATION DE GAUSS SUR LA LIGNE 8
      PIVOT = A(8,8)
      DO 85 I=1,11
         IF( I .NE. 4 .AND. I .NE. 8 ) THEN
            S = A(I,8) / PIVOT
            DO 82 J=1,11
               IF( J .NE. 4 .AND. J .NE. 8 ) THEN
                  A(I,J) = A(I,J) - S * A(8,J)
               ENDIF
 82         CONTINUE
         ENDIF
 85   CONTINUE
C
C     RECONSTRUCTION DANS AE AVEC SEULEMENT 9 DEGRES DE LIBERTE
C     REDUCTION DES DL VITESSE 4 8
      M  = 0
      I3 = 0
      DO 100 I = 1, 9
         J3 = 0
         DO 95 J = 1, I
            M = M + 1
            AE(M) = A(I+I3,J+J3)
            J3 = J / 3
 95      CONTINUE
         I3 = I / 3
 100  CONTINUE
C
      RETURN
      END
