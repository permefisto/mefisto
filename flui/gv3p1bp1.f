      SUBROUTINE GV3P1BP1( AE, AGAUSS )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    ELIMINATION PAR GAUSS DES DEGRES DE LIBERTE VITESSE 5 10 15
C -----    (VITESSE1 2 3 en moyenne sur e)
C          INTERNES ATTACHES AU TETRAEDRE BREZZI-FORTIN
C
C MODIFIES:
C --------
C AE     : MATRICE DE VISCOSITE (19,19) STOCKEE TRIANGULAIRE SUPERIEURE
C
C SORTIES:
C --------
C AE     : MATRICE DE VISCOSITE (16,16) STOCKEE TRIANGULAIRE SUPERIEURE
C AGAUSS : LES 3 LIGNES DE AE DES DEGRES DE LIBERTE ELIMINES
C          POUR LES CALCULER APRES RESOLUTION
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC Saint Pierre du Perray Decembre 2008
C23456---------------------------------------------------------------012
      DOUBLE PRECISION AE(190),  AGAUSS(19,3)
      DOUBLE PRECISION A(19,19), PIVOT, S
C
C     MISE SOUS FORME CARREE DE LA MATRICE ELEMENTAIRE
      M = 0
      DO 20 J=1,19
         DO 10 I=1,J
            M = M + 1
            A(I,J) = AE(M)
            A(J,I) = AE(M)
 10      CONTINUE
 20   CONTINUE
C
C     LA LIGNE DE AE DU DL 5 ELIMINE PAR GAUSS
C     POUR PERMETTRE LE CALCUL DES VITESSES MOYENNES APRES RESOLUTION
      DO 50 J=1,19
         AGAUSS(J,1) = A(5,J)
 50   CONTINUE
C
C     TRIANGULATION DE GAUSS SUR LA LIGNE 5
      PIVOT = A(5,5)
      DO 55 I=1,19
         IF( I .NE. 5 ) THEN
            S = A(I,5) / PIVOT
            DO 52 J=1,19
ccc               IF( J .NE. 5 ) THEN
               A(I,J) = A(I,J) - S * A(5,J)
ccc               ENDIF
 52         CONTINUE
         ENDIF
 55   CONTINUE
C
C     LA LIGNE DE AE DU DL 10 ELIMINE PAR GAUSS
C     POUR PERMETTRE LE CALCUL DES VITESSES MOYENNES APRES RESOLUTION
      DO 100 J=1,19
         AGAUSS(J,2) = A(10,J)
 100  CONTINUE
C     EN FAIT A(10,5)=0 SUITE A GAUSS MAIS VALEUR EXACTE IMPOSEE
      AGAUSS(5,2) = 0D0
C
C     TRIANGULATION DE GAUSS SUR LA LIGNE 10
      PIVOT = A(10,10)
      DO 105 I=1,19
         IF( I .NE. 5 .AND. I .NE. 10 ) THEN
            S = A(I,10) / PIVOT
            DO 102 J=1,19
ccc               IF( J .NE. 5 .AND. J .NE. 10 ) THEN
               A(I,J) = A(I,J) - S * A(10,J)
ccc               ENDIF
 102        CONTINUE
         ENDIF
 105  CONTINUE
C
C     LA LIGNE DE AE DU DL 15 ELIMINE PAR GAUSS
C     POUR PERMETTRE LE CALCUL DES VITESSES MOYENNES APRES RESOLUTION
      DO 150 J=1,19
         AGAUSS(J,3) = A(15,J)
 150  CONTINUE
C     EN FAIT A(15,5)=0  SUITE A GAUSS MAIS VALEUR EXACTE IMPOSEE
      AGAUSS(5,3) = 0D0
C     EN FAIT A(15,10)=0 SUITE A GAUSS MAIS VALEUR EXACTE IMPOSEE
      AGAUSS(10,3) = 0D0
C
C     TRIANGULATION DE GAUSS SUR LA LIGNE 15
      PIVOT = A(15,15)
      DO 155 I=1,19
         IF( I .NE. 5 .AND. I .NE. 10 .AND. I .NE. 15 ) THEN
            S = A(I,15) / PIVOT
            DO 152 J=1,19
ccc               IF( J .NE. 5 .AND. J .NE. 10 .AND. J .NE. 15 ) THEN
               A(I,J) = A(I,J) - S * A(15,J)
ccc               ENDIF
 152        CONTINUE
         ENDIF
 155   CONTINUE
C
C     RECONSTRUCTION DANS AE AVEC SEULEMENT 16 DEGRES DE LIBERTE
C     REDUCTION DES DL VITESSE 5 10 15
      M  = 0
      I4 = 0
      DO 200 I = 1, 16
         J4 = 0
         DO 190 J = 1, I
            M = M + 1
            AE(M) = A(I+I4,J+J4)
            J4 = J / 4
 190     CONTINUE
         I4 = I / 4
 200  CONTINUE
C
      RETURN
      END
