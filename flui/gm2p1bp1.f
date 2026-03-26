      SUBROUTINE GM2P1BP1( AE, AG )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    STOCKAGE DES LIGNES DE LA MATRICE ELEMENTAIRE DANS LA
C -----    MATRICE GLOBALE DES DL BARYCENTRES DES TRIANGLES
C          BREZZI-FORTIN
C
C ENTREE :
C --------
C AE     : MATRICE DE MASSE (11,11) STOCKEE TRIANGULAIRE SUPERIEURE
C
C SORTIE :
C --------
C AG     : LES 2 LIGNES DE AE DES DEGRES DE LIBERTE ELIMINES
C          SONT STOCKES DANS AG
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC Saint Pierre du Perray Decembre 2008
C23456---------------------------------------------------------------012
      DOUBLE PRECISION  AE(66),  AG(11,2)
C
C     LA LIGNE 4 DE LA MATRICE ELEMENTAIRE EST STOCKEE DANS AG(*,1)
      M = 6
      DO 10 K=1,4
         M = M + 1
         AG(K,1) = AE(M)
 10   CONTINUE
      DO 15 K=5,11
         M = M + K - 1
         AG(K,1) = AE(M)
 15   CONTINUE
C
C     LA LIGNE 8 DE LA MATRICE ELEMENTAIRE EST STOCKEE DANS AG(*,2)
      M = 28
      DO 20 K=1,8
         M = M + 1
         AG(K,2) = AE(M)
 20   CONTINUE
      DO 25 K=9,11
         M = M + K - 1
         AG(K,2) = AE(M)
 25   CONTINUE
C
C     COMPRESSION DANS AE DE 11 A SEULEMENT 9 DEGRES DE LIBERTE
C     REDUCTION DES DL VITESSE 4 8
      M  = 0
      M1 = 0
      DO 60 I = 1, 11
         IF( I .EQ. 4 ) THEN
            M = M + 4
         ELSE IF( I .EQ. 8 ) THEN
            M = M + 8
         ELSE
            DO 50 J = 1, I
               M = M + 1
               IF( J .NE. 4 .AND. J .NE. 8 ) THEN
                  M1 = M1 + 1
                  AE(M1) = AE(M)
               ENDIF
 50         CONTINUE
         ENDIF
 60   CONTINUE
C
      RETURN
      END
