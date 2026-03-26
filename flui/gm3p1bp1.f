      SUBROUTINE GM3P1BP1( AE, AG )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    STOCKAGE DES LIGNES DE LA MATRICE ELEMENTAIRE DANS LA
C -----    MATRICE GLOBALE DES DL BARYCENTRES DES TETRAEDRES
C          BREZZI-FORTIN
C
C ENTREE :
C --------
C AE     : MATRICE DE MASSE (19,19) STOCKEE TRIANGULAIRE SUPERIEURE
C
C SORTIE :
C --------
C AG     : LES 3 LIGNES DE AE DES DEGRES DE LIBERTE ELIMINES
C          SONT STOCKES DANS AG
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC Saint Pierre du Perray Decembre 2008
C23456---------------------------------------------------------------012
      DOUBLE PRECISION  AE(190),  AG(19,3)
C
C     LA LIGNE 5 DE LA MATRICE ELEMENTAIRE
      M = 10
      DO 10 K=1,5
         M = M + 1
         AG(K,1) = AE(M)
 10   CONTINUE
      DO 15 K=6,19
         M = M + K - 1
         AG(K,1) = AE(M)
 15   CONTINUE
C
C     LA LIGNE 10 DE LA MATRICE ELEMENTAIRE
      M = 45
      DO 20 K=1,10
         M = M + 1
         AG(K,2) = AE(M)
 20   CONTINUE
      DO 25 K=11,19
         M = M + K - 1
         AG(K,2) = AE(M)
 25   CONTINUE
C
C     LA LIGNE 15 DE LA MATRICE ELEMENTAIRE
      M = 105
      DO 30 K=1,15
         M = M + 1
         AG(K,3) = AE(M)
 30   CONTINUE
      DO 35 K=16,19
         M = M + K - 1
         AG(K,3) = AE(M)
 35   CONTINUE
C
C     COMPRESSION DANS AE DE 19 A SEULEMENT 16 DEGRES DE LIBERTE
C     REDUCTION DES DL VITESSE 5 10 15
      M  = 0
      M1 = 0
      DO 60 I = 1, 19
         IF( I .EQ. 5 ) THEN
            M = M + 5
         ELSE IF( I .EQ. 10 ) THEN
            M = M + 10
         ELSE IF( I .EQ. 15 ) THEN
            M = M + 15
         ELSE
            DO 50 J = 1, I
               M = M + 1
               IF( J .NE. 5 .AND. J .NE. 10 .AND. J .NE. 15 ) THEN
                  M1 = M1 + 1
                  AE(M1) = AE(M)
               ENDIF
 50         CONTINUE
         ENDIF
 60   CONTINUE
C
      RETURN
      END
