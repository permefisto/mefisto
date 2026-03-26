      SUBROUTINE GA3P1BP1( AE, BE, AGAUSS )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    ELIMINATION PAR GAUSS DES DEGRES DE LIBERTE 5 10 15 VITESSE
C -----    CAR INTERNES AU TETRAEDRE BREZZI-FORTIN
C
C MODIFIES:
C --------
C AE     : MATRICE DE VISCOSITE (19,19) STOCKEE TRIANGULAIRE SUPERIEURE
C BE     : SECOND MEMBRE ELEMENTAIRE (19)
C
C SORTIES:
C --------
C AE     : MATRICE DE VISCOSITE (16,16) STOCKEE TRIANGULAIRE SUPERIEURE
C BE     : SECOND MEMBRE ELEMENTAIRE (16)
C AGAUSS : LES 3 LIGNES DE AE et BE DES DEGRES DE LIBERTE ELIMINES
C          POUR LES CALCULER APRES RESOLUTION
CC+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC Saint Pierre du Perray Decembre 2008
C23456---------------------------------------------------------------012
      DOUBLE PRECISION AE(190),  BE(19), AGAUSS(20,3)
      DOUBLE PRECISION A(19,19), PIVOT
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
C     LA LIGNE DE AE et BE DU DL 5 ELIMINE PAR GAUSS
C     POUR PERMETTRE LE CALCUL DES VITESSES MOYENNES APRES RESOLUTION
      DO 50 J=1,19
         AGAUSS(J,1) = A(5,J)
 50   CONTINUE
      AGAUSS(20,1) = BE(5)
C
C     TRIANGULATION DE GAUSS SUR LA LIGNE 5
      PIVOT = A(5,5)
      DO 55 I=1,19
         IF( I .NE. 5 ) THEN
            DO 52 J=1,19
               IF( J .NE. 5 ) THEN
                  A(I,J) = A(I,J) - A(I,5) * A(5,J) / PIVOT
               ENDIF
 52         CONTINUE
            BE(I) = BE(I) - A(I,5) * BE(5) / PIVOT
         ENDIF
 55   CONTINUE
C
C     LA LIGNE DE AE et BE DU DL 10 ELIMINE PAR GAUSS
C     POUR PERMETTRE LE CALCUL DES VITESSES MOYENNES APRES RESOLUTION
      DO 60 J=1,19
         AGAUSS(J,2) = A(10,J)
 60   CONTINUE
      AGAUSS(20,2) = BE(10)
C     EN FAIT A(10,5)=0 SUITE A GAUSS MAIS CALCUL DANS A NON FAIT
      AGAUSS(5,2) = 0D0
C
C     TRIANGULATION DE GAUSS SUR LA LIGNE 10
      PIVOT = A(10,10)
      DO 68 I=1,19
         IF( I .NE. 5 .AND. I .NE. 10 ) THEN
            DO 65 J=1,19
               IF( J .NE. 5 .AND. J .NE. 10 ) THEN
                  A(I,J) = A(I,J) - A(I,10) * A(10,J) / PIVOT
               ENDIF
 65         CONTINUE
            BE(I) = BE(I) - A(I,10) * BE(10) / PIVOT
         ENDIF
 68   CONTINUE
C
C     LA LIGNE DE AE et BE DU DL 15 ELIMINE PAR GAUSS
C     POUR PERMETTRE LE CALCUL DES VITESSES MOYENNES APRES RESOLUTION
      DO 70 J=1,19
         AGAUSS(J,3) = A(15,J)
 70   CONTINUE
      AGAUSS(20,3) = BE(15)
C     EN FAIT A(15,5)=0 SUITE A GAUSS MAIS CALCUL DANS A NON FAIT
      AGAUSS(5,3) = 0D0
C     EN FAIT A(15,10)=0 SUITE A GAUSS MAIS CALCUL DANS A NON FAIT
      AGAUSS(10,3) = 0D0
C
C     TRIANGULATION DE GAUSS SUR LA LIGNE 15
      PIVOT = A(15,15)
      DO 85 I=1,19
         IF( I .NE. 5 .AND. I .NE. 10 .AND. I .NE. 15 ) THEN
            DO 80 J=1,19
               IF( J .NE. 5 .AND. J .NE. 10 .AND. J .NE. 15 ) THEN
                  A(I,J) = A(I,J) - A(I,15) * A(15,J) / PIVOT
               ENDIF
 80         CONTINUE
            BE(I) = BE(I) - A(I,15) * BE(15) / PIVOT
         ENDIF
 85   CONTINUE
C
ccc      print *,'AE et BE apres Gauss et avant reduction'
ccc      do 92 i=1,19
ccc         print 10110,(i,j,A(i,j),j=1,19)
ccc 92   continue
ccc      print *,('  be',i,'=',BE(i),i=1,19)
C
C     RECONSTRUCTION DANS AE et BE AVEC SEULEMENT 16 DEGRES DE LIBERTE
C     REDUCTION DES DL VITESSE 5 10 15
      M  = 0
      I4 = 0
      DO 100 I = 1, 16
         J4 = 0
         DO 95 J = 1, I
            M = M + 1
            AE(M) = A(I+I4,J+J4)
            J4 = J / 4
 95      CONTINUE
         BE(I) = BE(I+I4)
         I4 = I / 4
 100  CONTINUE
C
ccc      print *
ccc      print *,'AE apres Gauss dl 5 10 15 et apres compression'
ccc      m = 0
ccc      do 110 i=1,16
ccc         print 10110,(i,j,AE(m+j),j=1,i)
ccc         m = m + i
ccc 110  continue
ccc10110 FORMAT(5('  AE',I2,I2,'=',G14.6))
ccc      print *,'BE apres Gauss et apres compression'
ccc      print *,('  BE',i,'=',BE(i),i=1,16)
C
      RETURN
      END
