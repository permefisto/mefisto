      SUBROUTINE QUADCOUL2DBORD( XY, COUL, NC )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACE LE QUADRANGLE DE SOMMETS XY ET DE COULEURS COUL
C -----    SELON LES COULEURS INTERMEDIAIRES  (PALETTE 11 RECOMMANDEE)
C          AINSI QUE LES ARETES DU QUADRANGLE SELON LA COULEUR NC
C          ATTENTION (X,Y) EN COORDONNEES OBJETS 2D

C ENTREES:
C --------
C XY     : 2 COORDONNEES DES 4 SOMMETS
C COUL   : NUMERO REEL DE LA COULEUR AUX 4 SOMMETS DU QUADRANGLE
C          REEL DE VALEURS COMPRISES ENTRE N1COUL ET NDCOUL
C          (cf ~/incl/trvari.inc)
C NC     : COULEUR DE TRACE DES 4 ARETES
C          SI NC<0 PAS DE TRACE DU BORD
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC  PARIS        MAI 1999
C2345X7..............................................................012
      REAL           XY(2,4),  COUL(4)
      REAL           XYT(2,3), COULT(3)

C     LE QUADRANGLE EST DIVISE EN 4 A PARTIR DU BARYCENTRE B
      DO I=1,2
         XYT(I,1) = (XY(I,1) + XY(I,2) + XY(I,3) + XY(I,4) ) * 0.25
      ENDDO
      COULT(1) = ( COUL(1) + COUL(2) + COUL(3) + COUL(4) ) * 0.25

      KM1 = 4
      DO K=1,4

C        LES XY ET COULEURS DU TRIANGLE B K-1 K
         DO I=1,2
            XYT(I,2) = XY(I,KM1)
            XYT(I,3) = XY(I,K  )
         ENDDO
         COULT(2) = COUL(KM1)
         COULT(3) = COUL(K  )

C        TRACE DU TRIANGLE
         CALL TRIACOUL2DBORD( XYT, COULT, NC )

C        PASSAGE AU TRIANGLE SUIVANT
         KM1 = K

      ENDDO

C     SI QUADRANGLE CROISE IL FAUT TRACER B13 ET B24
      DO K=1,2

C        LES XYZ ET COULEURS DU TRIANGLE B K K+2
         DO I=1,2
            XYT(I,2) = XY(I,K)
            XYT(I,3) = XY(I,K+2)
         ENDDO
         COULT(2) = COUL(K)
         COULT(3) = COUL(K+2)

C        TRACE DU TRIANGLE
         CALL TRIACOUL2DBORD( XYT, COULT, NC )

      ENDDO

      RETURN
      END
