      SUBROUTINE BARYFACE( NOFACE, NOSOEF, XYZSOM, BARY )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCUL DES 3 COORDONNEES DU BARYCENTRE DE LA FACE NOFACE DE NOSOEF
C -----
C ENTREES:
C --------
C NOFACE : NUMERO NOSOEF DE LA FACE A AJOUTER ET TRACER
C NOSOEF : NUMERO XYZSOM DES 4 SOMMETS DE LA FACE
C          SI LE 4-EME EST NUL ALORS LA FACE EST UN TRIANGLE
C          SINON UN QUADRANGLE
C XYZSOM : XYZ DES SOMMETS DU MAILLAGE

C SORTIE :
C --------
C BARY   : XYZ DU BARYCENTRE DE LA FACE NOFACE DU TABLEAU NOSOEF
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  Saint Pierre du Perray            Avril 2020
C.......................................................................
      INTEGER  NOSOEF(4,*)
      REAL     XYZSOM(3,*), BARY(3)

C     LE BARYCENTRE DE LA FACE NOFACE DE NOSOEF
      DO K=1,3
         BARY( K ) = 0
      ENDDO

      NBS = 0
      DO N=1,4
C        LE NUMERO DU SOMMET N DE LA FACE NOFACE
         NS = NOSOEF( N, NOFACE )
         IF( NS .GT. 0 ) THEN
            NBS = NBS + 1
            DO K=1,3
               BARY( K ) = BARY( K ) + XYZSOM( K, NS )
            ENDDO
         ENDIF
      ENDDO

      DO K=1,3
         BARY( K ) =  BARY( K ) / NBS
      ENDDO

      RETURN
      END
