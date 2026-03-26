      SUBROUTINE BARYEF( NBSTEF, NOSOEL, XYZSOM, BARY )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCUL DES 3 COORDONNEES DU BARYCENTRE DE L'EF de SOMMETS NOSOEL
C -----
C ENTREES:
C --------
C NBSTEF : NOMBRE DE SOMMETS DE L'EF
C NOSOEL : NUMERO XYZSOM DES NBSTEF SOMMETS DE LA FACE
C XYZSOM : XYZ DES SOMMETS DU MAILLAGE

C SORTIE :
C --------
C BARY   : XYZ DU BARYCENTRE DE L'EF DE SOMMETS NOSOEL
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  Saint Pierre du Perray            Avril 2020
C.......................................................................
      INTEGER  NOSOEL(NBSTEF)
      REAL     XYZSOM(3,*), BARY(3)

C     LE BARYCENTRE DE L'EF DE SOMMETS NOSOEL
      DO K=1,3
         BARY( K ) = 0
      ENDDO

      NBS = 0
      DO N = 1, NBSTEF

C        LE NUMERO DU SOMMET N DE L'EF
         NS = NOSOEL( N )
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
