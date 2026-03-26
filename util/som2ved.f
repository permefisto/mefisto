      SUBROUTINE SOM2VED( NBVAR, B01, B02, B )
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! BUT :    B <= B01 + B02
! -----    B , B01 , B02 VECTEURS DE NBVAR VARIABLES DOUBLE PRECISION
!          VARIANTE OpenMP f95
!
! ENTREES:
! --------
! NBVAR  : NOMBRE DE VARIABLES DES VECTEURS B, B01, B02
! B01    : PREMIER VECTEUR REEL DOUBLE PRECISION DE LA SOMME
! B02    : SECOND  VECTEUR REEL DOUBLE PRECISION DE LA SOMME
!
! SORTIE :
! --------
! B      : VECTEUR SOMME DES 2 VECTEURS
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! AUTEUR : PERRONNET ALAIN LJLL UPMC & St Pierre du Perray Novembre 2012
!2345X7..............................................................012
      IMPLICIT  NONE
      DOUBLE PRECISION  B01(NBVAR), B02(NBVAR), B(NBVAR)
      INTEGER           NBVAR, K

      DO K = 1, NBVAR
         B(K) = B01(K) + B02(K)
      ENDDO


      RETURN
      END SUBROUTINE SOM2VED
