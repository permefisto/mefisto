      SUBROUTINE CHOIXFONTE( LHPX )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CHOISIR ET CHARGER UNE FONTE CORRECTE POUR MEFISTO
C -----    DE HAUTEUR MAXIMALE LHPX PIXELS
C
C ENTREE:
C --------
C LHPX  : LA HAUTEUR LA PLUS PROCHE EN PIXELS DE LA FONTE A CHARGER
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN LJLL UPMC & St PIERRE du PERRAY  Juillet 2012
C2345X7..............................................................012
      include"./incl/xvfontes.inc"
C
C     CHARGEMENT DE LA FONTE DE HAUTEUR LHPX
      IF( LHPX .LT. 0 .OR. LHPX .GT. MAXHPX ) THEN
         NUF = NUFOHPX( 7 )
      ELSE
         NUF = NUFOHPX( LHPX )
      ENDIF
C
      CALL CHARGEFONTE( NUF )
C
      RETURN
      END
