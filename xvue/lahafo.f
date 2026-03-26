      SUBROUTINE LAHAFO
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CHOIX DE LA FONTE COURANTE TEL QUE LE MENU TRACMAIL TIENNE
C ----- DANS LA FENETRE   NPHFCO est la HAUTEUR EN PIXELS D'UN CARACTERE
C       FONCTION de LHPXFE la HAUTEUR en PIXELS de la FENETRE
C       APPEL dans $MEFISTO/xvue/xvinit.f et luou.f
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  TIMS NTU TAIPEI TAIWAN          Octobre 2009
C2345X7..............................................................012
      include"./incl/trvari.inc"
      include"./incl/xvfontes.inc"

      IF( LHPXFE .LE. 600 ) THEN
         NPHFCO = 8
      ELSE IF( LHPXFE .LE. 700 ) THEN
         NPHFCO = 9
      ELSE IF( LHPXFE .LE. 800 ) THEN
         NPHFCO = 10
      ELSE IF( LHPXFE .LE. 900 ) THEN
         NPHFCO = 11
      ELSE
         NPHFCO = 12
      ENDIF

      RETURN
      END
