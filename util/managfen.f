      SUBROUTINE MANAGFEN
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : MANAGEMENT GESTION de la FENETRE GRAPHIQUE de MEFISTO
C -----
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : Alain PERRONNET  Saint PIERRE du PERRAY            AVRIL 2020
C2345X7..............................................................012
      include"./incl/mecoit.inc"
      include"./incl/trvari.inc"

C     LECTURE DU MOT-CLE OU OPTION A EXECUTER PARMI TOUTES LES OPTIONS
C     CONTENUES DANS LE MENU DE NOM 'managfen'
C     ----------------------------------------------------------------
 1    CALL LIMTCL( 'managfen', NMTCL )

C     TRAITEMENT DE L'OPTION NMTCL
      IF( NMTCL .LT.  0 ) GOTO 9999

      GOTO( 10, 20, 9999 ), NMTCL

C     NOMBRE PIXELS de la LARGEUR HAUTEUR de la FENETRE
 10   CALL XVPXFE
      GOTO 1

C     COULEUR du FOND de la FENETRE GRAPHIQUE Mefisto
 20   CALL INVITE( 17 )
      CALL LIMTCL( 'couleurs' , I )
      IF( I .LT. 0 ) THEN
         GOTO 1
      ELSE IF( I .EQ. 0 ) THEN
C        LA COULEUR NOIRE
         NCOFON = 0
      ELSE
         NCOFON = N1COEL + I
      ENDIF
C     LA COULEUR DU FOND EST IMPOSEE
      CALL XVFOND( NCOFON )
      CALL EFFACE
      GOTO 1

C     Sortie
 9999 RETURN
      END


