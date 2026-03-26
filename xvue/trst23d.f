      SUBROUTINE TRST23D( NCOUL, NOST, XYZSOM )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    en 2D ou 3D TRACE des CARACTERES '+NOST' DU SOMMET NOST
C -----    de XYZSOM le TABLEAU des XYZ des SOMMETS du MAILLAGE

C ENTREES:
C --------
C NCOUL  : NUMERO DE LA COULEUR DE TRACE
C NOST   : NUMERO XYZSOM du SOMMET NOST
C XYZSOM : XYZ DES SOMMETS DU MAILLAGE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : Alain PERRONNET  Saint PIERRE du PERRAY            Avril 2020
C2345X7..............................................................012
      include"./incl/trvari.inc"
      REAL           XYZSOM(3,*)

      IF( NDIMLI .EQ. 2 ) THEN
         CALL TRST2D( NCOUL, NOST, XYZSOM )
         GOTO 9999
      ENDIF

      IF( NDIMLI .EQ. 3 ) THEN
         CALL TRST3D( NCOUL, NOST, XYZSOM )
      ENDIF

 9999 RETURN
      END
