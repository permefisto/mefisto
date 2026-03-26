      SUBROUTINE TRAIT23D( NC, XYZ1, XYZ2 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : TRACE LE TRAIT en 2D ou 3D DE (XYZ1) A (XYZ2) EN COULEUR NC
C -----
C       LA TRANSFORMATION EN PIXELS EST ASSUREE DANS CE SP
C       LES AUTRES CARACTERISTIQUES DU TRACE SONT CELLES ACTUELLES
C       NDIMLI est UNE VARIABLE de trvari.inc 
C
C ENTREES:
C --------
C NC   : NUMERO DE LA COULEUR DU TRAIT A TRACER
C XYZ1 : 1-ERE EXTREMITE DU TRAIT
C XYZ2 : 2-EME EXTREMITE DU TRAIT
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : Alain PERRONNET  Saint PIERRE du PERRAY            Avril 2020
C2345X7..............................................................012
      include"./incl/trvari.inc"
      include"./incl/minint.inc"
      REAL   XYZ1(3), XYZ2(3)

C     SI COULEUR NEGATIVE PAS DE TRACE
      IF( NC .LT. 0 ) RETURN

      IF( NDIMLI .EQ. 2 ) THEN

C        TRACE EN 2D
         CALL TRAIT2D( NC, XYZ1(1),XYZ1(2), XYZ2(1),XYZ2(2) )
         GOTO 9999

      ENDIF

      IF( NDIMLI .EQ. 3 ) THEN

C        TRACE EN 3D
         CALL  TRAIT3D( NC, XYZ1, XYZ2 )

      ENDIF

 9999 RETURN
      END
