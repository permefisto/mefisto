      SUBROUTINE TRAXES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : TRACER DES AXES DU REPERE ORTHONORME SELON LE COMMON TRVARI
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC  PARIS       MARS 1994
C2345X7..............................................................012
      include"./incl/trvari.inc"
C
C     LE TRACE DES AXES
C     =================
      IF( NDIMLI .LE. 2 ) THEN
         CALL TRAXE2
      ELSE
         CALL TRAXE3
      ENDIF
      END
