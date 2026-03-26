      SUBROUTINE SO3C6C( NBSO3C, NB3C, NS3C6C )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    FOURNIR LE NUMERO DES SOMMETS DES 3CUBES D'UN 6-CUBE
C -----    POUR X4 X5 X6 FIXES AUX VALEURS 0 ET 1 => 8 premiers
C          POUR X1 X2 X3 FIXES AUX VALEURS 0 ET 1 => 8 derniers
C
C SORTIES:
C --------
C NBSO3C :  8=NOMBRE DE SOMMETS D'UN 3-CUBE D'UN 6-CUBE
C NB3C   : 16=NOMBRE DE 3CUBES 8 pour X4X5X6 FIXES D'UN 6-CUBE
C                              8 pour X1X2X3 FIXES D'UN 6-CUBE
C NS3C6C : (I,J) NO DU I-EME SOMMET DU 3-CUBE D'UN 6-CUBE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: : ALAIN PERRONNET  LJLL UPMC PARIS                OCTOBRE 2005
C2345X7..............................................................012
      INTEGER  NS3C6C(8,16)
      include"./incl/nusc3c6.inc"
C
      NBSO3C =  8
      NB3C   = 16
      DO 20 J=1,NB3C
         DO 10 I=1,NBSO3C
            NS3C6C(I,J) = NUS3C6C(I,J)
 10      CONTINUE
 20   CONTINUE
C
      RETURN
      END
