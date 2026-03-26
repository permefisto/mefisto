      SUBROUTINE AFCAR( CAR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     AFFICHER LES NBCAR CARACTERES DE CAR PRECEDE D'UN BLANC
C -----
C
C ENTREES :
C ---------
C NBCAR   : NOMBRE DE CARACTERES
C CAR     : LES CARACTERES A AFFICHER
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS   SEPTEMBRE 1988
C23456---------------------------------------------------------------012
      include"./incl/td.inc"
C.......................................................................
      CHARACTER*(*)      CAR
C
C     LA LIGNE EST BUFFERISEE A L'AIDE DU POINTEUR LCLIGN DERNIER
C     CARACTERE ENTRE DANS LA LIGNE
C
      IF( IMPRES * NOMUET .LE. 0 ) RETURN
C
C     LE NOMBRE DE CARACTERES DE LA CHAINE CAR
      NBCAR = LEN( CAR )
C
C     LE RESTE DE KLIGNE PEUT IL CONTENIR CAR ?
      IF( NCLIGN-LCLIGN .LT. NBCAR ) THEN
         CALL AFLIGN
      ENDIF
C
      KLIGNE(LCLIGN+1:LCLIGN+1+NBCAR) = ' ' // CAR(1:NBCAR)
      LCLIGN = LCLIGN + 1 + NBCAR
      END
