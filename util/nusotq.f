      INTEGER FUNCTION NUSOTQ( NBS1, NBS2, NBS4, I, J )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    NUMERO DU SOMMET (I,J) D'UN QUADRANGLE TRIANGULE-QUADRANGULE
C -----
C ENTREES :
C ---------
C NBS1   : NOMBRE DE SOMMETS DU COTE 1 ET 3 DU QUADRANGLE
C NBS2   : NOMBRE DE SOMMETS DU COTE 2 DU QUADRANGLE   NBS2 >= NBS4
C NBS4   : NOMBRE DE SOMMETS DU COTE 4 DU QUADRANGLE
C I      : INDICE SELON LE COTE 1 ET 3 (1<=I<=NBS1=NBS3)
C J      : INDICE SELON LE COTE 2 ET 4 (1<=J<=NBS2)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS     JUILLET 1993
C2345X7..............................................................012
      IF( NBS2 .EQ. NBS4 ) THEN
C
C        QUADRANGULATION PURE
         NUSOTQ = I + (J-1) * NBS1
C
      ELSE
C
C        TRIANGULATION-QUADRANGULATION
         IF( J .LE. NBS4 ) THEN
C
C           POINT DANS LA PARTIE INFERIEURE REGULIERE
            NUSOTQ = I + (J-1) * NBS1
C
         ELSE
C
C           POINT DANS LA PARTIE IRREGULIERE
            NUSOTQ = NBS1 * NBS4
            DO 10 JJ=NBS4+1,J
               NUSOTQ = NUSOTQ + NBS2 - JJ + 1
 10         CONTINUE
            NUSOTQ = NUSOTQ - NBS1 + I
         ENDIF
      ENDIF
      END
