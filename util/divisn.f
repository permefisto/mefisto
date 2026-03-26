      REAL FUNCTION DIVISN( Q , NR , NT )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCULER LA DIVISION D'UN INTERVALLE
C
C ENTREES :
C --------
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : BATOUCHE_EMMEL MARS DEA D'A.N. MAI 1989
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF (Q.LE.0) THEN
         STOP 'ERREUR DANS DIVISN : RAISON GEOMETRIQUE INCORRECT !'
      ENDIF
C
      IF ( Q .EQ. 1.0 ) THEN
         DIVISN = REAL(NR-1)/REAL(NT-1)
      GOTO 100
      ENDIF
C
      DIVISN = (1-Q**(NR-1))/(1-Q**(NT-1))
C
 100  RETURN
      END
