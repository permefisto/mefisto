      SUBROUTINE REORDD( NIVO, NBV, NBCOMP, VAL, VECT )
C **********************************************************************
C BUT :    ON REORDONNE PAR ORDRE CROISSANT LES VALEURS PROPRES
C -----    ET EVENTUELLEMENT LES VECTEURS PROPRES
C
C ENTREES:
C --------
C NIVO   : =0  LES VALEURS PROPRES SEULES      SONT ORDONNEES
C          =1  LES VALEURS ET VECTEURS PROPRES SONT ORDONNEES
C NBV    : NOMBRE DE VALEURS ET VECTEURS PROPRES
C NBCOMP : NOMBRE DE COMPOSANTES DES VECTEURS PROPRES
C
C MODIFIES:
C ---------
C VAL    : LES NBV VALEURS  PROPRES
C VECT   : LES NBV VECTEURS PROPRES
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LABO ANALYSE NUMERIQUE UPMC PARIS   AOUT 1998
C23456---------------------------------------------------------------012
      DOUBLE PRECISION VAL(NBV), VECT(NBV,NBCOMP), AUX
C
      NBV1 = NBV - 1
      IF( NIVO .EQ. 0 ) THEN
C
C        SEULES LES VALEURS PROPRES SONT REORDONNEES
         DO 10 I=1,NBV1
            I1=I+1
            DO 5 J=I1,NBV
               IF( VAL(J) .LT. VAL(I) ) THEN
                  AUX    = VAL(J)
                  VAL(J) = VAL(I)
                  VAL(I) = AUX
               ENDIF
   5        CONTINUE
  10     CONTINUE
C
      ELSE
C
C        LES VALEURS ET VECTEURS PROPRES SONT REORDONNES
         DO 40 I=1,NBV1
            I1=I+1
            DO 30 J=I1,NBV
               IF( VAL(J) .LT. VAL(I) ) THEN
                  AUX    = VAL(J)
                  VAL(J) = VAL(I)
                  VAL(I) = AUX
                  DO 25 K=1,NBCOMP
                     AUX       = VECT(J,K)
                     VECT(J,K) = VECT(I,K)
                     VECT(I,K) = AUX
   25             CONTINUE
               ENDIF
   30       CONTINUE
   40    CONTINUE
      ENDIF
      END
