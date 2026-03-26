      SUBROUTINE ENTMAX( A , B , MAXIMA )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     MAXIMA  INDIQUE LEQUEL DES 2 ENTIERS A ET B MULTI-MOTS
C -----     A LA PLUS GRANDE VALEUR ABSOLUE
C
C ENTREES :
C ---------
C A       : ENTIER MULTI-MOTS  CODE SOUS LA FORME
C                I=NA
C           A = SOMME ( A(I) * BASE ** I )   ET  A(-1) = SIGNE( A )
C                I=0                             A(-2) = NA
C
C B       : ENTIER MULTI-MOTS  CODE SOUS LA FORME
C                I=NB
C           B = SOMME ( B(I) * BASE ** I )   ET  B(-1) = SIGNE( B )
C                I=0                             B(-2) = NB
C
C SORTIES :
C ---------
C MAXIMA  : 1 SI A > B
C           0 SI A = B
C          -1 SI A < B
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS   FEVRIER 1988
C23456...............................................................012
      INTEGER A(-2:*),B(-2:*)
C
C     TEST SUR LE NOMBRE DE MOTS DES ENTIERS
      IF( A(-2) .GT. B(-2) ) THEN
C        A > B
         MAXIMA = 1
      ELSE IF( A(-2) .LT. B(-2) ) THEN
C        A < B
         MAXIMA = -1
      ELSE
C        MEME NOMBRE DE MOTS
         DO 10 I=A(-2),0,-1
C           COMPARAISON MOT A MOT
            IF( A(I) .GT. B(I) ) THEN
C              A > B
               MAXIMA = 1
               RETURN
            ELSE IF( A(I) .LT. B(I) ) THEN
C              A < B
               MAXIMA = -1
               RETURN
            ENDIF
 10      CONTINUE
C        A = B
         MAXIMA = 0
      ENDIF
      END
