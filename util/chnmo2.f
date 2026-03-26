      SUBROUTINE CHNMO2( NOMOPE, NCODE )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : RECHERCHE L'OPERATEUR NOMOPE DANS LA LISTE DES OPERATEURS BINAIRES
C ----- CE PEUT NE PAS ETRE UN TEL OPERATEUR
C
C ENTREE :
C --------
C NOMOPE : NOM DE L'OPERATEUR EVENTUEL (6 CARACTERES AU PLUS)
C
C SORTIE :
C --------
C NCODE  : =0 OPERATEUR BINAIRE NON RETROUVE
C          >0 NUMERO DE L'OPERATEUR
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS   DECEMBRE 1993
C23456---------------------------------------------------------------012
      CHARACTER*6    NOMOPE
      include"./incl/kmopbi.inc"
C
C     IDENTIFICATION DE L'OPERATEUR
      DO 10 I=N1OPBI,NBOPBI
         IF( NOMOPE .EQ. NMOPBI(I) ) THEN
            IF( I .GT. 0 ) THEN
C              OPERATEUR BINAIRE DIFFERENT DE 'OR ' 'AND '
               NCODE = I
            ELSE
C              EQUIVALENCE 'OR ' 'AND ' ET 'OU ' 'ET '
               NCODE = 2 - N1OPBI + I
            ENDIF
            RETURN
         ENDIF
 10   CONTINUE
C
C     NON RETROUVE
      NCODE = 0
      RETURN
      END
