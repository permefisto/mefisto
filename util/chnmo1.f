      SUBROUTINE CHNMO1( NOMOPE, NCODE )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : RECHERCHE L'OPERATEUR NOMOPE DANS LA LISTE DES OPERATEURS UNAIRES
C ----- CE PEUT NE PAS ETRE UN TEL OPERATEUR
C
C ENTREE :
C --------
C NOMOPE : NOM DE L'OPERATEUR EVENTUEL (5 CARACTERES AU PLUS)
C
C SORTIE :
C --------
C NCODE  : =0 OPERATEUR UNAIRE NON RETROUVE
C          >0 NUMERO DE L'OPERATEUR UNAIRE RETROUVE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS   DECEMBRE 1993
C23456---------------------------------------------------------------012
      CHARACTER*5    NOMOPE
      include"./incl/kmopun.inc"
C
C     IDENTIFICATION DE L'OPERATEUR
      DO 10 I=N1OPUN,NBOPUN
         IF( NOMOPE .EQ. NMOPUN(I) ) THEN
            IF( I .GT. 0 ) THEN
C              OPERATEUR UNAIRE DIFFERENT DE 'NOT'
               NCODE = I
            ELSE
C              EQUIVALENCE DU 'NOT' AVEC 'NON'
               NCODE = 1 - N1OPUN + I
            ENDIF
            RETURN
         ENDIF
 10   CONTINUE
C
C     NON RETROUVE
      NCODE = 0
      RETURN
      END
