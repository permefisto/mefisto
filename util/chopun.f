      SUBROUTINE CHOPUN( NLD , NCD , NCDF , NLF , NCF , NCODE )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : RECHERCHE A PARTIR DU CARACTERE (NLD,NCD) D'UN OPERATEUR UNAIRE
C ----- CE PEUT NE PAS ETRE UN TEL OPERATEUR
C
C ENTREES :
C ---------
C NLD,NCD : POSITION DANS KLG DU PREMIER CARACTERE A TRAITER
C NCDF    : SI >0 POSITION DU DERNIER CARACTERE DU NOM A COMPARER
C              =0 POSITION INCONNUE
C
C SORTIES :
C ---------
C NLF,NCF : SI NCODE>0 OPERATEUR RETROUVE ENTRE NLD,NCD ET NLF,NCF
C           SINON NLF=0
C NCODE   : 0  OPERATEUR UNAIRE NON RETROUVE
C           >0 NUMERO DE L'OPERATEUR
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS  OCTOBRE 1989
C23456---------------------------------------------------------------012
      include"./incl/lu.inc"
      include"./incl/kmopun.inc"
C
C     IDENTIFICATION DE L'OPERATEUR
      DO 100 I=N1OPUN,NBOPUN
C
C        LE NOMBRE DE CARACTERES DE L'OPERATEUR UNAIRE
         NBC = NCOPUN(I)
C
C        LA FIN DU MOT EST ELLE CONNUE?
         IF( NCDF .GT. 0 ) THEN
            NCF = NCDF
         ELSE
            NCF = NCD - 1 + NBC
         ENDIF
C
         IF( KLG(NLD)(NCD:NCF) .EQ. NMOPUN(I)(1:NBC) ) THEN
C           OPERATEUR RETROUVE
            IF( I .GT. 3 ) THEN
C              L'OPERATEUR N'EST PAS : 'NOT' , 'NON ' , '+'    ,'-'
C              EST IL SUIVI D'UN CARACTERE ( ?
               NL1 = NLD
               NC1 = NCF
               CALL CARPNB( NL1, NC1 )
               IF( KLG(NL1)(NC1:NC1) .NE. '(' ) GOTO 100
            ENDIF
            IF( I .GT. 0 ) THEN
C              OPERATEUR UNAIRE DIFFERENT DE 'NOT'
               NCODE = 100 + I
            ELSE
C              EQUIVALENCE DU 'NOT' AVEC 'NON'
               NCODE = 101 - N1OPUN + I
            ENDIF
            NLF = NLD
            RETURN
         ENDIF
 100  CONTINUE
C
C     OPERATEUR NON RETROUVE
      NCODE = 0
      NLF   = 0
      RETURN
      END
