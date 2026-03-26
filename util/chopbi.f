      SUBROUTINE CHOPBI( NLD , NCD , NCDF , NLF , NCF , NCODE )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : LE MOT KLG(NLD,NCD:NCDF) EST IL UN OPERATEUR BINAIRE ?
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
C NCODE   : 0  OPERATEUR BINAIRE NON RETROUVE
C           >0 NUMERO DE L'OPERATEUR
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS    OCTOBRE 1989
C23456---------------------------------------------------------------012
      include"./incl/lu.inc"
      include"./incl/kmopbi.inc"
C
C     IDENTIFICATION DE L'OPERATEUR
      DO 100 I=N1OPBI,NBOPBI
C
C        LE NOMBRE DE CARACTERES DE L'OPERATEUR BINAIRE
         NBC = NCOPBI(I)
C
C        LA FIN DU MOT EST ELLE CONNUE?
         IF( NCDF .GT. 0 ) THEN
            NCF = NCDF
         ELSE
            NCF = NCD - 1 + NBC
         ENDIF
C
         IF( KLG(NLD)(NCD:NCF) .EQ. NMOPBI(I)(1:NBC) ) THEN
C           OPERATEUR RETROUVE
            IF( I .GT. 14 ) THEN
C
C              L'OPERATEUR N'EST PAS:        'OR ' , 'AND ' ,
C                                    'OX ' , 'OU ' , 'ET '  , '<'  , '<=' ,
C                                    '='   , '<>'  , '>='   , '>'  , '+'  ,
C                                    '-'   , '*'   , '/'    , '**'
C              EST IL SUIVI D'UN CARACTERE ( ?
               NL1 = NLD
               NC1 = NCF
               CALL CARPNB( NL1, NC1 )
               IF( KLG(NL1)(NC1:NC1) .NE. '(' ) GOTO 100
            ENDIF
            IF( I .GT. 0 ) THEN
C              OPERATEUR BINAIRE DIFFERENT DE 'OR ' ET 'AND '
               NCODE = 200 + I
            ELSE
C              EQUIVALENCE 'OR ' ET 'AND ' AVEC 'OU ' ET 'ET '
               NCODE = 202 - N1OPBI + I
            ENDIF
C
            NLF = NLD
            RETURN
         ENDIF
 100  CONTINUE
C
C     OPERATEUR NON RETROUVE
      NCODE = 0
      NLF   = 0
C
      RETURN
      END
