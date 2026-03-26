      SUBROUTINE CHETIQ( NLD , NCD , NLF , NCF , NOETIQ )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : RECHERCHE A PARTIR DU CARACTERE (NLD,NCD) D'UNE ETIQUETTE
C ----- CE PEUT NE PAS ETRE UNE ETIQUETTE
C
C ENTREES :
C ---------
C NLD,NCD : POSITION DANS KLG DU PREMIER CARACTERE A TRAITER
C
C SORTIES :
C ---------
C NLF,NCF : SI NOETIQ>0 POSITION DANS KLG DU DERNIER CARACTERE DE
C                       L'ETIQUETTE
C           SINON NLF=0
C NOETIQ :  0 PAS D'ETIQUETTE RETROUVEE
C          >0 LE NUMERO DE L'ETIQUETTE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS    FEVRIER 1990
C23456---------------------------------------------------------------012
      include"./incl/lu.inc"
C
      DO 10 NOETIQ=NBETIQ,1,-1
         NBC = INDEX( KETIQ(NOETIQ) , ' ' ) - 1
         IF( NBC .LE. 0 ) NBC = NBCAET
         NBC = NBC - 1
         IF( KLG(NLD)(NCD:NCD+NBC) .EQ. KETIQ(NOETIQ)(1:1+NBC) ) THEN
C           L'ETIQUETTE NOETIQ EST RETROUVEE
            NLF = NLD
            NCF = NCD + NBC
            RETURN
         ENDIF
 10   CONTINUE
C
C     ETIQUETTE NON RETROUVEE
      NOETIQ = 0
      NLF    = 0
      END
