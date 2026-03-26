      SUBROUTINE CHITEM( NL , NC , NLD , NCD , NLF , NCF )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : RETOURNER (NLD,NCD) ET (NLF,NCF) DEBUT ET FIN DANS KTD
C ----- DU 1-ER CARACTERE ET DERNIER CARACTERE D'UN ITEM TERMINE PAR
C       PAR UN BLANC OU LA FIN DE LIGNE
C
C ENTREES :
C ---------
C NL , NC : POSITION DANS KTD DU CARACTERE QUI PRECEDE
C           LE PREMIER CARACTERE A ANALYSER
C           RESTENT INCHANGES EN SORTIE
C
C SORTIES :
C ---------
C NLD,NCD : POSITION DU PREMIER CARACTERE DE L'ITEM
C NLF,NCF : POSITION DU DERNIER CARACTERE DE L'ITEM
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS SEPTEMBRE 1988
C23456---------------------------------------------------------------012
C.......................................................................
      include"./incl/td.inc"
C.......................................................................
C     LE 1ER CARACTERE NON BLANC
      NLD = NL
      NCD = NC
      CALL CAR1NB( NLD , NCD )
C
C     RECHERCHE DU PREMIER BLANC
      NCF = INDEX( KTD(NLD)(NCD:NCKTD) , ' ' )
      IF( NCF .LE. 0 ) THEN
         NCF = NCKTD
      ELSE
         NCF = NCD + NCF - 2
      ENDIF
      NLF = NLD
      END
