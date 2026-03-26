      SUBROUTINE CHINEN( NL  , NC  ,
     %                   NLD , NCD ,  NLF , NCF , N1 , N2 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : RETOURNER LA VALEUR DE L'INDICE INITIAL ET FINAL
C ----- D'UN INTERVAL_ENT  DEFINI PAR
C         INTERVAL_ENT := ID_ENTIER .. ID_ENTIER
C         ID_ENTIER    := |ENTIER  |IDENT
C
C ENTREES :
C ---------
C NL , NC : POSITION DANS KTD DU CARACTERE QUI PRECEDE
C           LE PREMIER CARACTERE A ANALYSER
C           RESTENT INCHANGES EN SORTIE
C
C SORTIES :
C ---------
C NLD,NCD : POSITION DU PREMIER CARACTERE DE L'INTERVALLE
C NLF,NCF : POSITION DANS KTD DU DERNIER CARACTERE DE INTERVAL_ENT
C N1 N2   : VALEURS INITIALE ET FINALE DE L'INTERVALLE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS SEPTEMBRE 1988
C23456---------------------------------------------------------------012
C.......................................................................
      include"./incl/td.inc"
C.......................................................................
C     RECHERCHE DE LA PREMIERE VALEUR ENTIERE DE L'INTERVALLE
      NLD = NL
      NCD = NC
      CALL VAENTI( NL , NC , NLD , NCD , NLF , NCF , N1 )
C
C     RECHERCHE DE ..
 10   CALL CARAPR( NLF , NCF )
      IF( KTD(NLF)(NCF:NCF) .NE. '.' ) THEN
         GOTO 10
      ELSE
C        LE CARACTERE SUIVANT EST IL UN AUTRE . ?
         CALL CARAPR( NLF , NCF )
         IF( KTD(NLF)(NCF:NCF) .NE. '.' ) THEN
C           NON . ON PASSE AU CARACTERE SUIVANT
            GOTO 10
         ELSE
C           OUI . CALCUL DE LA VALEUR FINALE DE L'INTERVALLE
            CALL VAENTI( NLF , NCF , I , J , NLF , NCF , N2 )
         ENDIF
      ENDIF
      END
