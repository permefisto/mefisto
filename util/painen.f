      SUBROUTINE PAINEN( NL  , NC  ,
     %                   NLD , NCD ,  NLF , NCF , N1 , N2 , NRETOU )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : RETOURNER SOUS CONDITION LA VALEUR DE L'INDICE INITIAL ET FINAL
C ----- D'UN INTERVAL_ENT  DEFINI PAR
C         INTERVAL_ENT := ID_ENTIER .. ID_ENTIER
C         ID_ENTIER    := |ENTIER  |IDENT
C
C       SI ID_ENTIER = ENTIER    SEULEMENT
C                      ------    ---------
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
C           IINFO( 'GRAND' ) SI UN DES ID_ENTIER EST ERRONE
C NRETOU  : 0 SI N1 N2 SONT CONNUS   C-A-D  ID_ENTIER = ENTIER
C           1 SI N1 N2 SONT INCONNUS C-A-D  ID_ENTIER = IDENT
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS   OCTOBRE 1988
C23456---------------------------------------------------------------012
C.......................................................................
      include"./incl/td.inc"
C.......................................................................
C     RECHERCHE DE LA PREMIERE VALEUR ENTIERE DE L'INTERVALLE
      NLD = NL
      NCD = NC
      CALL VAENT2( NL , NC , NLD , NCD , NLF , NCF , N1 , NRETOU )
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
            CALL VAENT2( NLF , NCF , I , J , NLF , NCF , N2 , NRETO1 )
            NRETOU = NRETOU + NRETO1
         ENDIF
      ENDIF
      END
