      SUBROUTINE CHINT2( NL  , NC  ,
     %                   NLD , NCD ,  NLF , NCF )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : RETOURNER LA VALEUR DE L'INDICE INITIAL ET FINAL
C ----- D'UN INTERVALLE_E DEFINI PAR
C         INTERVALLE_E := |INTERVAL_ENT   |ID_ENTIER
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
C           NLD = 0 SI INTERVALLE_E INCORRECT
C NLF,NCF : POSITION DANS KTD DU DERNIER CARACTERE DE INTERVAL_ENT
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS SEPTEMBRE 1988
C23456---------------------------------------------------------------012
      include"./incl/td.inc"
      include"./incl/gsmenu.inc"
C
C     RECHERCHE DE LA PREMIERE VALEUR ENTIERE DE L'INTERVALLE
      NLD = NL
      NCD = NC
      CALL VAENT2( NL , NC , NLD , NCD , NLF , NCF , N1 , NR )
      IF( NLD .LE. 0 ) RETURN
      NL0 = NLF
      NC0 = NCF
C
C     RECHERCHE DE ..
      CALL CARAPR( NL0 , NC0 )
      IF( KTD(NL0)(NC0:NC0) .EQ. '.' ) THEN
C        LE CARACTERE SUIVANT EST IL UN AUTRE . ?
         CALL CARAPR( NL0 , NC0 )
         IF( KTD(NL0)(NC0:NC0) .NE. '.' ) THEN
C           ERREUR . ET NON PAS ..
            NBLGRC(NRERR) = 2
            KERR(1) = 'CHINT2:INTERVALLE .. INCORRECT DANS'
            KERR(2) = KTD(NL0)
            CALL LEREUR
            NLD = 0
            RETURN
         ENDIF
C
C        .. EST RETROUVE
C        CALCUL DE LA VALEUR FINALE DE L'INTERVALLE
         CALL VAENT2( NL0 , NC0 , I , J , NLF , NCF , N1 , NR )
         IF( I .LE. 0 ) THEN
            NLD = 0
            RETURN
         ENDIF
      ENDIF
      END
