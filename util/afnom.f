      SUBROUTINE AFNOM( NLD , NCD , NLF , NCF )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : AFFICHER LE NOM DEBUTANT EN (NLD,NCD) ET FINISSANT EN (NLF,NCF)
C ----- DU TABLEAU KTD
C
C ENTREES :
C ---------
C NLD, NCD: POSITION DANS KTD DU PREMIER CARACTERE DU NOM
C NLF, NCF: POSITION DANS KTD DU PREMIER CARACTERE DU NOM
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS   SEPTEMBRE 1988
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/td.inc"
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
C
C     LA LIGNE EST BUFFERISEE A L'AIDE DU POINTEUR LCLIGN DERNIER
C     CARACTERE ENTRE DANS LA LIGNE
C
      IF( IMPRES*NOMUET .LE. 0 ) RETURN
C
      IF( NLD .NE. NLF ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'AFNOM: NOM SUR PLUSIEURS LIGNES INCORRECT'
         ELSE
            KERR(1) = 'AFNOM: INCORRECT NAME on SEVERAL LINES'
         ENDIF
         CALL LEREUR
         RETURN
      ENDIF
C
      IF( NCD .GT. NCF ) THEN
         NBLGRC(NRERR) = 2
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'AFNOM: NOM INCORRECT dans'
         ELSE
            KERR(1) = 'AFNOM: INCORRECT NAME in'
         ENDIF
         KERR(2) = KTD(NLD)
         CALL LEREUR
         RETURN
      ENDIF
C
C     LE NOMBRE DE CARACTERES DU NOM
      NB = NCF - NCD + 1
C
C     Y A-T-IL SUFFISAMMENT DE PLACE DANS KLIGNE ?
      IF( NCLIGN - LCLIGN .LT. 1 + NB ) THEN
C        NON
         CALL AFLIGN
      ENDIF
C
C     STOCKAGE DU NOM
      KLIGNE(LCLIGN+1:NCLIGN) = ' ' // KTD(NLD)(NCD:NCF)
      LCLIGN = LCLIGN + 1 + NB
C
      RETURN
      END
