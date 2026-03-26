      SUBROUTINE LXNLOU( NTLX, NL, NOTAMS, MCTAMS )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : RETOURNER LE NUMERO NOTAMS D'UN TMS ET SON ADRESSE MCN MCTAMS
C ----- DU TMS NL DU LEXIQUE NTLX
C
C ENTREES :
C ---------
C NTLX   : NUMERO TMS DU LEXIQUE
C NL     : NUMERO LOCAL AU LEXIQUE DU TMS A OUVRIR
C
C SORTIES :
C ---------
C NOTAMS : NUMERO DU TMS
C MCTAMS : ADRESSE MCN DU TMS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     OCTOBRE 1988
C23456---------------------------------------------------------------012
      include"./incl/gsmenu.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
C
C     OUVERTURE DU LEXIQUE
      CALL TAMSOU( NTLX , MNLX )
C
      IF( MNLX .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'LXNLOU: LEXIQUE NON OUVRABLE'
         CALL LEREUR
         NOTAMS = 0
         MCTAMS = 0
         RETURN
      ENDIF
C
C     LE NUMERO DE TMS DU NL-EME NOM DANS LE LEXIQUE
      NOTAMS = MCN( MNLX + MCN(MNLX) * NL + MCN(MNLX+2) + 2 )
C
      IF( NOTAMS .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'LXNLOU: TMS du NL-EME NOM dans le LEXIQUE <=0'
         CALL LEREUR
         NOTAMS = 0
         MCTAMS = 0
         RETURN
      ENDIF
C
C     OUVERTURE DU TABLEAU NOTAMS
      CALL TAMSOU( NOTAMS , MCTAMS )
C
      RETURN
      END
