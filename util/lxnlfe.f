      SUBROUTINE LXNLFE( NTLX , NL )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    FERMER LE NOM DE NUMERO NL DU LEXIQUE NTLX
C -----
C
C ENTREES :
C ---------
C NTLX   : NUMERO TMS DU LEXIQUE
C NL     : NUMERO LOCAL AU LEXIQUE DU TMS A FERMER
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS  OCTOBRE 1988
C23456---------------------------------------------------------------012
      include"./incl/gsmenu.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
C
C     FERMETURE DU LEXIQUE
      CALL TAMSOU( NTLX , MNLX )
C
      IF( MNLX .LE. 0 ) THEN
          NBLGRC(NRERR) = 1
          KERR(1) = 'LXNLFE:LEXIQUE NON OUVRABLE'
          CALL LEREUR
          RETURN
      ENDIF
C
C     LE NUMERO DE TMS DU NL-EME NOM DANS LE LEXIQUE
      NOTAMS = MCN( MNLX + MCN(MNLX) * NL + MCN(MNLX+2) + 2 )
C
C     FERMETURE DU TABLEAU NOTAMS
      CALL TAMSFE( NOTAMS )
      END
