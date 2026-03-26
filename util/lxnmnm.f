      SUBROUTINE LXNMNM( NTLX , KNOM0 , KNOM1 )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CHANGER LE NOM KNOM0 EN LE NOM KNOM1 DANS LE LEXIQUE NTLX
C -----
C
C ENTREES :
C ---------
C NTLX    : NUMERO DU TMS DU LEXIQUE
C KNOM0   : NOM ANCIEN A REMPLACER
C KNOM1   : NOM NOUVEAU REMPLACANT
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : ALAIN PERRONNET ANALYSE NUMERIQUE  PARIS  DECEMBRE 1985
C..............................................................................
      include"./incl/gsmenu.inc"
      CHARACTER*(*)     KNOM0,KNOM1
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
C
C     RECHERCHE DU NOM KNOM0 DANS LE LEXIQUE NTLX
      CALL LXNMNO( NTLX , KNOM0 , NONOM , MNLX )
C
C     SI LE NOM N'EST PAS RETROUVE DIAGNOSTIC ET RETOUR
      IF( NONOM .LE. 0 ) THEN
          NBLGRC(NRERR) = 1
          KERR(1) =  'NOM INCONNU '//KNOM0
          CALL LEREUR
          CALL LXIM( NTLX )
          RETURN
      ENDIF
C
C     CONVERSION DU NOM KNOM1 EN ENTIERS ET ECRASEMENT
C     DE L'ANCIEN NOM
      MN = MNLX + MCN(MNLX) * NONOM
      CALL NOMENT( MCN(MNLX+2) , KNOM1 , MCN( MN ) )
      END
