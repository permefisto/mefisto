      SUBROUTINE LXRLFE( NTLX , KNOM )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : FERMER LE NOM KNOM RELATION DU LEXIQUE NTLX
C -----
C
C ENTREES :
C ---------
C NTLX   : NUMERO DU TABLEAU MS CONTENANT LE LEXIQUE
C KNOM   : CHAINE DE CARACTERES NOM A RETROUVER DANS LE LEXIQUE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS   OCTOBRE 1985
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      CHARACTER*(*)     KNOM
      CHARACTER*4       KTYPE,CHARX
C
C     RECHERCHE DU NOM KNOM PARMI LES NOMS DU LEXIQUE NTLX
      CALL LXNMNO( NTLX , KNOM , NO , MNLX )
C
C     LE NOM A T IL ETE RETROUVE ?
      IF( NO .LE. 0 ) THEN
C
C        NON
C        ===
         NBLGRC(NRERR) = 2
         KERR(1) = KNOM
         KERR(2) ='NOM INCONNU DANS LE LEXIQUE'
         CALL LEREUR
         CALL LXIM( NTLX )
         RETURN
      ENDIF
C
C     OUI
C     ===
C     ADRESSE MCN DU DEBUT DU NOM
      MN     = MNLX + MCN( MNLX ) * NO + MCN( MNLX+2 ) + 1
      KTYPE  = CHARX( MCN(MN) )
      IF( KTYPE .NE. 'RELA' ) THEN
C        TYPE INCORRECT =/ RELA
         NBLGRC(NRERR) = 2
         KERR(1) = KNOM
         KERR(2) = 'NOM DE TYPE '//KTYPE//' AU LIEU DE RELA'
         CALL LEREUR
         RETURN
      ENDIF
C
C     FERMETURE EFFECTIVE
      CALL TAMSFE( MCN( MN + 1 ) )
      END
