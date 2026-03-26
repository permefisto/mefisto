      SUBROUTINE LXTSFE( NTLX , KNOM )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : FERMER LE NOM KNOM TABLEAU MEMOIRE SECONDAIRE DU LEXIQUE NTLX
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
      include"./incl/pp.inc"
      COMMON   MCN(MOTMCN)
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      CHARACTER*(*) KNOM
      CHARACTER*4   KTYPE,CHARX
C
C     RECHERCHE DU NOM KNOM PARMI LES NOMS DU LEXIQUE NTLX
      CALL LXNMNO( NTLX , KNOM , NO , MNLX )
C
C     LE NOM A T IL ETE RETROUVE ?
      IF( NO .LE. 0 ) THEN
C        NON
         NBLGRC(NRERR) = 2
         KERR(1) = KNOM
         KERR(2) ='NOM INCONNU DANS LE LEXIQUE'
         CALL LEREUR
         CALL LXIM( NTLX )
         RETURN
      ENDIF
C
C     OUI
      MN     = MNLX + MCN( MNLX ) * NO + MCN( MNLX+2 ) + 1
      KTYPE  = CHARX( MCN(MN) )
      IF( KTYPE .NE. 'TAMS' ) THEN
C        TYPE INCORRECT =/ TAMS
         NBLGRC(NRERR) = 2
         KERR(1) = KNOM
         KERR(2) = 'NOM DE TYPE '//KTYPE//' AU LIEU DE TAMS'
         CALL LEREUR
         CALL ARRET(100)
      ENDIF
C
C     FERMETURE EFFECTIVE
      CALL TAMSFE( MCN( MN + 1 ) )
      END
