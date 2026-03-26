      SUBROUTINE NUOBNM( NMTOBJ, KNMOBJ, NUMOBJ )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    RETOURNER LE NUMERO D'UN OBJET CONNU PAR SON TYPE D'OBJET
C -----    ET SON NOM

C ENTREES:
C --------
C NMTOBJ : NOM DU TYPE DE L'OBJET
C KNMOBJ : NOM DE L'OBJET

C SORTIE :
C --------
C NUMOBJ : >0 NUMERO DE L'OBJET DANS SON LEXIQUE
C          =0 SI NON RETROUVE
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS  OCTOBRE 1988
C....................................................................
      include"./incl/langue.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      CHARACTER*(*)     NMTOBJ, KNMOBJ

C     LE NUMERO DE TMS DU LEXIQUE DU TYPE D'OBJET
      NOTYOB = NTYOBJ( NMTOBJ )
      IF( NOTYOB .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'NUOBNM: TYPE OBJET INCORRECT ' // NMTOBJ
         ELSE
            KERR(1) = 'NUOBNM: UNKNOWN TYPE of OBJECT ' // NMTOBJ
         ENDIF
         CALL LEREUR
         RETURN
      ENDIF

C     LE NUMERO DE L'OBJET DANS SON LEXIQUE
      CALL LXNMNO( NTMN( NOTYOB ), KNMOBJ, NUMOBJ, I )
      RETURN
      END
