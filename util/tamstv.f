      SUBROUTINE TAMSTV( NOTAMS, KTYPE, NBVARI )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    RETOURNER LE NOM DE TYPE KTYPE ET LE NOMBRE DE VARIABLES
C -----    NBVARI DU TABLEAU MS DE NUMERO NOTAMS
C
C ENTREE :
C --------
C NOTAMS : NUMERO DU TABLEAU MS DE TYPE ET DE NOMBRE DE VARIABLES A DECRIRE
C
C SORTIES :
C ---------
C KTYPE  : NOM DU TYPE DU TABLEAU MS DE NUMERO NOTAMS 9 CARACTERES
C          LOGIQUE , CARACTERE, ENTIER/2 , ENTIER ,  REEL,
C          REEL2   , REEL4    , COMPLEXE , COMPLEXE2, MOTS,
C          ^LEXIQUE, XYZ
C NBVARI : NOMBRE DE VARIABLES DE CE TABLEAU MS NOTAMS
C          0 SI LE TABLEAU N'EST PAS DECLARE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : PERRONNET ALAIN ANALYSE NUMERIQUE PARIS    DECEMBRE 1983
C.......................................................................
      include"./incl/motmcg.inc"
      CHARACTER*(*)  KTYPE
      CHARACTER*9    TYPNUM, KVIDE
      DATA           KVIDE/ 'VIDE' /
C
C     RECHERCHE DU DESCRIPTEUR DU TABLEAU MS NOTAMS
      CALL TAMSRE( NOTAMS, MGTAMS )
C
C     LE TYPE DU TABLEAU
      NUMTYV = ABS( MCG(MGTAMS) )
      IF( NUMTYV .EQ. 0 ) THEN
         KTYPE = KVIDE
      ELSE
         KTYPE = TYPNUM( NUMTYV )
      ENDIF
C
C     LE NOMBRE DE SES VARIABLES
      NBVARI = MCG( MGTAMS + 1 )
C
ccc      print *,'tamstv: numtyv=',numtyv,' mgtams=',mgtams,
ccc     %        ' ktype=',ktype,' nbvari=',nbvari
C
      RETURN
      END
