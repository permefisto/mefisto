      SUBROUTINE NMOBNU( NMTOBJ , NUOBJE , KNMOBJ )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    RETOURNER LE NOM D'UN OBJET CONNU PAR SON TYPE D'OBJET
C -----    ET SON NUMERO
C
C ENTREES:
C --------
C NMTOBJ : NOM DU TYPE DE L'OBJET
C NUOBJE : NUMERO DE L'OBJET DANS SON LEXIQUE
C
C SORTIE :
C --------
C KNMOBJ : NOM DE L'OBJET
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS  OCTOBRE 1988
C....................................................................
      include"./incl/ntmnlt.inc"
      include"./incl/gsmenu.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      CHARACTER*(*)     NMTOBJ , KNMOBJ
      INTEGER           NUOBJE
      CHARACTER*24      KNOM
C
C     LE NUMERO DE TMS DU LEXIQUE DU TYPE D'OBJET
      NOTYOB = NTYOBJ( NMTOBJ )
      IF( NOTYOB .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'NMOBNU: TYPE OBJET INCORRECT '//NMTOBJ
         CALL LEREUR
         RETURN
      ENDIF
C
C     L'ADRESSE DU LEXIQUE DE CE TYPE D'OBJET
      CALL TAMSOU( NTMN(NOTYOB) , MNLX )
C
C     L'ADRESSE DU 1-ER MOT DE CET OBJET DANS LE LEXIQUE
      MN = MNLX + MCN( MNLX ) * NUOBJE
C
C     LE NOM CONVERTI D'ENTIERS EN CARACTERES
C     KNOM OFFRE UNE PROTECTION DE DEBORDEMENT
      CALL ENTNOM( MCN(MNLX+2) , MCN(MN) , KNOM )
C
C     LE NOM EST REPORTE DANS KNMOBJ
      KNMOBJ = KNOM
C
      RETURN
      END
