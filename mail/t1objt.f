      SUBROUTINE T1OBJT( NOMOBJ )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACER LES OBJETS PREMIERS PLSV DE L'OBJET NOMOBJ
C------    SANS INTERPOLATION
C
C ENTREES:
C --------
C NOMOBJ : NOM DE L'OBJET
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS  NOVEMBRE 1989
C.......................................................................
      IMPLICIT         INTEGER (W)
      include"./incl/gsmenu.inc"
      COMMON /UNITES / LECTEU , IMPRIM , NUNITE(30)
      include"./incl/pp.inc"
      COMMON           MCN(MOTMCN)
      include"./incl/ntmnlt.inc"
      include"./incl/a_objet__definition.inc"
      CHARACTER*(*)    NOMOBJ
      CHARACTER*24     KNOMOB
      CHARACTER*10     NMTYOB,KNOMTY
C
C     RECHERCHE DE L'OBJET
      CALL LXLXOU( NTOBJE , NOMOBJ , NTLXOB , MNLXOB )
      IF( NTLXOB .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'OBJET INCONNU: '//NOMOBJ
         CALL LEREUR
         RETURN
      ENDIF
C
C     OUVERTURE DU TABLEAU  'DEFINITION' DE L'OBJET
      CALL LXTSOU( NTLXOB , 'DEFINITION' , NTDFOB , MNDFOB )
      IF( NTDFOB .LE. 0 ) THEN
C        TABLEAU 'DEFINITION' INEXISTANT
         NBLGRC(NRERR) = 1
         KERR(1) = 'OBJET NON DEFINI:'//NOMOBJ
         CALL LEREUR
         RETURN
      ENDIF
C
C     LES OBJETS DE L'OBJET SONT EXPRIMES EN TERMES DE
C     POINTS LIGNES SURFACES VOLUMES  EXCLUSIVEMENT
C     C'EST A DIRE D'OBJETS 'PREMIERS'
      MNOBPR = 0
      CALL OBJPRE( NOMOBJ , NBOBPR , MOOBPR , MNOBPR ,  IERR )
      IF( IERR .NE. 0 ) RETURN
C
C     CADRE AUTOMATIQUE DE L'OBJET
CCC      CALL VISEE0
C
C     LE TRACE DES OBJETS DE L'OBJET
      MN = MNOBPR
      DO 10 I=1,NBOBPR
C        LE NO DU TYPE ET LE NUMERO DE L'OBJET
         NTY = MCN( MN     )
         NUO = MCN( MN + 1 )
C        LE NOM DU TYPE
         KNOMTY = NMTYOB( NTY )
C        LE NOM KNOMOB DE L'OBJET
         CALL NMOBNU( KNOMTY , NUO , KNOMOB )
         CALL T1MOBJ( KNOMTY , KNOMOB , NUO )
         MN = MN + 2
 10   CONTINUE
C
C     TRACE DU TITRE ET FERMETURE
      CALL TRFINS( NOMOBJ )
      END
