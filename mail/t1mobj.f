      SUBROUTINE T1MOBJ( NMTOBJ, NOMOBJ, NUMOBJ )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACER LE PLSV NOMOBJ DE TYPE NMTOBJ
C------

C ENTREES:
C --------
C NMTOBJ : NOM DU TYPE DE L'OBJET 'POINT' 'LIGNE' ... 'VOLUME'
C NOMOBJ : NOM DU PLSV
C NUMOBJ : NUMERO DE PLSV DANS LE LEXIQUE DE CE TYPE DE PLSV
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     OCTOBRE 1988
C....................................................................012
      include"./incl/langue.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/gsmenu.inc"
      COMMON /UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      CHARACTER*(*)    NMTOBJ, NOMOBJ

C     NUMERO DU TYPE D'OBJET A TRACER
C     -------------------------------
      NOTYOB = NTYOBJ( NMTOBJ )
      IF( NOTYOB .LE. 0 .OR. NOTYOB .GE. 5 ) RETURN

C     RECHERCHE OU GENERATION DU MAILLAGE DE L'OBJET
C     ----------------------------------------------
      NBLGRC(NRERR) = 0
      CALL MAILEX( NMTOBJ, NOMOBJ, NTTSMA, MNTSMA,
     %             NTSOMM, MNSOMM, IERR )
C     SI SORTIE DE LA DESTRUCTION DU PLSV (IERR=82)
      IF( IERR .EQ. 82 ) RETURN
      IF( IERR .NE. 0  ) THEN
         NBLGRC(NRERR) = NBLGRC(NRERR) + 2
         IF( LANGAG .EQ. 0 ) THEN
          KERR(NBLGRC(NRERR)-1) = NMTOBJ // ' SANS MAILLAGE: ' // NOMOBJ
          KERR(NBLGRC(NRERR)  ) = NOMOBJ // ' EST SUPPRIME'
         ELSE
           KERR(NBLGRC(NRERR)-1) = NMTOBJ // ' WITHOUT MESH: ' // NOMOBJ
           KERR(NBLGRC(NRERR)  ) = NOMOBJ // ' IS DELETED'
         ENDIF
         CALL LEREUR
C        DESTRUCTION DU PLSV NON MAILLE
         CALL LXLXDS( NTMN(NOTYOB), NOMOBJ )
         GOTO 9999
      ENDIF

C     TRACE EFFECTIF DE L'OBJET
C     -------------------------
      IF( INTERA .GE. 1 .AND. NOTYOB .GT. 0 .AND. MNSOMM .GT. 0 ) THEN
         CALL T1SOBJ( NOTYOB, NOMOBJ, NUMOBJ,
     %                NTTSMA, MNTSMA, NTSOMM, MNSOMM )
      ENDIF

 9999 RETURN
      END
