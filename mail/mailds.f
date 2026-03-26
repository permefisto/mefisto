      SUBROUTINE MAILDS( NMTOBJ , NOMOBJ )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    DETRUIRE LE MAILLAGE D'UN OBJET
C------
C
C ENTREES:
C --------
C NMTOBJ : NOM DU TYPE DE L'OBJET 'POINT' 'LIGNE' ... 'OBJET'
C NOMOBJ : NOM DE L'OBJET
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS  OCTOBRE 1988
C.......................................................................
      CHARACTER*(*)     NMTOBJ,NOMOBJ
      include"./incl/ntmnlt.inc"
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (RMCN(1),MCN(1))
C
C     NUMERO DU TYPE D'OBJET A MAILLER 1:POINT 2:LIGNE ...
      NOTYOB = NTYOBJ( NMTOBJ )
C
C     NUMERO TS ET ADRESSE DU LEXIQUE DE CE TYPE D'OBJET
      NTLX   = NTMN( NOTYOB )
C
C     NUMERO TS DU LEXIQUE DE CET OBJET
      CALL LXLXOU( NTLX , NOMOBJ , NTLXOB , I )
      IF( NTLXOB .LE. 0 ) RETURN
C
C     DESTRUCTION DU TABLEAU 'NSEF'
      CALL LXTSOU( NTLXOB ,  'NSEF' , NT , MN )
      IF( NT .GT. 0 ) THEN
         CALL LXTSDS( NTLXOB ,  'NSEF' )
      ENDIF
C
C     DESTRUCTION DU TABLEAU 'XYZSOMMET'
      CALL LXTSOU( NTLXOB ,  'XYZSOMMET' , NT , MN )
      IF( NT .GT. 0 ) THEN
         CALL LXTSDS( NTLXOB ,  'XYZSOMMET' )
      ENDIF
C
C     CET OBJET EST-IL UNE UNION D'OBJETS ?
      CALL LXTSOU( NTLXOB , 'UNION' , NT , MN )
      IF( NT .GT. 0 ) THEN
C        OUI : DESTRUCTION DU TABLEAU 'UNION'
         CALL LXTSDS( NTLXOB , 'UNION' )
      ENDIF
      END
