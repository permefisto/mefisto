        SUBROUTINE PLSV82( NUTYOB , NUPLSV , NTPLSV , LADEFI ,
     %                     IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    DETRUIRE UN POINT, LIGNE, SURFACE, VOLUME
C -----
C          ATTENTION : EN SORTIE LE LEXIQUE DE CE PLSV EST DETRUIT
C                                NTPLSV=0
C
C ENTREES:
C --------
C NUTYOB : NUMERO DU TYPE DU P L S V   1:POINT, 2:LIGNE, ...
C NUPLSV : NUMERO DU P L S V DANS SON LEXIQUE
C LADEFI : TABLEAU DE DEFINITION DU P L S V
C          CF '~TD/D/A_..._DEFINITION'
C
C MODIFIE :
C ---------
C NTPLSV : NUMERO 0 DU TABLEAU TS DU LEXIQUE DU P L S V
C          EN ENTREE LEXIQUE DU P L S V
C          EN SORTIE =0
C SORTIE :
C --------
C IERR   : 82>0 NECESSAIRE POUR NE PAS TRACER UN PLSV DETRUIT
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : A.PERRONNET LABORATOIRE ANALYSE NUMERIQUE PARIS DECEMBRE 1990
C23456...............................................................012
      IMPLICIT INTEGER (W)
      include"./incl/gsmenu.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/a_point__definition.inc"
C
      INTEGER           LADEFI(0:*)
      CHARACTER*24      KNOM
      CHARACTER*10      NMTOBJ,NMTYOB
C
C     LE NUMERO DE TMS DU LEXIQUE DE TYPE NUTYOB
      IF( NUTYOB .LE. 0 .OR. NUTYOB .GT. 5 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'PLSV82:TYPE DIFFERENT DE P L S V'
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF
C
C     LE NOM DU TYPE DU P L S V
      NMTOBJ = NMTYOB( NUTYOB )
      NTLXLX = NTMN( NUTYOB )
      IF( NTLXLX .LE. 0 ) THEN
         NBLGRC(NRERR) = 2
         KERR(1) = 'PLSV82:' // NMTOBJ
         KERR(2) = 'DE LEXIQUE INCONNU'
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF
C
C     VERIFICATION DU TYPE DU P L S V
      IF( LADEFI(WUTYPO) .NE. 82 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'PLSV82:TYPE DE P L S V DIFFERENT DE 82'
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF
C
C     LE  P L S V A DETRUIRE EXISTE-T-IL ?
      CALL LXNLOU( NTLXLX , NUPLSV , NT0 , MN0 )
      IF( NT0 .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IERR = NUDCNB( NMTOBJ )
         IF( NUTYOB .EQ. 2 .OR. NUTYOB .EQ. 3 ) THEN
            KERR(1) = NMTOBJ(1:IERR-1) // ' ANCIENNE INCONNUE'
         ELSE
            KERR(1) = NMTOBJ(1:IERR-1) // ' ANCIEN INCONNU'
         ENDIF
         CALL LEREUR
         IERR = 2
         RETURN
      ENDIF
C
C     DESTRUCTION DU LEXIQUE DU P L S V
      CALL NMOBNU( NMTOBJ , NUPLSV , KNOM )
C
C     DESTRUCTION DU POINT NOUVEAU
      CALL LXLXDS( NTLXLX , KNOM )
      NTPLSV = 0
      IERR   = 82
      END
