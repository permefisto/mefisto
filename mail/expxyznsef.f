      SUBROUTINE EXPXYZNSEF
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   EXPORTER SUR UN FICHIER de  CARACTERES ASCII
C -----   les COORDONNEES DES SOMMETS et TANGENTES
C         les NUMEROS DES SOMMETS DU MAILLAGE D'UN PLSV
C         C-A-D les TMS XYZSOMMET et NSEF Mefisto

C SORTIE :
C --------
C         LE FICHIER ASCII xyznsef.PLSV.KNMPLSV est ECRIT dans le
C         REPERTOIRE du PROJET
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     OCTOBRE 1997
C MODIFS : Alain PERRONNET Saint PIERRE du PERRAY             Avril 2020
C23456---------------------------------------------------------------012
      include"./incl/gsmenu.inc"
      include"./incl/langue.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/pp.inc"
      COMMON             MCN(MOTMCN)
      REAL               RMCN(1)
      EQUIVALENCE       (MCN(1),RMCN(1))
      CHARACTER*10       NMTYOB
      CHARACTER*24       KNMPLSV
      CHARACTER*1        KSUFIX(4)
      DATA               KSUFIX / 'p', 'l', 's', 'v' /

C     NOM DU TYPE DU PLSV A TRAITER
 10   CALL LIMTCL( 'typ_plsv', NUTYOB )
      IF( NUTYOB .LT. 0 ) GOTO 9999

C     NOM DU PLSV A TRAITER
      GOTO( 11, 12, 13, 14 ),NUTYOB
 11   CALL INVITE( 51 )
      GOTO 18
 12   CALL INVITE( 40 )
      GOTO 18
 13   CALL INVITE( 42 )
      GOTO 18
 14   CALL INVITE( 60 )
C
 18   IERR   = 0
      NCVALS = 0
      CALL LIRCAR( NCVALS, KNMPLSV )
      IF( NCVALS .EQ. -1 ) GOTO 9999

C     RECHERCHE DU NOM DE L'OBJET DANS LE LEXIQUE DES OBJETS
      CALL LXLXOU( NTMN(NUTYOB), KNMPLSV, NTLXOB, MNLXOB )
C     S'IL N'EXISTE PAS RETOUR A LA DEMANDE DU NOM DE L'OBJET
      IF( NTLXOB .LE. 0 ) THEN
         NBLGRC(NRERR) = 2
         KERR(1) = KNMPLSV
         KERR(2) ='ERREUR:' // NMTYOB(NUTYOB) // 'INCONNU'
         CALL LEREUR
         CALL LXIM( NTMN(NUTYOB) )
         GOTO 10
      ENDIF
      PRINT*,'expxyznsef: EXPORT ',NMTYOB(NUTYOB),' ', KNMPLSV

C     RECHERCHE DU TABLEAU XYZSOMMET
      CALL LXTSOU( NTLXOB, 'XYZSOMMET', NTXYZS, MNXYZS )
C     S'IL N'EXISTE PAS RETOUR A LA DEMANDE DU NOM DE L'OBJET
      IF( NTXYZS .LE. 0 ) THEN
         NBLGRC(NRERR) = 2
         KERR(1) = KNMPLSV
         IF( LANGAG .EQ. 0 ) THEN
            KERR(2) = 'PAS DE TMS XYZSOMMET'
         ELSE
            KERR(2) = 'TMS XYZSOMMET UNKNOWN'
         ENDIF
         CALL LEREUR
         GOTO 10
      ENDIF

C     RECHERCHE DU TABLEAU NSEF
      IF( NUTYOB .GT. 1 ) THEN
         CALL LXTSOU( NTLXOB, 'NSEF', NTNSEF, MNNSEF )
C        S'IL N'EXISTE PAS RETOUR A LA DEMANDE DU NOM DE L'OBJET
         IF( NTNSEF .LE. 0 ) THEN
            NBLGRC(NRERR) = 2
            KERR(1) = KNMPLSV
            IF( LANGAG .EQ. 0 ) THEN
               KERR(2) = 'PAS DE TMS NSEF'
            ELSE
               KERR(2) = 'TMS NSEF UNKNOWN'
            ENDIF
            CALL LEREUR
            GOTO 10
         ENDIF
      ELSE
C        POINT
         NTNSEF = 0
         MNNSEF = 0
      ENDIF

C     SANS ou AVEC les TANGENTES dans le FICHIER?
C     ===========================================
      CALL LIMTCL( 'tgoupas', NOTG )
      IF( NOTG .LT. 0 ) GOTO 9999

C     DECLARATION OUVERTURE ECRITURE du FICHIER xyznsef.NOM_PLSV
C     ==========================================================
      CALL SAVXYZNSEF( NUTYOB, KNMPLSV, MNXYZS, MNNSEF, IERR )
      GOTO 10

 9999 RETURN
      END
