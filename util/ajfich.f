      SUBROUTINE AJFICH( NCVALS )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   AJOUTER UN FICHIER D ACCES DIRECT A LA MEMOIRE SECONDAIRE
C -----   DECLARATION ET OUVERTURE
C
C SORTIE :
C --------
C NMTCL : -1  SI ABANDON DE L'AJOUT DU FICHIER
C         >=0 SI PAS D'ERREUR
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS  NOVEMBRE 1989
C23456---------------------------------------------------------------012
      include"./incl/motmcg.inc"
      include"./incl/gsmenu.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      COMMON / MSSFTA / NOFISF,NBPASF,MOPASF,MGBUSF,NSFLIB,
     %                  M1FIMS,M2FIMS,MGFIMS,NSFIMS,LPFIMS,
     %                  M1TAMS,M2TAMS,MGBUTA,NBBUTA,NPTAMS,NATAMS,
     %                  NBCTMS,LLTAMS,LFTAMS,MGNPSF,NSFNPS,NPSNPS,
     %                  MGZLMG,MGZLMK,MGZLMN,MOTSMG,MOTSMK,MOTSMN,NTADAM
      NBPAGE = 64
      MOPAGE = 256
      NCVALS = 0
C
C     LECTURE DES DONNEES
 1    CALL LIMTCL( 'fichier' , NMTCL )
C     TRAITEMENT DE L'ABANDON EVENTUEL DE LA LECTURE
      IF( NMTCL .EQ. -1 ) THEN
         NCVALS = -1
         GOTO 9000
      ENDIF
C
      GOTO( 10 , 20 , 30 ) , NMTCL
C
C     NOMBRE_PAGES
 10   NCVALS = 0
      CALL INVITE( 70 )
      CALL LIRENT( NCVALS , NBPAGE )
      IF( NCVALS .EQ. -1 ) RETURN
      NBPAGE = MAX( NBPAGE , 16 )
      GOTO 1
C
C     MOTS_UNE_PAGE
 20   NCVALS = 0
      CALL INVITE( 69 )
      CALL LIRENT( NCVALS , MOPAGE )
      IF( NCVALS .EQ. -1 ) RETURN
      MOPAGE = MAX( MOPAGE , 16 )
      GOTO 1
C
C     EXECUTION
 30   IF( NBPAGE .LE. 0 .OR. MOPAGE .LE. 0 ) THEN
         KERR(1) = 'NOMBRE DE  PAGES OU MOTS PAR PAGE INCORRECT'
         NBLGRC(NRERR) = 1
         CALL LEREUR
         GOTO 1
      ENDIF
C
C     RECHERCHE DU PREMIER NUMERO DE FICHIER DISPONIBLE
C     LE PARCOURS DES FICHIERS NUMERIQUES
      IG = MGFIMS
      DO 40 NOFICH=1,M2FIMS
         IF( MCG( IG ) .GT. 0 ) THEN
            IG = IG + M1FIMS
            GOTO 40
         ELSE
C           LE FICHIER EST DECLARE ET OUVERT
            CALL MSFI( NOFICH , NBPAGE , MOPAGE )
            GOTO 9000
         ENDIF
 40   CONTINUE
      KERR(1) = 'TROP DE FICHIERS DECLARES'
      KERR(2) = 'AUGMENTER LE NOMBRE DE FICHIERS DECLARABLES'
      NBLGRC(NRERR) = 2
      CALL LEREUR
      NCVALS = -1
C
9000  RETURN
      END
