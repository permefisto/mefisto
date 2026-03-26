      SUBROUTINE LIRTPLSV( NOTYPE, KNOM, NUPLSV )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     LIRE LE TYPE ET NOM D'UN PLSV
C -----     RETOURNER SON TYPE SON NOM ET SON NUMERO

C SORTIES :
C ---------
C NOTYPE : 1:POINT 2:LIGNE 3:SURFACE 4:VOLUME 5:OBJET
C         -1 ABANDON DE LA LECTURE PAR ENTREE DE @ ou ECHAPPEMENT
C KNOM   : NOM DU PLSVO en MAJUSCULES
C NUPLSV : NUMERO DANS SON LEXIQUE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  LJLL ANALYSE NUMERIQUE UPMC PARIS  JUIN 2004
C23456---------------------------------------------------------------012
      include"./incl/td.inc"
      include"./incl/typobj.inc"
      include"./incl/gsmenu.inc"
      include"./incl/pilect.inc"
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      COMMON / MSSFTA / MSSF(28),NTADAM
      CHARACTER*10      NMTYOB, KNOMTY
      CHARACTER*(*)     KNOM

C     LECTURE DU TYPE D'OBJET
 10   NCVALS = 0
      CALL INVITE( 91 )
      CALL LIRCAR( NCVALS, KNOM )
      IF( NCVALS .EQ. -1 ) THEN
         NOTYPE = -1
         RETURN
      ENDIF

      CALL MAJUSC( KNOM )
      NOTYPE = NTYOBJ( KNOM )
      IF( NOTYPE .EQ. 0 ) THEN
C        ERREUR A REDONNER
         GOTO 10
      ENDIF

C     OUVERTURE DU LEXIQUE DANS LE LEXIQUE ADAM
      KNOMTY = NMTYOB( NOTYPE )
      CALL LXLXOU( NTADAM, KNOMTY, NTLX, MNLX )

C     LECTURE DU NOM DE L'OBJET DANS SON LEXIQUE
      GOTO( 11, 12, 13, 14, 15 ),NOTYPE
C     POINT
 11   CALL INVITE( 51 )
      GOTO 30
C     LIGNE
 12   CALL INVITE( 40 )
      GOTO 30
C     SURFACE
 13   CALL INVITE( 42 )
      GOTO 30
C     VOLUME
 14   CALL INVITE( 60 )
      GOTO 30
C     OBJET
 15   CALL INVITE( 45 )

C     LECTURE DU NOM DU PLSVO
 30   NCVALS = 0
      CALL LIRCAR( NCVALS, KNOM )
      IF( NCVALS .EQ. -1 ) THEN
         NOTYPE = -1
         RETURN
      ENDIF

C     OUVERTURE DANS LE LEXIQUE APRES MISE EN MAJUSCULES du NOM
      CALL MAJUSC( KNOM )
      CALL LXNMNO( NTLX , KNOM , NUPLSV , N2 )
      IF( NUPLSV .LE. 0 ) THEN
C        NOM INCORRECT NON RETROUVE DANS LE LEXIQUE
         IF (LHLECT.GT. 1 .AND. LHLECT.LE.MXLECT ) THEN
C           LECTURE ACTUELLE SUR UN FICHIER AVEC UNE ERREUR
C           AFFICHAGE DU RETOUR AU CLAVIER
            CALL QUIFIC
         ELSE
            NBLGRC(NRERR) = 0
         ENDIF
         N2 = NUDCNB( KNOM )
         KERR(NBLGRC(NRERR)+1) = 'NOM INCONNU: ' // KNOM(1:N2)
         KERR(NBLGRC(NRERR)+2) = 'A choisir PARMI:'
         NBLGRC(NRERR) = NBLGRC(NRERR) + 2
         CALL LXIM0( MNLX )
         GOTO 30
      ENDIF

      RETURN
      END
