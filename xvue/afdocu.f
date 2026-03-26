      SUBROUTINE AFDOCU
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     AFFICHER LA DOCUMENTATION ASSOCIE A L'IDENTIFICATEUR COURANT
C -----     SON NUMERO NUIDEN SE TROUVE DANS LE COMMON / ITDESC /
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS        MARS 1996
C23456---------------------------------------------------------------012
      include"./incl/td.inc"
      include"./incl/homdir.inc"
      include"./incl/gsmenu.inc"
      include"./incl/langue.inc"
      include"./incl/lu.inc"
C.......................................................................
      COMMON /UNITES/ LECTEU,IMPRIM,INTERA,NFDOCU,NUNITE(28)
      CHARACTER*160   KNOM
      CHARACTER*24    KKIDE
      CHARACTER*1     CARLU, CHAR
C
      IF( INTERA .LE. 0 ) GOTO 9999
C
C     INTERACTIF
C     SAUVEGARDE DE LA FENETRE ACTUELLE
      IF( INTERA .GE. 3 ) CALL SAUVEFENETRE
C
C     EXISTE_T-IL UN IDENTIFICATEUR POUR CETTE DEMANDE DE DOC?
      CARLU  = ' '
      LIGLU0 = 0
      LIGLUE = 0
      IFIN   = 1
      IF( NUIDEN .LE. 0 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            KDOCU(1) = 'DESOLE: ICI PAS DE DOCUMENTATION'
         ELSE
            KDOCU(1) = 'SORRY: NO DOCUMENTATION HERE'
         ENDIF
         NBLGRC(NRDOCU) = 1
         GOTO 45
      ENDIF
C
C     LA LECTURE DU FICHIER DOCUMENTATION ASSOCIE A L'IDENTIFICATEUR
C     ==============================================================
      KKIDE = KIDENT( NUIDEN )
      CALL MINUSC( KKIDE )
      KNOM  = HOMDIR // '/doc/' // KKIDE
      OPEN( FILE=KNOM, UNIT=NFDOCU , STATUS='OLD' , IOSTAT=IOERR )
      IF( IOERR .NE. 0 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            KDOCU(1)='DESOLE: '//KIDENT(NUIDEN)//' SANS DOCUMENTATION'
         ELSE
            KDOCU(1)='SORRY: NO DOCUMENTATION FROM '//KIDENT(NUIDEN)
         ENDIF
         NBLGRC(NRDOCU) = 1
         GOTO 45
      ENDIF
C
C     LECTURE ET AFFICHAGE DE LA DOCUMENTATION
 10   LIGLU0 = LIGLUE
      DO 30 I = 1 , MXLGDO-1
         READ(UNIT=NFDOCU, FMT='(A)', END=40, IOSTAT=IOERR) KDOCU(I)
         IF( IOERR .NE. 0 ) GOTO 9900
 30   CONTINUE
C     IL RESTE DE LA DOCUMENTATION A AFFICHER
      IFIN   = 0
      LIGLUE = MXLGDO - 1
      NBLGRC(NRDOCU) = LIGLUE
      GOTO 45
C
C     C'EST LA DERNIERE PAGE DE LA DOCUMENTATION
C     AFFICHAGE DES NBLGRC(NRDOCU) LIGNES DE LA DOCUMENTATION
 40   LIGLUE = I - 1
      IF( LIGLUE .LE. 0 ) GOTO 9900
      NBLGRC(NRDOCU) = LIGLUE
      IFIN   = 1
C
 45   IF( INTERA .LE. 1 ) THEN
C
C        AFFICHAGE A L'ECRAN STANDARD
C        ============================
         DO 50 I = 1 , NBLGRC(NRDOCU)
            WRITE(IMPRIM,'(A)' ) KDOCU(I)
 50      CONTINUE
C
      ELSE
C
C        TRACE SUR L'ECRAN
C        =================
C        AFFICHAGE DU MOYEN DE SORTIR
         NBLGRC(NRDOCU) = MIN( NBLGRC(NRDOCU)+1 , MXLGDO )
         IF( LANGAG .EQ. 0 ) THEN
          KDOCU(NBLGRC(NRDOCU))='Pour SORTIR CLIQUER le BOUTON 2 ou TAPE
     %R Echappement'
         ELSE
            KDOCU(NBLGRC(NRDOCU)) = 'TO EXIT CLICK the BUTTON 2 or TYPE
     %Escape'
         ENDIF
C
C        LE TRACE DE LA DOCUMENTATION
         CALL RECTTX( NRDOCU , KDOCU , 0 , NA )
C
C        ATTENTE D'UN CLIC D'UN BOUTON DE LA SOURIS OU SAISIE D'UN CARACTERE
C        ===================================================================
 300     CALL SAIPTC ( NOCODE, NX, NY, NASCII )
         IF( NOCODE .EQ. 0 ) THEN
C
C           ABANDON
            GOTO 9900
C
         ELSE IF( NOCODE .GT. 0 ) THEN
C
C           CLIC DU BOUTON NOCODE DE LA SOURIS DANS L'ECRAN GRAPHIQUE EN X,Y
C           ----------------------------------------------------------------
C           GESTION DU POINT CLIQUE
            IF( (NX.GE.XRECT(NRDOCU)) .AND.
     &          (NX.LE.XRECT(NRDOCU)+DXRECT(NRDOCU)) .AND.
     &          (NY.GE.YRECT(NRDOCU)) .AND.
     &          (NY.LE.YRECT(NRDOCU)+DYRECT(NRDOCU)) ) THEN
C
C              CLIC DANS LE RECTANGLE DE DOCUMENTATION : LA SUITE ?
               IF( NOCODE .EQ. 1 ) THEN
C
C                 BOUTON 1 DE LA SOURIS => LECTURE DE LA PAGE SUIVANTE DE LA DOC
                  IF( IFIN .EQ. 0 ) THEN
C                    IL RESTE DE LA DOCUMENTATION A LIRE  => PASSAGE A LA PAGE
                     CALL RECTEF( NRDOCU )
                     GOTO 10
                  ELSE
C                    PLUS DE DOCUMENTATION => RETOUR EN SAISIE
                     GOTO 300
                  ENDIF
C
               ELSE IF( NOCODE .EQ. 2 ) THEN
C
C                 BOUTON 2 => SORTIE
                  GOTO 9900
C
               ELSE IF( NOCODE .EQ. 3 ) THEN
C
C                 BOUTON 3 => RETOUR EN ARRIERE DANS LA DOCUMENTATION
                  NASCII = 303
                  GOTO 310
C
               ENDIF
            ENDIF
C
C           CLIC EN DEHORS DES LIMITES DU RECTANGLE DE DOCUMENTATION => RETOUR
            GOTO 9900
C
         ELSE IF( NOCODE .EQ. 0 ) THEN
C
C           ERREUR DANS L'EVENEMENT SAISI
            GOTO 300
C
         ENDIF
C
C        ENTREE DU CARACTERE CARLU A L'AIDE DU CLAVIER PHYSIQUE
C        ------------------------------------------------------
         CARLU = CHAR( NASCII )
C
C        ABANDON AVEC LE CARACTERE @ OU '  ' OU  ECHAPPEMENT OU SUPPRESSION
         IF( CARLU  .EQ. '@' .OR. CARLU  .EQ. ' ' .OR.
     %       NASCII .EQ.  27 .OR. NASCII .EQ. 127 ) GOTO 9900
C
 310     IF( CARLU  .EQ. '^' .OR. CARLU  .EQ. '-' .OR.
     %       NASCII .EQ.  8  .OR. NASCII .EQ. 303 ) THEN
C
C           FLECHES VERS LE HAUT(^) OU - OU BACKSPACE(8) OU
C           BOUTON 3 SOURIS(303)
C           => REMONTEE DANS LA DOCUMENTATION D'UNE PAGE
C           NBLP = NOMBRE DE LIGNES A REMONTER SUR LE FICHIER
            NBLP = LIGLU0 + LIGLUE
C
            IF( NBLP .GT. 0 ) THEN
C
C              REMONTEE D'UNE PAGE DANS LA DOCUMENTATION
C              -----------------------------------------
               CALL RECTEF( NRDOCU )
               DO 320 I = 1 , NBLP
                  BACKSPACE( UNIT=NFDOCU , IOSTAT=IOERR )
                  IF( IOERR .NE. 0 ) THEN
                     REWIND( UNIT=NFDOCU , IOSTAT=IOERR)
                     GOTO 10
                  ENDIF
 320           CONTINUE
               LIGLUE = MXLGDO-1
               GOTO 10
C
            ENDIF
         ENDIF
C
C        AJOUTE POUR LECTURE D'UNE SEULE PAGE ET SORTIE
         IFIN = 1
      ENDIF
C
C     PASSAGE A LA PAGE SUIVANTE OU SORTIE
      IF( IFIN .EQ. 0 ) THEN
C        IL RESTE DE LA DOCUMENTATION A LIRE  => PASSAGE A LA PAGE SUIVANTE
         CALL RECTEF( NRDOCU )
         GOTO 10
      ENDIF
C
C     FIN DE LECTURE DE LA LISTE DES MOTS CLES
 9900 CLOSE( UNIT=NFDOCU )
      IF( INTERA .GE. 3 ) THEN
C        RESTAURATION DE LA FENETRE
         CALL RESTAUREFENETRE
      ENDIF
C
C     MISE A BLANC DE LA LIGNE DE SAISIE
 9999 KLG(1) = '  '
      LHKLG  = 0

      RETURN
      END
