      SUBROUTINE DONNMF( NOTYPE, NOTYPS, NCVALS, DBLVAL, KCHAIN )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : LIRE ET RETOURNER UNE VALEUR DONNEE PAR L'UTILISATEUR
C ----- CETTE VALEUR EST REELLE DOUBLE PRECISION  SI NCVALS=1
C       CETTE VALEUR EST UNE CHAINE DE CARACTERES SI NCVALS=2
C       CETTE VALEUR EST INCORRECTE               SI NCVALS=0
C
C       CF ~LU/GRAMMAIRE DE DEFINITION DU LANGAGE UTILISATEUR
C
C ENTREES :
C ---------
C NOTYPE  : =0  LE TYPE DE LA VARIABLE N'EST PAS INITIALISE EN ENTREE
C           >0  LE TYPE DE LA VARIABLE INITIALISEE ET CODEE DANS
C               DBLVAL SI TYPE NUMERIQUE
C               KCHAIN SI TYPE CARACTERE
C           =-2 UN FICHIER DE DONNEES EST FOURNI POUR UNE EXECUTION EN BATCH
C               ET KLG(1) = 'READF ' // NMFILEDO(1:N) // ' ; '
C
C           LES VALEURS DES DIFFERENTS TYPES POSSIBLES :
C              LOGIQUE  => 1 CARACTERE=> 2  ENTIER/2 => 3 ENTIER   => 4
C              REEL     => 5 REEL2    => 6  REEL4    => 7 COMPLEXE => 8
C              COMPLEXE2=> 9 MOTS     => 10 ^LEXIQUE =>11 XYZ      =>12
C              TYPEOBJET=>13                ( CF ~/incl/msvaau.inc )
C
C NOTYPS  : =0 AUCUN TYPE PRECIS N'EST DEMANDE EN SORTIE
C           >0 LE TYPE EST IMPOSE EN SORTIE
C
C ENTREES ET SORTIES :
C --------------------
C NCVALS : 0 DBLVAL ET KCHAIN NE SONT PAS INITIALISES
C          1 DBLVAL EST INITIALISEE
C          2 KCHAIN EST INITIALISEE
C         -1 ABANDON DE LA LECTURE D'UNE VALEUR
C DBLVAL : VALEUR REELLE DOUBLE PRECISION
C KCHAIN : CHAINE DE CARACTERES

C INTERA : 0:BATCH      PAS d' ECRAN GRAPHIQUE PAS de CLAVIER PAS de SOURIS
C          1:INTERACTIF AVEC   ECRAN GRAPHIQUE PAS de CLAVIER PAS de SOURIS
C          3:INTERACTIF AVEC   ECRAN GRAPHIQUE AVEC   CLAVIER AVEC   SOURIS
C
C ATTENTION : ? ET = DEMANDE AUX PARAMETRES D'ETRE INITIALISES
C                            EN ENTREE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS    JANVIER 1990
C MODIFS : PERRONNET ALAIN LJLL UPMC & St PIERRE DU PERRAY  AVRIL   2013
C23456---------------------------------------------------------------012
      include"./incl/lu.inc"
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/pilect.inc"
      include"./incl/trvari.inc"
      COMMON / UNITES /  LECTEU, IMPRIM, INTERA, NUNITE(29)
      DOUBLE PRECISION   DBLVAL
      CHARACTER*(*)      KCHAIN
      LOGICAL            LOGIQ, LEXIST
      CHARACTER*160      KNOMFI
C
C     VALEUR NON INITIALISEE
      NCVALS = 0
C
C     ============================================================
C     TRAITEMENT DE LA CHAINE DE CARACTERES STOCKEE OU A STOCKER
C     A PARTIR DE KLG(NLPTV1)(NCPTV1+1) JUSQU'A KLG(LHKLG)(NBCALI)
C     ============================================================
 10   IF( NLPTV1 .GT. 0 ) GOTO 18
C
C     RIEN DANS LE BUFFER LIGNES KLG
C     ------------------------------
 15   LHKLG  = 0
      NLPTV1 = 1
      NCPTV1 = 0
C
C     LECTURE DE LA LIGNE SUIVANTE DES DONNEES UTILISATEUR
 16   CALL LIRLIG( I )
C     SI @ => ABANDON
 17   IF( I .EQ. -1 ) GOTO 9500
      IF( I .NE.  0 ) RETURN
C
C     ICI : LHKLG EST LA DERNIERE LIGNE CONTENANT DES DONNEES LUES
      IF( LHKLG .LE. 0 ) GOTO 15
C
C     L'ANALYSE DES CARACTERES AU DELA DE NLPTV1,NCPT1 DEBUT
 18   NL0 = NLPTV1
      NC0 = NCPTV1
C
C     RECHERCHE DU PROCHAIN CARACTERE NON BLANC
      CALL CARPNB( NL0 , NC0 )
      IF( NL0 .LE. 0 ) GOTO 15
C
C     PROTECTION DU DEBUT DE CHAINE
      NL00 = NL0
      NC00 = NC0
C
C     ==================================
C     EFFACEMENT DE KLG DES COMMENTAIRES
C     ==================================
 20   CALL LICMMT( NL0 , NC0 )
      IF( NL0 .EQ.  0 ) GOTO 9000
ccc      IF( NL0 .EQ. -1 ) GOTO 9500
      IF( NUDCNB( KLG(NL0) ) .EQ. 0 ) GOTO 16
C
C     RESTAURATION DU DEBUT DE CHAINE
      NL0 = NL00
      NC0 = NC00
C
C     ===============================================================
C     REPERAGE D'UN EVENTUEL @ CARACTERE D'ABANDON DE LA LECTURE
C     ===============================================================
 60   NCD = INDEX( KLG(NL0)(NC0:NBCALI) , '@' )
      IF( NCD .GT. 0 ) THEN
C        LA LECTURE DOIT ETRE ABANDONNEE PUISQUE @ A ETE LU
         GOTO 9500
      ENDIF
C
C     ===============================================================
C     RECHERCHE DE LA POSITION DU PREMIER ;  FIN DE LA DONNEE
C     ===============================================================
      NCPTV2 = INDEX( KLG(NL0)(NC0:NBCALI) , ';' )
      IF( NCPTV2 .EQ. 0 ) THEN
         IF( NL0 .GE. LHKLG ) THEN
C           DERNIERE LIGNE LUE SANS ;
            IF( LHLECT .EQ. 1 ) THEN
               IF( LHKLG .LT. MXLGER .AND. NL0 .LT. MXLGER ) THEN
                  DO 66 I=1,LHKLG
                     KERR(I) = KLG(I)
 66               CONTINUE
                  I = LHKLG + 1
               ELSE
                  KERR(1) = KLG(NL0)
                  I = 2
               ENDIF
               NBLGRC(NRERR) = I
               WRITE(KERR(MXLGER)(1:5),'(I5)') NL0
               IF( LANGAG .EQ. 0 ) THEN
                  KERR(I) = 'LIGNE'//KERR(MXLGER)(1:5)
     %                    //' NON TERMINEE PAR '';'' '
               ELSE
                  KERR(I) ='INPUT LINE'//KERR(MXLGER)(1:5)
     %                    //' UNFINISHED. LACK of '';'' '
               ENDIF
            ENDIF
C           LECTURE D'UNE LIGNE SUPPLEMENTAIRE
            CALL LIRLIG( I )
C           SI @ => ABANDON
            IF( I .EQ. -1 ) GOTO 9500
            IF( I .NE.  0 ) GOTO 9000
            NL0 = LHKLG
            NC0 = 1
            GOTO 20
         ELSE
C           PASSAGE A LA LIGNE SUIVANTE DE KLG
            NL0 = NL0 + 1
            NC0 = 1
            GOTO 60
         ENDIF
      ELSE
C        LE ; EST RETROUVE
C        NLPTV2,NCPTV2 EST LE DERNIER CARACTERE DE LA DONNEE A TRAITER
         NLPTV2 = NL0
         NCPTV2 = NC0 - 1 + NCPTV2
      ENDIF
C
C     ==========================================================
C     L'INTERPRETATION EST FAITE DE KLG(NL00  )(NC00:...   )  A
C                                   KLG(NLPTV2)(... :NCPTV2)
C     ==========================================================
      NL = NL00
      NC = NC00 - 1
C
      NLF = NLPTV2
      NCF = NCPTV2
C
C     RECHERCHE DU PREMIER MOT CLE
C     ============================
C     LECTURE DE L'ITEM SUIVANT
      IF( NL .GT. NLPTV2 .OR.
     %  ( NL .EQ. NLPTV2 .AND. NC .GE. NCPTV2 ) ) GOTO 8000
C
C     RECHERCHE DU PREMIER CARACTERE NON BLANC
      CALL CARPNB( NL , NC )
      IF( NL .EQ. 0 ) GOTO 9000
C
      IF( KLG(NL)(NC:NC) .EQ. ';' ) THEN
C        LE ; FINAL EST ATTEINT
         GOTO 8000
C
      ELSE IF( KLG(NL)(NC:NC) .EQ. '?' ) THEN
C        DEMANDE DE DOCUMENTATION
         CALL AFDOCU
C        DEMANDE D'AFFICHAGE DE LA VALEUR A DONNER
         IF( NOTYPE .EQ. 2 ) THEN
            NBLGRC(NRERR) = 1
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'NOM INITIAL='//KCHAIN
            ELSE
               KERR(1) = 'INITIAL NAME='//KCHAIN
            ENDIF
         ELSE IF( NOTYPE .GT. 0 .AND. NOTYPE .NE. 2 ) THEN
            NBLGRC(NRERR) = 1
            WRITE(KERR(2)(1:15),'(G15.7)') DBLVAL
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'VALEUR INITIALE='//KERR(2)(1:15)
            ELSE
               KERR(1) = 'INITIAL VALUE='//KERR(2)(1:15)
            ENDIF
         ELSE
            NBLGRC(NRERR) = 1
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'PAS DE VALEUR INITIALE'
            ELSE
               KERR(1) = 'NO INITIAL VALUE'
            ENDIF
         ENDIF
         GOTO 8000
C
      ELSE IF( KLG(NL)(NC:NC) .EQ. '=' ) THEN
C        LA VALEUR ENTREE NE DOIT PAS ETRE MODIFIEE
         IF( NOTYPE .EQ. 2 ) THEN
            IF( LHLECT .EQ. 1 ) THEN
               NBLGRC(NRERR) = 1
               IF( LANGAG .EQ. 0 ) THEN
                  KERR(1) = 'NOM EN RETOUR='//KCHAIN
               ELSE
                  KERR(1) = 'OUTPUT NAME='//KCHAIN
               ENDIF
            ENDIF
            NCVALS = 2
         ELSE IF( NOTYPE .GT. 0 .AND. NOTYPE .NE. 2 ) THEN
            IF( LHLECT .EQ. 1 ) THEN
               NBLGRC(NRERR) = 1
               WRITE(KERR(2)(1:15),'(G15.7)') DBLVAL
               IF( LANGAG .EQ. 0 ) THEN
                  KERR(1) = 'VALEUR EN RETOUR='//KERR(2)(1:15)
               ELSE
                  KERR(1) = 'OUTPUT VALUE='//KERR(2)(1:15)
               ENDIF
            ENDIF
            NCVALS = 1
         ELSE
            NBLGRC(NRERR) = 2
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'PAS DE VALEUR INITIALE PREVUE'
               KERR(2) = 'VOUS DEVEZ DONNER UNE VALEUR'
            ELSE
               KERR(1) = 'NO INITIAL VALUE'
               KERR(2) = 'YOU HAVE TO GIVE A VALUE'
            ENDIF
            CALL LEREUR
         ENDIF
         GOTO 8000
C
      ELSE IF( KLG(NL)(NC:NC+8) .EQ. 'INITIERPS' ) THEN
C        INITIALISER LE FICHIER DE RECUEIL DES INSTRUCTIONS
C        POSTSCRIPT POUR IMPRIMER PLUS TARD CE QUI VA ETRE TRACE SUR L'ECRAN
C        FERME EVENTUELLEMENT LE FICHIER ACTUEL PS
         CALL XVINITIERPS( 1 )
         LASOPS = 1
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'INITIALISATION DU FICHIER POSTSCRIPT'
         ELSE
            KERR(1) = 'POSTSCRIPT FILE INITIATION'
         ENDIF
         GOTO 15
C
      ELSE IF( KLG(NL)(NC:NC+9) .EQ. 'AVECMENUPS' ) THEN
         IF( LASOPS .EQ. 0 ) THEN
            NBLGRC(NRERR) = 1
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'INITIERPS doit preceder AVECMENUPS'
            ELSE
               KERR(1) = 'INITIERPS must precede AVECMENUPS'
            ENDIF
            GOTO 15
         ENDIF
C        AJOUTER LES TRACES DES MENUS EN POSTSCRIPT
         LASOPS = 2
         CALL XVPOSTSCRIPT(LASOPS)
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'AJOUT DU TRACE DES MENUS EN PS'
         ELSE
            KERR(1) = 'DRAWING OF MENUS IS ADDED TO PS'
         ENDIF
         GOTO 15
C
      ELSE IF( KLG(NL)(NC:NC+9) .EQ. 'SANSMENUPS' ) THEN
C        ARRETER LE TRACE DES MENUS EN POSTSCRIPT
         LASOPS = 1
         CALL XVPOSTSCRIPT(LASOPS)
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'ARRET DU TRACE DES MENUS EN PS'
         ELSE
            KERR(1) = 'DRAWING OF MENUS TO PS IS STOPPED'
         ENDIF
         GOTO 15
C
      ELSE IF( KLG(NL)(NC:NC+7) .EQ. 'SAUVERPS' ) THEN
C        SAUVERPS 'NOM_FICHIER_PS'
C        FERMER LE FICHIER NOM_FICHIER_PS CONTENANT TOUTES
C        LES INSTRUCTIONS POSTSCRIPT DEPUIS LE DERNIER APPEL A INITPS
C        LIRE A PARTIR DE MAINTENANT LE NOM DU FICHIER PS
         IF( LASOPS .EQ. 0 ) THEN
            NBLGRC(NRERR) = 1
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'ERREUR: INITIERPS N''A PAS PRECEDE SAUVERPS'
            ELSE
               KERR(1) = 'ERROR: INITIERPS DO''NT PRECEDE SAUVERPS'
            ENDIF
            CALL LEREUR
            GOTO 9500
         ENDIF
         NC  = NC + 8
         CALL LIUTCH( NL , NC , NLD , NCD , NLF , NCF )
         IF( NLF .LE. 0 ) GOTO 9500
C
C        RECHERCHE DU 1-ER CARACTERE NON BLANC A PARTIR DE (NL,NC)
         DO 75 I=NC,NBCALI
            IF( KLG(NLD)(I:I) .EQ. ';' ) THEN
               DO 74 J=NC,I-1
                  IF( KLG(NLD)(J:J) .NE. ' ' ) GOTO 77
 74            CONTINUE
C              CHAINE VIDE
               NBLGRC(NRERR) = 1
               IF( LANGAG .EQ. 0 ) THEN
                  KERR(1) ='Taper SAUVERPS NOM_DU_FICHIER ;'
               ELSE
                  KERR(1) ='Type SAUVERPS FILE_NAME ;'
               ENDIF
               CALL LEREUR
               GOTO 9500
            ENDIF
 75      CONTINUE
C        PAS DE ; . SA DEMANDE EST IMPOSEE
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) ='Taper SAUVERPS NOM_DU_FICHIER ;'
         ELSE
            KERR(1) ='Type SAUVERPS FILE_NAME ;'
         ENDIF
         CALL LEREUR
         GOTO 9500
C
 77      CALL MINUSC( KLG(NLD)(NCD:NCF) )
         KNOMFI = KLG(NLD)(NCD:NCF)  // '.eps'
         INQUIRE( FILE=KNOMFI , OPENED=LOGIQ , NUMBER=NF , IOSTAT=N )
         IF( LOGIQ .OR. N .NE. 0 ) THEN
            NBLGRC(NRERR) = 3
            KERR(2) = KLG(NLD)(NCD:NCF)
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) ='FICHIER DEJA OUVERT'
               KERR(3) ='Retaper SAUVERPS NOM_CORRECT_DU_FICHIER ;'
            ELSE
               KERR(1) ='ALREADY OPENED FILE'
               KERR(3) ='Type again SAUVERPS CORRECT_FILE_NAME'
            ENDIF
            CALL LEREUR
            GOTO 9500
         ENDIF
         LONGCH = NCF-NCD+1
         CALL XVSAUVERPS( KLG(NLD)(NCD:NCF) , LONGCH )
         IF ( LONGCH.NE.0 ) THEN
           IF ( LONGCH.EQ.-2 ) THEN
             NBLGRC(NRERR) = 2
             IF( LANGAG .EQ. 0 ) THEN
                KERR(1)='SAUVEGARDE DU FICHIER PS ' // KLG(NLD)(NCD:NCF)
     S                //'.eps IMPOSSIBLE'
                KERR(2) ='SANS DOUTE DEJA OUVERT ... A REFERMER'
             ELSE
                KERR(1)='IMPOSSIBILITY TO SAVE THE PS FILE '
     S                // KLG(NLD)(NCD:NCF) // '.eps'
                KERR(2) ='IF ALREADY OPENED => CLOSE IT'
             ENDIF
             CALL LEREUR
           ELSE
             LASOPS = 0
             NBLGRC(NRERR) = 2
             IF( LANGAG .EQ. 0 ) THEN
                KERR(1) ='SAUVEGARDE DU FICHIER DES QUALITES '
     S                 // KLG(NLD)(NCD:NCF) // '.qua IMPOSSIBLE'
                KERR(2) ='LES QUALITES SONT DANS TEMPORAIRE.QUA'
                KERR(3) ='POUR EVITER QU''ELLES SOIENT PERDUES'
                KERR(4) ='SAUVEGARDER CE FICHIER SOUS '
     S                 // KLG(NLD)(NCD:NCF) // '.qua'
             ELSE
                KERR(1) ='IMPOSSIBILITY TO SAVE THE QUALITY FILE'
     S                 // KLG(NLD)(NCD:NCF) // '.qua'
                KERR(2)='QUALITIES ARE STORED IN TEMPORAIRE.QUA'
                KERR(3)='TO AVOID THE LOSS SAVE RENAME THIS FILE'
                KERR(4)= KLG(NLD)(NCD:NCF) // '.qua'
             ENDIF
             CALL LEREUR
           ENDIF
         ELSE
           LASOPS = 0
           NBLGRC(NRERR) = 2
           IF( LANGAG .EQ. 0 ) THEN
              KERR(1) = 'SAUVEGARDE DU FICHIER POSTSCRIPT DE NOM'
           ELSE
              KERR(1) = 'SAFEGUARD OF THE PS FILE'
           ENDIF
           KERR(2) =  KLG(NLD)(NCD:NCF)
         ENDIF
         GOTO 15
C
      ELSE IF( KLG(NL)(NC:NC+9) .EQ. 'IMPRIMERPS' ) THEN
C        IMPRIMERPS 'NOM_FICHIER_PS'
C        ENVOYER SUR L'IMPRIMANTE LES INSTRUCTIONS POSTSCRIPT
C        STOCKEES DANS LE FICHIER NOM_FICHIER_PS
C        LIRE A PARTIR DE MAINTENANT LE NOM DU FICHIER PS
         NC  = NC + 10
         CALL LIUTCH( NL , NC , NLD , NCD , NLF , NCF )
         IF( NLF .LE. 0 ) GOTO 9500
C
C        RECHERCHE DU 1-ER CARACTERE NON BLANC A PARTIR DE (NL,NC)
         DO 85 I=NC,NBCALI
            IF( KLG(NLD)(I:I) .EQ. ';' ) THEN
               DO 84 J=NC,I-1
                  IF( KLG(NLD)(J:J) .NE. ' ' ) GOTO 87
 84            CONTINUE
C              CHAINE VIDE
               NBLGRC(NRERR) = 1
               IF( LANGAG .EQ. 0 ) THEN
                  KERR(1) ='Taper IMPRIMERPS NOM_DU_FICHIER ;'
               ELSE
                  KERR(1) ='Type IMPRIMERPS FILE_NAME ;'
               ENDIF
               CALL LEREUR
               GOTO 9500
            ENDIF
 85      CONTINUE
C        PAS DE ; . SA DEMANDE EST IMPOSEE
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) ='Taper IMPRIMERPS NOM_DU_FICHIER ;'
         ELSE
            KERR(1) ='Type IMPRIMERPS FILE_NAME ;'
         ENDIF
         CALL LEREUR
         GOTO 9500
C
 87      CALL MINUSC( KLG(NLD)(NCD:NCF) )
         KNOMFI = KLG(NLD)(NCD:NCF)  // '.eps'
         INQUIRE( FILE=KNOMFI , OPENED=LOGIQ , NUMBER=NF ,
     %            EXIST=LEXIST, IOSTAT=N )
         IF( .NOT. LEXIST .OR. LOGIQ .OR. N .NE. 0 ) THEN
            NBLGRC(NRERR) = 3
            KERR(2) = KLG(NLD)(NCD:NCF)
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) ='FICHIER INCONNU ou OUVERT ou avec PB'
               KERR(3) ='RETAPER IMPRIMERPS NOM_CORRECT_DU_FICHIER ;'
            ELSE
               KERR(1) ='UNKNOWN or OPENED or PROBLEM FILE'
               KERR(3) ='Type again IMPRIMERPS CORRECT_FILE_NAME ;'
            ENDIF
            CALL LEREUR
            GOTO 9500
         ENDIF
         IF( NLF .EQ. 0 ) GOTO 9500
         CALL MINUSC( KLG(NLD)(NCD:NCF) )
         IF( N .NE. 0 ) THEN
            NBLGRC(NRERR) = 3
            KERR(2) = KLG(NLD)(NCD:NCF)
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) ='FICHIER INCONNU OU DEJA OUVERT'
               KERR(3) ='RETAPER IMPRIMERPS NOM_CORRECT_DU_FICHIER ;'
            ELSE
               KERR(1) ='UNKNOWN or ALREADY OPENED FILE'
               KERR(3) ='Type again IMPRIMERPS CORRECT_FILE_NAME ;'
            ENDIF
            CALL LEREUR
            GOTO 9500
         ENDIF
         LONGCH = NCF-NCD+1
         CALL XVIMPRIMERPS( KLG(NLD)(NCD:NCF) , LONGCH )
         IF ( LONGCH .NE. 0 ) THEN
           IF ( LONGCH .LT. 0 ) THEN
             NBLGRC(NRERR) = 4
             IF( LANGAG .EQ. 0 ) THEN
             KERR(1) ='IMPRESSION DU FICHIER ' // KLG(NLD)(NCD:NCF)
     S                //'.eps'
             KERR(2) ='IMPOSSIBLE DE CETTE CONSOLE'
             KERR(3) ='TAPER DIRECTEMENT SOUS UNIX :'
             KERR(4) ='lpr -PNomImprimante '//KLG(NLD)(NCD:NCF)//'.eps'
             ELSE
             KERR(1) ='PRINTING of FILE ' // KLG(NLD)(NCD:NCF)
     S                //'.eps'
             KERR(2) ='IS NOT POSSIBLE FROM THIS DEVICE'
             KERR(3) ='Type DIRECTLY FROM UNIX :'
             KERR(4) ='lpr -PPrinterName '//KLG(NLD)(NCD:NCF)//'.eps'
             ENDIF
             CALL LEREUR
           ELSE
             NBLGRC(NRERR) = 4
             IF( LANGAG .EQ. 0 ) THEN
                KERR(1) ='IMPRESSION DU FICHIER ' // KLG(NLD)(NCD:NCF)
     S                   //'.eps IMPOSSIBLE'
                KERR(2) ='MAUVAIS ORDRE D''IMPRESSION'
                KERR(3) ='A CHANGER DANS LA PROCEDURE xvimprimerps'
                KERR(4) ='DU FICHIER xvue/xvuelc.c'
             ELSE
                KERR(1) ='NO PRINTING of FILE ' // KLG(NLD)(NCD:NCF)
     S                   //'.eps IMPOSSIBLE'
                KERR(2) ='WRONG INSTRUCTION OF PRINTING'
                KERR(3) ='CHANGE IT IN PROCEDURE xvimprimerps'
                KERR(4) ='OF FILE xvue/xvuelc.c'
             ENDIF
             CALL LEREUR
           ENDIF
         ELSE
           NBLGRC(NRERR) = 2
           IF( LANGAG .EQ. 0 ) THEN
              KERR(1) = 'IMPRESSION DU FICHIER POSTSCRIPT DE NOM'
              KERR(2) =  KLG(NLD)(NCD:NCF)
           ELSE
              KERR(1) = 'PRINTING OF PS FILE'
              KERR(2) =  KLG(NLD)(NCD:NCF)
           ENDIF
         ENDIF
         GOTO 15
C
      ELSE IF( KLG(NL)(NC:NC+6) .EQ. 'QUITTER' .OR.
     &         KLG(NL)(NC:NC+3) .EQ. 'QUIT'    .OR.
     &         KLG(NL)(NC:NC+3) .EQ. 'EXIT'    .OR.
     &         KLG(NL)(NC:NC+5) .EQ. 'FERMER'  .OR.
     &         KLG(NL)(NC:NC+4) .EQ. 'CLOSE' ) THEN
C        QUITTER OU FERMER LE FICHIER ACTUEL DE LECTURE
C        RETOUR AU LECTEUR PRECEDENT DANS LA PILE
         IF ( LHLECT.LE.0 .OR. LHLECT.GT.MXLECT ) THEN
            NBLGRC(NRERR) = 2
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'AUCUN AUTRE LECTEUR DANS LA PILE'
               KERR(2) = 'RETOUR AU LECTEUR INITIAL'
            ELSE
               KERR(1) = 'NO OTHER INPUT FILE IN THE STACK'
               KERR(2) = 'RETURN TO INITIAL INPUT DEVICE'
            ENDIF
            CALL LEREUR
C           INTERA = 0:BATCH      0 ECRAN GRAPHIQUE 0 CLAVIER 0 SOURIS
C                    1:INTERACTIF 1 ECRAN GRAPHIQUE 0 CLAVIER 0 SOURIS
C                    3:INTERACTIF 1ECRAN GRAPHIQUE  1 CLAVIER 1 SOURIS
            IF( INTERA .LE. 1 ) THEN
C              ERREUR EN BATCH OU LECTURE FICHIER => ARRET DU CALCUL
               CALL ARRET( 100 )
            ENDIF
            LHLECT    = 1
            LECTEU    = IINFO('LECTEUR INITIAL')
            LPLECT(1) = LECTEU
            INTERA    = IINFO('INTERACTIVITE INITIALE')
            INTERB(1) = INTERA
            GOTO 9000
         ENDIF

         IF( KLG(NL)(NC:NC+5) .EQ. 'FERMER' .OR.
     %       KLG(NL)(NC:NC+4) .EQ. 'CLOSE' ) THEN
C           FERMETURE DU FICHIER
            CLOSE( LPLECT(LHLECT) )
         ENDIF
C        QUITTER => PAS DE FERMETURE DU FICHIER POUR ASSURER
C                   UNE EVENTUELLE POURSUITE DE LA LECTURE

C        DESCENTE AU LECTEUR DE HAUTEUR INFERIEURE
         LHLECT = LHLECT - 1
         IF( LHLECT .EQ. 1 ) THEN
            IF( INTERB(LHLECT) .LE. 0 ) THEN
C              RETOUR AU NIVEAU CLAVIER EN MODE BATCH => ARRET
C              EN BATCH ARRET DU CALCUL
               CALL ARRET( 100 )
            ENDIF
         ENDIF

         LECTEU = LPLECT( LHLECT )
C        LE NIVEAU D'INTER-ACTION EST REGENERE
         INTERA = INTERB( LHLECT )
C
C        AFFICHAGE DU RETOUR AU FICHIER PRECEDENT
         IF( INTERA .GE. 1 ) THEN
            CALL RECTEF( NRERR  )
            CALL RECTEF( NRINVI )
            CALL RECTEF( NRMENU )
         ENDIF
C        REMISE A ZERO DES BUFFERS
         NBLGRC(NRMENU) = 0
         NBLGRC(NRLGLU) = 0
C
C        LE RESULTAT
         NBLGRC(NRLGSA) = 0
         NBLGRC(NRERR ) = 2
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'FIN de LECTURE du FICHIER UTILISATEUR'
            KERR(2) = 'Pour CONTINUER frapper @ ou ABANDON'
         ELSE
            KERR(1) = 'END of READING a USER INPUT FILE'
            KERR(2) = 'To CONTINUE TYPE @ or ESCAPE'
         ENDIF
         GOTO 15
C
      ELSE IF( KLG(NL)(NC:NC+4) .EQ. 'LIREF' .OR.
     %         KLG(NL)(NC:NC+4) .EQ. 'READF' ) THEN
C        LIRE A PARTIR DE MAINTENANT LE FICHIER DONT LE NOM SUIT
C        LECTURE DU NOM DU FICHIER
         NC  = NC + 5
         CALL INVITE( 49 )
         CALL LIUTCH( NL , NC , NLD , NCD , NLF , NCF )
         IF( NLF .LE. 0 ) GOTO 9500
C
         CALL MINUSC( KLG(NLD)(NCD:NCF) )
         INQUIRE( FILE=KLG(NLD)(NCD:NCF) , OPENED=LOGIQ , NUMBER=NF ,
     %            IOSTAT=N )
         IF( LOGIQ ) GOTO 88
C        RECHERCHE D'UN NUMERO D'UNITE LIBRE
         CALL TRUNIT( NF )
         OPEN(FILE=KLG(NLD)(NCD:NCF),UNIT=NF,ACCESS='SEQUENTIAL',
     %        FORM='FORMATTED',STATUS='OLD',
     %        IOSTAT= N )
 88      IF( N .NE. 0 ) THEN
            NBLGRC(NRERR) = 3
            KERR(2) = KLG(NLD)(NCD:NCF)
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) ='FICHIER INCONNU ou DEJA OUVERT'
               KERR(3) ='Retaper LIREF NOM_CORRECT_DU_FICHIER ;'
            ELSE
               KERR(1) ='UNKNOWN FILE or ALREADY OPENED'
               KERR(3) ='Type again READF CORRECT_NAME_FILE ;'
            ENDIF
            CALL LEREUR
            IF( INTERA .LE. 1 ) THEN
C              ERREUR EN BATCH ou LECTURE FICHIER L'EXECUTION EST ARRETEE
               CALL ARRET( 100 )
            ENDIF
            NLPTV1 = 0
            GOTO 9000
         ENDIF
         IF (LHLECT.GE.MXLECT) THEN
             NBLGRC(NRERR) = 2
             IF( LANGAG .EQ. 0 ) THEN
                KERR(1) ='PILE SATUREE DES LECTEURS'
                KERR(2) ='Augmenter MXLECT dans ./incl/pilect.inc'
             ELSE
                KERR(1) ='SATURATED STACK of INPUT FILES'
                KERR(2) ='Augment MXLECT in ./incl/pilect.inc'
             ENDIF
             CALL LEREUR
             IF( INTERA .LE. 1 ) THEN
C               ERREUR EN BATCH L'EXECUTION EST ARRETEE
                CALL ARRET( 100 )
             ENDIF
             GOTO 9500
         ENDIF
C        AFFICHAGE DE LA REDIRECTION DE LA LECTURE
         IF( INTERA .GE. 1 ) CALL RECTEF( NRERR )
         WRITE(IMPRIM,*)
         NBLGRC(NRERR) = 2
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'REDIRECTION de la LECTURE sur le FICHIER'
         ELSE
            KERR(1) = 'DATA INPUT is REDIRECTED on the FILE'
         ENDIF
         KERR(2) = KLG(NLD)(NCD:NCF)
         CALL EFFACEMEMPX
         CALL LERESU
         WRITE(IMPRIM,*)
         LECTEU = NF
         LHLECT = LHLECT + 1
         LPLECT(LHLECT) = LECTEU
C        LECTURE SUR FICHIER DE DONNEES
C        => PAS DE VISUALISATION GRAPHIQUE PERMISE
C        INTERA = 0:BATCH      PAS d' ECRAN GRAPHIQUE PAS de CLAVIER PAS de SOURIS
C                 1:INTERACTIF AVEC   ECRAN GRAPHIQUE PAS de CLAVIER PAS de SOURIS
C                 3:INTERACTIF AVEC   ECRAN GRAPHIQUE AVEC   CLAVIER AVEC   SOURIS
         INTERA = MIN( 1, INTERA )
         INTERB(LHLECT) = INTERA
         IF( NOTYPE .EQ. -2 ) THEN
            I = -1
            GOTO 17
         ELSE
            GOTO 15
ccc         GOTO 8000
         ENDIF
C
      ELSE IF( KLG(NL)(NC:NC+6) .EQ. 'DEFFONC' .OR.
     %         KLG(NL)(NC:NC+6) .EQ. 'DEFFUNC' ) THEN
C        DEFINITION D'UNE FONCTION_UTILISATEUR :
C        DEFFONC <IDENT_FONC> ( <PARAM> {,<PARAM>} ) = <VALEUR_D> ;
         CALL LIDFFO( NL , NC , NLF , NCF , NRETOU )
         IF( NRETOU .EQ. -1 ) GOTO 9500
         IF( NRETOU .NE.  0 ) GOTO 9000
         GOTO 8000
C
      ELSE IF( KLG(NL)(NC:NC+5) .EQ. 'DEFVAR'  ) THEN
C        DEFINITION D'UNE VARIABLE_UTILISATEUR : DEFVAR <IDENT_VAR>
         CALL LIDFVA( NL , NC , NLF , NCF , NRETOU )
         IF( NRETOU .NE. 0 ) GOTO 9000
         GOTO 8000
C
      ELSE IF( KLG(NL)(NC:NC) .EQ. '''' ) THEN
C        CHAINE DE CARACTERES ENCADREE PAR DES ' NON DOUBLEES
         CALL LICHAI( NL , NC , NLF , NCF , NCVALS , KCHAIN )
         GOTO 8000
C
      ELSE
C        RECHERCHE D'UNE EVENTUELLE VARIABLE UTILISATEUR SUIVIE
C        DU SIGNE = ET D'UNE AFFECTATION
C        <IDENT_VAR> = <EXPRESS_D>
         NL0 = NL
         NC0 = NC
         CALL CHVARU( NL , NC , NLF , NCF , NOVARU )
         IF( NOVARU .GT. 0 ) THEN
C           L'IDENTIFICATEUR EST IL IMMEDIATEMENT SUIVI DE =
            NL = NLF
            NC = NCF
            CALL CARPNB( NL , NC )
            IF( NL .EQ. 0 ) GOTO 9000
            IF( KLG(NL)(NC:NC) .EQ. '=' ) THEN
C               OUI : TRAITEMENT DE <IDENT_VAR> = <EXPRESS_D>
C               RECHERCHE DU 1-ER CARACTERE NON BLANC DERRIERE =
                CALL CARPNB( NL , NC )
                IF( NL .LE. 0 ) GOTO 9000
                CALL LIEXPD( NL, NC , NLPTV2 , NCPTV2,
     %                       NCVALS , DBLVAL , NRETOU )
                IF( NRETOU .NE. 0 ) GOTO 9000
C               LA POSITION DU DENIER CARACTERE TRAITE
                NLF = NLPTV2
                NCF = NCPTV2
                IF( NCVALS .EQ. 1 ) THEN
C                  LA VARIABLE EST INITIALISEE
                   DVARU( NOVARU ) = DBLVAL
                   I = INDEX(KVARU(NOVARU),' ')
                   IF( I .EQ. 0 ) THEN
                      I = NBCAVA
                   ELSE
                      I = I - 1
                   ENDIF
                   IF( LHLECT .EQ. 1 ) THEN
                      NBLGRC(NRERR) = 1
                      WRITE(KERR(2)(1:15),'(G15.7)') DBLVAL
                      KERR(1) = KVARU(NOVARU)(1:I) //
     %                          '='//KERR(2)(1:15)
                   ENDIF
                ENDIF
                NCVALS = 0
                GOTO 8000
            ELSE
C               LA VARIABLE N'EST PAS SUIVIE DE =
C               REMISE EN L'ETAT AVANT LA DECOUVERTE DE LA VARIABLE
                NL = NL0
                NC = NC0
            ENDIF
         ENDIF
C
         IF( NOTYPS .EQ. 2 ) THEN
C           LECTURE FORCEE D'UNE CHAINE DE CARACTERES NON ENCADREE DE '
            CALL LIUTCH( NL,NC, NLD,NCD, NLF,NCF )
            IF( NLF .EQ. -1 ) GOTO 9500
            IF( NLF .LE.  0 ) GOTO 8000
            CALL MOTRES( KLG(NLD)(NCD:NCF), NONOUI )
            IF( NONOUI .NE. 0 ) GOTO 7000
C           LA CHAINE N'EST PAS UN MOT RESERVE => RETOUR CORRECT
            KCHAIN = KLG(NLD)(NCD:NCF)
            NCVALS = 2
            GOTO 8000
         ENDIF
C
C        ICI CE NE PEUT ETRE QU'UNE <EXPRESS_D>
 7000    CALL LIEXPD( NL, NC , NLPTV2 , NCPTV2,
     %                NCVALS , DBLVAL , NRETOU )
         IF( NRETOU .NE.  0 ) GOTO 9000
         NLF = NLPTV2
         NCF = NCPTV2
C
C        IMPRESSION DU RESULTAT
C        ======================
 8000    IF( LUIMPR .GT. 0 ) THEN
            IF( LHLECT .EQ. 1 ) THEN
               IF( NCVALS .EQ. 1 ) THEN
                  IF( INTERA .GE. 1 ) CALL RECTEF( NRERR )
                  NBLGRC(NRERR) = 1
                  WRITE(KERR(2)(1:15),'(G15.7)') DBLVAL
                  IF( LANGAG .EQ. 0 ) THEN
                     KERR(1) = 'Valeur Lue='//KERR(2)(1:15)
                  ELSE
                     KERR(1) = 'Read Value='//KERR(2)(1:15)
                  ENDIF
               ELSE IF( NCVALS .EQ. 2 ) THEN
                  IF( INTERA .GE. 1 ) CALL RECTEF( NRERR )
                  NBLGRC(NRERR) = 1
                  IF( LANGAG .EQ. 0 ) THEN
                     KERR(1) = 'Chaine lue='// KCHAIN
                  ELSE
                     KERR(1) = 'Read String='// KCHAIN
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
C
C        OU REPRENDRE L'ANALYSE DES DONNEES ?
C        ====================================
         IF( NLF .GT. 0 .AND. (NLF .LT. NLPTV2 .OR.
     %      (NLF .EQ. NLPTV2 .AND. NCF .LT. NCPTV2)) ) THEN
            NL = NLF
            NC = NCF
         ELSE
            NL = NLPTV2
            NC = NCPTV2
         ENDIF
C
C        RECHERCHE DU PROCHAIN ;
         IF( KLG(NL)(NC:NC) .NE. ';' ) CALL CARPNB( NL , NC )
C        PROTECTION SI PAS DE ;
         IF( NL .LE. 0 .OR. NC .LE. 0 ) GOTO 15
         IF( KLG(NL)(NC:NC) .EQ. ';' ) THEN
C
C           LE ; A ETE ATTEINT
C           LES LIGNES DE KLG(NLPTV2:LHKLG) SONT TRANSLATEES EN DEBUT DE KLG
            IF( NLPTV2 .GT. 1 ) THEN
C              LES LIGNES LUES SONT EFFACEES
               CALL RECTEF( NRLGLU )
               DO 8010 I=NLPTV2,LHKLG
                  KLG(I-NLPTV2+1) = KLG(I)
 8010          CONTINUE
               LHKLG  = LHKLG - NLPTV2 + 1
               NLPTV2 = 1
            ENDIF
C
            IF(LHKLG.EQ.1 .AND. KLG(NLPTV2)(NCPTV2+1:NBCALI).EQ.' ')THEN
C              LE RESTE DE LA LIGNE EST BLANC.LE BUFFER EST REMIS A BLANC
               NLPTV2 = 0
               NCPTV2 = 0
            ENDIF
         ELSE
            NBLGRC(NRERR) = 2
            IF( LANGAG .EQ. 0 ) THEN
              KERR(1)=
     %    'DONNMF: INCORRECTE FIN SANS ; FINAL DANS LA LIGNE SUIVANTE:'
            ELSE
              KERR(1)=
     %    'DONNMF: INCORRECT END WITHOUT FINAL ; IN THE FOLLOWING LINE:'
            ENDIF
            KERR(2) = KLG(NL)
            CALL LEREUR
            NCVALS = -1
            NLPTV2 = 0
            NCPTV2 = 0
         ENDIF
C
C        LES VALEURS DE DEPART DANS KLG POUR LA SUITE
         NLPTV1 = NLPTV2
         NCPTV1 = NCPTV2
C
         IF( NCVALS .EQ. 0 ) GOTO 10
         GOTO 9900
      ENDIF
C
C     ERREUR
C     ======
 9000 NCVALS = 0
      NLPTV1 = 0
      IF ( LHLECT .GT. 1 .AND. LHLECT .LE. MXLECT ) THEN
C        LECTURE ACTUELLE SUR UN FICHIER AVEC UNE ERREUR
C        AFFICHAGE DU RETOUR AU CLAVIER
         CALL QUIFIC
         CALL LEREUR
      ENDIF
      GOTO 10
C
C     ABANDON DE LA LECTURE
C     =====================
 9500 NCVALS = -1
      NLPTV1 =  0
      IF( LHLECT .EQ. 1 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'Frappe de @ => ABANDON de la LECTURE'
         ELSE
            KERR(1) = 'Type of @ => ESCAPE of READING'
         ENDIF
      ENDIF
C
C     L'INVITE EST EFFACEE
CCC 9900 CALL RECTEF( NRINVI )
 9900 RETURN
      END
