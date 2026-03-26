      SUBROUTINE LIRLIG( IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   LIRE DANS KLG(LHKLG) UNE LIGNE SUPPLEMENTAIRE DE DONNEES
C -----   SAUVEGARDER SUR LE FICHIER FRAPPE CETTE LIGNE
C
C SORTIE :
C --------
C IERR : 1 SI KLG EST SATURE
C        0 SI LECTURE CORRECTE
C       -1 SI LE CARACTERE '@' APPARTIENT A LA DERNIERE LIGNE LUE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS    JANVIER 1990
C MODIFS : PERRONNET ALAIN LJLL UPMC & St PIERRE DU PERRAY  AVRIL   2013
C23456---------------------------------------------------------------012
      include"./incl/lu.inc"
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/pilect.inc"
      include"./incl/nisafr.inc"
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NFDOCU,NFFRAP,NUNITE(27)

C     EXISTE-T-IL UNE LIGNE DISPONIBLE DANS KLG ?
 1    IF( LHKLG .GT. 0 ) THEN
         IF( NUDCNB( KLG(LHKLG) ) .EQ. 0 ) THEN
            LHKLG = LHKLG - 1
            GOTO 1
         ENDIF
      ENDIF

      LHKLG = LHKLG + 1
      IF( LHKLG .GT. MXKLG ) THEN
          NBLGRC(NRERR) = 2
          IF( LANGAG .EQ. 0 ) THEN
             KERR(1)='TROP DE LIGNES DE DONNEES. A REDUIRE ou'
             KERR(2)='PARENTHESE OUVRANTE { NON TERMINEE PAR  }'
          ELSE
             KERR(1)='TOO INPUT LINES to be REDUCED or'
             KERR(2)='OPENED { NOT FINISHED BY  }'
          ENDIF
          CALL LEREUR
          LHKLG = 0
C         LE CODE DE RETOUR
          IERR  = 1
          RETURN
      ENDIF

C     MISE A BLANC ET LECTURE DE LA LIGNE
 3    KLG(LHKLG) = ' '

C     LECTURE DE LA LIGNE DE DONNEE SELON LE MODE : FICHIER ou CLAVIER ou SOURIS
C     --------------------------------------------------------------------------
 5    IF( INTERA .LE. 1 ) THEN

C         LECTURE SUR UN FICHIER DE DONNEES
C         =================================
cccC        L'INVITE A ENTRER UNE DONNEE
ccc         IF( LHLECT .EQ. 1 ) THEN
ccc            IF( LANGAG .EQ. 0 ) THEN
ccc               WRITE(IMPRIM,10000)
ccc            ELSE
ccc               WRITE(IMPRIM,10001)
ccc            ENDIF
ccc         ENDIF
ccc10000    FORMAT('Frappez la LIGNE de DONNEE: ')
ccc10001    FORMAT('Type the INPUT DATA LINE: ')

C        LECTURE DE LA LIGNE DE DONNEE SUR UN FICHIER DE NO LECTEU
         READ (LECTEU,'(A)',IOSTAT=IERR) KLG(LHKLG)
         IF( IERR .NE. 0 ) THEN
            IF( LANGAG .EQ. 0 ) THEN
               WRITE(IMPRIM,*)
               WRITE(IMPRIM,*) 'FIN du FICHIER de LECTURE ',LECTEU
               WRITE(IMPRIM,*) 'PASSAGE au LECTEUR PRECEDENT'
            ELSE
               WRITE(IMPRIM,*)
               WRITE(IMPRIM,*) 'END of INPUT FILE ',LECTEU
               WRITE(IMPRIM,*) 'BACK to the PREVIOUS INPUT FILE'
            ENDIF
            WRITE(IMPRIM,*)

ccc            IF( INTERA .LE. 0 ) THEN
cccC              EN BATCH L'EXECUTION EST ARRETEE
ccc               CALL ARRET( 100 )
ccc            ENDIF

            IF( LHLECT .LE. 1 .OR. LHLECT .GT. MXLECT ) THEN
               IF( LANGAG .EQ. 0 ) THEN
                  WRITE(IMPRIM,*)'AUCUN LECTEUR dans la PILE'
                  WRITE(IMPRIM,*)'CORRIGER LE FICHIER DE DONNEES'
               ELSE
                  WRITE(IMPRIM,*)'NO INPUT FILE in the STACK'
                  WRITE(IMPRIM,*)'CORRECT THE DATA FILE'
               ENDIF
               WRITE(IMPRIM,*)
               CALL ARRET( 100 )
            ENDIF
            LHLECT = LHLECT - 1
            LECTEU = LPLECT( LHLECT )
C           LE NIVEAU D'INTERACTIVITE
            INTERA = INTERB( LHLECT )
            GOTO 3
         ENDIF

      ELSE

C        ENTREE A LA SOURIS ou AU CLAVIER
C        ================================
         CALL EVMENU

      ENDIF

C     SUPPRESSION DES TABULATIONS
 56   NDC = INDEX( KLG(LHKLG), '\t' )
      IF( NDC .GT. 0 ) THEN
         KLG(LHKLG)(NDC:NDC) = ' '
         GOTO 56
      ENDIF

C     SUPPRESSION DES RETOUR CHARIOT
 57   NDC = INDEX( KLG(LHKLG), '\r' )
      IF( NDC .GT. 0 ) THEN
         KLG(LHKLG)(NDC:NDC) = ' '
         GOTO 57
      ENDIF

C     SUPPRESSION DES BEEP SONORES
 58   NDC = INDEX( KLG(LHKLG), '\a' )
      IF( NDC .GT. 0 ) THEN
         KLG(LHKLG)(NDC:NDC) = ' '
         GOTO 58
      ENDIF

C     SUPPRESSION DES SAUTS A LA LIGNE
 59   NDC = INDEX( KLG(LHKLG), '\n' )
      IF( NDC .GT. 0 ) THEN
         KLG(LHKLG)(NDC:NDC) = ' '
         GOTO 59
      ENDIF

C     SUPPRESSION DES SAUTS A LA PAGE SUIVANTE
 60   NDC = INDEX( KLG(LHKLG), '\f' )
      IF( NDC .GT. 0 ) THEN
         KLG(LHKLG)(NDC:NDC) = ' '
         GOTO 60
      ENDIF

C     SUPPRESSION DES TABULATIONS VERTICALES
 61   NDC = INDEX( KLG(LHKLG), '\v' )
      IF( NDC .GT. 0 ) THEN
         KLG(LHKLG)(NDC:NDC) = ' '
         GOTO 61
      ENDIF

C     SUPPRESSION DES SAUTS DES BACKSPACES
 62   NDC = INDEX( KLG(LHKLG), '\b' )
      IF( NDC .GT. 0 ) THEN
         KLG(LHKLG)(NDC:NDC) = ' '
         GOTO 62
      ENDIF

C     SI LA LIGNE EST BLANCHE RETOUR A LA DEMANDE DE DONNEES
      NDC = NUDCNB( KLG(LHKLG) )
      IF( NDC .LE. 0 ) GOTO 5
C     AFFICHAGE DE LA LIGNE LUE
      IF( LANGAG .EQ. 0 ) THEN
         IF( LHLECT .EQ. 1 ) THEN
            WRITE(IMPRIM,10057) KLG(LHKLG)
ccc         ELSE
ccc            WRITE(IMPRIM,10058) LHLECT,KLG(LHKLG)  26/11/2015
         ENDIF
      ELSE
         IF( LHLECT .EQ. 1 ) THEN
            WRITE(IMPRIM,20057) KLG(LHKLG)
ccc         ELSE
ccc            WRITE(IMPRIM,20058) LHLECT,KLG(LHKLG)  26/11/2015
         ENDIF
      ENDIF
10057 FORMAT(' Clavier:',A)
ccc10058 FORMAT(' FichierLu',I3,':',A)
20057 FORMAT(' Keyboard:',A)
ccc20058 FORMAT(' Read File',I3,':',A)

C     SAUVEGARDE SUR LE FICHIER FRAPPE DE LA LIGNE DE DONNEES
      IF( LHLECT .EQ. 1 .AND. NISAFR .NE. 0 ) THEN
         CALL SANSDBL( KLG(LHKLG), IERR )
         WRITE(NFFRAP,*) KLG(LHKLG)(1:IERR)
      ENDIF

C     TRANSFORMATION DES LETTRES MINUSCULES EN MAJUSCULES
      CALL MAJUSC( KLG(LHKLG) )

C     LE CARACTERE D'ABANDON @ EXISTE T IL DANS CETTE LIGNE?
      IF( INDEX( KLG(LHKLG)(1:NBCALI), '@' ) .GT. 0 ) THEN
C        OUI
         IERR = -1
      ELSE
C        NON
         IERR = 0
      ENDIF

      RETURN
      END
