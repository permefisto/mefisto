      SUBROUTINE TEXTSB( TEXT, NUDCNONB )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    SUPPRIMER LES BLANCS DE DEBUT D'UN TEXTE ET INTERMEDIAIRES
C -----    RETOURNER LE NUMERO DU DERNIER CARACTERE NON BLANC DU TEXTE
C
C MODIFIE:
C --------
C TEXT   : LA CHAINE DE CARACTERES
C
C SORTIE :
C --------
C NUDCNONB: NUMERO DU DERNIER CARACTERE NON BLANC APRES SUPPRESSION
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC &St PIERRE DU PERRAY Septembre 2012
C23456---------------------------------------------------------------012
      CHARACTER*(*)  TEXT
      CHARACTER*128  TEXTE
C
C     NOMBRE DE CARACTERES DE TEXT
      NBC = LEN( TEXT )
C
      TEXTE  = TEXT
      NBLANC = 0
      DO L=1,NBC
C        NBLANC NOMBRE DE BLANCS EN DEBUT DE TEXT
         IF( TEXTE(L:L) .EQ. ' ' ) THEN
            NBLANC = NBLANC+1
         ELSE
            GOTO 10
         ENDIF
      ENDDO
C
C     SI DES BLANCS ALORS TRANSLATION
 10   IF( NBLANC .GT. 0 ) THEN
         DO L=NBLANC+1,NBC
            TEXTE(L-NBLANC:L-NBLANC) = TEXTE(L:L)
         ENDDO
C        AJOUT DES NBLANC BLANCS A LA FIN DU TEXTE
         DO L=NBC-NBLANC+1,NBC
            TEXTE(L-NBLANC:L-NBLANC) = ' '
         ENDDO
      ENDIF
C
C     NUMERO DU DERNIER CARACTERE NON BLANC APRES SUPPRESSION
      NUDCNONB = NUDCNB( TEXTE )
      NB = 1
C
C     SUPPRESSION DES BLANCS INTERMEDIAIRES
 20   CONTINUE
      DO L=NB,NUDCNONB
         IF( TEXTE(L:L) .EQ. ' ' ) THEN
            DO LL=L+1,NUDCNONB
               TEXTE(LL-1:LL-1) = TEXTE(LL:LL)
            ENDDO
            TEXTE(NUDCNONB:NUDCNONB) = ' '
            NUDCNONB = NUDCNONB - 1
            NB = L
            GOTO 20
         ENDIF
      ENDDO
C
      TEXT(1:NUDCNONB) = TEXTE(1:NUDCNONB)
C
      RETURN
      END
