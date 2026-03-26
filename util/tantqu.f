      SUBROUTINE TANTQU( NLMIN  , NCMIN  , NLMAX  , NCMAX ,  NRETOU )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : REMPLACER LES INSTRUCTIONS TANTQUE EXP_LOGIQ ; ... FINTANT ;
C ----- PAR LEUR EQUIVALENCE AVEC SI ALORS SINON
C
C       TANTQUE EXP_LOGIQ ;    ou    WHILE EXP_LOGIQ ;
C          INSTRUCTION                  INSTRUCTION
C       FINTANT ;                    ENDWHILE ;
C  <=>
C       ETIQUETTE:
C       SI  EXP_LOGIQ  ALORS
C          INSTRUCTION
C          ALLER ETIQUETTE ;
C       FINSI;
C
C CF $MEFISTO/td/g/grammaire_lu DE DEFINITION DU LANGAGE UTILISATEUR
C
C ENTREES :
C ---------
C NLMIN,NCMIN : POSITION DANS KLG DU PREMIER CARACTERE A TRAITER
C
C MODIFIES :
C ----------
C NLMAX,NCMAX : POSITION DANS KLG DU DERNIER CARACTERE A TRAITER
C
C SORTIES :
C ---------
C NRETOU :  0 SI PAS D'ERREUR
C           9 SI TABLE DES INSTRUCTIONS SATUREE
C          >0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS SEPTEMBRE 1999
C23456---------------------------------------------------------------012
      include"./incl/lu.inc"
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      CHARACTER*7       ETIQUE
      CHARACTER*4       KFORMA
C
C     RECHERCHE DU PREMIER 'TANTQUE' OU 'WHILE'
 1    NC1 = NCMIN
      NC2 = NBCALI
      DO 20 NL1=NLMIN,NLMAX
         IF( NL1 .NE. NLMIN ) NC1=1
         IF( NL1 .EQ. NLMAX ) NC2=NCMAX
         NC0 = INDEX( KLG(NL1)(NC1:NC2), 'TANTQUE ' )
         NCC = INDEX( KLG(NL1)(NC1:NC2), 'WHILE ' )
         IF( NC0 .LE. 0 .AND. NCC .LE. 0 ) GOTO 20
C        UN DES 2 MOTS AU MOINS A ETE RETROUVE
         IF( NC0 .GT. 0 ) THEN
C           TANTQUE RETROUVE
            IF( NCC .GT. 0 .AND. NCC .LT. NC0 ) THEN
C              WHILE EST AVANT TANTQUE
               GOTO 12
            ELSE
C              TANTQUE EST AVANT
               GOTO 11
            ENDIF
         ELSE
C           PAS DE TANTQUE
            GOTO 12
         ENDIF
C
C        TANTQUE EST LE PREMIER
 11      NBC1 = 7
         GOTO 21
C
C        WHILE EST LE PREMIER
 12      NC0  = NCC
         NBC1 = 5
C        2 BLANCS AJOUTES DERRIERE WHILE
         CALL DEPKLG( NL1, NC0+NC1+4, NLMAX, NCMAX, 2 )
         GOTO 21
 20   CONTINUE
C
C     PLUS DE TANTQUE
      GOTO 8000
C
C     TANTQUE OU WHILE DEBUTE EN NL1,NC0
 21   NC0 = NC1 - 1 + NC0
C
C     RECHERCHE DU DERNIER FINTANT OU ENDWHILE
      NC1 = 1
      NC2 = NCMAX
      DO 22 NL2=NLMAX,NL1,-1
         IF( NL2 .EQ. NL1   ) NC1=NC0+NBC1
         IF( NL2 .NE. NLMAX ) NC2=NBCALI
         NC3 = INDEX( KLG(NL2)(NC1:NC2) , 'FINTANT' )
         IF( NC3 .GT. 0 ) THEN
C           EXISTE T IL UN AUTRE FINTANT SUR CETTE LIGNE?
            NC3 = NC1 - 1 + NC3
 15         NQ  = INDEX( KLG(NL2)(NC3+1:NC2) , 'FINTANT' )
            IF( NQ .GT. 0 ) THEN
               NC3 = NC3 + NQ
               GOTO 15
            ENDIF
         ENDIF
         NCC = INDEX( KLG(NL2)(NC1:NC2) , 'ENDWHILE' )
         IF( NCC .GT. 0 ) THEN
C           EXISTE T IL UN AUTRE ENDWHILE SUR CETTE LIGNE?
            NCC = NC1 - 1 + NCC
 17         NQ  = INDEX( KLG(NL2)(NCC+1:NC2) , 'ENDWHILE' )
            IF( NQ .GT. 0 ) THEN
               NCC = NCC + NQ
               GOTO 17
            ENDIF
         ENDIF
         IF( NC3 .LE. 0 .AND. NCC .LE. 0 ) GOTO 22
C
C        RECHERCHE DU DERNIER DES 2
         IF( NC3 .LT. NCC ) THEN
C           ENDWHILE AU DELA DE FINTANT
            NC3  = NCC
            NBC2 = 8
C           E FINAL DE ENDWHILE DEVIENT UN BLANC
            KLG(NL2)(NC3+7:NC3+7) = ' '
         ELSE
C           FINTANT AU DELA DE ENDWHILE
            NBC2 = 7
         ENDIF
         GOTO 23
 22   CONTINUE
C
      NBLGRC(NRERR) = 1
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = 'LU: TANTQUE sans ''FINTANT'' '
      ELSE
         KERR(1) = 'LU: WHILE without ''ENDWHILE'' '
      ENDIF
      CALL LEREUR
      GOTO 9900
 23   NC3 = NC1 - 1 + NC3
C
C     REMPLACEMENT DE TANTQUE PAR UNE ETIQUETTE: SI
      IF( NBETIQ .GE. MXETIQ ) GOTO 9500
C     AJOUT DE L'ETIQUETTE
      NBETIQ = NBETIQ + 1
C     LE NOMBRE DE CHIFFRES DE NBETIQ
      NC1    = NBCHIF( NBETIQ )
      ETIQUE = 'TWKZYO0'
      KFORMA = '(I0)'
      WRITE(KFORMA(3:3),'(I1)') NC1
      WRITE(ETIQUE(5-NC1:4),KFORMA) NBETIQ
      KETIQ(NBETIQ) = ETIQUE(1:4)
CCC      NBLGRC(NRERR) = 1
CCC      KERR(1) = 'TANTQU: AJOUT DE L''ETIQUETTE ' // KETIQ(NBETIQ)
CCC      CALL LERESU
C
C     TANTQUE OU WHILE EST ECRASE PAR CETTE ETIQUETTE
      ETIQUE(5:7) = ':SI'
      KLG(NL1)(NC0:NC0+6) = ETIQUE
C
C     SI ; DERRIERE IL EST SUPPRIME
      NC0 = NC0 + 7
      CALL CARPNB( NL1, NC0 )
      IF( KLG(NL1)(NC0:NC0) .EQ. ';' ) KLG(NL1)(NC0:NC0)=' '
C
C     RECHERCHE DU ; DERRIERE L'EXPRESSION LOGIQUE
      NC1 = NBCALI
      DO 25 NL = NL1 , NL2
         IF( NL .NE. NL1 ) NC0 = 1
         IF( NL .EQ. NL2 ) NC1 = NC3
         NC2 = INDEX( KLG(NL)(NC0:NC1) , ';' )
         IF( NC2 .GT. 0 ) GOTO 27
 25   CONTINUE
      NBLGRC(NRERR) = 1
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = 'LU: EXPRESSION LOGIQUE DERRIERE TANTQUE SANS '';'' '
      ELSE
         KERR(1) = 'LU: NO '';'' after the LOGICAL EXPRESSION of WHILE'
      ENDIF
      CALL LEREUR
      GOTO 9900
C
C     TRANSLATION DES CARACTERES POUR AJOUTER ALORS
 27   NC0 = NC0 - 1 + NC2
C     LE ; DISPARAIT
      KLG(NL)(NC0:NC0) = ' '
      NC0 = NC0 + 1
      CALL DEPKLG( NL , NC0 , NLMAX , NCMAX , 6 )
      IF( NL .LE. 0 ) GOTO 9900
      KLG(NL)(NC0:NC0+5) = 'ALORS '
C
C     RECHERCHE DU DERNIER FINTANT OU ENDWHILE QUI A PU ETRE TRANSLATE
      NC1 = 1
      NC2 = NCMAX
      DO 42 NL2=NLMAX,NL,-1
         IF( NL2 .EQ. NL    ) NC1=NC0+6
         IF( NL2 .NE. NLMAX ) NC2=NBCALI
         NC3 = INDEX( KLG(NL2)(NC1:NC2) , 'FINTANT' )
         IF( NC3 .GT. 0 ) THEN
C           EXISTE T IL UN AUTRE FINTANT SUR CETTE LIGNE?
            NC3 = NC1 - 1 + NC3
 45         NQ  = INDEX( KLG(NL2)(NC3+1:NC2) , 'FINTANT' )
            IF( NQ .GT. 0 ) THEN
               NC3 = NC3 + NQ
               GOTO 45
            ENDIF
         ENDIF
         NCC = INDEX( KLG(NL2)(NC1:NC2) , 'ENDWHIL' )
         IF( NCC .GT. 0 ) THEN
C           EXISTE T IL UN AUTRE ENDWHILE SUR CETTE LIGNE?
            NCC = NC1 - 1 + NCC
 47         NQ  = INDEX( KLG(NL2)(NCC+1:NC2) , 'ENDWHIL' )
            IF( NQ .GT. 0 ) THEN
               NCC = NCC + NQ
               GOTO 47
            ENDIF
         ENDIF
         IF( NC3 .LE. 0 .AND. NCC .LE. 0 ) GOTO 42
C
C        RECHERCHE DU DERNIER DES 2
         IF( NC3 .LT. NCC ) THEN
C           ENDWHILE AU DELA DE FINTANT
            NC3  = NCC
         ENDIF
         GOTO 50
 42   CONTINUE
C
C     REMPLACEMENT DE FINTANT PAR 'FINSI'
 50   KLG(NL2)(NC3:NC3+6) = ' FINSI '
C
C     DEPLACER LA SUITE DE FACON A LIBERER DE LA PLACE DANS KLG
C     POUR AJOUTER ' ALORS ALLER ETIQUETTE ;'
      NC0 = NC3
      CALL DEPKLG( NL2 , NC0 , NLMAX , NCMAX , 11 )
      IF( NL2 .LE. 0 ) GOTO 9900
      KLG(NL2)(NC0:NC0+10) = 'ALLER ' // ETIQUE(1:4) // ';'
C
C     EXISTE T IL D'AUTRES TANTQUE ...
      GOTO 1
C
C     FIN DE TRAITEMENT
 8000 NRETOU = 0
      RETURN
C
 9500 NBLGRC(NRERR) = 1
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = 'LU: PILE SATUREE DES ETIQUETTES'
      ELSE
         KERR(1) = 'LU: SATURATED STACK of LABELS'
      ENDIF
      CALL LEREUR
C
 9900 NRETOU = 1
      RETURN
      END
