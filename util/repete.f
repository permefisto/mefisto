      SUBROUTINE REPETE( NLMIN  , NCMIN  , NLMAX  , NCMAX ,  NRETOU )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : REMPLACER LES INSTRUCTIONS REPETER ... JUSQUE EXP_LOGIQ ; PAR
C ----- LEUR EQUIVALENCE AVEC SI ALORS SINON
C
C       REPETER ;              ou     REPEAT
C          INSTRUCTION                   INSTRUCTION
C       JUSQUE  EXP_LOGIQ ;           UNTIL EXP_LOGIQ ;
C  <=>
C       ETIQUETTE:
C          INSTRUCTION
C       SI NON EXP_LOGIQ ALORS ALLER ETIQUETTE ;
C       FINSI ;
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
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS    FEVRIER 1990
C23456---------------------------------------------------------------012
      include"./incl/lu.inc"
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      CHARACTER*7       ETIQUE
      CHARACTER*4       KFORMA
C
C     RECHERCHE DU PREMIER 'REPETER' OU 'REPEAT'
 1    NC1 = NCMIN
      NC2 = NBCALI
      DO 20 NL1=NLMIN,NLMAX
         IF( NL1 .NE. NLMIN ) NC1=1
         IF( NL1 .EQ. NLMAX ) NC2=NCMAX
         NC0 = INDEX( KLG(NL1)(NC1:NC2) , 'REPETER' )
         NCC = INDEX( KLG(NL1)(NC1:NC2) , 'REPEAT' )
         IF( NC0 .LE. 0 .AND. NCC .LE. 0 ) GOTO 20
C        UN DES 2 MOTS AU MOINS A ETE RETROUVE
         IF( NC0 .GT. 0 ) THEN
C           REPETER RETROUVE
            IF( NCC .GT. 0 .AND. NCC .LT. NC0 ) THEN
C              REPEAT EST AVANT REPETER
               GOTO 12
            ELSE
C              REPETER EST AVANT
               GOTO 11
            ENDIF
         ELSE
C           PAS DE REPETER
            GOTO 12
         ENDIF
C
C        REPETER EST LE PREMIER
 11      NBC1 = 7
         GOTO 21
C
C        REPEAT EST LE PREMIER
 12      NC0  = NCC
         NBC1 = 6
         GOTO 21
 20   CONTINUE
C
C     PLUS DE REPETER
      GOTO 8000
C
C     REPETER OU REPEAT DEBUTE EN NL1,NC0
 21   NC0 = NC1 - 1 + NC0
C
C     RECHERCHE DU DERNIER JUSQUE OU UNTIL
      NC1 = 1
      NC2 = NCMAX
      DO 22 NL2=NLMAX,NL1,-1
         IF( NL2 .EQ. NL1   ) NC1=NC0+NBC1
         IF( NL2 .NE. NLMAX ) NC2=NBCALI
         NC3 = INDEX( KLG(NL2)(NC1:NC2) , 'JUSQUE ' )
         IF( NC3 .GT. 0 ) THEN
C           EXISTE T IL UN AUTRE JUSQUE SUR CETTE LIGNE?
            NC3 = NC1 - 1 + NC3
 15         NQ  = INDEX( KLG(NL2)(NC3+1:NC2) , 'JUSQUE ' )
            IF( NQ .GT. 0 ) THEN
               NC3 = NC3 + NQ
               GOTO 15
            ENDIF
         ENDIF
         NCC = INDEX( KLG(NL2)(NC1:NC2) , 'UNTIL ' )
         IF( NCC .GT. 0 ) THEN
C           EXISTE T IL UN AUTRE UNTIL SUR CETTE LIGNE?
            NCC = NC1 - 1 + NCC
 17         NQ  = INDEX( KLG(NL2)(NCC+1:NC2) , 'UNTIL ' )
            IF( NQ .GT. 0 ) THEN
               NCC = NCC + NQ
               GOTO 17
            ENDIF
         ENDIF
         IF( NC3 .LE. 0 .AND. NCC .LE. 0 ) GOTO 22
C
C        RECHERCHE DU DERNIER DES 2
         IF( NC3 .LT. NCC ) THEN
C           UNTIL AU DELA DE JUSQUE
            NC3  = NCC
            NBC2 = 6
         ELSE
C           JUSQUE AU DELA DE UNTIL
            NBC2 = 7
         ENDIF
         GOTO 23
 22   CONTINUE
      NBLGRC(NRERR) = 1
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = 'LU: REPETER sans ''JUSQUE '' '
      ELSE
         KERR(1) = 'LU: REPEAT without ''UNTIL '' '
      ENDIF
      CALL LEREUR
      GOTO 9900
C
C     REMPLACEMENT DE REPETER OU REPEAT PAR UNE ETIQUETTE
 23   IF( NBETIQ .GE. MXETIQ ) GOTO 9500
C     AJOUT DE L'ETIQUETTE
      NBETIQ = NBETIQ + 1
C     LE NOMBRE DE CHIFFRES DE NBETIQ
      NC1    = NBCHIF( NBETIQ )
      ETIQUE = 'RO0O0O'
      KFORMA = '(I0)'
      WRITE(KFORMA(3:3),'(I1)') NC1
      WRITE(ETIQUE(6-NC1:5),KFORMA) NBETIQ
      ETIQUE(6:6) = ':'
      KETIQ(NBETIQ) = ETIQUE(1:5)
C
C     REPETER OU REPEAT EST ECRASE PAR CETTE ETIQUETTE
      KLG(NL1)(NC0:NC0+5) = ETIQUE(1:6)
C     LE R FINAL DE REPETER EST REMPLACE PAR UN BLANC
      IF( NBC1 .EQ. 7 ) KLG(NL1)(NC0+6:NC0+6) = ' '
C
C     SI ; DERRIERE IL EST SUPPRIME
      NC0 = NC0 + 5
      CALL CARPNB( NL1, NC0 )
      IF( KLG(NL1)(NC0:NC0) .EQ. ';' ) KLG(NL1)(NC0:NC0)=' '
C
C     REMPLACEMENT DE 'JUSQUE ' OU 'UNTIL ' PAR 'SI NON '
      IF( NBC2 .EQ. 6 ) CALL DEPKLG( NL2, NC3+6, NL2, NBCALI, 1 )
      KLG(NL2)(NC3:NC3+6) = 'SI NON '
C
C     RECHERCHE DU ; APRES L'EXPR_LOGIQ
      NC1 = NC3+6
      NC2 = NBCALI
 30   NC0 = INDEX( KLG(NL2)(NC1:NC2) , ';' )
      IF( NC0 .LE. 0 ) THEN
C        ; NON RETROUVE
         NL2 = NL2 + 1
         IF( NL2 .GT. LHKLG ) THEN
            NBLGRC(NRERR) = 1
            IF( LANGAG .EQ. 0 ) THEN
            KERR(1)='LU: '';'' OUBLIEE DERRIERE L''EXPRESSION LOGIQUE'
            ELSE
            KERR(1)='LU: FORGOTTEN '';'' AFTER THE LOGICAL EXPRESSION'
            ENDIF
            CALL LEREUR
            GOTO 9900
         ENDIF
         NC1 = 1
         GOTO 30
      ENDIF
      NC0 = NC1 - 1 + NC0
C
C     DEPLACER LA SUITE DE FACON A LIBERER DE LA PLACE DANS KLG
C     POUR AJOUTER ' ALORS ALLER ETIQUETTE ;'
      KLG(NL2)(NC0:NC0) = ' '
      CALL DEPKLG( NL2 , NC0 , NLMAX , NCMAX , 28 )
      IF( NL2 .LE. 0 ) GOTO 9900
      KLG(NL2)(NC0:NC0+28)=' ALORS ALLER ' // ETIQUE(1:5) // '; FINSI;'
C
C     EXISTE T IL UN AUTRE REPETER OU REPEAT?
      GOTO 1
C
 8000 NRETOU = 0
      RETURN
C
 9500 NBLGRC(NRERR) = 1
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = 'LU: PILE SATUREE DES ETIQUETTES'
      ELSE
         KERR(1) = 'LU: SATURATED STACK OF LABELS'
      ENDIF
      CALL LEREUR
C
 9900 NRETOU = 1
      RETURN
      END
