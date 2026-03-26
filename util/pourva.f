      SUBROUTINE POURVA( NLMIN  , NCMIN  , NLMAX  , NCMAX ,  NRETOU )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : REMPLACER LES INSTRUCTIONS POUR ... FINPOUR ; PAR
C ----- LEUR EQUIVALENCE AVEC SI ALORS SINON
C
C       POUR NOM_VAR = EXP_ARITH1  A  EXP_ARITH2 ;  (; TRES IMPORTANT)
C          INSTRUCTION
C       FINPOUR ;
C
C  ou
C       FOR NOM_VAR = EXP_ARITH1  TO  EXP_ARITH2 ;  (; TRES IMPORTANT)
C          INSTRUCTION
C       ENDFOR ;
C  <=>
C       NOM_VAR = EXP_ARITH1 ;
C       ETIQUETTE: SI NOM_VAR <= EXP_ARITH2 ALORS
C          INSTRUCTION
C          NOM_VAR = NOM_VAR + 1 ;
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
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS  SEPTEMBRE 1999
C23456---------------------------------------------------------------012
      include"./incl/lu.inc"
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      CHARACTER*7       ETIQUE
      CHARACTER*4       KFORMA
C
C     RECHERCHE DU PREMIER 'POUR' OU 'FOR'
 1    NC1 = NCMIN
      NC2 = NBCALI
      DO 20 NL1=NLMIN,NLMAX
         IF( NL1 .NE. NLMIN ) NC1=1
         IF( NL1 .EQ. NLMAX ) NC2=NCMAX
         NC0 = INDEX( KLG(NL1)(NC1:NC2), 'POUR ' )
         NCC = INDEX( KLG(NL1)(NC1:NC2), 'FOR ' )
         IF( NC0 .LE. 0 .AND. NCC .LE. 0 ) GOTO 20
C        UN DES 2 MOTS AU MOINS A ETE RETROUVE
         IF( NC0 .GT. 0 ) THEN
C           POUR RETROUVE
            IF( NCC .GT. 0 .AND. NCC .LT. NC0 ) THEN
C              FOR EST AVANT POUR
               GOTO 12
            ELSE
C              POUR EST AVANT
               GOTO 11
            ENDIF
         ELSE
C           PAS DE POUR
            GOTO 12
         ENDIF
C
C        POUR EST LE PREMIER
 11      NBC1 = 4
         GOTO 21
C
C        FOR EST LE PREMIER
 12      NC0  = NCC
         NBC1 = 3
C        1 BLANC AJOUTE DERRIERE FOR
         CALL DEPKLG( NL1, NC0+NC1+2, NLMAX, NCMAX, 1 )
         GOTO 21
 20   CONTINUE
C
C     PLUS DE POUR OU FOR
      GOTO 8000
C
C     POUR OU WHILE DEBUTE EN NL1,NC0
 21   NC0 = NC1 - 1 + NC0
C     AU DELA DE POUR ou FOR
      NC00 = NC0 + 4
C     POUR ou FOR EST EFFACE
      KLG(NL1)(NC0:NC00) = '    '
      NL0 = NL1
C
C     RECHERCHE DU DERNIER FINPOUR OU ENDFOR
      NC1 = 1
      NC2 = NCMAX
      DO 22 NL2=NLMAX,NL1,-1
         IF( NL2 .EQ. NL1   ) NC1=NC0+NBC1
         IF( NL2 .NE. NLMAX ) NC2=NBCALI
         NC3 = INDEX( KLG(NL2)(NC1:NC2) , 'FINPOUR' )
         IF( NC3 .GT. 0 ) THEN
C           EXISTE T IL UN AUTRE FINPOUR SUR CETTE LIGNE?
            NC3 = NC1 - 1 + NC3
            NL3 = NL2
 15         NQ  = INDEX( KLG(NL2)(NC3+1:NC2) , 'FINPOUR' )
            IF( NQ .GT. 0 ) THEN
               NC3 = NC3 + NQ
               GOTO 15
            ENDIF
         ENDIF
         NCC = INDEX( KLG(NL2)(NC1:NC2) , 'ENDFOR' )
         IF( NCC .GT. 0 ) THEN
C           EXISTE T IL UN AUTRE ENDFOR SUR CETTE LIGNE?
            NCC = NC1 - 1 + NCC
            NL3 = NL2
 17         NQ  = INDEX( KLG(NL2)(NCC+1:NC2) , 'ENDFOR' )
            IF( NQ .GT. 0 ) THEN
               NCC = NCC + NQ
               GOTO 17
            ENDIF
         ENDIF
         IF( NC3 .LE. 0 .AND. NCC .LE. 0 ) GOTO 22
C
C        RECHERCHE DU DERNIER DES 2
         IF( NC3 .LT. NCC ) THEN
C           ENDFOR AU DELA DE FINPOUR
            NC3  = NCC
            NBC2 = 6
C           1 BLANC AJOUTE DERRIERE FOR
            CALL DEPKLG( NL2, NC1+NC3+5, NLMAX, NCMAX, 1 )
         ELSE
C           FINPOUR AU DELA DE ENDFOR
            NBC2 = 7
         ENDIF
         GOTO 23
 22   CONTINUE
C
      NBLGRC(NRERR) = 1
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = 'LU: POUR sans ''FINPOUR'' '
      ELSE
         KERR(1) = 'LU: FOR without ''ENDFOR'' '
      ENDIF
      CALL LEREUR
      GOTO 9900
C
 23   NC3 = NC1 - 1 + NC3
C
C     RECHERCHE DE = DERRIERE POUR ou FOR
      NC0 = NC00
      NC1 = NBCALI
      DO 25 NL = NL1, NL3
         IF( NL .NE. NL1 ) NC0 = 1
         IF( NL .EQ. NL3 ) NC1 = NC3
         NC2 = INDEX( KLG(NL)(NC0:NC1) , '=' )
         IF( NC2 .GT. 0 ) GOTO 27
 25   CONTINUE
C
      NBLGRC(NRERR) = 1
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = 'LU: POUR sans ''='' DERRIERE la VARIABLE'
      ELSE
         KERR(1) = 'LU: FOR without ''='' after the VARIABLE'
      ENDIF
      CALL LEREUR
      GOTO 9900
C
C     = EST EN NL,NC0
 27   NC0 = NC0 - 1 + NC2
C
C     RECHERCHE DE NOM_VAR ENTRE NL1,NC00 ET NL,NC0
      NC1 = NC00 - 1
      CALL CARPNB( NL1 , NC1 )
      NC2 = NC1
      NL2 = NL1
      CALL CARPBL( NL2 , NC2 )
      IF( NL2 .GT. NL1 ) NC2=NBCALI+1
C     PEUT ETRE PAS DE BLANC ENTRE NOM_VAR ET =
      IF( NC2 .GT. NC0 ) THEN
C        IL EXISTE = COLLE AU NOM DE LA VARIABLE
         NL2 = NL1
         NC2 = NC0
      ENDIF
C     LE NOMBRE DE CARACTERES  DE NOM_VARIABLE
      NBCAR = NC2 - NC1
      NC2   = NC2 - 1
C     NOM_VARIABLE SE TROUVE ENTRE NC1 ET NC2 DE LA LIGNE NL1
C
C     RECHERCHE DE 'A' OU 'TO' ENTRE EXP_ARITH1 ET EXP_ARITH2
      NC4 = NC0
      NC5 = NBCALI
      DO 35 NL = NL2 , NL3
         IF( NL .NE. NL2 ) NC4 = 1
         IF( NL .EQ. NL3 ) NC5 = NC3
         NC  = INDEX( KLG(NL)(NC4:NC5) , ' A ' )
         NCC = INDEX( KLG(NL)(NC4:NC5) , ' TO ' )
         NC  = MAX( NC, NCC )
         IF( NC .GT. 0 ) GOTO 37
 35   CONTINUE
C
      NBLGRC(NRERR) = 1
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1)='LU: POUR SANS ''a'' ENTRE EXPRESSIONS ARITHMETIQUES'
      ELSE
         KERR(1)='LU: FOR without ''to'' between ARITHMETIC EXPRESSIONS'
      ENDIF
      CALL LEREUR
      GOTO 9900
C
C     'A ' OU 'TO' EST REMPLACE PAR '; '
 37   NC = NC4 + NC
      KLG(NL)(NC:NC+1) = '; '
C
C     CREATION D'UNE ETIQUETTE
      IF( NBETIQ .GE. MXETIQ ) GOTO 9500
C     AJOUT DE L'ETIQUETTE
      NBETIQ = NBETIQ + 1
C     LE NOMBRE DE CHIFFRES DE NBETIQ
      NC4    = NBCHIF( NBETIQ )
      ETIQUE = 'PWKZ0O0'
      KFORMA = '(I0)'
      WRITE(KFORMA(3:3),'(I1)') NC4
      WRITE(ETIQUE(7-NC4:6),KFORMA) NBETIQ
      ETIQUE(7:7) = ':'
      KETIQ(NBETIQ) = ETIQUE(1:6)
C
CCC      NBLGRC(NRERR) = 1
CCC      KERR(1) = 'POUR: AJOUT DE L''ETIQUETTE '// KETIQ(NBETIQ)
CCC      CALL LERESU
C
C     TRANSLATION DES CARACTERES POUR AJOUTER
C     ETIQ:SI NOM_VAR <= EXP_ARITH2 ALORS
      NC4 = NBCAR + 15
      NC0 = NC + 2
      NL0 = NL
      CALL DEPKLG( NL0 , NC0 , NLMAX , NCMAX , NC4 )
      IF( NL0 .LE. 0 ) GOTO 9900
      KLG(NL0)(NC0:NC0+NC4-1) = ETIQUE // ' SI ' //
     %                          KLG(NL1)(NC1:NC2) // ' <= '
C
C     RECHERCHE DU ; DERRIERE L'EXPRESSION_ARITHMETIQUE2 FINALE
      NC0 = NC0 + NC4
      NC5 = NBCALI
      DO 45 NL = NL0 , NLMAX
         IF( NL .NE. NL0 ) NC0 = 1
         IF( NL .EQ. NLMAX ) NC5 = NCMAX
         NC = INDEX( KLG(NL)(NC0:NC5) , ';' )
         IF( NC .GT. 0 ) GOTO 47
 45   CONTINUE
C
      NBLGRC(NRERR) = 1
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1)='LU: POUR SANS '';'' DERRIERE EXP_ARITH FINALE'
      ELSE
         KERR(1)='LU: FOR without '';'' AFTER the FINAL ARITH_EXP'
      ENDIF
      CALL LEREUR
      GOTO 9900
C
 47   NC0 = NC0 - 1 + NC
C     ; EST EN NL,NC0
C
C     TRANSLATION DES CARACTERES POUR AJOUTER  ALORS
      KLG(NL)(NC0:NC0) = ' '
      CALL DEPKLG( NL , NC0 , NLMAX , NCMAX , 7 )
      IF( NL .LE. 0 ) GOTO 9900
      KLG(NL)(NC0:NC0+6) = ' ALORS '
C
C     RECHERCHE DU DERNIER FINPOUR OU ENDFOR QUI A PU ETRE TRANSLATE
      NC4 = 1
      NC5 = NCMAX
      DO 60 NL3=NLMAX,NL,-1
         IF( NL3 .EQ. NL    ) NC4=NC0+6
         IF( NL3 .NE. NLMAX ) NC5=NBCALI
         NC3 = INDEX( KLG(NL3)(NC4:NC5) , 'FINPOUR' )
         IF( NC3 .GT. 0 ) THEN
C           EXISTE T IL UN AUTRE FINPOUR SUR CETTE LIGNE?
            NC3 = NC4 - 1 + NC3
 55         NQ  = INDEX( KLG(NL3)(NC3+1:NC5) , 'FINPOUR' )
            IF( NQ .GT. 0 ) THEN
               NC3 = NC3 + NQ
               GOTO 55
            ENDIF
         ENDIF
         NCC = INDEX( KLG(NL3)(NC4:NC5) , 'ENDFOR' )
         IF( NCC .GT. 0 ) THEN
C           EXISTE T IL UN AUTRE ENDFOR SUR CETTE LIGNE?
            NCC = NC4 - 1 + NCC
 57         NQ  = INDEX( KLG(NL3)(NCC+1:NC5) , 'ENDFOR' )
            IF( NQ .GT. 0 ) THEN
               NCC = NCC + NQ
               GOTO 57
            ENDIF
         ENDIF
         IF( NC3 .LE. 0 .AND. NCC .LE. 0 ) GOTO 60
C
C        RECHERCHE DU DERNIER DES 2
         IF( NC3 .LT. NCC ) THEN
C           ENDFOR AU DELA DE FINPOUR
            NC3  = NCC
         ENDIF
         GOTO 63
 60   CONTINUE
C
C     TRANSLATION DES CARACTERES POUR AJOUTER
C     NOM_VAR=NOM_VAR+1 ; ALLER ETIQ ; FINSI
C     A LA PLACE DE FINPOUR QUI DEBUTE EN NL3,NC3
C     FINPOUR EST BLANCHI
 63   IF( KLG(NL3)(NC4:NC4) .EQ. 'E' ) THEN
C        CAS 'ENDFOR'
         NCC = 5
      ELSE
C        CAS 'FINPOUR'
         NCC = 6
      ENDIF
      KLG(NL3)(NC3:NC3+NCC) = ' '
      NBCAR = NBCAR + NBCAR + 26
      CALL DEPKLG( NL3 , NC3 , NLMAX , NCMAX , NBCAR )
      IF( NL3 .LE. 0 ) GOTO 9900
      KLG(NL3)(NC3:NC3+NBCAR-1) = KLG(NL1)(NC1:NC2) // '='
     %// KLG(NL1)(NC1:NC2) // '+1 ; ALLER ' // ETIQUE(1:6) // ' ; FINSI'
C
C     EXISTE T IL D'AUTRES POUR OU FOR?
      GOTO 1
C
C     SORTIE
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
