      SUBROUTINE AFTAR2( NBIND, MIMX, MOT, T )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : AFFICHER LE TABLEAU REEL DOUBLE PRECISION DE NBIND INDICES
C -----
C
C ENTREES:
C --------
C NBIND  : NOMBRE D'INDICES
C MIMX   : MINIMUM ET MAXIMUM DE CHACUN DES NBIND INDICES DU TABLEAU
C MOT    : NOMBRE DE MOTS D'UNE VARIABLE DU TABLEAU
C T      : TABLEAU A AFFICHER
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 1994
C23456---------------------------------------------------------------012
      PARAMETER  (MXIND=5)
      include"./incl/td.inc"
      COMMON / UNITES / LECTEU, IMPRIM,  NUNITE(30)
      INTEGER           MIMX( 1:2, 1:NBIND )
      DOUBLE PRECISION  T( 0:* )
C
      CHARACTER*(MXIND) KNBC
      CHARACTER*66      KFORMA
C
C     UNE VARIABLE OU UN TABLEAU?
 5    IF( NBIND .LE. 0 ) THEN
C
C        UNE VARIABLE SEULEMENT
         NBV = 1
         CALL AFREE2( T )
         GOTO 9000
C
      ENDIF
C
C     UN TABLEAU AVEC AU MOINS UN INDICE
C     ==================================
C     LE NOMBRE DE CHIFFRES DE LA VARIATION DE CHACUN DES INDICES
C     LE NOMBRE TOTAL DE VARIABLES
      NBV = 1
      DO 10 I3 = 1, NBIND
         NBV = NBV * ( MIMX(2,I3) - MIMX(1,I3) + 1 )
         I1  = NBCHIF( MIMX(1,I3) )
         I2  = NBCHIF( MIMX(2,I3) )
         I1  = MAX( I1, I2 )
         WRITE( KNBC(I3:I3), '(I1)' ) I1
 10   CONTINUE
C
      IF( NBV .EQ. 1 ) THEN
C         TABLEAU D'UNE VARIABLE TRAITE COMME UNE VARIABLE
          NBIND = 0
          GOTO 5
      ENDIF
C
C     LA LIGNE ACTUELLE EST AFFICHEE
      CALL AFLIGN
C
C     LE NOMBRE DE VARIABLES POUR LE PREMIER INDICE
      NBV1   = MIMX(2,1) - MIMX(1,1) + 1
      KFORMA = ' '
C
C     AFFICHAGE SELON LE NOMBRE D'INDICES
C     ===================================
      IF( NBIND .EQ. 1 ) THEN
C
C        TABLEAU A 1 INDICE
C        ------------------
         KFORMA = '(5('' ('',I' // KNBC(1:1) // ','':'',G24.16,'')'')) '
         WRITE(IMPRIM,KFORMA)(I1,T(I1-MIMX(1,1)),I1=MIMX(1,1),MIMX(2,1))
C
      ELSE IF( NBIND .EQ. 2 ) THEN
C
C        TABLEAU A 2 INDICES
C        -------------------
         KFORMA='(5('' ('',I' // KNBC(1:1) // ',''_'',I' // KNBC(2:2) //
     %          ','':'',G24.16,'')'')) '
         I = 0
         DO 20 I2 = MIMX(1,2), MIMX(2,2)
            WRITE(IMPRIM,KFORMA) ( I1, I2, T(I+I1-MIMX(1,1)),
     %                             I1 = MIMX(1,1), MIMX(2,1) )
            I = I + NBV1
 20      CONTINUE
C
      ELSE IF( NBIND .EQ. 3 ) THEN
C
C        TABLEAU A 3 INDICES
C        -------------------
         KFORMA='(5('' ('',I' // KNBC(1:1) // ',''_'',I' // KNBC(2:2) //
     %                                        ',''_'',I' // KNBC(3:3) //
     %          ','':'',G24.16,'')'')) '
         I = 0
         DO 35 I3 = MIMX(1,3), MIMX(2,3)
            DO 30 I2 = MIMX(1,2), MIMX(2,2)
               WRITE(IMPRIM,KFORMA) ( I1, I2, I3, T(I+I1-MIMX(1,1)),
     %                                I1 = MIMX(1,1), MIMX(2,1) )
               I = I + NBV1
 30         CONTINUE
 35      CONTINUE
C
C
      ELSE IF( NBIND .EQ. 4 ) THEN
C
C        TABLEAU A 4 INDICES
C        -------------------
         KFORMA='(5('' ('',I' // KNBC(1:1) // ',''_'',I' // KNBC(2:2) //
     %                                        ',''_'',I' // KNBC(3:3) //
     %                                        ',''_'',I' // KNBC(4:4) //
     %            ','':'',G24.16,'')'')) '
         I = 0
         DO 48 I4 = MIMX(1,4), MIMX(2,4)
            DO 44 I3 = MIMX(1,3), MIMX(2,3)
               DO 40 I2 = MIMX(1,2), MIMX(2,2)
                  WRITE(IMPRIM,KFORMA) ( I1,I2,I3,I4,T(I+I1-MIMX(1,1)),
     %                                   I1 = MIMX(1,1), MIMX(2,1) )
                  I = I + NBV1
 40            CONTINUE
 44         CONTINUE
 48      CONTINUE
C
C
      ELSE IF( NBIND .EQ. 5 ) THEN
C
C        TABLEAU A 5 INDICES
C        -------------------
         KFORMA='(5('' ('',I' // KNBC(1:1) // ',''_'',I' // KNBC(2:2) //
     %                                        ',''_'',I' // KNBC(3:3) //
     %                                        ',''_'',I' // KNBC(4:4) //
     %                                        ',''_'',I' // KNBC(5:5) //
     %          ','':'',G24.16,'')'')) '
         I = 0
         DO 59 I5 = MIMX(1,5), MIMX(2,5)
         DO 58 I4 = MIMX(1,4), MIMX(2,4)
            DO 54 I3 = MIMX(1,3), MIMX(2,3)
               DO 50 I2 = MIMX(1,2), MIMX(2,2)
                  WRITE(IMPRIM,KFORMA)(I1,I2,I3,I4,I5,T(I+I1-MIMX(1,1)),
     %                                 I1 = MIMX(1,1), MIMX(2,1) )
                  I = I + NBV1
 50            CONTINUE
 54         CONTINUE
 58      CONTINUE
 59      CONTINUE
C
      ELSE
C
C        PLUS DE MXIND INDICES
         WRITE(IMPRIM,*) 'ERREUR AFTAR2: PLUS DE ',MXIND,' INDICES'
         CALL XVPAUSE
      ENDIF
C
C     DECALAGE DANS LE TABLEAU AU DELA POUR ATTEINDRE LA DERNIERE VARIABLE
 9000 LDTS( LHTMS ) = LDTS( LHTMS ) + MOT * NBV
      END
