      SUBROUTINE VAOPUN( NOPERA , D , NCODEV , DBLVAL )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCUL DE DBLVAL = ( OPERATEUR DE NUMERO NOPERA) D
C ----- POUR UN OPERATEUR UNAIRE
C
C ENTREES :
C ---------
C NOPERA : NUMERO DE L'OPERATION
C      SI  101 ALORS L'OPERATION EST NON
C      SI  102                       +
C      SI  103                       -
C      SI  104                       ABS
C      SI  105                       SQR         CF FORTRAN 77
C      SI  106                       LOG
C      SI  107                       LOG10
C      SI  108                       EXP
C      SI  109                       ASIN
C      SI  110                       ACOS
C      SI  111                       ATAN
C      SI  112                       SIN
C      SI  113                       COS
C      SI  114                       TAN
C      SI  115                       SINH
C      SI  116                       COSH
C      SI  117                       TANH
C      SI  118                       NINT
C      SI  119                       LLIGNE
C D     : L' OPERANDE REEL DOUBLE PRECISION
C
C SORTIES :
C ---------
C NCODEV : 0 DBLVAL N'EST PAS INITIALISEE
C          1 DBLVAL EST INITIALISEE
C DBLVAL : VALEUR REELLE DOUBLE PRECISION
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS  OCTOBRE 1989
C23456---------------------------------------------------------------012
      include"./incl/lu.inc"
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      DOUBLE PRECISION  DBLVAL,D
      LOGICAL           L
C
      GOTO ( 1 , 2 , 3 , 4 , 5 , 6 , 7 , 8 , 9 , 10 ,
     %      11 , 12 , 13 , 14 , 15 , 16 , 17 , 18 , 19 ),NOPERA-100
C
C     NON
 1    L = .NOT. ( D .NE. 0D0 )
      GOTO 1000
C
C     + UNAIRE
 2    DBLVAL = D
      GOTO 2000
C
C     - UNAIRE
 3    DBLVAL = -D
      GOTO 2000
C
C     ABS
 4    DBLVAL = ABS( D )
      GOTO 2000
C
C     SQRT
 5    IF( D .LT. 0D0 ) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(2)(1:25),'(G25.17)') D
         KERR(1) =  'LU: RACINE CARREE DE (' // KERR(2)(1:25)
     %           // ' ) IMPOSSIBLE'
         CALL LEREUR
         GOTO 9000
      ENDIF
      DBLVAL = SQRT( D )
      GOTO 2000
C
C     LOG
 6    IF( D .LE. 0D0 ) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(2)(1:25),'(G25.17)') D
         KERR(1) =  'LU: LOG DE ('// KERR(2)(1:25) //') IMPOSSIBLE'
         CALL LEREUR
         GOTO 9000
      ENDIF
      DBLVAL = LOG( D )
      GOTO 2000
C
C     LOG10
 7    IF( D .LE. 0D0 ) THEN
         WRITE(KERR(2)(1:25),'(G25.17)') D
         NBLGRC(NRERR) = 1
         KERR(1) =  'LU: LOG10 DE (' // KERR(2)(1:25) // ' ) IMPOSSIBLE'
         CALL LEREUR
         GOTO 9000
      ENDIF
      DBLVAL = LOG10( D )
      GOTO 2000
C
C     EXP
 8    DBLVAL = EXP( D )
      GOTO 2000
C
C     ASIN
 9    IF( ABS(D) .GT. 1D0 ) THEN
         WRITE(KERR(2)(1:25),'(G25.17)') D
         NBLGRC(NRERR) = 2
         KERR(1) =  'LU: ERREUR CALCUL DE ASIN('//KERR(2)(1:25)
         KERR(2) =  'LU: ABS(OPERANDE)>1'
         CALL LEREUR
         GOTO 9000
      ENDIF
      DBLVAL = ASIN( D )
      GOTO 2000
C
C     ACOS
 10   IF( ABS(D) .GT. 1D0 ) THEN
         WRITE(KERR(2)(1:25),'(G25.17)') D
         NBLGRC(NRERR) = 2
         KERR(1) =  'LU: ERREUR CALCUL DE ACOS('//KERR(2)(1:25)
         KERR(2) =  'LU: ABS(OPERANDE)>1'
         CALL LEREUR
         GOTO 9000
      ENDIF
      DBLVAL = ACOS( D )
      GOTO 2000
C
C     ATAN
 11   IF( ABS(D) .GT. RINFO( 'GRAND' ) ) THEN
C         VALEUR FORCEE A PI/2
          DBLVAL = ATAN( 1D0 ) * 2D0
          IF( D .LT. 0 ) DBLVAL = -DBLVAL
      ELSE
         DBLVAL = ATAN( D )
      ENDIF
      GOTO 2000
C
C     SIN
 12   DBLVAL = SIN( D )
      GOTO 2000
C
C     COS
 13   DBLVAL = COS( D )
      GOTO 2000
C
C     TAN
 14   DBLVAL = TAN( D )
      GOTO 2000
C
C     SINH
 15   DBLVAL = SINH( D )
      GOTO 2000
C
C     COSH
 16   DBLVAL = COSH( D )
      GOTO 2000
C
C     TANH
 17   DBLVAL = TANH( D )
      GOTO 2000
C
C     NINT
 18   DBLVAL = NINT( D )
      GOTO 2000
C
C     LLIGNE
 19   N = NINT( D )
C     N POSITION DANS KPACH DU NOM DE LA LIGNE
      CALL LLIGNE( KPACH(N) , DBLVAL )
      GOTO 2000
C
 1000 IF( L ) THEN
         DBLVAL = 1D0
      ELSE
         DBLVAL = 0D0
      ENDIF
C
 2000 NCODEV = 1
      RETURN
C
C     ERREUR
 9000 NCODEV = 0
      END
