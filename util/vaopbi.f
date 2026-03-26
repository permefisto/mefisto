      SUBROUTINE VAOPBI( NOPERA , D1 , D2 , NCODEV , DBLVAL )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCUL DE DBLVAL = D1 ( OPERATEUR DE NUMERO NOPERA) D2
C ----- POUR UN OPERATEUR BINAIRE
C
C ENTREES :
C ---------
C NOPERA : NUMERO DE L'OPERATION
C      SI  201 => OX
C      SI  202    OU
C      SI  203    ET
C      SI  204    <
C      SI  205    <=            CF FORTRAN 77
C      SI  206    =
C      SI  207    <>
C      SI  208    >=
C      SI  209    >
C      SI  210    +
C      SI  211    -
C      SI  212    *
C      SI  213    /
C      SI  214    **
C      SI  215    MIN
C      SI  216    MAX
C      SI  217    MOD
C      SI  218    ISIGN ISIGN(X1,X2) = + |X1| SI X2 >=0
C                                      - |X1| SI X2 < 0
C      SI  219    NUOBJT( 'TYPE_OBJET' , 'NOM_OBJET' )   AVEC
C                         'TYPE_OBJET' = 'POINT'    'LIGNE'
C                                        'SURFACE'  'VOLUME'
C                                        'OBJET'
C                                        'FONCTION' 'TRANSFO'
C                  RETOURNE LE NUMERO DANS LE LEXIQUE DE CET OBJET
C                           0 SI NOM NON RETROUVE
C
C      SI  220     DIST2P( 'NOM_POINT_1' , 'NOM_POINT_2' )
C                  RETOURNE LA DISTANCE ENTRE LES 2 POINTS
C                          -1 SI l'UN DES POINTS N'EST PAS RETROUVE
C
C D1 , D2 : LES OPERANDES REEL DOUBLE PRECISION  OU BIEN
C           LES REELS DOUBLE PRECISION NUMEROS DES CONSTANTES CHAINES
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
      DOUBLE PRECISION  DBLVAL,D1,D2
      INTEGER           N1,N2
      LOGICAL           L
C
      GOTO ( 1 ,  2 ,  3 ,  4 ,  5 ,  6 ,  7 ,  8 ,  9 , 10 ,
     %      11 , 12 , 13 , 14 , 15 , 16 , 17 , 18 , 19 , 20 ),
     %      NOPERA-200
C
C     OX
 1    L = D1 .EQ. D2
      GOTO 1000
C
C     OU
 2    L = D1.NE.0D0 .OR. D2.NE.0D0
      GOTO 1000
C
C     ET
 3    L = D1.NE.0D0 .AND. D2.NE.0D0
      GOTO 1000
C
C     <
 4    L = D1 .LT. D2
      GOTO 1000
C
C     <
 5    L = D1 .LE. D2
      GOTO 1000
C
C     =
 6    L = D1 .EQ. D2
      GOTO 1000
C
C     <>
 7    L = D1 .NE. D2
      GOTO 1000
C
C     >=
 8    L = D1 .GE. D2
      GOTO 1000
C
C     >
 9    L = D1 .GT. D2
      GOTO 1000
C
C     +
 10   DBLVAL = D1 + D2
      GOTO 2000
C
C     -
 11   DBLVAL = D1 - D2
      GOTO 2000
C
C     *
 12   DBLVAL = D1 * D2
      GOTO 2000
C
C     /
 13   IF( D2 .EQ. 0D0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'LU: ERREUR /0 INTERDITE'
         CALL LEREUR
         GOTO 9000
      ENDIF
      DBLVAL = D1 / D2
      GOTO 2000
C
C     **
 14   N2 = NINT( D2 )
      IF( N2 .EQ. D2 ) THEN
         IF( D1 .GE. 0 ) THEN
C           TEST SUITE A DES PROBLEMES APOLLO
            DBLVAL = D1 ** N2
         ELSE
            IF( MOD(N2,2) .EQ. 0 ) THEN
               DBLVAL = (-D1) ** N2
            ELSE
               DBLVAL = -( (-D1) ** N2 )
            ENDIF
         ENDIF
      ELSE
         IF( D1 .LT. 0D0 ) THEN
            NBLGRC(NRERR) = 2
            WRITE(KERR(3)(1:25),'(G25.17)') D1
            WRITE(KERR(4)(1:25),'(G25.17)') D2
            KERR(1) =  'LU: ERREUR CALCUL DE '//KERR(3)(1:25)//
     %                 ' ** ' // KERR(4)(1:25)
            KERR(2) =  'LU: LE PREMIER OPERANDE DOIT ETRE >=0'
            CALL LEREUR
            GOTO 9000
         ENDIF
         DBLVAL = D1 ** D2
      ENDIF
      GOTO 2000
C
C     MIN
 15   DBLVAL = MIN(D1,D2)
      GOTO 2000
C
C     MAX
 16   DBLVAL = MAX(D1,D2)
      GOTO 2000
C
C     MOD
 17   N1 = NINT( D1 )
      N2 = NINT( D2 )
      IF( N1 .LT. 0 .OR. N2 .LE. 0 ) THEN
         NBLGRC(NRERR) = 2
         WRITE(KERR(3)(1:12),'(I12)') N1
         WRITE(KERR(4)(1:12),'(I12)') N2
         KERR(1) =  'LU: ERREUR CALCUL DE '//KERR(3)(1:12)//
     %              ' MODULO ' // KERR(4)(1:12)
         KERR(2) =  'LU: LE PREMIER OPERANDE DOIT ETRE >=0'
         CALL LEREUR
         GOTO 9000
      ENDIF
      DBLVAL = MOD(N1,N2)
      GOTO 2000
C
C     ISIGN
 18   N1 = NINT( D1 )
      N2 = NINT( D2 )
      DBLVAL = ISIGN( N1 , N2 )
      GOTO 2000
C
C     NUOBJT
 19   N1 = NINT( D1 )
      N2 = NINT( D2 )
C     N1 N2 SONT LES POSITIONS DANS KPACH
C     DU TYPE DE L'OBJET ET DU NOM DE L'OBJET
      CALL NUOBNM( KPACH(N1) , KPACH(N2) , NCODEV )
      DBLVAL = NCODEV
      GOTO 2000
C
C     DIST2P
 20   N1 = NINT( D1 )
      N2 = NINT( D2 )
C     N1 N2 SONT LES POSITIONS DANS KPACH
C     DES 2 NOMS DES POINTS
      CALL DI2PUT( KPACH(N1) , KPACH(N2) , DBLVAL )
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
