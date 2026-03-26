      SUBROUTINE INTERP( NOINTE, X, Y, Z,   NBP, P )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : RETOURNER SELON LE NO DE L INTERPOLATION ET
C ----- LE POINT (X,Y,Z) DE L ELEMENT DE REFERENCE
C       LE NOMBRE DE FONCTIONS DE BASE
C       LES VALEURS DE CES FONCTIONS AU POINT (X,Y,Z)
C
C PARAMETRES D ENTREE :
C ---------------------
C NOINTE : NO DE L INTERPOLATION
C          NOINTE = NDIM * 1000 + (J - 1) * 30 + I
C                   NDIM DIMENSION DE L ESPACE ( 1 OU 2 OU 3 )
C                   J NO DU TYPE D ELEMENT (NDIM=1 J=1 SEGMENT
C                                           NDIM=2 J=1 TRIANGLE
C                                                  J=2 QUADRANGLE
C                                           NDIM=3 J=1 TETRAEDRE
C                                                  J=2 PENTAEDRE
C                                                  J=3 HEXAEDRE
C                                                  J=4 PYRAMIDE
C                                           NDIM=6 J=1 6-CUBE )
C                   I NO DANS LE TYPE DE L ELEMENT FINI
C X,Y,Z  : COORDONNEES DU POINT DE L ELEMENT DE REFERENCE OU LES
C          VALEURS DES FONCTIONS DE BASE DOIVENT ETRE CALCULEES
C
C PARAMETRES RESULTATS :
C ----------------------
C NBP    : NOMBRE DE FONCTIONS DE BASE
C P      : NBP VALEURS DES FONCTIONS DE BASE AU POINT (X,Y,Z)
C          VALEURS REELLES DOUBLE PRECISION
C
C ATTENTION : LES FONCTIONS ET POINTS SONT DEFINIS SUR L ELEMENT DE
C ----------- REFERENCE RECTANGLE UNITE
C             POUR LES ELEMENTS DE TYPE SERENDIP LE NOMBRE DE MONOMES
C             DES POLYNOMES PEUT EXCEDER LE NOMBRE DE POLYNOMES
C             EXEMPLE : HEXA 3Q2C  NBP=20 NBRE MONOMES=27
C             LA PROGRAMMATION D UNE NOUVELLE INTERPOLATION DOIT VEILLER
C             A NE PAS DEPASSER LE NOMBRE DE POLYNOMES
C
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS        AOUT 1981
C ......................................................................
      include"./incl/gsmenu.inc"
      DOUBLE PRECISION  D1,D2,D3,D4,S1,S2,S3,S4,S5,S6,X,Y,Z,P(*)
      DOUBLE PRECISION  X1,Y1,Z1,X2,Y2,Z2
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
C
C     DIMENSION DE L ESPACE DE L INTERPOLATION  (1 OU 2 OU 3 )
C     --------------------------------------------------------
      NDIM  = NOINTE / 1000
      NUM   = NOINTE - NDIM * 1000
      IF( NDIM .LE. 0 .OR. NDIM .GT. 3 ) GOTO 1
C
C     ******************************************************************
      GOTO ( 1000 , 2000 , 3000 ) , NDIM
C     ******************************************************************
 1    NBLGRC(NRERR) = 1
      WRITE(KERR(MXLGER)(1:4),'(I4)') NOINTE
      KERR(1) = 'INTERP:INTERPOLATION '//KERR(MXLGER)(1:4)//' INCONNUE'
      CALL LEREUR
      CALL ARRET(100)
      RETURN
C
C     &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     INTERPOLATION DANS R ** 1
C     &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
 1000 GOTO ( 1001, 1002, 1 ) , NUM
C     ==================================================================
C     SEGM 1P1D
C     ==================================================================
 1001 NBP   = 2
      P( 1) = 1D0 - X
      P( 2) = X
      RETURN
C
C     ==================================================================
C     SEGM 1P2D
C     ==================================================================
 1002 NBP   = 3
      P( 1) = ( 1D0 - X ) * ( 1D0 - 2D0 * X )
      P( 2) =      X      * ( 2D0 * X - 1D0 )
      P( 3) = 4D0 * X * ( 1D0 - X )
      RETURN
C
C     &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     INTERPOLATION DANS R ** 2
C     &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
 2000 IF( NUM .GT. 30 ) GOTO 2030
C
C     TRIANGLES
C     .........
C
      GOTO ( 2001 , 2002 , 2003 , 2004 , 2005 , 2006 ) , NUM
C
C     ==================================================================
C     TRIA AP1D ET TRIA 2P1D
C     ==================================================================
C
 2001 NBP   = 3
      P( 1) = 1.D0 - X - Y
      P( 2) = X
      P( 3) = Y
      RETURN
C
C     ==================================================================
C     TRIA AP2C ET TRIA 2P2C
C     ==================================================================
C
 2002 NBP   = 6
      D1    = 1.D0 - X - Y
      P( 1) = D1   * ( 2.D0 * D1 - 1.D0 )
      P( 2) = 2.D0 * X * X - X
      P( 3) = 2.D0 * Y * Y - Y
      D1    = 4.D0 * D1
      P( 4) = D1   * X
      P( 5) = 4.D0 * X * Y
      P( 6) = D1   * Y
      RETURN
C
C     ==================================================================
C     TRIANGLE BREZZI-FORTIN  P1 + BULLE CB1 CB2 CB3
C     ==================================================================
C     Version avec le degre de liberte v(BARYCENTRE)
 2003 NBP = 4
      D1 = 1D0 - X - Y
      D3 = 9D0 * D1 * X * Y
      P( 1) =  D1 - D3
      P( 2) =  X  - D3
      P( 3) =  Y  - D3
      P( 4) = 3D0 * D3
      RETURN
CCC      DOUBLE PRECISION  P(4,4,4)
ccc      Version avec le degre de liberte Int v dX / mes e
CCCC     LES 4 POLYNOMES DE BASE DE BREZZI-FORTIN SUR LE TRIANGLE REFERENCE
CCC      DATA P/ 1D0,-1D0,0D0,0D0, -1D0,-20D0,20D0,0D0, 0D0,20D0,0D0,0D0,
CCC     %        0D0,0D0,0D0,0D0,
CCC     %        0D0,1D0,0D0,0D0,   0D0,-20D0,20D0,0D0, 0D0,20D0,0D0,0D0,
CCC     %        0D0,0D0,0D0,0D0,
CCC     %        0D0,0D0,0D0,0D0,   1D0,-20D0,20D0,0D0, 0D0,20D0,0D0,0D0,
CCC     %        0D0,0D0,0D0,0D0,
CCC     %        0D0,0D0,0D0,0D0,   0D0,60D0,-60D0,0D0, 0D0,-60D0,0D0,0D0,
CCC     %        0D0,0D0,0D0,0D0 /
C
C     ==================================================================
C     TRIA MT10
C     ==================================================================
C
 2004 GOTO 1
C
C     ==================================================================
C     TRIA MT21
C     ==================================================================
C
 2005 GOTO 1
C
C     ==================================================================
C     TRIA EQ06 : INTERPOLATION SUR CHAQUE ARETE INDEPENDAMMENT
C     ==================================================================
C
 2006 NBP = 6
      D1  = DSQRT(3.D0)
      D2  = (1.D0+D1)/2.D0
      D3  = (1.D0-D1)/2.D0
      P( 1) = -D1*X+D2
      P( 2) =  D1*X+D3
      P( 3) = (-D1*(X-Y)+1.D0)/2.D0
      P( 4) = (D1*(X-Y)+1.D0)/2.D0
      P( 5) =  D1*Y+D3
      P( 6) = -D1*Y+D2
      RETURN
C
C     QUADRANGLES
C     ...........
C
 2030 NUM  = NUM - 30
      GOTO ( 2031 , 2032 , 1 , 2034 , 2035) , NUM
C
C     ==================================================================
C     QUAD AQ1C ET QUAD 2Q1C
C     ==================================================================
C
 2031 NBP   = 4
      P( 1) = 1.D0 - X - Y + X * Y
      P( 2) = X - X * Y
      P( 3) = X * Y
      P( 4) = Y - X * Y
      RETURN
C
C     ==================================================================
C     QUAD AQ2C ET QUAD 2Q2C
C     ==================================================================
C
C     QUADRANGLE Serendipity Valeur imposee au barycentre
C     Cf Cours Thomas Raviart 1972  page X-22
C     P(Barycentre) = - SOM 1 a 4 P(Sommet) / 4 + Som 5 a 8 P(Milieu) / 2
 2032 NBP   = 8
      D1    = 1.D0 - X
      D2    = 1.D0 - Y
      P( 1) = ( D1 - Y + X * Y ) * ( 1.D0 - 2.D0 * ( X + Y ) )
      P( 2) = X * D2 * ( -1.D0 + 2.D0 * ( X - Y ) )
      P( 3) = X * Y * ( 2.D0 * ( X + Y ) - 3.D0 )
      P( 4) = Y * D1 * (-1.D0 + 2.D0 * ( Y - X ) )
C
      D1    = 4.D0 * D1
      P( 5) = D1 * X * D2
      P( 6) = 4.D0 * X * Y * D2
      P( 7) = X * Y * D1
      P( 8) = Y * D1 * D2
      RETURN
C
C     ==================================================================
C     QUAD MQ10
C     ==================================================================
C
 2034 GOTO 1
C
C     ==================================================================
C     QUAD MQ21
C     ==================================================================
C
 2035 GOTO 1
C
C     &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     INTERPOLATION DANS R ** 3
C     &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
 3000 IF( NUM .GT. 30 ) GOTO 3030
C
C     LES TETRAEDRES
C     ..............
C
C     ******************************************************************
      GOTO ( 3001 , 3002 , 3003 , 1 , 1 , 1 ) , NUM
C     ******************************************************************
C
C     ==================================================================
C     TETR 3P1D
C     ==================================================================
C
 3001 NBP  = 4
      P(1) = 1.D0 - X - Y - Z
      P(2) = X
      P(3) = Y
      P(4) = Z
      RETURN
C
C     ==================================================================
C     TETR 3P2C
C     ==================================================================
C
 3002 NBP   = 10
      D1    = 1.D0 - X - Y - Z
      P( 1) =  D1 * (2D0 * D1 - 1D0)
      P( 2) =   X * (2D0 *  X - 1D0)
      P( 3) =   Y * (2D0 *  Y - 1D0)
      P( 4) =   Z * (2D0 *  Z - 1D0)
      P( 5) = 4D0 * D1 * X
      P( 6) = 4D0 *  X * Y
      P( 7) = 4D0 *  Y * D1
      P( 8) = 4D0 * D1 * Z
      P( 9) = 4D0 *  X * Z
      P(10) = 4D0 *  Y * Z
      RETURN
C
C     ==================================================================
C     TETRAEDRE BREZZI-FORTIN  P1 + BULLE CB1 CB2 CB3 CB4
C     ==================================================================
C     Version avec le degre de liberte v(BARYCENTRE)
 3003 NBP = 5
      D1 =  1D0 - X - Y - Z
      D4 = 64D0 * D1 * X * Y * Z
      P( 1) = D1  - D4
      P( 2) =  X  - D4
      P( 3) =  Y  - D4
      P( 4) =  Z  - D4
      P( 5) = 4D0 * D4
C
CCCC  LES 5 POLYNOMES DE BASE DE BREZZI-FORTIN SUR LE TETRAEDRE REFERENCE
ccc   Version avec le degre de liberte Int v dX / mes e
CCC   DOUBLE PRECISION P(5,5,5, 5)
CCC      DATA P/ 1D0,-1D0,3*0D0,  -1D0, 19*0D0,
CCC     %       -1D0, 4*0D0, 0D0,-210D0,210D0,2*0D0, 0D0,210D0,13*0D0,
CCC     %        6*0D0,210D0,68*0D0,
CCC
CCC     %        0D0,1D0,23*0D0,
CCC     %        6*0D0,-210D0,210D0,3*0D0,210D0,13*0D0,
CCC     %        6*0D0,210D0,68*0D0,
CCC
CCC     %        5*0D0,1D0,19*0D0,
CCC     %        6*0D0,-210D0,210D0,3*0D0,210D0,13*0D0,
CCC     %        6*0D0,210D0,68*0D0,
CCC
CCC     %        25*0D0,
CCC     %        1D0,5*0D0,-210D0,210D0,3*0D0,210D0,13*0D0,
CCC     %        6*0D0,210D0,68*0D0,
CCC
CCC     %        31*0D0, 840D0,-840D0,3*0D0,-840D0,13*0D0,
CCC     %        6*0D0,-840D0,68*0D0 /
C
      RETURN
C
 3030 NUM    = NUM - 30
      IF(NUM .GT. 30 ) GOTO 3060
C
C
C     LES PENTAEDRES
C     ..............
C
C     ******************************************************************
      GOTO ( 3031 , 3032 , 1 , 1 , 1 , 1 , 1 ) , NUM
C     ******************************************************************
C
C     ==================================================================
C     PENT 3R1C
C     ==================================================================
C
 3031 NBP   = 6
      D1    = 1D0 - X - Y
      D2    = 1D0 - Z
      P( 1) = D1 * D2
      P( 2) =  X * D2
      P( 3) =  Y * D2
      P( 4) = D1 * Z
      P( 5) =  X * Z
      P( 6) =  Y * Z
      RETURN
C
C     ==================================================================
C     PENT 3R2C
C     ==================================================================
C
 3032 NBP   = 15
C
C     LA VALEUR DES POLYNOMES DE BASE DU TRIANGLE P2 UNITE
C
      D1    = 1D0 - X - Y
      D2    = X * (2D0 * X - 1D0)
      D3    = Y * (2D0 * Y - 1D0)
      S1    = 4D0 * D1 * X
      S2    = 4D0 *  X * Y
      S3    = 4D0 *  Y * D1
      D1    = D1  * (2D0 * D1 - 1D0)
C
C     LA VALEUR DES POLYNOMES R2 COMPLET
C
      D4    = (1D0 - Z) * (1D0 - 2D0 * Z)
      P( 1) = D1 * D4
      P( 2) = D2 * D4
      P( 3) = D3 * D4
      P( 7) = S1 * D4
      P( 8) = S2 * D4
      P( 9) = S3 * D4
      D4    = 4D0 * Z * (1D0 - Z)
      P(10) = D1 * D4
      P(11) = D2 * D4
      P(12) = D3 * D4
      Y1    = S1 * D4
      Z1    = S2 * D4
      X1    = S3 * D4
      D4    =  Z * (2D0 * Z - 1D0)
      P( 4) = D1 * D4
      P( 5) = D2 * D4
      P( 6) = D3 * D4
      P(13) = S1 * D4
      P(14) = S2 * D4
      P(15) = S3 * D4
C
C     ELIMINATION DU BARYCENTRE DES 3 FACES CARREES
C
      D1    = ( X1 + Y1 ) * 0.25D0
      D2    = ( Y1 + Z1 ) * 0.25D0
      D3    = ( Z1 + X1 ) * 0.25D0
      P( 1) = P( 1) - D1
      P( 2) = P( 2) - D2
      P( 3) = P( 3) - D3
      P( 4) = P( 4) - D1
      P( 5) = P( 5) - D2
      P( 6) = P( 6) - D3
      P( 7) = P( 7) + Y1 * 0.5D0
      P( 8) = P( 8) + Z1 * 0.5D0
      P( 9) = P( 9) + X1 * 0.5D0
      P(10) = P(10) + D1 * 2.0D0
      P(11) = P(11) + D2 * 2.0D0
      P(12) = P(12) + D3 * 2.0D0
      P(13) = P(13) + Y1 * 0.5D0
      P(14) = P(14) + Z1 * 0.5D0
      P(15) = P(15) + X1 * 0.5D0
      RETURN
C
C     LES HEXAEDRES
C     .............
C
 3060 NUM    = NUM - 30
      IF(NUM .GT. 30 ) GOTO 3090
C
C     ******************************************************************
      GOTO ( 3061 , 3062 , 1 , 1 , 1 , 1 ) , NUM
C     ******************************************************************
C
C     ==================================================================
C     HEXA 3Q1C
C     ==================================================================
C
 3061 NBP    = 8
      D1     = 1.D0 - X
      D2     = 1.D0 - Y
      D3     = 1.D0 - Z
C
      P(1) = D1 * D2 * D3
      P(2) =  X * D2 * D3
      P(3) =  X *  Y * D3
      P(4) = D1 *  Y * D3
      P(5) = D1 * D2 *  Z
      P(6) =  X * D2 *  Z
      P(7) =  X *  Y *  Z
      P(8) = D1 *  Y *  Z
      RETURN
C
C     ==================================================================
C     HEXA 3Q2C
C     ==================================================================
C
 3062 X1    = 1D0 - X
      Y1    = 1D0 - Y
      Z1    = 1D0 - Z
      X2    = 1D0 - 2D0 * X
      Y2    = 1D0 - 2D0 * Y
      Z2    = 1D0 - 2D0 * Z
C
      S4    =  X2 * Y2 * Z2
      D4    =  Z1 * S4
      P( 1) =  X1 * Y1 * D4
      P( 2) = - X * Y1 * D4
      P( 3) =   X *  Y * D4
      P( 4) = -X1 *  Y * D4
C
      D4    =   Z * S4
      P( 5) = -X1 * Y1 * D4
      P( 6) =   X * Y1 * D4
      P( 7) = - X *  Y * D4
      P( 8) =  X1 *  Y * D4
C
      S4    =  4D0 * X
      D1    =   S4 * X1 * Y1 * Y2
      D2    = - S4 * X2 * Y  * Y1
      D3    = - S4 * X1 * Y  * Y2
      D4    =   X1 * X2 * Y  * Y1 * 4D0
C
      S4    =  Z1 * Z2
      P( 9) =  D1 * S4
      P(10) =  D2 * S4
      P(11) =  D3 * S4
      P(12) =  D4 * S4
C
C
      S4    = -Z  * Z2
      P(17) =  D1 * S4
      P(18) =  D2 * S4
      P(19) =  D3 * S4
      P(20) =  D4 * S4
C
      D1    =   Z * Z1 * 4D0
      S4    =  Y1 * Y2 * D1
      P(13) =  X1 * X2 * S4
      P(14) = -X  * X2 * S4
      S4    =  Y  * Y2 * D1
      P(15) =  X  * X2 * S4
      P(16) = -X1 * X2 * S4
C
      D1 =   X * X1 * Y  * Y1 * 16D0
C     S1=P(21),S2=P(22),S3=P(23),S4=P(24),S5=P(25),S6=P(26)
      S1 =  D1 * Z1 * Z2
      S4 = -D1 * Z  * Z2
      D1 =  Y  * Y1 * Z  * Z1 * 16D0
      S2 =  X1 * X2 * D1
      S5 = -X  * X2 * D1
      D1 =  X  * X1 * Z  * Z1 * 16D0
      S3 =  Y1 * Y2 * D1
      S6 = -Y  * Y2 * D1
C
C     D1 = P(27)
      D1 =  X  * X1 * Y  * Y1 * Z  * Z1 * 16D0
C
C     ELIMINATION DU BARYCENTRE DES 6 FACES CARREES
      NBP = 20
      P( 1) = P( 1) - (S1 + S2 + S3) * 0.25D0 - D1
      P( 2) = P( 2) - (S1 + S3 + S5) * 0.25D0 - D1
      P( 3) = P( 3) - (S1 + S5 + S6) * 0.25D0 - D1
      P( 4) = P( 4) - (S1 + S2 + S6) * 0.25D0 - D1
      P( 5) = P( 5) - (S2 + S3 + S4) * 0.25D0 - D1
      P( 6) = P( 6) - (S3 + S4 + S5) * 0.25D0 - D1
      P( 7) = P( 7) - (S4 + S5 + S6) * 0.25D0 - D1
      P( 8) = P( 8) - (S2 + S4 + S6) * 0.25D0 - D1
      P( 9) = P( 9) + (S1 + S3) * 0.5D0 + D1
      P(10) = P(10) + (S1 + S5) * 0.5D0 + D1
      P(11) = P(11) + (S1 + S6) * 0.5D0 + D1
      P(12) = P(12) + (S1 + S2) * 0.5D0 + D1
      P(13) = P(13) + (S2 + S3) * 0.5D0 + D1
      P(14) = P(14) + (S3 + S5) * 0.5D0 + D1
      P(15) = P(15) + (S5 + S6) * 0.5D0 + D1
      P(16) = P(16) + (S2 + S6) * 0.5D0 + D1
      P(17) = P(17) + (S3 + S4) * 0.5D0 + D1
      P(18) = P(18) + (S4 + S5) * 0.5D0 + D1
      P(19) = P(19) + (S4 + S6) * 0.5D0 + D1
      P(20) = P(20) + (S2 + S4) * 0.5D0 + D1
      RETURN
C
C     LES PYRAMIDES
C     .............
C
 3090 NUM = NUM - 30
      IF(NUM .GT. 30 ) GOTO 3090
C
C     ******************************************************************
      GOTO ( 3091 , 3092 , 1 , 1 , 1 , 1 ) , NUM
C     ******************************************************************
C
C     ==================================================================
C     PYRA 3PY1
C     ==================================================================
 3091 NBP= 5
      X1 = 1.D0 - X
      Y1 = 1.D0 - Y
      Z1 = 1.D0 - Z
      P(1) = X1 * Y1 * Z1
      P(2) = X  * Y1 * Z1
      P(3) = X  * Y  * Z1
      P(4) = X1 * Y  * Z1
      P(5) = Z
      RETURN
C
C     ==================================================================
C     PYRA 3PY2
C     ==================================================================
 3092 NBP = 13
C
C     LA BASE DE LA PYRAMIDE EN Z=0
      X1 = 1.D0 - X
      Y1 = 1.D0 - Y
      Z1 = 1.D0 - Z
      S3 = ( 1.D0 - Z ) * ( 1.D0 - 2D0 * Z )
C     LES 4 SOMMETS DU BAS DE LA PYRAMIDE
      P( 1) = ( X1 - Y + X * Y ) * ( 1.D0 - 2.D0 * ( X + Y ) ) * S3
      P( 2) =   X * Y1 * ( -1.D0 + 2.D0 * ( X - Y ) ) * S3
      P( 3) =   X * Y  * (  2.D0 * ( X + Y ) - 3.D0 ) * S3
      P( 4) =   Y * X1 * ( -1.D0 + 2.D0 * ( Y - X ) ) * S3
C
C     LES 4 MILIEUX D'ARETE DE LA BASE DE LA PYRAMIDE
      D3    = 4.D0 * X1
      P( 6) = D3 * X * Y1 * S3
      P( 7) = 4.D0 * X * Y * Y1 * S3
      P( 8) = X *  Y * D3  * S3
      P( 9) = Y * D3 * Y1  * S3
C
C     LE PLAN MEDIAN EN Z=1/2
      S3 = 4D0 * Z * Z1
      X2 = 1D0 - 2D0 * X
      Y2 = 1D0 - 2D0 * Y
C     LES 4 MILIEUX DU PLAN MEDIAN DE LA PYRAMIDE
      P(10) =  X2 * Y2  * S3
      P(11) = 2D0 * X   * Y2 * S3
      P(12) = 4D0 * X   * Y  * S3
      P(13) =  X2 * 2D0 * Y  * S3
C
C     LE SOMMET 5 EN HAUT DE LA PYRAMIDE
C     CE NUMERO POUR ATTEINDRE DE 1 A 5 LES 5 SOMMETS DE LA PYRAMIDE
C     PARMI LES 13 POINTS=NOEUDS
      P( 5) = Z * ( 2D0 * Z - 1D0 )
C
      RETURN
      END
