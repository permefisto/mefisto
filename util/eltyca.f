      SUBROUTINE ELTYCA ( NO )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : FOURNIR SELON LE NUMERO DE TYPE DE L'EF
C ----- SES CARACTERISTIQUES GEOMETRIQUES ET D'INTERPOLATION NECESSAIRES
C       A LA CONSTRUCTION DES TABLEAUX POUR CET EF
C
C PARAMETRES D ENTREE :
C ---------------------
C NO     : NO DE L EF DANS LES SP UTILITAIRES
C          L INTERPOLATION A ETE PRISE EN COMPTE
C
C PARAMETRES RESULTATS CONTENUS DANS LE COMMON / PONOEL / :
C ---------------------------------------------------------
C NBPOE  : NOMBRE DE POINTS DE L EF NO
C NBNOE  : NOMBRE DE NOEUDS DE L EF NO
C NOTRAE : CODE TRAITEMENT  DE L EF NO
C NBNSOM : NOMBRE DE NOEUDS-SOMMETS DE L EF
C NBPOIN : NOMBRE DE POINTS INTERNES DE L EF
C NBNOIN : NOMBRE DE NOEUDS INTERNES DE L EF
C NARET  : NOMBRE DE SES ARETES
C NOSOAR : NO DES 2 SOMMETS DE CHACUNE DE SES ARETES
C NOTYAR : NO DU TYPE       DE CHACUNE DE SES ARETES
C          ( CF SP TYARCP )
C NBPOAR : NOMBRE DE POINTS-NON SOMMETS DE CHACUNE DE SES ARETES
C NOPOAR : NO ELEMENTAIRE DES POINTS NON SOMMETS DE CHACUNE DE SES
C          ARETES
C NBNOAR : NOMBRE DE NOEUDS-NON SOMMETS DE CHACUNE DE SES ARETES
C NONOAR : NO ELEMENTAIRE DES NOEUDS NON SOMMETS DE CHACUNE DE SES
C          ARETES
C NFACE  : NOMBRE DE SES FACES
C NBSOFA : NOMBRE DE SOMMETS DE CHACUNE DE SES FACES
C NOSOFA : NO ELEMENTAIRE DES SOMMETS DE CHACUNE DE SES FACES
C NOTYFA : NO DU TYPE       DE CHACUNE DE SES FACES
C NBPOFA : NOMBRE DE POINTS-NON SUR LES ARETES DE CHACUNE DE SES FACES
C NOPOFA : NO ELEMENTAIRE DES POINTS-NON SUR LES ARETES DE CHACUNE DE
C          SES FACES
C NBNOFA : NOMBRE DE NOEUDS-NON SUR LES ARETES DE CHACUNE DE SES FACES
C NONOFA : NO ELEMENTAIRE DES NOEUDS-NON SUR LES ARETES DE CHACUNE DE
C          SES FACES
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE UPMC PARIS    OCTOBRE 1990
C2345X7..............................................................012
      include"./incl/ponoel.inc"
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
C
C     ATTENTION : SI DECLARATION INSUFFISANTE MODIFIER CES VALEURS ICI
C              ET DANS LE SP DEFELE OU ILS SONT EFFECTIVEMENT DECLARES
C     =================================================================
C     NCOGEL LE CODE GEOMETRIQUE EST CALCULE
      CALL ELNUCG( NO, NCOGEL )
C
C     LE NO DES SOMMETS DES ARETES ET DES FACES
C     -----------------------------------------
      CALL SOARFA( NCOGEL, NOSOAR, NOSOFA )
C
C     LE NOMBRE DE POINTS,NOEUDS ET CODE TRAITEMENT
C     ---------------------------------------------
      CALL ELPONO( NO, NBPOE, NBNOE, NOTRAE, I )
C
C     A PRIORI PAS DE POINT NI NOEUD INTERNE A L'EF
      NBPOIN = 0
      NBNOIN = 0
C
C     LE NOMBRE D'ARETES DE L'EF
C     --------------------------
      NARET = NBARET( NCOGEL )
C
C     LE NOMBRE DE FACES DE L'EF
C     --------------------------
      NFACE = NBFACE( NCOGEL )
C
C     LE NOMBRE ET NO DES ARETES DES FACES
C     ------------------------------------
      CALL NUARFA( NCOGEL, NBARFA, NOARFA )
C
C     NOEUD=POINT=SOMMET
C     ==================
      IF( NCOGEL .EQ. 1 ) THEN
         NBNSOM = 1
         NOSOAR(1,1) = 1
         RETURN
      ENDIF
C
C     TRAITEMENT SELON LE TYPE D'EF : PARTIE SPECIFIQUE DE CHACUN D EUX
C     ******************************************************************
      GOTO ( 110 , 120 , 130 , 140 , 150 , 170 , 150 , 170 , 190 , 200 ,
     &       210 , 220 , 110 , 240 , 120 , 110 , 240 , 120 , 290 , 300 ,
     &       310 , 320 , 330 , 340 , 240 , 170 , 105 , 380 , 110 , 600 ,
     &       410 , 420 , 430 , 440 ) , NO
C     ******************************************************************
C
  105 NBLGRC(NRERR) = 2
      KERR(1) = 'TYPE D''EF NON PROGRAMME'
      WRITE( KERR(3)(1:10), '(I10)' ) NO
      KERR(2) = 'DE NUMERO ' // KERR(3)(1:10)
      CALL LEREUR
      RETURN
C
C     ==================================================================
C     TRIA AP1D  ET  TRIA 2P1D  ET  TRIA 2P1C
C     ==================================================================
  110 NOTFA  = 6
  112 NBNSOM = NARET
      DO 115 I=1,NARET
           NBPOAR(I) = 0
C          NOMBRE DE POINTS-NON SOMMETS DE CHAQUE ARETE
           NBNOAR(I) = 0
C          NOMBRE DE NOEUDS-NON SOMMETS DE CHAQUE ARETE
           NOTYAR(I) = 1
C          NO DU TYPE DE CHAQUE ARETE
  115 CONTINUE
C
C     LE TYPE NOMBRE DE SOMMETS POINTS NOEUDS DE LA FACE
  117 NOTYFA(1) = NOTFA
      NBSOFA(1) = NBNSOM
      NBPOFA(1) = 0
      NBNOFA(1) = 0
      RETURN
C
C     ==================================================================
C     TRIA AP2C  ET  TRIA 2P2C
C     ==================================================================
  120 NOTFA  = 7
  122 NBNSOM = NARET
      DO 125 I=1,NARET
           NBPOAR(I) = 1
C          NOMBRE DE POINTS-NON SOMMETS DE L ARETE I
           NOPOAR(1,I) = I + NARET
C          NO ELEMENTAIRE DU POINT MILIEU DE L ARETE I
           NBNOAR(I) = 1
C          NOMBRE DE NOEUDS-NON SOMMETS DE L ARETE I
           NONOAR(1,I) = I + NARET
C          NO ELEMENTAIRE DU NOEUD MILIEU DE L ARETE I
           NOTYAR(I) = 2
C          NO DU TYPE DE L ARETE I (TYPE P2)
 125  CONTINUE
C     LA FACE
      GOTO 117
C
C     ==================================================================
C     QUAD AQ1C  ET  QUAD 2Q1C
C     ==================================================================
 130  NOTFA = 4
      GOTO 112
C
C     ==================================================================
C     QUAD AQ2C  ET  QUAD 2Q2C    LE BARYCENTRE N EST NI POINT NI NOEUD
C     ==================================================================
 140  NOTFA = 5
      GOTO 122
C
C     ==================================================================
C     TRIA MT10  ET  QUAD MQ10
C     ==================================================================
  150 NBNSOM = 0
      DO 155 I=1,NARET
           NBPOAR(I)   = 0
           NBNOAR(I)   = 1
           NONOAR(1,I) = I
           NOTYAR(I)   = 3
  155 CONTINUE
C     LA FACE EST A METTRE A JOUR
      RETURN
C
C     ==================================================================
C     TRIA MT21  ET  QUAD MQ21 ET TRIA EQ06
C     ==================================================================
  170 NBNSOM = 0
      DO 175 I=1,NARET
           NBPOAR(I)   = 0
           NBNOAR(I)   = 2
           NONOAR(1,I) = 2 * I - 1
           NONOAR(2,I) = 2 * I
           NOTYAR(I)   = 4
  175 CONTINUE
C     LA FACE EST A METTRE A JOUR
      RETURN
C
C     ==================================================================
C     TETR M3T1
C     ==================================================================
  190 NBNFA  = 1
      NOTFA  = 1
  191 NBSFA  = 3
C
  193 NBNSOM = 0
      DO 195 I=1,NARET
           NBPOAR(I) = 0
           NBNOAR(I) = 0
           NOTYAR(I) = 1
  195 CONTINUE
      L      = 0
      DO 197 I=1,NFACE
           NBSOFA(I) = NBSFA
           NBPOFA(I) = 0
           NOTYFA(I) = NOTFA
           NBNOFA(I) = NBNFA
           DO 196 J=1,NBNFA
                L       = L + 1
                NONOFA(J,I) = L
  196      CONTINUE
  197 CONTINUE
      RETURN
C
C     ==================================================================
C     TETR M3T2
C     ==================================================================
  200 NBNFA  = 3
      NOTFA  = 2
      GOTO 191
C
C     ==================================================================
C     HEXA M3H1
C     ==================================================================
  210 NBNFA  = 1
      NOTFA  = 4
  211 NBSFA  = 4
      GOTO 193
C
C     ==================================================================
C     HEXA M3H2
C     ==================================================================
  220 NBNFA  = 4
      NOTFA  = 5
      GOTO 211
C
C     ==================================================================
C     TRIA 2P2D   ET    QUAD 2P2D  ET  TRIA HD06
C     ==================================================================
  240 NBNSOM = NARET
      DO 245 I=1,NARET
           NBPOAR(I) = 0
           NBNOAR(I) = 1
           NONOAR(1,I) = I + NARET
           NOTYAR(I) = 1
  245 CONTINUE
C     LE TYPE NOMBRE DE SOMMETS POINTS NOEUDS DE LA FACE
      IF( NARET .EQ. 3 ) THEN
         NOTFA = 7
      ELSE
         NOTFA = 4
      ENDIF
      NOTYFA(1) = NOTFA
      NBSOFA(1) = NBNSOM
      NBPOFA(1) = 0
      NBNOFA(1) = 0
      RETURN
C
C     ==================================================================
C     TETR 3P1D
C     ==================================================================
  290 NBNSOM = 4
C
C     TYPE ARETE P1
C     -------------
      NOTAR  = 1
      NBPAR  = 0
      NBNAR  = 0
C
C     TYPE FACE TRIANGULAIRE P1
C     -------------------------
      NBSFA  = 3
      NOTFA  = 6
C
C     LES ARETES
C     ----------
  292 DO 294 I=1,NARET
           NOTYAR(I) = NOTAR
           NBPOAR(I) = NBPAR
           NBNOAR(I) = NBNAR
  294 CONTINUE
C
C     LES FACES
C     ---------
      DO 296 I=1,NFACE
           NBSOFA(I) = NBSFA
           NOTYFA(I) = NOTFA
           NBPOFA(I) = 0
           NBNOFA(I) = 0
  296 CONTINUE
      RETURN
C
C     ==================================================================
C     TETR 3P2C
C     ==================================================================
  300 NBNSOM = 4
C
C     TYPE FACE TRIANGULAIRE P2
C     -------------------------
      NBSFA = 3
      NOTFA = 7
C
C     TYPE ARETE P2
C     -------------
  301 NOTAR = 2
      NBPAR = 1
      NBNAR = 1
      DO 304 I=1,NARET
           NOPOAR(1,I) = NBNSOM + I
           NONOAR(1,I) = NBNSOM + I
  304 CONTINUE
      GOTO 292
C
C     ==================================================================
C     PENT 3R1C
C     ==================================================================
  310 NBNSOM = 6
C
C     ARETE TYPE P1
C     -------------
      NOTAR  = 1
      NBPAR  = 0
      NBNAR  = 0
C
C     2 TYPES DE FACE : TRIANGLE P1  OU  QUADRANGLE Q1
C     -----------------------------------------------
      NOTYFA(1) = 6
      NOTYFA(2) = 8
      NOTYFA(3) = 8
      NOTYFA(4) = 6
      NOTYFA(5) = 8
C
C     LES ARETES
C     ----------
  312 DO 314 I=1,NARET
           NOTYAR(I) = NOTAR
           NBPOAR(I) = NBPAR
           NBNOAR(I) = NBNAR
  314 CONTINUE
C
C     LES FACES
C     ---------
      NBSOFA(1) = 3
      NBSOFA(2) = 4
      NBSOFA(3) = 4
      NBSOFA(4) = 3
      NBSOFA(5) = 4
C
      DO 318 I=1,NFACE
           NBPOFA(I) = 0
           NBNOFA(I) = 0
  318 CONTINUE
      RETURN
C
C     ==================================================================
C     PENT 3R2C
C     ==================================================================
  320 NBNSOM = 6
C
C     ARETE DE TYPE P2
C     ----------------
      NOTAR  = 2
      NBPAR  = 1
      NBNAR  = 1
C
C     2 TYPES DE FACE TRIANGLE P2 :7 ET QUADRANGLE Q2 :9
C     --------------------------------------------------
      NOTYFA(1) = 7
      NOTYFA(2) = 9
      NOTYFA(3) = 9
      NOTYFA(4) = 7
      NOTYFA(5) = 9
C
C     LE NO ELEMENTAIRE DES POINTS=NOEUDS DES ARETES
C     ----------------------------------------------
      DO 322 I=1,NARET
           NOPOAR(1,I) = NBNSOM + I
           NONOAR(1,I) = NBNSOM + I
  322 CONTINUE
      GOTO 312
C
C     ==================================================================
C     HEXA 3Q1C
C     ==================================================================
  330 NBNSOM = 8
C
C     ARETE DE TYPE P1
C     ----------------
      NOTAR  = 1
      NBPAR  = 0
      NBNAR  = 0
C
C     FACE DE TYPE QUADRANGULAIRE Q1
C     -------------------------------
      NBSFA  = 4
      NOTFA  = 8
      GOTO 292
C
C     ==================================================================
C     HEXA 3Q2C
C     ==================================================================
 340  NBNSOM = 8
C
C     ARETE DE TYPE P2
C     ----------------
C     FACES DE TYPE QUADRANGLE Q2
C     ---------------------------
      NBSFA  = 4
      NOTFA  = 9
      GOTO 301
C
C     ==================================================================
C     6CUB 6Q1C     reduit a un HEXA 3Q1C
C     ==================================================================
 600  NBNSOM = 64
C
C     ARETE DE TYPE P1
C     ----------------
      NOTAR  = 1
      NBPAR  = 0
      NBNAR  = 0
C
C     FACE DE TYPE QUADRANGULAIRE Q1
C     ------------------------------
C     NBSFA = 2*5 = 32 DANS LA REALITE DU 6-CUBE!
      NBSFA  = 4
      NOTFA  = 8
      GOTO 292
C
C     ==================================================================
C     PYRA 3PY1
C     ==================================================================
 410  NBNSOM = 5
C
C     ARETE TYPE P1
C     -------------
      NOTAR  = 1
      NBPAR  = 0
      NBNAR  = 0
C
C     2 TYPES DE FACE : TRIANGLE P1:6  OU  QUADRANGLE Q1:8
C     ----------------------------------------------------
      NOTYFA(1) = 8
      NOTYFA(2) = 6
      NOTYFA(3) = 6
      NOTYFA(4) = 6
      NOTYFA(5) = 6
C
C     LES ARETES
C     ----------
 412  DO 414 I=1,NARET
         NOTYAR(I) = NOTAR
         NBPOAR(I) = NBPAR
         NBNOAR(I) = NBNAR
 414  CONTINUE
C
C     LES FACES
C     ---------
      NBSOFA(1) = 4
      NBSOFA(2) = 3
      NBSOFA(3) = 3
      NBSOFA(4) = 3
      NBSOFA(5) = 3
C
C     LES POINTS=NOEUDS INTERNES A CHAQUE FACE
      DO 418 I=1,NFACE
         NBPOFA(I) = 0
         NBNOFA(I) = 0
 418  CONTINUE
      RETURN
C
C     ==================================================================
C     PYRA 3PY2
C     ==================================================================
 420  NBNSOM = 5
C
C     ARETE DE TYPE P2
C     ----------------
      NOTAR  = 2
      NBPAR  = 1
      NBNAR  = 1
C
C     2 TYPES DE FACE TRIANGLE P2 :7 ET QUADRANGLE Q2 :9
C     --------------------------------------------------
      NOTYFA(1) = 9
      NOTYFA(2) = 7
      NOTYFA(3) = 7
      NOTYFA(4) = 7
      NOTYFA(5) = 7
C
C     LE NO ELEMENTAIRE DES POINTS=NOEUDS DES ARETES
C     ----------------------------------------------
      DO 422 I=1,NARET
         NOPOAR(1,I) = NBNSOM + I
         NONOAR(1,I) = NBNSOM + I
 422  CONTINUE
      GOTO 412
C
C     ==================================================================
C     SEGM 1P1D
C     ==================================================================
C     NOMBRE DE NOEUDS SOMMETS
 380  NOTFA  = 0
      NBNSOM = 2
C     NOMBRE DE POINTS-NON SOMMETS DE CHAQUE ARETE
      NBPOAR(1) = 0
C     NOMBRE DE NOEUDS-NON SOMMETS DE CHAQUE ARETE
      NBNOAR(1) = 0
C     NO DU TYPE DE CHAQUE ARETE  P1-LAGRANGE
      NOTYAR(1) = 1
      GOTO 445
C
C     ==================================================================
C     SEGM 1P2D
C     ==================================================================
C     NOMBRE DE NOEUDS SOMMETS
 430  NOTFA  = 0
      NBNSOM = 2
C     NOMBRE DE POINTS-NON SOMMETS DE CHAQUE ARETE
      NBPOAR(1) = 1
C     NOMBRE DE NOEUDS-NON SOMMETS DE CHAQUE ARETE
      NBNOAR(1) = 1
C     NO DU TYPE DE CHAQUE ARETE  P2-LAGRANGE
      NOTYAR(1) = 2
C     NO LOCAL DU NOEUD MILIEU DE L'ARETE
      NOPOAR(1,1) = 3
C     NO LOCAL DU POINT MILIEU DE L'ARETE
      NONOAR(1,1) = 3
      GOTO 445
C
C     ==================================================================
C     SEGM 1P3D
C     ==================================================================
C     NOMBRE DE NOEUDS SOMMETS
 440  NOTFA  = 0
      NBNSOM = 2
C     NOMBRE DE POINTS-NON SOMMETS DE CHAQUE ARETE
      NBPOAR(1) = 0
C     NOMBRE DE NOEUDS-NON SOMMETS DE CHAQUE ARETE
      NBNOAR(1) = 0
C     NO DU TYPE DE CHAQUE ARETE  P3-HERMITE
      NOTYAR(1) = 3
C
C     LE TYPE NOMBRE DE SOMMETS POINTS NOEUDS DE LA FACE
 445  NOTYFA(1) = 0
      NBSOFA(1) = 0
      NBPOFA(1) = 0
      NBNOFA(1) = 0
C
      RETURN
      END
