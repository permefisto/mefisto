      SUBROUTINE GFACE ( NDIM, NCOGEL, K,  ALFA )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : ETABLIR LE PARAMETRAGE DE LA K-EME FACE D UN ELEMENT ET LE
C ----- TABLEAU ALFA(NDIM,NDIM) TEL QUE LA RESTRICTION DE G DE
C       L APPLICATION F DE L ELEMENT DE REFERENCE
C       TRIDIMENSIONNEL A LA FACE K S ECRIVE
C       G( PSI , ETA ) = F( ALFA(1,1) + ALFA(1,2) * PSI + ALFA(1,3)*ETA,
C                           ALFA(2,1) + ALFA(2,2) * PSI + ALFA(2,3)*ETA,
C                           ALFA(3,1) + ALFA(3,2) * PSI + ALFA(3,3)*ETA)
C       BIDIMENSIONNEL A L'ARETE K S ECRIVE
C       G(    PSI    ) = F( ALFA(1,1) + ALFA(1,2) * PSI ,
C                           ALFA(2,1) + ALFA(2,2) * PSI )
C
C PARAMETRES D ENTREE :
C ---------------------
C NDIM   : DIMENSION DE L ESPACE DE L ELEMENT (R ** 2 OU R ** 3)
C NCOGEL : CODE GEOMETRIQUE DE L ELEMENT (3:TRIANGLE,...)
C K      : NO DE LA FACE DANS L ELEMENT
C
C PARAMETRE RESULTAT :
C --------------------
C ALFA   : MATRICE (NDIM,NDIM)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     OCTOBRE 1981
C ......................................................................
      include"./incl/gsmenu.inc"
      DOUBLE PRECISION  ALFA(NDIM,NDIM)
      COMMON / UNITES / LECTEU , IMPRIM , NUNITE(30)
C
C     MISE A ZERO DE ALFA
C     ===================
      CALL AZEROD( NDIM*NDIM, ALFA )
C
      GOTO( 10 , 10 , 300 , 400 , 500 , 600 , 700 , 10, 900 ) , NCOGEL
C
C     ELEMENT NON TRAITE
C     ==================
 10   NBLGRC(NRERR) = 1
      WRITE(KERR(MXLGER)(1:4),'(I4)') NCOGEL
      KERR(1) = 'GFACE:ELEMENT DE TYPE'//KERR(MXLGER)(1:4)
     %           //' NON TRAITE'
      CALL LEREUR
      RETURN
C
C     TRIANGLE ARETE (12) (23) (31)
C     =============================
 300  GOTO( 410 , 320 , 440 ) , K
C
C     G(PSI) = F(1-PSI,PSI)
C     ---------------------
 320  ALFA(1,1) =  1D0
      ALFA(1,2) = -1D0
      ALFA(2,2) =  1D0
      GOTO 1000
C
C     QUADRANGLE ARETE (12) (23) (34) (41)
C     ====================================
 400  GOTO( 410 , 420 , 430 , 440 ) , K
C
C     G(PSI) = F(PSI,0)
C     -----------------
 410  ALFA(1,2) =  1D0
      GOTO 1000
C
C     G(PSI) = F(1,PSI)
C     -----------------
 420  ALFA(1,1) =  1D0
      ALFA(2,2) =  1D0
      GOTO 1000
C
C     G(PSI) = F(1-PSI,1)
C     -------------------
 430  ALFA(1,1) =  1D0
      ALFA(1,2) = -1D0
      ALFA(2,1) =  1D0
      GOTO 1000
C
C     G(PSI) = F(0,1-PSI)
C     -------------------
 440  ALFA(2,1) =  1D0
      ALFA(2,2) = -1D0
      GOTO 1000
C
C     TETRAEDRE  FACES (132) (143) (124) (234)
C     ========================================
 500  GOTO ( 710 , 720 , 730 , 540 ) , K
C
C     G(PSI,ETA) = F(1-PSI-ETA,PSI,ETA)
C     ---------------------------------
 540  ALFA(1,1) =  1D0
      ALFA(1,2) = -1D0
      ALFA(1,3) = -1D0
      ALFA(2,2) =  1D0
      ALFA(3,3) =  1D0
      GOTO 1000
C
C     PENTAEDRE  FACES (132) (1463) (1254) (456) (2365)
C     =================================================
 600  GOTO ( 710 , 720 , 730 , 740 , 650 ) , K
C
C     G(PSI,ETA) = F(1-PSI,PSI,ETA)
C     -----------------------------
 650  ALFA(1,1) =  1D0
      ALFA(1,2) = -1D0
      GOTO 755
C
C     HEXAEDRE  FACES (1432) (1584) (1265) (5678) (2376) (4873)
C     =========================================================
 700  GOTO ( 710 , 720 , 730 , 740 , 750 , 760 ) , K
C
C     G(PSI,ETA) = F(ETA,PSI,0)
C     -------------------------
 710  ALFA(1,3) =  1D0
      ALFA(2,2) =  1D0
      GOTO 1000
C
C     G(PSI,ETA) = F(0,ETA,PSI)
C     -------------------------
 720  ALFA(2,3) =  1D0
      ALFA(3,2) =  1D0
      GOTO 1000
C
C     G(PSI,ETA) = F(PSI,0,ETA)
C     -------------------------
 730  ALFA(1,2) =  1D0
      ALFA(3,3) =  1D0
      GOTO 1000
C
C     G(PSI,ETA) = F(PSI,ETA,1)
C     -------------------------
 740  ALFA(1,2) =  1D0
      ALFA(2,3) =  1D0
      ALFA(3,1) =  1D0
      GOTO 1000
C
C     G(PSI,ETA) = F(1,PSI,ETA)
C     -------------------------
 750  ALFA(1,1) =  1D0
 755  ALFA(2,2) =  1D0
      ALFA(3,3) =  1D0
      GOTO 1000
C
C     G(PSI,ETA) = F(ETA,1,PSI)
C     -------------------------
 760  ALFA(1,3) =  1D0
      ALFA(2,1) =  1D0
      ALFA(3,2) =  1D0
      GOTO 1000
C
C     PYRAMIDE  FACES (1432) (125) (235) (345) (154)
C     ==============================================
 900  GOTO ( 910, 920, 930, 940, 950 ) , K
C
C     G(PSI,ETA) = F(ETA,PSI,0)  FACE 1432
C     -------------------------
 910  ALFA(1,3) =  1D0
      ALFA(2,2) =  1D0
      GOTO 1000
C
C     G(PSI,ETA) = F(PSI,0,ETA)  FACE 125
C     -------------------------
 920  ALFA(1,2) =  1D0
      ALFA(3,3) =  1D0
      GOTO 1000
C
C     G(PSI,ETA) = F(1-ETA,PSI,ETA)  FACE 235
C     -----------------------------
 930  ALFA(1,1) =  1D0
      ALFA(1,3) = -1D0
      ALFA(2,2) =  1D0
      ALFA(3,3) =  1D0
      GOTO 1000
C
C     G(PSI,ETA) = F(ETA,1-PSI,PSI)  FACE 453
C     -----------------------------
 940  ALFA(1,3) =  1D0
      ALFA(2,1) =  1D0
      ALFA(2,2) = -1D0
      ALFA(3,2) =  1D0
      GOTO 1000
C
C     G(PSI,ETA) = F(0,ETA,PSI)  FACE 154
C     -------------------------
 950  ALFA(2,3) =  1D0
      ALFA(3,2) =  1D0
      GOTO 1000
C
 1000 RETURN
      END
