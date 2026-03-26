      SUBROUTINE ELCONO( NO,
     &                   NDIM, NBNOE, DCOORN )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : RETOURNER SELON LE NOM DE L ELEMENT FINI
C ----- LE NOMBRE DE NOEUDS DE L ELEMENT FINI
C       LE NOMBRE DE COORDONNEES DES NOEUDS (1, 2, 3 ou 6)
C       LES COORDONNEES DE CHAQUE NOEUD DE L ELEMENT FINI DE REFERENCE
C
C ENTREE :
C --------
C NO     : NO DE TYPE DE L ELEMENT FINI
C
C SORTIES:
C --------
C NDIM   : NOMBRE DE COMPOSANTES DES COORDONNEES
C NBNOE  : NOMBRE DE NOEUDS DE L ELEMENT FINI
C DCOORN : DCOORN(6,NBNOE) LES COORDONNEES REELLES DOUBLE PRECISION
C                          DES NBNOE NOEUDS
C         (SEULES LES NDIM-ERES SONT ICI INITIALISEES)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     JANVIER 1991
C23456---------------------------------------------------------------012
      include"./incl/gsmenu.inc"
      DOUBLE PRECISION  DCOORN(6,*),D1,D2,DSQRT
      CHARACTER*4       NOMELE(2)
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
C
C     ******************************************************************
      GOTO ( 10,  20,  30,  40, 50, 60,  70,  60,  90, 100,
     &      110, 120,  10,  20, 20, 30,  40,  40, 190, 200,
     &      210, 220, 230, 240, 20, 60,   5, 280,  10, 600,
     &      310, 320, 330, 340 ), NO
C     ******************************************************************
C
    5 NBLGRC(NRERR) = 2
      CALL ELNUNM( NO , NOMELE )
      KERR(1) = 'ELCONO: TYPE EF ' // NOMELE(1) // NOMELE(2)
      KERR(2) = 'COORDONNEES DES NOEUDS NON PROGRAMMEES'
      CALL LEREUR
      NO = 0
      RETURN
C
C     ==================================================================
C     TRIA AP1D ET TRIA 2P1D ET TRIA 2P1C
C     ==================================================================
   10 NBNOE = 3
C
   11 DCOORN(1,3) = 0.0D0
      DCOORN(2,3) = 1.0D0
C
   12 DCOORN(1,2) = 1.0D0
      DCOORN(2,2) = 0.0D0
      DCOORN(1,1) = 0.0D0
      DCOORN(2,1) = 0.0D0
      NDIM       = 2
      RETURN
C
C     ==================================================================
C     TRIA AP2C ET TRIA 2P2C ET TRIA HD06
C     ==================================================================
C
   20 NBNOE      = 6
      DCOORN(1,6) = 0.0D0
      DCOORN(2,6) = 0.5D0
      DCOORN(1,5) = 0.5D0
      DCOORN(2,5) = 0.5D0
C
      DCOORN(1,4) = 0.5D0
      DCOORN(2,4) = 0.0D0
      GOTO 11
C
C     ==================================================================
C     QUAD AQ1C ET QUAD 2Q1C
C     ==================================================================
C
   30 NBNOE      = 4
C
   31 DCOORN(1,4) = 0.0D0
      DCOORN(2,4) = 1.0D0
      DCOORN(1,3) = 1.0D0
      DCOORN(2,3) = 1.0D0
      GOTO 12
C
C     ==================================================================
C     QUAD AQ2C ET QUAD 2Q2C
C     ==================================================================
C
   40 NBNOE      = 8
      DCOORN(1,8) = 0.0D0
      DCOORN(2,8) = 0.5D0
      DCOORN(1,7) = 0.5D0
      DCOORN(2,7) = 1.0D0
      DCOORN(1,6) = 1.0D0
      DCOORN(2,6) = 0.5D0
      DCOORN(1,5) = 0.5D0
      DCOORN(2,5) = 0.0D0
      GOTO 31
C
C     ==================================================================
C     TRIA MT10
C     ==================================================================
C
   50 NBNOE      = 3
      DCOORN(1,3) = 0.0D0
      DCOORN(2,3) = 0.5D0
      DCOORN(1,2) = 0.5D0
      DCOORN(2,2) = 0.5D0
   51 DCOORN(1,1) = 0.5D0
      DCOORN(2,1) = 0.0D0
      NDIM       = 2
      RETURN
C
C     ==================================================================
C     TRIA MT21 ET TRIA EQ06
C     ==================================================================
C
   60 D2         = 0.5D0 / DSQRT(3.D0)
      D1         = 0.5D0 - D2
      D2         = 0.5D0 + D2
      IF(NO .EQ. 8) GOTO 80
C
      NBNOE      = 6
      DCOORN(1,6) = 0.0D0
      DCOORN(2,6) = D1
      DCOORN(1,5) = 0.0D0
      DCOORN(2,5) = D2
      DCOORN(1,4) = D1
      DCOORN(2,4) = D2
      DCOORN(1,3) = D2
      DCOORN(2,3) = D1
   61 DCOORN(1,2) = D2
      DCOORN(2,2) = 0.0D0
      DCOORN(1,1) = D1
      DCOORN(2,1) = 0.0D0
      NDIM       = 2
      RETURN
C
C     ==================================================================
C     QUAD MQ10
C     ==================================================================
C
   70 NBNOE      = 4
      DCOORN(1,4) = 0.0D0
      DCOORN(2,4) = 0.5D0
      DCOORN(1,3) = 0.5D0
      DCOORN(2,3) = 1.0D0
      DCOORN(1,2) = 1.0D0
      DCOORN(2,2) = 0.0D0
      GOTO 51
C
C     ==================================================================
C     QUAD MQ21
C     ==================================================================
C
   80 NBNOE      = 8
      DCOORN(1,8) = 0.0D0
      DCOORN(2,8) = D1
      DCOORN(1,7) = 0.0D0
      DCOORN(2,7) = D2
      DCOORN(1,6) = D1
      DCOORN(2,6) = 1.0D0
      DCOORN(1,5) = D2
      DCOORN(2,5) = 1.0D0
      DCOORN(1,4) = 1.0D0
      DCOORN(2,4) = D2
      DCOORN(1,3) = 1.0D0
      DCOORN(2,3) = D1
      GOTO 61
C     ==================================================================
C     TETR M3T1
C     ==================================================================
C
   90 NBNOE      = 4
      NDIM       = 3
      D1         = 1.D0/3.D0
      DCOORN(1,1) = D1
      DCOORN(2,1) = D1
      DCOORN(3,1) = 0D0
      DCOORN(1,2) = 0D0
      DCOORN(2,2) = D1
      DCOORN(3,2) = D1
      DCOORN(1,3) = D1
      DCOORN(2,3) = 0D0
      DCOORN(3,3) = D1
      DCOORN(1,4) = D1
      DCOORN(2,4) = D1
      DCOORN(3,4) = D1
      RETURN
C     ==================================================================
C     TETR M3T2
C     ==================================================================
C
  100 NBNOE      = 12
      NDIM       = 3
      D1         = .5D0
      DCOORN(1,1) = 0D0
      DCOORN(2,1) = D1
      DCOORN(3,1) = 0D0
      DCOORN(1,2) = D1
      DCOORN(2,2) = D1
      DCOORN(3,2) = 0D0
      DCOORN(1,3) = D1
      DCOORN(2,3) = 0D0
      DCOORN(3,3) = 0D0
      DCOORN(1,4) = 0D0
      DCOORN(2,4) = 0D0
      DCOORN(3,4) = D1
      DCOORN(1,5) = 0D0
      DCOORN(2,5) = D1
      DCOORN(3,5) = D1
      DCOORN(1,6) = 0D0
      DCOORN(2,6) = D1
      DCOORN(3,6) = 0D0
      DCOORN(1,7) = D1
      DCOORN(2,7) = 0D0
      DCOORN(3,7) = 0D0
      DCOORN(1,8) = D1
      DCOORN(2,8) = 0D0
      DCOORN(3,8) = D1
      DCOORN(1,9) = 0D0
      DCOORN(2,9) = 0D0
      DCOORN(3,9) = D1
      DCOORN(1,10) = D1
      DCOORN(2,10) = D1
      DCOORN(3,10) = 0D0
      DCOORN(1,11) = 0D0
      DCOORN(2,11) = D1
      DCOORN(3,11) = D1
      DCOORN(1,12) = D1
      DCOORN(2,12) = 0D0
      DCOORN(3,12) = D1
      RETURN
C     ==================================================================
C     HEXA M3H1
C     ==================================================================
C
  110 NBNOE      = 6
      NDIM       = 3
      D1         = .5D0
      DCOORN(1,1) = D1
      DCOORN(2,1) = D1
      DCOORN(3,1) = 0D0
      DCOORN(1,2) = 0D0
      DCOORN(2,2) = D1
      DCOORN(3,2) = D1
      DCOORN(1,3) = D1
      DCOORN(2,3) = 0D0
      DCOORN(3,3) = D1
      DCOORN(1,4) = D1
      DCOORN(2,4) = D1
      DCOORN(3,4) = 1D0
      DCOORN(1,5) = 1D0
      DCOORN(2,5) = D1
      DCOORN(3,5) = D1
      DCOORN(1,6) = D1
      DCOORN(2,6) = 1D0
      DCOORN(3,6) = D1
      RETURN
C     ==================================================================
C     HEXA M3H2
C     ==================================================================
C
  120 NBNOE      = 24
      NDIM       = 3
      D2         = .5D0 / DSQRT(3.D0)
      D1         = .5D0 - D2
      D2         = .5D0 + D2
      DCOORN(1,1) = D1
      DCOORN(2,1) = D1
      DCOORN(3,1) = 0D0
      DCOORN(1,2) = D1
      DCOORN(2,2) = D2
      DCOORN(3,2) = 0D0
      DCOORN(1,3) = D2
      DCOORN(2,3) = D2
      DCOORN(3,3) = 0D0
      DCOORN(1,4) = D2
      DCOORN(2,4) = D1
      DCOORN(3,4) = 0D0
      DCOORN(1,5) = 0D0
      DCOORN(2,5) = D1
      DCOORN(3,5) = D1
      DCOORN(1,6) = 0D0
      DCOORN(2,6) = D1
      DCOORN(3,6) = D2
      DCOORN(1,7) = 0D0
      DCOORN(2,7) = D2
      DCOORN(3,7) = D2
      DCOORN(1,8) = 0D0
      DCOORN(2,8) = D2
      DCOORN(3,8) = D1
      DCOORN(1,9) = D1
      DCOORN(2,9) = 0D0
      DCOORN(3,9) = D1
      DCOORN(1,10) = D2
      DCOORN(2,10) = 0D0
      DCOORN(3,10) = D1
      DCOORN(1,11) = D2
      DCOORN(2,11) = 0D0
      DCOORN(3,11) = D2
      DCOORN(1,12) = D1
      DCOORN(2,12) = 0D0
      DCOORN(3,12) = D2
      DCOORN(1,13) = D1
      DCOORN(2,13) = D1
      DCOORN(3,13) = 1D0
      DCOORN(1,14) = D2
      DCOORN(2,14) = D1
      DCOORN(3,14) = 1D0
      DCOORN(1,15) = D2
      DCOORN(2,15) = D2
      DCOORN(3,15) = 1D0
      DCOORN(1,16) = D1
      DCOORN(2,16) = D2
      DCOORN(3,16) = 1D0
      DCOORN(1,17) = 1D0
      DCOORN(2,17) = D1
      DCOORN(3,17) = D1
      DCOORN(1,18) = 1D0
      DCOORN(2,18) = D2
      DCOORN(3,18) = D1
      DCOORN(1,19) = 1D0
      DCOORN(2,19) = D2
      DCOORN(3,19) = D2
      DCOORN(1,20) = 1D0
      DCOORN(2,20) = D1
      DCOORN(3,20) = D2
      DCOORN(1,21) = D1
      DCOORN(2,21) = 1D0
      DCOORN(3,21) = D1
      DCOORN(1,22) = D1
      DCOORN(2,22) = 1D0
      DCOORN(3,22) = D2
      DCOORN(1,23) = D2
      DCOORN(2,23) = 1D0
      DCOORN(3,23) = D2
      DCOORN(1,24) = D2
      DCOORN(2,24) = 1D0
      DCOORN(3,24) = D1
      RETURN
C
C     ==================================================================
C     TETR 3P1D
C     ==================================================================
C
  190 NBNOE  = 4
C
  195 NDIM   = 3
C
      DCOORN(1,1) = 0D0
      DCOORN(2,1) = 0D0
      DCOORN(3,1) = 0D0
C
      DCOORN(1,2) = 1D0
      DCOORN(2,2) = 0D0
      DCOORN(3,2) = 0D0
C
      DCOORN(1,3) = 0D0
      DCOORN(2,3) = 1D0
      DCOORN(3,3) = 0D0
C
      DCOORN(1,4) = 0D0
      DCOORN(2,4) = 0D0
      DCOORN(3,4) = 1D0
      RETURN
C
C     ==================================================================
C     TETR 3P2C
C     ==================================================================
C
  200 NBNOE  = 10
C
      DCOORN(1, 5) = 0.5D0
      DCOORN(2, 5) =   0D0
      DCOORN(3, 5) =   0D0
C
      DCOORN(1, 6) = 0.5D0
      DCOORN(2, 6) = 0.5D0
      DCOORN(3, 6) =   0D0
      DCOORN(1, 7) =   0D0
      DCOORN(2, 7) = 0.5D0
      DCOORN(3, 7) =   0D0
C
      DCOORN(1, 8) =   0D0
      DCOORN(2, 8) =   0D0
      DCOORN(3, 8) = 0.5D0
C
      DCOORN(1, 9) = 0.5D0
      DCOORN(2, 9) =   0D0
      DCOORN(3, 9) = 0.5D0
C
      DCOORN(1,10) =   0D0
      DCOORN(2,10) = 0.5D0
      DCOORN(3,10) = 0.5D0
      GOTO 195
C
C     ==================================================================
C     PENT 3R1C
C     ==================================================================
C
  210 NBNOE  = 6
C
  215 NDIM   = 3
C
      DCOORN(1,1) = 0D0
      DCOORN(2,1) = 0D0
      DCOORN(3,1) = 0D0
C
      DCOORN(1,2) = 1D0
      DCOORN(2,2) = 0D0
      DCOORN(3,2) = 0D0
C
      DCOORN(1,3) = 0D0
      DCOORN(2,3) = 1D0
      DCOORN(3,3) = 0D0
C
      DCOORN(1,4) = 0D0
      DCOORN(2,4) = 0D0
      DCOORN(3,4) = 1D0
C
      DCOORN(1,5) = 1D0
      DCOORN(2,5) = 0D0
      DCOORN(3,5) = 1D0
C
      DCOORN(1,6) = 0D0
      DCOORN(2,6) = 1D0
      DCOORN(3,6) = 1D0
C
      RETURN
C
C     ==================================================================
C     PENT 3R2C
C     ==================================================================
C
  220 NBNOE  = 15
C
      DCOORN(1, 7) = 0.5D0
      DCOORN(2, 7) =   0D0
      DCOORN(3, 7) =   0D0
C
      DCOORN(1, 8) = 0.5D0
      DCOORN(2, 8) = 0.5D0
      DCOORN(3, 8) =   0D0
C
      DCOORN(1, 9) =   0D0
      DCOORN(2, 9) = 0.5D0
      DCOORN(3, 9) =   0D0
C
      DCOORN(1,10) =   0D0
      DCOORN(2,10) =   0D0
      DCOORN(3,10) = 0.5D0
C
      DCOORN(1,11) = 1.0D0
      DCOORN(2,11) =   0D0
      DCOORN(3,11) = 0.5D0
C
      DCOORN(1,12) =   0D0
      DCOORN(2,12) = 1.0D0
      DCOORN(3,12) = 0.5D0
      DCOORN(1,13) = 0.5D0
      DCOORN(2,13) =   0D0
      DCOORN(3,13) = 1.0D0
C
      DCOORN(1,14) = 0.5D0
      DCOORN(2,14) = 0.5D0
      DCOORN(3,14) = 1.0D0
C
      DCOORN(1,15) =   0D0
      DCOORN(2,15) = 0.5D0
      DCOORN(3,15) = 1.0D0
      GOTO 215
C
C     ==================================================================
C     HEXA 3Q1C
C     ==================================================================
C
  230 NBNOE  = 8
C
  235 NDIM   = 3
C
      DCOORN(1,1) = 0D0
      DCOORN(2,1) = 0D0
      DCOORN(3,1) = 0D0
C
      DCOORN(1,2) = 1D0
      DCOORN(2,2) = 0D0
      DCOORN(3,2) = 0D0
C
      DCOORN(1,3) = 1D0
      DCOORN(2,3) = 1D0
      DCOORN(3,3) = 0D0
C
      DCOORN(1,4) = 0D0
      DCOORN(2,4) = 1D0
      DCOORN(3,4) = 0D0
C
      DCOORN(1,5) = 0D0
      DCOORN(2,5) = 0D0
      DCOORN(3,5) = 1D0
C
      DCOORN(1,6) = 1D0
      DCOORN(2,6) = 0D0
      DCOORN(3,6) = 1D0
C
      DCOORN(1,7) = 1D0
      DCOORN(2,7) = 1D0
      DCOORN(3,7) = 1D0
C
      DCOORN(1,8) = 0D0
      DCOORN(2,8) = 1D0
      DCOORN(3,8) = 1D0
C
      RETURN
C     ==================================================================
C     HEXA 3Q2C
C     ==================================================================
C
  240 NBNOE  = 20
C
      DCOORN(1, 9) = 0.5D0
      DCOORN(2, 9) =   0D0
      DCOORN(3, 9) =   0D0
C
      DCOORN(1,10) = 1.0D0
      DCOORN(2,10) = 0.5D0
      DCOORN(3,10) =   0D0
C
      DCOORN(1,11) = 0.5D0
      DCOORN(2,11) = 1.0D0
      DCOORN(3,11) =   0D0
C
      DCOORN(1,12) =   0D0
      DCOORN(2,12) = 0.5D0
      DCOORN(3,12) =   0D0
C
      DCOORN(1,13) =   0D0
      DCOORN(2,13) =   0D0
      DCOORN(3,13) = 0.5D0
C
      DCOORN(1,14) = 1.0D0
      DCOORN(2,14) =   0D0
      DCOORN(3,14) = 0.5D0
C
      DCOORN(1,15) = 1.0D0
      DCOORN(2,15) = 1.0D0
      DCOORN(3,15) = 0.5D0
C
      DCOORN(1,16) =   0D0
      DCOORN(2,16) = 1.0D0
      DCOORN(3,16) = 0.5D0
C
      DCOORN(1,17) = 0.5D0
      DCOORN(2,17) =   0D0
      DCOORN(3,17) = 1.0D0
C
      DCOORN(1,18) = 1.0D0
      DCOORN(2,18) = 0.5D0
      DCOORN(3,18) = 1.0D0
C
      DCOORN(1,19) = 0.5D0
      DCOORN(2,19) = 1.0D0
      DCOORN(3,19) = 1.0D0
C
      DCOORN(1,20) =   0D0
      DCOORN(2,20) = 0.5D0
      DCOORN(3,20) = 1.0D0
      GOTO 235
C
C     ==================================================================
C     6CUB 6Q1C  GENERALISATION A 6 DIMENSIONS DU Q1C
C     ==================================================================
C
 600  NBNOE = 64
      NDIM  = 6
C
C     CALCUL DES 6 COORDONNEES DES NBSOM SOMMETS=NOEUDS
      NO = 0
      DO 660 N=0,1
         DO 650 M=0,1
            DO 640 L=0,1
               DO 630 K=0,1
                  DO 620 J=0,1
                     DO 610 I=0,1
                        NO = NO + 1
                        DCOORN(1,NO) = I
                        DCOORN(2,NO) = J
                        DCOORN(3,NO) = K
                        DCOORN(4,NO) = L
                        DCOORN(5,NO) = M
                        DCOORN(6,NO) = N
 610                 CONTINUE
 620              CONTINUE
 630           CONTINUE
 640        CONTINUE
 650     CONTINUE
 660  CONTINUE
      RETURN
C
C
C     ==================================================================
C     PYRA 3PY1
C     ==================================================================
C
  310 NBNOE  = 5
C
  315 NDIM   = 3
C
      DCOORN(1,1) = 0D0
      DCOORN(2,1) = 0D0
      DCOORN(3,1) = 0D0
C
      DCOORN(1,2) = 1D0
      DCOORN(2,2) = 0D0
      DCOORN(3,2) = 0D0
C
      DCOORN(1,3) = 1D0
      DCOORN(2,3) = 1D0
      DCOORN(3,3) = 0D0
C
      DCOORN(1,4) = 0D0
      DCOORN(2,4) = 1D0
      DCOORN(3,4) = 0D0
C
      DCOORN(1,5) = 0D0
      DCOORN(2,5) = 0D0
      DCOORN(3,5) = 1D0
C
      RETURN
C
C     ==================================================================
C     PYRA 3PY2
C     ==================================================================
C
  320 NBNOE  = 13
C
      DCOORN(1, 6) = 0.5D0
      DCOORN(2, 6) =   0D0
      DCOORN(3, 6) =   0D0
C
      DCOORN(1, 7) =   1D0
      DCOORN(2, 7) = 0.5D0
      DCOORN(3, 7) =   0D0
C
      DCOORN(1, 8) = 0.5D0
      DCOORN(2, 8) =   1D0
      DCOORN(3, 8) =   0D0
C
      DCOORN(1, 9) =   0D0
      DCOORN(2, 9) = 0.5D0
      DCOORN(3, 9) =   0D0
C
      DCOORN(1,10) =   0D0
      DCOORN(2,10) =   0D0
      DCOORN(3,10) = 0.5D0
C
      DCOORN(1,11) = 0.5D0
      DCOORN(2,11) =   0D0
      DCOORN(3,11) = 0.5D0
C
      DCOORN(1,12) = 0.5D0
      DCOORN(2,12) = 0.5D0
      DCOORN(3,12) = 0.5D0
C
      DCOORN(1,13) =   0D0
      DCOORN(2,13) = 0.5D0
      DCOORN(3,13) = 0.5D0
C
      GOTO 315
C
C     ==================================================================
C     SEGM 1P1D
C     ==================================================================
 280  NBNOE = 2
      GOTO 335
C
C     ==================================================================
C     SEGM 1P2D
C     ==================================================================
 330  NBNOE = 3
      DCOORN(1,3) = 0.5D0
C
 335  DCOORN(1,1) = 0.0D0
      DCOORN(1,2) = 1.0D0
      NDIM  = 1
      RETURN
C
C     ==================================================================
C     SEGM 1P3D
C     ==================================================================
 340  NBNOE = 2
      GOTO 335
      END
