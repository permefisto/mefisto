      SUBROUTINE PIVMTR( NOTYEF, NOINTI, NBPIEX, NOPIEX, COPIEX )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT     : RECUPERER LE NUMERO D'INTERPOLATION, LES POINTS
C --------- D'INTEGRATION UTILISES POUR EXTRAPOLER LE CRITERE DE
C           VON MISES-TRESCA AU NIVEAU DES NOEUDS EN FONCTION
C           DU TYPE D'EF
C
C ENTREES :
C ---------
C NOTYEF  : NUMERO DU TYPE D'ELEMENT FINI
C
C SORTIES:
C --------
C NOINTI  : NUMERO D'INTERPOLATION
C NBPIEX  : NOMBRE DE  POINTS D'INTEGRATION NECESSAIRES
C NOPIEX  : NUMERO DES POINTS D'INTEGRATION UTILISES
C COPIEX  : COORDONNEES DES NBPIEX POINTS D'INTEGRATION
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : Alain PERRONNET Laboratoire J-L. LIONS UPMC PARIS    MAI 2007
C23456---------------------------------------------------------------012
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      DOUBLE PRECISION  COPIEX(3,8), RAC, D, DD
      INTEGER           NOPIEX(30)
C
      GOTO ( 1001, 1002, 1030, 1040,   10,   10, 10,   10,   10,   10,
     &         10,   10, 1013,   10, 1002, 1030, 10, 1040, 1019, 1020,
     &       1021, 1022, 1023, 1024,   10,   10, 10,   10, 1001,   10)
     &                                                           ,NOTYEF
C
C     ==================================================================
C     TRIA AP1D
C     ==================================================================
 1001 NOINTI=2001
      NBPIEX=3
      NOPIEX(1)=1
      NOPIEX(2)=1
      NOPIEX(3)=1
C
C     EN FAIT LE POINT D'INTEGRATION 1 de TRP2
      COPIEX(1,1)  = 1.D0 / 6.D0
      COPIEX(2,1)  = COPIEX(1,1)
      COPIEX(3,1)  = 0D0
C
C     EN FAIT LE POINT D'INTEGRATION 2 de TRP2
      COPIEX(1,2)  = 2.D0 / 3.D0
      COPIEX(2,2)  = COPIEX(1,1)
      COPIEX(3,2)  = 0D0
C
C     EN FAIT LE POINT D'INTEGRATION 3 de TRP2
      COPIEX(1,3)  = COPIEX(1,1)
      COPIEX(2,3)  = COPIEX(1,2)
      COPIEX(3,3)  = 0D0
      RETURN
C
C     ==================================================================
C     TRIA AP2C ET TRIA 2P2C
C     ==================================================================
 1002 NOINTI=2001
      NBPIEX=3
      NOPIEX(1)=1
      NOPIEX(2)=2
      NOPIEX(3)=3
C
C     EN FAIT LE POINT D'INTEGRATION 1 de TRP5
      RAC = DSQRT(15.D0)
      COPIEX(1,1)  = (6.D0 - RAC) / 21.D0
      COPIEX(2,1)  = COPIEX(1,1)
      COPIEX(3,1)  = 0D0
C
C     EN FAIT LE POINT D'INTEGRATION 2 de TRP5
      COPIEX(1,2)  = (9.D0 + 2.D0 * RAC) / 21.D0
      COPIEX(2,2)  = COPIEX(1,1)
      COPIEX(3,2)  = 0D0
C
C     EN FAIT LE POINT D'INTEGRATION 3 de TRP5
      COPIEX(1,3)  = COPIEX(1,1)
      COPIEX(2,3)  = COPIEX(1,2)
      COPIEX(3,3)  = 0D0
      RETURN
C
C     ==================================================================
C     TRIA 2P1D
C     ==================================================================
 1013 NOINTI=2001
      NBPIEX=3
      NOPIEX(1)=1
      NOPIEX(2)=1
      NOPIEX(3)=1
C
C     EN FAIT LE SOMMET 1
      COPIEX(1,1)  = 0D0
      COPIEX(2,1)  = 0D0
      COPIEX(3,1)  = 0D0
C
C     EN FAIT LE SOMMET 2
      COPIEX(1,2)  = 1D0
      COPIEX(2,2)  = 0D0
      COPIEX(3,2)  = 0D0
C
C     EN FAIT LE SOMMET 3
      COPIEX(1,3)  = 0D0
      COPIEX(2,3)  = 1D0
      COPIEX(3,3)  = 0D0
      RETURN
C
C     ==================================================================
C     QUAD AQ1C ET QUAD 2Q1C
C     ==================================================================
 1030 NOINTI=2031
      NBPIEX=4
C     ATTENTION: ORDRE DE L'INTERPOLATION UTILISEE ENSUITE
      NOPIEX(1)=1
      NOPIEX(2)=2
      NOPIEX(3)=4
      NOPIEX(4)=3
C
      RAC = 0.5D0 / DSQRT(3.D0)
      COPIEX(1,1)  = 0.5D0 - RAC
      COPIEX(2,1)  = COPIEX(1,1)
      COPIEX(3,1)  = 0D0
C
      COPIEX(1,2)  = 0.5D0 + RAC
      COPIEX(2,2)  = COPIEX(1,1)
      COPIEX(3,2)  = 0D0
C
      COPIEX(1,3)  = COPIEX(1,1)
      COPIEX(2,3)  = COPIEX(1,2)
      COPIEX(3,3)  = 0D0
C
      COPIEX(1,4)  = COPIEX(1,2)
      COPIEX(2,4)  = COPIEX(1,2)
      COPIEX(3,4)  = 0D0
      RETURN
C
C     ==================================================================
C     QUAD AQ2C ET QUAD 2Q2C
C     ==================================================================
 1040 NOINTI=2031
      NBPIEX=4
C     ATTENTION: ORDRE DE L'INTERPOLATION UTILISEE ENSUITE
      NOPIEX(1)=1
      NOPIEX(2)=3
      NOPIEX(3)=7
      NOPIEX(4)=9
C
      RAC = 0.5D0 * DSQRT(3.D0 / 5.D0)
      COPIEX(1,1) = 0.5D0 - RAC
      COPIEX(2,1) = COPIEX(1,1)
      COPIEX(3,1) = 0D0
C
      COPIEX(1,2) = 0.5D0 + RAC
      COPIEX(2,2) = COPIEX(1,1)
      COPIEX(3,2) = 0D0
C
      COPIEX(1,3) = COPIEX(1,1)
      COPIEX(2,3) = COPIEX(1,2)
      COPIEX(3,3) = 0D0
C
      COPIEX(1,4) = COPIEX(1,2)
      COPIEX(2,4) = COPIEX(1,2)
      COPIEX(3,4) = 0D0
      RETURN
C
C     ==================================================================
C     'TETR','3P1D'
C     ==================================================================
 1019 NOINTI=3001
      NBPIEX=4
      NOPIEX(1)=1
      NOPIEX(2)=1
      NOPIEX(3)=1
      NOPIEX(4)=1
C
C     EN FAIT LE POINT D'INTEGRATION 1 de TEP1
      COPIEX(1,1) = 0.D0
      COPIEX(2,1) = 0.D0
      COPIEX(3,1) = 0.D0
C
C     EN FAIT LE POINT D'INTEGRATION 2 de TEP1
      COPIEX(1,2) = 1.D0
      COPIEX(2,2) = 0.D0
      COPIEX(3,2) = 0.D0
C
C     EN FAIT LE POINT D'INTEGRATION 3 de TEP1
      COPIEX(1,3) = 0.D0
      COPIEX(2,3) = 1.D0
      COPIEX(3,3) = 0.D0
C
C     EN FAIT LE POINT D'INTEGRATION 4 de TEP1
      COPIEX(1,4) = 0.D0
      COPIEX(2,4) = 0.D0
      COPIEX(3,4) = 1.D0
      RETURN
CC
C     ==================================================================
C     'TETR','3P2C'
C     ==================================================================
 1020 NOINTI=3001
      NBPIEX=4
      NOPIEX(1)=2
      NOPIEX(2)=3
      NOPIEX(3)=4
      NOPIEX(4)=5
C
      RAC = DSQRT(15D0)
      D   = (7D0 - RAC)/34D0
      DD  = (13D0 + 3D0 * RAC)/34D0
C
C     EN FAIT LE POINT D'INTEGRATION 2 de TEP5
      COPIEX(1,1)   = D
      COPIEX(2,1)   = D
      COPIEX(3,1)   = D
C
C     EN FAIT LE POINT D'INTEGRATION 3 de TEP5
      COPIEX(1,2)   = DD
      COPIEX(2,2)   = D
      COPIEX(3,2)   = D
C
C     EN FAIT LE POINT D'INTEGRATION 4 de TEP5
      COPIEX(1,3)   = D
      COPIEX(2,3)   = DD
      COPIEX(3,3)   = D
C
C     EN FAIT LE POINT D'INTEGRATION 5 de TEP5
      COPIEX(1,4)   = D
      COPIEX(2,4)   = D
      COPIEX(3,4)   = DD
      RETURN
C
C     ==================================================================
C     'PENT','3R1C'
C     ==================================================================
 1021 NOINTI=3031
      NBPIEX=6
      NOPIEX(1)=1
      NOPIEX(2)=2
      NOPIEX(3)=3
      NOPIEX(4)=4
      NOPIEX(5)=5
      NOPIEX(6)=6
C
C     EN FAIT LE POINT D'INTEGRATION 1 de PER1
      COPIEX(1,1) = 0D0
      COPIEX(2,1) = 0D0
      COPIEX(3,1) = 0D0
C
      COPIEX(1,2) = 1D0
      COPIEX(2,2) = 0D0
      COPIEX(3,2) = 0D0
C
      COPIEX(1,3) = 0D0
      COPIEX(2,3) = 1D0
      COPIEX(3,3) = 0D0
C
      COPIEX(1,4) = 0D0
      COPIEX(2,4) = 0D0
      COPIEX(3,4) = 1D0
C
      COPIEX(1,5) = 1D0
      COPIEX(2,5) = 0D0
      COPIEX(3,5) = 1D0
C
      COPIEX(1,6) = 0D0
      COPIEX(2,6) = 1D0
      COPIEX(3,6) = 1D0
      RETURN
C
C     ==================================================================
C     'PENT','3R2C'
C     ==================================================================
 1022 NOINTI=3031
      NBPIEX=6
      NOPIEX(1)=1
      NOPIEX(2)=2
      NOPIEX(3)=3
      NOPIEX(4)=15
      NOPIEX(5)=16
      NOPIEX(6)=17
C
C     EN FAIT LE POINT D'INTEGRATION 1 de TRP5 * PT 1 de SEP5
      RAC = DSQRT(15.D0)
      COPIEX(1,1)  = (6.D0 - RAC) / 21.D0
      COPIEX(2,1)  = COPIEX(1,1)
      DD = 0.5D0 - 0.5D0 * DSQRT(3.D0 / 5.D0)
      COPIEX(3,1)  = DD
C
C     EN FAIT LE POINT D'INTEGRATION 2 de TRP5 * PT 1 de SEP5
      COPIEX(1,2)  = (9.D0 + 2.D0 * RAC) / 21.D0
      COPIEX(2,2)  = COPIEX(1,1)
      COPIEX(3,2)  = DD
C
C     EN FAIT LE POINT D'INTEGRATION 3 de TRP5 * PT 1 de SEP5
      COPIEX(1,3)  = COPIEX(1,1)
      COPIEX(2,3)  = COPIEX(1,2)
      COPIEX(3,3)  = DD
C
C     EN FAIT LE POINT D'INTEGRATION 1 de TRP5 * PT 3 de SEP5
      COPIEX(1,4)  = (6.D0 - RAC) / 21.D0
      COPIEX(2,4)  = COPIEX(1,1)
      DD = 0.5D0 + 0.5D0 * DSQRT(3.D0 / 5.D0)
      COPIEX(3,4)  = DD
C
C     EN FAIT LE POINT D'INTEGRATION 2 de TRP5 * PT 3 de SEP5
      COPIEX(1,5)  = (9.D0 + 2.D0 * RAC) / 21.D0
      COPIEX(2,5)  = COPIEX(1,1)
      COPIEX(3,5)  = DD
C
C     EN FAIT LE POINT D'INTEGRATION 3 de TRP5 * PT 3 de SEP5
      COPIEX(1,6)  = COPIEX(1,1)
      COPIEX(2,6)  = COPIEX(1,2)
      COPIEX(3,6)  = DD
      RETURN
C
C     ==================================================================
C     'HEXA','3Q1C'
C     ==================================================================
 1023 NOINTI=3061
      NBPIEX=8
      NOPIEX(1)=1
      NOPIEX(2)=2
      NOPIEX(3)=3
      NOPIEX(4)=4
      NOPIEX(5)=5
      NOPIEX(6)=6
      NOPIEX(7)=7
      NOPIEX(8)=8
C
C     EN FAIT LE POINT D'INTEGRATION i de HEQ1 = SOMMET i
      COPIEX(1,1) = 0D0
      COPIEX(2,1) = 0D0
      COPIEX(3,1) = 0D0
C
      COPIEX(1,2) = 1D0
      COPIEX(2,2) = 0D0
      COPIEX(3,2) = 0D0
C
      COPIEX(1,3) = 0D0
      COPIEX(2,3) = 1D0
      COPIEX(3,3) = 0D0
C
      COPIEX(1,4) = 1D0
      COPIEX(2,4) = 1D0
      COPIEX(3,4) = 0D0
C
      COPIEX(1,5) = 0D0
      COPIEX(2,5) = 0D0
      COPIEX(3,5) = 1D0
C
      COPIEX(1,6) = 1D0
      COPIEX(2,6) = 0D0
      COPIEX(3,6) = 1D0
C
      COPIEX(1,7) = 0D0
      COPIEX(2,7) = 1D0
      COPIEX(3,7) = 1D0
C
      COPIEX(1,8) = 1D0
      COPIEX(2,8) = 1D0
      COPIEX(3,8) = 1D0
      RETURN
C
C     ==================================================================
C     'HEXA','3Q2C'
C     ==================================================================
 1024 NOINTI=3061
      NBPIEX=8
      NOPIEX(1)=1
      NOPIEX(2)=3
      NOPIEX(3)=7
      NOPIEX(4)=9
      NOPIEX(5)=19
      NOPIEX(6)=21
      NOPIEX(7)=25
      NOPIEX(8)=27
C
      RAC = 0.5D0 * DSQRT(3.D0 / 5.D0)
      D   = 0.5D0 - RAC
      DD  = 0.5D0 + RAC
C
C     EN FAIT LE POINT D'INTEGRATION 1 de HEQ5
      COPIEX(1,1) = D
      COPIEX(2,1) = D
      COPIEX(3,1) = D
C
C     EN FAIT LE POINT D'INTEGRATION 3 de HEQ5
      COPIEX(1,2) = DD
      COPIEX(2,2) = D
      COPIEX(3,2) = D
C
C     EN FAIT LE POINT D'INTEGRATION 7 de HEQ5
      COPIEX(1,3) = D
      COPIEX(2,3) = DD
      COPIEX(3,3) = D
C
C     EN FAIT LE POINT D'INTEGRATION 9 de HEQ5
      COPIEX(1,4) = DD
      COPIEX(2,4) = DD
      COPIEX(3,4) = D
C
C     EN FAIT LE POINT D'INTEGRATION 19 de HEQ5
      COPIEX(1,5) = D
      COPIEX(2,5) = D
      COPIEX(3,5) = DD
C
C     EN FAIT LE POINT D'INTEGRATION 21 de HEQ5
      COPIEX(1,6) = DD
      COPIEX(2,6) = D
      COPIEX(3,6) = DD
C
C     EN FAIT LE POINT D'INTEGRATION 25 de HEQ5
      COPIEX(1,7) = D
      COPIEX(2,7) = DD
      COPIEX(3,7) = DD
C
C     EN FAIT LE POINT D'INTEGRATION 27 de HEQ5
      COPIEX(1,8) = DD
      COPIEX(2,8) = DD
      COPIEX(3,8) = DD
C
 10   RETURN
      END
