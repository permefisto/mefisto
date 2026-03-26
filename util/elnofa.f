      SUBROUTINE ELNOFA ( NO, K, NBNOFK, NONOFK )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    FOURNIR LE NO DE TOUS LES NOEUDS DE L EF DE NUMERO NO
C -----    DE LA FACE K SI L EF EST TRIDIMENSIONNEL
C          DE L ARETE K SI L EF EST  BIDIMENSIONNEL
C
C ATTENTION : NE PAS CONFONDRE AVEC LES RESULTATS DU SP ELTYCA
C ----------- ICI LE NO DE TOUS LES NOEUDS DE LA FACE EST FOURNI
C             DANS ELTYCA SEULS LES NOEUDS INTERNES SONT FOURNIS
C
C ENTREES:
C --------
C NO     : NO DE L EF DANS LES SP UTILITAIRES
C K      : NO DE LA FACE DONT LE NO DES NOEUDS EST DEMANDE
C
C SORTIES:
C --------
C NBNOFK : NOMBRE DE NOEUDS DE LA FACE K
C NONOFK : NO DES NOEUDS DE LA FACE K
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     JANVIER 1991
C23456---------------------------------------------------------------012
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      INTEGER           NONOFK(*)
      CHARACTER*4       NOMELE(2)
C
      INTEGER  NOTEP2(6,4)
      INTEGER  NOPER2(8,5)
      INTEGER  NOPYR2(8,5)
      INTEGER  NOHEQ2(8,6)
C
      DATA     NOTEP2 / 1,3,2,7,6,5,   1,4,3,8,10,7,   1,2,4,5,9,8,
     &                  2,3,4,6,10,9/
C
      DATA     NOPER2 / 1,3,2,9,8,7,0,0,      1,4,6,3,10,15,12,9,
     &                  1,2,5,4,7,11,13,10,   4,5,6,13,14,15,0,0,
     &                  2,3,6,5,8,12,14,11/
C
      DATA     NOHEQ2 / 1,4,3,2,12,11,10, 9,  1,5,8,4,13,20,16,12,
     &                  1,2,6,5, 9,14,17,13,  5,6,7,8,17,18,19,20,
     &                  2,3,7,6,10,15,18,14,  4,8,7,3,16,19,15,11/
C
      DATA     NOPYR2 / 1,4,3,2,9,8,7,6,
     &                  1,2,5,6,11,10,0,0,   2,3,5,7,12,11,0,0,
     &                  3,4,5,8,13,12,0,0,   4,1,5,9,10,13,0,0 /
C
C     TRAITEMENT SELON L EF   PARTIE SPECIFIQUE DE CHACUN D EUX
C     ******************************************************************
      GOTO ( 130, 140, 160, 170,   1,   1,   1,   1,   1,   1 ,
     &         1,   1, 130, 140, 140, 160, 170, 170, 190, 200 ,
     &       210, 220, 230, 240, 140, 250,   1, 330, 130,   1 ,
     &       310, 320, 330, 330 ), NO
C     ******************************************************************
 1    IF( 0 .LT. NO .AND. NO .LT. 30 ) THEN
         CALL ELNUNM( NO, NOMELE )
      ELSE
         NOMELE(1) = 'INCO'
         NOMELE(2) = 'NNU '
      ENDIF
      NBLGRC(NRERR) = 3
      KERR(1) = 'ELNOFA: EF AVEC DES FACES NON DEFINIES'
      KERR(2) = NOMELE(1) // NOMELE(2)
      KERR(3) = 'TYPE D''EF A PROGRAMMER'
      CALL LEREUR
      RETURN
C
C     ==================================================================
C     TRIA 2P1D   TRIA AP1D   TRIA 2P1C
C     ==================================================================
  130 NBARE  = 3
C
  132 NBNOFK = 2
C
  134 K1     = K + 1
      IF( K1 .GT. NBARE ) K1 = 1
      NONOFK( 1 ) = K
      NONOFK( 2 ) = K1
      RETURN
C
C     ==================================================================
C     TRIA 2P2D TRIA 2P2C ET TRIA HD06
C     ==================================================================
  140 NBARE  = 3
C
  142 NBNOFK = 3
      NONOFK ( 3 ) = K + NBARE
      GOTO 134
C
C     ==================================================================
C     QUAD 2Q1C
C     ==================================================================
  160 NBARE  = 4
      GOTO 132
C
C     ==================================================================
C     QUAD 2Q2C QUAD 2Q2D
C     ==================================================================
  170 NBARE  = 4
      GOTO 142
C
C     ==================================================================
C     TETR 3P1D
C     ==================================================================
  190 NBNOFK = 3
      GOTO 202
C
C     ==================================================================
C     TETR 3P2C
C     ==================================================================
C     TYPE FACE TRIANGULAIRE P2
C     -------------------------
  200 NBNOFK = 6
  202 DO 205 I=1,NBNOFK
         NONOFK(I) = NOTEP2(I,K)
  205 CONTINUE
      RETURN
C
C     ==================================================================
C     PENT 3R1C
C     ==================================================================
  210 GOTO ( 211 , 212 , 212 , 211 , 212 ) , K
C
C     TYPE FACE TRIANGULAIRE P1
C     -------------------------
  211 NBNOFK = 3
      GOTO 224
C
C     TYPE FACE QUADRANGULAIRE Q1
C     ---------------------------
  212 NBNOFK = 4
      GOTO 224
C
C     ==================================================================
C     PENT 3R2C
C     ==================================================================
  220 GOTO ( 221 , 222 , 222 , 221 , 222 ) , K
C
C     TYPE FACE TRIANGULAIRE P2
C     -------------------------
  221 NBNOFK = 6
      GOTO 224
C
C     TYPE FACE QUADRANGULAIRE Q2 SANS BARYCENTRE
C     -------------------------------------------
  222 NBNOFK = 8
  224 DO 225 I=1,NBNOFK
         NONOFK(I) = NOPER2(I,K)
  225 CONTINUE
      RETURN
C
C     ==================================================================
C     HEXA 3Q1C
C     ==================================================================
  230 NBNOFK = 4
      GOTO 242
C
C     ==================================================================
C     HEXA 3Q2C
C     ==================================================================
C
C     TYPE FACE QUADRANGULAIRE Q2 SANS BARYCENTRE
C     -------------------------------------------
  240 NBNOFK = 8
  242 DO 245 I=1,NBNOFK
         NONOFK(I) = NOHEQ2(I,K)
  245 CONTINUE
      RETURN
C
C     ==================================================================
C     TRIA EQ06
C     ==================================================================
  250 NBNOFK=2
      NONOFK(1)=2*K-1
      NONOFK(2)=2*K
      RETURN
C
C     ==================================================================
C     PYRA 3PY1
C     ==================================================================
  310 GOTO ( 312, 311, 311, 311, 311 ) , K
C
C     TYPE FACE TRIANGULAIRE P1
C     -------------------------
  311 NBNOFK = 3
      GOTO 324
C
C     TYPE FACE QUADRANGULAIRE Q1
C     ---------------------------
  312 NBNOFK = 4
      GOTO 324
C
C     ==================================================================
C     PYRA 3PY2
C     ==================================================================
  320 GOTO ( 322, 321, 321, 321, 321 ) , K
C
C     TYPE FACE TRIANGULAIRE P2
C     -------------------------
  321 NBNOFK = 6
      GOTO 324
C
C     TYPE FACE QUADRANGULAIRE Q2 SANS BARYCENTRE
C     -------------------------------------------
  322 NBNOFK = 8
  324 DO 325 I=1,NBNOFK
         NONOFK(I) = NOPYR2(I,K)
  325 CONTINUE
      RETURN
C
C     ==================================================================
C     SEGM 1P1D 1P2D 1P3D
C     ==================================================================
  330 NBARE  = 1
      NBNOFK = 0
      NONOFK( 1 ) = 0
      RETURN
      END
