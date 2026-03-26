      SUBROUTINE HACREN( NCOGEL, NBS, NBT0, NUSOM )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : PERMUTER LES NUMEROS DES NBS SOMMETS DE L'EF DE CODE GEOMETRIQUE
C ----- NCOGEL EN METTANT EN PREMIER LE PLUS PETIT NUMERO DE SOMMET ET
C       EN CONSERVANT UN VOLUME POSITIF
C       LES EVENTUELLES TANGENTES SUBISSENT LA MEME PERMUTATION
C
C ENTREES :
C ---------
C NCOGEL : CODE DE GEOMETRIE DE L'ELEMENT
C          1:NOEUD 2:SEGMENT 3:TRIANGLE 4:QUADRANGLE 5:TETRAEDRE
C          6:PENTAEDRE 7:HEXAEDRE 9:PYRAMIDE
C NBS    : NOMBRE DE SOMMETS A TRAITER
C NBT0   : NUMERO DANS NUSOM DE LA PREMIERE TANGENTE A TRAITER
C          0 SI PAS DE TANGENTE A TRAITER
C          LE NOMBRE DE TANGENTES TRAITEES EST IMPOSE SELON NCOGEL
C           2 PAR SEGMENT,   6 PAR TRIANGLE,   8 PAR QUADRANGLE,
C          12 PAR TETRAEDRE,18 PAR PENTAEDRE, 24 PAR HEXAEDRE
C           0 PAR PYRAMIDE
C
C ENTREE ET SORTIE :
C ------------------
C NUSOM  : NUMERO DES NBS SOMMETS ET DES TANGENTES A PARTIR DE NBT0
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS       AVRIL 1989
C ......................................................................
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      INTEGER  NUSOM(1:*)
C
C     TRI SELON LE CODE GEOMETRIQUE
C     =============================
      GOTO( 9000, 20, 30, 30, 30, 30, 30, 30, 30 ) , NCOGEL
C
C     SEGMENT
C     -------
 20   IF( NUSOM(1) .GT. NUSOM(2) ) THEN
         K        = NUSOM(1)
         NUSOM(1) = NUSOM(2)
         NUSOM(2) = K
         IF( NBT0 .GT. 0 ) THEN
C           LES 2 TANGENTES SONT A PERMUTER
            NBT1        = NBT0 + 1
            K           = NUSOM(NBT0)
            NUSOM(NBT0) = NUSOM(NBT1)
            NUSOM(NBT1) = K
         ENDIF
      ENDIF
      GOTO 9000
c
C     RECHERCHE DU PLUS PETIT DES NBS NUMEROS DE SOMMETS
 30   NOSENS = 1
      DO 310 J=2,NBS
         IF( NUSOM(NOSENS) .GT. NUSOM(J) ) THEN
             NOSENS = J
         ENDIF
 310  CONTINUE
C
C     LE PLUS PETIT NUMERO DE SOMMET EST POSITIONNE EN TETE
      GOTO( 9000, 9000, 300, 400, 500, 600, 700, 700, 900 ) , NCOGEL
C
C     TRIANGLE
C     --------
 300  GOTO( 9000, 312, 313 ) , NOSENS
C
C     LE 2-EME SOMMET DEVIENT LE 1-ER : 231 => 123
 312  L        = NUSOM(1)
      NUSOM(1) = NUSOM(2)
      NUSOM(2) = NUSOM(3)
      NUSOM(3) = L
      IF( NBT0 .GT. 0 ) THEN
C        6 TANGENTES A PERMUTER CIRCULAIREMENT A PARTIR DE NBT0
         NT1 = NUSOM(NBT0  )
         NT2 = NUSOM(NBT0+1)
         NUSOM(NBT0  ) = NUSOM(NBT0+2)
         NUSOM(NBT0+1) = NUSOM(NBT0+3)
         NUSOM(NBT0+2) = NUSOM(NBT0+4)
         NUSOM(NBT0+3) = NUSOM(NBT0+5)
         NUSOM(NBT0+4) = NT1
         NUSOM(NBT0+5) = NT2
      ENDIF
      GOTO 9000
C
C     LE 3-EME SOMMET DEVIENT LE 1-ER : 312 => 123
 313  L        = NUSOM(1)
      NUSOM(1) = NUSOM(3)
      NUSOM(3) = NUSOM(2)
      NUSOM(2) = L
      IF( NBT0 .GT. 0 ) THEN
C        6 TANGENTES A PERMUTER CIRCULAIREMENT A PARTIR DE NBT0
         NT1 = NUSOM(NBT0  )
         NT2 = NUSOM(NBT0+1)
         NUSOM(NBT0  ) = NUSOM(NBT0+4)
         NUSOM(NBT0+1) = NUSOM(NBT0+5)
         NUSOM(NBT0+4) = NUSOM(NBT0+2)
         NUSOM(NBT0+5) = NUSOM(NBT0+3)
         NUSOM(NBT0+2) = NT1
         NUSOM(NBT0+3) = NT2
      ENDIF
      GOTO 9000
C
C     QUADRANGLE
C     ----------
 400  GOTO( 9000 , 412 , 413 , 414 ) , NOSENS
C
C     LE 2-EME SOMMET DEVIENT LE 1-ER : 2341 => 1234
 412  L        = NUSOM(4)
      NUSOM(4) = NUSOM(1)
      NUSOM(1) = NUSOM(2)
      NUSOM(2) = NUSOM(3)
      NUSOM(3) = L
      IF( NBT0 .GT. 0 ) THEN
C        8 TANGENTES A PERMUTER CIRCULAIREMENT A PARTIR DE NBT0
         NT1 = NUSOM(NBT0  )
         NT2 = NUSOM(NBT0+1)
         NUSOM(NBT0  ) = NUSOM(NBT0+2)
         NUSOM(NBT0+1) = NUSOM(NBT0+3)
         NUSOM(NBT0+2) = NUSOM(NBT0+4)
         NUSOM(NBT0+3) = NUSOM(NBT0+5)
         NUSOM(NBT0+4) = NUSOM(NBT0+6)
         NUSOM(NBT0+5) = NUSOM(NBT0+7)
         NUSOM(NBT0+6) = NT1
         NUSOM(NBT0+7) = NT2
      ENDIF
      GOTO 9000
C
C     LE 3-EME SOMMET DEVIENT LE 1-ER : 3412 => 1234
 413  L        = NUSOM(1)
      NUSOM(1) = NUSOM(3)
      NUSOM(3) = L
      L        = NUSOM(2)
      NUSOM(2) = NUSOM(4)
      NUSOM(4) = L
      IF( NBT0 .GT. 0 ) THEN
C        8 TANGENTES A PERMUTER CIRCULAIREMENT A PARTIR DE NBT0
         NT1 = NUSOM(NBT0  )
         NT2 = NUSOM(NBT0+1)
         NUSOM(NBT0  ) = NUSOM(NBT0+4)
         NUSOM(NBT0+1) = NUSOM(NBT0+5)
         NUSOM(NBT0+4) = NT1
         NUSOM(NBT0+5) = NT2
         NT1 = NUSOM(NBT0+2)
         NT2 = NUSOM(NBT0+3)
         NUSOM(NBT0+2) = NUSOM(NBT0+6)
         NUSOM(NBT0+3) = NUSOM(NBT0+7)
         NUSOM(NBT0+6) = NT1
         NUSOM(NBT0+7) = NT2
      ENDIF
      GOTO 9000
C
C     LE 4-EME SOMMET DEVIENT LE 1-ER : 4123 => 1234
 414  L        = NUSOM(1)
      NUSOM(1) = NUSOM(4)
      NUSOM(4) = NUSOM(3)
      NUSOM(3) = NUSOM(2)
      NUSOM(2) = L
      IF( NBT0 .GT. 0 ) THEN
C        8 TANGENTES A PERMUTER CIRCULAIREMENT A PARTIR DE NBT0
         NT1 = NUSOM(NBT0  )
         NT2 = NUSOM(NBT0+1)
         NUSOM(NBT0  ) = NUSOM(NBT0+6)
         NUSOM(NBT0+1) = NUSOM(NBT0+7)
         NUSOM(NBT0+6) = NUSOM(NBT0+4)
         NUSOM(NBT0+7) = NUSOM(NBT0+5)
         NUSOM(NBT0+4) = NUSOM(NBT0+2)
         NUSOM(NBT0+5) = NUSOM(NBT0+3)
         NUSOM(NBT0+2) = NT1
         NUSOM(NBT0+3) = NT2
      ENDIF
      GOTO 9000
C
C     TETRAEDRE
C     ---------
C     LE PLUS PETIT NUMERO DE SOMMET EST PORTE EN TETE
 500  GOTO( 590 , 520 , 530 , 540 ) , NOSENS
C
C     TETRAEDRE : NO TANGENTE1(S1S2), NT2(S1S3),  NT3(S1S4),
C                 NO TANGENTE4(S2S3), NT5(S2S1),  NT6(S2S4),
C                 NO TANGENTE7(S3S1), NT8(S3S2),  NT9(S3S4),
C                 NO TANGENT10(S4S1), NT11(S4S2), NT12(S4S3)
C
C     LE 2-EME SOMMET DEVIENT LE 1-ER : 2314 => 1234
 520  L        = NUSOM(3)
      NUSOM(3) = NUSOM(1)
      NUSOM(1) = NUSOM(2)
      NUSOM(2) = L
      IF( NBT0 .GT. 0 ) THEN
C        12 TANGENTES A PERMUTER CIRCULAIREMENT A PARTIR DE NBT0
         NT1 = NUSOM(NBT0  )
         NT2 = NUSOM(NBT0+1)
         NT3 = NUSOM(NBT0+2)
         NUSOM(NBT0  ) = NUSOM(NBT0+3)
         NUSOM(NBT0+1) = NUSOM(NBT0+4)
         NUSOM(NBT0+2) = NUSOM(NBT0+5)
         NUSOM(NBT0+3) = NUSOM(NBT0+6)
         NUSOM(NBT0+4) = NUSOM(NBT0+7)
         NUSOM(NBT0+5) = NUSOM(NBT0+8)
         NUSOM(NBT0+6) = NT1
         NUSOM(NBT0+7) = NT2
         NUSOM(NBT0+8) = NT3
         NT3 = NUSOM(NBT0+9)
         NUSOM(NBT0+ 9) = NUSOM(NBT0+10)
         NUSOM(NBT0+10) = NUSOM(NBT0+11)
         NUSOM(NBT0+11) = NT3
      ENDIF
      GOTO 590
C
C     LE 3-EME SOMMET DEVIENT LE 1-ER : 3124 => 1234
 530  L        = NUSOM(2)
      NUSOM(2) = NUSOM(1)
      NUSOM(1) = NUSOM(3)
      NUSOM(3) = L
      IF( NBT0 .GT. 0 ) THEN
C        12 TANGENTES A PERMUTER CIRCULAIREMENT A PARTIR DE NBT0
         NT1 = NUSOM(NBT0  )
         NT2 = NUSOM(NBT0+1)
         NT3 = NUSOM(NBT0+2)
         NUSOM(NBT0  ) = NUSOM(NBT0+6)
         NUSOM(NBT0+1) = NUSOM(NBT0+7)
         NUSOM(NBT0+2) = NUSOM(NBT0+8)
         NUSOM(NBT0+3) = NUSOM(NBT0+3)
         NUSOM(NBT0+4) = NUSOM(NBT0+4)
         NUSOM(NBT0+5) = NUSOM(NBT0+5)
         NUSOM(NBT0+6) = NT1
         NUSOM(NBT0+7) = NT2
         NUSOM(NBT0+8) = NT3
         NT3 = NUSOM(NBT0+9)
         NUSOM(NBT0+ 9) = NUSOM(NBT0+11)
         NUSOM(NBT0+11) = NUSOM(NBT0+10)
         NUSOM(NBT0+10) = NT3
      ENDIF
      GOTO 590
C
C     LE 4-EME SOMMET DEVIENT LE 1-ER : 4132 => 1234
 540  L        = NUSOM(2)
      NUSOM(2) = NUSOM(1)
      NUSOM(1) = NUSOM(4)
      NUSOM(4) = L
      IF( NBT0 .GT. 0 ) THEN
C        12 TANGENTES A PERMUTER CIRCULAIREMENT A PARTIR DE NBT0
         NT1 = NUSOM(NBT0  )
         NT2 = NUSOM(NBT0+1)
         NT3 = NUSOM(NBT0+2)
         NUSOM(NBT0   ) = NUSOM(NBT0+ 9)
         NUSOM(NBT0+ 1) = NUSOM(NBT0+10)
         NUSOM(NBT0+ 2) = NUSOM(NBT0+11)
         NUSOM(NBT0+ 9) = NUSOM(NBT0+5)
         NUSOM(NBT0+10) = NUSOM(NBT0+4)
         NUSOM(NBT0+11) = NUSOM(NBT0+3)
         NUSOM(NBT0+ 3) = NT2
         NUSOM(NBT0+ 4) = NT3
         NUSOM(NBT0+ 5) = NT1
         NT3 = NUSOM(NBT0+6)
         NUSOM(NBT0+ 6) = NUSOM(NBT0+8)
         NUSOM(NBT0+ 8) = NUSOM(NBT0+7)
         NUSOM(NBT0+ 7) = NT3
      ENDIF
C
C     PARMI LES SOMMETS 2 3 4 RECHERCHE DU PLUS PETIT
C     ET POSITIONNEMENT EN SECONDE POSITION TOUT EN CONSERVANT
C     LA REGLE DE NUMEROTATION ( VOLUME >0 )
 590  IF( NUSOM(2) .LT. NUSOM(3) ) THEN
         IF( NUSOM(4) .LT. NUSOM(2) ) THEN
C           S1 < S4 < S2 < S3 : 1342 => 1234
            L        = NUSOM(4)
            NUSOM(4) = NUSOM(3)
            NUSOM(3) = NUSOM(2)
            NUSOM(2) = L
            IF( NBT0 .GT. 0 ) THEN
C              12 TANGENTES A PERMUTER CIRCULAIREMENT A PARTIR DE NBT0
               NT1 = NUSOM(NBT0)
               NUSOM(NBT0  ) = NUSOM(NBT0+2)
               NUSOM(NBT0+2) = NUSOM(NBT0+1)
               NUSOM(NBT0+1) = NT1
               NT1 = NUSOM(NBT0+3)
               NT2 = NUSOM(NBT0+4)
               NT3 = NUSOM(NBT0+5)
               NUSOM(NBT0+3) = NUSOM(NBT0+10)
               NUSOM(NBT0+4) = NUSOM(NBT0+9)
               NUSOM(NBT0+5) = NUSOM(NBT0+11)
               NUSOM(NBT0+ 9) = NUSOM(NBT0+6)
               NUSOM(NBT0+10) = NUSOM(NBT0+8)
               NUSOM(NBT0+11) = NUSOM(NBT0+7)
               NUSOM(NBT0+6) = NT2
               NUSOM(NBT0+7) = NT3
               NUSOM(NBT0+8) = NT1
            ENDIF
C        ELSE
C           S1 < S2 < S3 ET S4 : 1234 => 1234
         ENDIF
      ELSE
         IF( NUSOM(3) .LT. NUSOM(4) ) THEN
C           S1 < S3 < S2 ET S4 : 1342 => 1234
            L        = NUSOM(3)
            NUSOM(3) = NUSOM(4)
            NUSOM(4) = NUSOM(2)
            NUSOM(2) = L
            IF( NBT0 .GT. 0 ) THEN
C              12 TANGENTES A PERMUTER CIRCULAIREMENT A PARTIR DE NBT0
               NT1 = NUSOM(NBT0)
               NUSOM(NBT0  ) = NUSOM(NBT0+1)
               NUSOM(NBT0+1) = NUSOM(NBT0+2)
               NUSOM(NBT0+2) = NT1
               NT1 = NUSOM(NBT0+3)
               NT2 = NUSOM(NBT0+4)
               NT3 = NUSOM(NBT0+5)
               NUSOM(NBT0+3) = NUSOM(NBT0+8)
               NUSOM(NBT0+4) = NUSOM(NBT0+6)
               NUSOM(NBT0+5) = NUSOM(NBT0+7)
               NUSOM(NBT0+6) = NUSOM(NBT0+9)
               NUSOM(NBT0+7) = NUSOM(NBT0+11)
               NUSOM(NBT0+8) = NUSOM(NBT0+10)
               NUSOM(NBT0+ 9) = NT2
               NUSOM(NBT0+10) = NT1
               NUSOM(NBT0+11) = NT3
            ENDIF
         ELSE
C           S1 < S4 < S3 < S2 : 1423 => 1234
            L        = NUSOM(4)
            NUSOM(4) = NUSOM(3)
            NUSOM(3) = NUSOM(2)
            NUSOM(2) = L
            IF( NBT0 .GT. 0 ) THEN
C              12 TANGENTES A PERMUTER CIRCULAIREMENT A PARTIR DE NBT0
               NT1 = NUSOM(NBT0)
               NUSOM(NBT0  ) = NUSOM(NBT0+2)
               NUSOM(NBT0+2) = NUSOM(NBT0+1)
               NUSOM(NBT0+1) = NT1
               NT1 = NUSOM(NBT0+3)
               NT2 = NUSOM(NBT0+4)
               NT3 = NUSOM(NBT0+5)
               NUSOM(NBT0+3) = NUSOM(NBT0+10)
               NUSOM(NBT0+4) = NUSOM(NBT0+ 9)
               NUSOM(NBT0+5) = NUSOM(NBT0+11)
               NUSOM(NBT0+ 9) = NUSOM(NBT0+6)
               NUSOM(NBT0+10) = NUSOM(NBT0+8)
               NUSOM(NBT0+11) = NUSOM(NBT0+7)
               NUSOM(NBT0+6) = NT2
               NUSOM(NBT0+7) = NT3
               NUSOM(NBT0+8) = NT1
            ENDIF
         ENDIF
      ENDIF
      GOTO 9000
C
C     PENTAEDRE
C     ---------
C     LE PLUS PETIT NUMERO DE SOMMET EST PORTE EN TETE
 600  GOTO( 9000 , 620 , 630 , 640 , 650 , 660 ) , NOSENS
C
C     PENTAEDRE : NO TANGENTE1(S1S2), NT2(S1S3),  NT3(S1S4),
C                 NO TANGENTE4(S2S3), NT5(S2S1),  NT6(S2S5),
C                 NO TANGENTE7(S3S1), NT8(S3S2),  NT9(S3S6),
C                 NO TANGENT10(S4S5), NT11(S4S6), NT12(S4S1),
C                 NO TANGENT13(S5S6), NT14(S5S4), NT15(S5S2),
C                 NO TANGENT16(S6S4), NT17(S6S5), NT18(S6S3)
C
C     LE 2-EME SOMMET DEVIENT LE 1-ER  : 312645 => 123456
 620  L        = NUSOM(1)
      NUSOM(1) = NUSOM(2)
      NUSOM(2) = NUSOM(3)
      NUSOM(3) = L
      L        = NUSOM(4)
      NUSOM(4) = NUSOM(5)
      NUSOM(5) = NUSOM(6)
      NUSOM(6) = L
      IF( NBT0 .GT. 0 ) THEN
C        18 TANGENTES A PERMUTER CIRCULAIREMENT A PARTIR DE NBT0
         NT1 = NUSOM(NBT0  )
         NT2 = NUSOM(NBT0+1)
         NT3 = NUSOM(NBT0+2)
         NUSOM(NBT0  ) = NUSOM(NBT0+3)
         NUSOM(NBT0+1) = NUSOM(NBT0+4)
         NUSOM(NBT0+2) = NUSOM(NBT0+5)
         NUSOM(NBT0+3) = NUSOM(NBT0+6)
         NUSOM(NBT0+4) = NUSOM(NBT0+7)
         NUSOM(NBT0+5) = NUSOM(NBT0+8)
         NUSOM(NBT0+6) = NT1
         NUSOM(NBT0+7) = NT2
         NUSOM(NBT0+8) = NT3
         NT1 = NUSOM(NBT0+ 9)
         NT2 = NUSOM(NBT0+10)
         NT3 = NUSOM(NBT0+11)
         NUSOM(NBT0+ 9) = NUSOM(NBT0+12)
         NUSOM(NBT0+10) = NUSOM(NBT0+13)
         NUSOM(NBT0+11) = NUSOM(NBT0+14)
         NUSOM(NBT0+12) = NUSOM(NBT0+15)
         NUSOM(NBT0+13) = NUSOM(NBT0+16)
         NUSOM(NBT0+14) = NUSOM(NBT0+17)
         NUSOM(NBT0+15) = NT1
         NUSOM(NBT0+16) = NT2
         NUSOM(NBT0+17) = NT3
      ENDIF
      GOTO 9000
C
C     LE 3-EME SOMMET DEVIENT LE 1-ER : 231564 => 123456
 630  L        = NUSOM(1)
      NUSOM(1) = NUSOM(3)
      NUSOM(3) = NUSOM(2)
      NUSOM(2) = L
      L        = NUSOM(4)
      NUSOM(4) = NUSOM(6)
      NUSOM(6) = NUSOM(5)
      NUSOM(5) = L
      IF( NBT0 .GT. 0 ) THEN
C        18 TANGENTES A PERMUTER CIRCULAIREMENT A PARTIR DE NBT0
         NT1 = NUSOM(NBT0  )
         NT2 = NUSOM(NBT0+1)
         NT3 = NUSOM(NBT0+2)
         NUSOM(NBT0  ) = NUSOM(NBT0+6)
         NUSOM(NBT0+1) = NUSOM(NBT0+7)
         NUSOM(NBT0+2) = NUSOM(NBT0+8)
         NUSOM(NBT0+6) = NUSOM(NBT0+3)
         NUSOM(NBT0+7) = NUSOM(NBT0+4)
         NUSOM(NBT0+8) = NUSOM(NBT0+5)
         NUSOM(NBT0+3) = NT1
         NUSOM(NBT0+4) = NT2
         NUSOM(NBT0+5) = NT3
         NT1 = NUSOM(NBT0+ 9)
         NT2 = NUSOM(NBT0+10)
         NT3 = NUSOM(NBT0+11)
         NUSOM(NBT0+ 9) = NUSOM(NBT0+15)
         NUSOM(NBT0+10) = NUSOM(NBT0+16)
         NUSOM(NBT0+11) = NUSOM(NBT0+17)
         NUSOM(NBT0+15) = NUSOM(NBT0+12)
         NUSOM(NBT0+16) = NUSOM(NBT0+13)
         NUSOM(NBT0+17) = NUSOM(NBT0+14)
         NUSOM(NBT0+12) = NT1
         NUSOM(NBT0+13) = NT2
         NUSOM(NBT0+14) = NT3
      ENDIF
      GOTO 9000
C
C     LE 4-EME SOMMET DEVIENT LE 1-ER : 465132 => 123456
 640  L        = NUSOM(1)
      NUSOM(1) = NUSOM(4)
      NUSOM(4) = L
      L        = NUSOM(2)
      NUSOM(2) = NUSOM(6)
      NUSOM(6) = L
      L        = NUSOM(3)
      NUSOM(3) = NUSOM(5)
      NUSOM(5) = L
      IF( NBT0 .GT. 0 ) THEN
C        18 TANGENTES A PERMUTER CIRCULAIREMENT A PARTIR DE NBT0
         NT1 = NUSOM(NBT0  )
         NT2 = NUSOM(NBT0+1)
         NT3 = NUSOM(NBT0+2)
         NUSOM(NBT0  ) = NUSOM(NBT0+10)
         NUSOM(NBT0+1) = NUSOM(NBT0+ 9)
         NUSOM(NBT0+2) = NUSOM(NBT0+11)
         NUSOM(NBT0+ 9) = NT2
         NUSOM(NBT0+10) = NT1
         NUSOM(NBT0+11) = NT3
         NT1 = NUSOM(NBT0+3)
         NT2 = NUSOM(NBT0+4)
         NT3 = NUSOM(NBT0+5)
         NUSOM(NBT0+3) = NUSOM(NBT0+16)
         NUSOM(NBT0+4) = NUSOM(NBT0+15)
         NUSOM(NBT0+5) = NUSOM(NBT0+17)
         NUSOM(NBT0+15) = NT2
         NUSOM(NBT0+16) = NT1
         NUSOM(NBT0+17) = NT3
         NT1 = NUSOM(NBT0+6)
         NT2 = NUSOM(NBT0+7)
         NT3 = NUSOM(NBT0+8)
         NUSOM(NBT0+6) = NUSOM(NBT0+13)
         NUSOM(NBT0+7) = NUSOM(NBT0+12)
         NUSOM(NBT0+8) = NUSOM(NBT0+14)
         NUSOM(NBT0+12) = NT2
         NUSOM(NBT0+13) = NT1
         NUSOM(NBT0+14) = NT3
      ENDIF
      GOTO 9000
C
C     LE 5-EME SOMMET DEVIENT LE 1-ER : 546213 => 123456
 650  L        = NUSOM(1)
      NUSOM(1) = NUSOM(5)
      NUSOM(5) = L
      L        = NUSOM(2)
      NUSOM(2) = NUSOM(4)
      NUSOM(4) = L
      L        = NUSOM(3)
      NUSOM(3) = NUSOM(6)
      NUSOM(6) = L
      IF( NBT0 .GT. 0 ) THEN
C        18 TANGENTES A PERMUTER CIRCULAIREMENT A PARTIR DE NBT0
         NT1 = NUSOM(NBT0  )
         NT2 = NUSOM(NBT0+1)
         NT3 = NUSOM(NBT0+2)
         NUSOM(NBT0   ) = NUSOM(NBT0+12)
         NUSOM(NBT0+ 1) = NUSOM(NBT0+13)
         NUSOM(NBT0+ 2) = NUSOM(NBT0+14)
         NUSOM(NBT0+12) = NUSOM(NBT0+6)
         NUSOM(NBT0+13) = NUSOM(NBT0+7)
         NUSOM(NBT0+14) = NUSOM(NBT0+8)
         NUSOM(NBT0+ 6) = NUSOM(NBT0+ 9)
         NUSOM(NBT0+ 7) = NUSOM(NBT0+10)
         NUSOM(NBT0+ 8) = NUSOM(NBT0+11)
         NUSOM(NBT0+ 9) = NUSOM(NBT0+3)
         NUSOM(NBT0+10) = NUSOM(NBT0+4)
         NUSOM(NBT0+11) = NUSOM(NBT0+5)
         NUSOM(NBT0+ 3) = NUSOM(NBT0+15)
         NUSOM(NBT0+ 4) = NUSOM(NBT0+16)
         NUSOM(NBT0+ 5) = NUSOM(NBT0+17)
         NUSOM(NBT0+15) = NT1
         NUSOM(NBT0+16) = NT2
         NUSOM(NBT0+17) = NT3
      ENDIF
      GOTO 9000
C
C     LE 6-EME SOMMET DEVIENT LE 1-ER : 654321 => 123456
 660  L        = NUSOM(1)
      NUSOM(1) = NUSOM(6)
      NUSOM(6) = L
      L        = NUSOM(2)
      NUSOM(2) = NUSOM(5)
      NUSOM(5) = L
      L        = NUSOM(3)
      NUSOM(3) = NUSOM(4)
      NUSOM(4) = L
      IF( NBT0 .GT. 0 ) THEN
C        18 TANGENTES A PERMUTER CIRCULAIREMENT A PARTIR DE NBT0
         NT1 = NUSOM(NBT0  )
         NT2 = NUSOM(NBT0+1)
         NT3 = NUSOM(NBT0+2)
         NUSOM(NBT0   ) = NUSOM(NBT0+15)
         NUSOM(NBT0+ 1) = NUSOM(NBT0+16)
         NUSOM(NBT0+ 2) = NUSOM(NBT0+17)
         NUSOM(NBT0+15) = NUSOM(NBT0+3)
         NUSOM(NBT0+16) = NUSOM(NBT0+4)
         NUSOM(NBT0+17) = NUSOM(NBT0+5)
         NUSOM(NBT0+ 3) = NUSOM(NBT0+ 9)
         NUSOM(NBT0+ 4) = NUSOM(NBT0+10)
         NUSOM(NBT0+ 5) = NUSOM(NBT0+11)
         NUSOM(NBT0+ 9) = NUSOM(NBT0+6)
         NUSOM(NBT0+10) = NUSOM(NBT0+7)
         NUSOM(NBT0+11) = NUSOM(NBT0+8)
         NUSOM(NBT0+ 6) = NUSOM(NBT0+12)
         NUSOM(NBT0+ 7) = NUSOM(NBT0+13)
         NUSOM(NBT0+ 8) = NUSOM(NBT0+14)
         NUSOM(NBT0+12) = NT1
         NUSOM(NBT0+13) = NT2
         NUSOM(NBT0+14) = NT3
      ENDIF
      GOTO 9000
C
C     HEXAEDRE
C     --------
C     LE PLUS PETIT NUMERO DE SOMMET EST PORTE EN TETE
 700  GOTO( 790 , 720 , 730 , 740 , 750 , 760 , 770 , 780 ) , NOSENS
C
C     HEXAEDRE  : NO TANGENTE1(S1S2), NT2(S1S4),  NT3(S1S5),
C                 NO TANGENTE4(S2S3), NT5(S2S1),  NT6(S2S6),
C                 NO TANGENTE7(S3S4), NT8(S3S2),  NT9(S3S7),
C                 NO TANGENT10(S4S1), NT11(S4S3), NT12(S4S8),
C                 NO TANGENT13(S5S6), NT14(S5S8), NT15(S5S1),
C                 NO TANGENT16(S6S7), NT17(S6S5), NT18(S6S2),
C                 NO TANGENT19(S7S8), NT20(S7S6), NT21(S7S3),
C                 NO TANGENT22(S8S5), NT23(S8S7), NT24(S8S4)
C
C     LE 2-EME SOMMET DEVIENT LE 1-ER : 41238567 => 12345678
 720  L        = NUSOM(1)
      NUSOM(1) = NUSOM(2)
      NUSOM(2) = NUSOM(3)
      NUSOM(3) = NUSOM(4)
      NUSOM(4) = L
      L        = NUSOM(5)
      NUSOM(5) = NUSOM(6)
      NUSOM(6) = NUSOM(7)
      NUSOM(7) = NUSOM(8)
      NUSOM(8) = L
      IF( NBT0 .GT. 0 ) THEN
C        24 TANGENTES A PERMUTER CIRCULAIREMENT A PARTIR DE NBT0
         NT1 = NUSOM(NBT0  )
         NT2 = NUSOM(NBT0+1)
         NT3 = NUSOM(NBT0+2)
         NUSOM(NBT0  ) = NUSOM(NBT0+3)
         NUSOM(NBT0+1) = NUSOM(NBT0+4)
         NUSOM(NBT0+2) = NUSOM(NBT0+5)
         NUSOM(NBT0+3) = NUSOM(NBT0+6)
         NUSOM(NBT0+4) = NUSOM(NBT0+7)
         NUSOM(NBT0+5) = NUSOM(NBT0+8)
         NUSOM(NBT0+6) = NUSOM(NBT0+9)
         NUSOM(NBT0+7) = NUSOM(NBT0+10)
         NUSOM(NBT0+8) = NUSOM(NBT0+11)
         NUSOM(NBT0+ 9) = NT1
         NUSOM(NBT0+10) = NT2
         NUSOM(NBT0+11) = NT3
         NT1 = NUSOM(NBT0+12)
         NT2 = NUSOM(NBT0+13)
         NT3 = NUSOM(NBT0+14)
         NUSOM(NBT0+12) = NUSOM(NBT0+15)
         NUSOM(NBT0+13) = NUSOM(NBT0+16)
         NUSOM(NBT0+14) = NUSOM(NBT0+17)
         NUSOM(NBT0+15) = NUSOM(NBT0+18)
         NUSOM(NBT0+16) = NUSOM(NBT0+19)
         NUSOM(NBT0+17) = NUSOM(NBT0+20)
         NUSOM(NBT0+18) = NUSOM(NBT0+21)
         NUSOM(NBT0+19) = NUSOM(NBT0+22)
         NUSOM(NBT0+20) = NUSOM(NBT0+23)
         NUSOM(NBT0+21) = NT1
         NUSOM(NBT0+22) = NT2
         NUSOM(NBT0+23) = NT3
      ENDIF
      GOTO 790
C
C     LE 3-EME SOMMET DEVIENT LE 1-ER : 34127856 => 12345678
 730  L        = NUSOM(1)
      NUSOM(1) = NUSOM(3)
      NUSOM(3) = L
      L        = NUSOM(2)
      NUSOM(2) = NUSOM(4)
      NUSOM(4) = L
      L        = NUSOM(5)
      NUSOM(5) = NUSOM(7)
      NUSOM(7) = L
      L        = NUSOM(6)
      NUSOM(6) = NUSOM(8)
      NUSOM(8) = L
      IF( NBT0 .GT. 0 ) THEN
C        24 TANGENTES A PERMUTER CIRCULAIREMENT A PARTIR DE NBT0
         NT1 = NUSOM(NBT0  )
         NT2 = NUSOM(NBT0+1)
         NT3 = NUSOM(NBT0+2)
         NUSOM(NBT0  ) = NUSOM(NBT0+6)
         NUSOM(NBT0+1) = NUSOM(NBT0+7)
         NUSOM(NBT0+2) = NUSOM(NBT0+8)
         NUSOM(NBT0+6) = NT1
         NUSOM(NBT0+7) = NT2
         NUSOM(NBT0+8) = NT3
         NT1 = NUSOM(NBT0+3)
         NT2 = NUSOM(NBT0+4)
         NT3 = NUSOM(NBT0+5)
         NUSOM(NBT0+3) = NUSOM(NBT0+ 9)
         NUSOM(NBT0+4) = NUSOM(NBT0+10)
         NUSOM(NBT0+5) = NUSOM(NBT0+11)
         NUSOM(NBT0+ 9) = NT1
         NUSOM(NBT0+10) = NT2
         NUSOM(NBT0+11) = NT3
         NT1 = NUSOM(NBT0+12)
         NT2 = NUSOM(NBT0+13)
         NT3 = NUSOM(NBT0+14)
         NUSOM(NBT0+12) = NUSOM(NBT0+18)
         NUSOM(NBT0+13) = NUSOM(NBT0+19)
         NUSOM(NBT0+14) = NUSOM(NBT0+20)
         NUSOM(NBT0+18) = NT1
         NUSOM(NBT0+19) = NT2
         NUSOM(NBT0+20) = NT3
         NT1 = NUSOM(NBT0+15)
         NT2 = NUSOM(NBT0+16)
         NT3 = NUSOM(NBT0+17)
         NUSOM(NBT0+15) = NUSOM(NBT0+21)
         NUSOM(NBT0+16) = NUSOM(NBT0+22)
         NUSOM(NBT0+17) = NUSOM(NBT0+23)
         NUSOM(NBT0+21) = NT1
         NUSOM(NBT0+22) = NT2
         NUSOM(NBT0+23) = NT3
      ENDIF
      GOTO 790
C
C     LE 4-EME SOMMET DEVIENT LE 1-ER : 23416785 => 12345678
 740  L        = NUSOM(1)
      NUSOM(1) = NUSOM(4)
      NUSOM(4) = NUSOM(3)
      NUSOM(3) = NUSOM(2)
      NUSOM(2) = L
      L        = NUSOM(5)
      NUSOM(5) = NUSOM(8)
      NUSOM(8) = NUSOM(7)
      NUSOM(7) = NUSOM(6)
      NUSOM(6) = L
      IF( NBT0 .GT. 0 ) THEN
C        24 TANGENTES A PERMUTER CIRCULAIREMENT A PARTIR DE NBT0
         NT1 = NUSOM(NBT0  )
         NT2 = NUSOM(NBT0+1)
         NT3 = NUSOM(NBT0+2)
         NUSOM(NBT0  ) = NUSOM(NBT0+ 9)
         NUSOM(NBT0+1) = NUSOM(NBT0+10)
         NUSOM(NBT0+2) = NUSOM(NBT0+11)
         NUSOM(NBT0+ 9) = NUSOM(NBT0+6)
         NUSOM(NBT0+10) = NUSOM(NBT0+7)
         NUSOM(NBT0+11) = NUSOM(NBT0+8)
         NUSOM(NBT0+6) = NUSOM(NBT0+3)
         NUSOM(NBT0+7) = NUSOM(NBT0+4)
         NUSOM(NBT0+8) = NUSOM(NBT0+5)
         NUSOM(NBT0+ 9) = NT1
         NUSOM(NBT0+10) = NT2
         NUSOM(NBT0+11) = NT3
         NT1 = NUSOM(NBT0+12)
         NT2 = NUSOM(NBT0+13)
         NT3 = NUSOM(NBT0+14)
         NUSOM(NBT0+12) = NUSOM(NBT0+21)
         NUSOM(NBT0+13) = NUSOM(NBT0+22)
         NUSOM(NBT0+14) = NUSOM(NBT0+23)
         NUSOM(NBT0+21) = NUSOM(NBT0+18)
         NUSOM(NBT0+22) = NUSOM(NBT0+19)
         NUSOM(NBT0+23) = NUSOM(NBT0+20)
         NUSOM(NBT0+18) = NUSOM(NBT0+15)
         NUSOM(NBT0+19) = NUSOM(NBT0+16)
         NUSOM(NBT0+20) = NUSOM(NBT0+17)
         NUSOM(NBT0+15) = NT1
         NUSOM(NBT0+16) = NT2
         NUSOM(NBT0+17) = NT3
      ENDIF
      GOTO 790
C
C     LE 5-EME SOMMET DEVIENT LE 1-ER : 58761432 => 12345678
 750  L        = NUSOM(1)
      NUSOM(1) = NUSOM(5)
      NUSOM(5) = L
      L        = NUSOM(2)
      NUSOM(2) = NUSOM(8)
      NUSOM(8) = L
      L        = NUSOM(3)
      NUSOM(3) = NUSOM(7)
      NUSOM(7) = L
      L        = NUSOM(4)
      NUSOM(4) = NUSOM(6)
      NUSOM(6) = L
      IF( NBT0 .GT. 0 ) THEN
C        24 TANGENTES A PERMUTER CIRCULAIREMENT A PARTIR DE NBT0
         NT1 = NUSOM(NBT0  )
         NT2 = NUSOM(NBT0+1)
         NT3 = NUSOM(NBT0+2)
         NUSOM(NBT0  ) = NUSOM(NBT0+13)
         NUSOM(NBT0+1) = NUSOM(NBT0+12)
         NUSOM(NBT0+2) = NUSOM(NBT0+14)
         NUSOM(NBT0+12) = NT2
         NUSOM(NBT0+13) = NT1
         NUSOM(NBT0+14) = NT3
         NT1 = NUSOM(NBT0+3)
         NT2 = NUSOM(NBT0+4)
         NT3 = NUSOM(NBT0+5)
         NUSOM(NBT0+3) = NUSOM(NBT0+22)
         NUSOM(NBT0+4) = NUSOM(NBT0+21)
         NUSOM(NBT0+5) = NUSOM(NBT0+23)
         NUSOM(NBT0+21) = NT2
         NUSOM(NBT0+22) = NT1
         NUSOM(NBT0+23) = NT3
         NT1 = NUSOM(NBT0+6)
         NT2 = NUSOM(NBT0+7)
         NT3 = NUSOM(NBT0+8)
         NUSOM(NBT0+6) = NUSOM(NBT0+19)
         NUSOM(NBT0+7) = NUSOM(NBT0+18)
         NUSOM(NBT0+8) = NUSOM(NBT0+20)
         NUSOM(NBT0+18) = NT3
         NUSOM(NBT0+19) = NT1
         NUSOM(NBT0+20) = NT2
         NT1 = NUSOM(NBT0+ 9)
         NT2 = NUSOM(NBT0+10)
         NT3 = NUSOM(NBT0+11)
         NUSOM(NBT0+ 9) = NUSOM(NBT0+16)
         NUSOM(NBT0+10) = NUSOM(NBT0+15)
         NUSOM(NBT0+11) = NUSOM(NBT0+17)
         NUSOM(NBT0+15) = NT2
         NUSOM(NBT0+16) = NT1
         NUSOM(NBT0+17) = NT3
      ENDIF
      GOTO 790
C
C     LE 6-EME SOMMET DEVIENT LE 1-ER : 65872143 => 12345678
 760  L        = NUSOM(1)
      NUSOM(1) = NUSOM(6)
      NUSOM(6) = L
      L        = NUSOM(2)
      NUSOM(2) = NUSOM(5)
      NUSOM(5) = L
      L        = NUSOM(3)
      NUSOM(3) = NUSOM(8)
      NUSOM(8) = L
      L        = NUSOM(4)
      NUSOM(4) = NUSOM(7)
      NUSOM(7) = L
      IF( NBT0 .GT. 0 ) THEN
C        24 TANGENTES A PERMUTER CIRCULAIREMENT A PARTIR DE NBT0
         NT1 = NUSOM(NBT0  )
         NT2 = NUSOM(NBT0+1)
         NT3 = NUSOM(NBT0+2)
         NUSOM(NBT0  ) = NUSOM(NBT0+16)
         NUSOM(NBT0+1) = NUSOM(NBT0+15)
         NUSOM(NBT0+2) = NUSOM(NBT0+17)
         NUSOM(NBT0+15) = NT2
         NUSOM(NBT0+16) = NT1
         NUSOM(NBT0+17) = NT3
         NT1 = NUSOM(NBT0+3)
         NT2 = NUSOM(NBT0+4)
         NT3 = NUSOM(NBT0+5)
         NUSOM(NBT0+3) = NUSOM(NBT0+13)
         NUSOM(NBT0+4) = NUSOM(NBT0+12)
         NUSOM(NBT0+5) = NUSOM(NBT0+14)
         NUSOM(NBT0+12) = NT2
         NUSOM(NBT0+13) = NT1
         NUSOM(NBT0+14) = NT3
         NT1 = NUSOM(NBT0+6)
         NT2 = NUSOM(NBT0+7)
         NT3 = NUSOM(NBT0+8)
         NUSOM(NBT0+6) = NUSOM(NBT0+22)
         NUSOM(NBT0+7) = NUSOM(NBT0+21)
         NUSOM(NBT0+8) = NUSOM(NBT0+23)
         NUSOM(NBT0+21) = NT3
         NUSOM(NBT0+22) = NT1
         NUSOM(NBT0+23) = NT2
         NT1 = NUSOM(NBT0+ 9)
         NT2 = NUSOM(NBT0+10)
         NT3 = NUSOM(NBT0+11)
         NUSOM(NBT0+ 9) = NUSOM(NBT0+19)
         NUSOM(NBT0+10) = NUSOM(NBT0+18)
         NUSOM(NBT0+11) = NUSOM(NBT0+20)
         NUSOM(NBT0+18) = NT2
         NUSOM(NBT0+19) = NT1
         NUSOM(NBT0+20) = NT3
      ENDIF
      GOTO 790
C
C     LE 7-EME SOMMET DEVIENT LE 1-ER : 76583214 => 12345678
 770  L        = NUSOM(1)
      NUSOM(1) = NUSOM(7)
      NUSOM(7) = L
      L        = NUSOM(2)
      NUSOM(2) = NUSOM(6)
      NUSOM(6) = L
      L        = NUSOM(3)
      NUSOM(3) = NUSOM(5)
      NUSOM(5) = L
      L        = NUSOM(4)
      NUSOM(4) = NUSOM(8)
      NUSOM(8) = L
      IF( NBT0 .GT. 0 ) THEN
C        24 TANGENTES A PERMUTER CIRCULAIREMENT A PARTIR DE NBT0
         NT1 = NUSOM(NBT0  )
         NT2 = NUSOM(NBT0+1)
         NT3 = NUSOM(NBT0+2)
         NUSOM(NBT0  ) = NUSOM(NBT0+19)
         NUSOM(NBT0+1) = NUSOM(NBT0+18)
         NUSOM(NBT0+2) = NUSOM(NBT0+20)
         NUSOM(NBT0+18) = NT2
         NUSOM(NBT0+19) = NT1
         NUSOM(NBT0+20) = NT3
         NT1 = NUSOM(NBT0+3)
         NT2 = NUSOM(NBT0+4)
         NT3 = NUSOM(NBT0+5)
         NUSOM(NBT0+3) = NUSOM(NBT0+16)
         NUSOM(NBT0+4) = NUSOM(NBT0+15)
         NUSOM(NBT0+5) = NUSOM(NBT0+17)
         NUSOM(NBT0+15) = NT2
         NUSOM(NBT0+16) = NT1
         NUSOM(NBT0+17) = NT3
         NT1 = NUSOM(NBT0+6)
         NT2 = NUSOM(NBT0+7)
         NT3 = NUSOM(NBT0+8)
         NUSOM(NBT0+6) = NUSOM(NBT0+13)
         NUSOM(NBT0+7) = NUSOM(NBT0+12)
         NUSOM(NBT0+8) = NUSOM(NBT0+14)
         NUSOM(NBT0+12) = NT3
         NUSOM(NBT0+13) = NT1
         NUSOM(NBT0+14) = NT2
         NT1 = NUSOM(NBT0+ 9)
         NT2 = NUSOM(NBT0+10)
         NT3 = NUSOM(NBT0+11)
         NUSOM(NBT0+ 9) = NUSOM(NBT0+22)
         NUSOM(NBT0+10) = NUSOM(NBT0+21)
         NUSOM(NBT0+11) = NUSOM(NBT0+23)
         NUSOM(NBT0+21) = NT2
         NUSOM(NBT0+22) = NT1
         NUSOM(NBT0+23) = NT3
      ENDIF
      GOTO 790
C
C     LE 8-EME SOMMET DEVIENT LE 1-ER : 87654321 => 12345678
 780  L        = NUSOM(1)
      NUSOM(1) = NUSOM(8)
      NUSOM(8) = L
      L        = NUSOM(2)
      NUSOM(2) = NUSOM(7)
      NUSOM(7) = L
      L        = NUSOM(3)
      NUSOM(3) = NUSOM(6)
      NUSOM(6) = L
      L        = NUSOM(4)
      NUSOM(4) = NUSOM(5)
      NUSOM(5) = L
      IF( NBT0 .GT. 0 ) THEN
C        24 TANGENTES A PERMUTER CIRCULAIREMENT A PARTIR DE NBT0
         NT1 = NUSOM(NBT0  )
         NT2 = NUSOM(NBT0+1)
         NT3 = NUSOM(NBT0+2)
         NUSOM(NBT0  ) = NUSOM(NBT0+22)
         NUSOM(NBT0+1) = NUSOM(NBT0+21)
         NUSOM(NBT0+2) = NUSOM(NBT0+23)
         NUSOM(NBT0+21) = NT2
         NUSOM(NBT0+22) = NT1
         NUSOM(NBT0+23) = NT3
         NT1 = NUSOM(NBT0+3)
         NT2 = NUSOM(NBT0+4)
         NT3 = NUSOM(NBT0+5)
         NUSOM(NBT0+3) = NUSOM(NBT0+19)
         NUSOM(NBT0+4) = NUSOM(NBT0+18)
         NUSOM(NBT0+5) = NUSOM(NBT0+20)
         NUSOM(NBT0+18) = NT2
         NUSOM(NBT0+19) = NT1
         NUSOM(NBT0+20) = NT3
         NT1 = NUSOM(NBT0+6)
         NT2 = NUSOM(NBT0+7)
         NT3 = NUSOM(NBT0+8)
         NUSOM(NBT0+6) = NUSOM(NBT0+16)
         NUSOM(NBT0+7) = NUSOM(NBT0+15)
         NUSOM(NBT0+8) = NUSOM(NBT0+17)
         NUSOM(NBT0+15) = NT3
         NUSOM(NBT0+16) = NT1
         NUSOM(NBT0+17) = NT2
         NT1 = NUSOM(NBT0+ 9)
         NT2 = NUSOM(NBT0+10)
         NT3 = NUSOM(NBT0+11)
         NUSOM(NBT0+ 9) = NUSOM(NBT0+13)
         NUSOM(NBT0+10) = NUSOM(NBT0+12)
         NUSOM(NBT0+11) = NUSOM(NBT0+14)
         NUSOM(NBT0+12) = NT2
         NUSOM(NBT0+13) = NT1
         NUSOM(NBT0+14) = NT3
      ENDIF
C
C     PARMI LES SOMMETS 2 4 5 RECHERCHE DU PLUS PETIT
C     ET POSITIONNEMENT EN SECONDE POSITION TOUT EN CONSERVANT
C     LA REGLE DE NUMEROTATION ( VOLUME >0 )
 790  IF( NUSOM(2) .LT. NUSOM(4) ) THEN
         IF( NUSOM(5) .LT. NUSOM(2) ) THEN
C           S5 < S2 < S4 :  14852376 => 12345678
            L        = NUSOM(2)
            NUSOM(2) = NUSOM(5)
            NUSOM(5) = NUSOM(4)
            NUSOM(4) = L
            L        = NUSOM(3)
            NUSOM(3) = NUSOM(6)
            NUSOM(6) = NUSOM(8)
            NUSOM(8) = L
            IF( NBT0 .GT. 0 ) THEN
C              24 TANGENTES A PERMUTER CIRCULAIREMENT A PARTIR DE NBT0
               NT1 = NUSOM(NBT0  )
               NT2 = NUSOM(NBT0+1)
               NT3 = NUSOM(NBT0+2)
               NUSOM(NBT0  ) = NT3
               NUSOM(NBT0+1) = NT1
               NUSOM(NBT0+2) = NT2
               NT1 = NUSOM(NBT0+3)
               NT2 = NUSOM(NBT0+4)
               NT3 = NUSOM(NBT0+5)
               NUSOM(NBT0+3) = NUSOM(NBT0+12)
               NUSOM(NBT0+4) = NUSOM(NBT0+14)
               NUSOM(NBT0+5) = NUSOM(NBT0+13)
               NUSOM(NBT0+12) = NUSOM(NBT0+11)
               NUSOM(NBT0+13) = NUSOM(NBT0+10)
               NUSOM(NBT0+14) = NUSOM(NBT0+ 9)
               NUSOM(NBT0+ 9) = NT2
               NUSOM(NBT0+10) = NT3
               NUSOM(NBT0+11) = NT1
               NT1 = NUSOM(NBT0+6)
               NT2 = NUSOM(NBT0+7)
               NT3 = NUSOM(NBT0+8)
               NUSOM(NBT0+6) = NUSOM(NBT0+17)
               NUSOM(NBT0+7) = NUSOM(NBT0+16)
               NUSOM(NBT0+8) = NUSOM(NBT0+15)
               NUSOM(NBT0+15) = NUSOM(NBT0+22)
               NUSOM(NBT0+16) = NUSOM(NBT0+23)
               NUSOM(NBT0+17) = NUSOM(NBT0+21)
               NUSOM(NBT0+21) = NT1
               NUSOM(NBT0+22) = NT3
               NUSOM(NBT0+23) = NT2
            ENDIF
C        ELSE
C           S2 < S4 ET S5 : 12345678 => 12345678
         ENDIF
      ELSE
         IF( NUSOM(4) .LT. NUSOM(5) ) THEN
C           S4 < S2 ET S5 : 15624873 => 12345678
            L        = NUSOM(2)
            NUSOM(2) = NUSOM(4)
            NUSOM(4) = NUSOM(5)
            NUSOM(5) = L
            L        = NUSOM(3)
            NUSOM(3) = NUSOM(8)
            NUSOM(8) = NUSOM(6)
            NUSOM(6) = L
            IF( NBT0 .GT. 0 ) THEN
C              24 TANGENTES A PERMUTER CIRCULAIREMENT A PARTIR DE NBT0
               NT1 = NUSOM(NBT0  )
               NT2 = NUSOM(NBT0+1)
               NT3 = NUSOM(NBT0+2)
               NUSOM(NBT0  ) = NT2
               NUSOM(NBT0+1) = NT3
               NUSOM(NBT0+2) = NT1
               NT1 = NUSOM(NBT0+3)
               NT2 = NUSOM(NBT0+4)
               NT3 = NUSOM(NBT0+5)
               NUSOM(NBT0+3) = NUSOM(NBT0+11)
               NUSOM(NBT0+4) = NUSOM(NBT0+ 9)
               NUSOM(NBT0+5) = NUSOM(NBT0+10)
               NUSOM(NBT0+ 9) = NUSOM(NBT0+14)
               NUSOM(NBT0+10) = NUSOM(NBT0+13)
               NUSOM(NBT0+11) = NUSOM(NBT0+12)
               NUSOM(NBT0+12) = NT1
               NUSOM(NBT0+13) = NT3
               NUSOM(NBT0+14) = NT2
               NT1 = NUSOM(NBT0+6)
               NT2 = NUSOM(NBT0+7)
               NT3 = NUSOM(NBT0+8)
               NUSOM(NBT0+6) = NUSOM(NBT0+21)
               NUSOM(NBT0+7) = NUSOM(NBT0+23)
               NUSOM(NBT0+8) = NUSOM(NBT0+22)
               NUSOM(NBT0+21) = NUSOM(NBT0+17)
               NUSOM(NBT0+22) = NUSOM(NBT0+15)
               NUSOM(NBT0+23) = NUSOM(NBT0+16)
               NUSOM(NBT0+15) = NT3
               NUSOM(NBT0+16) = NT2
               NUSOM(NBT0+17) = NT1
            ENDIF
         ELSE
C           S5 < S4 < S2 :  14852376 => 12345678
            L        = NUSOM(2)
            NUSOM(2) = NUSOM(5)
            NUSOM(5) = NUSOM(4)
            NUSOM(4) = L
            L        = NUSOM(3)
            NUSOM(3) = NUSOM(6)
            NUSOM(6) = NUSOM(8)
            NUSOM(8) = L
            IF( NBT0 .GT. 0 ) THEN
C              24 TANGENTES A PERMUTER CIRCULAIREMENT A PARTIR DE NBT0
               NT1 = NUSOM(NBT0  )
               NT2 = NUSOM(NBT0+1)
               NT3 = NUSOM(NBT0+2)
               NUSOM(NBT0  ) = NT3
               NUSOM(NBT0+1) = NT1
               NUSOM(NBT0+2) = NT2
               NT1 = NUSOM(NBT0+3)
               NT2 = NUSOM(NBT0+4)
               NT3 = NUSOM(NBT0+5)
               NUSOM(NBT0+3) = NUSOM(NBT0+12)
               NUSOM(NBT0+4) = NUSOM(NBT0+14)
               NUSOM(NBT0+5) = NUSOM(NBT0+13)
               NUSOM(NBT0+12) = NUSOM(NBT0+11)
               NUSOM(NBT0+13) = NUSOM(NBT0+10)
               NUSOM(NBT0+14) = NUSOM(NBT0+ 9)
               NUSOM(NBT0+ 9) = NT2
               NUSOM(NBT0+10) = NT3
               NUSOM(NBT0+11) = NT1
               NT1 = NUSOM(NBT0+6)
               NT2 = NUSOM(NBT0+7)
               NT3 = NUSOM(NBT0+8)
               NUSOM(NBT0+6) = NUSOM(NBT0+17)
               NUSOM(NBT0+7) = NUSOM(NBT0+16)
               NUSOM(NBT0+8) = NUSOM(NBT0+15)
               NUSOM(NBT0+15) = NUSOM(NBT0+22)
               NUSOM(NBT0+16) = NUSOM(NBT0+23)
               NUSOM(NBT0+17) = NUSOM(NBT0+21)
               NUSOM(NBT0+21) = NT1
               NUSOM(NBT0+22) = NT3
               NUSOM(NBT0+23) = NT2
            ENDIF
         ENDIF
      ENDIF
      GOTO 9000
C
C     PYRAMIDE
C     --------
 900  GOTO( 9000, 912, 913, 914, 915 ) , NOSENS
C
C     LE 2-EME SOMMET DEVIENT LE 1-ER : 23415 => 12345
 912  L        = NUSOM(4)
      NUSOM(4) = NUSOM(1)
      NUSOM(1) = NUSOM(2)
      NUSOM(2) = NUSOM(3)
      NUSOM(3) = L
      IF( NBT0 .GT. 0 )  GOTO 9900
C     TANGENTES A PERMUTER CIRCULAIREMENT A PARTIR DE NBT0
      GOTO 9000
C
C     LE 3-EME SOMMET DEVIENT LE 1-ER : 34125 => 12345
 913  L        = NUSOM(1)
      NUSOM(1) = NUSOM(3)
      NUSOM(3) = L
      L        = NUSOM(2)
      NUSOM(2) = NUSOM(4)
      NUSOM(4) = L
      IF( NBT0 .GT. 0 )  GOTO 9900
C     TANGENTES A PERMUTER CIRCULAIREMENT A PARTIR DE NBT0
      GOTO 9000
C
C     LE 4-EME SOMMET DEVIENT LE 1-ER : 41235 => 12345
 914  L        = NUSOM(1)
      NUSOM(1) = NUSOM(4)
      NUSOM(4) = NUSOM(3)
      NUSOM(3) = NUSOM(2)
      NUSOM(2) = L
      IF( NBT0 .GT. 0 ) GOTO 9900
      GOTO 9000
C
C     LE 5-EME SOMMET NE PEUT PAS DEVENIR LE PREMIER
C     CAR LA PREMIERE FACE EST 1234 et les AUTRES SONT TRIANGULAIRES
C     ET CETTE STRUCTURE NE PEUT ETRE MODIFIEE
C     LA CORRECTION EST FAITE DANS HACHAG
 915  RETURN
C
C     SORTIE NORMALE
 9000 RETURN
C
 9900 IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = 'HACREN: PYRAMIDE SANS TANGENTES PROGRAMMEES'
      ELSE
         KERR(1) = 'HACREN: PYRAMID WITHOUT PROGRAMMED TANGENTS'
      ENDIF
      CALL LEREUR
      RETURN
      END
