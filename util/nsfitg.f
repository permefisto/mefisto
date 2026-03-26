      FUNCTION NSFITG( NS, NTG, NCOGEL )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  RETOURNER LE NUMERO DU SOMMET PROCHE DE L'EXTREMITE DE
C -----  LA TANGENTE NTG AU SOMMET NS DE L'EF DE CODE GEOMETRIQUE NCOGEL
C
C ENTREES:
C --------
C NS     : NUMERO DU SOMMET DE L'EF
C NUTG   : NUMERO DE LA TANGENTE AU SOMMET
C NCOGEL : CODE GEOMETRIQUE DE L ELEMENT
C          1:POINT 2:SEGMENT 3:TRIANGLE 4:QUADRANGLE
C          5:TETRAEDRE 6:PENTAEDRE 7:HEXAEDRE >7:ERREUR
C
C SORTIE :
C --------
C NSFITG : NUMERO DU SOMMET PROCHE DE L'EXTREMITE DE LA TANGENTE NT
C          AU SOMMET NS DE L'EF DE CODE GEOMETRIQUE NCOGEL
C
C    TRIANGLE  :
C                NO TANGENTE1(S1S2), NT2(S1S3),   NT3(S2S3), NT4(S2S1),
C                NO TANGENTE5(S3S1), NT6(S3S2)
C    QUADRANGLE:
C                NO TANGENTE1(S1S2), NT2(S1S3),   NT3(S2S3), NT4(S2S1),
C                NO TANGENTE5(S3S4), NT6(S3S2),   NT7(S4S1), NT8(S4S3)
C    TETRAEDRE :
C                NO TANGENTE1(S1S2), NT2(S1S3),  NT3(S1S4),
C                NO TANGENTE4(S2S3), NT5(S2S1),  NT6(S2S4),
C                NO TANGENTE7(S3S1), NT8(S3S2),  NT9(S3S4),
C                NO TANGENT10(S4S1), NT11(S4S2), NT12(S4S3)
C    PENTAEDRE :
C                NO TANGENTE1(S1S2), NT2(S1S3),  NT3(S1S4),
C                NO TANGENTE4(S2S3), NT5(S2S1),  NT6(S2S5),
C                NO TANGENTE7(S3S1), NT8(S3S2),  NT9(S3S6),
C                NO TANGENT10(S4S5), NT11(S4S6), NT12(S4S1),
C                NO TANGENT13(S5S6), NT14(S5S4), NT15(S5S2),
C                NO TANGENT16(S6S4), NT17(S6S5), NT18(S6S3)
C    HEXAEDRE  :
C                NO TANGENTE1(S1S2), NT2(S1S4),  NT3(S1S5),
C                NO TANGENTE4(S2S3), NT5(S2S1),  NT6(S2S6),
C                NO TANGENTE7(S3S4), NT8(S3S2),  NT9(S3S7),
C                NO TANGENT10(S4S1), NT11(S4S3), NT12(S4S8),
C                NO TANGENT13(S5S6), NT14(S5S8), NT15(S5S1),
C                NO TANGENT16(S6S7), NT17(S6S5), NT18(S6S2),
C                NO TANGENT19(S7S8), NT20(S7S6), NT21(S7S3),
C                NO TANGENT22(S8S5), NT23(S8S7), NT24(S8S4)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC PARIS    DECEMBRE 1996
C2345X7..............................................................012
      INTEGER  NSTTETR(3,4), NSTPENT(3,6), NSTHEXA(3,8)
      DATA     NSTTETR / 2,3,4, 3,1,4, 1,2,4, 1,2,3 /
      DATA     NSTPENT / 2,3,4, 3,1,5, 1,2,6, 5,6,1, 6,4,2, 4,5,3 /
      DATA     NSTHEXA / 2,4,5, 3,1,6, 4,2,7, 1,3,8,
     %                   6,8,1, 7,5,2, 8,6,3, 1,7,4  /
C
      GOTO( 10, 20, 30, 30, 50, 60, 70 ) , NCOGEL
C
C     POINT
C     =====
 10   NSFITG = 0
      RETURN
C
C     SEGMENT
C     =======
 20   IF( NS .EQ. 1 ) THEN
         NSFITG = 2
      ELSE
         NSFITG = 1
      ENDIF
      RETURN
C
C     TRIANGLE OU QUADRANGLE
C     ======================
 30   IF( NTG .EQ. 1 ) THEN
C
C        PREMIERE TANGENTE AU SOMMET
         IF( NS .NE. NCOGEL ) THEN
            NSFITG = NS + 1
         ELSE
            NSFITG = 1
         ENDIF
C
      ELSE
C
C        SECONDE TANGENTE AU SOMMET
         IF( NS .NE. 1 ) THEN
            NSFITG = NS - 1
         ELSE
            NSFITG = NCOGEL
         ENDIF
      ENDIF
      RETURN
C
C     TETRAEDRE
C     =========
 50   NSFITG = NSTTETR( NTG, NS )
      RETURN
C
C     PENTAEDRE
C     =========
 60   NSFITG = NSTPENT( NTG, NS )
      RETURN
C
C     HEXAEDRE
C     ========
 70   NSFITG = NSTHEXA( NTG, NS )
      RETURN
      END
