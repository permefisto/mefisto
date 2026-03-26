      SUBROUTINE COMORS( NTDL , LPLIGN , LPCOLO , A )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT:     SUPPRIMER LES ZEROS DANS LA MATRICE A DE STOCKAGE MORSE
C ----
C
C ENTREE :
C --------
C NTDL   : ORDRE DE LA MATRICE AG
C
C MODIFIES:
C ---------
C LPLIGN : LE POINTEUR SUR LE DERNIER COEFFICIENT DE CHAQUE LIGNE DE A
C LPCOLO : LE NUMERO DE COLONNE DES COEFFICIENTS  DE LA MATRICE MORSE A
C A      : LA MATRICE MORSE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY  ANALYSE NUMERIQUE UPMC PARIS            MAI 1990
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      COMMON / EPSSSS / EPZERO,EPSXYZ
      INTEGER           LPLIGN(NTDL+1),LPCOLO(1:*)
      DOUBLE PRECISION  A(1:*),PIV
C
       LONG0 = LPLIGN(NTDL+1)
       LP    = 0
       NPIV  = 0
       J1    = 1
       DO 100 I=1,NTDL
          J2 = LPLIGN(I+1)
          PIV = A(J2)
          IF (DABS(PIV).LT.EPZERO) THEN
              NPIV = NPIV + 1
              A(J2) = EPZERO
              GO TO 100
          END IF
          DO 101 J=J1,J2-1
             IF ( DABS( A(J) ) / PIV .GT. EPZERO ) THEN
                LP = LP + 1
                LPCOLO(LP) = LPCOLO(J)
                A(LP) = A(J)
             END IF
 101      CONTINUE
          LP = LP + 1
          LPCOLO(LP) = LPCOLO(J2)
          A(LP) = A(J2)
          LPLIGN(I+1)=LP
          J1 = J2 + 1
 100   CONTINUE
       LONG1 = LPLIGN(NTDL+1)
       IF( LANGAG .EQ. 0 ) THEN
          WRITE (IMPRIM,10000) NTDL,LONG0,LONG1
       ELSE
          WRITE (IMPRIM,20000) NTDL,LONG0,LONG1
       ENDIF
10000 FORMAT(' COMPRESSION DE LA MATRICE MORSE',/,
     %        ' NOMBRE D''INCONNUES',T40,I10/,
     %        ' NOMBRE INITIAL DE COEFFICIENTS ',T40,I10/,
     %        ' NOMBRE FINAL   DE COEFFICIENTS ',T40,I10)
20000 FORMAT(' SPARSE MATRIX COMPRESSION',/,
     %        ' DEGREES of FREEDOM NUMBER',T40,I10/,
     %        ' INITIAL COEFFICIENT NUMBER ',T40,I10/,
     %        ' FINAL   COEFFICIENT NUMBER ',T40,I10)
      IF( NPIV .GT. 0 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            WRITE (IMPRIM,10001) NPIV,EPZERO
         ELSE
            WRITE (IMPRIM,20001) NPIV,EPZERO
         ENDIF
      ENDIF
10001 FORMAT(' ATTENTION :',I8,' COEFFICIENTS DIAGONAUX <',G13.6)
20001 FORMAT(' ATTENTION :',I8,' DIAGONAL COEFFICIENTS <',G13.6)
C
      RETURN
      END
