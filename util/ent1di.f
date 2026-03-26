      SUBROUTINE ENT1DI( BASE , A , B , R )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     R EST L'ENTIER MULTI-MOTS DIFFERENCE DE 2 ENTIERS A>0 ET B>0
C -----
C
C ENTREES :
C ---------
C BASE    : LA BASE DE TRAITEMENT
C           ATTENTION ELLE DOIT ETRE CALCULEE DE TELLE SORTE QUE
C                     BASE + (BASE-1)  NE DEBORDE PAS UN ENTIER MACHINE
C
C A       : ENTIER
C B       : ENTIER
C
C SORTIES :
C ---------
C R       : ENTIER MULTI-MOTS  CODE SOUS LA FORME
C                I=NR
C           R = SOMME ( R(I) * BASE ** I )   ET  R(-1) = SIGNE( R )
C                I=0                             R(-2) = NR
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS         MAI 1991
C23456...............................................................012
      INTEGER A, B, R(-2:0)
      INTEGER BASE, DIFF
C
      print *,'ENT1DI: BASE=',BASE
C
      DIFF = A - B
C
C     STOCKAGE EN MULTI-MOTS DE 1 MOT ICI
C     LE NOMBRE DE MOTS - 1 DE L'ENTIER
      R(-2) = 0
      IF( DIFF .LT. 0 ) THEN
C        L'ENTIER EST NEGATIF
         R(-1)  = -1
         R(0)   = -DIFF
      ELSE
C        L'ENTIER EST POSITIF OU NUL
         R(-1)  = 1
         R( 0)  = DIFF
      ENDIF
C
      RETURN
      END
