      SUBROUTINE ENT1AD( BASE , A , B , R )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     R EST L'ENTIER MULTI-MOTS   SOMME DE 2 ENTIERS A>0 ET B>0
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
      INTEGER A, B, R(-2:*)
      INTEGER BASE,SOM
C
C     0 + 0 + 0 =< SOM =< (BASE-1) + (BASE-1) + 1
C                         (BASE-1) + BASE
C     LA RETENUE EST AU PLUS EGALE A 1
C
      SOM = A + B
      IF( SOM .GE. BASE ) THEN
C        LE NOMBRE DE MOTS - 1 DE L'ENTIER
         R(-2)  = 1
C        L'ENTIER EST TOUJOURS POSITIF
         R(-1)  = 1
         R( 0)  = SOM - BASE
         R( 1)  = 1
      ELSE
C        LE NOMBRE DE MOTS - 1 DE L'ENTIER
         R(-2)  = 0
C        L'ENTIER EST TOUJOURS POSITIF
         R(-1)  = 1
         R( 0)  = SOM
      ENDIF
      END
