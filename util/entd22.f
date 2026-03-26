      SUBROUTINE ENTD22( BASE ,
     %                   A11 , A12 ,
     %                   A21 , A22 ,
     %                   DET22 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :      DET22 EST LE DETERMINANT EXACT DE LA MATRICE A 2 * 2
C -----
C
C ENTREES :
C ---------
C BASE    : LA BASE DE TRAITEMENT
C AIJ     : ENTIER MULTI-MOTS EN POSITION I J DE LA MATRICE A
C
C SORTIES :
C ---------
C DET22   : ENTIER MULTI-MOTS  CODE SOUS LA FORME
C                   I=NDET
C           DET22 = SOMME ( DET22(I) * BASE ** I )   ET
C                   I=0
C
C           DET22(-1) = SIGNE( DET )
C           DET22(-2) = NDET   LE NOMBRE DE MOTS ENTIERS - 1  DE DET22
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS   FEVRIER 1988
C23456...............................................................012
      INTEGER BASE
      INTEGER A11(-2:*),A12(-2:*),
     %        A21(-2:*),A22(-2:*),
     %        DET22(-2:*),AUX(-2:10)
C
C     DET22 = A11 * A22
      CALL ENTMUL( BASE , A11 , A22 , DET22 )
C
C     AUX = A21 * A12
      CALL ENTMUL( BASE , A21 , A12 , AUX   )
C
C     AUX = - A21 * A12
      AUX(-1) = - AUX(-1)
C
C     DET22 = A11 * A22 - A21 * A12
      CALL ENTADD( BASE , DET22  , AUX , DET22 )
      END
