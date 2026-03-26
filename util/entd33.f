      SUBROUTINE ENTD33( BASE ,
     %                   A11 , A12 , A13 ,
     %                   A21 , A22 , A23 ,
     %                   A31 , A32 , A33 ,
     %                   DET33 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     DET33 EST LE DETERMINANT EXACT DE LA MATRICE A 3*3
C -----
C
C ENTREES :
C ---------
C BASE    : LA BASE DE TRAITEMENT
C AIJ     : ENTIER MULTI-MOTS EN POSITION I J DE LA MATRICE A
C
C SORTIES :
C ---------
C DET33   : ENTIER MULTI-MOTS  CODE SOUS LA FORME
C                   I=NDET33
C           DET33 = SOMME ( DET33(I) * BASE ** I )   ET
C                   I=0
C
C           DET33(-1) = SIGNE( DET33 )
C           DET33(-2) = NDET33 LE NOMBRE DE MOTS ENTIERS - 1  DE DET33
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS   FEVRIER 1988
C23456...............................................................012
      INTEGER BASE
      INTEGER A11(-2:*),A12(-2:*),A13(-2:*),
     %        A21(-2:*),A22(-2:*),A23(-2:*),
     %        A31(-2:*),A32(-2:*),A33(-2:*),
     %        DET33(-2:*)
      INTEGER DET22(-2:10),PROD(-2:10)
C
C                      (A22     A23)            (A21     A23
C     DET33 = A11 * DET(           ) - A12 * DET(           )  +
C                      (A32     A33)            (A31     A33)
C
C                      (A21     A22)
C           + A13 * DET(           )
C                      (A31     A32)
C
      CALL ENTD22( BASE ,
     %             A22 , A23 ,
     %             A32 , A33 , DET22 )
C
      CALL ENTMUL( BASE , A11 , DET22  , DET33 )
C
C     LE SIGNE - EST PRIS EN COMPTE PAR PERMUTATION DES 2 COLONNES
      CALL ENTD22( BASE ,
     %             A23 , A21 ,
     %             A33 , A31 , DET22 )
C
      CALL ENTMUL( BASE , A12    , DET22 , PROD  )
      CALL ENTADD( BASE , DET33  , PROD  , DET33 )
C
      CALL ENTD22( BASE ,
     %             A21 , A22 ,
     %             A31 , A32 , DET22 )
C
      CALL ENTMUL( BASE , A13 , DET22 , PROD   )
C
C     DET33 = A11 * D1 + A12 * D2 + A13 * D3
      CALL ENTADD( BASE , DET33 , PROD , DET33 )
      END
