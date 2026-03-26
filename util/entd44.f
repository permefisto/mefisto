      SUBROUTINE ENTD44( BASE ,
     %                   A11 , A12 , A13 , A14 ,
     %                   A21 , A22 , A23 , A24 ,
     %                   A31 , A32 , A33 , A34 ,
     %                   A41 , A42 , A43 , A44 ,
     %                   DET44 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     DET44 EST LE DETERMINANT EXACT DE LA MATRICE A 4 * 4
C -----
C
C ENTREES :
C ---------
C BASE    : LA BASE DE TRAITEMENT
C AIJ     : ENTIER MULTI-MOTS EN POSITION I J DE LA MATRICE A
C
C SORTIES :
C ---------
C DET44   : ENTIER MULTI-MOTS  CODE SOUS LA FORME
C                   I=NDET44
C           DET44 = SOMME ( DET44(I) * BASE ** I )   ET
C                   I=0
C
C           DET44(-1) = SIGNE( DET44 )
C           DET44(-2) = NDET44 LE NOMBRE DE MOTS ENTIERS - 1  DE DET44
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS   FEVRIER 1988
C23456...............................................................012
      INTEGER BASE
      INTEGER A11(-2:*),A12(-2:*),A13(-2:*),A14(-2:*),
     %        A21(-2:*),A22(-2:*),A23(-2:*),A24(-2:*),
     %        A31(-2:*),A32(-2:*),A33(-2:*),A34(-2:*),
     %        A41(-2:*),A42(-2:*),A43(-2:*),A44(-2:*),
     %        DET44(-2:*)
C     LES TABLEAUX AUXILIAIRES
      INTEGER DET33(-2:10),PROD(-2:10)
C
C                      (A22 A23 A24)            (A21 A23 A24)
C     DET44 = A11 * DET(A32 A33 A34) - A12 * DET(A31 A33 A34)  +
C                      (A42 A43 A44)            (A41 A43 A44)
C
C                      (A21 A22 A24)            (A21 A22 A23)
C           + A13 * DET(A31 A32 A34) - A14 * DET(A31 A32 A33)
C                      (A41 A42 A44)            (A41 A42 A43)
C
      CALL ENTD33( BASE ,
     %             A22 , A23 , A24 ,
     %             A32 , A33 , A34 ,
     %             A42 , A43 , A44 ,    DET33 )
      CALL ENTMUL( BASE , A11 , DET33 , DET44 )
C
C     LE SIGNE - EST PRIS EN COMPTE PAR LA PERMUTATION DE 2 COLONNES
      CALL ENTD33( BASE ,
     %             A23 , A21 , A24 ,
     %             A33 , A31 , A34 ,
     %             A43 , A41 , A44 ,       DET33 )
      CALL ENTMUL( BASE , A12    , DET33 , PROD  )
      CALL ENTADD( BASE , DET44  , PROD  , DET44 )
C
      CALL ENTD33( BASE ,
     %             A21 , A22 , A24 ,
     %             A31 , A32 , A34 ,
     %             A41 , A42 , A44 ,       DET33 )
      CALL ENTMUL( BASE , A13    , DET33 , PROD  )
      CALL ENTADD( BASE , DET44  , PROD  , DET44 )
C
C     LE SIGNE - EST PRIS EN COMPTE PAR LA PERMUTATION DE 2 COLONNES
      CALL ENTD33( BASE ,
     %             A22 , A21 , A23 ,
     %             A32 , A31 , A33 ,
     %             A42 , A41 , A43 ,    DET33 )
      CALL ENTMUL( BASE , A14 , DET33 , PROD  )
C
C     DET44 = A11 * D1 + A12 * D2 + A13 * D3 + A14 * D4
      CALL ENTADD( BASE , DET44 , PROD , DET44 )
      END
