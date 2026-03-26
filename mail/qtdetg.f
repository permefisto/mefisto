      SUBROUTINE QTDETG( MNSTCA, MN2DER, MNSOFA, MNFASU )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRANSFORMER DES DERIVEES ds/du ET ds/dv EN CHAQUE SOMMET
C -----    DES EF (TRIANGLE OU QUADRANGLE) DU CARRE OU TRIANGLE UNITE
C          EN LES 6 TANGENTES SUR CHAQUE TRIANGLE
C          EN LES 8 TANGENTES SUR CHAQUE QUADRANGLE
C          DU MAILLAGE DU QUADRANGLE OU TRIANGLE COURBE
C
C ENTREES:
C --------
C MNSTCA : ADRESSE MCN DU TABLEAU DES 2 COORDONNEES DES SOMMETS DES EF
C         (TRIANGLES OU QUADRANGLES) MAILLAGE DU CARRE OU TRIANGLE UNITE
C MN2DER : ADRESSE MCN DU TABLEAU DES DERIVEES ds/du ET ds/dv EN CHACUN
C          DES SOMMETS DU MAILLAGE DU CARRE OU TRIANGLE UNITE
C MNFASU : ADRESSE MCN DU TMS 'NSEF' DES NUMEROS DES SOMMETS DES EF
C          CF $MEFISTO/td/d/a___nsef
C MNSOFA : ADRESSE MCN DU TMS 'XYZSOMMET' DE LA SURFACE
C          CF $MEFISTO/td/d/a___xyzsommet
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS           AVRIL 1998
C234567..............................................................012
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/pp.inc"
      COMMON             MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE       (MCN(1),RMCN(1))
C
      INTEGER            NOSOEL(1:64),MN(0:5)
      DOUBLE PRECISION   U,V,DU1,DU2,DV1,DV2,DRU,DRV
C
C     LES PARAMETRES DES NO SOMMET DU MAILLAGE
      CALL NSEFPA( MCN(MNFASU) ,
     %             NUTYMA, NBSOEL, NBSOEF, NBTGEF,
     %             LDAPEF, LDNGEF, LDTGEF, NBEFOB,
     %             NX    , NY    , NZ    ,
     %             IERR   )
      IF( IERR .NE. 0 ) RETURN
C
      MNTG = MNSOFA + WYZSOM + 3 * MCN(MNSOFA+WNBSOM) - 1
      MNC  = MNSTCA - 2
C
      DO 100 N=1,NBEFOB
C
C        LE NUMERO DES NBSOEF SOMMETS DE L'EF N
         CALL NSEFNS( N     , NUTYMA, NBSOEF, NBTGEF,
     %                LDAPEF, LDNGEF, LDTGEF,
     %                MNFASU, NX, NY, NZ,
     %                NCOGEL, NUGEEF, NUEFTG, NOSOEL, IERR )
C        LE NOMBRE DE SOMMETS DE CET ELEMENT = NCOGEL 3 OU 4
C
C        ADRESSE DE LA PREMIERE COORDONNEE DES NCOGEL SOMMETS DE L'EF
         MN(1) = MNC + 2 * NOSOEL(1)
         MN(2) = MNC + 2 * NOSOEL(2)
         MN(3) = MNC + 2 * NOSOEL(3)
         MN(4) = MNC + 2 * NOSOEL(4)
C
C        POUR ASSURER LA PERMUTATION CIRCULAIRE SANS TEST
         MN(NCOGEL+1) = MN(1)
         MN(0)        = MN(NCOGEL)
C
         DO 30 I=1,NCOGEL
C
C           AU SOMMET I DU TRIANGLE N DU PLAN
            U = RMCN(MN(I)  )
            V = RMCN(MN(I)+1)
C
C           LA DIRECTION Si Si+1
            DU1 = RMCN(MN(I+1)  ) - U
            DV1 = RMCN(MN(I+1)+1) - V
C
C           LA DIRECTION Si Si-1
            DU2 = RMCN(MN(I-1)  ) - U
            DV2 = RMCN(MN(I-1)+1) - V
C           L'ADRESSE DES 2 DERIVEES /du ET /dv DU SOMMET
            MD1 = MN2DER + 6 * NOSOEL(I) - 7
            MD2 = MD1 + 3
C
C           LE CALCUL DES 3 COMPOSANTES DES 2 TANGENTES AUX ARETES
            DO 20 K=1,3
C
C              LA DERIVEE SELON LE 1-ER  PARAMETRE
               DRU = RMCN(MD1+K)
C              LA DERIVEE SELON LE 2-EME PARAMETRE
               DRV = RMCN(MD2+K)
C
C              LA DERIVEE SELON LA DIRECTION DE L'ARETE
               RMCN(MNTG  +K) = REAL( DRU * DU1 + DRV * DV1 )
               RMCN(MNTG+3+K) = REAL( DRU * DU2 + DRV * DV2 )
 20         CONTINUE
C
C           PASSAGE A LA TANGENTE SUIVANTE
            MNTG = MNTG + 6
C
 30      CONTINUE
C
 100  CONTINUE
      END
