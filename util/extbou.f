      SUBROUTINE EXTBOU( BASE , P , P1 , P2 , P3 , P4 ,  NSIGNE )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     LE POINT P EST IL EXTERIEUR A LA BOULE CIRCONSCRITE
C -----     AU TETRAEDRE DE SOMMETS P1,P2,P3,P4 ?
C
C           OUI => NSIGNE=-1    ( P HORS DE LA BOULE )
C           NON => NSIGNE=+1    ( P SUR OU DANS LA BOULE)
C           P1P2P3P4 DEGENERE => NSIGNE= 0
C
C           POUR ETRE CERTAIN DU RESULTAT
C           EMPLOI D'ENTIERS MULTI-MOTS  CODES SOUS LA FORME
C                 I=NENT
C           ENT = SOMME ( ENT(I) * BASE ** I )
C                 I=0
C           ENT(-1) = SIGNE( ENT )
C           ENT(-2) = NENT   LE NOMBRE DE MOTS ENTIERS - 1  DE ENT
C
C ENTREES :
C ---------
C BASE    : LA BASE DE TRAITEMENT
C P       : LES 3 COORDONNEES X Y Z DU POINT P
C P1 P2 P3 P4 : LES 4 POINTS DEFINIS PAR LEURS 3 COORDONNEES
C
C SORTIES :
C ---------
C D OU DELTA : ENTIER MULTI-MOTS  CODE SOUS LA FORME
C                I=ND
C           D = SOMME ( D(I) * BASE ** I )   ET
C                I=0
C
C           D(-1) = SIGNE( D )
C           D(-2) = ND   LE NOMBRE DE MOTS ENTIERS - 1  DE D
C D       : RESULTAT INTERMEDIAIRE
C DELTA   : 6 FOIS LE VOLUME DU TETRAEDRE P1 P2 P3 P4
C
C NSIGNE  :  0 SI LE TETRAEDRE P1 P2 P3 P4 EST DEGENERE
C           -1 SI LE POINT P EST DANS OU SUR LA BOULE CIRCONSCRITE
C                  AU TETRAEDRE DE SOMMETS P1 P2 P3 P4
C            1 SI LE POINT P EST EXTERIEUR A LA BOULE CIRCONSCRITE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS         MAI 1991
C23456...............................................................012
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      INTEGER           P(3),P1(3),P2(3),P3(3),P4(3)
      INTEGER           BASE,DELTA(-2:20),NSIGNE
C
      INTEGER           IX12(-2:1),IY12(-2:1),IZ12(-2:1),
     %                  IX13(-2:1),IY13(-2:1),IZ13(-2:1),
     %                  IX14(-2:1),IY14(-2:1),IZ14(-2:1),
     %                  JX  (-2:1),JY  (-2:1),JZ  (-2:1),
     %                  A(-2:20),B(-2:20),C(-2:20),
     %                  D(-2:20),E(-2:20),F(-2:20),
     %                  G(-2:20),H(-2:20),I(-2:20),
     %                  AUX1(-2:20),AUX2(-2:20)
ccc      DOUBLE PRECISION  VOLTET
C
C     LES DIFFERENCES ENTRE COORDONNEES
C     =================================
C     L'ENTIER SUR UN MOT CAR < BASE
      CALL ENT1DI( BASE , P1(1) , P2(1) , IX12 )
      CALL ENT1DI( BASE , P1(1) , P3(1) , IX13 )
      CALL ENT1DI( BASE , P1(1) , P4(1) , IX14 )
C
      CALL ENT1DI( BASE , P1(2) , P2(2) , IY12 )
      CALL ENT1DI( BASE , P1(2) , P3(2) , IY13 )
      CALL ENT1DI( BASE , P1(2) , P4(2) , IY14 )
C
      CALL ENT1DI( BASE , P1(3) , P2(3) , IZ12 )
      CALL ENT1DI( BASE , P1(3) , P3(3) , IZ13 )
      CALL ENT1DI( BASE , P1(3) , P4(3) , IZ14 )
C
C     LE SECOND MEMBRE A B C DU SYSTEME POUR TROUVER LE CENTRE DU CERCLE
C     ==================================================================
C     CALCUL DE A
      CALL ENT1AD( BASE , P1(1) , P2(1) , JX )
      CALL ENT1AD( BASE , P1(2) , P2(2) , JY )
      CALL ENT1AD( BASE , P1(3) , P2(3) , JZ )
      CALL ENTMUL( BASE , IX12 , JX , A )
      CALL ENTMUL( BASE , IY12 , JY , AUX1 )
      CALL ENTADD( BASE , A , AUX1 , A )
      CALL ENTMUL( BASE , IZ12 , JZ , AUX1 )
      CALL ENTADD( BASE , A , AUX1 , A )
C
C     CALCUL DE B
      CALL ENT1AD( BASE , P1(1) , P3(1) , JX )
      CALL ENT1AD( BASE , P1(2) , P3(2) , JY )
      CALL ENT1AD( BASE , P1(3) , P3(3) , JZ )
      CALL ENTMUL( BASE , IX13 , JX , B )
      CALL ENTMUL( BASE , IY13 , JY , AUX1 )
      CALL ENTADD( BASE , B , AUX1 , B )
      CALL ENTMUL( BASE , IZ13 , JZ , AUX1 )
      CALL ENTADD( BASE , B , AUX1 , B )
C
C     CALCUL DE C
      CALL ENT1AD( BASE , P1(1) , P4(1) , JX )
      CALL ENT1AD( BASE , P1(2) , P4(2) , JY )
      CALL ENT1AD( BASE , P1(3) , P4(3) , JZ )
      CALL ENTMUL( BASE , IX14 , JX , C )
      CALL ENTMUL( BASE , IY14 , JY , AUX1 )
      CALL ENTADD( BASE , C , AUX1 , C )
      CALL ENTMUL( BASE , IZ14 , JZ , AUX1 )
      CALL ENTADD( BASE , C , AUX1 , C )
C
C     DELTA = 6 * VOLUME DU TETRAEDRE
C     ===============================
      CALL ENTMUL( BASE , IY13 , IZ14 , AUX1 )
      CALL ENTMUL( BASE , IY14 , IZ13 , AUX2 )
C
C     LA DIFFERENCE SE FAIT PAR SOMME APRES CHANGEMENT DE SIGNE
      AUX2(-1) = -AUX2(-1)
      CALL ENTADD( BASE , AUX1 , AUX2 , D )
C
      CALL ENTMUL( BASE , IY14 , IZ12 , AUX1 )
      CALL ENTMUL( BASE , IY12 , IZ14 , AUX2 )
      AUX2(-1) = -AUX2(-1)
      CALL ENTADD( BASE , AUX1 , AUX2 , E )
C
      CALL ENTMUL( BASE , IY12 , IZ13 , AUX1 )
      CALL ENTMUL( BASE , IY13 , IZ12 , AUX2 )
      AUX2(-1) = -AUX2(-1)
      CALL ENTADD( BASE , AUX1 , AUX2 , F )
C
      CALL ENTMUL( BASE , IX12 , D    , AUX1  )
      CALL ENTMUL( BASE , IX13 , E    , AUX2  )
      CALL ENTADD( BASE , AUX1 , AUX2 , DELTA )
      CALL ENTMUL( BASE , IX14 , F    , AUX2  )
      CALL ENTADD( BASE , DELTA, AUX2 , DELTA )
C
C     LE TETRAEDRE EST IL DEGENERE ?
C     ==============================
      IF( DELTA(-2) .EQ. 0 .AND. DELTA(0) .EQ. 0 ) THEN
C        OUI
         NSIGNE = 0
ccc         WRITE(IMPRIM,10010) VOLTET( P1 , P2 , P3 , P4 )
ccc10010 FORMAT(' EXTBOU:VOLUME DU TETRAEDRE=',G15.6)
         RETURN
      ENDIF
C
C     TETRAEDRE NON DEGENERE
C     ======================
C     CALCUL DE G
      CALL ENTMUL( BASE , A    , D    , AUX1  )
      CALL ENTMUL( BASE , B    , E    , AUX2  )
      CALL ENTADD( BASE , AUX1 , AUX2 , G     )
      CALL ENTMUL( BASE , C    , F    , AUX2  )
      CALL ENTADD( BASE , G    , AUX2 , G     )
C
C     CALCUL DE H
      CALL ENTMUL( BASE , IX14 , IZ13 , AUX1 )
      CALL ENTMUL( BASE , IX13 , IZ14 , AUX2 )
      AUX2(-1) = -AUX2(-1)
      CALL ENTADD( BASE , AUX1 , AUX2 , D )
C
      CALL ENTMUL( BASE , IX12 , IZ14 , AUX1 )
      CALL ENTMUL( BASE , IX14 , IZ12 , AUX2 )
      AUX2(-1) = -AUX2(-1)
      CALL ENTADD( BASE , AUX1 , AUX2 , E )
C
      CALL ENTMUL( BASE , IX13 , IZ12 , AUX1 )
      CALL ENTMUL( BASE , IX12 , IZ13 , AUX2 )
      AUX2(-1) = -AUX2(-1)
      CALL ENTADD( BASE , AUX1 , AUX2 , F )
C
      CALL ENTMUL( BASE , A    , D    , AUX1  )
      CALL ENTMUL( BASE , B    , E    , AUX2  )
      CALL ENTADD( BASE , AUX1 , AUX2 , H     )
      CALL ENTMUL( BASE , C    , F    , AUX2  )
      CALL ENTADD( BASE , H    , AUX2 , H     )
C
C     CALCUL DE I
      CALL ENTMUL( BASE , IX13 , IY14 , AUX1 )
      CALL ENTMUL( BASE , IX14 , IY13 , AUX2 )
      AUX2(-1) = -AUX2(-1)
      CALL ENTADD( BASE , AUX1 , AUX2 , D )
C
      CALL ENTMUL( BASE , IX14 , IY12 , AUX1 )
      CALL ENTMUL( BASE , IX12 , IY14 , AUX2 )
      AUX2(-1) = -AUX2(-1)
      CALL ENTADD( BASE , AUX1 , AUX2 , E )
C
      CALL ENTMUL( BASE , IX12 , IY13 , AUX1 )
      CALL ENTMUL( BASE , IX13 , IY12 , AUX2 )
      AUX2(-1) = -AUX2(-1)
      CALL ENTADD( BASE , AUX1 , AUX2 , F )
C
      CALL ENTMUL( BASE , A    , D    , AUX1  )
      CALL ENTMUL( BASE , B    , E    , AUX2  )
      CALL ENTADD( BASE , AUX1 , AUX2 , I     )
      CALL ENTMUL( BASE , C    , F    , AUX2  )
      CALL ENTADD( BASE , I    , AUX2 , I     )
C
C     LES PRODUITS FINAUX
C     ===================
      CALL ENT1AD( BASE , P1(1) , P(1) , IX12 )
      CALL ENTMUL( BASE , DELTA , IX12 , AUX1 )
      G(-1) = -G(-1)
      CALL ENTADD( BASE , AUX1  , G    , A    )
C
      CALL ENT1AD( BASE , P1(2) , P(2) , IY12 )
      CALL ENTMUL( BASE , DELTA , IY12 , AUX1 )
      H(-1) = -H(-1)
      CALL ENTADD( BASE , AUX1  , H    , B    )
C
      CALL ENT1AD( BASE , P1(3) , P(3) , IZ12 )
      CALL ENTMUL( BASE , DELTA , IZ12 , AUX1 )
      I(-1) = -I(-1)
      CALL ENTADD( BASE , AUX1  , I    , C    )
C
      CALL ENT1DI( BASE , P1(1) , P(1) , IX12 )
      CALL ENTMUL( BASE , IX12  , A   , AUX1 )
C
      CALL ENT1DI( BASE , P1(2) , P(2) , IY12 )
      CALL ENTMUL( BASE , IY12 , B    , AUX2 )
      CALL ENTADD( BASE , AUX1 , AUX2 , D    )
C
      CALL ENT1DI( BASE , P1(3) , P(3) , IZ12 )
      CALL ENTMUL( BASE , IZ12 , C    , AUX2 )
      CALL ENTADD( BASE , D    , AUX2 , D    )
C
C     LE SIGNE FINAL
C     ==============
C     LE POINT P EST EXTERIEUR A LA BOULE  SI SIGNE(D)*SIGNE(DELTA)<0
      IF( D(-2) .EQ. 0 .AND. D(0) .EQ. 0 ) THEN
C        LE POINT P EST SUR LA BOULE
         NSIGNE = -1
      ELSE
C        LE PRODUIT DES SIGNES DECIDE
         NSIGNE = - D(-1)  * DELTA(-1)
      ENDIF
      END
