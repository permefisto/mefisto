      SUBROUTINE ENTDET( BASE ,
     %                   P , P1 , P2 , P3 , P4 ,
     %                   DET55 , DET44 , NSIGNE )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     CALCUL DU DETERMINANT DE LA MATRICE 5 * 5
C -----
C           ( L**2  L1**2  L2**2  L3**2  L4**2 )
C           (   1      1      1      1      1  )
C           (   X     X1     X2     X3     X4  )
C           (   Y     Y1     Y2     Y3     Y4  )
C           (   Z     Z1     Z2     Z3     Z4  )
C
C           AVEC L **2 = X  * X  + Y  * Y  + Z  * Z
C                LI**2 = XI * XI + YI * YI + ZI * ZI
C                PI = XI YI ZI COORDONNEES DU POINT PI I=1,2,3,4
C
C           ET DU SIGNE DE DET(55) * DET(44 INFERIEURE DROITE)
C           SI CE SIGNE EST - , LE POINT P EST INTERIEUR
C           A LA BOULE CIRCONSCRITE AU TETRAEDRE P1 P2 P3 P4
C
C ENTREES :
C ---------
C BASE    : LA BASE DE TRAITEMENT
C P       : LES 3 COORDONNEES X Y Z DU POINT P
C P1 P2 P3 P4 : LES 4 POINTS DEFINIS PAR LEURS 3 COORDONNEES
C
C SORTIES :
C ---------
C DET     : ENTIER MULTI-MOTS  CODE SOUS LA FORME
C                I=NDET
C           DET = SOMME ( DET(I) * BASE ** I )   ET
C                I=0
C
C           DET(-1) = SIGNE( DET )
C           DET(-2) = NDET   LE NOMBRE DE MOTS ENTIERS - 1  DE DET
C
C NSIGNE  :  0 SI LE TETRAEDRE EST DEGENERE
C           -1 SI LE POINT P EST DANS OU SUR LA BOULE CIRCONSCRITE
C                  AU TETRAEDRE DE SOMMETS P1 P2 P3 P4
C            1 SI LE POINT P EST EXTERIEUR A LA BOULE CIRCONSCRITE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS   FEVRIER 1988
C23456...............................................................012
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      INTEGER           P(3),P1(3),P2(3),P3(3),P4(3)
      INTEGER           BASE,DET55(-2:10),DET44(-2:10),NSIGNE
C
ccc      DOUBLE PRECISION  VOLTET
      INTEGER      NP(1:3,0:4),MP(-2:0,1:3,0:4),ML(-2:2,0:4)
      INTEGER      UN(-2:0)
      INTEGER      X (-2:0),Y (-2:0),Z (-2:0),L (-2:2),
     %             X1(-2:0),Y1(-2:0),Z1(-2:0),L1(-2:2),
     %             X2(-2:0),Y2(-2:0),Z2(-2:0),L2(-2:2),
     %             X3(-2:0),Y3(-2:0),Z3(-2:0),L3(-2:2),
     %             X4(-2:0),Y4(-2:0),Z4(-2:0),L4(-2:2)
      EQUIVALENCE (MP(-2,1,0),X (-2)) , (MP(-2,1,1),X1(-2)) ,
     %            (MP(-2,1,2),X2(-2)) , (MP(-2,1,3),X3(-2)) ,
     %            (MP(-2,1,4),X4(-2))
      EQUIVALENCE (MP(-2,2,0),Y (-2)) , (MP(-2,2,1),Y1(-2)) ,
     %            (MP(-2,2,2),Y2(-2)) , (MP(-2,2,3),Y3(-2)) ,
     %            (MP(-2,2,4),Y4(-2))
      EQUIVALENCE (MP(-2,3,0),Z (-2)) , (MP(-2,3,1),Z1(-2)) ,
     %            (MP(-2,3,2),Z2(-2)) , (MP(-2,3,3),Z3(-2)) ,
     %            (MP(-2,3,4),Z4(-2))
C
      EQUIVALENCE (ML(-2,0),L (-2)) , (ML(-2,1),L1(-2)) ,
     %            (ML(-2,2),L2(-2)) , (ML(-2,3),L3(-2)) ,
     %            (ML(-2,4),L4(-2))
C
C     L'ENTIER 1 ECRIT EN ENTIER MULTI-MOTS
      DATA     UN/0,1,1/
C
C
C
C     TRANSFORMATION DES 5 POINTS EN COORDONNEES ENTIERES
C     REEL => ENTIER SANS MODIFICATION
      CALL ENTIDE( P  , NP(1,0) )
      CALL ENTIDE( P1 , NP(1,1) )
      CALL ENTIDE( P2 , NP(1,2) )
      CALL ENTIDE( P3 , NP(1,3) )
      CALL ENTIDE( P4 , NP(1,4) )
C
C     TRANSFORMATION DE CHAQUE COORDONNEE ENTIERE EN
C     COORDONNEE ENTIERE MULTI-MOTS
      DO 30 J=0,4
         DO 20 I=1,3
C           LE NOMBRE DE MOTS - 1 DE L'ENTIER
            MP(-2,I,J) = 0
C           L'ENTIER EST TOUJOURS POSITIF
            MP(-1,I,J) = 1
C           L'ENTIER SUR UN MOT CAR < BASE
            MP( 0,I,J) = NP(I,J)
 20      CONTINUE
 30   CONTINUE
C
C     LES CARRES DES DISTANCES DES POINTS A L'ORIGINE
      CALL ENTDIS( BASE , NP(1,0) , ML(-2,0) )
      CALL ENTDIS( BASE , NP(1,1) , ML(-2,1) )
      CALL ENTDIS( BASE , NP(1,2) , ML(-2,2) )
      CALL ENTDIS( BASE , NP(1,3) , ML(-2,3) )
      CALL ENTDIS( BASE , NP(1,4) , ML(-2,4) )
C
C     DEVELOPPEMENT DU DETERMINANT SUIVANT LA LIGNE DE UN
C     DET55 = - DET1 + DET2 - DET3 + DET4 - DET5
      CALL ENTD44( BASE ,
     %             X1 , X2 , X3 , X4 ,
     %             L1 , L2 , L3 , L4 ,
     %             Y1 , Y2 , Y3 , Y4 ,
     %             Z1 , Z2 , Z3 , Z4 ,
     %             DET55 )
C
      CALL ENTD44( BASE ,
     %             L  , L2 , L3 , L4 ,
     %             X  , X2 , X3 , X4 ,
     %             Y  , Y2 , Y3 , Y4 ,
     %             Z  , Z2 , Z3 , Z4 ,
     %             DET44 )
      CALL ENTADD( BASE , DET55 , DET44 , DET55 )
C
      CALL ENTD44( BASE ,
     %             X  , X1 , X3 , X4 ,
     %             L  , L1 , L3 , L4 ,
     %             Y  , Y1 , Y3 , Y4 ,
     %             Z  , Z1 , Z3 , Z4 ,
     %             DET44 )
      CALL ENTADD( BASE , DET55 , DET44 , DET55 )
C
      CALL ENTD44( BASE ,
     %             L  , L1 , L2 , L4 ,
     %             X  , X1 , X2 , X4 ,
     %             Y  , Y1 , Y2 , Y4 ,
     %             Z  , Z1 , Z2 , Z4 ,
     %             DET44 )
      CALL ENTADD( BASE , DET55 , DET44 , DET55 )
C
      CALL ENTD44( BASE ,
     %             X  , X1 , X2 , X3 ,
     %             L  , L1 , L2 , L3 ,
     %             Y  , Y1 , Y2 , Y3 ,
     %             Z  , Z1 , Z2 , Z3 ,
     %             DET44 )
      CALL ENTADD( BASE , DET55 , DET44 , DET55 )
CCC      WRITE(IMPRIM,10001) DET55(-2),DET55(-1),
CCC     %                   (DET55(I),I=DET55(-2),0,-1)
CCC10001 FORMAT(' DET55:',I4,' MOTS SIGNE:',I2/(5I14))
C
C     LE DETERMINANT DE LA MATRICE 4 4 INFERIEURE DROITE
      CALL ENTD44( BASE ,
     %             UN , UN , UN , UN ,
     %             X1 , X2 , X3 , X4 ,
     %             Y1 , Y2 , Y3 , Y4 ,
     %             Z1 , Z2 , Z3 , Z4 ,
     %             DET44 )
CCC      WRITE(IMPRIM,10002) DET44(-2),DET44(-1),
CCC     %                   (DET44(I),I=DET44(-2),0,-1)
CCC10002 FORMAT(' DET44:',I4,' MOTS SIGNE:',I2/(5I14))
C
C     LE TETRAEDRE EST IL DEGENERE ?
      IF( DET44(-2) .EQ. 0 .AND. DET44(0) .EQ. 0 ) THEN
C        OUI
         NSIGNE = 0
ccc         WRITE(IMPRIM,10010) VOLTET( P1 , P2 , P3 , P4 )
ccc10010 FORMAT(' ENTDET:VOLUME DU TETRAEDRE=',G15.6)
      ELSE
C        LE POINT P EST INTERIEUR A LA BOULE
C        SI SIGNE(DET55)*SIGNE(DET44)<0
         IF( DET55(-2) .EQ. 0 .AND. DET55(0) .EQ. 0 ) THEN
C           LE POINT P EST SUR LA BOULE
            NSIGNE = -1
         ELSE
            NSIGNE = DET55(-1)  * DET44(-1)
         ENDIF
      ENDIF
CCC      WRITE(IMPRIM,10100) P,NSIGNE,1,P1,2,P2,3,P3,4,P4
CCC10100 FORMAT(' P :',3I15,'  RESULTAT ENTDET=',I2/
CCC     %4(' P',I1,':',3I15/))
      END
