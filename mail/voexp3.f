         SUBROUTINE  VOEXP3( M, N, XYSTR1, XYSTR5, XYSCA1,
     %                       COSO )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    MAILLER UN PENTAEDRE DEFINI PAR SES FACES
C -----    PAR INTERPOLATION TRANSFINIE DEGRE 1
C          (cf CRAS A.PERRONNET)
C
C ENTREES:
C --------
C M      : NOMBRE DE SOMMETS PAR ARETE DES FACES TRIANGULAIRES
C N      : NOMBRE DE SOMMETS PAR ARETE EN Z DES FACES QUADRILATERES
C XYSTR1 : LES 2 COORDONNEES DES POINTS SUR LE TRIANGLE UNITE
C          DE LA FACE 1 DU PENTAEDRE UNITE
C XYSTR5 : LES 2 COORDONNEES DES POINTS SUR LE TRIANGLE UNITE
C          DE LA FACE 5 DU PENTAEDRE UNITE
C XYSCA1 : LES 2 COORDONNEES DES POINTS SUR LE CARRE UNITE
C          DES FACES 2 3 4 DU PENTAEDRE UNITE
C COSO   : LES COORDONNEES DES SOMMETS SUR LES FACES DU PENTAEDRE COURBE
C
C ENTREES ET SORTIES :
C --------------------
C COSO   : LES COORDONNEES DES SOMMETS DU PENTAEDRE COURBE APRES MAILLAGE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 1997
C23456---------------------------------------------------------------012
      REAL         XYSTR1(2,*), XYSTR5(2,*)
      REAL         XYSCA1(2,M,N,2:4)
      REAL         COSO(3,*)
      REAL         XYZ(3,0:8), XYZF(3,8), XYZA(3,15)
      EQUIVALENCE (CBA2,XYZ(1,0)), (CBA3,XYZ(2,0)), (Z,XYZ(3,0))
C
C     LES FONCTIONS FORMULES DES NUMEROS DES SOMMETS DU TRIANGLE ET PENTAEDRE
      NUSOTR(I,J)   = ( I * I - I ) / 2 + J
      NUSOPE(I,J,K) = J + ( I * I - I + (K-1) * ( M * M + M ) ) / 2
C
C     LE NUMERO DES SOMMETS DU PENTAEDRE
      NBST = NUSOTR(M,M)
      N2   = NUSOPE(M,1,1)
      N3   = NUSOPE(M,M,1)
      N4   = NUSOPE(1,1,N)
      N5   = NUSOPE(M,1,N)
      N6   = NUSOPE(M,M,N)
C
C     LA BOUCLE SUR LES SOMMETS INTERNES
      DO 90 K = 2, N-1
         DO 80 I = 3, M-1
            DO 70 J = 2, I-1
C
C              LA DROITE D1: SF2(I,K) SF4(M+1-I,K)
C              LE SOMMET SF2(I,K)
               XYZ(1,1) = XYSCA1(1,I,K,2)
               XYZ(2,1) = 0
               XYZ(3,1) = XYSCA1(2,I,K,2)
C              LE SOMMET SF4(M+1-I,K)
               XYZ(1,5) = 0
               XYZ(2,5) = 1 - XYSCA1(1,M+1-I,K,4)
               XYZ(3,5) = XYSCA1(2,M+1-I,K,4)
C
C              LA DROITE D2: SF4(M+1-J,K) SF3(J,K)
C              LE SOMMET SF4(M+1-J,K)
               XYZ(1,2) = 0
               XYZ(2,2) = 1 - XYSCA1(1,M+1-J,K,4)
               XYZ(3,2) = XYSCA1(2,M+1-J,K,4)
C              LE SOMMET SF3(J,K)
               XYZ(1,6) = 1 - XYSCA1(1,J,K,3)
               XYZ(2,6) = XYSCA1(1,J,K,3)
               XYZ(3,6) = XYSCA1(2,J,K,3)
C
C              LA DROITE D3: SF3(M-I+J,K) SF2(I-J+1,K)
C              LE SOMMET SF3(M-I+J,K)
               XYZ(1,3) = 1 - XYSCA1(1,M-I+J,K,3)
               XYZ(2,3) = XYSCA1(1,M-I+J,K,3)
               XYZ(3,3) = XYSCA1(2,M-I+J,K,3)
C              LE SOMMET SF2(I-J+1,K)
               XYZ(1,7) = XYSCA1(1,I-J+1,K,2)
               XYZ(2,7) = 0
               XYZ(3,7) = XYSCA1(2,I-J+1,K,2)
C
C              LA DROITE D4: SF1(I,J) SF5(I,J)
C              LE SOMMET SF1(I,J)
               XYZ(1,4) = XYSTR1(1,NUSOTR(I,J))
               XYZ(2,4) = XYSTR1(2,NUSOTR(I,J))
               XYZ(3,4) = 0
C              LE SOMMET SF5(I,J)
               XYZ(1,8) = XYSTR5(1,NUSOTR(I,J))
               XYZ(2,8) = XYSTR5(2,NUSOTR(I,J))
               XYZ(3,8) = 1
C
C              RECHERCHE DU POINT A DISTANCE MINIMALE DE CES 4 DROITES
C              DANS LE PENTAEDRE UNITE
               CALL MINDIS( 4, XYZ )
C
C              RECHERCHE DES POINTS SUR LES FACES DU PENTAEDRE COURBE
               CALL PEN5FA( XYZ, I, J, K, M, N,
     %                      XYSTR1, XYSTR5, XYSCA1, COSO,  XYZF )
C
C              RECHERCHE DES POINTS SUR LES ARETES DU PENTAEDRE COURBE
               CALL PEN9AR( XYZ, M, N, XYSCA1, COSO,  XYZA )
C
C              LA FORMULE D'INTERPOLATION TRANSFINIE DE DEGRE 1
C              SUR LE PENTAEDRE UNITE EST CALCULEE EN CE POINT XYZ
C              CF L'EQUIVALENCE POUR CBA2, CBA3 ET Z
               Z1   = 1 - Z
               CBA1 = 1 - CBA2 - CBA3
               NUSO = NUSOPE(I,J,K)
C
               DO 20 L=1,3
                  COSO(L,NUSO) = Z1 * XYZF(L,1) + Z * XYZF(L,2)
     %           +CBA1*( XYZF(L,3)+XYZF(L,8)-XYZA(L,13)
     %                  -( Z1*( XYZA(L,1)+XYZA(L,11)-COSO(L, 1) )
     %                    +Z *( XYZA(L,2)+XYZA(L,12)-COSO(L,N4) ) ) )
     %           +CBA2*( XYZF(L,4)+XYZF(L,5)-XYZA(L,14)
     %                  -( Z1*( XYZA(L,3)+XYZA(L, 5)-COSO(L,N2) )
     %                    +Z *( XYZA(L,4)+XYZA(L, 6)-COSO(L,N5) ) ) )
     %           +CBA3*( XYZF(L,6)+XYZF(L,7)-XYZA(L,15)
     %                  -( Z1*( XYZA(L,7)+XYZA(L, 9)-COSO(L,N3) )
     %                    +Z *( XYZA(L,8)+XYZA(L,10)-COSO(L,N6) ) ) )
 20            CONTINUE
C
 70         CONTINUE
 80      CONTINUE
 90   CONTINUE
      RETURN
      END
