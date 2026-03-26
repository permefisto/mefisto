         SUBROUTINE  VOEXT6( M, NBSTTR, XYSFTR, COSO )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    MAILLER UN TETRAEDRE COURBE DEFINI PAR SES FACES
C -----    PAR INTERPOLATION TRANSFINIE DE DEGRE 1
C          (cf CRAS A. PERRONNET)
C
C ENTREES:
C --------
C M      : LE NOMBRE CONSTANT DE SOMMETS PAR ARETE DU TETRAEDRE
C NBSTTR : LE NOMBRE DE SOMMETS PAR TRIANGLE NBSTTR=M*(M+1)/2
C XYSFTR : LES 3 COORDONNEES DES SOMMETS DES 4 FACES DU TETRAEDRE UNITE
C
C ENTREE ET SORTIE :
C ------------------
C COSO   : LES 3 COORDONNEES DES SOMMETS DU TETRAEDRE COURBE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 1997
C23456---------------------------------------------------------------012
      REAL         XYSFTR(2,NBSTTR,4),
     %             COSO(3,*)
      REAL         XYZ(3,0:12), XYZF(3,12), XYZA(3,12)
      EQUIVALENCE (CBA2,XYZ(1,0)), (CBA3,XYZ(2,0)), (CBA4,XYZ(3,0))
C
C     LA FONCTION FORMULE DU NUMERO DES SOMMETS DU TRIANGLE
      NUSOTR(I,J)   = (I*I-I)/2 + J
C     LA FONCTION FORMULE DU NUMERO DES SOMMETS DU TETRAEDRE
      NUSOTE(I,J,K) = (I-1)*I*(I+1)/6 + (J+K-2)*(J+K-1)/2 + K
C
C     LE NOMBRE DE SOMMETS D'UNE FACE TRIANGULAIRE
      NBSTTR = NUSOTR(M,M)
C
C     LE NUMERO DES SOMMETS DU TETRAEDRE
      N2 = NUSOTE(M,1,1)
      N3 = NUSOTE(M,M,1)
      N4 = NUSOTE(M,1,M)
C
C     LA BOUCLE SUR LES SOMMETS INTERNES
      DO 90 I = 4, M-1
         DO 80 J = 2, I-2
            DO 70 K = 2, I-J
C
C              LA DROITE D1: SF3(I,I+1-K)-SF2(I,K)
C              LE SOMMET SF3(I,I+1-K)
               L = NUSOTR(I,I+1-K)
               XYZ(1,1) = XYSFTR(2,L,3)
               XYZ(2,1) = 0
               XYZ(3,1) = XYSFTR(1,L,3)
C              LE SOMMET SF2(I,K)
               L = NUSOTR(I,K)
               XYZ(1,7) = 0
               XYZ(2,7) = XYSFTR(1,L,2)
               XYZ(3,7) = XYSFTR(2,L,2)
C
C              LA DROITE D2: SF1(I,J)-SF2(I,I+1-J)
C              LE SOMMET SF1(I,J)
               L = NUSOTR(I,J)
               XYZ(1,2) = XYSFTR(1,L,1)
               XYZ(2,2) = XYSFTR(2,L,1)
               XYZ(3,2) = 0
C              LE SOMMET SF2(I,I+1-J)
               L = NUSOTR(I,I+1-J)
               XYZ(1,8) = 0
               XYZ(2,8) = XYSFTR(1,L,2)
               XYZ(3,8) = XYSFTR(2,L,2)
C
C              LA DROITE D3: SF2(K+J-1,K)-SF4(J+K-1,K)
C              LE SOMMET SF2(K+J-1,K)
               L = NUSOTR(K+J-1,K)
               XYZ(1,3) = 0
               XYZ(2,3) = XYSFTR(1,L,2)
               XYZ(3,3) = XYSFTR(2,L,2)
C              LE SOMMET SF4(J+K-1,K)
               L = NUSOTR(J+K-1,K)
               XYZ(1,9) = 1-XYSFTR(1,L,4)-XYSFTR(2,L,4)
               XYZ(2,9) = XYSFTR(1,L,4)
               XYZ(3,9) = XYSFTR(2,L,4)
C
C              LA DROITE D4: SF1(I,J+K-1)-SF3(I,I+2-J-K)
C              LE SOMMET SF1(I,J+K-1)
               L = NUSOTR(I,J+K-1)
               XYZ(1,4) = XYSFTR(1,L,1)
               XYZ(2,4) = XYSFTR(2,L,1)
               XYZ(3,4) = 0
C              LE SOMMET SF3(I,I+2-J-K)
               L = NUSOTR(I,I+2-J-K)
               XYZ(1,10) = XYSFTR(2,L,3)
               XYZ(2,10) = 0
               XYZ(3,10) = XYSFTR(1,L,3)
C
C              LA DROITE D5: SF3(I-J+1,I+2-K-J)-SF4(M-1-I+K+J,K)
C              LE SOMMET SF3(I-J+1,I+2-K-J)
               L = NUSOTR(I-J+1,I+2-K-J)
               XYZ(1,5) = XYSFTR(2,L,3)
               XYZ(2,5) = 0
               XYZ(3,5) = XYSFTR(1,L,3)
C              LE SOMMET SF4(M-1-I+K+J,K)
               L = NUSOTR(M-1-I+K+J,K)
               XYZ(1,11) = 1-XYSFTR(1,L,4)-XYSFTR(2,L,4)
               XYZ(2,11) = XYSFTR(1,L,4)
               XYZ(3,11) = XYSFTR(2,L,4)
C
C              LA DROITE D6: SF1(I-K+1,J)-SF4(M-I+K+J-1,M-I+K)
C              LE SOMMET SF1(I-K+1,J)
               L = NUSOTR(I-K+1,J)
               XYZ(1,6) = XYSFTR(1,L,1)
               XYZ(2,6) = XYSFTR(2,L,1)
               XYZ(3,6) = 0
C              LE SOMMET SF4(M-I+K+J-1,M-I+K)
               L = NUSOTR(M-I+K+J-1,M-I+K)
               XYZ(1,12) = 1-XYSFTR(1,L,4)-XYSFTR(2,L,4)
               XYZ(2,12) = XYSFTR(1,L,4)
               XYZ(3,12) = XYSFTR(2,L,4)
C
C              RECHERCHE DU POINT A DISTANCE MINIMALE DE CES 6 DROITES
               CALL MINDIS( 6, XYZ )
C
C              RECHERCHE DES POINTS SUR LES FACES DU TETRAEDRE COURBE
               CALL TET4FA( XYZ, M, NBSTTR, XYSFTR, COSO,
     %                      XYZF )
C
C              RECHERCHE DES POINTS SUR LES ARETES DU TETRAEDRE COURBE
               CALL TET6AR( XYZ, M, NBSTTR, XYSFTR, COSO,
     %                      XYZA )
C
C              LA FORMULE D'INTERPOLATION TRANSFINIE DE DEGRE 1
C              SUR LE TETRAEDRE UNITE EST CALCULEE EN CE POINT XYZ
C              CF L'EQUIVALENCE POUR CBA2, CBA3 ET CB4
               CBA1 = 1 - CBA2 - CBA3 - CBA4
               NUSO = NUSOTE(I,J,K)
C
               DO 20 L=1,3
               COSO(L,NUSO) =
     %         CBA1*( XYZF(L,1) + XYZF(L, 9) + XYZF(L,11)
     %               -XYZA(L,1) - XYZA(L, 6) - XYZA(L, 7) + COSO(L,1) )
     %        +CBA2*( XYZF(L,2) + XYZF(L, 4) + XYZF(L,12)
     %               -XYZA(L,2) - XYZA(L, 3) - XYZA(L, 9) + COSO(L,N2) )
     %        +CBA3*( XYZF(L,3) + XYZF(L, 5) + XYZF(L, 7)
     %               -XYZA(L,4) - XYZA(L, 5) - XYZA(L,11) + COSO(L,N3) )
     %        +CBA4*( XYZF(L,6) + XYZF(L, 8) + XYZF(L,10)
     %               -XYZA(L,8) - XYZA(L,10) - XYZA(L,12) + COSO(L,N4) )
 20            CONTINUE
C
 70         CONTINUE
 80      CONTINUE
 90   CONTINUE
      RETURN
      END
