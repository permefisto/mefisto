      PROGRAM PPBREZFORT3D
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCULER LES MATRICES tP P et tDP DP pour les POLYNOMES 
C ----- DE BREZZI-FORTIN integrees sur le tetraedre unite
C       les AFFICHER pour etre integrees dans f3mp1bp1 f3rp1bp1 f3sp1bp1
C
C       Version avec le degre de liberte Int v dX / mes e
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: Alain PERRONNET LJLL UPMC et St Pierre du Perray Novembre 2008
C2345X7..............................................................012
      DOUBLE PRECISION  La(2,2,2,4),   P(5,5,5, 5),
     %                  PP(9,9,9),    DP(5,5,5, 5,3),
     %                  LaDP(7,7,7),  PL(6,6,6),
     %                  TPP(5,5),     TDPDP(5,3,5,3),
     %                  TDPLa(5,3,4), TPLa(5,4)
      double precision a,b,c,d
C
C     LES COORDONNEES BARYCENTIQUES SUR LE TETRAEDRE DE REFERENCE
      DATA La/ 1D0, -1D0, -1D0, 0D0,   -1D0,  0D0, 0D0, 0D0,
     %         0D0,  1D0,  0D0, 0D0,    0D0,  0D0, 0D0, 0D0,
     %         0D0,  0D0,  1D0, 0D0,    0D0,  0D0, 0D0, 0D0,
     %         0D0,  0D0,  0D0, 0D0,    1D0,  0D0, 0D0, 0D0 /
C
C     LES 5 POLYNOMES DE BASE DE BREZZI-FORTIN SUR LE TETRAEDRE REFERENCE
      DATA P/ 1D0,-1D0,3*0D0,  -1D0, 19*0D0,
     %       -1D0, 4*0D0, 0D0,-210D0,210D0,2*0D0, 0D0,210D0,13*0D0,
     %        6*0D0,210D0,68*0D0,

     %        0D0,1D0,23*0D0,
     %        6*0D0,-210D0,210D0,3*0D0,210D0,13*0D0,
     %        6*0D0,210D0,68*0D0,

     %        5*0D0,1D0,19*0D0,
     %        6*0D0,-210D0,210D0,3*0D0,210D0,13*0D0,
     %        6*0D0,210D0,68*0D0,

     %        25*0D0,
     %        1D0,5*0D0,-210D0,210D0,3*0D0,210D0,13*0D0,
     %        6*0D0,210D0,68*0D0,

     %        31*0D0, 840D0,-840D0,3*0D0,-840D0,13*0D0,
     %        6*0D0,-840D0,68*0D0 /
C
C---------- INITIALISATION DES DERIVEES DES POLYNOMES DE BREZZI FORTIN
      DO 20 K=1,3
         DO 10 I=1,5
            CALL PN3DDE( K, 4+1, P(1,1,1,I), DP(1,1,1,I,K) )
 10      CONTINUE
 20   CONTINUE
C
C-----------  CALCUL de Integrale( tP P ) __________________________
      DO 35 J=1,5
         DO 30 I=1,5
            CALL PN3DPR( 0, 5, P(1,1,1,I),  0, 5, P(1,1,1,J),  9, PP )
            CALL PN3DIN( 0, 9, PP,  TPP(I,J) )
 30      CONTINUE
 35   CONTINUE
      print *
      print 10035,((i,j,TPP(I,J),j=1,5),i=1,5)
10035 FORMAT(5('  tP',I1,'P',I1,'=',D25.17))
      print 10056,((TPP(i,j),i=1,5),j=1,5)
C
C-----------  CALCUL de Integrale( tDP DP ) ________________________
      DO 50 J=1,5
         DO 45 I=1,5
            DO 43 L=1,3
               DO 41 K=1,3
            CALL PN3DPR( 0, 5, DP(1,1,1,I,K),  0, 5, DP(1,1,1,J,L),
     %                   9, PP )
            CALL PN3DIN( 0, 9, PP,  TDPDP(I,K,J,L) )
 41            CONTINUE
 43         CONTINUE
 45      CONTINUE
 50   CONTINUE
      print *
      print 10055,((((i,k,j,l,TDPDP(i,k,j,l),i=1,5),k=1,3),j=1,5),l=1,3)
10055 FORMAT(5('  tDP',I1,I1,'DP',I1,I1,'=',D25.17))
      print *
      print 10056,((((TDPDP(i,k,j,l),i=1,5),k=1,3),j=1,5),l=1,3)
10056 FORMAT('     %',D25.17,',  ',D25.17,',')
C
C-----------  CALCUL de Integrale( tDP Lambda ) _____________________
      DO 70 J=1,4
         DO 65 I=1,5
            DO 63 K=1,3
            CALL PN3DPR( 0, 5, DP(1,1,1,I,K),  0, 2, La(1,1,1,J),
     %                   7, LaDP )
            CALL PN3DIN( 0, 7, LaDP,  TDPLa(I,K,J) )
            IF( ABS(TDPLa(I,K,J)) .LT. 1D-14 ) TDPLa(I,K,J)=0D0
 63         CONTINUE
 65      CONTINUE
 70   CONTINUE
      print *
      print 10070,(((i,k,j,TDPLa(i,k,j),i=1,5),k=1,3),j=1,4)
10070 FORMAT(5('  tDP',I1,I1,'La',I1,'=',D25.17))
      print *
      print 10056,(((TDPLa(i,k,j),i=1,5),k=1,3),j=1,4)
C
C-----------  CALCUL de Integrale( tP Lambda ) ______________________
      DO 90 J=1,4
         DO 80 I=1,5
            CALL PN3DPR( 0, 5, P(1,1,1,I),  0, 2, La(1,1,1,J),  6, PL )
            CALL PN3DIN( 0, 6, PL,  TPLa(I,J) )
 80      CONTINUE
 90      CONTINUE
      print *
      print 10090,((i,j,TPLa(I,J),i=1,5),j=1,4)
10090 FORMAT(5('  tP',I1,'La',I1,'=',D25.17))
      print 10056,((TPLa(i,j),i=1,5),j=1,4)
C
      A    =  107D0 / 7920D0
      B    =   41D0 / 7920D0
      C    =   23D0 /  792D0
      D    =   28D0 /   99D0
      print 100,a,b,c,d
 100  format('A=',D25.17,' B=',D25.17,' C=',D25.17,' D=',D25.17)
C
      STOP
      END
