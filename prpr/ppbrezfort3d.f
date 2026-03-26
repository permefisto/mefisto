      PROGRAM PPBREZFORT3D
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCULER LES MATRICES tP P et tDP DP pour les POLYNOMES 
C ----- DE BREZZI-FORTIN integrees sur le tetraedre unite
C       les AFFICHER pour etre integrees dans f3mp1bp1 f3rp1bp1 f3sp1bp1
C
C       Version avec le degre de liberte v(BARYCENTRE)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: Alain PERRONNET LJLL UPMC et St Pierre du Perray Novembre 2008
C MODIF : Alain PERRONNET LJLL UPMC et St Pierre du Perray      Mai 2010
C2345X7..............................................................012
      DOUBLE PRECISION  La(2,2,2,4),   P(3,3,3, 5),
     %                  PP(5,5,5),    DP(3,3,3, 5,3),
     %                  LaDP(4,4,4),  PL(4,4,4),      TP(5),
     %                  TPP(5,5),     TDPDP(5,3,5,3), TPDP(5,3,5),
     %                  DP1bLa(5,3,4), TPLa(5,4)
      double precision a,b,c,d
C
C     LES COORDONNEES BARYCENTIQUES SUR LE TETRAEDRE DE REFERENCE
      DATA La/ 1D0, -1D0, -1D0, 0D0,   -1D0,  0D0, 0D0, 0D0,
     %         0D0,  1D0,  0D0, 0D0,    0D0,  0D0, 0D0, 0D0,
     %         0D0,  0D0,  1D0, 0D0,    0D0,  0D0, 0D0, 0D0,
     %         0D0,  0D0,  0D0, 0D0,    1D0,  0D0, 0D0, 0D0 /
C
C     LES 5 POLYNOMES DE BASE DE BREZZI-FORTIN SUR LE TETRAEDRE REFERENCE
      DATA P/ 1D0,-1D0,0D0,   -1D0,  0D0, 0D0,   0D0, 0D0,0D0,
     %       -1D0, 0D0,0D0,    0D0,-64D0,64D0,   0D0,64D0,0D0,
     %        0D0, 0D0,0D0,    0D0, 64D0, 0D0,   0D0, 0D0,0D0,

     %        0D0, 1D0,0D0,    0D0, 0D0,  0D0,   0D0, 0D0,0D0,
     %        0D0, 0D0,0D0,    0D0,-64D0,64D0,   0D0,64D0,0D0,
     %        0D0, 0D0,0D0,    0D0, 64D0, 0D0,   0D0, 0D0,0D0,

     %        0D0, 0D0,0D0,    1D0,  0D0, 0D0,   0D0, 0D0,0D0,
     %        0D0, 0D0,0D0,    0D0,-64D0,64D0,   0D0,64D0,0D0,
     %        0D0, 0D0,0D0,    0D0, 64D0, 0D0,   0D0, 0D0,0D0,

     %        0D0, 0D0,0D0,    0D0,  0D0, 0D0,   0D0, 0D0,0D0,
     %        1D0, 0D0, 0D0,   0D0,-64D0,64D0,   0D0,64D0,0D0,
     %        0D0, 0D0, 0D0,   0D0, 64D0, 0D0,   0D0, 0D0,0D0,

     %        0D0, 0D0,0D0,    0D0,  0D0,   0D0, 0D0,   0D0,0D0,
     %        0D0, 0D0,0D0,    0D0,256D0,-256D0, 0D0,-256D0,0D0,
     %        0D0, 0D0,0D0,    0D0,-256D0,  0D0, 0D0,   0D0,0D0 /
C
C---------- INITIALISATION DES DERIVEES DES POLYNOMES DE BREZZI FORTIN
      DO 20 K=1,3
         DO 10 I=1,5
            CALL PN3DDE( K, 2+1, P(1,1,1,I), DP(1,1,1,I,K) )
 10      CONTINUE
 20   CONTINUE
C
C-----------  CALCUL de Integrale( P ) _____________________________
      DO I=1,5
         CALL PN3DIN( 0, 2+1, P(1,1,1,I), TP(I) )
      ENDDO
      print *
      print 10005,(I,TP(I),i=1,5)
10005 FORMAT(5('  tP1b',I1,'=',D26.18))
      print 10056,(TP(i),i=1,5)
      print *,' 73.D0 / 2520.D0=',73.D0 / 2520.D0,
     %        '  16.D0 / 315.D0=',16.D0 / 315.D0
C
C-----------  CALCUL de Integrale( tP P ) __________________________
      DO 35 J=1,5
         DO 30 I=1,5
            CALL PN3DPR( 1, 3, P(1,1,1,I),  1, 3, P(1,1,1,J),  5, PP )
            CALL PN3DIN( 0, 5, PP,  TPP(I,J) )
 30      CONTINUE
 35   CONTINUE
      print *
      print 10035,((i,j,TPP(I,J),j=1,5),i=1,5)
10035 FORMAT(5('  tP1b',I1,' P1b',I1,'=',D26.18))
      print 10056,((TPP(i,j),i=1,5),j=1,5)
C
C     LES 4 COEFFICIENTS DE LA MATRICE TPP
      A = 29836D0 / 2494800D0
      B =  9046D0 / 2494800D0
      C =   956D0 /  155925D0
      D =  4096D0 /  155925D0
      print 100,a,b,c,d
 100  format(/'TP1bP1b: A=',D26.18,' B=',D26.18,' C=',D26.18,
     %        ' D=',D26.18)
C
C-----------  CALCUL de Integrale( tP DP ) ________________________
      DO J=1,5
         DO I=1,5
            DO K=1,3
            CALL PN3DPR( 1, 3, P(1,1,1,I),  1, 3, DP(1,1,1,J,K),
     %                   5, PP )
            CALL PN3DIN( 0, 5, PP,  TPDP(I,K,J) )
            IF( ABS(TPDP(I,K,J)) .LT. 1D-14 ) TPDP(I,K,J)=0D0
            ENDDO
         ENDDO
      ENDDO
      print *
      print 10037,(((i,k,j,TPDP(i,k,j),i=1,5),k=1,3),j=1,5)
10037 FORMAT(5('  tP1b',I1,' DP1b',I1,I1,'=',D26.18))
      print *
      print 10038,(((TPDP(i,k,j),i=1,5),k=1,3),j=1,5)
10038 FORMAT('     %',D26.18,',  ',D26.18,',')
C
C-----------  CALCUL de Integrale( tDP DP ) ________________________
      DO 50 J=1,5
         DO 45 I=1,5
            DO 43 L=1,3
               DO 41 K=1,3
            CALL PN3DPR( 1, 3, DP(1,1,1,I,K),  1, 3, DP(1,1,1,J,L),
     %                   5, PP )
            CALL PN3DIN( 0, 5, PP,  TDPDP(I,K,J,L) )
 41            CONTINUE
 43         CONTINUE
 45      CONTINUE
 50   CONTINUE
      print *
      print 10055,((((i,k,j,l,TDPDP(i,k,j,l),i=1,5),k=1,3),j=1,5),l=1,3)
10055 FORMAT(5('  tDP1b',I1,I1,' DP1b',I1,I1,'=',D26.18))
      print *
      print 10056,((((TDPDP(i,k,j,l),i=1,5),k=1,3),j=1,5),l=1,3)
10056 FORMAT('     %',D26.18,',  ',D26.18,',')
C
C-----------  CALCUL de Integrale( tDP Lambda ) _____________________
      DO 70 J=1,4
         DO 65 I=1,5
            DO 63 K=1,3
            CALL PN3DPR( 1, 3, DP(1,1,1,I,K),  1, 2, La(1,1,1,J),
     %                   4, LaDP )
            CALL PN3DIN( 0, 4, LaDP,  DP1bLa(I,K,J) )
            IF( ABS(DP1bLa(I,K,J)) .LT. 1D-14 ) DP1bLa(I,K,J)=0D0
 63         CONTINUE
 65      CONTINUE
 70   CONTINUE
      print *
      print 10070,(((i,k,j,DP1bLa(i,k,j),i=1,5),k=1,3),j=1,4)
10070 FORMAT(5('  DP1b',I1,I1,' La',I1,'=',D26.18))
      print *
      print 10056,(((DP1bLa(i,k,j),i=1,5),k=1,3),j=1,4)
C
C-----------  CALCUL de Integrale( tP Lambda ) ______________________
      DO 90 J=1,4
         DO 80 I=1,5
            CALL PN3DPR( 1, 3, P(1,1,1,I),  1, 2, La(1,1,1,J),  4, PL )
            CALL PN3DIN( 0, 4, PL,  TPLa(I,J) )
 80      CONTINUE
 90      CONTINUE
      print *
      print 10090,((i,j,TPLa(I,J),i=1,5),j=1,4)
10090 FORMAT(5('  tP1b',I1,' La',I1,'=',D26.18))
      print 10056,((TPLa(i,j),i=1,5),j=1,4)
      a = 17D0 / 1260D0
      b = 13D0 / 2520D0
      c =  4D0 /  315D0
      print 10091,a,b,c
10091 format(/'TP1bLa: A=',D26.18,' B=',D26.18,' C=',D26.18)
C
      STOP
      END
