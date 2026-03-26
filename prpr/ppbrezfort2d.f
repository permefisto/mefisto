      PROGRAM PPBREZFORT2D
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCULER LES MATRICES tP P et tDP DP pour les POLYNOMES 
C ----- DE BREZZI-FORTIN integrees sur le triangle unite
C       les AFFICHER pour etre integrees dans f2mp1bp1 f2rp1bp1 f2sp1bp1
C
C       version avec le dl v(BARYCENTRE)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: Alain PERRONNET LJLL UPMC et St Pierre du Perray Novembre 2008
C MODIF : Alain PERRONNET LJLL UPMC et St Pierre du Perray      Mai 2010
C2345X7..............................................................012
      DOUBLE PRECISION  PP(7,7),       La(2,2,3),
     %                  P(4,4,4),      DP(4,4,4,2),
     %                  LaDP(5,5),     PL(5,5),
     %                  TPP(4,4),      TDPDP(4,2,4,2), TPDP(4,2,4),
     %                  DP1BLa(4,2,3), TPLa(4,3),
     %                  A, B, C, D
C
C     LES COORDONNEES BARYCENTIQUES SUR LE TRIANGLE DE REFERENCE
      DATA La/ 1D0, -1D0, -1D0, 0D0,
     %         0D0,  1D0,  0D0, 0D0,
     %         0D0,  0D0,  1D0, 0D0 /
C
C     LES 4 POLYNOMES DE BASE DE BREZZI-FORTIN SUR LE TRIANGLE REFERENCE
      DATA P/ 1D0,-1D0,0D0,0D0, -1D0,-9D0,  9D0,0D0,  0D0,  9D0,0D0,0D0,
     %        0D0, 0D0,0D0,0D0,
     %        0D0, 1D0,0D0,0D0,  0D0,-9D0,  9D0,0D0,  0D0,  9D0,0D0,0D0,
     %        0D0, 0D0,0D0,0D0,
     %        0D0, 0D0,0D0,0D0,  1D0,-9D0,  9D0,0D0,  0D0,  9D0,0D0,0D0,
     %        0D0, 0D0,0D0,0D0,
     %        0D0, 0D0,0D0,0D0,  0D0,27D0,-27D0,0D0,  0D0,-27D0,0D0,0D0,
     %        0D0, 0D0,0D0,0D0 /
C
C---------- INITIALISATION DES DERIVEES DES POLYNOMES DE BREZZI FORTIN
      DO 20 K=1,2
         DO 10 I=1,4
            CALL PN2DDE( K, 3+1, P(1,1,I), DP(1,1,I,K) )
 10      CONTINUE
 20   CONTINUE
C
C-----------  CALCUL des Integrale( tP P ) __________________________
      DO 35 J=1,4
         DO 30 I=1,4
            CALL PN2DPR( 0, 4, P(1,1,I),  0, 4, P(1,1,J),  7, PP )
            CALL PN2DIN( 0, 7, PP,  TPP(I,J) )
 30      CONTINUE
 35   CONTINUE
      print *
      print 10035,((i,j,TPP(I,J),j=1,4),i=1,4)
10035 FORMAT(4('  tP',I1,' P',I1,'=',D26.18))
      A    = 83D0 / 1680D0
      B    = 13D0 / 1680D0
      C    =  3D0 /  112D0
      D    = 81D0 /  560D0
      print 10036,A,B,C,D
10036 format('TP4P4: A=',D26.18,' B=',D26.18,' C=',D26.18,' D+',D26.18)
C
C-----------  CALCUL des Integrale( tP DP ) _________________________
      DO J=1,4
         DO I=1,4
            DO K=1,2
            CALL PN2DPR( 0, 4, P(1,1,I),  0, 4, DP(1,1,J,K),  7, PP )
            CALL PN2DIN( 0, 7, PP,  TPDP(I,K,J) )
            IF( ABS(TPDP(I,K,J)) .LT. 1D-14 ) TPDP(I,K,J)=0D0
            ENDDO
         ENDDO
      ENDDO
      print *
      print 10037,(((i,k,j,TPDP(i,k,j),i=1,4),k=1,2),j=1,4)
10037 FORMAT(4('  tP',I1,' DP',I1,I1,'=',D26.18))
      print *
      print 10038,(((TPDP(i,k,j),i=1,4),k=1,2),j=1,4)
10038 FORMAT('     %',D26.18,',  ',D26.18,',')
C
C-----------  CALCUL des Integrale( tDP DP ) ________________________
      DO 50 J=1,4
         DO 45 I=1,4
            DO 43 L=1,2
               DO 41 K=1,2
            CALL PN2DPR( 0, 4, DP(1,1,I,K),  0, 4, DP(1,1,J,L),  7, PP )
            CALL PN2DIN( 0, 7, PP,  TDPDP(I,K,J,L) )
 41            CONTINUE
 43         CONTINUE
 45      CONTINUE
 50   CONTINUE
      print *
      print 10055,((((i,k,j,l,TDPDP(i,k,j,l),i=1,4),k=1,2),j=1,4),l=1,2)
10055 FORMAT(4('  tDP',I1,I1,' DP',I1,I1,'=',D26.18))
      print *
      print 10056,((((TDPDP(i,k,j,l),i=1,4),k=1,2),j=1,4),l=1,2)
10056 FORMAT('     %',D26.18,',  ',D26.18,',')
C
C-----------  CALCUL des Integrale( tDP Lambda ) _____________________
      DO 70 J=1,3
         DO 65 I=1,4
            DO 63 K=1,2
            CALL PN2DPR( 0, 4, DP(1,1,I,K),  0, 2, La(1,1,J),  5, LaDP )
            CALL PN2DIN( 0, 5, LaDP,  DP1BLa(I,K,J) )
            IF( ABS(DP1BLa(I,K,J)) .LT. 1D-14 ) DP1BLa(I,K,J)=0D0
 63         CONTINUE
 65      CONTINUE
 70   CONTINUE
      print *
      print 10070,(((i,k,j,DP1BLa(i,k,j),i=1,4),k=1,2),j=1,3)
10070 FORMAT(4('  tDP1B',I1,I1,' La',I1,'=',D26.18))
      print *
      print 10056,(((DP1BLa(i,k,j),i=1,4),k=1,2),j=1,3)
C
C-----------  CALCUL des Integrale( tP Lambda ) ______________________
      DO 90 J=1,3
         DO 80 I=1,4
            CALL PN2DPR( 0, 4, P(1,1,I),  0, 2, La(1,1,J),  5, PL )
            CALL PN2DIN( 0, 5, PL,  TPLa(I,J) )
 80      CONTINUE
 90      CONTINUE
      print *
      print 10090,((i,j,TPLa(I,J),i=1,4),j=1,3)
10090 FORMAT(4('  tP',I1,' La',I1,'=',D26.18))
      print 10056,((TPLa(i,j),i=1,4),j=1,3)
C
      STOP
      END
