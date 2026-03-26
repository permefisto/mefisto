      PROGRAM PPDVSR
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCULER LES 10 COEFFICIENTS DES 3 POLYNOMES P3 DES 3 TRIANGLES
C ----- DU CARRE UNITE POUR CHACUNE DES 12 FONCTIONS DE BASE
C       DU CARRE DE de VEUBKE SANDER REDUIT
C       (DERIVEE NORMALE P1 SUR LES 4 COTES)
C
C             S4                     S3
C              X----->------------<---X
C              | \     p1+p2+p3     / |
C              |   \              /   |
C              -     \          /     -
C              |       \      /       |
C              |          \ /         |
C              | p1+p3   /   \  p1+p2 |
C              |       /       \      |
C              ^     /           \    ^
C              |   /               \  |
C              | /        p1         \|
C              X----->----------<-----X
C             S1                     S2
C
C SORTIE :
C --------
C FdVeSa : FdVeSa(10,3,12) TABLEAU DES COEFFICIENTS DES FONCTIONS DE BASE
C          FdVeSa(I,J,K) I=NUMERO DU COEFFICIENT DU POLYNOME (1 a 10)
C                        J=NUMERO DU POLYNOME pj
C                        K=NUMERO DU DEGRE DE LIBERTE
C
C AVEC L'ORDRE SUIVANT POUR CHAQUE INDICE :
C P(U,V) = SOMME pkl U**k V**l = SOMME pi U**k V**l
C I=1 p00, I=2 p10, I=3 p01,  I=4 p11, I=5 p20, I=6 p02,
C I=7 p21, I=8 p12, I=9 p30, I=10 p03
C
C J=1 2 3 COMME INDIQUE SUR LA FIGURE CI DESSUS
C J=1 => P1 AGIT SUR LES 4 TRIANGLES
C J=2 => P2 AGIT SUR LE TRIANGLE S2S3S4 (P2=0 SUR ARETE S2S4)
C J=3 => P3 AGIT SUR LE TRIANGLE S1S3S4 (P3=0 SUR ARETE S1S3)
C
C K=1  F(S1),        K=2   F(S2),        K=3   F(S3),        K=4   F(S4),
C K=5 DF(S1)(S2-S1), K=6  DF(S1)(S4-S1), K=7  DF(S2)(S3-S2), K=8  DF(S2)(S1-S2),
C K=9 DF(S3)(S4-S3), K=10 DF(S3)(S2-S3), K=11 DF(S4)(S1-S4), K=12 DF(S4)(S3-S4)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : OUACHTAOUI DEA  ANALYSE NUMERIQUE UPMC PARIS     FEVRIER 1996
C MODIFS : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS         MAI 1996
C2345X7..............................................................012
      DOUBLE PRECISION  A(30,42),
     %                  FdVeSa(10,3,12)
      EQUIVALENCE      (A(1,31),FdVeSa(1,1,1))
      REAL              DUDVSR(12), DUHCT(9)
C
C---------- INITIALISATION DE LA MATRICE ET DU SECOND MEMBRE --
C
      DO 20 J=1,42
         DO 10 I=1,30
            A(I,J)=0.D0
 10      CONTINUE
 20   CONTINUE
C
C-----------  REMPLISSAGE DES 12 SECONDS MEMBRES ---------------
C
      DO 30 I=1,12
         A(I,30+I)=1D0
 30   CONTINUE
C
C-----------  REMPLISSAGE DE LA MATRICE  -----------------------
C
C     LIGNE 1 = x1,y1,z1
      A(1,1)=1.D0
C								
C     LIGNE 2 = x2,y2,z2
      A(2,1)=1.D0
      A(2,2)=1.D0
      A(2,5)=1.D0
      A(2,9)=1.D0
C	
C     LIGNE 3 = x3,y3,z3
      DO 2 I=1,20
         A(3,I)=1.D0
 2    CONTINUE
C
C     LIGNE 4 = x4,y4,z4
      A(4, 1)=1.D0
      A(4, 3)=1.D0
      A(4, 6)=1.D0
      A(4,10)=1.D0
      A(4,21)=1.D0
      A(4,23)=1.D0
      A(4,26)=1.D0
      A(4,30)=1.D0
C	
C     LIGNE 5  = a la composante de D(s1)(s2-s1) en x1,y1,z1
      A(5,2)=1.D0
C
C     LIGNE 6  = a la composante de D(s1)(s4-s1) en x1,y1,z1
      A(6,3)=1.D0
C      	
C     LIGNE 7  = a la composante de D(s2)(s3-s2) en x2,y2,z2
      A(7,3)=1.D0
      A(7,4)=1.D0
      A(7,7)=1.D0
C
C     LIGNE 8  = a la composante de D(s2)(s1-s2) en x2,y2,z2
      A(8,2)=-1.D0
      A(8,5)=-2.D0
      A(8,9)=-3.D0
C
C     LIGNE 9  = a la composante de D(s3)(s4-s3) en x3,y3,z3
      A(9, 2)=-1.D0
      A(9, 4)=-1.D0
      A(9, 5)=-2.D0
      A(9, 7)=-2.D0
      A(9, 8)=-1.D0
      A(9, 9)=-3.D0
      A(9,12)=-1.D0
      A(9,14)=-1.D0
      A(9,15)=-2.D0
      A(9,17)=-2.D0	
      A(9,18)=-1.D0	
      A(9,19)=-3.D0
C
C     LIGNE 10 = a la composante de D(s3)(s2-s3) en x3,y3,z3
      A(10, 3)=-1.D0
      A(10, 4)=-1.D0
      A(10, 6)=-2.D0
      A(10, 7)=-1.D0
      A(10, 8)=-2.D0
      A(10,10)=-3.D0
      A(10,13)=-1.D0
      A(10,14)=-1.D0
      A(10,16)=-2.D0
      A(10,17)=-1.D0
      A(10,18)=-2.D0
      A(10,20)=-3.D0
C
C     LIGNE 11 = a la composante de D(s4)(s1-s4) en x4,y4,z4
      A(11, 3)=-1.D0
      A(11, 6)=-2.D0
      A(11,10)=-3.D0
      A(11,23)=-1.D0
      A(11,26)=-2.D0
      A(11,30)=-3.D0
C
C     LIGNE 12 = a la composante de D(s4)(s3-s4) en x4,y4,z4
      A(12, 2)=1.D0
      A(12, 4)=1.D0
      A(12, 8)=1.D0
      A(12,22)=1.D0
      A(12,24)=1.D0
      A(12,28)=1.D0
C
C     LIGNE 13
      A(13,11)=1.D0
      A(13,12)=1.D0
      A(13,15)=1.D0	
      A(13,19)=1.D0
C
C     LIGNE 14
      A(14,11)=1.D0	
      A(14,13)=1.D0	
      A(14,16)=1.D0	
      A(14,20)=1.D0
C	
C     LIGNE 15
      DO 3 I=11,20
         IF(I.LE.11) A(15,I)=8.D0
         IF((I.GT.11).AND.(I.LE.13)) A(15,I)=4.D0
         IF((I.GT.13).AND.(I.LE.16)) A(15,I)=2.D0
         IF(I.GT.16) A(15,I)=1.D0
 3    CONTINUE
C
C     LIGNE 16
      A(16,11)=64.D0	
      A(16,12)=48.D0
      A(16,13)=16.D0	
      A(16,14)=12.D0
      A(16,15)=36.D0
      A(16,16)=4.D0
      A(16,17)=9.D0
      A(16,18)=3.D0
      A(16,19)=27.D0
      A(16,20)=1.D0
C
C     LIGNE 17
      A(17,12)=1.D0	
      A(17,13)=1.D0	
      A(17,14)=1.D0	
      A(17,15)=2.D0	
      A(17,17)=1.D0	
      A(17,19)=3.D0
C	
C     LIGNE 18
      A(18,12)=1.D0
      A(18,13)=1.D0
      A(18,14)=1.D0
      A(18,16)=2.D0
      A(18,18)=1.D0
      A(18,20)=3.D0
C
C     LIGNE 19
      A(19,12)=4.D0
      A(19,13)=4.D0
      A(19,14)=4.D0
      A(19,15)=4.D0
      A(19,16)=4.D0
      A(19,17)=3.D0
      A(19,18)=3.D0
      A(19,19)=3.D0
      A(19,20)=3.D0
C
C     LIGNE 20
      A(20,21)=1.D0
C	
C     LIGNE 21
      DO 4 I=21,30
         A(21,I)=1.D0
 4    CONTINUE
C
C     LIGNE 22
      DO 5 I=21,30
         IF(I.LE.21) A(22,I)=8.D0
         IF((I.GT.21).AND.(I.LE.23)) A(22,I)=4.D0
         IF((I.GT.23).AND.(I.LE.26)) A(22,I)=2.D0
         IF(I.GT.26) A(22,I)=1.D0
 5    CONTINUE
C
C     LIGNE 23
      DO 6 I=21,30
         IF(I.LE.21) A(23,I)=64.D0
         IF((I.GT.21).AND.(I.LE.23)) A(23,I)=16.D0
         IF((I.GT.23).AND.(I.LE.26)) A(23,I)=4.D0
         IF(I.GT.26) A(23,I)=1.D0
 6    CONTINUE
C
C     LIGNE 24
      A(24,22)=1.D0
      A(24,23)=-1.D0
C
C     LIGNE 25
      A(25,22)=1.D0
      A(25,23)=-1.D0	
      A(25,25)=2.D0
      A(25,26)=-2.D0
      A(25,27)=1.D0
      A(25,28)=-1.D0
      A(25,29)=3.D0
      A(25,30)=-3.D0
C
C     LIGNE 26
      A(26,22)=4.D0
      A(26,23)=-4.D0
      A(26,25)=4.D0
      A(26,26)=-4.D0
      A(26,27)=1.D0
      A(26,28)=-1.D0
      A(26,29)=3.D0
      A(26,30)=-3.D0
C
C     LIGNE 27
      A(27,7)=1.D0
C
C     LIGNE 28
      A(28,8)=1.D0
      A(28,18)=1.D0
C
C     LIGNE 29
      A(29,7)=1.D0
      A(29,17)=1.D0
      A(29,27)=1.D0
C
C     LIGNE 30
      A(30,8)=1.D0
      A(30,28)=1.D0
C
C     RESOLUTION DU SYSTEME AVEC 12 SECONDS MEMBRES
      CALL GAUSPT(30,12,A,IER)
C
C     RESULTAT MIS SUR UN FICHIER
      OPEN(UNIT=20,FILE='incl/fdvesa.inc')
C
      WRITE(20,10000)
      WRITE(* ,10000)
10000 FORMAT('       DOUBLE PRECISION  FDVESA(10,3,12)')
      WRITE(20,10001)
      WRITE(* ,10001)
10001 FORMAT('       COMMON / FDBDVS / FDVESA')
C
C     AFFICHAGE ET STOCKAGE
      DO 310 K=1,3
C        LE SOUS-TRIANGLE K
         DO 300 J=1,12
C           LE DEGRE DE LIBERTE J
10300 FORMAT('      FBDVSR(',I2,')= FBDVSR(',I2,')')
10301 FORMAT('     & +',D25.17)
10302 FORMAT('     & +',D25.17,' * U' )
10303 FORMAT('     & +',D25.17,' * V')
10304 FORMAT('     & +',D25.17,' * U * V')
10305 FORMAT('     & +',D25.17,' * U * U')
10306 FORMAT('     & +',D25.17,' * V * V')
10307 FORMAT('     & +',D25.17,' * U * U * V')
10308 FORMAT('     & +',D25.17,' * U * V * V')
10309 FORMAT('     & +',D25.17,' * U * U * U')
10310 FORMAT('     & +',D25.17,' * V * V * V')
C
      WRITE(20,10300) J,J
      IF( ABS(FdVeSa(1,K,J)) .GT. 1D-14 )
     %    WRITE(20,10301) FdVeSa(1,K,J)
      IF( ABS(FdVeSa(2,K,J)) .GT. 1D-14 )
     %    WRITE(20,10302) FdVeSa(2,K,J)
      IF( ABS(FdVeSa(3,K,J)) .GT. 1D-14 )
     %    WRITE(20,10303) FdVeSa(3,K,J)
      IF( ABS(FdVeSa(4,K,J)) .GT. 1D-14 )
     %    WRITE(20,10304) FdVeSa(4,K,J)
      IF( ABS(FdVeSa(5,K,J)) .GT. 1D-14 )
     %    WRITE(20,10305) FdVeSa(5,K,J)
      IF( ABS(FdVeSa(6,K,J)) .GT. 1D-14 )
     %    WRITE(20,10306) FdVeSa(6,K,J)
      IF( ABS(FdVeSa(7,K,J)) .GT. 1D-14 )
     %    WRITE(20,10307) FdVeSa(7,K,J)
      IF( ABS(FdVeSa(8,K,J)) .GT. 1D-14 )
     %    WRITE(20,10308) FdVeSa(8,K,J)
      IF( ABS(FdVeSa(9,K,J)) .GT. 1D-14 )
     %    WRITE(20,10309) FdVeSa(9,K,J)
      IF( ABS(FdVeSa(10,K,J)) .GT. 1D-14 )
     %    WRITE(20,10310) FdVeSa(10,K,J)
C
 300     CONTINUE
 310  CONTINUE
C
      CLOSE( UNIT=20 )
C
C     VERIFICATION DU DVS REDUIT
10400 FORMAT('AU POINT',F3.1,' ',F3.1,'  DVS   ='/12F6.2/)
      CALL VFBDVS( 0., 0., DUDVSR )
      PRINT 10400,0.,0.,DUDVSR
      CALL VFBDVS( 1., 0., DUDVSR )
      PRINT 10400,1.,0.,DUDVSR
      CALL VFBDVS( 1., 1., DUDVSR )
      PRINT 10400,1.,1.,DUDVSR
      CALL VFBDVS( 0., 1., DUDVSR )
      PRINT 10400,0.,1.,DUDVSR
      PRINT *
      PRINT *
C
10410 FORMAT('AU POINT',F3.1,' ',F3.1,' DDVS/DU ='/12F6.2/)
      CALL DDUDVS( 0., 0., DUDVSR )
      PRINT 10410,0.,0.,DUDVSR
      CALL DDUDVS( 1., 0., DUDVSR )
      PRINT 10410,1.,0.,DUDVSR
      CALL DDUDVS( 1., 1., DUDVSR )
      PRINT 10410,1.,1.,DUDVSR
      CALL DDUDVS( 0., 1., DUDVSR )
      PRINT 10410,0.,1.,DUDVSR
      PRINT *
C
10420 FORMAT('AU POINT',F3.1,' ',F3.1,' DDVS/DV ='/12F6.2/)
      CALL DDVDVS( 0., 0., DUDVSR )
      PRINT 10420,0.,0.,DUDVSR
      CALL DDVDVS( 1., 0., DUDVSR )
      PRINT 10420,1.,0.,DUDVSR
      CALL DDVDVS( 1., 1., DUDVSR )
      PRINT 10420,1.,1.,DUDVSR
      CALL DDVDVS( 0., 1., DUDVSR )
      PRINT 10420,0.,1.,DUDVSR
      PRINT *
      PRINT *
      PRINT *
      PRINT *
C
C
C
C     VERIFICATION DU HCT REDUIT
10500 FORMAT('AU POINT',F3.1,' ',F3.1,'  HCT   ='/9F6.2/)
      CALL VFBHCT( 0., 0., DUHCT )
      PRINT 10500,0.,0.,DUHCT
      CALL VFBHCT( 1., 0., DUHCT )
      PRINT 10500,1.,0.,DUHCT
      CALL VFBHCT( 0., 1., DUHCT )
      PRINT 10500,0.,1.,DUHCT
      PRINT *
      PRINT *
C
10510 FORMAT('AU POINT',F3.1,' ',F3.1,' DHCT/DU ='/9F6.2/)
      CALL DDUHCT( 0., 0., DUHCT )
      PRINT 10510,0.,0.,DUHCT
      CALL DDUHCT( 1., 0., DUHCT )
      PRINT 10510,1.,0.,DUHCT
      CALL DDUHCT( 0., 1., DUHCT )
      PRINT 10510,0.,1.,DUHCT
      PRINT *
C
      CALL DDUHCT( 0.5, 0., DUHCT )
      PRINT 10510,0.5,0.,DUHCT
      CALL DDUHCT( 0.5, 0.5, DUHCT )
      PRINT 10510,0.5, 0.5,DUHCT
      CALL DDUHCT( 0., 0.5, DUHCT )
      PRINT 10510,0.,0.5,DUHCT
      PRINT *
C
10520 FORMAT('AU POINT',F3.1,' ',F3.1,' DHCT/DV ='/9F6.2/)
      CALL DDVHCT( 0., 0., DUHCT )
      PRINT 10520,0.,0.,DUHCT
      CALL DDVHCT( 1., 0., DUHCT )
      PRINT 10520,1.,0.,DUHCT
      CALL DDVHCT( 0., 1., DUHCT )
      PRINT 10520,0.,1.,DUHCT
      PRINT *
C
      CALL DDVHCT( 0.5, 0., DUHCT )
      PRINT 10520,0.5,0.,DUHCT
      CALL DDVHCT( 0.5, 0.5, DUHCT )
      PRINT 10520,0.5, 0.5,DUHCT
      CALL DDVHCT( 0., 0.5, DUHCT )
      PRINT 10520,0.,0.5,DUHCT
      PRINT *
      END
