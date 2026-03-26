      SUBROUTINE PN3DPR( NPOUQ1,N1,P1, NPOUQ2,N2,P2, N3,P3 )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : RECHERCHE DES COEFFICIENTS DU POLYNOME PRODUIT P3 DE P1 PAR P2
C -----
C
C PARAMETRES D ENTREE :
C ---------------------
C NPOUQI : 0 SI LE POLYNOME PI EST DE DEGRE(NI-1) EN (X,Y,Z)
C        : 1 SI LE POLYNOME PI EST DE DEGRE(NI-1) EN X ET Y ET Z
C        : 2 SI LE POLYNOME PI EST DE DEGRE(NI-1) EN (X,Y) ET Z
C NI     : DEGRE+1 DU POLYNOME PI AU SENS CI DESSUS
C PI     : TABLEAU A 3 INDICES CONTENANT LES COEFFICIENTS DU POLYNOME PI
C          PI(I,J,K)=COEFFICIENT DE X**(I-1) Y**(J-1) Z**(K-1)
C          I= 1 OU 2
C N3     : DIMENSION DU POLYNOME PRODUIT P3
C
C PARAMETRE-RESULTAT :
C --------------------
C P3     : LE POLYNOME PRODUIT DE P1 PAR P2
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : A.PERRONNET LAN189 PARIS ET IRIA  OCTOBRE 1979
C ......................................................................
      DOUBLE PRECISION P1(N1,N1,N1), P2(N2,N2,N2), P3(N3,N3,N3)
C
C     INITIALISATION A ZERO DE P3
C
      DO I=1,N3
         DO J=1,N3
            DO K=1,N3
               P3(I,J,K) = 0.D0
            ENDDO
         ENDDO
      ENDDO
C
C     QUELQUES VARIABLES AUXILIAIRES
C
      N12 = N1  + 1
      N13 = N12 + 1
      N22 = N2  + 1
      N23 = N22 + 1
      NPOU1 = NPOUQ1 - 1
      NPOU2 = NPOUQ2 - 1
C
      DO 10 I=1,N1
         J1 = N1
         IF(NPOU1 .NE. 0) J1 = N12 - I
         II1 = I - 1
         DO 20 J=1,J1
            K1 = N1
            IF(NPOUQ1 .EQ. 0) K1 = N13 - I - J
            JJ1 = J - 1
            DO 30 K=1,K1
               KK1 = K - 1
               DO 40 IP=1,N2
                  J1P = N2
                  IF(NPOU2 .NE. 0) J1P = N22 - IP
                  IIP = II1 + IP
                  DO 50 JP=1,J1P
                     K1P = N2
                     IF(NPOUQ2 .EQ. 0) K1P = N23 - IP - JP
                     JJP = JJ1 + JP
                     DO 60 KP=1,K1P
                        KKP = KK1 + KP
                        P3(IIP,JJP,KKP) = P3(IIP,JJP,KKP)
     &                                  + P1(I,J,K) * P2(IP,JP,KP)
 60                  CONTINUE
 50               CONTINUE
 40            CONTINUE
 30         CONTINUE
 20      CONTINUE
 10   CONTINUE
C
      RETURN
      END
