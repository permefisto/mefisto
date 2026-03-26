      SUBROUTINE PN2DPR( NPOUQ1,N1,P1, NPOUQ2,N2,P2, N3,P3 )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : RECHERCHE DES COEFFICIENTS DU POLYNOME PRODUIT P3 DE P1 PAR P2
C -----

C PARAMETRES D ENTREE :
C ---------------------
C NPOUQI : 0 SI LE POLYNOME PI EST DE DEGRE(N1-1) EN (X,Y)
C        : 1 SI LE POLYNOME PI EST DE DEGRE(N1-1) EN X ET Y
C NI     : DEGRE+1 DU POLYNOME PI AU SENS CI DESSUS
C PI     : TABLEAU A 3 INDICES CONTENANT LES COEFFICIENTS DU POLYNOME PI
C          PI(I,J)=COEFFICIENT DE X**(I-1) Y**(J-1)
C          I= 1 OU 2
C N3     : DIMENSION DU POLYNOME PRODUIT P3

C PARAMETRE-RESULTAT :
C --------------------
C P3     : LE POLYNOME PRODUIT DE P1 PAR P2
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : A.PERRONNET LAN189 PARIS ET IRIA  OCTOBRE 1979
C ......................................................................
      DOUBLE PRECISION P1(N1,N1), P2(N2,N2), P3(N3,N3)

C     INITIALISATION A ZERO DE P3
      DO I=1,N3
         DO J=1,N3
            P3(I,J) = 0.D0
         ENDDO
      ENDDO

C     QUELQUES VARIABLES AUXILIAIRES
      N12 = N1  + 1
      N22 = N2  + 1
      NPOU1 = NPOUQ1 - 1
      NPOU2 = NPOUQ2 - 1

      DO I=1,N1
         J1 = N1
         IF(NPOU1 .NE. 0) J1 = N12 - I
         II1 = I - 1
         DO J=1,J1
            JJ1 = J - 1
            DO IP=1,N2
               J1P = N2
               IF(NPOU2 .NE. 0) J1P = N22 - IP
               IIP = II1 + IP
               DO JP=1,J1P
                  JJP = JJ1 + JP
                  P3(IIP,JJP) = P3(IIP,JJP)
     %                        + P1(I,J) * P2(IP,JP)
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      RETURN
      END
