      SUBROUTINE AEBFV3( Rho, DtMhu, NONOEF, NBSOM, XYZSOM,
     %                   AE,  AEGAUSS )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DE LA MATRICE ELEMENTAIRE DU TETRAEDRE BREZZI FORTIN
C -----    Integrale ( Rho P1B P1B + dt Mhu Grad P1B Grad P1B ) dX
C          P1+BULLE CONTINU POUR UNE COMPOSANTE DE LA VITESSE
C          ELIMINATION DE GAUSS DE LA COMPOSANTE DE LA VITESSE AU BARYCENTRE
C
C ENTREES:	
C --------
C Rho    : DENSITE DE MASSE DU FLUIDE
C DtMhu  : Pas de TEMPS * VISCOSITE DYNAMIQUE DU FLUIDE
C NONOEF : 4 NUMEROS DES SOMMETS DU TETRAEDRE ACTIF
C NBSOM  : NOMBRE DE SOMMETS DE LA TETRAEDRISATION
C XYZSOM : 3 COORDONNEES DES SOMMETS DE LA TETRAEDRISATION
C
C SORTIES:
C --------
C AE     : MATRICE ELEMENTAIRE DE 4-EME LIGNE ELIMINEE PAR GAUSS
C AEGAUSS: 5 COEFFICIENTS DE LA MATRICE ELEMENTAIRE ELIMINES PAR GAUSS
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY Decembre 2010
C MODIFS : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY  Octobre 2012
C23456---------------------------------------------------------------012
      IMPLICIT NONE
      INTEGER           NONOEF(5), NBSOM, I, J, K, L, M
      REAL              XYZSOM(3,NBSOM)
      DOUBLE PRECISION  Rho, DtMhu, AEGAUSS(5)
      DOUBLE PRECISION  AE(15), DELTA, X1, Y1, Z1, A, B, C, D, S
      DOUBLE PRECISION  DF(3,3), TDFDF(3,3), DFM1(3,3), PIVOT
      INTRINSIC         ABS
C
C     DOUBLE PRECISION  TDPDP(5,3,5,3)
C     TDPDP(i,k,j,l) = integrale dPi/dxk dPj/dxl dX
      include"./incl/tdp53dp53.inc"
C
C     CALCUL DE LA MATRICE JACOBIENNE ET DE SON INVERSE
      I  = NONOEF(1)
      X1 = XYZSOM(1,I)
      Y1 = XYZSOM(2,I)
      Z1 = XYZSOM(3,I)
C
      J = NONOEF(2)
      DF(1,1) = XYZSOM(1,J) - X1
      DF(1,2) = XYZSOM(2,J) - Y1
      DF(1,3) = XYZSOM(3,J) - Z1
C
      K = NONOEF(3)
      DF(2,1) = XYZSOM(1,K) - X1
      DF(2,2) = XYZSOM(2,K) - Y1
      DF(2,3) = XYZSOM(3,K) - Z1
C
      L = NONOEF(4)
      DF(3,1) = XYZSOM(1,L) - X1
      DF(3,2) = XYZSOM(2,L) - Y1
      DF(3,3) = XYZSOM(3,L) - Z1
C
C     LE DETERMINANT DE DF
      DELTA = ABS( DF(1,1) * ( DF(2,2) * DF(3,3) - DF(3,2) * DF(2,3) )
     %           + DF(2,1) * ( DF(3,2) * DF(1,3) - DF(1,2) * DF(3,3) )
     %           + DF(3,1) * ( DF(1,2) * DF(2,3) - DF(2,2) * DF(1,3) ) )
C     LE TETRAEDRE EST SUPPOSE DE VOLUME NON NUL
C
C     Integrale Rho P1B P1B dX
C     ------------------------
C     LES 4 COEFFICIENTS GENERIQUES DE LA MATRICE DE MASSE
      S = Rho * DELTA
      A = 29836D0 / 2494800D0 * S
      B =  9046D0 / 2494800D0 * S
      C =   956D0 /  155925D0 * S
      D =  4096D0 /  155925D0 * S
C
CCC     A=  0.11959275292608627D-01
CCC     B=  0.36259419592752926D-02
CCC     C=  0.61311527978194641D-02
CCC     D=  0.26269039602372937D-01
C
C     AFFECTATION DES COEFFICIENTS DE LA MATRICE TRIANGULAIRE INFERIEURE
C     LIGNE PAR LIGNE ET DE LA GAUCHE VERS LA DROITE
      AE(1) = A
C
      AE(2) = B
      AE(3) = A
C
      AE(4) = B
      AE(5) = B
      AE(6) = A
C
      AE(7) = B
      AE(8) = B
      AE(9) = B
      AE(10)= A
C
      AE(11) = C
      AE(12) = C
      AE(13) = C
      AE(14) = C
      AE(15) = D
C
C     Integrale dt Mhu Grad P1B Grad P1B dX
C     -------------------------------------
C     CALCUL DE LA SOUS MATRICE tDF-1 * DF-1 * dt * Mhu / DELTA
      S = DtMhu / DELTA
C
C     LES 9 COEFFICIENTS DE LA MATRICE INVERSE DFM1 SANS / DELTA
      DFM1(1,1) = ( DF(2,2) * DF(3,3) - DF(3,2) * DF(2,3) )
      DFM1(2,1) = ( DF(2,3) * DF(3,1) - DF(3,3) * DF(2,1) )
      DFM1(3,1) = ( DF(2,1) * DF(3,2) - DF(3,1) * DF(2,2) )
C
      DFM1(1,2) = ( DF(1,3) * DF(3,2) - DF(1,2) * DF(3,3) )
      DFM1(2,2) = ( DF(1,1) * DF(3,3) - DF(1,3) * DF(3,1) )
      DFM1(3,2) = ( DF(1,2) * DF(3,1) - DF(1,1) * DF(3,2) )
C
      DFM1(1,3) = ( DF(1,2) * DF(2,3) - DF(2,2) * DF(1,3) )
      DFM1(2,3) = ( DF(1,3) * DF(2,1) - DF(2,3) * DF(1,1) )
      DFM1(3,3) = ( DF(1,1) * DF(2,2) - DF(2,1) * DF(1,2) )
C
C     t[DF]-1 [DF]-1  EST UNE MATRICE SYMETRIQUE
      TDFDF(1,1) = ( DFM1(1,1) * DFM1(1,1)
     %             + DFM1(2,1) * DFM1(2,1)
     %             + DFM1(3,1) * DFM1(3,1) ) * S
C
      TDFDF(2,1) = ( DFM1(1,2) * DFM1(1,1)
     %             + DFM1(2,2) * DFM1(2,1)
     %             + DFM1(3,2) * DFM1(3,1) ) * S
      TDFDF(1,2) = TDFDF(2,1)
C
      TDFDF(3,1) = ( DFM1(1,3) * DFM1(1,1)
     %             + DFM1(2,3) * DFM1(2,1)
     %             + DFM1(3,3) * DFM1(3,1) ) * S
      TDFDF(1,3) = TDFDF(3,1)
C
      TDFDF(2,2) = ( DFM1(1,2) * DFM1(1,2)
     %             + DFM1(2,2) * DFM1(2,2)
     %             + DFM1(3,2) * DFM1(3,2) ) * S
C
      TDFDF(3,2) = ( DFM1(1,3) * DFM1(1,2)
     %             + DFM1(2,3) * DFM1(2,2)
     %             + DFM1(3,3) * DFM1(3,2) ) * S
      TDFDF(2,3) = TDFDF(3,2)
C
      TDFDF(3,3) = ( DFM1(1,3) * DFM1(1,3)
     %             + DFM1(2,3) * DFM1(2,3)
     %             + DFM1(3,3) * DFM1(3,3) ) * S
C
      M = 0
      DO J=1,5
         DO I=1,J
            M = M + 1
            DO K=1,3
               DO L=1,3
                  AE(M) = AE(M) + TDPDP(I,K,J,L) * TDFDF(K,L)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
C     ELIMINATION DE GAUSS DU DL 5 DE LA VITESSE P1B AU BARYCENTRE
C     LA LIGNE 5 DE AE EST ELIMINEE PAR GAUSS ET STOCKEE DANS AEGAUSS
C     POUR PERMETTRE LE CALCUL DES VITESSES APRES RESOLUTION
C     ---------------------------------------------------------------
      DO J=1,5
         AEGAUSS(J) = AE(10+J)
      ENDDO
C
C     TRIANGULATION DE GAUSS DE LA LIGNE 5 AVEC PIVOT=AE(5,5)
      PIVOT = AE(15)
      M = 0
      DO I = 1, 4
         S = AE(10+I) / PIVOT
         DO J = 1, I
C
C           COEFFICIENT AE(I,J)=AE(M)
            M = M + 1
            AE(M) = AE(M) - AE(10+J) * S
C
         ENDDO
      ENDDO
C
      RETURN
      END
