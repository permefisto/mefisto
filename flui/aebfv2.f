      SUBROUTINE AEBFV2( Rho, DtMhu, NONOEF, NBSOM, XYZSOM,
     %                   AE,  AEGAUSS )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DE LA MATRICE ELEMENTAIRE DU TRIANGLE BREZZI FORTIN
C -----    Integrale ( Rho P1B P1B + dt Mhu Grad P1B Grad P1B ) dX
C          P1+BULLE CONTINU POUR UNE COMPOSANTE DE LA VITESSE
C          ELIMINATION DE GAUSS DE LA COMPOSANTE DE LA VITESSE AU BARYCENTRE
C
C ENTREES:	
C --------
C Rho    : DENSITE DE MASSE DU FLUIDE
C DtMhu  : Pas de TEMPS * VISCOSITE DYNAMIQUE DU FLUIDE
C NONOEF : 3 NUMEROS DES SOMMETS DU TRIANGLE ACTIF
C NBSOM  : NOMBRE DE SOMMETS DE LA TRIANGULATION
C XYZSOM : 3 COORDONNEES DES SOMMETS DE LA TRIANGULATION
C
C SORTIES:
C --------
C AE     : MATRICE ELEMENTAIRE DE 4-EME LIGNE ELIMINEE PAR GAUSS
C AEGAUSS: 4 COEFFICIENTS DE LA MATRICE ELEMENTAIRE ELIMINES PAR GAUSS
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY Decembre 2010
C MODIFS : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY  Octobre 2012
C23456---------------------------------------------------------------012
      IMPLICIT NONE
      INTEGER           NONOEF(4), NBSOM, I, J, K, L, M
      REAL              XYZSOM(3,NBSOM)
      DOUBLE PRECISION  Rho, DtMhu, AEGAUSS(4)
      DOUBLE PRECISION  AE(10), DELTA, X21, Y21, X31, Y31, A, B, C, D, S
      DOUBLE PRECISION  TDFDF(2,2), PIVOT
      INTRINSIC         ABS
C
      DOUBLE PRECISION  TDP1BDP1B(4,2,4,2)
C     TDP1BDP1B(i,k,j,l) = integrale dPi/dxk DPj/dxl dX
      DATA TDP1BDP1B/
     %  0.95D+00,   -0.5D-01,
     %  0.45D+00,   -0.135D+01,
     %  0.725D+00,   0.225D+00,
     % -0.275D+00,  -0.675D+00,
     % -0.5D-01,     0.95D+00,
     %  0.45D+00,   -0.135D+01,
     % -0.275D+00,   0.225D+00,
     %  0.725D+00,  -0.675D+00,
     %  0.45D+00,    0.45D+00,
     %  0.45D+00,   -0.135D+01,
     %  0.225D+00,   0.225D+00,
     %  0.225D+00,  -0.675D+00,
     % -0.135D+01,  -0.135D+01,
     % -0.135D+01,   0.405D+01,
     % -0.675D+00,  -0.675D+00,
     % -0.675D+00,   0.2025D+01,
     %  0.725D+00,  -0.275D+00,
     %  0.225D+00,  -0.675D+00,
     %  0.95D+00,    0.45D+00,
     % -0.5D-01,    -0.135D+01,
     %  0.225D+00,   0.225D+00,
     %  0.225D+00,  -0.675D+00,
     %  0.45D+00,    0.45D+00,
     %  0.45D+00,   -0.135D+01,
     % -0.275D+00,   0.725D+00,
     %  0.225D+00,  -0.675D+00,
     % -0.5D-01,     0.45D+00,
     %  0.95D+00,   -0.135D+01,
     % -0.675D+00,  -0.675D+00,
     % -0.675D+00,   0.2025D+01,
     % -0.135D+01,  -0.135D+01,
     % -0.135D+01,   0.405D+01 /
C
C     CALCULS INITIAUX SUR Fe: e ref -> e
      I = NONOEF(1)
      J = NONOEF(2)
      K = NONOEF(3)
C
      A = XYZSOM(1,I)
      X21 = XYZSOM(1,J) - A
      X31 = XYZSOM(1,K) - A
C
      B = XYZSOM(2,I)
      Y21 = XYZSOM(2,J) - B
      Y31 = XYZSOM(2,K) - B
C
C     CALCUL DU DETERMINANT DE LA JACOBIENNE DFe
      DELTA = ABS( X21*Y31 - X31*Y21 )
C
C     Integrale Rho P1B P1B dX
C     ------------------------
C     LES 4 COEFFICIENTS GENERIQUES DE LA MATRICE DE MASSE
      S = Rho * DELTA
      A = S * 83D0 / 1680D0
      B = S * 13D0 / 1680D0
      C = S *  3D0 /  112D0
      D = S * 81D0 /  560D0
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
      AE(7)  = C
      AE(8)  = C
      AE(9)  = C
      AE(10) = D
C
C     Integrale dt Mhu Grad P1B Grad P1B dX
C     -------------------------------------
C     CALCUL DE LA SOUS MATRICE dt * tDF-1 * DF-1 * Mhu / DELTA
      S = DtMhu / DELTA
      TDFDF(1,1) = S * (  X31 * X31 + Y31 * Y31 )
      TDFDF(1,2) = S * ( -X21 * X31 - Y21 * Y31 )
      TDFDF(2,1) = TDFDF(1,2)
      TDFDF(2,2) = S * (  X21 * X21 + Y21 * Y21 )
C
      M = 0
      DO J=1,4
         DO I=1,J
            M = M + 1
            DO K=1,2
               DO L=1,2
                  AE(M) = AE(M) + TDP1BDP1B(I,K,J,L) * TDFDF(K,L)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
C     ELIMINATION DE GAUSS DU DL 4 DE LA VITESSE P1B AU BARYCENTRE
C     LA LIGNE 4 DE AE EST ELIMINEE PAR GAUSS ET STOCKEE DANS AEGAUSS
C     POUR PERMETTRE LE CALCUL DES VITESSES APRES RESOLUTION
C     ---------------------------------------------------------------
      DO I=1,4
         AEGAUSS(I) = AE(6+I)
      ENDDO
C
C     TRIANGULATION DE GAUSS DE LA LIGNE 4 SUR AE
      PIVOT = AE(10)
      M = 0
      DO I = 1, 3
         S = AE(6+I) / PIVOT
         DO J = 1, I
C           COEFFICIENT AE(I,J)=AE(M)
            M = M + 1
            AE(M) = AE(M) - AE(6+J) * S
         ENDDO
      ENDDO
C
      RETURN
      END
