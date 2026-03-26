      SUBROUTINE F2EX2P1BP1( XYZSOM, NUSOTR,
     %                       NOOBLA, NUMILI, NUMALI, LTDELI,
     %                       NOOBSF, NUMISU, NUMASU, LTDESU,
     %                       DELTAT, CoGrPr, NBSOM,  PRESS0,
     %                       VE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DU SECOND MEMBRE DU TRIANGLE BREZZI FORTIN
C -----    P1+BULLE CONTINU POUR LES 2 COMPOSANTES DE LA VITESSE
C          P1       CONTINU POUR LA PRESSION
C          INTERPOLATION DE FOMEGA AUX 3 SOMMETS+BARYCENTRE DU TRIANGLE
C          INTERPOLATION DE FGAMMA AUX 2 SOMMETS DE L'ARETE
C          Integrale  dt CoGrPr div V  P(tn) + dt  V FOmega(tn) dX
C ENTREES:
C --------
C XYZSOM : 3 COORDONNEES DES 3 SOMMETS DE LA TRIANGULATION
C          ICI LA COMPOSANTE Z=XYZSOM(3,.) N'EST PAS UTILISEE
C NUSOTR : NUSOTR(I) NO GLOBAL DU I-EME SOMMET DU TRIANGLE I=1,...,3
C          NUSOTR(4) = NBSOMMET + NUMERO DU TRIANGLE = NO BARYCENTRE
C
C NOOBLA : NUMERO DES OBJETS LIGNES DES ARETES DE L'ELEMENT FINI
C NUMILI : NUMERO MINIMAL DES OBJETS LIGNES
C NUMALI : NUMERO MAXIMAL DES OBJETS LIGNES
C LTDELI : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES DU FLUIDE
C          DES OBJETS LIGNES
C
C NOOBSF : NUMERO D'OBJET SURFACE DE L'ELEMENT FINI
C NUMISU : NUMERO MINIMAL DES OBJETS SURFACES
C NUMASU : NUMERO MAXIMAL DES OBJETS SURFACES
C LTDESU : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES DU FLUIDE
C          DES OBJETS SURFACES
C
C DELTAT : PAS DE TEMPS POUR L'INTEGRATION EN TEMPS
C CoGrPr : COEFFICIENT DU GRADIENT DE PRESSION DANS LES EDP
C
C NTDLVI : NOMBRE DE DL EN VITESSES = 2 * (NBSOMMETS+NBTRIANGLES)
C NTDLPR : NOMBRE DE DL EN PRESSION = NBSOMMETS
C NBSOM  : NOMBRE DE SOMMETS DE LA TETRAEDRISATION
C          = NOMBRE DE DL DE LA PRESSION
C PRESS0 : DL DE LA PRESSION AUX SOMMETS DE LA TETRAEDRISATION
C
C SORTIES:
C --------
C VE     : DL ELEMENTAIRES DE LA VITESSE A L'INSTANT temps+deltat
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC & St Pierre du Perray  Decembre 2010
C23456---------------------------------------------------------------012
      IMPLICIT   NONE
      include"./incl/donflu.inc"
      include"./incl/ctemps.inc"
      INTEGER         LECTEU, IMPRIM, NUNITE
      COMMON /UNITES/ LECTEU, IMPRIM, NUNITE(30)
C
      REAL              XYZSOM(3,NBSOM)
      DOUBLE PRECISION  PRESS0(NBSOM)
      INTEGER           NOOBLA(3), NUMILI, NUMALI
      INTEGER           LTDELI(1:MXDOFL,NUMILI:NUMALI)
      INTEGER           NOOBSF(1), NUMISU, NUMASU
      INTEGER           LTDESU(1:MXDOFL,NUMISU:NUMASU)
      INTEGER           NBSOM
      INTEGER           NUSOTR(4), I, J, K, KK, M, NL
      DOUBLE PRECISION  CoGrPr, DELTAT
      DOUBLE PRECISION  VE(4,2), FORCE(3,4)
      DOUBLE PRECISION  TP1BP1(4,4)
      DOUBLE PRECISION  DELTA, X21, X31, Y21, Y31, VN(2)
      DOUBLE PRECISION  X1, Y1, A, B, C, D, S
      INTRINSIC         ABS
C
C     L'INTEGRALE  + dt CoGrPr  t(Div Vitesse) Pression dx
C     DELTAT CogrPr integrale t( DF-1 ligne1  DP1B
C                               +DF-1 ligne2  DP1B ) P1 dX {dl Pression}
C     ------------------------------------------------------------------
C     LE TRIANGLE EST SUPPOSE DE SURFACE NON NULLE
C     LES 4 COEFFICIENTS DE LA MATRICE INVERSE DFM1 SANS / DELTA
C     CAR COMPENSE PAR LA MULTIPLICATION PAR DELTA DE L'INTEGRATION
      CALL F2EX4P1BP1( XYZSOM, NUSOTR, DELTAT, CoGrPr,
     %                 NBSOM,  PRESS0, VE )
C
C     CONTRIBUTION DES EFFORTS SURFACIQUES INTERPOLES BREZZI-FORTIN
C     -------------------------------------------------------------
      IF( LTDESU(LPFORC,NOOBSF(1)) .GT. 0 ) THEN
C
C        VALEUR DES EFFORTS SURFACIQUES AUX 3 SOMMETS DU TRIANGLE
         DO J=1,3
            K = NUSOTR(J)
            X1 = XYZSOM(1,K)
            Y1 = XYZSOM(2,K)
            CALL REFORC( 3, NOOBSF(1), 2, X1,Y1,0D0, 0D0,0D0,0D0,
     %                   LTDESU(LPFORC,NOOBSF(1)), FORCE(1,J) )
         ENDDO
C
C        VALEUR DES EFFORTS SURFACIQUES AU BARYCENTRE
         I = NUSOTR(1)
         J = NUSOTR(2)
         K = NUSOTR(3)

         X1 = (XYZSOM(1,I) + XYZSOM(1,J) + XYZSOM(1,K) ) / 3D0
         Y1 = (XYZSOM(2,I) + XYZSOM(2,J) + XYZSOM(2,K) ) / 3D0
         CALL REFORC( 3, NOOBSF(1), 2, X1,Y1,0D0, 0D0,0D0,0D0,
     %                LTDESU(LPFORC,NOOBSF(1)), FORCE(1,4) )
C
C        CALCUL DE LA MATRICE JACOBIENNE
         X21 = XYZSOM(1,J) - XYZSOM(1,I)
         X31 = XYZSOM(1,K) - XYZSOM(1,I)
C
         Y21 = XYZSOM(2,J) - XYZSOM(2,I)
         Y31 = XYZSOM(2,K) - XYZSOM(2,I)
C
C        CALCUL DU DETERMINANT DE LA JACOBIENNE DFe
         DELTA = ABS( X21*Y31 - X31*Y21 )
C
C        LES 4 COEFFICIENTS GENERIQUES DE LA MATRICE DE MASSE tP1B P1B
         S = DELTAT * DELTA
         A = 83D0 / 1680D0 * S
         B = 13D0 / 1680D0 * S
         C =  3D0 /  112D0 * S
         D = 81D0 /  560D0 * S
C
         TP1BP1(1,1) = A
         TP1BP1(2,1) = B
         TP1BP1(3,1) = B
         TP1BP1(4,1) = C
C
         TP1BP1(1,2) = B
         TP1BP1(2,2) = A
         TP1BP1(3,2) = B
         TP1BP1(4,2) = C
C
         TP1BP1(1,3) = B
         TP1BP1(2,3) = B
         TP1BP1(3,3) = A
         TP1BP1(4,3) = C
C
         TP1BP1(1,4) = C
         TP1BP1(2,4) = C
         TP1BP1(3,4) = C
         TP1BP1(4,4) = D
C
C        COEFFICIENTS DU SECOND MEMBRE LIES AUX EFFORTS SURFACIQUES
         DO I=1,4
            DO J=1,4
               D = TP1BP1(I,J)
               VE(I,1) = VE(I,1) + D * FORCE(1,J)
               VE(I,2) = VE(I,2) + D * FORCE(2,J)
            ENDDO
         ENDDO
C
      ENDIF
C
C     CONTRIBUTION DES EFFORTS  LINEIQUES SUR LES ARETES DU TRIANGLE
C     --------------------------------------------------------------
      DO K=1,3
C
C        NO DE LA LIGNE DE L'ARETE K
         NL = NOOBLA(K)
         IF( NL .GT. 0 ) THEN
            IF( LTDELI(LPFORC,NL) .GT. 0 ) THEN
C
C              IL EXISTE UNE FORCE DONNEE SUR CETTE LIGNE
C              LE NUMERO DES SOMMETS DU COTE K DU TRIANGLE
               IF( K .NE. 3 ) THEN
                  KK = K+1
               ELSE
                  KK = 1
               ENDIF
C
C              LE VECTEUR ORTHOGONAL AU VECTEUR TANGENT
               I = NUSOTR(K)
               J = NUSOTR(KK)
               VN(1) = XYZSOM(2,J) - XYZSOM(2,I)
               VN(2) = XYZSOM(1,I) - XYZSOM(1,J)
C
C              LE JACOBIEN
               D = SQRT( VN(1)**2 + VN(2)**2 )
C
C              LE VECTEUR NORMAL UNITAIRE
               VN(1) = VN(1) / D
               VN(2) = VN(2) / D
C
C              LE PRODUIT DU JACOBIEN PAR LE POIDS et DELTAT
               D = D / 6D0 * DELTAT
C
C              CALCUL DES EFFORTS LIES AU SOMMET K
               X1 = XYZSOM(1,I)
               Y1 = XYZSOM(2,I)
               CALL REFORC( 2, NL, 2,
     %                      X1,Y1,0D0, VN(1),VN(2),0D0,
     %                      LTDELI(LPFORC,NL), FORCE(1,1) )
C
C              CALCUL DES EFFORTS LIES AU SOMMET KK
               X1 = XYZSOM(1,J)
               Y1 = XYZSOM(2,J)
               CALL REFORC( 2, NL, 2,
     %                      X1,Y1,0D0, VN(1),VN(2),0D0,
     %                      LTDELI(LPFORC,NL), FORCE(1,2) )
C
C              PRISE EN COMPTE DES COEFF Integrale( tLa La ) Force
               DO M=1,2
C
C                 LA COMPOSANTE M DE LA FORCE ELEMENTAIRE AU SOMMET K
C                 COEFFICIENT DU SECOND MEMBRE ELEMENTAIRE DE LA COMPOSANTE M
                  VE(K,M) = VE(K,M) + D * (2D0*FORCE(M,1)+FORCE(M,2))
C
C                 LA COMPOSANTE M DE LA FORCE ELEMENTAIRE AU SOMMET KK
C                 COEFFICIENT DU SECOND MEMBRE ELEMENTAIRE DE LA COMPOSANTE M
                  VE(KK,M) = VE(KK,M) + D * (FORCE(M,1)+2D0*FORCE(M,2))
C
               ENDDO
C
            ENDIF
         ENDIF
C
      ENDDO
C
      RETURN
      END
