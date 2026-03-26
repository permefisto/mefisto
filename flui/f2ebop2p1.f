      SUBROUTINE F2EB0P2P1( NBNOVI, XYZNOE, NONOTR, NONOSO,
     %                      NOOBLA, NUMILI, NUMALI, LTDELI,
     %                      NOOBSF, NUMISU, NUMASU, LTDESU,
     %                      DELTAT, CoGrPr, TP2P2,  NBSOM,  PRESS0,
     %                      VE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DU SECOND MEMBRE DU TRIANGLE TAYLOR HOOD
C -----    P2 CONTINU POUR LES 2 COMPOSANTES DE LA VITESSE
C          P1 CONTINU POUR LA PRESSION
C          INTERPOLATION DE FOMEGA AUX 3 SOMMETS + 3 MILIEUX DES ARETES
C          INTERPOLATION DE FGAMMA AUX 2 SOMMETS + MILIEU DE L'ARETE
C          Integrale  dt CoGrPr div V  P(tn) + dt  V FOmega(tn) dX
C ENTREES:
C --------
C NBNOVI : NOMBRE DE NOEUDS DE LA TRIANGULATION
C XYZNOE : 3 COORDONNEES DES NBNOVI NOEUDS DU TRIANGLE
C NONOTR : NONOTR(I) NO GLOBAL DU I-EME NOEUD DU TRIANGLE I=1,...,6
C NONOSO : NONOSO(I) = NUMERO DE SOMMET DE 1 A NBSOM DU I-EME NOEUD GLOBAL
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
C TP2P2  : Integrale tP2 P2 dX sur LE TRIANGLE UNITE
C
C NTDLPR : NOMBRE DE DL EN PRESSION = NBSOMMETS
C NBSOM  : NOMBRE DE SOMMETS DE LA TETRAEDRISATION
C          = NOMBRE DE DL DE LA PRESSION
C PRESS0 : DL DE LA PRESSION AUX SOMMETS DE LA TETRAEDRISATION
C
C SORTIES:
C --------
C VE     : DL ELEMENTAIRES DE LA VITESSE A L'INSTANT temps+deltat
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC & St Pierre du Perray     Avril 2011
C23456---------------------------------------------------------------012
!$    USE OMP_LIB
      IMPLICIT   NONE
      include"./incl/langue.inc"
      include"./incl/donflu.inc"
      include"./incl/ctemps.inc"
      INTEGER         LECTEU, IMPRIM, NUNITE
      COMMON /UNITES/ LECTEU, IMPRIM, NUNITE(30)
C
      REAL              XYZNOE(3,NBNOVI)
      DOUBLE PRECISION  PRESS0(NBSOM), TP2P2(6,6)
      INTEGER           NOOBLA(3), NUMILI, NUMALI
      INTEGER           LTDELI(1:MXDOFL,NUMILI:NUMALI)
      INTEGER           NOOBSF(1), NUMISU, NUMASU
      INTEGER           LTDESU(1:MXDOFL,NUMISU:NUMASU)
      INTEGER           NBSOM, NBNOVI
      INTEGER           NONOTR(6), NONOSO(*), I, J, K, KK, M, NL
      DOUBLE PRECISION  CoGrPr, DELTAT
      DOUBLE PRECISION  VE(6,2), FORCE(3,6)
      DOUBLE PRECISION  DELTAe, X21, X31, Y21, Y31, VN(2)
      DOUBLE PRECISION  X1, Y1, D, S
      INTRINSIC         ABS
C
C     CONTRIBUTION Integrale t(div Vitesse P2) CoefP Pression dX
C     DELTAT CogrPr integrale t( DF-1 ligne1  DP2
C                               +DF-1 ligne2  DP2 ) P1 dX {dl Pression}
C     -----------------------------------------------------------------
      CALL F2EX4P2P1( DELTAT*CoGrPr, NBNOVI, XYZNOE, NONOTR, NONOSO,
     %                NBSOM, PRESS0, VE )
C
C     CONTRIBUTION DES EFFORTS SURFACIQUES INTERPOLES TAYLOR-HOOD
C     -----------------------------------------------------------
      IF( LTDESU(LPFORC,NOOBSF(1)) .GT. 0 ) THEN
C
C        VALEUR DES EFFORTS SURFACIQUES AUX 6 NOEUDS DU TRIANGLE
         DO J=1,6
            M = NONOTR(J)
            X1 = XYZNOE(1,M)
            Y1 = XYZNOE(2,M)
            CALL REFORC( 3, NOOBSF(1), 2, X1,Y1,0D0, 0D0,0D0,0D0,
     %                   LTDESU(LPFORC,NOOBSF(1)), FORCE(1,J) )
         ENDDO
C
C        COEFFICIENTS DU SECOND MEMBRE LIES AUX EFFORTS SURFACIQUES
C        CALCUL DE LA MATRICE JACOBIENNE ET DE SON INVERSE
         M  = NONOTR(1)
         X1 = XYZNOE(1,M)
         X21 = XYZNOE(1,NONOTR(2)) - X1
         X31 = XYZNOE(1,NONOTR(3)) - X1
C
         Y1 = XYZNOE(2,M)
         Y21 = XYZNOE(2,NONOTR(2)) - Y1
         Y31 = XYZNOE(2,NONOTR(3)) - Y1
C
C        CALCUL DU DETERMINANT DE LA JACOBIENNE DFe
         DELTAe = ABS( X21*Y31 - X31*Y21 )
         IF( DELTAe .LE. 0D0 ) THEN
            IF( LANGAG .EQ. 0 ) THEN
               WRITE(IMPRIM,*) 'F2EX2P2P1: ATTENTION EF',
     %                     ' de SURFACE*2=',DELTAe,' NON PRIS EN COMPTE'
            ELSE
               WRITE(IMPRIM,*) 'F2EX2P2P1: ATTENTION FE',
     %                       ' of SURFACE*2=',DELTAe,' is NOT COMPUTED'
            ENDIF
            GOTO 90
         ENDIF
C
C        COEFFICIENTS DU SECOND MEMBRE LIES AUX EFFORTS SURFACIQUES
C        dt Integrale tP2 P2 dX {Force(t,ne)}
         D = DELTAT * DELTAe
         DO I=1,6
            DO J=1,6
               S = TP2P2(I,J) * D
               VE(I,1) = VE(I,1) + S * FORCE(1,J)
               VE(I,2) = VE(I,2) + S * FORCE(2,J)
            ENDDO
         ENDDO
C
      ENDIF
C
C     CONTRIBUTION DES EFFORTS LINEIQUES SUR LES ARETES DU TRIANGLE
C     -------------------------------------------------------------
 90   DO K=1,3
C        NO DE LA LIGNE DE L'ARETE K
         NL = NOOBLA(K)
         IF( NL .GT. 0 ) THEN
            IF( LTDELI(LPFORC,NL) .GT. 0 ) THEN
C
C              IL EXISTE UNE FORCE DONNEE SUR CETTE LIGNE
C              INTERPOLEE AUX 2 SOMMETS ET LE MILIEU DE L'ARETE
C              LE NUMERO DES 2 SOMMETS DU COTE K DU TRIANGLE
               IF( K .NE. 3 ) THEN
                  KK = K+1
               ELSE
                  KK = 1
               ENDIF
C
C              LE VECTEUR ORTHOGONAL AU VECTEUR TANGENT DE L'ARETE K
               VN(1) = XYZNOE(2,NONOTR(KK)) - XYZNOE(2,NONOTR(K))
               VN(2) = XYZNOE(1,NONOTR(K)) - XYZNOE(1,NONOTR(KK))
C
C              LE JACOBIEN
               D = SQRT( VN(1)**2 + VN(2)**2 )
C
C              LE VECTEUR NORMAL UNITAIRE
               VN(1) = VN(1) / D
               VN(2) = VN(2) / D
C
C              CALCUL DES EFFORTS LIES AU SOMMET K
               X1 = XYZNOE(1,NONOTR(K))
               Y1 = XYZNOE(2,NONOTR(K))
               CALL REFORC( 2, NL, 2,
     %                      X1,Y1,0D0, VN(1),VN(2),0D0,
     %                      LTDELI(LPFORC,NL), FORCE(1,1) )
C
C              CALCUL DES EFFORTS LIES AU SOMMET KK
               X1 = XYZNOE(1,NONOTR(KK))
               Y1 = XYZNOE(2,NONOTR(KK))
               CALL REFORC( 2, NL, 2,
     %                      X1,Y1,0D0, VN(1),VN(2),0D0,
     %                      LTDELI(LPFORC,NL), FORCE(1,2) )
C
C              CALCUL DES EFFORTS LIES AU MILIEU DE L'ARETE
               X1 = ( X1 + XYZNOE(1,NONOTR(K)) ) / 2D0
               Y1 = ( Y1 + XYZNOE(2,NONOTR(K)) ) / 2D0
               CALL REFORC( 2, NL, 2,
     %                      X1,Y1,0D0, VN(1),VN(2),0D0,
     %                      LTDELI(LPFORC,NL), FORCE(1,3) )
C
C              LE PRODUIT DU JACOBIEN PAR DELTAT PAR POIDS
               D = D * DELTAT / 30D0
C
C              PRISE EN COMPTE DES COEFF Integrale( tP2 P2 ) Force
               DO M=1,2
C
C                 LA COMPOSANTE M DE LA FORCE ELEMENTAIRE AU SOMMET K
C                 COEFFICIENT DU SECOND MEMBRE ELEMENTAIRE DE LA COMPOSANTE M
                  VE(K,M) = VE(K,M) + D *
     %                        (4D0*FORCE(M,1)-FORCE(M,2)+2D0*FORCE(M,3))
C
C                 LA COMPOSANTE M DE LA FORCE ELEMENTAIRE AU SOMMET KK
C                 COEFFICIENT DU SECOND MEMBRE ELEMENTAIRE DE LA COMPOSANTE M
                  VE(KK,M) = VE(KK,M) + D *
     %                       (-FORCE(M,1)+4D0*FORCE(M,2)+2D0*FORCE(M,3))
C
C                 LA COMPOSANTE M DE LA FORCE ELEMENTAIRE AU MILIEU DE l'ARETE
C                 COEFFICIENT DU SECOND MEMBRE ELEMENTAIRE DE LA COMPOSANTE M
                  VE(K+3,M) = VE(K+3,M) + D *
     %                      2D0 * (FORCE(M,1)+FORCE(M,2)+8D0*FORCE(M,3))
               ENDDO
C
            ENDIF
         ENDIF
C
      ENDDO
C
      RETURN
      END
