      SUBROUTINE F2EX4P2P1GRAD( XYZNOE, NONOTR, NONOSO, DELTAT, CoGrPr,
     %                          NBNOVI, NBSOM,  PRESP1,  BE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  CALCUL DU SECOND MEMBRE DU TRIANGLE TAYLOR YOUNG
C -----  BEk = -DELTAT CoGrPr Integrale sur ef  tV(x) gradk PRESP1(x) dx
C        P2 CONTINU POUR LES 2 COMPOSANTES DE LA VITESSE k=1,2
C        P1 CONTINU POUR LA PRESSION

C ENTREES:
C --------
C XYZNOE : XYZNOE(3,NBNOVI)  3 COORDONNEES DES NOEUDS DU MAILLAGE P2
C          ICI LA COMPOSANTE Z=XYZNOE(3,.) N'EST PAS UTILISEE
C NONOTR : NONOTR(I) NO GLOBAL DU I-EME NOEUD DU TRIANGLE P2 I=1,...,6
C NONOSO : NONOSO(I) = NUMERO DE SOMMET DE 1 A NBSOM DU I-EME NOEUD GLOBAL

C DELTAT : PAS DE TEMPS
C CoGrPr : COEFFICIENT DU GRADIENT DE PRESSION DANS LES EDP
C NBNOVI : NOMBRE DE NOEUDS VITESSE P2
C NBSOM  : NOMBRE DE SOMMETS DE LA TRIANGULATION P2
C          = NOMBRE DE DL DE LA PRESSION
C PRESP1 : DL DE LA PRESSION P1 AUX SOMMETS DE LA TRIANGULATION P2

C SORTIES:
C --------
C BE     : SECOND MEMBRE ELEMENTAIRE DE LA VITESSE A L'INSTANT temps+deltat
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC & St Pierre du Perray   Octobre 2022
C23456---------------------------------------------------------------012
      IMPLICIT NONE
      REAL              XYZNOE(3,NBNOVI)
      INTEGER           NBNOVI, NBSOM, NONOTR(6), NONOSO(*),
     %                  I, J, K, NSJ
      DOUBLE PRECISION  PRESP1(NBSOM), DELTAT, CoGrPr, BE(6,2)
      DOUBLE PRECISION  DFM1DLa(2,3), S, COEF

Cccc     DOUBLE PRECISION  INTP2(6)
Cccc     INTP2(i) = Integrale sur e chapeau P2 dx dy
Cccc     INTP2(1:3)=0, INTP2(4:6)=1/6
Cccc     DATA  INTP2 / 0, 0, 0, 1/6, 1/6, 1/6 /  INTEGRE AU CALCUL

C     DFM1 DLambda(2,3)
C     -----------------
      I = NONOTR(1)
      J = NONOTR(2)
      K = NONOTR(3)

      DFM1DLa(1,1) = XYZNOE(2,J) - XYZNOE(2,K)
      DFM1DLa(2,1) = XYZNOE(1,K) - XYZNOE(1,J)

      DFM1DLa(1,2) = XYZNOE(2,K) - XYZNOE(2,I)
      DFM1DLa(2,2) = XYZNOE(1,I) - XYZNOE(1,K)

      DFM1DLa(1,3) = XYZNOE(2,I) - XYZNOE(2,J)
      DFM1DLa(2,3) = XYZNOE(1,J) - XYZNOE(1,I)

C     - DELTAT * CoGrPr INTEGRALE  tVitesse  Grad Pression dX
C     -------------------------------------------------------
      COEF = - DELTAT * CoGrPr / 6D0

C     K EST LA COMPOSANTE DE LA VITESSE
      DO K=1,2

C        L'INTEGRALE Lambda i ( 2 Lambda i - 1 ) AUX 3 SOMMETS EST NULLE
         DO I=1,3
C           LE COEFFICIENT BE(I,K) EST NUL
            BE(I,K) = 0D0
         ENDDO

C        - DELTAT * CoGrPr Integrale tP2 dX  (=1/6 AUX NOEUDS MILIEUX DES COTES)
         S = 0D0
         DO J=1,3
C           NONOTR(J) NO GLOBAL DU MILIEU DU COTE J DU TRIANGLE P2
C           NSJ       NO GLOBAL DU SOMMET J DU TRIANGLE P1
            NSJ = NONOSO( NONOTR(J) )
            S = S + DFM1DLa(K,J) * PRESP1( NSJ )
         ENDDO

         DO I=4,6
C           LE COEFFICIENT BE(I,K) EST INITIALISE
            BE(I,K) = S * COEF
         ENDDO

      ENDDO

      RETURN
      END
