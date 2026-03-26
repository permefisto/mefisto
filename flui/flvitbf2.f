      SUBROUTINE FLVITBF2( NBPOE,  X,      SURFTR, NONOEF,
     %                     NTDLVP, NBVECT, VXYPR,  NDDLNV,
     %                     INTVIT, NOFOVI, TIMES,  INTVER )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DE L'INTEGRALE L2 DES NBVECT VITESSE
C -----    SUR UN TRIANGLE DE BREZZI-FORTIN DE VITESSE AVEC 2 COMPOSANTES
C          POLYNOME DE LAGRANGE DE DEGRE 1 + BULLE AU BARYCENTRE
C
C ENTREES:
C --------
C NBPOE  : NOMBRE DE NOEUDS=POINTS DE L'EF
C X      : 2 COORDONNEES DES NBPOE NOEUDS DU TRIANGLE
C SURFTR : SURFACE DU TRIANGLE
C NONOEF : NUMEROS DES 6 NOEUDS DU TRIANGLE COURANT
C NTDLVP : NOMBRE TOTAL DE DEGRES DE LIBERTE DU MAILLAGE DE L'OBJET
C NBVECT : NOMBRE DE VECTEURS SOLUTION VITESSEPRESSION
C VXYPR  : NBVECT VECTEURS VITESSEPRESSION SUR LE MAILLAGE
C NDDLNV : TABLEAU DES POINTEURS SUR LE DERNIER DL POUR CHAQUE NOEUD VITESSE
C          NDDLNV(I)= NO DERNIER DL DU NOEUD (SOMMET ou MILIEU ou BARYCENTRE)
C          CE TABLEAU EST DIMENSIONNE A 0:NBNOVI
C NOFOVI : NUMERO DE LA FONCTION UTILISATEUR VITESSE_EXACTE(t,x,y,z,nocomp)
C          0 SI ELLE N'EXISTE PAS => INTVER N'EST PAS CALCULE
C
C SORTIE :
C --------                                                      (Vke1)
C INTVIT : Som  Jacobien Som (Vke1,...,Vke4) [int P1Bi P1Bj DX] (  . )
C          e de E        k=1 ... 4                              (Vke4)
C                                                                    (Vk1E-Vk1C)
C INTVER : Som  Som Jacobien (Vk1E-Vk1C,...,Vk4E-Vk4C) [int Pi Pj dX](   ...   )
C       e de E  k=1,...,2                                            (Vk4E-Vk4C)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET TEXAS A & M UNIVERSITY at QATAR   Janvier 2012
C23456---------------------------------------------------------------012
C     LES 4 COEFFICIENTS GENERIQUES DE LA MATRICE Integrale tP1B P1B dX
      DOUBLE PRECISION  A, B, C, D
      PARAMETER        (A=83D0/1680D0, B=13D0/1680D0,
     %                  C= 3D0/ 112D0, D=81D0/ 560D0 )
C
      INTEGER           NONOEF(NBPOE), NDDLNV(0:*)
      REAL              X(NBPOE,2)
C           NBPOE=3? et pour X?
      DOUBLE PRECISION  VXYPR(NTDLVP,NBVECT),
     %                  INTVIT(NBVECT), INTVER(NBVECT),
     %                  DELTA, SURFTR, INTVEF, INTVERR, PARAMF(5),
     %                  VitExkI, VitExkJ, XBAR, YBAR
      REAL              TIMES(NBVECT)
C
C     INTEGRALES P1B P1B SUR LE TRIANGLE P1B de REFERENCE
      DOUBLE PRECISION  P1BP1B2D(4,4)
C                       P1BP1B2D(i,j) = integrale P1Bi P1Bj dX
      DATA   P1BP1B2D / A, B, B, C,
     %                  B, A, B, C,
     %                  B, B, A, C,
     %                  C, C, C, D /
C
C     JACOBIEN = 2 x (SURFACE DU TRIANGLE)
      DELTA = 2D0 * SURFTR
C
C     COORDONNEES DU BARYCENTRE
      IF( NOFOVI .GT. 0 ) THEN
         XBAR = ( X(1,1) + X(2,1) + X(3,1) ) / 3D0
         YBAR = ( X(1,2) + X(2,2) + X(3,2) ) / 3D0
         PARAMF(4) = 0D0
      ELSE
         XBAR = 0D0
         YBAR = 0D0
      ENDIF
C
      DO NV=1, NBVECT
C
         INTVEF  = 0D0
         INTVERR = 0D0
C
         DO J=1,4
C
C           NUMERO DU NOEUD J DU TRIANGLE
            NSJ = NONOEF( J )
C
C           NUMERO DU DL DE VITESSE AU NOEUD AVANT NSJ
            NDLJ = NDDLNV( NSJ-1 )
C
            DO I=1,4
C
C              NUMERO DU NOEUD I DU TRIANGLE
               NSI = NONOEF( I )
C
C              NUMERO DU DL DE VITESSE AU NOEUD AVANT NSI
               NDLI = NDDLNV( NSI-1 )
C
C              NORME L1 DE VITESSE**2 INTEGREE EXACTEMENT
               DO K=1,2
C                 COMPOSANTE K DE LA VITESSE
                  INTVEF = INTVEF + P1BP1B2D(I,J) * VXYPR(NDLI+K,NV)
     %                                            * VXYPR(NDLJ+K,NV)
                  IF( NOFOVI .GT. 0 ) THEN
C                    CALCUL DE LA VITESSE_EXACTE(temps,x,y,z,k) AU NOEUD J
                     PARAMF(1) = TIMES(NV)
                     IF( J .EQ. 4 ) THEN
                        PARAMF(2) = XBAR
                        PARAMF(3) = YBAR
                     ELSE
                        PARAMF(2) = X(J,1)
                        PARAMF(3) = X(J,2)
                     ENDIF
                     PARAMF(5) = K
                     CALL FONVAL( NOFOVI, 5, PARAMF, NCODEV, VitExkJ )
                     IF( NCODEV .LE. 0 ) RETURN
C
C                    CALCUL DE LA VITESSE_EXACTE(temps,x,y,z,k) AU NOEUD I
C                    PARAMF(1) = TIMES(NV)
                     IF( I .EQ. 4 ) THEN
                        PARAMF(2) = XBAR
                        PARAMF(3) = YBAR
                     ELSE
                        PARAMF(2) = X(I,1)
                        PARAMF(3) = X(I,2)
                     ENDIF
C                    PARAMF(5) = K
                     CALL FONVAL( NOFOVI, 5, PARAMF, NCODEV, VitExkI )
                     IF( NCODEV .LE. 0 ) RETURN
C
C                    COMPOSANTE K DE LA VITESSE
                     INTVERR = INTVERR + P1BP1B2D(I,J)
     %                                 * ( VitExkI-VXYPR(NDLI+K,NV) )
     %                                 * ( VitExkJ-VXYPR(NDLJ+K,NV) )
                  ENDIF
               ENDDO
C
            ENDDO
C
         ENDDO
C
C        CARRE DE LA NORME L2 DE LA VITESSE INTEGREE EXACTEMENT
         INTVIT(NV) = INTVIT(NV) + INTVEF  * DELTA
         INTVER(NV) = INTVER(NV) + INTVERR * DELTA
C
      ENDDO
C
      RETURN
      END
