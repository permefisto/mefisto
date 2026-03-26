      SUBROUTINE FL2VITTH3( NBPOE,  VOLTET, NONOEF,
     %                      NTDLVP, NBVECT, VXYZPR, NDDLNV,
     %                      INTVIT, NOFOVI, NBNOVI, VITEXA, INTVER )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DE L'INTEGRALE DES NBVECT |VITESSE|**2
C -----    SUR UN TETRAEDRE DE TAYLOR-HOOD DE VITESSE AVEC 3 COMPOSANTES
C          POLYNOME DE LAGRANGE DE DEGRE 2
C
C ENTREES:
C --------
C NBPOE  : NOMBRE DE NOEUDS=POINTS DE L'EF
C VOLTET : VOLUME DU TETRAEDRE
C NONOEF : NUMEROS DES NBNOEF SOMMETS DU TETRAEDRE COURANT
C NTDLVP : NOMBRE TOTAL DE DEGRES DE LIBERTE DU MAILLAGE DE L'OBJET
C NBVECT : NOMBRE DE VECTEURS SOLUTION VITESSEPRESSION
C VXYZPR : NBVECT VECTEURS VITESSEPRESSION SUR LE MAILLAGE
C NDDLNV : TABLEAU DES POINTEURS SUR LE DERNIER DL POUR CHAQUE NOEUD VITESSE
C          NDDLNV(I)= NO DERNIER DL DU NOEUD (SOMMET ou MILIEU ou BARYCENTRE)
C          CE TABLEAU EST DIMENSIONNE A 0:NBNOVI
C NOFOVI : NUMERO DE LA FONCTION UTILISATEUR VITESSE_EXACTE(t,x,y,z,nocomp)
C          0 SI ELLE N'EXISTE PAS => INTVER N'EST PAS CALCULE
C NBNOVI : NOMBRE DE NOEUDS SUPPORT DE LA VITESSE SUR LE MAILLAGE
C VITEXA : VITEXA(NBNOVI,NBVECT,NDIM) VITESSE EXACTE EN CHAQUE NOEUD TEMPS
C
C SORTIES:
C --------
C                                                              (Vke1 )
C INTVIT : Som  Jacobien Som (Vke1,...,Vke10) [int P2i P2j DX] (  .  )
C          e de E        k=1 ... 10                            (Vke10)
C                                                                       (Vk1E -
C INTVER : Som  Som Jacobien (Vk1E-Vk1C,...,Vk10E-Vk10C) [int Pi Pj dX] (   ...
C       e de E  k=1,...,3                                               (Vk10E-V
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET TEXAS A & M UNIVERSITY at QATAR   Janvier 2012
C23456---------------------------------------------------------------012
C     INTEGRALES P2 P2 SUR LE TETRAEDRE P2 de REFERENCE
      include"./incl/p2p23d.inc"
CCC   DOUBLE PRECISION  P2P23D(10,10)
CCCC                    P2P23D(i,j) = integrale P2i P2j dX
C
      INTEGER           NONOEF(NBPOE), NDDLNV(0:NBNOVI)
      DOUBLE PRECISION  VXYZPR(NTDLVP,NBVECT),
     %                  VITEXA(NBNOVI,NBVECT,3),
     %                  INTVIT(NBVECT), INTVER(NBVECT),
     %                  VOLTET, DELTA,  INTVEF, INTVERR,
     %                  VitExkI, VitExkJ(3)
C
C     JACOBIEN = 6 x (VOLUME DU TETRAEDRE)
      DELTA = VOLTET * 6D0
C
      DO NV=1, NBVECT
C
         INTVEF  = 0D0
         INTVERR = 0D0
         DO J=1,10
C
C           NUMERO DU NOEUD J DU TETRAEDRE
            NSJ = NONOEF( J )
C
C           NUMERO DU DL DE VITESSE AU NOEUD AVANT NSJ
            NDLJ = NDDLNV( NSJ-1 )
C
            IF( NOFOVI .GT. 0 ) THEN
C              LA VITESSE_EXACTE(temps,x,y,z,k) AU NOEUD J
               DO K=1,3
                  VitExkJ(K) = VITEXA(NSJ,NV,K)
               ENDDO
            ENDIF
C
            DO I=1,10
C
C              NUMERO DU NOEUD I DU TETRAEDRE
               NSI = NONOEF( I )
C
C              NUMERO DU DL DE VITESSE AU NOEUD AVANT NSI
               NDLI = NDDLNV( NSI-1 )
C
C              CARRE DE LA NORME L2 DE |VITESSE| INTEGREE EXACTEMENT
               DO K=1,3
C
C                 COMPOSANTE K DE LA VITESSE
                  INTVEF = INTVEF + P2P23D(I,J) * VXYZPR(NDLI+K,NV)
     %                                          * VXYZPR(NDLJ+K,NV)
C
                  IF( NOFOVI .GT. 0 ) THEN
C
C                    LA VITESSE_EXACTE(temps,x,y,z,k) AU NOEUD I
                     VitExkI = VITEXA(NSI,NV,K)
C
C                    COMPOSANTE K DE LA VITESSE
                     INTVERR = INTVERR+ P2P23D(I,J)
     %                                * ( VitExkI   -VXYZPR(NDLI+K,NV) )
     %                                * ( VitExkJ(K)-VXYZPR(NDLJ+K,NV) )
                  ENDIF
C
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
