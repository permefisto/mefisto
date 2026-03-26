      SUBROUTINE FLUXVITSF( NBNOVI, XYZNOE,
     %                      NUTYEL, NBNOEF, NBELEM, NUNOEF,
     %                      NUDDL,  NBSFEF, NUSFEF,
     %                      VXYZPN, NUMISF, NUMASF,
     %                      FLUNEG, FLUPOS, FLUVSF )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LE FLUX 3D DE LA VITESSE A TRAVERS LES FACES
C -----    FRONTALIERES DECLAREES DANS LA DEFINITION DU FLUIDE
C
C ENTREES:
C --------
C NBNOVI : NOMBRE TOTAL DE NOEUDS VITESSE
C XYZNOE : XYZNOE(3,NBNOVI) 3 COORDONNEES DES NOEUDS DE LA TETREDRISATION
C NBNOEF : NOMBRE DE NOEUDS D'UN EF DE CE TYPE
C NBELEM : NOMBRE D'EF DU MAILLAGE
C NUNOEF : NUNOEF(NBELEM,NBNOEF) NO DES NOEUDS DES NBELEM EF
C NUDDL  : NUMERO DU DERNIER DL DE CHAQUE NOEUD
C VXYZPN : VITESSE PRESSION EN TOUS LES NOEUDS DU MAILLAGE
C NUMISF : NUMERO MINIMAL DES SURFACES DANS LA DEFINITION DU FLUIDE
C NUMASF : NUMERO MAXIMAL DES SURFACES DANS LA DEFINITION DU FLUIDE
C
C SORTIES:
C --------
C FLUNEG : FLUX NEGATIF DE LA VITESSE A TRAVERS LES FACES FRONTALIERES
C FLUPOS : FLUX POSITIF DE LA VITESSE A TRAVERS LES FACES FRONTALIERES
C FLUVSF : FLUX DE LA VITESSE A TRAVERS CHACUNE DES SURFACES DU FLUIDE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY     Juin 2010
C23456---------------------------------------------------------------012
      IMPLICIT NONE
      include"./incl/langue.inc"
      INTEGER           LECTEU, IMPRIM, NUNITE
      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)
      REAL              XYZNOE(3,NBNOVI)
      INTEGER           NUNOEF(NBELEM,NBNOEF), NUSFEF(NBSFEF,NBELEM)
      INTEGER           NUDDL(0:NBNOVI),
     %                  NBNOVI, NUTYEL, NBNOEF, NBELEM, NBSFEF,
     %                  NUMISF, NUMASF
      DOUBLE PRECISION  VXYZPN(*)
      DOUBLE PRECISION  DGL(2,3), DN(3), FLUPOS, FLUNEG, FLUXV, VN
      DOUBLE PRECISION  FLUVSF(NUMISF:NUMASF)
      INTEGER           K,I,LEDECA,NN,N1,N2,N3,NODL,NUELEM,NUDCNB,NSF
      CHARACTER*24      KNOMSF

C     NO DES 3 SOMMETS ET 3 MILIEUX DES ARETES DES 4 FACES DU TETRAEDRE
C     CE SONT AUSSI LES POINTS D'INTEGRATION
      INTEGER  NOSFTE(6,4)
      DATA     NOSFTE / 1,3,2, 7, 6,5,
     %                  1,4,3, 8,10,7,
     %                  1,2,4, 5, 9,8,
     %                  2,3,4, 6,10,9 /

C     INITIALISATION DES FLUX A ZERO
      FLUNEG = 0D0
      FLUPOS = 0D0
      DO K = NUMISF, NUMASF
         FLUVSF( K ) = 0D0
      ENDDO

C     COEFFICIENTS DE L'INTEGRALE DES POLYNOMES SUR L'EF DE REFERENCE
      IF( NUTYEL .EQ. 19 ) THEN
C        TETRAEDRE BREZZI-FORTIN   INTEGRATION EXACTE P1 SUR LE TRIANGLE
         LEDECA = 0
      ELSE
C        TETRAEDRE TAYLOR-HOOD  INTEGRATION EXACTE P2 SUR LE TRIANGLE
         LEDECA = 3
      ENDIF

C     LA BOUCLE SUR LES ELEMENTS FINIS POUR CALCULER LE FLUX A LA FRONTIERE
C     =====================================================================
      DO NUELEM = 1, NBELEM

C        CALCUL DU VECTEUR NORMAL A LA FACE K
         DO K = 1, 4

            NSF = NUSFEF(K,NUELEM)
            IF( NSF .GT. 0 ) THEN

C              RECHERCHE DU JACOBIEN DE G
               N1 = NUNOEF( NUELEM, NOSFTE(1,K) )
               N2 = NUNOEF( NUELEM, NOSFTE(2,K) )
               N3 = NUNOEF( NUELEM, NOSFTE(3,K) )

               DGL(1,1) = XYZNOE( 1, N2 ) - XYZNOE( 1, N1 )
               DGL(2,1) = XYZNOE( 1, N3 ) - XYZNOE( 1, N1 )

               DGL(1,2) = XYZNOE( 2, N2 ) - XYZNOE( 2, N1 )
               DGL(2,2) = XYZNOE( 2, N3 ) - XYZNOE( 2, N1 )

               DGL(1,3) = XYZNOE( 3, N2 ) - XYZNOE( 3, N1 )
               DGL(2,3) = XYZNOE( 3, N3 ) - XYZNOE( 3, N1 )

C              CALCUL DE LA NORMALE * DELTAK: PRODUIT VECTORIEL(DG/DX1,DG/DX2)
               DN(1) = DGL(1,2) * DGL(2,3) - DGL(2,2) * DGL(1,3)
               DN(2) = DGL(1,3) * DGL(2,1) - DGL(2,3) * DGL(1,1)
               DN(3) = DGL(1,1) * DGL(2,2) - DGL(2,1) * DGL(1,2)

C              FLUX NORMAL DE LA VITESSE A TRAVERS LA FACE K
               FLUXV = 0D0
               DO I = 1, 3
C                 NUMERO DU DL AVANT CEUX DE LA VITESSE AU NOEUD I
                  NODL = NUDDL( NUNOEF(NUELEM,NOSFTE(LEDECA+I,K)) - 1 )
C                 Vitesse Normale  DELTAK  au NOEUD I de la FACE K
                  VN = VXYZPN(NODL+1) * DN(1)
     %               + VXYZPN(NODL+2) * DN(2)
     %               + VXYZPN(NODL+3) * DN(3)
                  FLUXV = FLUXV + VN
               ENDDO
               FLUXV = FLUXV / 6D0

C              REPARTITION EN FLUX NEGATIF ET POSITIF ET POIDS
               IF( FLUXV .LE. 0D0 ) THEN
                  FLUNEG = FLUNEG + FLUXV
               ELSE
                  FLUPOS = FLUPOS + FLUXV
               ENDIF

C              FLUX DE LA VITESSE SUR LA SURFACE NOOBSF(K)
               FLUVSF(NSF) = FLUVSF(NSF) + FLUXV

            ENDIF

         ENDDO

      ENDDO

C     AFFICHAGE DU FLUX DE LA VITESSE SUR LES SURFACES DU FLUIDE
C     ===========================================================
      DO K = NUMISF, NUMASF

C        LE FLUX DE LA VITESSE SUR LA SURFACE K A-T-ELLE ETE CALCULEE
         FLUXV = FLUVSF( K )

         IF( FLUXV .NE. 0.D0 ) THEN
C           LE NOM DE LA SURFACE
            CALL NMOBNU( 'SURFACE', K, KNOMSF )
            NN = NUDCNB( KNOMSF )
            IF( LANGAG .EQ. 0 ) THEN
               WRITE(IMPRIM,*)'FLUX NORMAL de la VITESSE=',FLUXV,
     %                     ' sur la SURFACE ',KNOMSF(1:NN)
            ELSE
               WRITE(IMPRIM,*)'VELOCITY NORMAL FLUX=',FLUXV,
     %                     ' through the SURFACE ',KNOMSF(1:NN)
            ENDIF
         ENDIF

      ENDDO


C     AFFICHAGE DU TOTAL DES FLUX POSITIFS et NEGATIFS DE LA VITESSE
C     SUR LES SURFACES DU FLUIDE
      FLUXV = FLUNEG + FLUPOS
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10100) FLUNEG, FLUPOS, FLUXV
      ELSE
         WRITE(IMPRIM,20100) FLUNEG, FLUPOS, FLUXV
      ENDIF

10100 FORMAT(' Totaux sur les surfaces: FLUX<0 de la VITESSE=',
     %G14.6,'  FLUX>0=',G14.6,'  FLUX TOTAL=',G14.6)
20100 FORMAT(' Total on all Surfaces: VELOCITY NEGATIVE FLUX=',G14.6,
     %'  POSITIVE FLUX=',G14.6,'  TOTAL FLUX=',G14.6)

      RETURN
      END
