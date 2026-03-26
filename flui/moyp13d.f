      SUBROUTINE MOYP13D( NBNOVI, XYZNOE, NBELEM, NBNOEF, NUNOEF,NUNOSO,
     %                    NBSOM,  P1SOLU,  VOLUME, INTP1 )
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! BUT:  CALCULER LE VOLUME DU MAILLAGE ET
! ----  L'INTEGRALE DE L'INTERPOLATION P1 de la PRESSION SUR LE MAILLAGE

! ENTREES:
! --------
! NBNOVI : NOMBRE DE NOEUDS VITESSE SOMMETS et MILIEUX DES ARETES DES EF
! XYZNOE : XYZNOE(3,NBNOVI)  3 COORDONNEES DES NOEUDS DU MAILLAGE
! NBELEM : NOMBRE DE TETRAEDRES DU MAILLAGE
! NBNOEF : NOMBRE DE NOEUDS D'UN TETRAEDRE (P1B=>4, P2=>10)
! NUNOEF : NO DES NBNOEF NOEUDS DES NBELEM TETRAEDRES
! NUNOSO : NUNOSO(I) = NUMERO DU SOMMET 1 A NBSOM DU NOEUD GLOBAL I
! P1SOLU : INTERPOLATION P1 AUX SOMMETS DE LA TETRAEDRISATION

! SORTIES:
! --------
! VOLUME : VOLUME DU MAILLAGE 3D
! INTP1  : INTEGRALE DE P1SOLU SUR LE MAILLAGE
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! AUTEUR : ALAIN PERRONNET LJLL UPMC & St Pierre du Perray Decembre 2012
!23456---------------------------------------------------------------012
      IMPLICIT NONE
      INTEGER           NBSOM, NBNOVI, NBELEM, NBNOEF,
     %                  NUNOEF(NBELEM,NBNOEF), NUNOSO(NBNOVI)
      REAL              XYZNOE(3,NBNOVI)
      DOUBLE PRECISION  P1SOLU(NBSOM), VOLUME, INTP1
      DOUBLE PRECISION  DF(3,3), DET, X, Y, Z
      INTEGER           NUELEM, I, J, K, L

      VOLUME = 0D0
      INTP1  = 0D0

!     BOUCLE SUR LES TETRAEDRES
      DO NUELEM = 1, NBELEM

!        NUMERO DES 4 NOEUDS-SOMMETS DU TETRAEDRE NUELEM
         I = NUNOEF(NUELEM,1)
         J = NUNOEF(NUELEM,2)
         K = NUNOEF(NUELEM,3)
         L = NUNOEF(NUELEM,4)

!        CONSTRUCTION DE LA MATRICE DF
         X = XYZNOE(1,I)
         Y = XYZNOE(2,I)
         Z = XYZNOE(3,I)

         DF(1,1) = XYZNOE(1,J) - X
         DF(2,1) = XYZNOE(1,K) - X
         DF(3,1) = XYZNOE(1,L) - X

         DF(1,2) = XYZNOE(2,J) - Y
         DF(2,2) = XYZNOE(2,K) - Y
         DF(3,2) = XYZNOE(2,L) - Y

         DF(1,3) = XYZNOE(3,J) - Z
         DF(2,3) = XYZNOE(3,K) - Z
         DF(3,3) = XYZNOE(3,L) - Z

!        LE DETERMINANT DE DF
         DET = DF(1,1) * ( DF(2,2) * DF(3,3) - DF(3,2) * DF(2,3) )
     %       + DF(2,1) * ( DF(3,2) * DF(1,3) - DF(1,2) * DF(3,3) )
     %       + DF(3,1) * ( DF(1,2) * DF(2,3) - DF(2,2) * DF(1,3) )

!        VOLUME DU TETRAEDRE * 6
         VOLUME = VOLUME + DET

         X = 0D0
         DO I=1,4
            K = NUNOEF(NUELEM,I)
            IF( NBNOEF .EQ. 10 ) K = NUNOSO( K )
            X = X + P1SOLU( K )
         ENDDO
         INTP1 = INTP1 + DET * X / 24D0

      ENDDO   ! FIN DE LA BOUCLE DES TETRAEDRES DU MAILLAGE

      VOLUME = VOLUME / 6D0

      print *,'moyp13d.f: VOLUME du MAILLAGE=',VOLUME,
     %' INTEGRALE Pression=',INTP1,' INTEGRALE Pression/VOLUME=',
     %INTP1/VOLUME

      RETURN
      END
