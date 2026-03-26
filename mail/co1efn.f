      SUBROUTINE CO1EFN( NBNOEF, NBELEM,  NUNOEF, NBNOVI, NU1EFN )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CONSTRUIRE LE TABLEAU NU1EFN No 1EF pour chaque NOEUD
C -----
C ENTREES:
C --------
C NBNOEF : NOMBRE DE NOEUDS D'UN TETRAEDRE = 10
C NBELEM : NOMBRE DE TETRAEDRES DU MAILLAGE
C NUNOEF : NUMERO DES 4 SOMMETS ET 6 MILIEUX DES ARETES DES TETRAEDRES
C NBNOVI : NOMBRE DE NOEUDS DU MAILLAGE

C SORTIES:
C --------
C NU1EFN : NUMERO D'UN EF CONTENANT LE NOEUD
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St Pierre du Perray Decembre 2012
C23456---------------------------------------------------------------012
      IMPLICIT  NONE
      INTEGER   NBNOVI, NBNOEF, NBELEM, NUNOEF(NBELEM,NBNOEF),
     %          NU1EFN(NBNOVI), NUELEM, K

      DO NUELEM=1,NBELEM
         DO K=1,NBNOEF
            NU1EFN( NUNOEF(NUELEM,K) ) = NUELEM
         ENDDO
      ENDDO

      RETURN
      END
