      SUBROUTINE ASMEEX2( NDIM, NBNOEF, NONOEF, BE1, BE2, NBNOVI, VG )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    ASSEMBLAGE DES 2 VECTEURS ELEMENTAIRES DANS LE VECTEUR GLOBAL
C -----
C
C ENTREES:
C --------
C NDIM   : NOMBRE DE COLONNES DES VECTEURS STOCKES DANS VG
C NBNOEF : NOMBRE DE DL (4 TRIANGLE BF, 5 TETRAEDRE BF, ...)
C NONOEF : NUMERO GLOBAL DES NBNOEF NOEUDS DE L'EF
C BE1,BE2: 2 VECTEURS ELEMENTAIRES
C NBNOVI : NOMBRE DE DL D'UNE COMPOSANTE DU VECTEUR VG
C
C SORTIES:
C --------
C VG     : NDIM COMPOSANTES MODIFIEES DU VECTEUR GLOBAL
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC & St Pierre du Perray  Decembre 2010
C23456---------------------------------------------------------------012
C$ USE OMP_LIB
      DOUBLE PRECISION  BE1(NBNOEF,NDIM),
     %                  BE2(NBNOEF,NDIM),
     %                  VG(NBNOVI,NDIM), S
      INTEGER           NONOEF(NBNOEF)
C
      DO I=1,NBNOEF
C
C        NUMERO DU I-EME NOEUD DE L'ELEMENT FINI
         NS = NONOEF(I)
C
C        ASSEMBLAGE DU COEFFICIENT ELEMENTAIRE DANS LE COEFFICIENT GLOBAL
         DO K=1,NDIM
            S =  BE1(I,K) + BE2(I,K)
!$OMP ATOMIC
            VG( NS, K ) = VG( NS, K ) + S
         ENDDO
C
      ENDDO
C
      RETURN
      END
