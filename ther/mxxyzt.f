      SUBROUTINE MXXYZT( MOARET, MXARET, LARETE,
     %                   NBPTAF, CPNFAF, COIN )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL LE MIN MAX DES COORDONNEES DES POINTS DES ARETES
C -----
C
C ENTREES :
C ---------
C MOARET : NOMBRE DE MOTS STOCKES POUR CHAQUE ARETE DU TABLEAU LARETE
C MXARET : NOMBRE MAXIMAL D'ARETES DECLARABLES DANS LARETE
C LARETE : LE TABLEAU DES ARETES SELON LE HACHAGE
C NBPTAF : NOMBRE DE POINTS PAR ARETE OU LE FLUX EST CALCULE
C CPNFAF : 3 COORDONNEES DES POINTS D'INTEGRATION ET NORMALES SUR LES ARETES
C
C SORTIES :
C ---------
C COIN   : MIN ET MAX DES COORDONNEES DES POINTS DES ARETES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS         MAI 1995
C23456---------------------------------------------------------------012
      include"./incl/trvari.inc"
      INTEGER     LARETE(1:MOARET,1:MXARET)
      REAL        CPNFAF(1:3,1:NBPTAF,1:MXARET,1:2),
     %            COIN(6,2)
C
C     INITIALISATIONS
      RMAX = RINFO( 'GRAND' )
      DO 10 I=1,3
C        LE MINIMUM
         COIN(I,1) =  RMAX
C        LE MAXIMUM
         COIN(I,2) = -RMAX
 10   CONTINUE
C
      DO 50 N=1,MXARET
         IF( LARETE(1,N) .NE. 0 ) THEN
C           L'ARETE EXISTE
            DO 40 J=1,NBPTAF
               DO 30 I=1,3
C                 LE MINIMUM DES COORDONNEES
                  S = CPNFAF( I, J, N, 1 )
                  COIN(I,1) = MIN( COIN(I,1) , S )
C                 LE MAXIMUM DES COORDONNEES
                  COIN(I,2) = MAX( COIN(I,2) , S )
 30            CONTINUE
 40         CONTINUE
         ENDIF
 50   CONTINUE
C
C     RECHERCHE DE LA DIMENSION MAXIMALE ECAMAX DE L'OBJET
      ECAMAX = 0.
      DO 60  I=1,NDIMLI
         ECAMAX = MAX( ECAMAX, COIN(I,2) - COIN(I,1) )
 60   CONTINUE
      IF( ECAMAX .LE. 0. ) ECAMAX = 1.
C     UNE MARGE DE 2%
      ECAMAX = ECAMAX * 0.51
C
      DO 70 I=1,NDIMLI
         RMAX = ( COIN(I,1) + COIN(I,2) ) * 0.5
         COIN(I,1) = RMAX - ECAMAX
         COIN(I,2) = RMAX + ECAMAX
 70   CONTINUE
C
      IF( NDIMLI .EQ. 2 ) THEN
         COIN(3,1) = 0
         COIN(3,2) = 0
      ENDIF
      END
