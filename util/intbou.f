      SUBROUTINE INTBOU( NCENTR , NPS , NP , D )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    LE POINT NP EST IL DANS , SUR OU EXTERIEUR A LA BOULE
C -----    DE CENTRE NCENTR ?
C
C ENTREES:
C --------
C NCENTR : LES 3 COORDONNEES ENTIERES DU CENTRE DU CERCLE
C NPS    : LES 3 COORDONNEES D'UN POINT SUR LA BOULE
C NP     : LES 3 COORDONNEES DU POINT A  TRAITER
C
C SORTIE :
C --------
C D      : <0 DEDANS      LA BOULE
C          =0 SUR         LA BOULE
C          >0 EXTERIEUR A LA BOULE
C
C EN FAIT :
C D = (X-XS)(X+XS-2*XC)+ (Y-YS)(Y+YS-2*YC)+ (Z-ZS)(Z+ZS-2*ZC)
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LABORATOIRE ANALYSE NUMERIQUE PARIS NOVEMBRE 1988
C...............................................................................
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      INTEGER           NPS(3),NP(3),NCENTR(3)
      DOUBLE PRECISION  D1,D2,D
C
      D1 = NP(1) - NPS(1)
      D2 = NP(1) - NCENTR(1) + NPS(1) - NCENTR(1)
      D  = D1 * D2
C
      D1 = NP(2) - NPS(2)
      D2 = NP(2) - NCENTR(2) + NPS(2) - NCENTR(2)
      D  = D + D1 * D2
C
      D1 = NP(3) - NPS(3)
      D2 = NP(3) - NCENTR(3) + NPS(3) - NCENTR(3)
      D  = D + D1 * D2
C
C     ATTENTION : D PEUT ETRE > MAXIMUM( ENTIERS )
      END
