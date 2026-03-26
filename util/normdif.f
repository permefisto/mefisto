      SUBROUTINE NORMDIF( NBC, VEC0, VEC1,  NORDIF, NORVEC )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  CALCUL DE LA NORME DE LA DIFFERENCE DES 2 VECTEURS VEC0 et VEC1
C -----  / Nb COMPOSANTES
C        CALCUL DE LA NORME DU VECTEUR VEC0 / Nb COMPOSANTES
C
C ENTREES :
C ---------
C NBC      : NOMBRE DE COMPOSANTES DE CHACUN DES 2 VECTEURS
C VEC0,VEC1: 2 VECTEURS DE NBC COMPOSANTES
C
C SORTIES :
C ---------
C NORDIF  : SOMME I=1,NBC  de  |VEC1(I) - VEC0(I)|  / NBC
C NORVEC  : SOMME I=1,NBC  de  |VEC0(I)| / NBC
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET TEXAS A & M UNIVERSITY at QATAR  JANVIER 2012
C23456---------------------------------------------------------------012
      DOUBLE PRECISION  VEC0(NBC), VEC1(NBC), NORDIF, NORVEC
C
      NORDIF = 0D0
      NORVEC = 0D0
      DO I = 1, NBC
         NORDIF = NORDIF + ABS( VEC1(I) - VEC0(I) )
         NORVEC = NORVEC + ABS( VEC0(I) )
      ENDDO
C
      NORDIF = NORDIF / NBC
      NORVEC = NORVEC / NBC
C
      RETURN
      END
