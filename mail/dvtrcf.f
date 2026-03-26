      SUBROUTINE DVTRCF( N1ARCF, NSUIV, CF, PXYD )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACE DU CONTOUR FERME CF
C -----
C ENTREES:
C --------
C N1ARCF : NUMERO DANS CF DE LA PREMIERE ARETE DU CF
C NSUIV  : NUMERO DU PREMIER INDICE DONNANT L'ARETE SUIVANTE
C CF     : CF(1,NA) NUMERO DU SOMMET DE L'ARETE NA DU CF
C          CF(2,NA) NUMERO DE L'ARETE SUIVANTE DANS CF
C          NA=1 EST LA PREMIERE ARETE DU CF QUI EST PARCOURU SELON
C          LE SENS DES AIGUILLES D'UNE MONTRE
C PXYD   : TABLEAU DES COORDONNEES 2D DES POINTS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC       JUIN 1994
C....................................................................012
      include"./incl/trvari.inc"
      INTEGER           CF(1:2,*)
      DOUBLE PRECISION  PXYD(3,*)
C
C     TRACE OU PAS DE TRACE ?
      LOGICAL           TRATRI
      COMMON / DV2DCO / TRATRI
C
      IF( .NOT. TRATRI ) RETURN
      IF( N1ARCF .LE. 0 ) RETURN
C
      NA0 = N1ARCF
      NA1 = N1ARCF
      NA2 = CF( NSUIV, NA1 )
      IF( NA2 .LE. 0 ) RETURN
      NS1 = CF( 1, NA1 )
      X1  = REAL( PXYD(1,NS1) )
      Y1  = REAL( PXYD(2,NS1) )
C
C     LA COULEUR DE TRACE DU CF
      NCOUL = CF( 1, NA2 )
      IF( NCOUL .GT. NDCOUL ) NCOUL = MOD( NCOUL, NDCOUL )
      IF( NCOUL .LE. NDCORE ) NCOUL = NDCORE + 1
C
C     TRACE DE L'ARETE
 10   NS2 = CF( 1, NA2 )
      X2  = REAL( PXYD(1,NS2) )
      Y2  = REAL( PXYD(2,NS2) )
      CALL TRAIT2D( NCOUL, X1,Y1, X2,Y2 )
C
C     LE NUMERO DU SOMMET
      CALL ENTIER2D( NCOUL, X1,Y1, NS1 )
C
C     LE NUMERO DU SOMMET
      CALL ENTIER2D( NCOUL, X2,Y2, NS2 )
C
C     L'ARETE SUIVANTE
      NS1 = NS2
      X1  = X2
      Y1  = Y2
      NA2 = CF( NSUIV, NA2 )
      IF( NA2 .NE. NA0 ) GOTO 10
C
C     TRACE DE LA DERNIERE ARETE
      NS2 = CF( 1, NA2 )
      X2  = REAL( PXYD(1,NS2) )
      Y2  = REAL( PXYD(2,NS2) )
      CALL TRAIT2D(  NCOUL, X1,Y1, X2,Y2 )
      CALL ENTIER2D( NCOUL, X1,Y1, NS1 )
      CALL ENTIER2D( NCOUL, X2,Y2, NS2 )
C
      RETURN
      END
