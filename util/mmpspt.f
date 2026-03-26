      SUBROUTINE MMPSPT( NBCOOR, NBPOIN , POINTS, XYZPOI , HMMX)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    METTRE LE MINIMUM ET MAXIMUM DANS LE TABLEAU HMMX
C -----    DES COORDONNEES DES POINTS PROJETES SUR LA DROITE DEFINIE
C          PAR UN POINT ET UN VECTEUR
C
C
C ENTREES :
C ---------
C NBCOOR : NOMBRE DE COORDONNEES PAR POINT (3 ou 6)
C NBPOIN : NOMBRE DE POINTS
C XYZPOI : 3 COORDONNEES DES NBPOIN POINTS
C POINTS : POINTS(1:3,1) : LE POINT
C        : POINTS(1:3,2) : LE VECTEUR
C
C SORTIE :
C --------
C HMMX   : HMMX(1) MIN DES COORDONNEES DES POINTS
C          HMMX(2) MAX DES COORDONNEES DES POINTS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL HAVE DEA ANALYSE NUMERIQUE UPMC PARIS     JANVIER 2000
C23456---------------------------------------------------------------012
      REAL  HMMX(2), XYZPOI(NBCOOR,NBPOIN)
      REAL  POINTS(3,2)
C
C     L'INITIALISATION
      HMMX(1) = (XYZPOI(1,1)-POINTS(1,1)) * POINTS(1,2)
     %        + (XYZPOI(2,1)-POINTS(2,1)) * POINTS(2,2)
     %        + (XYZPOI(3,1)-POINTS(3,1)) * POINTS(3,2)
      HMMX(2) = HMMX(1)
C
C     LE MIN ET MAX
      DO 20 I=2,NBPOIN
            HMMX(1) = MIN( (XYZPOI(1,I)-POINTS(1,1)) * POINTS(1,2)
     %                   + (XYZPOI(2,I)-POINTS(2,1)) * POINTS(2,2)
     %                   + (XYZPOI(3,I)-POINTS(3,1)) * POINTS(3,2)
     %                   , HMMX(1) )
            HMMX(2) = MAX( (XYZPOI(1,I)-POINTS(1,1)) * POINTS(1,2)
     %                   + (XYZPOI(2,I)-POINTS(2,1)) * POINTS(2,2)
     %                   + (XYZPOI(3,I)-POINTS(3,1)) * POINTS(3,2)
     %                   , HMMX(2) )
 20   CONTINUE
C
      END
