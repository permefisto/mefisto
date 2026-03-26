      SUBROUTINE MIMXPT( NBCOOR, NBPOIN, XYZPOI, COIN )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LE MINIMUM ET MAXIMUM DES 3-ERES COORDONNEES
C -----    DES POINTS DANS LE TABLEAU COIN

C ENTREES :
C ---------
C NBCOOR : NOMBRE DE COORDONNEES D'UN POINT
C NBPOIN : NOMBRE DE POINTS
C XYZPOI : LES NBCOOR COORDONNEES DES NBPOIN POINTS

C SORTIE :
C --------
C COIN   : COIN(.,1) MIN DES NBCOOR COORDONNEES DES POINTS
C          COIN(.,2) MAX DES NBCOOR COORDONNEES DES POINTS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE  UPMC PARIS   NOVEMBRE 1989
C MODIFS : ALAIN PERRONNET  TEXAS A & M UNIVERSITY          JUILLET 2005
C.......................................................................
      REAL  COIN(6,2), XYZPOI(NBCOOR,NBPOIN)

C     L'INITIALISATION
      DO I=1,NBCOOR
         COIN(I,1) = XYZPOI(I,1)
         COIN(I,2) = XYZPOI(I,1)
      ENDDO

C     LE MIN ET MAX
      DO J=2,NBPOIN
         DO I=1,NBCOOR
            COIN(I,1) = MIN( XYZPOI(I,J) , COIN(I,1) )
            COIN(I,2) = MAX( XYZPOI(I,J) , COIN(I,2) )
         ENDDO
      ENDDO

      RETURN
      END
