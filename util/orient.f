      SUBROUTINE ORIENT( X , Y )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    ORIENTER LE PARCOURS DES SOMMETS D'UN QUADRILATERE
C -----
C
C ENTREES :
C ---------
C X Y  : COORDONNEES DES SOMMETS DU QUADRILATERE
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY     ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 1990
C23456---------------------------------------------------------------012
      REAL  X(4), Y(4)
C
C     ORIENTATION
C     -----------
C     LES SOMMETS 3 ET 4 SONT-ILS DU MEME COTE DE 1 - 2 ?
      SIGNE3 = ( X(2) - X(1) ) * ( Y(3) - Y(1) )
     %       - ( X(3) - X(1) ) * ( Y(2) - Y(1) )
      SIGNE4 = ( X(2) - X(1) ) * ( Y(4) - Y(1) )
     %       - ( X(4) - X(1) ) * ( Y(2) - Y(1) )
      IF ( SIGNE3 * SIGNE4 .LE. 0.0 ) THEN
C        NON : ECHANGE DE 1 ET 4
         X4 = X(1)
         Y4 = Y(1)
         X(1) = X(4)
         Y(1) = Y(4)
         X(4) = X4
         Y(4) = Y4
      ENDIF
C     LE CENTRE DE GRAVITE DU QUADRILATERE
      XG = ( X(1)+X(2)+X(3)+X(4) ) / 4
      YG = ( Y(1)+Y(2)+Y(3)+Y(4) ) / 4
C     LE TRIANGLE 1 - 2 - G
      SENS1 = ( X(2) - XG ) * ( Y(1) - YG )
     %      - ( X(1) - XG ) * ( Y(2) - YG )
C     LE TRIANGLE 3 - 4 - G
      SENS2 = ( X(4) - XG ) * ( Y(3) - YG )
     %      - ( X(3) - XG ) * ( Y(4) - YG )
C     LE QUADRILATERE EST-IL CONVEXE ?
      IF ( SENS1 * SENS2 .LE. 0.0 ) THEN
C        NON : ECHANGE DE 3 ET 4
         X3 = X(3)
         Y3 = Y(3)
         X(3) = X(4)
         Y(3) = Y(4)
         X(4) = X3
         Y(4) = Y3
      ENDIF
C
      RETURN
      END
