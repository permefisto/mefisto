      SUBROUTINE SGIPERM( L1CH, LP1SGI, LCHSGI )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  PERMUTER CIRCULAIREMENT LES SGI POUR AMENER EN TETE LE SGI L1CH
C -----
C MODIFIES:
C----------
C LP1SGI : POINTEUR SUR LE 1-ER SEGMENT INTERSECTION (SGI)
C LCHSGI : NUMERO DU SGI et CHAINAGE SUR LE SUIVANT
C          LCHSGI(2,K) = SUIVANT DE K DANS LCHSGI
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET Alain LJLL UPMC & St Pierre du Perray  Octobre 2011
C2345X7..............................................................012
      INTEGER  LP1SGI, LCHSGI(2,*)
C
C     RECHERCHE DE L1CH DANS LE CHAINAGE DES SGI
C     L1CH EST SUPPOSE ETRE UN DES SGI
C     ------------------------------------------
C     LE PREMIER SGI
      LCH = LP1SGI
      IF( LCH .LE. 0 ) GOTO 9999
C
 10   IF( LCH .NE. L1CH ) THEN
C        LE SUIVANT
         LCH = LCHSGI( 2, LCH )
         IF( LCH .GT. 0 ) GOTO 10
      ENDIF
C
C     MODIFICATION DU CHAINAGE
C     LE DERNIER ACTUEL VA POINTER SUR LE PREMIER
C     -------------------------------------------
C     RECHERCHE DU DERNIER SGI DU CHAINAGE
C     LE SUIVANT DE LCH
 20   LCH0 = LCH
      LCH  = LCHSGI( 2, LCH  )
      IF( LCH .LE. 0 ) GOTO 20
C
C     LCH0 EST L'ANCIEN DERNIER, DE SUIVANT L'ANCIEN PREMIER
      LCHSGI( 2, LCH0 ) = LP1SGI
C
C     LE NOUVEAU PREMIER EST L1CH
C     ---------------------------
      LP1SGI = L1CH
C
 9999 RETURN
      END
