      SUBROUTINE MIMXCO( COIN )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    METTRE LE MINIMUM ET MAXIMUM DES COORDONNEES ACTUELLES
C -----    DU TABLEAU COIN
C
C SORTIE :
C --------
C COIN   : COIN(.,1) MIN DES COORDONNEES DES SOMMETS ACTUELS
C          COIN(.,2) MAX DES COORDONNEES DES SOMMETS ACTUELS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE  UPMC PARIS  OCTOBRE 1988
C.......................................................................
      include"./incl/xyzext.inc"
      REAL COIN(6,2)
C
C     LE TABLEAU COINCO DU COMMON / XYZEXT / CONTIENT LES EXTREMES
C     DES COORDONNEES ACTUELLES
      DO 20 J=1,2
         DO 10 I=1,6
            COIN( I , J ) = COOEXT( I , J )
 10      CONTINUE
 20   CONTINUE
      END
