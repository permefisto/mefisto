      SUBROUTINE AMPXYZ
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  AMPLIFIER LES COORDONNEES POUR AVOIR UN TRACE HOMOGENE EN X Y Z
C -----
C
C MODIFS:   Dans $MEFISTO/incl/xyzext.inc
C -------
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR:ALAIN PERRONNET LJLL UPMC PARIS & St Pierre du Perray JUIN 2009
C23456---------------------------------------------------------------012
      include"./incl/xyzext.inc"
C
      IF( INIEXT .EQ. 0 ) RETURN
C
C     ECARTEMENT DES ABSCISSES
      ECARTX = COOEXT(1,2) - COOEXT(1,1)
      IF( ECARTX .LE. 0 ) RETURN
C
C     LES ECARTS EN Y ET Z SONT AMPLIFIES POUR ATTEINDRE ECARTX
      MOAXYZ = 1
      XYZAMPLI(1) = 1.0
C
      ECART = COOEXT(2,2) - COOEXT(2,1)
      IF( ECART .GT. 0 ) THEN
         XYZAMPLI(2) = ECARTX / ECART
      ELSE
         XYZAMPLI(2) = 1.0
      ENDIF
      COOEXT(2,1) = COOEXT(2,1) * XYZAMPLI(2)
      COOEXT(2,2) = COOEXT(2,2) * XYZAMPLI(2)
C
      ECART = COOEXT(3,2) - COOEXT(3,1)
      IF( ECART .GT. 0 ) THEN
         XYZAMPLI(3) = ECARTX / ECART
      ELSE
         XYZAMPLI(3) = 1.0
      ENDIF
      COOEXT(3,1) = COOEXT(3,1) * XYZAMPLI(3)
      COOEXT(3,2) = COOEXT(3,2) * XYZAMPLI(3)
C
      RETURN
      END
