      SUBROUTINE EMPIL3 ( LHPILE , MXPILE , NPILE , K1 , K2 , K3 )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : EMPILER LES 3 ENTIERS K1 K2 K3 DANS LA PILE NPILE
C -----
C
C PARAMETRES D ENTREE :
C ---------------------
C MXPILE : NOMBRE MAXIMAL DE TERMES DANS LA PILE
C K1     : 1-ER TERME A EMPILER
C K2     : 2-ME TERME A EMPILER
C K3     : 3-ME TERME A EMPILER
C
C PARAMETRES MODIFIES :
C ----------------------
C LHPILE : POINTEUR SUR LE SOMMET DE LA PILE
C NPILE  : PILE ( 3 , MXPILE ) LA PILE D ENTIERS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 1990
C23456---------------------------------------------------------------012
      include"./incl/gsmenu.inc"
      INTEGER        NPILE(3,MXPILE)
C
C     LA PILE EST-ELLE SATUREE ?
C     ==========================
      IF( LHPILE .LT. MXPILE ) THEN
C
C        K1,K2,K3 SONT EMPILES EN POSITION LHPILE + 1
         LHPILE = LHPILE + 1
         NPILE( 1 , LHPILE ) = K1
         NPILE( 2 , LHPILE ) = K2
         NPILE( 3 , LHPILE ) = K3
C
      ELSE
C
C        OUI.DIAGNOSTIC ET ARRET
C        =======================
         NBLGRC(NRERR) = 2
         KERR(1) = 'EMPIL3: PILE SATUREE'
         KERR(2) = 'AUGMENTER MXPILE'
         CALL LEREUR
      ENDIF
      END
