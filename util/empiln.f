      SUBROUTINE EMPILN ( IAPILE , MXPILE , NPILE , N  , NOS )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    EMPILER LES N ENTIERS DANS LA PILE NPILE
C -----
C
C ENTREES:
C --------
C MXPILE : NOMBRE MAXIMAL DE N - UPLETS DANS LA PILE
C N      : NOMBRE D ENTIERS A EMPILER
C NOS    : TABLEAU DE CES ENTIERS
C
C MODIFIES:
C ---------
C IAPILE : POINTEUR SUR LE SOMMET DE LA PILE
C          -1 EN SORTIE EN CAS DE SATURATION DE LA PILE
C NPILE  : PILE ( N , MXPILE ) LA PILE DE N - UPLETS ENTIERS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     OCTOBRE 1994
C23456---------------------------------------------------------------012
      INTEGER           NPILE(N,1)
      INTEGER           NOS(N)
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
10010 FORMAT('ERREUR EMPILN: PILE SATUREE MXPILE=',I13)
C
C     LA PILE EST-ELLE SATUREE ?
C     ==========================
      IF( IAPILE .GE. MXPILE ) THEN
C
C        OUI.DIAGNOSTIC ET ARRET
         WRITE(IMPRIM,10010) MXPILE
         IAPILE = -1
C
      ELSE
C
C        NON : NOS( . ) SONT EMPILES EN POSITION IAPILE + 1
         IAPILE = IAPILE + 1
         DO 11 I = 1 , N
            NPILE( I , IAPILE ) = NOS( I )
11       CONTINUE
      ENDIF
      END
