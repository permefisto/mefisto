      SUBROUTINE MXTRIE( NBPT, PT,
     %                   NP0, NP1, NP2, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     RECHERCHE DES MEILLEURS 3 POINTS DES NBPT POINTS PT
C -----     POUR DEFINIR UN TRIEDRE DONT LE PLAN DE POINTS
C           NP0 NP1 NP2 SOIT LE PLAN DES NBPT POINTS PT
C
C ENTREES :
C ---------
C NBPT   : NOMBRE DE POINTS
C PT     : LES 3 COORDONNEES DES NBPT POINTS
C
C SORTIES :
C ---------
C NP0, NP1, NP2 : LES 3 POINTS QUI DONNENT LA VALEUR ABSOLUE MINIMALE
C                 DE COSINUS( NP0-NP1, NP0-NP2 )
C IERR   : 0 S'IL EXISTE CES 3 POINTS
C          1 SI LES POINTS SONT ALIGNES OU CONFONDUS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     JANVIER 1993
C23456...............................................................012
      PARAMETER ( COSEUI=0.9962 )
C     COSEUI : SEUIL DU COSINUS DE L'ANGLE
C              ( 0.99     => 8.11 DEGRES )
C              ( 0.9962   => 5    DEGRES )
C              ( 0.99756  => 4    DEGRES )
C              ( 0.99863  => 3    DEGRES )
C              ( 0.999    => 2.56 DEGRES )
C              ( 0.9999   => 0.8  DEGRES )
      REAL   PT(3,NBPT)
C
      COSMIN = 2
      DO 50 N0 = 1, NBPT/2
         N1 = N0 + 1
         DO 40 N2=N1+1,NBPT
C           RECHERCHE DU COSINUS MINIMAL EN VALEUR ABSOLUE
            COSI = ABS( COS3PT( PT(1,N0), PT(1,N1), PT(1,N2) ) )
            IF( COSI .LT. COSMIN ) THEN
C              MEILLEUR MINIMUM
               COSMIN = COSI
               NP0 = N0
               NP1 = N1
               NP2 = N2
               IF( COSMIN .LT. 0.2 ) GOTO 90
            ENDIF
 40      CONTINUE
 50   CONTINUE
C
 90   IF( COSMIN .GT. COSEUI ) THEN
         IERR = 1
      ELSE
         IERR = 0
      ENDIF
      END
