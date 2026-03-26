      SUBROUTINE POKISO ( K , NOP , X , COPOTR , NOPOTR , IERR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  AJOUTER LE POINT NOP DE COORDONNEES X(2) SUR LA K-EME ISOVALEUR
C -----
C
C PARAMETRES D ENTREE :
C ---------------------
C K      : NUMERO DE L ISOVALEUR
C X      : COORDONNEES DU POINT NOP
C
C PARAMETRES MODIFIES :
C ----------------------
C NOP    : NUMERO DU DERNIER POINT AJOUTE (AVANT ET APRES AJOUT)
C COPOTR : COORDONNEES DES POINTS D INTERSECTION DES ARETES DES
C          TRIANGLES ET DES ISOVALEURS
C NOPOTR : NOPOTR(1,K) NOMBRE DE POINTS D INTERSECTION AVEC L ISO K
C          NOPOTR(2,K) NO DU 1-ER POINT D INTERSECTION AVEC L ISO K
C          NOPOTR(3,K) NO DU 2-ME POINT D INTERSECTION AVEC L ISO K
C
C PARAMETRE RESULTAT :
C --------------------
C IERR   : 0 SI PAS D ANOMALIE RENCONTREE
C          1 PLUS DE 2 POINTS D INTERSECTION ENTRE L ISOVALEUR K ET
C            LES ARETES DU TRIANGLE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 1990
C23456---------------------------------------------------------------012
      include"./incl/gsmenu.inc"
      REAL              COPOTR(2,*),X(2)
      INTEGER           NOPOTR(3,*)
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
C
C     N NOMBRE DE POINTS D INTERSECTION ISO K ET COTES DU TRIANGLE
      IERR   = 0
      N      = NOPOTR( 1 , K )
      IF( N .LE. 0 ) GOTO 30
C
C     LE NOUVEAU POINT A T IL DEJA ETE STOCKE ?
C     -----------------------------------------
      DO 20 I=1,N
         L = NOPOTR( I + 1 , K )
         DO 10 J=1,2
            IF( X( J ) .NE. COPOTR( J , L ) ) GOTO 20
   10    CONTINUE
C
C        LE POINT EST DEJA RANGE
         GOTO 100
   20 CONTINUE
C
C     SI POSSIBLE LE POINT EST AJOUTE
C     -------------------------------
      IF( N .GE. 2 ) THEN
C        LA PLACE EST TOTALEMENT PRISE
         NBLGRC(NRERR) = 1
         KERR(1) = 'PLUS DE 2 INTERSECTIONS ARETES-ISO'
         CALL LEREUR
         IERR = 1
         GOTO 100
      ENDIF
C
C     ADJONCTION DU POINT
   30 NOP = NOP + 1
      N   = N   + 1
      NOPOTR( 1 ,  K  ) = N
      NOPOTR(N+1,  K  ) = NOP
      COPOTR( 1 , NOP ) = X( 1 )
      COPOTR( 2 , NOP ) = X( 2 )
C
  100 RETURN
      END
