      SUBROUTINE MJN1TEVI( NOMSP,  MXTETR, NOTETR, IVOLTE, NVOLTE,
     %                     N1TEVI, NUDTETR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT:     METTRE A JOUR N1TEVI POINTEUR SUR LE 1-ER TETRAEDRE VIDE
C -----
C ENTREES:
C --------
C NOMSP  : NOM DU SOUS-PROGRAMME D'APPEL
C MXTETR : NOMBRE DE TETRAEDRES DECLARABLES DANS NOTETR
C NOTETR : LISTE DES TETRAEDRES
C          SOMMET1,    SOMMET2,    SOMMET3,    SOMMET4,
C          TETRAEDRE1, TETRAEDRE2, TETRAEDRE3, TETRAEDRE4
C          DE L'AUTRE COTE DE LA FACE
C          1: 123      2: 234      3: 341      4: 412
C IVOLTE : 0 PAS       DE TABLEAU NVOLTE A L'APPEL
C          1 EXISTENCE DU TABLEAU NVOLTE A L'APPEL
C NVOLTE : >0 NUMERO DU VOLUME DE 1 A NBVOPA DE CHAQUE TETRAEDRE ACTIF
C          -1 SI TETRAEDRE INACTIF
C SORTIES:
C --------
C N1TEVI : NO DU PREMIER TETRAEDRE VIDE  DANS NOTETR
C NUDTETR: NO DU DERNIER TETRAEDRE ACTIF DANS NOTETR
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & VEULETTES SUR MER    Octobre 2014
C2345X7..............................................................012
      CHARACTER*(*)  NOMSP
      INTEGER        NOTETR(8,MXTETR), NVOLTE(MXTETR)

C     LE TETRAEDRE SUIVANT A UN INDICE >0
      PRINT*
      PRINT*, NOMSP,': Debut mjn1tevi avec N1TEVI=',N1TEVI,' et',
     %        MXTETR,' tetraedres declares'

C     VERIFIER SI TOUT TETRAEDRE VIDE A SON 1-ER NUMERO NUL
      CALL VETEVIDE( N1TEVI, MXTETR, NOTETR, IERR )

      NBTVID = 0
      N1TEVI = 0
      NUDTETR = 0
      DO NTE = MXTETR, 1, -1

         IF( NOTETR(5,NTE) .LT. 0 ) THEN
            PRINT *,'mjn1tevi: NOTETR(',NTE,')=',(NOTETR(K,NTE),K=1,8),
     %              ' ???...'
         ENDIF

         IF( NOTETR(1,NTE) .LE. 0 ) THEN

C           TETRAEDRE VIDE
            NOTETR(1,NTE) = 0

C           MISE A JOUR DE SON SUIVANT
            NOTETR(5,NTE) = N1TEVI
            N1TEVI = NTE

C           NUMERO DE VOLUME INCONNU
            IF( IVOLTE .NE. 0 ) NVOLTE(NTE) = -1

            NBTVID = NBTVID + 1

         ELSE

C           NO DU DERNIER TETRAEDRE ACTIF DANS NOTETR
            IF( NUDTETR .EQ. 0 ) NUDTETR = NTE

         ENDIF

      ENDDO

      NBTETR = MXTETR - NBTVID
      PRINT *, NOMSP,': Fin  mjn1tevi avec N1TEVI=',N1TEVI,
     %        ' NBTVID=',NBTVID,' tetraedres VIDES,',
     %        ' NBTETR=',NBTETR,' tetraedres OCCUPES pour',
     %          MXTETR,' tetraedres declares,',
     %          NUDTETR,' DERNIER tetraedre ACTIF,',
     %          N1TEVI,' No du PREMIER TETRAEDRE VIDE'
      PRINT *

      RETURN
      END
