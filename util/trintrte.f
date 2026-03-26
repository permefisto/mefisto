      SUBROUTINE TRINTRTE( STR, STE, LINTER, XYZINT )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACE L'INTERSECTION D'UN TRIANGLE ET D'UN TETRAEDRE
C -----    VIDE SI AUCUNE ARETE DU TRIANGLE INTERSECTE LES FACES DU TETRAEDRE
C          +  AUCUNE ARETE DU TETRAEDRE INTERSECTE LE TRIANGLE
C          +  AUCUN DES 3 SOMMETS DU TRIANGLE EST INTERNE AU TETRAEDRE
C ENTREES:
C --------
C STR    : 3 COORDONNEES DES 3 SOMMETS DU TRIANGLE
C STE    : 3 COORDONNEES DES 4 SOMMETS DU TETRAEDRE
C LINTER : 0 PAS D'INTERSECTION
C          1 SI UNE ARETE DU TRIANGLE INTERSECTE UNE FACE  DU TETRAEDRE
C          2 SI LE TRIANGLE EST INTERSECTE PAR   UNE ARETE DU TETRAEDRE
C          3 SI UN DES SOMMETS DU TRIANGLE EST STRICTEMENT INTERNE
C            AU TETRAEDRE
C XYZINT : XYZ DU POINT D'INTERSECTION SAUF SI LINTER=0
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Saint Pierre du Perray           Janvier 2020
C23456...............................................................012
      include"./incl/trvari.inc"
      include"./incl/xyzext.inc"
      COMMON / TRTETR / STOPTE, TRACTE
      LOGICAL           STOPTE, TRACTE
C     TRACTE = FAUX ==> PAS DE TRACE DES TETRAEDRES
C              VRAI ==> TRACE DES TETRAEDRES DE L'ETOILE

      DOUBLE PRECISION  STR(3,3), STE(3,4), XYZINT(3)
      REAL              XYZ1(3), XYZ2(3)
      CHARACTER*88      KTITRE
      INTEGER           NOSOARTE( 2, 6 )
      DATA              NOSOARTE/ 1,2, 2,3, 3,1, 1,4, 2,4, 3,4 /

      IF( .NOT. TRACTE ) RETURN

C     MINIMUM ET MAXIMUM DES COORDONNEES DES FACES PERDUES
C     ====================================================
      R = RINFO( 'GRAND' )
      DO K=1,3
         COOEXT(K,1) =  R
         COOEXT(K,2) = -R
      ENDDO

C     XYZ MIN MAX DES 3 SOMMETS DU TRIANGLE
      DO N = 1, 3
         DO K=1,3
            R = REAL( STR(K,N) )
            IF( R .LT. COOEXT(K,1) )  COOEXT(K,1)=R
            IF( R .GT. COOEXT(K,2) )  COOEXT(K,2)=R
         ENDDO
      ENDDO

C     XYZ MIN MAX DES 4 SOMMETS DU TETRAEDRE
      DO N = 1, 4
         DO K=1,3
            R = REAL( STE(K,N) )
            IF( R .LT. COOEXT(K,1) )  COOEXT(K,1)=R
            IF( R .GT. COOEXT(K,2) )  COOEXT(K,2)=R
         ENDDO
      ENDDO

      CALL VISEE0

C     TRACE EN MODE 3 BOUTONS POUR DEPLACEMENT ROTATIONS ZOOM
      NORBITE = 0
      LORBITE = 1

      IF( LORBITE .EQ. 0 ) GOTO 20

C     INITIALISATION DE L'ORBITE
C     ==========================
      CALL ORBITE0( NOTYEV )
      GOTO 20

C     TRACE SELON L'ORBITE OU ZOOM OU TRANSLATION ACTIFS
C     ==================================================
 10   CALL ORBITE1( NOTYEV )
      IF( NOTYEV .EQ. 0 ) GOTO 9000
C
C     TRACE DES AXES 3D
 20   CALL TRAXE3

C     LE TRACE DES 6 ARETES DU TETRAEDRE
      DO N = 1, 6

         NS1 = NOSOARTE( 1, N )
         XYZ1(1) = REAL( STE( 1, NS1 ) )
         XYZ1(2) = REAL( STE( 2, NS1 ) )
         XYZ1(3) = REAL( STE( 3, NS1 ) )
  
         NS2 = NOSOARTE( 2, N )
         XYZ2(1) = REAL( STE( 1, NS2 ) )
         XYZ2(2) = REAL( STE( 2, NS2 ) )
         XYZ2(3) = REAL( STE( 3, NS2 ) )

         CALL TRAIT3D( NCROUG, XYZ1, XYZ2 )

      ENDDO

C     LE TRACE DES 3 ARETES DU TRIANGLE
      DO N = 1, 3

C        NOSOARTR = NOSOARTE POUR LES 3 PREMIERES ARETES
         NS1 = NOSOARTE( 1, N )
         XYZ1(1) = REAL( STR( 1, NS1 ) )
         XYZ1(2) = REAL( STR( 2, NS1 ) )
         XYZ1(3) = REAL( STR( 3, NS1 ) )
  
         NS2 = NOSOARTE( 2, N )
         XYZ2(1) = REAL( STR( 1, NS2 ) )
         XYZ2(2) = REAL( STR( 2, NS2 ) )
         XYZ2(3) = REAL( STR( 3, NS2 ) )

         CALL TRAIT3D( NCBLEU, XYZ1, XYZ2 )

      ENDDO

C     TRACE EVENTUEL DU POINT D'INTERSECTION
      IF( LINTER .GT. 0 ) THEN
         XYZ1(1) = REAL( XYZINT(1) )
         XYZ1(2) = REAL( XYZINT(2) )
         XYZ1(3) = REAL( XYZINT(3) )
         CALL SYMBOLE3D( NCNOIR, XYZ1, '*' )
      ENDIF

C     TRACE DU TITRE
      KTITRE='INTERSECTION TRIANGLE-TETRAEDRE LINTER=    '
      WRITE(KTITRE(40:40),'(I1)') LINTER

      IF( LINTER .EQ. 0 ) THEN
         KTITRE = KTITRE(1:40) // '  PAS d''INTERSECTION'
      ENDIF

      IF( LINTER .EQ. 1 ) THEN
         KTITRE = KTITRE(1:40) // '  ARETE TRIANGLE -> FACE TETRAEDRE'
      ENDIF

      IF( LINTER .EQ. 2 ) THEN
         KTITRE = KTITRE(1:40) // '  ARETE TETRAEDRE -> TRIANGLE'
      ENDIF

      IF( LINTER .EQ. 3 ) THEN
         KTITRE = KTITRE(1:40) //
     %            '  UN SOMMET TRIANGLE INTERNE au TETRAEDRE'
      ENDIF

      CALL TRFINS( KTITRE )

C     REPRISE DE L'ORBITE
C     ===================
      IF( LORBITE .NE. 0 ) GOTO 10

 9000 RETURN
      END
