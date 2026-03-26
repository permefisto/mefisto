      SUBROUTINE EX3T1TPL( COSMAX, XYZSOM, NOTRIA, NT, K )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    LE TRIANGLE NT A T IL UN SOMMET K CENTRE D'UN GRAND TRIANGLE
C -----    FORME DE NT ET DE 2 TRIANGLES ADJACENTS S'ENROULANT AUTOUR
C          ET CES 3 TRIANGLES SONT ILS COPLANAIRES?

C ENTREES:
C --------
C COSMAX : COSINUS DE L'ANGLE ENTRE LES 2 PLANS AU DESSOUS DUQUEL
C          2 TRIANGLES SONT JUGES NON COPLANAIRES
C XYZSOM : 3 COORDONNEES DES SOMMETS
C NOTRIA : NUMERO DES 3 SOMMETS ET 3 TRIANGLES ADJACENTS PAR LES ARETES
C NT     : NUMERO DANS NOTRIA DU TRIANGLE A TRAITER

C SORTIE :
C --------
C K      : =0 AUCUN DES 3 SOMMETS DE NT EST LE CENTRE DE 3 TRIANGLES
C          >0 LE SOMMET K DE NT         EST LE CENTRE DE 3 TRIANGLES
C             QUI FORMENT UN GRAND TRIANGLE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC et St PIERRE DU PERRAY Janvier 2016
C2345X7..............................................................012
      INTEGER   NOTRIA(6,*)
      REAL      XYZSOM(3,*), C, COSMAX

      IF( NT .EQ. 0 ) GOTO 30

      DO 20 K=1,3

C        SOMMET K DU TRIANGLE NT
ccc      NS1 = NOTRIA( K, NT )

C        LE TRIANGLE ADJACENT PAR L'ARETE K
         NT12 = NOTRIA( 3+K, NT )
         IF( NT12 .EQ. 0 ) GOTO 20

C        NUMERO K0 DE L'ARETE QUI PRECEDE L'ARETE K DANS LE TRIANGLE NT
         IF( K .EQ. 1 ) THEN
            K0 = 3
         ELSE
            K0 = K-1
         ENDIF

C        LE TRIANGLE ADJACENT PAR L'ARETE K0
         NT31 = NOTRIA( 3+K0, NT )
         IF( NT31 .EQ. 0 ) GOTO 20

         DO I=4,6
            IF( NOTRIA(I,NT12) .EQ. NT31 ) THEN
C              LE SOMMET K DE NT EST LE CENTRE DES 3 TRIANGLES NT NT12 NT13
               GOTO 50
            ENDIF
         ENDDO

 20   ENDDO

C     AUCUN DES 3 SOMMETS DE NT EST LE CENTRE DE 3 TRIANGLES
C     OU LES 3 TRIANGLES NE SONT PAS COPLANAIRES
C     ------------------------------------------------------
 30   K = 0
      RETURN

C     LES 3 TRIANGLES NT NT12 NT31 SONT ILS COPLANAIRES?
C     --------------------------------------------------
C     LES 3 SOMMETS DU TRIANGLE NT
 50   NS1 = NOTRIA( K,  NT )
      IF( K .EQ. 3 ) THEN
         K1 = 1
      ELSE
         K1 = K+1
      ENDIF
      NS2 = NOTRIA( K1, NT )
      NS3 = NOTRIA( K0, NT )

C     NS4 SOMMET DE NT12 DIFFERENT DE NS1 ET NS2
      DO L=1,3
         NS4 = NOTRIA( L, NT12 )
         IF( NS4 .NE. NS1 .AND. NS4 .NE. NS2 ) GOTO 60
      ENDDO

C     CALCUL  DE L'ANGLE ENTRE LES TRIANGLES NT ET NT12
C     COSINUS DE L'ANGLE ENTRE LES NORMALES AUX 2 TRIANGLES
C     ORIENTES DANS LE MEME SENS
 60   CALL COS2TR( XYZSOM(1,NS1),
     %             XYZSOM(1,NS2),
     %             XYZSOM(1,NS3),
     %             XYZSOM(1,NS1),
     %             XYZSOM(1,NS4),
     %             XYZSOM(1,NS2),
     %             C, IERR1, IERR2 )
C     ( 0.95     => 18.2  DEGRES )
C     ( 0.96     => 16.3  DEGRES )
C     ( 0.97     => 14.1  DEGRES )
C     ( 0.98     => 11.5  DEGRES )
C     ( 0.99     =>  8.11 DEGRES )
C     ( 0.9962   =>  5    DEGRES )
C     ( 0.99756  =>  4    DEGRES )
C     ( 0.99863  =>  3    DEGRES )
C     ( 0.999    =>  2.56 DEGRES )
C     ( 0.9999   =>  0.8  DEGRES )

      IF( C .LT. COSMAX ) GOTO 30

C     CALCUL  DE L'ANGLE ENTRE LES TRIANGLES NT ET NT31
      CALL COS2TR( XYZSOM(1,NS1),
     %             XYZSOM(1,NS2),
     %             XYZSOM(1,NS3),
     %             XYZSOM(1,NS1),
     %             XYZSOM(1,NS3),
     %             XYZSOM(1,NS4),
     %             C, IERR1, IERR2 )

      IF( C .LT. COSMAX ) GOTO 30

      RETURN
      END
