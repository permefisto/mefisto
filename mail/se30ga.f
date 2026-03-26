      SUBROUTINE SE30GA( NT0,   IMAX,   COSMAXPL,
     %                   NBSOM, XYZSOM, NBTRIA, NOTRIA, NOSTFR, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     POUR AMELIORER LA QUALITE D'UNE TRIANGULATION
C -----     TRAITER LE TRIANGLE NT0 AYANT UN TROP GRAND ANGLE IMAX
C           NT1 LE TRIANGLE OPPOSE A L'ARETE OPPOSEE AU GRAND ANGLE IMAX
C           SI ANGLE 2 PLANS (NT0,NT1) est PETIT ALORS
C            LE SOMMET DE NT0 EST PLACE AU MILIEU DE L'ARETE OPPOSEE
C            LEQUEL EST JOINT AU SOMMET OPPOSE DU TRIANGLE ADJACENT.
C            2 TRIANGLES => 2 TRIANGLES
C           SINON
C            ANGLE 2 PLANS (NT0,NT1) GRAND => CREATION DU SOMMET NS5 MILIEU
C            DE L'ARETE 2 DE NT0 ET DECOUPAGE EN 2 DE NT0 ET NT1
C            2 TRIANGLES => 4 TRIANGLES DONT 2 AJOUTES
C           ET SEULEMENT SI LA QUALITE MINIMALE DES TRIANGLES EST AMELIOREE

C ENTREES:
C --------
C NT0    : NUMERO NOTRIA DU TRIANGLE D'ANGLE IMAX TROP GRAND
C IMAX   : NUMERO DU SOMMET D'ANGLE TROP GRAND
C COSMAXPL:COSINUS DE L'ANGLE ENTRE LES 2 PLANS AU DESSOUS DUQUEL
C          IL N'Y A PAS COPLANEARITE DES 2 TRIANGLES

C NBSOM  : NOMBRE DE SOMMETS DU MAILLAGE
C NBTRIA : NOMBRE DE TRIANGLES DU MAILLAGE

C MODIFIES:
C ---------
C XYZSOM : 3 COORDONNEES DES NBSOM SOMMETS
C NOTRIA : NUMERO DES 3 SOMMETS ET 3 TRIANGLES ADJACENTS PAR LES ARETES
C NOSTFR : 0 OU 1 SELON QUE LE SOMMET EST INTERNE OU FRONTALIER

C SORTIES:
C --------
C NBSOM  : NOMBRE DE SOMMETS   DU MAILLAGE
C NBTRIA : NOMBRE DE TRIANGLES DU MAILLAGE
C IERR   : 0 SI LE PLUS GRAND ANGLE DU TRIANGLE NT0 A ETE REDUIT
C          1 SI IMPOSSIBLE DE REDUIRE LE PLUS GRAND ANGLE DU TRIANGLE NT0
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : A.PERRONNET LABORATOIRE ANALYSE NUMERIQUE PARIS NOVEMBRE 1993
C MODIFS : A.PERRONNET LJLL UPMC & St Pierre du Perray      JANVIER 2016
C2345X7..............................................................012
      COMMON / TRTETR / STOPTE,TRACTE
      LOGICAL           STOPTE,TRACTE,TRACTE0

      INTEGER           NOTRIA(6,*), NOSTFR(*), LITRIA(8), NOSOTR(3)
      REAL              XYZSOM(3,NBSOM), XYZM23(3), COS2FA

      tracte0 = tracte

C     TRACE DU TRIANGLE NT0 ET DES 3 TRIANGLES ADJACENTS
      LITRIA(1) = NT0
      CALL TRTRIAN( 'se30ga', XYZSOM, 6, 1, LITRIA, NOTRIA )

C     CALCUL DE LA QUALITE Q DU TRIANGLE NT0
      CALL QUATRI( NOTRIA(1,NT0), XYZSOM, Q0 )

C     SI 2 ARETES DE NT0 SONT FRONTALIERES
C     ALORS LE TRIANGLE NE PEUT ETRE RESORBE CAR PERTE DE FRONTIERE
      NT1 = NOTRIA(6,NT0)
      DO I=4,6
         NT3 = NOTRIA(I,NT0)
         IF( NT1 .EQ. 0 .AND. NT3 .EQ. 0 ) THEN
            IERR = 1
            RETURN
         ENDIF
         NT1 = NT3
      ENDDO

      IERR  = 0
      NSENS = 0
      Q1    = 0

C     PAR PERMUTATION => 1 SEUL CAS IMAX=1 NUMERO DU SOMMET DE
C                                          TROP GRAND ANGLE
      IF( IMAX .EQ. 2 ) THEN
         I                = NOTRIA( 1, NT0 )
         NOTRIA( 1, NT0 ) = NOTRIA( 2, NT0 )
         NOTRIA( 2, NT0 ) = NOTRIA( 3, NT0 )
         NOTRIA( 3, NT0 ) = I
         I                = NOTRIA( 4, NT0 )
         NOTRIA( 4, NT0 ) = NOTRIA( 5, NT0 )
         NOTRIA( 5, NT0 ) = NOTRIA( 6, NT0 )
         NOTRIA( 6, NT0 ) = I
      ELSE IF( IMAX .EQ. 3 ) THEN
         I                = NOTRIA( 1, NT0 )
         NOTRIA( 1, NT0 ) = NOTRIA( 3, NT0 )
         NOTRIA( 3, NT0 ) = NOTRIA( 2, NT0 )
         NOTRIA( 2, NT0 ) = I
         I                = NOTRIA( 4, NT0 )
         NOTRIA( 4, NT0 ) = NOTRIA( 6, NT0 )
         NOTRIA( 6, NT0 ) = NOTRIA( 5, NT0 )
         NOTRIA( 5, NT0 ) = I
      ENDIF

C     LES 3 SOMMETS DE NT0
      NS1 = NOTRIA( 1, NT0 )
      NS2 = NOTRIA( 2, NT0 )
      NS3 = NOTRIA( 3, NT0 )

C     LE TRIANGLE NT1 OPPOSE PAR L'ARETE 2 OPPOSE AU SOMMET 1
C     ( DE PLUS GRAND ANGLE ) DE NT0 SOIT ENCORE
C     LE TRIANGLE ADJACENT A L'ARETE NOTRIA(2:3,NT0)
      NT1 = NOTRIA( 5, NT0 )

C     LES 3 COORDONNEES DU MILIEU DE L'ARETE NS2-NS3 OPPOSEE
C     AU PLUS GRAND ANGLE DU TRIANGLE NT0 DE SOMMET NS1
      DO I=1,3
         XYZM23(I) = ( XYZSOM(I,NS2) + XYZSOM(I,NS3) ) * 0.5
      ENDDO

C     LE MILIEU DE NS2-NS3 EST IDENTIFIABLE AU SOMMET OPPOSE NS1?
      CALL XYZIDE( XYZM23, XYZSOM(1,NS1), IDENT )
      IF( IDENT .NE. 0 ) THEN
C        OUI: NS1 = MILIEU NS2-NS3
         IF( NT1 .GT. 0 ) THEN

C           NT1 EST LE TRIANGLE OPPOSE A NT0 PAR L'ARETE NS2-NS3
C           -> CAS 1: NT0 SUPPRIME et NT1 DECOUPE EN 2
C           ----------------------------------------------------
            GOTO 22

         ELSE

C           NT0 EST QUASI-PLAT SANS TRIANGLE OPPOSE
C            -> NT0 EST SUPPRIME et SES VOISINS LE PERDENT EN OPPOSE
C           --------------------------------------------------------
            DO I=1,3
               NT1 = NOTRIA( 3+I, NT0 )
               IF( NT1 .GT. 0 ) THEN
                  DO J=1,3
                     IF( NOTRIA( 3+J, NT1 ) .EQ. NT0 ) THEN
                         NOTRIA( 3+J, NT1 ) = 0
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
            DO I=1,6
               NOTRIA( I, NT0 ) = 0
            ENDDO
            GOTO 9999

         ENDIF
      ENDIF

C     PAS D'IDENTIFICATION DU MILIEU DE NS2-NS3 A NS1
C     LES TRIANGLES ADJACENTS AUX ARETES DE NT0
      NT12 = NOTRIA( 4, NT0 )
      NT31 = NOTRIA( 6, NT0 )

      IF( NT1 .LE. 0 ) THEN
C        -> CAS 3
         GOTO 300
      ENDIF


C     NS4 3-EME SOMMET DE NT1 N'APPARTENANT PAS A L'ARETE NS2-NS3 DE NT0
      DO I=1,3
         NS4 = NOTRIA(I,NT1)
         IF( NS4 .NE. NS2 .AND. NS4 .NE. NS3 ) GOTO 15
      ENDDO

C     CALCUL  DE L'ANGLE ENTRE LES TRIANGLES NT0 ET NT1
C     COSINUS DE L'ANGLE ENTRE LES NORMALES AUX 2 FACES NT0 NT1
C     ORIENTEES DANS LE MEME SENS
 15   CALL COS2TR( XYZSOM(1,NS2), XYZSOM(1,NS3), XYZSOM(1,NS1),
     %             XYZSOM(1,NS3), XYZSOM(1,NS2), XYZSOM(1,NS4),
     %             COS2FA, IERR1, IERR2 )
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

      IF( IERR1 .NE. 0 .OR. IERR2 .NE. 0 ) THEN
C        SI UN TRIANGLE EST SANS NORMALE (IL EST PLAT)
C        ALORS L'ECHANGE DES 2 DIAGONALES EST FORCE
         CALL EC2DIA( COSMAXPL, NT0, 2, NOTRIA, NBSOM, XYZSOM,  NT3 )
         IF( NT3 .GT. 0 ) GOTO 9999
      ENDIF

      IF( ABS( COS2FA ) .LT. 0.99 ) GOTO 100


C     -------------------------------------------------------------
C     CAS 1: ANGLE(NT0,NT1) PETIT => LE SOMMET 1 (PLUS GRAND ANGLE)
C            L'ARETE NOTRIA(2:3,NT0) DISPARAIT
C            LE TRIANGLE NT0 DISPARAIT
C            LE TRIANGLE OPPOSE NT1 EST DECOUPE EN 2 -> NT0 et NT1
C            2T -> 2T
C     -------------------------------------------------------------

C     NS1 EST IL FRONTALIER ?
      IF( NOSTFR(NS1) .NE. 0 .OR.
     %   (NOSTFR(NS2) .NE. 0 .AND. NOSTFR(NS3) .NE. 0) ) THEN
C        NS1 EST FRONTALIER
         NOSTFR(NS1) = 1
      ENDIF

C     L'ARETE DANS NT1 DE SOMMET NS2 NS3
 22   DO I=1,3
         IF( I .NE. 3 ) THEN
            I1 = I + 1
         ELSE
            I1 = 1
         ENDIF
         IF( NOTRIA(I ,NT1) .EQ. NS3 .AND.
     %       NOTRIA(I1,NT1) .EQ. NS2 ) THEN
            NSENS = 1
            GOTO 24
         ELSE IF( NOTRIA(I ,NT1) .EQ. NS2 .AND.
     %            NOTRIA(I1,NT1) .EQ. NS3 ) THEN
            NSENS = -1
            GOTO 24
         ENDIF
      ENDDO
C
C     LE 3-EME SOMMET DU TRIANGLE NT1
 24   IF( I1 .NE. 3 ) THEN
         I2 = I1 + 1
      ELSE
         I2 = 1
      ENDIF
      NS4 = NOTRIA(I2,NT1)

C     CALCUL DE LA QUALITE Q1 DU TRIANGLE NT1
      CALL QUATRI( NOTRIA(1,NT1), XYZSOM, Q1 )

C     CALCUL DE LA QUALITE DES 2 TRIANGLES A CREER
      NOSOTR(1) = NS1
      NOSOTR(2) = NS4
      NOSOTR(3) = NS3
      CALL QUATRI( NOSOTR, XYZSOM, Q2 )

      NOSOTR(1) = NS1
      NOSOTR(2) = NS2
      NOSOTR(3) = NS4
      CALL QUATRI( NOSOTR, XYZSOM, Q3 )

      IF( MIN(Q2,Q3) .LE. MIN(Q0,Q1) ) THEN
C        LA QUALITE CREEE DETERIORAIT LA QUALITE INITIALE
         IERR = 1
         GOTO 9999
      ENDIF

      PRINT*,'se30ga: CAS 1  2T->2T NT0=',NT0,' NT1=',NOTRIA(5,NT0),
     %       ' NS1=',NS1,' DEVIENT MILIEU de l''ARETE OPPOSEE'

C     LES TRIANGLES VOISINS
      IF( NSENS .GT. 0 ) THEN

         NT24 = NOTRIA(3+I1,NT1)
         NT43 = NOTRIA(3+I2,NT1)

         NOTRIA(1,NT0) = NS1
         NOTRIA(2,NT0) = NS4
         NOTRIA(3,NT0) = NS3
         NOTRIA(4,NT0) = NT1
         NOTRIA(5,NT0) = NT43
         NOTRIA(6,NT0) = NT31

         NOTRIA(1,NT1) = NS1
         NOTRIA(2,NT1) = NS2
         NOTRIA(3,NT1) = NS4
         NOTRIA(4,NT1) = NT12
         NOTRIA(5,NT1) = NT24
         NOTRIA(6,NT1) = NT0

      ELSE

         NT24 = NOTRIA(3+I2,NT1)
         NT43 = NOTRIA(3+I1,NT1)

         NOTRIA(1,NT0) = NS1
         NOTRIA(2,NT0) = NS3
         NOTRIA(3,NT0) = NS4
         NOTRIA(4,NT0) = NT31
         NOTRIA(5,NT0) = NT43
         NOTRIA(6,NT0) = NT1

         NOTRIA(1,NT1) = NS1
         NOTRIA(2,NT1) = NS4
         NOTRIA(3,NT1) = NS2
         NOTRIA(4,NT1) = NT0
         NOTRIA(5,NT1) = NT24
         NOTRIA(6,NT1) = NT12

      ENDIF

C     MISE A JOUR DANS LES TRIANGLES ADJACENTS DES TRIANGLES MODIFIES
      IF( NT43 .GT. 0 ) THEN
         DO I=4,6
            IF( NOTRIA(I,NT43) .EQ. NT1 ) GOTO 30
         ENDDO
         PRINT*,'PB se30ga: TRIANGLE NT1=',NT1,
     % ' NON ADJACENT AU TRIANGLE NT43=',NT43
         IERR = 3
         GOTO 33
 30      NOTRIA(I,NT43) = NT0
      ENDIF

 33   IF( NT12 .GT. 0 ) THEN
         DO I=4,6
            IF( NOTRIA(I,NT12) .EQ. NT0 ) GOTO 40
         ENDDO
         PRINT*,'PB se30ga: TRIANGLE NT0=',NT0,
     % ' NON ADJACENT AU TRIANGLE NT12=',NT12
         IERR = 4
         GOTO 9000
 40      NOTRIA(I,NT12) = NT1
      ENDIF

C     ESSAI D'ECHANGER LES 2 DIAGONALES
      IF( NT1 .GT. 0 ) THEN
C        LES 2 TRIANGLES NT0 ET NT1 EXISTENT
         CALL EC2DIA( COSMAXPL, NT0, 2, NOTRIA, NBSOM, XYZSOM,  NT3 )
         CALL EC2DIA( COSMAXPL, NT1, 2, NOTRIA, NBSOM, XYZSOM,  NT3 )
      ENDIF

      GOTO 9999


C     --------------------------------------------------------------
C     CAS 2: ANGLE(NT0,NT1) GRAND => CREATION DU SOMMET NS5 MILIEU
C            DE L'ARETE 2 DE NT0 ET DECOUPAGE EN 2 DE NT0 ET NT1
C            2T -> 4T TRIANGLES DE SOMMET COMMUN NS5
C     --------------------------------------------------------------
C     CREATION DU MILIEU NS5 DE L'ARETE NS2-NS3
C     DANS se30.f LA DECLARATION DE XYZSOM ET NSEF A ETE DOUBLEE
C     CE QUI EXPLIQUE LE NON CONTROLE DE DEBORDEMENT DE CES 2 TABLEAUX
 100  NBSOM = NBSOM + 1
      NS5   = NBSOM
      DO I=1,3
         XYZSOM(I,NS5) = XYZM23(I)
      ENDDO

C     NS5 EST IL FRONTALIER ?
      IF( NOSTFR(NS2) .NE. 0 .AND. NOSTFR(NS3) .NE. 0 ) THEN
C        NS5 EST FRONTALIER
         NOSTFR(NS5) = 1
      ENDIF
C
C     LES TRIANGLES ADJACENTS AUX ARETES DE NT0
      NT12 = NOTRIA( 4, NT0 )
      NT31 = NOTRIA( 6, NT0 )

      IF( NT1 .LE. 0 ) GOTO 108

C     L'ARETE DANS NT1 DE SOMMET NS2 NS3
      DO I=1,3
         IF( I .NE. 3 ) THEN
            I1 = I + 1
         ELSE
            I1 = 1
         ENDIF
         IF( NOTRIA(I ,NT1) .EQ. NS3 .AND.
     %       NOTRIA(I1,NT1) .EQ. NS2 ) THEN
            NSENS = 1
            GOTO 104
         ELSE IF( NOTRIA(I ,NT1) .EQ. NS2 .AND.
     %            NOTRIA(I1,NT1) .EQ. NS3 ) THEN
            NSENS = -1
            GOTO 104
         ENDIF
      ENDDO

C     LE 3-EME SOMMET DU TRIANGLE NT1
 104  IF( I1 .NE. 3 ) THEN
         I2 = I1 + 1
      ELSE
         I2 = 1
      ENDIF
      NS4 = NOTRIA(I2,NT1)

C     LES TRIANGLES VOISINS AVEC 2 TRIANGLES AJOUTES
 108  IF( NT1 .GT. 0 ) THEN

C        LE TRIANGLE NT1 OPPOSE EXISTE et EST DECOUPE EN 2 TRIANGLES
C        SI LA QUALITE MINIMALE EST AMELIOREE
C        CALCUL DE LA QUALITE DES TRIANGLES A CREER
         NOSOTR(1) = NS1
         NOSOTR(2) = NS2
         NOSOTR(3) = NS5
         CALL QUATRI( NOSOTR, XYZSOM, Q2 )

         NOSOTR(1) = NS2
         NOSOTR(2) = NS4
         NOSOTR(3) = NS5
         CALL QUATRI( NOSOTR, XYZSOM, Q3 )

         NOSOTR(1) = NS1
         NOSOTR(2) = NS5
         NOSOTR(3) = NS3
         CALL QUATRI( NOSOTR, XYZSOM, Q4 )

         NOSOTR(1) = NS5
         NOSOTR(2) = NS4
         NOSOTR(3) = NS3
         CALL QUATRI( NOSOTR, XYZSOM, Q5 )

         IF( MIN(Q2,Q3,Q4,Q5) .LE. MIN(Q0,Q1) ) THEN
C           LA QUALITE CREEE DETERIORAIT LA QUALITE INITIALE
            IERR = 1
            GOTO 9999
         ENDIF

         PRINT*,'se30ga: CAS 2  2T->4T NT0=',NT0,' NT1=',NT1,
     %       ' NS5=',NS5,' NOUVEAU SOMMET MILIEU de l''ARETE',NS2,NS3

         IF( NSENS .GT. 0 ) THEN
            NT24 = NOTRIA(3+I1,NT1)
            NT43 = NOTRIA(3+I2,NT1)
         ELSE
            NT24 = NOTRIA(3+I2,NT1)
            NT43 = NOTRIA(3+I1,NT1)
         ENDIF

         NBTRIA = NBTRIA + 1
         NT2    = NBTRIA

         NBTRIA = NBTRIA + 1
         NT3    = NBTRIA

         NOTRIA(1,NT0) = NS1
         NOTRIA(2,NT0) = NS2
         NOTRIA(3,NT0) = NS5
         NOTRIA(4,NT0) = NT12
         NOTRIA(5,NT0) = NT1
         NOTRIA(6,NT0) = NT2

         NOTRIA(1,NT1) = NS2
         NOTRIA(2,NT1) = NS4
         NOTRIA(3,NT1) = NS5
         NOTRIA(4,NT1) = NT24
         NOTRIA(5,NT1) = NT3
         NOTRIA(6,NT1) = NT0

         NOTRIA(1,NT2) = NS1
         NOTRIA(2,NT2) = NS5
         NOTRIA(3,NT2) = NS3
         NOTRIA(4,NT2) = NT0
         NOTRIA(5,NT2) = NT3
         NOTRIA(6,NT2) = NT31

         NOTRIA(1,NT3) = NS5
         NOTRIA(2,NT3) = NS4
         NOTRIA(3,NT3) = NS3
         NOTRIA(4,NT3) = NT1
         NOTRIA(5,NT3) = NT43
         NOTRIA(6,NT3) = NT2
      
C        MISE A JOUR DANS LES TRIANGLES ADJACENTS DES TRIANGLES MODIFIES
         IF( NT43 .GT. 0 ) THEN
            DO I=4,6
               IF( NOTRIA(I,NT43) .EQ. NT1 ) GOTO 110
            ENDDO
            PRINT*,'PB se30ga: TRIANGLE NT1=',NT1,
     %             ' NON ADJACENT AU TRIANGLE NT43=',NT43
            IERR = 3
            GOTO 112

 110        NOTRIA(I,NT43) = NT3
         ENDIF
C
 112     IF( NT31 .GT. 0 ) THEN
            DO I=4,6
               IF( NOTRIA(I,NT31) .EQ. NT0 ) GOTO 115
            ENDDO
            PRINT*,'PB se30ga: TRIANGLE NT0=',NT0,
     %             ' NON ADJACENT AU TRIANGLE NT31=',NT31
            IERR = 4
            GOTO 9000

 115        NOTRIA(I,NT31) = NT2
         ENDIF
C
C        ESSAI D'ECHANGER LES 2 DIAGONALES POUR AMELIORER LA QUALITE
         CALL EC2DIA( COSMAXPL, NT0, 1, NOTRIA, NBSOM, XYZSOM,  N )
         CALL EC2DIA( COSMAXPL, NT1, 1, NOTRIA, NBSOM, XYZSOM,  N )
         CALL EC2DIA( COSMAXPL, NT2, 3, NOTRIA, NBSOM, XYZSOM,  N )
         CALL EC2DIA( COSMAXPL, NT3, 2, NOTRIA, NBSOM, XYZSOM,  N )

         GOTO 9000

      ENDIF

C     ---------------------------------------------
C            PAS DE TRIANGLE NT1=0 OPPOSE
C     CAS 3: NT0 EST DECOUPE EN 2 TRIANGLES NT0 NT2
C            SI LA QUALITE MINIMALE EST AMELIOREE
C     ---------------------------------------------
C     CREATION MOMENTANEE DU MILIEU DE L'ARETE NS2-NS3
 300  NBSOM = NBSOM + 1
      NS5   = NBSOM
      DO I=1,3
         XYZSOM(I,NS5) = XYZM23(I)
      ENDDO

C     CALCUL DE LA QUALITE DES 2 TRIANGLES A CREER
      NOSOTR(1) = NS1
      NOSOTR(2) = NS2
      NOSOTR(3) = NS5
      CALL QUATRI( NOSOTR, XYZSOM, Q2 )

      NOSOTR(1) = NS1
      NOSOTR(2) = NS5
      NOSOTR(3) = NS3
      CALL QUATRI( NOSOTR, XYZSOM, Q3 )

      IF( MIN(Q2,Q3) .LE. Q0 ) THEN
C        LA QUALITE CREEE DETERIORAIT LA QUALITE INITIALE
         NBSOM = NBSOM - 1
         IERR = 1
         GOTO 9999
      ENDIF

      PRINT*,'se30ga: CAS 3  1T->2T NT0=',NT0,' NT1=',NT1,
     %       ' NS5=',NS5,' NOUVEAU SOMMET MILIEU de l''ARETE',NS2,NS3

      NT24 = 0
      NT43 = 0
      NT3  = 0

      NBTRIA = NBTRIA + 1
      NT2    = NBTRIA

      NOTRIA(1,NT0) = NS1
      NOTRIA(2,NT0) = NS2
      NOTRIA(3,NT0) = NS5
      NOTRIA(4,NT0) = NT12
      NOTRIA(5,NT0) = 0
      NOTRIA(6,NT0) = NT2

      NOTRIA(1,NT2) = NS1
      NOTRIA(2,NT2) = NS5
      NOTRIA(3,NT2) = NS3
      NOTRIA(4,NT2) = NT0
      NOTRIA(5,NT2) = 0
      NOTRIA(6,NT2) = NT31

C     ESSAI D'ECHANGER LES 2 DIAGONALES POUR AMELIORER LA QUALITE
      CALL EC2DIA( COSMAXPL, NT0, 1, NOTRIA, NBSOM, XYZSOM,  N )
      CALL EC2DIA( COSMAXPL, NT2, 3, NOTRIA, NBSOM, XYZSOM,  N )


C     TRACE DES TRIANGLES MODIFIES ET DE LEURS 3 TRIANGLES ADJACENTS
C     ==============================================================
 9000 LITRIA(1) = NT0
      LITRIA(2) = NT1
      LITRIA(3) = NT12
      LITRIA(4) = NT31
      LITRIA(5) = NT24
      LITRIA(6) = NT43
      N = 6
      IF( NT2 .GT. 0 ) THEN
         LITRIA(7) = NT2
         LITRIA(8) = NT3
         N = 8
      ENDIF

C     COMPRESSION DU TABLEAU LITRIA PAR SUPPRESSION DES ZEROS
      M = 0
      DO K=1,N
         L = LITRIA(K)
         IF( L .GT. 0 ) THEN
            M = M+1
            LITRIA(M) = L
         ENDIF
      ENDDO

C     TRACE DES M TRIANGLES
      CALL TRTRIAN( 'se30ga F', XYZSOM, 6, M, LITRIA, NOTRIA )
C
 9999 tracte = tracte0
      RETURN
      END
