        SUBROUTINE EC2D2TRI( COSMAX, XYZSOM, NT1, NT2, NOTRIA,  MODIF )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CHOISIR LA DIAGONALE DES TRIANGLES NT1 NT2 COPLANAIRES
C -----    QUI MAXIMISE LE MINIMUM DES 2 QUALITES
C          ATTENTION: DANS R3 LES 4 SOMMETS DES 2 TRIANGLES
C                     FORMENT UN TETRAEDRE QUI DISPARAIT OU
C                     APPARAIT SELON LE CHOIX DE LA DIAGONALE

C ENTREE :
C --------
C COSMAX : COSINUS DE L'ANGLE ENTRE LES 2 PLANS AU DESSOUS DUQUEL
C          LES 2 TRIANGLES SONT JUGES NON COPLANAIRES
C XYZSOM : 3 COORDONNEES DES SOMMETS DU MAILLAGE
C NT1    : NUMERO DANS NOTRIA DU 1-ER  TRIANGLE
C NT2    : NUMERO DANS NOTRIA DU 2-EME TRIANGLE

C MODIFIES :
C ----------
C NOTRIA : LISTE CHAINEE DES TRIANGLES
C                         ------- ------- ------- ---
C          PAR TRIANGLE : SOMMET1 SOMMET2 SOMMET3  0
C                         ------- ------- ------- ---
C                         SOMMET EST SON NUMERO DANS XYZSOM
C
C SORTIE :
C --------
C MODIF  : 1 SI LA DIAGONALE DE NT1+NT2 A CHANGE
C          0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC PARIS & VEULETTES SUR MER Avril 2017
C....................................................................012
      REAL     XYZSOM(3,*)
      INTEGER  NOTRIA(4,*), NOSOTR(3)

      IF( NT1 .LE. 0 .OR. NT2 .LE. 0 ) GOTO 9000
      IF( NOTRIA(1,NT1) .LE. 0 ) GOTO 9000
      IF( NOTRIA(1,NT2) .LE. 0 ) GOTO 9000

C     RECHERCHE DE L'ARETE COMMUNE DE NT1 et NT2
C     ------------------------------------------
      DO N1 = 1, 3

         IF( N1 .EQ. 3 ) THEN
            N2 = 1
         ELSE
            N2 = N1 + 1
         ENDIF
         NS1 = NOTRIA( N1, NT1 )
         NS2 = NOTRIA( N2, NT1 )

         DO NN1 = 1, 3
            IF( NN1 .EQ. 3 ) THEN
               NN2 = 1
            ELSE
               NN2 = NN1 + 1
            ENDIF
            NNS1 = NOTRIA( NN1, NT2 )
            NNS2 = NOTRIA( NN2, NT2 )

            IF( ( NS1 .EQ. NNS2 .AND. NS2 .EQ. NNS1 ) .OR.
     %          ( NS1 .EQ. NNS1 .AND. NS2 .EQ. NNS2 ) ) THEN
               GOTO 10
            ENDIF
         ENDDO

      ENDDO

ccc      NT1 et NT2 N'ONT PAS D'ARETE COMMUNE
ccc      PRINT *,'ec2d2tri: TRIANGLES',NT1,NT2,' SANS ARETE COMMUNE'
ccc      PRINT *,'ec2d2tri: NOTRIA(',NT1,')=',(NOTRIA(L,NT1),L=1,4)
ccc      PRINT *,'ec2d2tri: NOTRIA(',NT2,')=',(NOTRIA(L,NT2),L=1,4)

      GOTO 9000

C     LE 3-EME SOMMET DU TRIANGLE NT1
C     -------------------------------
 10   IF( N2 .EQ. 3 ) THEN
         N3 = 1
      ELSE
         N3 = N2 + 1
      ENDIF
      NS3 = NOTRIA( N3, NT1 )
      
C     LE 3-EME SOMMET DU TRIANGLE NT2
C     -------------------------------
      IF( NN2 .EQ. 3 ) THEN
         NN3 = 1
      ELSE
         NN3 = NN2 + 1
      ENDIF
      NS4 = NOTRIA( NN3, NT2 )

C     LE COSINUS DE L'ANGLE DES NORMALES AUX TRIANGLES NT1 NT2
C     EST IL INFERIEUR A COSMAX?
C     --------------------------------------------------------
      CALL COS2TR( XYZSOM(1,NS1), XYZSOM(1,NS2), XYZSOM(1,NS3),
     %             XYZSOM(1,NS2), XYZSOM(1,NS1), XYZSOM(1,NS4),
     %             COS2PL, IERR1, IERR2 )

      IF( IERR1 .NE. 0 .OR. IERR2 .NE. 0 ) THEN
C        SI UN TRIANGLE EST SANS NORMALE
C        ALORS L'ECHANGE EST FORCE INDEPENDAMMENT DE COSMAX
       print*,'ec2d2tri: COS2PL=',COS2PL,' IERR1=',IERR1,' IERR2=',IERR2
         GOTO 20
      ENDIF

      IF( COS2PL .LT. COSMAX ) THEN
         GOTO 9000
      ENDIF

C     LES 2 TRIANGLES FORMENT ILS UN QUADRANGLE 1423 CONVEXE?
C     -------------------------------------------------------
 20   CALL QUADCXR( XYZSOM(1,NS1), XYZSOM(1,NS4),
     %              XYZSOM(1,NS2), XYZSOM(1,NS3),  NONOUI )
      IF( NONOUI .EQ. 0 ) GOTO 9000

C     OUI: QUALITE DES 2 CHOIX POSSIBLES DE TRIANGULATIONS
C     ----------------------------------------------------
C     TRIANGLE 123
      CALL QUATRI( NOTRIA(1,NT1), XYZSOM, Q12 )

C     TRIANGLE 142
      CALL QUATRI( NOTRIA(1,NT2), XYZSOM, Q   )
      Q12 = MIN( Q12, Q )

C     TRIANGLE 143
      NOSOTR(1) = NS1
      NOSOTR(2) = NS4
      NOSOTR(3) = NS3
      CALL QUATRI( NOSOTR, XYZSOM, Q )
C     COMPARAISON DES QUALITES  MIN( Q123, Q142 ) ET Q143
      IF( Q .LE. Q12 ) GOTO 9000

C     TRIANGLE 234
      NOSOTR(1) = NS2
      NOSOTR(2) = NS3
      NOSOTR(3) = NS4
      CALL QUATRI( NOSOTR, XYZSOM, Q )
C     COMPARAISON DES QUALITES  MIN( Q123, Q142 ) ET Q234
      IF( Q .LE. Q12 ) GOTO 9000

C     L'ECHANGE DE L'ARETE NS1-NS2 PAR NS3-NS4 EST POSSIBLE
C     SI L'ARETE NS3-NS4 N'EST PAS DEJA UNE ARETE DE LA TRIANGULATION
C     ---------------------------------------------------------------

cccccccccccc      NS3-NS4 DANS NARFA?






C     ECHANGE EFFECTIF DE LA DIAGONALE DES 2 TRIANGLES NT1 NT2
C     ========================================================
C     ICI NT1 NT2 SONT REMPLACES SUR EUX MEMES
      NOTRIA( 1, NT1 ) = NS1
      NOTRIA( 2, NT1 ) = NS4
      NOTRIA( 3, NT1 ) = NS3
      NOTRIA( 4, NT1 ) = 0

      NOTRIA( 1, NT2 ) = NS2
      NOTRIA( 2, NT2 ) = NS3
      NOTRIA( 3, NT2 ) = NS4
      NOTRIA( 4, NT2 ) = 0

      MODIF = 1
      GOTO 9999


C     PAS D'ECHANGE DES 2 DIAGONALES
C     ==============================
 9000 MODIF = 0

 9999 RETURN
      END
