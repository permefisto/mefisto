        SUBROUTINE EC2D2TR( COSMAX, XYZSOM, NQT1,   NQT2,  NOSTQT,
     %                      LIBREA, L1ARFA, L2ARFA, NARFA,
     %                      MODIF )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CHOISIR LA DIAGONALE DES TRIANGLES NQT1 NQT2 COPLANAIRES
C -----    QUI MAXIMISE LE MINIMUM DES 2 QUALITES
C          ATTENTION: DANS R3 LES 4 SOMMETS DES 2 TRIANGLES
C                     FORMENT UN TETRAEDRE QUI DISPARAIT OU
C                     APPARAIT SELON LE CHOIX DE LA DIAGONALE

C ENTREE :
C --------
C COSMAX : COSINUS DE L'ANGLE ENTRE LES 2 PLANS AU DESSOUS DUQUEL
C          LES 2 TRIANGLES SONT JUGES NON COPLANAIRES
C XYZSOM : 3 COORDONNEES DES SOMMETS DU MAILLAGE
C NQT1   : NUMERO DANS NOSTQT DU 1-ER  TRIANGLE
C NQT2   : NUMERO DANS NOSTQT DU 2-EME TRIANGLE
C L1ARFA : NOMBRE DE MOTS PAR ARFA DU TABLEAU NARFA
C L2ARFA : NOMBRE DE QTANGLES DU TABLEAU NARFA

C MODIFIES :
C ----------
C NOSTQT : LISTE CHAINEE DES TRIANGLES
C                         ------- ------- ------- ---
C          PAR TRIANGLE : SOMMET1 SOMMET2 SOMMET3  0
C                         ------- ------- ------- ---
C                         SOMMET EST SON NUMERO DANS XYZSOM
C LIBREA : LA 1-ERE ARFA LIBRE A PARTIR DES DERNIERES PAR VALEUR DECROISSANTE
C NARFA  : TABLEAU DES ARETES DES QTANGLES DU MAILLAGE
C          NARFA(1,I)= NO DU 1-ER  SOMMET DE L'ARETE
C          NARFA(2,I)= NO DU 2-EME SOMMET > 1-ER  SOMMET
C          NARFA(3,I)= CHAINAGE HACHAGE SUR LE QTANGLE SUIVANT
C          NARFA(4:L1ARFA,I)= NO NOSTQT DU QTANGLE CONTENANT L'ARETE
C          SI UNE ARETE APPARTIENT A PLUS DE L1ARFA-3 QTANGLES, 
C          LE NUMERO DE QTANGLE en POSITION L1ARFA EST RENDU NEGATIF
C          POUR INDIQUER QUE LA LISTE DES QTANGLES EST INCOMPLETE

C SORTIE :
C --------
C MODIF  : 1 SI LA DIAGONALE DE NQT1+NQT2 A ETE CHANGEE
C          0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET St Pierre du Perray               Avril   2017
C MODIFS: ALAIN PERRONNET Veulettes sur mer                 Octobre 2019
C....................................................................012
      REAL     XYZSOM(3,*)
      INTEGER  NOSTQT(4,*), NARFA(L1ARFA,L2ARFA), NOSOTR(3),
     %         NS(2), NOQTAR(16)

      IF( NQT1 .LE. 0  .OR.  NQT2 .LE. 0 ) GOTO 9000
      IF( NOSTQT(1,NQT1) .LE. 0 ) GOTO 9000
      IF( NOSTQT(1,NQT2) .LE. 0 ) GOTO 9000

C     RECHERCHE DE L'ARETE COMMUNE DE NQT1 et NQT2
C     --------------------------------------------
      DO N1 = 1, 3

         IF( N1 .EQ. 3 ) THEN
            N2 = 1
         ELSE
            N2 = N1 + 1
         ENDIF
         NS1 = NOSTQT( N1, NQT1 )
         NS2 = NOSTQT( N2, NQT1 )

         DO NN1 = 1, 3
            IF( NN1 .EQ. 3 ) THEN
               NN2 = 1
            ELSE
               NN2 = NN1 + 1
            ENDIF
            NNS1 = NOSTQT( NN1, NQT2 )
            NNS2 = NOSTQT( NN2, NQT2 )

            IF( ( NS1 .EQ. NNS2 .AND. NS2 .EQ. NNS1 ) .OR.
     %          ( NS1 .EQ. NNS1 .AND. NS2 .EQ. NNS2 ) ) THEN
               GOTO 10
            ENDIF
         ENDDO

      ENDDO

ccc      NQT1 et NQT2 N'ONT PAS D'ARETE COMMUNE
ccc      PRINT *,'ec2d2tr: TRIANGLES',NQT1,NQT2,' SANS ARETE COMMUNE'
ccc      PRINT *,'ec2d2tr: NOSTQT(',NQT1,')=',(NOSTQT(L,NQT1),L=1,4)
ccc      PRINT *,'ec2d2tr: NOSTQT(',NQT2,')=',(NOSTQT(L,NQT2),L=1,4)

      GOTO 9000

C     VERIFICATION ET LIMITATION
C     LISTE DES NBQTAR QTANGLES DE L'ARETE NS1-NS2
C     --------------------------------------------
 10   CALL LISTQTAR( L1ARFA,   L2ARFA, NARFA,
     %               NS1, NS2, NOARET, NBQTAR, NOQTAR )
      IF( NBQTAR .NE. 2 ) THEN
C        ARETE INEXISTANTE ou O ou >2 QTANGLES POUR CETTE ARETE
C        => PAS DE MODIFICATION
         GOTO 9000
      ENDIF

C     LE 3-EME SOMMET DU TRIANGLE NQT1
C     --------------------------------
      IF( N2 .EQ. 3 ) THEN
         N3 = 1
      ELSE
         N3 = N2 + 1
      ENDIF
      NS3 = NOSTQT( N3, NQT1 )
      
C     LE 3-EME SOMMET DU TRIANGLE NQT2
C     --------------------------------
      IF( NN2 .EQ. 3 ) THEN
         NN3 = 1
      ELSE
         NN3 = NN2 + 1
      ENDIF
      NS4 = NOSTQT( NN3, NQT2 )

C     LE COSINUS DE L'ANGLE DES NORMALES AUX TRIANGLES NQT1 NQT2
C     EST IL INFERIEUR A COSMAX?
C     ----------------------------------------------------------
      CALL COS2TR( XYZSOM(1,NS1), XYZSOM(1,NS2), XYZSOM(1,NS3),
     %             XYZSOM(1,NS2), XYZSOM(1,NS1), XYZSOM(1,NS4),
     %             COS2PL, IERR1, IERR2 )

      IF( IERR1 .NE. 0 .OR. IERR2 .NE. 0 ) THEN
C        SI UN TRIANGLE EST SANS NORMALE
C        ALORS L'ECHANGE EST FORCE INDEPENDAMMENT DE COSMAX
        print*,'ec2d2tr: COS2PL=',COS2PL,' IERR1=',IERR1,' IERR2=',IERR2
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
      CALL QUATRI( NOSTQT(1,NQT1), XYZSOM, Q12 )

C     TRIANGLE 142
      CALL QUATRI( NOSTQT(1,NQT2), XYZSOM, Q   )
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
C     SI L'ARETE NS3-NS4 N'EST PAS DEJA UNE ARETE DE LA QTRIANGULATION
C     ----------------------------------------------------------------
C     NS3-NS4 DANS NARFA?
      IF( NS3 .LE. NS4 ) THEN
         NS(1) = NS3
         NS(2) = NS4
      ELSE
         NS(1) = NS4
         NS(2) = NS3
      ENDIF
      CALL HACHAR( 2, NS,  L1ARFA, L2ARFA, NARFA, 3, NAR )
      IF( NAR .GT. 0 ) THEN
C        L'ARETE NS EST DEJA DANS LE TABLEAU NARFA DES ARETES
         GOTO 9000
      ENDIF


C     ECHANGE EFFECTIF DE LA DIAGONALE DES 2 TRIANGLES NQT1 NQT2
C     ==========================================================
C     SUPPRESSION DES 3 ARETES DE NQT1 DU TABLEAU NARFA
      CALL SU3ANARF( L1ARFA, L2ARFA, NARFA, NQT1, NOSTQT, IERR )
      IF( IERR .NE. 0 ) GOTO 9000

C     SUPPRESSION DES 3 ARETES DE NQT2 DU TABLEAU NARFA
      CALL SU3ANARF( L1ARFA, L2ARFA, NARFA, NQT2, NOSTQT, IERR )
      IF( IERR .NE. 0 ) GOTO 9000

C     ICI NQT1 NQT2 SONT REMPLACES SUR EUX MEMES
      NOSTQT( 1, NQT1 ) = NS1
      NOSTQT( 2, NQT1 ) = NS4
      NOSTQT( 3, NQT1 ) = NS3
      NOSTQT( 4, NQT1 ) = 0

      NOSTQT( 1, NQT2 ) = NS2
      NOSTQT( 2, NQT2 ) = NS3
      NOSTQT( 3, NQT2 ) = NS4
      NOSTQT( 4, NQT2 ) = 0

C     AJOUT DES 3 ARETES DE NQT1 DANS LE TABLEAU NARFA
      CALL AJ3ANARF( LIBREA, L1ARFA, L2ARFA, NARFA,
     %               NQT1,   NOSTQT, IERR )

C     AJOUT DES 3 ARETES DE NQT2 DANS LE TABLEAU NARFA
      CALL AJ3ANARF( LIBREA, L1ARFA, L2ARFA, NARFA,
     %               NQT2,   NOSTQT, IERR )

C     ECHANGE DE NS1-NS2 PAR NS3-NS4 EFFECTUE
      MODIF = 1
      GOTO 9999


C     PAS D'ECHANGE DES 2 DIAGONALES
C     ==============================
 9000 MODIF = 0


 9999 RETURN
      END
