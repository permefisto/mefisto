      SUBROUTINE SE30PA( ITRQM,  IMIN,   NBTRQM, NOTRQM,
     %                     NBTRIA, NOTRIA, NOSTFR, XYZSOM, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    POUR AMELIORER LA QUALITE D'UNE TRIANGULATION
C -----    TRAITER LE TRIANGLE NT AYANT UN TROP PETIT ANGLE
C          SON ARETE OPPOSEE EST REDUITE EN UN POINT
C          SOIT MILIEU DE L'ARETE, SOIT L'UN DES 2 SOMMETS DE L'ARETE
C          ET L'AUTRE SOMMET EST ALORS SUPPRIME
C
C ATTENTION: UNE ARETE DOIT APPARTENIR AU PLUS A 2 TRIANGLES
C
C ENTREES:
C --------
C ITRQM  : NUMERO DANS NOTRQM DU TRIANGLE D'ANGLE IMIN TROP PETIT
C IMIN   : NUMERO DU SOMMET D'ANGLE TROP PETIT
C NBTRQM : NOMBRE DE TRIANGLES A TRAITER
C
C MODIFIES :
C ----------
C NBTRIA : NOMBRE DE TRIANGLES DE LA TRIANGULATION AVANT ET APRES
C NOTRIA : NUMERO DES 3 SOMMETS ET 3 TRIANGLES ADJACENTS PAR LES ARETES
C XYZSOM : 3 COORDONNEES DES SOMMETS
C NOSTFR : NUMERO 0 OU 1 SELON QUE LE SOMMET APPARTIENT A UNE ARETE
C          APPARTENANT ELLE-MEME A UN SEUL TRIANGLE
C NOTRQM : NUMERO DES NBTRQM TRIANGLES A TRAITER
C IERR   : 0 SI TRIANGLE TRAITE, 1 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : A.PERRONNET LABORATOIRE ANALYSE NUMERIQUE PARIS NOVEMBRE 1993
C MODIFS : A.PERRONNET LJLL UPMC et St PIERRE du PERRAY     Octobre 2015
C2345X7..............................................................012
      COMMON / TRTETR / STOPTE,TRACTE
      LOGICAL           STOPTE,TRACTE,TRACTE0

      INTEGER   NOTRIA(6,NBTRIA),
     %          NOSTFR(*),
     %          NOTRQM(NBTRQM),
     %          LITRIA(8)
      REAL      XYZSOM(3,*)

      tracte0 = tracte
ccc      tracte = .true.

C     TRACE DU TRIANGLE NT DE TRES PETIT ANGLE ET 3 TRIANGLES ADJACENTS
      IERR = 0
      NT = NOTRQM( ITRQM )

      IF( NOTRIA(1,NT) .LE. 0 ) THEN
         PRINT*, 'se30pa: ANOMALIE dans NOTRIA(',NT,')=',
     %           (NOTRIA(I1,NT),I1=1,6)
         IERR = 1
         GOTO 9900
      ENDIF

      LITRIA(1) = NT
      CALL TRTRIAN( 'se30pa', XYZSOM, 6, 1, LITRIA, NOTRIA )

      I1 = 0
      I2 = 0
C
C     TRAITEMENT DE L'ANGLE IMIN DU TRIANGLE NT
C     -----------------------------------------
C     PAR PERMUTATION => 1 SEUL CAS IMIN=1
      IF( IMIN .EQ. 2 ) THEN
         I               = NOTRIA( 1, NT )
         NOTRIA( 1, NT ) = NOTRIA( 2, NT )
         NOTRIA( 2, NT ) = NOTRIA( 3, NT )
         NOTRIA( 3, NT ) = I
         I               = NOTRIA( 4, NT )
         NOTRIA( 4, NT ) = NOTRIA( 5, NT )
         NOTRIA( 5, NT ) = NOTRIA( 6, NT )
         NOTRIA( 6, NT ) = I
      ELSE IF( IMIN .EQ. 3 ) THEN
         I               = NOTRIA( 1, NT )
         NOTRIA( 1, NT ) = NOTRIA( 3, NT )
         NOTRIA( 3, NT ) = NOTRIA( 2, NT )
         NOTRIA( 2, NT ) = I
         I               = NOTRIA( 4, NT )
         NOTRIA( 4, NT ) = NOTRIA( 6, NT )
         NOTRIA( 6, NT ) = NOTRIA( 5, NT )
         NOTRIA( 5, NT ) = I
      ENDIF

C     NUMERO DU SOMMET 1 DU TRIANGLE NT D'ANGLE MIN
      NS1  = NOTRIA( 1, NT )

C     L'ARETE NOTRIA(2:3,NT)
      NS2  = NOTRIA( 2, NT )
      NS3  = NOTRIA( 3, NT )

C     LES 3 TRIANGLES ADJACENTS PAR LES COTES
      NT12 = NOTRIA( 4, NT )
      NT23 = NOTRIA( 5, NT )
      NT31 = NOTRIA( 6, NT )

C     LE SOMMET NS1=NOTRIA(1,NT) EST IL LE CENTRE DE 3 TRIANGLES?
C     DANS CE CAS INTERDIT => LE TRIANGLE NT N'EST PAS MODIFIE
C     -----------------------------------------------------------
      IF( NT12 .EQ. 0 .AND. NT31 .EQ. 0 ) GOTO 10

      IF( NT12 .GT. 0 ) THEN
         DO I=4,6
            IF( NOTRIA(I,NT12) .EQ. NT31 ) THEN
C              CAS INTERDIT : SOMMET CENTRE DE 3 TRIANGLES
               NOTRQM( ITRQM ) = 0
               IERR = 1
               GOTO 9900
            ENDIF
         ENDDO
      ENDIF

C     RECHERCHE DU SOMMET OPPOSE A NT DANS NT23
      IF( NT23 .GT. 0 ) THEN
C        RECHERCHE DE L'ARETE NS2-NS3 DANS NT23
         DO I=1,3
            IF( I .EQ. 3 ) THEN
               I1 = 1
            ELSE
               I1 = I + 1
            ENDIF
            IF( (NOTRIA(I ,NT23) .EQ. NS2 .AND.
     %           NOTRIA(I1,NT23) .EQ. NS3) .OR.
     %          (NOTRIA(I ,NT23) .EQ. NS3 .AND.
     %           NOTRIA(I1,NT23) .EQ. NS2) ) GOTO 4
         ENDDO
 4       IF( I1 .EQ. 3 ) THEN
            I2 = 1
         ELSE
            I2 = I1 + 1
         ENDIF
         NT3 = NOTRIA(3+I1,NT23)
         NT2 = NOTRIA(3+I2,NT23)
         DO I=4,6
            IF( NOTRIA(I,NT2) .EQ. NT3 ) THEN
C              CAS INTERDIT : SOMMET CENTRE DE 3 TRIANGLES
               NOTRQM( ITRQM ) = 0
               DO I1=ITRQM+1, NBTRQM
                  IF( NOTRQM(I1) .EQ. NT23 ) THEN
                     NOTRQM(I1) = 0
                     IERR = 1
                     GOTO 9900
                  ENDIF
               ENDDO
               IERR = 1
               GOTO 9900
            ENDIF
         ENDDO
      ENDIF
C
C     ICI: AUCUN DES SOMMETS OPPOSES N'EST LE CENTRE DE 3 TRIANGLES

C     RECHERCHE DU MEILLEUR POINT NSD: MILIEU DE NS2-NS3 ou NS2 ou NS3?
C     -------------------------------- --------------------------------
C     AUTOUR DU SOMMET NS2 CALCUL DE L'ANGLE ENTRE 2 TRIANGLES ADJACENTS
C     PAR UNE ARETE DE SOMMET NS2, DEPART DE NT
 10   COSMI2 = 2.0
C     1-ERE ARETE D'ADJACENCE NS1-NS2 DES 2 TRIANGLES NTT1 NTT2
      NBA2 = 0
      NSS1 = NS1
      NSS3 = NS3
      NTT1 = NT
      NTT2 = NT12

 20   IF( NTT2 .EQ. 0 ) GOTO 29

C     RECHERCHE DU SOMMET NSS4 DU TRIANGLE NTT2 NON NSS1 NS2
      NBA2 = NBA2 + 1
      IF( NBA2 .GT. 32 ) THEN
C        BOUCLE -> ABANDON
         IERR = 1
         GOTO 9900
      ENDIF

      DO K=1,3
         NSS4 = NOTRIA(K,NTT2)
         IF( NSS4 .LE. 0 ) THEN
            IERR = 1
            GOTO 9900
         ENDIF
         IF( NSS4 .NE. NSS1 .AND. NSS4 .NE. NS2 ) GOTO 23
      ENDDO

C     CALCUL DU COSINUS DE L'ANGLE ENTRE LES NORMALES AUX 2 TRIANGLES
C     ORIENTES DANS LE MEME SENS
 23   CALL COS2TR( XYZSOM(1,NSS1), XYZSOM(1,NS2),  XYZSOM(1,NSS3),
     %             XYZSOM(1,NSS1), XYZSOM(1,NSS4), XYZSOM(1,NS2),
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

      COSMI2 = MIN( COSMI2, C )

      IF( NTT2 .EQ. NT ) GOTO 29

C     RECHERCHE DU TRIANGLE NTT3 ADJACENT A NTT2 PAR L'ARETE NS2-NSS4
      NTT3 = NOTRIA(3+K,NTT2)
      IF( NTT3 .EQ. 0 ) GOTO 29

C     LE TRIANGLE NTT3 A T IL NS2 POUR SOMMET?
 25   NTT30 = NTT3
      DO L=1,3
         IF( NOTRIA(L,NTT3) .EQ. NS2 ) GOTO 27
      ENDDO
C     LE TRIANGLE NTT3 EST SANS SOMMET NS2
      IF( K .EQ. 1 ) THEN
         K0 = 3
      ELSE
         K0 = K-1
      ENDIF
      NTT3 = NOTRIA(3+K0,NTT2)
      IF( NTT3 .GT. 0 ) THEN
         IF( NTT3 .NE. NTT30 ) THEN
            K = K0
            GOTO 25
         ENDIF
      ENDIF
      GOTO 29

C     LE SOMMET L DU TRIANGLE NTT3 EST NS2. OU EST NSS4?
C     ARETE SUIVANT L'ARETE L DANS NTT3
 27   IF( L .EQ. 3 ) THEN
         L1 = 1
      ELSE
         L1 = L + 1
      ENDIF

C     ARETE PRECEDANT L'ARETE L DANS NTT3
      IF( L .EQ. 1 ) THEN
         L0 = 3
      ELSE
         L0 = L - 1
      ENDIF

C     TRIANGLE ADJACENT PAR ARETE NS2-NSS4?
      NTT1 = NTT2
      NTT2 = NTT3
      NSS3 = NSS1
      NSS1 = NSS4
      IF( NOTRIA(L1,NTT3) .EQ. NSS4 ) THEN
C        ARETE NS2-NSS4 EST L'ARETE L DE NTT3
         NSS4 = NOTRIA(L0,NTT3)
      ELSE
C        ARETE NS2-NSS4 EST L'ARETE L DE NTT3
         NSS4 = NOTRIA(L1,NTT3)
      ENDIF
      GOTO 20

C     AUTOUR DU SOMMET NS3 CALCUL DE L'ANGLE ENTRE 2 TRIANGLES ADJACENTS
C     PAR UNE ARETE DE SOMMET NS3
C     ------------------------------------------------------------------
 29   COSMI3 = 2.0
C     1-ERE ARETE D'ADJACENCE NS1-NS3 DES 2 TRIANGLES NTT1 NTT2
      NBA3  = 0
      NSS1 = NS1
      NSS3 = NS2
      NTT1 = NT
      NTT2 = NT31

 30   IF( NTT2 .LE. 0 ) GOTO 39

C     RECHERCHE DU SOMMET NSS4 DU TRIANGLE NTT2 NON NSS1 NS3
      NBA3 = NBA3 + 1
      IF( NBA3 .GT. 32 ) THEN
C        BOUCLE -> ABANDON
         IERR = 1
         GOTO 9900
      ENDIF

      DO K=1,3
         NSS4 = NOTRIA(K,NTT2)
         IF( NSS4 .LE. 0 ) THEN
            IERR = 1
            GOTO 9900
         ENDIF
         IF( NSS4 .NE. NSS1 .AND. NSS4 .NE. NS3 ) GOTO 33
      ENDDO
C     CALCUL DU COSINUS DE L'ANGLE ENTRE LES NORMALES AUX 2 TRIANGLES
C     NTT1 NTT2 ORIENTES DANS LE MEME SENS
 33   IF( NSS3 .EQ. NSS4 ) THEN
         COSMI2 = 1.
         COSMI3 = 1.
         GOTO 39
      ENDIF
      CALL COS2TR( XYZSOM(1,NSS1), XYZSOM(1,NS3),  XYZSOM(1,NSS3),
     %             XYZSOM(1,NSS1), XYZSOM(1,NSS4), XYZSOM(1,NS3),
     %             C, IERR1, IERR2 )
      COSMI3 = MIN( COSMI3, C )

      IF( NTT2 .EQ. NT ) GOTO 39

C     RECHERCHE DU TRIANGLE NTT3 ADJACENT A NTT2 PAR L'ARETE NS3-NSS4
      NTT3 = NOTRIA(3+K,NTT2)
      IF( NTT3 .LE. 0 ) GOTO 39

C     LE TRIANGLE NTT3 A T IL NS3 POUR SOMMET?
 35   NTT30 = NTT3
      DO L=1,3
         IF( NOTRIA(L,NTT3) .EQ. NS3 ) GOTO 37
      ENDDO
C     LE TRIANGLE NTT3 EST SANS SOMMET NS3
      IF( K .EQ. 1 ) THEN
         K0 = 3
      ELSE
         K0 = K-1
      ENDIF
      NTT3 = NOTRIA(3+K0,NTT2)
      IF( NTT3 .GT. 0 ) THEN
         IF( NTT3 .NE. NTT30 ) THEN
            GOTO 35
         ENDIF
      ENDIF
      GOTO 39

C     LE SOMMET L DU TRIANGLE NTT3 EST NS3. OU EST NSS4?
C     ARETE SUIVANT L'ARETE L DANS NTT3
 37   IF( L .EQ. 3 ) THEN
         L1 = 1
      ELSE
         L1 = L + 1
      ENDIF

C     ARETE PRECEDANT L'ARETE L DANS NTT3
      IF( L .EQ. 1 ) THEN
         L0 = 3
      ELSE
         L0 = L - 1
      ENDIF

C     TRIANGLE ADJACENT PAR ARETE NS3-NSS4?
      NTT1 = NTT2
      NTT2 = NTT3
      NSS3 = NSS1
      NSS1 = NSS4
      IF( NOTRIA(L1,NTT3) .EQ. NSS4 ) THEN
C        ARETE NS3-NSS4 EST L'ARETE L DE NTT3
         NSS4 = NOTRIA(L0,NTT3)
      ELSE
C        ARETE NS3-NSS4 EST L'ARETE L DE NTT3
         NSS4 = NOTRIA(L1,NTT3)
      ENDIF
      GOTO 30


C     BILAN SUR LES ANGLES AUTOUR DE NS2 et NS3
C     -----------------------------------------
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

 39   IF( COSMI2 .GE. 0.99 .AND. COSMI3 .GE. 0.99 ) THEN

C        LE SOMMET NSD=MIN(NS2,NS3) DEVIENT LE MILIEU
C        DE L'ARETE NOTRIA(2:3,NT)
C        LE SOMMET DOUBLE ET LE SOMMET PERDU
         MILIEU = 1
         IF( NS2 .LT. NS3 ) THEN
            NSD = NS2
            NSP = NS3
         ELSE
            NSD = NS3
            NSP = NS2
         ENDIF

      ELSE

         IF( COSMI2 .GE. COSMI3 ) THEN
C           NS2 GLISSE JUSQU'A NS3 QUI NE BOUGE PAS
            NSD = NS3
            NSP = NS2
            MILIEU = 3
         ELSE
C           NS3 GLISSE JUSQU'A NS2 QUI NE BOUGE PAS
            NSD = NS2
            NSP = NS3
            MILIEU = 2
         ENDIF

      ENDIF

C     MISE A JOUR DES TRIANGLES ADJACENTS AU TRIANGLE NT
      NT12 = NOTRIA(4,NT)
      NT31 = NOTRIA(6,NT)
      PRINT*,'se30pa: NT modifie=',NT,' NT23 detruit=',NT23,
     %       ' NoStDouble=',NSD,' NoStPerdu=',NSP,' Identifie a',NSD,
     %       ' COSMI2=',COSMI2,' COSMI3=',COSMI3,' MILIEU=',MILIEU

C     RECHERCHE DE L'ARETE NS1-NS2 DANS NT12
      IF( NT12 .GT. 0 ) THEN
         DO I=4,6
            IF( NOTRIA(I,NT12) .EQ. NT ) GOTO 65
         ENDDO
      ENDIF

C     RECHERCHE DE L'ARETE NS3-NS1 DANS NT31
 65   IF( NT31 .GT. 0 ) THEN
         DO J=4,6
            IF( NOTRIA(J,NT31) .EQ. NT ) GOTO 69
         ENDDO
      ENDIF

 69   IF( NT12 .GT. 0 ) NOTRIA(I,NT12) = NT31
      IF( NT31 .GT. 0 ) NOTRIA(J,NT31) = NT12


C     TRAITEMENT DU TRIANGLE OPPOSE A L'ANGLE IMIN DU TRIANGLE NT
C     -----------------------------------------------------------
      IF( NT23 .GT. 0 ) THEN

C        LES TRIANGLES OPPOSES AUX ARETES 31 ET 12 DE NT23
C        RECHERCHE DE L'ARETE NS2-NS3 DANS NT23
         DO I=1,3
            IF( I .EQ. 3 ) THEN
               I1 = 1
            ELSE
               I1 = I + 1
            ENDIF
            IF( (NOTRIA(I ,NT23) .EQ. NS2 .AND.
     %           NOTRIA(I1,NT23) .EQ. NS3) .OR.
     %          (NOTRIA(I ,NT23) .EQ. NS3 .AND.
     %           NOTRIA(I1,NT23) .EQ. NS2) ) GOTO 81
         ENDDO
C        L'ARETE I DE NT23 EST NS2-NS3
 81      IF( I1 .EQ. 3 ) THEN
            I2 = 1
         ELSE
            I2 = I1 + 1
         ENDIF
         NT31 = NOTRIA(3+I1,NT23)
         NT12 = NOTRIA(3+I2,NT23)

C        RECHERCHE DE L'ARETE NS1-NS2 DANS NT12
         IF( NT12 .GT. 0 ) THEN
            DO I=4,6
               IF( NOTRIA(I,NT12) .EQ. NT23 ) GOTO 84
            ENDDO
         ENDIF
C        RECHERCHE DE L'ARETE NS3-NS1 DANS NT31
 84      IF( NT31 .GT. 0 ) THEN
            DO J=4,6
               IF( NOTRIA(J,NT31) .EQ. NT23 ) GOTO 86
            ENDDO
         ENDIF

 86      IF( NT12 .GT. 0 ) NOTRIA(I,NT12) = NT31
         IF( NT31 .GT. 0 ) NOTRIA(J,NT31) = NT12

C        LE TRIANGLE NT23 EST DETRUIT
         NOTRIA(1,NT23) = 0

      ENDIF

C     LE TRIANGLE NT EST DETRUIT
      NOTRIA(1,NT) = 0

C     LES COORDONNEES DU SOMMET NSD SUR LA FRONTIERE OU NON
      IF( NOSTFR(NS2).NE.0 .AND. NOSTFR(NS3).EQ.0 ) NS3 = NS2
      IF( NOSTFR(NS3).NE.0 .AND. NOSTFR(NS2).EQ.0 ) NS2 = NS3

C     COORDONNEES DU MEILLEUR POINT NSD: MILIEU DE NS2-NS3 ou NS2 ou NS3?
C     ----------------------------------
      IF( MILIEU .EQ. 1 ) THEN
C        NSD EST LE MILIEU DE NS2-NS3 ET NSP EST PERDU
         DO I=1,3
            XYZSOM(I,NSD) = ( XYZSOM(I,NS2) + XYZSOM(I,NS3) ) * 0.5
         ENDDO
      ENDIF

C     POSITION DE NSD PAR RAPPORT A LA FRONTIERE
      IF( NOSTFR(NS2).NE.0 .OR. NOSTFR(NS3).NE.0 ) NOSTFR(NSD)=1

C     LE SOMMET NSP DEVIENT NSD DANS TOUS LES TRIANGLES
      DO J = 1, NBTRIA
         IF( NOTRIA(1,J) .GT. 0 ) THEN
C           LE TRIANGLE EXISTE
            DO I = 1, 3
               IF( NOTRIA(I,J) .EQ. NSP ) NOTRIA(I,J)=NSD
            ENDDO
         ENDIF
      ENDDO

C     LISTE DES TRIANGLES MODIFIES ET DE LEURS 3 TRIANGLES ADJACENTS
      LITRIA(1) = NT12
      LITRIA(2) = NT31
      NBT = 2
      IF( NT12 .GT. 0 ) THEN
         LITRIA(3) = NOTRIA(4,NT12)
         LITRIA(4) = NOTRIA(5,NT12)
         LITRIA(5) = NOTRIA(6,NT12)
         NBT = 5
      ENDIF

      IF( NT31 .GT. 0 ) THEN
         LITRIA(NBT+1) = NOTRIA(4,NT31)
         LITRIA(NBT+2) = NOTRIA(5,NT31)
         LITRIA(NBT+3) = NOTRIA(6,NT31)
         NBT = NBT + 3
      ENDIF

C     TRACE DES TRIANGLES MODIFIES ET DE LEURS 3 TRIANGLES ADJACENTS
      CALL TRTRIAN( 'se30pa F', XYZSOM, 6, NBT, LITRIA, NOTRIA )

 9900 IF( IERR .EQ. 0 ) THEN

      ENDIF
      tracte = tracte0
      RETURN
      END
