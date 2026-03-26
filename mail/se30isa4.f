      SUBROUTINE SE30ISA4( NBSOM,  XYZSOM, NOSOM1, NOSOM2,
     %                     MXQTAR, L1ARET, L2ARET, MNARET,
     %                     NBTRIA, NOTRIA, NBSTID )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  IDENTIFICATION DES SOMMETS PROCHES I.E. des ARETES TRES PETITES
C -----  I.E. TRAITEMENT DU CAS D'UN ANGLE TRES PETIT DU TRIANGLE
C        LES 2 SOMMETS DE LA TRES PETITE ARETE SONT IDENTIFIES ET
C        LES TRIANGLES AYANT CETTE ARETE SONT DETRUITS

C ENTREE :
C --------
C NOSOM1 : NUMERO XYZSOM DU PREMIER SOMMET A TENTER D'IDENTIFIER
C NOSOM2 : NUMERO XYZSOM DU DERNIER SOMMET A TENTER D'IDENTIFIER
C MXQTAR : NOMBRE MAXIMAL DE QT PAR ARETE

C MODIFIES:
C ---------
C L1ARET : NOMBRE DE MOTS PAR ARETE DU TABLEAU LARETE (3+MXQTAR)
C L2ARET : NOMBRE DE ARETES         DU TABLEAU LARETE
C MNARET : ADRESSE MCN DU TABLEAU LARETE DES ARETES DU MAILLAGE
C NBSOM  : NOMBRE DE SOMMETS   DU MAILLAGE AVANT et APRES
C XYZSOM : 3 COORDONNEES DES NBSOM SOMMETS AVANT et APRES
C NBTRIA : NOMBRE DE TRIANGLES DU MAILLAGE AVANT et APRES
C NOTRIA : NUMERO DES (3 SOMMETS ET 0) DES NBTRIA TRIANGLES

C SORTIE :
C --------
C NBSTID : NOMBRE DE SOMMETS IDENTIFIES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : A.PERRONNET  Saint Pierre du Perray                 Mars 2020
C2345X7..............................................................012
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      COMMON / TRTETR / STOPTE,TRACTE
      LOGICAL           STOPTE,TRACTE,TRACTE0
      INTEGER           NOTRIA(4,*)
      REAL              XYZSOM(3,NBSOM)

      TRACTE0 = TRACTE
      NBSTID  = 0
      NBSOM0  = NBSOM

C     IDENTIFICATION DES SOMMETS PROCHES I.E. des ARETES TRES PETITES
C     I.E. TRAITEMENT DU CAS D'UN ANGLE TRES PETIT DU TRIANGLE
C     ===============================================================
C     NOQT1S EST UN TABLEAU AUXILIAIRE
      MXQT1S = NBTRIA
      CALL TNMCDC( 'ENTIER', MXQT1S, MNQT1S )
      IF( MNQT1S .LE. 0 ) GOTO 9999

C     NEWS EST UN TABLEAU AUXILIAIRE
      CALL TNMCDC( 'ENTIER', 1+2*NBSOM0, MNNEWS )
      IF( MNNEWS .LE. 0 ) GOTO 9999

C     BOUCLE SUR LES SOMMETS NoSt DE LA TRIANGULATION
      DO 10 NoSt = NOSOM2, NOSOM1, -1

C        LISTE des TRIANGLES DE SOMMET NoSt
         CALL LIEF1ST( NoSt,   4,      NBTRIA,      NOTRIA,
     %                 MXQT1S, NBQT1S, MCN(MNQT1S), IERR )
         IF( NBQT1S .LE. 0 ) GOTO 10

C        LONGUEUR MOYENNE DES ARETES DES NBQT1S TRIANGLES DE SOMMET NoSt
C        RECHERCHE DE L'ARETE DE SOMMET NoSt DE LONGUEUR MINIMALE
         QMIN    = 2
         NTRMIN  = 0
         ARPEMIN = 1E28
         DO K = 1, NBQT1S

C           ADRESSE DU TRIANGLE K DE SOMMET NoSt
            NTR = MCN( MNQT1S-1+K )

C           QUALITE DU TRIANGLE
            CALL QUATRI( NOTRIA(1,NTR), XYZSOM, Q )
            QMIN = MIN( QMIN, Q )

            PERIME  = 0.
            ARETMIN = 1E28
            NS1     = NOTRIA(3,NTR)
            DO L=1,3

               NS2 = NOTRIA(L,NTR)
               ARETE = ( XYZSOM(1,NS2) - XYZSOM(1,NS1) ) **2
     %               + ( XYZSOM(2,NS2) - XYZSOM(2,NS1) ) **2
     %               + ( XYZSOM(3,NS2) - XYZSOM(3,NS1) ) **2
               ARETE = SQRT( ARETE )

C              PERIMETRE DU TRIANGLE
               PERIME = PERIME + ARETE

               IF( NS1 .NE. NoSt .AND. NS2 .NE. NoSt ) THEN
C                 ARETE SANS LE SOMMET NoSt
                  GOTO 5
               ENDIF

               IF( NS1 .EQ. NoSt .AND. NS2 .EQ. NoSt ) THEN
C                 NS1 ou NS2 DEJA IDENTIFIE A NoSt
                  GOTO 5
               ENDIF

               IF( ARETE .LT. ARETMIN ) THEN
C                 L'ARETE MIN DE SOMMET NoSt DU TRIANGLE 
                  ARETMIN = ARETE
               ENDIF

 5             NS1 = NS2

            ENDDO

C           RECHERCHE DU MIN DU RAPPORT ARETE Min/PERIMETRE du TRIANGLE
            ARPE = ARETMIN / PERIME
            IF( ARPE .LT. ARPEMIN ) THEN
               ARPEMIN = ARPE
C              NUMERO DU TRIANGLE DU MIN des RAPPORTS ARETE/PERIMETRE
               NTRMIN  = NTR
            ENDIF

         ENDDO

ccc      IF( ARPEMIN .LE. 0.015 ) THEN    IDENTIFIE PAS ASSEZ
ccc      IF( ARPEMIN .LE. 0.01  ) THEN    IDENTIFIE PAS ASSEZ
ccc      IF( ARPEMIN .LE. 0.03  ) THEN    IDENTIFIE TROP
ccc      IF( ARPEMIN .LE. 0.026 ) THEN

         IF( ARPEMIN .LE. 0.028 ) THEN

C           LE TRIANGLE NTRMIN DE SOMMET NoSt A UNE PETITE ARETE
C           NoSt-NSMIN DE SOMMET NSMIN A IDENTIFIER
C           ----------------------------------------------------
            ARETMIN = 1E28
            NS1     = NOTRIA(3,NTRMIN)
            DO L=1,3

               NS2 = NOTRIA(L,NTRMIN)
               IF( NS1 .NE. NoSt .AND. NS2 .NE. NoSt ) THEN
C                 ARETE SANS LE SOMMET NoSt
                  GOTO 7
               ENDIF

               ARETE = ( XYZSOM(1,NS2) - XYZSOM(1,NS1) ) **2
     %               + ( XYZSOM(2,NS2) - XYZSOM(2,NS1) ) **2
     %               + ( XYZSOM(3,NS2) - XYZSOM(3,NS1) ) **2
               ARETE = SQRT( ARETE )

               IF( ARETE .LT. ARETMIN ) THEN
C                 L'ARETE MIN DU TRIANGLE ISSUE DE NoSt
                  ARETMIN = ARETE
                  IF( NS1 .EQ. NoSt ) THEN
                     NSMIN = NS2
                  ELSE
                     NSMIN = NS1
                  ENDIF
               ENDIF

 7             NS1 = NS2

            ENDDO

            IF( NoSt .EQ. NSMIN ) THEN
C              NoSt DEJA IDENTIFIE
               GOTO 10
            ENDIF

C           IDENTIFICATION DU SOMMET NoSt -> NSMIN DANS LES TRIANGLES
C           ---------------------------------------------------------
            NBSTID = NBSTID + 1
            PRINT*
            PRINT*,'se30isa4 0: le SOMMET',NoSt,' DEVIENT le SOMMET',
     %             NSMIN,' car ARETE Min/PERIMETRE du TRIANGLE=',
     %             ARPEMIN,' de Qualite MIN=',QMIN

C           TRACE DES TRIANGLES DE SOMMET NoSt ET DE LEURS VOISINS
            tracte = .true.
            CALL TRTRIAR( 'se30isa4 0', XYZSOM,
     %                     MXQTAR, L1ARET, L2ARET, MNARET,
     %                     NBTRIA, NOTRIA,
     %                     MXQT1S,NBQT1S,NBTRPVV,NBTRPVVV,MCN(MNQT1S) )

            DO K=1,3
               XYZSOM( K, NoSt ) = XYZSOM( K, NSMIN )
            ENDDO

            DO K = 1, NBQT1S

C              LE K-EME TRIANGLE K DE SOMMET NoSt
               NTR = MCN( MNQT1S-1+K )
               DO L=1,3
                  IF( NOTRIA( L, NTR ) .EQ. NoSt ) THEN
C                     LE SOMMET NoSt DEVIENT NSMIN
                      NOTRIA( L, NTR ) = NSMIN
                  ENDIF
               ENDDO

C              COMBIEN DE FOIS LE TRIANGLE NTR A T IL LE SOMMET NSMIN?
               NB = 0
               DO L=1,3
                  IF( NOTRIA( L, NTR ) .EQ. NSMIN ) THEN
                     NB = NB + 1
                  ENDIF
               ENDDO

               IF( NB .GE. 2 ) THEN
C                 LE TRIANGLE NTR DEGENERE EST SUPPRIME
                  DO L=1,4
                     NOTRIA( L, NTR ) = 0
                  ENDDO
                  MCN( MNQT1S-1+K ) = 0
               ENDIF

            ENDDO

C           COMPRESSION DU TABLEAU QT1S ELARGI AUX VOISINS
            NBQT1S = NBTRPVVV
            NB = 0
            DO K = 1, NBQT1S
C              LE K-EME TRIANGLE K DE SOMMET NoSt
               NTR = MCN( MNQT1S-1+K )
               IF( NTR .GT. 0 ) THEN
                  NB = NB + 1
                  MCN( MNQT1S-1+NB ) = NTR
               ENDIF
            ENDDO
            NBQT1S = NB

C           TRACE DES TRIANGLES DE FAIBLE QUALITE ET DE LEURS VOISINS
            CALL TRTRIAR( 'se30isa4 1', XYZSOM,
     %                     MXQTAR, L1ARET, L2ARET, MNARET,
     %                     NBTRIA, NOTRIA,
     %                     MXQT1S,NBQT1S,NBTRPVV,NBTRPVVV,MCN(MNQT1S) )

C           RECHERCHE ET SUPPRESSION DES TRIANGLES AYANT 2 MEMES SOMMETS
C           NSMIN APRES L'IDENTIFICATION NoSt -> NSMIN
C           MISE A JOUR DU TABLEAU XYZSOM ET NOSOEF EN RENUMEROTANT
C           LES SOMMETS ET ELIMINANT LES EF DESACTIVES
            CALL MAJXYZNSE( NBSOM, XYZSOM, MCN(MNNEWS),
     %                          4, NBTRIA, NOTRIA )

C           VERIFICATION BILAN PAR ARETE DU NOMBRE DE QTANGLES ADJACENTS
C           CONSTRUCTION DU TABLEAU DES ARETES DE LA QTRIANGULATION
C           DESTRUCTION DU TABLEAU DES ARETES DU MAILLAGE
            IF(MNARET.GT.0) CALL TNMCDS('ENTIER', L1ARET*L2ARET, MNARET)
            CALL GEARSU( 4,      NBTRIA, NOTRIA, MXQTAR,
     %                   L1ARET, L2ARET, MNARET, IERR )
C           ARETE(1,I)= NO DU 1-ER  SOMMET DE L'ARETE
C           ARETE(2,I)= NO DU 2-EME SOMMET > 1-ER  SOMMET
C           ARETE(3,I)= CHAINAGE HACHAGE SUR ARETE SUIVANTE
C           ARETE(4:L1ARET,I)= NUMERO DU QTANGLE CONTENANT CETTE ARETE

            IF( NSMIN .GT. NoSt ) THEN
C              DANS LA MISE A JOUR NSMIN A NoSt PLACE AVANT et
C              DISPARAISSANT SON NO DIMINUE DE UN
               NSMIN = NSMIN - 1
            ENDIF

C           LISTE des TRIANGLES QT1S DE SOMMET NSMIN
            CALL LIEF1ST( NSMIN,  4,      NBTRIA,      NOTRIA,
     %                    MXQT1S, NBQT1S, MCN(MNQT1S), IERR )
            IF( NBQT1S .LE. 0 ) GOTO 10

C           TRACE DES TRIANGLES DE FAIBLE QUALITE ET DE LEURS VOISINS
            CALL TRTRIAR( 'se30isa4 2', XYZSOM,
     %                     MXQTAR, L1ARET, L2ARET, MNARET,
     %                     NBTRIA, NOTRIA,
     %                     MXQT1S,NBQT1S,NBTRPVV,NBTRPVVV,MCN(MNQT1S) )
            tracte = tracte0

         ENDIF

 10   ENDDO

      IF( MNNEWS .GT. 0 ) CALL TNMCDS( 'ENTIER', 1+2*NBSOM0, MNNEWS )
      IF( MNQT1S .GT. 0 ) CALL TNMCDS( 'ENTIER', MXQT1S,     MNQT1S )

      PRINT*,'se30isa4:',NBSTID,' SOMMETS IDENTIFIES pour ELIMINER les P
     %ETITES ARETES'

 9999 RETURN
      END
