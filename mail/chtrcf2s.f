      SUBROUTINE CHTRCF2S( KTITRE, PTXYZD,
     %                     MXTETR, NOTETR, NUDTETR, N1TEVI, N1TETS,
     %                     MXFACO, LEFACO, NO0FAR,
     %                     NBTRCF, NOTRCF, NBSTIS, NOSTIS,
     %                     NBSTCF, NOSTCF,
     %                     MXTE1S, NOTE1S,
     %                     MXTECF, NBTECF, NOTECF, MODIF, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    EN CAS D'ETOILE DE TETRAEDRES DU CONTOUR FERME
C -----    AYANT 3 SOMMETS DU CF + 2 SOMMETS 1 AU DESSOUS et 1 AU DESSUS
C          LES NBTRCF TRIANGLES PERDUS SONT JOINTS AUX
C          2 SOMMETS POUR FORMER 2 NBTRCF TETRAEDRES et
C          LES NBTRCF TRIANGLES PERDUS SONT RETROUVES
C          DANS LA TETRAEDRISATION

C ENTREES:
C --------
C KTITRE : TITRE D'UN TRACE
C PTXYZD : TABLEAU DES COORDONNEES DES POINTS
C          PAR POINT : X  Y  Z DISTANCE_SOUHAITEE
C MXTETR : NOMBRE DE TETRAEDRES DECLARABLES DANS NOTETR
C NOTETR : LISTE DES TETRAEDRES
C          SOMMET1,    SOMMET2,    SOMMET3,    SOMMET4,
C          TETRAEDRE1, TETRAEDRE2, TETRAEDRE3, TETRAEDRE4
C          DE L'AUTRE COTE DE LA FACE
C NUDTETR: PLUS GRAND NUMERO DANS NOTETR DE TETRAEDRE OCCUPE
C N1TEVI : NO NOTETR DU PREMIER TETRAEDRE VIDE
C N1TETS : N1TETS(NS) NUMERO D'UN TETRAEDRE AYANT POUR SOMMET NS

C MXFACO : NOMBRE MAXIMAL DECLARABLE DE FACES DU CONTOUR
C LEFACO : FACE DU CONTOUR OU INTERFACES ENTRE VOLUMES
C          IL CONTIENT DANS CET ORDRE
C          1:   =0 POUR UNE FACE VIDE
C          123: NO (DANS PTXYZD) DU SOMMET 1, SOMMET 2, SOMMET 3
C          45:  NO (DANS NUVOPA 0 SINON) DU VOLUME1 , VOLUME2 DE LA FACE
C          678: NO (DANS LEFACO) DE LA FACE ADJACENTE PAR L'ARETE 1 2 3
C          9: ATTENTION: UNE ARETE PEUT APPARTENIR A PLUS DE 2 FACES
C             => CHAINAGE CIRCULAIRE DE CES FACES DANS LEFACO
C             LEFACO(9,*) -> FACE SUIVANTE (*=0:VIDE, *<>0:NON VIDE)
C          10: HACHAGE AVEC LA SOMME DES 3 SOMMETS MODULO MXFACO
C              LF = MOD( NOSOFA(1)+NOSOFA(2)+NOSOFA(3) , MXFACO ) + 1
C              NF = LEFACO( 10, LF ) LE NUMERO DE LA 1-ERE FACE DANS LEFACO
C              SI LA FACE NE CONVIENT PAS. PASSAGE A LA SUIVANTE
C              NF = LEFACO( 9, NF )  ...
C          11: >0  NO NOTETR D'UN TETRAEDRE AYANT CETTE FACE,
C              =0  SINON
cccC          12: = NO FACEOC DE 1 A NBFACES D'OpenCascade
C NO0FAR : NUMERO DES 3 SOMMETS DES FACES AJOUTEES AU CF

C NBTRCF : NOMBRE DE FACES DE NOTRCF
C NOTRCF : >0 NUMERO DANS LEFACO DES TRIANGLES PERDUS  DU CF
C          <0 NUMERO DANS NO0FAR DES TRIANGLES AJOUTES AU CF
C NBSTIS : NOMBRE DE SOMMETS ISOLES DU CF
C NOSTIS : NUMERO PTXYZD DES NBSTIS SOMMETS ISOLES
C NBSTCF : NOMBRE DE SOMMETS DES ARETES PERIPHERIQUES DU CF
C NOSTCF : NUMERO PTXYZD DES NBSTCF SOMMETS DES ARETES PERIPHERIQUES
C MXTE1S : NOMBRE D'ENTIERS DU TABLEAU NOTE1S
C NOTE1S : TABLEAU AUXILIAIRE DE MXTE1S ENTIERS

C SORTIES:
C --------
C NBTECF : NOMBRE DE TETRAEDRES DE L'ETOILE AVANT et APRES
C NOTECF : NUMERO NOTETR DES NBTECF TETRAEDRES DE L'ETOILE
C MODIF  : =0 PAS DE MODIFICATION DES NBTECF TETRAEDRES DE NOTECF
C          =1 REMPLACEMENT DES NBTECF TETRAEDRES PAR DES TETRAEDRES
C             DE FACE NOTRCF PERDUE ET DONC RETROUVEE
C IERR   : =0 SI PAS D'ERREUR DETECTEE >0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  Saint PIERRE du PERRAY          Janvier 2020
C2345X7..............................................................012
      CHARACTER*(*)     KTITRE
      DOUBLE PRECISION  PTXYZD(4,*), V, VOLTET
      INTEGER           NOTETR(8,*), LEFACO(11,0:MXFACO), NO0FAR(3,*),
     %                  NOTRCF(NBTRCF), NOSTIS(NBSTIS), NOSTCF(NBSTCF),
     %                  NOTE1S(MXTE1S), NOTECF(MXTECF), N1TETS(*)

      IERR   = 0
      MODIF  = 0
      NBTECF = 0

C     AJOUT DES SOMMETS ISOLES AUX SOMMETS DES TRIANGLES DU CF
      NBSTICF = NBSTIS + NBSTCF
      DO N = 1, NBSTIS
         NOSTCF( NBSTCF + N ) = NOSTIS( N )
      ENDDO

C     LES SOMMETS NON CF DES TETRAEDRES NOTECF SONT STOKES
C     A LA SUITE DES SOMMETS CF ou ISOLES
      NST4 = 0
      NBSTET = NBSTICF
      DO N = 1, NBSTICF

C        LE TRIANGLE N PERDU DU CF
         NSTCF0 = NOSTCF( N )

C        RECHERCHE DES TETRAEDRES DE SOMMET NSTCF0
         CALL TETR1S( NSTCF0, N1TETS, NOTETR,
     %                NBTE1S, MXTE1S, NOTE1S, IERR )
         IF( IERR .NE. 0 ) GOTO 9999

         DO 20 K = 1, NBTE1S

C           LE TETRAEDRE K DE SOMMET NSTCF0
            NTE = NOTE1S( K )

C           LE TETRAEDRE NTE DE SOMMET NSTCF0 DU CF A T IL
C           2 AUTRES SOMMETS DU CF?
            NBSTECF = 0
            DO I = 1, 4

C              LE SOMMET I DU TETRAEDRE NTE EST IL UN SOMMET CF?
               NSTE = NOTETR(I,NTE)

               DO 10 NN = 1, NBSTICF
C                 LE NO DU SOMMET NN DU CF ou ISOLE
                  IF( NSTE .EQ. NOSTCF(NN) ) THEN 
C                    LE SOMMET I DU TETRAEDRE NTE EST UN SOMMET CF
                     NBSTECF = NBSTECF + 1
C                    MARQUAGE NEGATIF DU SOMMET CF DE NTE
                     NOTETR(I,NTE) = -NSTE
                     IF( NBSTECF .EQ. 3 ) THEN
C                       LE TETRAEDRE NTE A 3 SOMMETS CF
C                       QUEL EST LE NO DU SOMMET NON CF DANS NTE?
                        DO I4=1,4
                           IF( NOTETR(I4,NTE) .LT. 0 ) THEN
C                             LE SIGNE REDEVIENT >0
                              NOTETR(I4,NTE) = ABS( NOTETR(I4,NTE) )
                           ELSE
C                             LE 4-EME SOMMET DE NTE NON CF
                              NST4 = NOTETR(I4,NTE)
                           ENDIF
                        ENDDO
                        GOTO 15
                     ENDIF
                  ENDIF
 10            ENDDO

            ENDDO

C           RETOUR AUX VALEURS POSITIVES DES SOMMETS DE NTE
            DO I4=1,4
               IF( NOTETR(I4,NTE) .LT. 0 ) THEN
C                 LE SIGNE REDEVIENT >0
                  NOTETR(I4,NTE) = ABS( NOTETR(I4,NTE) )
               ENDIF
            ENDDO
C           PASSAGE AU TETRAEDRE SUIVANT DE SOMMET NSTCF0
            GOTO 20

C           AJOUT DU TETRAEDRE NTE DANS NOTECF S'IL N'Y EST PAS DEJA
 15         DO M=1,NBTECF
               IF( NTE .EQ. NOTECF(M) ) GOTO 20
            ENDDO
            IF( NBTECF .GE. MXTECF ) THEN
C              SATURATION DU TABLEAU NOTECF
               PRINT*,'chtrcf2s: SATURATION DU TABLEAU NOTECF MXTECF=',
     %                 MXTECF
               IERR = 1
               GOTO 9999
            ENDIF
            NBTECF = NBTECF + 1
            NOTECF( NBTECF ) = NTE

C           AJOUT DU 4-EME SOMMET DE NTE NON CF AU TABLEAU NOSTCF?
            DO I = 1, NBSTET
               IF( NST4 .EQ. NOSTCF(I) ) GOTO 20
            ENDDO
            NBSTET = NBSTET + 1
            NOSTCF( NBSTET ) = NST4

 20      ENDDO

      ENDDO

      IF( NBTECF .LE. 0 ) GOTO 9999

      IF( NBTECF .EQ. 2*NBTRCF .AND. NBSTET .EQ. NBSTICF+2 ) THEN

C        CONFIGURATION DE NBTRCF TRIANGLES PERDUS AVEC
C        UN SOMMET AU DESSUS ET UN SOMMET AU DESSOUS
C        => REMPLACEMENT DES 2*NBTRCF TETRAEDRES NOTECF
C           PAR LES TETRAEDRES TRANGLE CF + 1 SOMMET DESSUS
C                              TRANGLE CF + 1 SOMMET DESSOUS
C        ---------------------------------------------------

C        TRACE DES NBTECF TETRAEDRES ET OPPOSES
      KTITRE='chtrcf2s:         TETRAEDRES AVEC 3 SOMMETS CF ET 2 POINTS
     % DESSOUS DESSUS'
         WRITE( KTITRE(11:17), '(I7)' ) NBTECF
         CALL SANSDBL( KTITRE, L )
         CALL TRFETO4( KTITRE(1:L), PTXYZD, 0, NFETOI,
     %                 NBTRCF, NOTRCF, LEFACO, NO0FAR,
     %                 NBTECF, NOTECF, NOTETR )
ccc         CALL TRFETO13( KTITRE, PTXYZD, NBTECF, NOTECF, NOTETR )

         NBTECR = 0
         DO M=1,2

C           LE SOMMET INFERIEUR ou SUPERIEUR
            NST4 = NOSTCF( NBSTICF+M )

            DO N = 1, NBTRCF

C              NO LEFACO DU TRIANGLE PERDU
               NTRCF = NOTRCF( N )

C              LE TETRAEDRE TRIANGLE N + NST4
               IF( N1TEVI .LE. 0 ) THEN
C                 SATURATION DU TABLEAU NOTETR
                PRINT*,'chtrcf2s: SATURATION DU TABLEAU NOTETR MXTETR=',
     %                  MXTETR
                  IERR = 2
                  GOTO 9999
               ENDIF
               NTE = N1TEVI
               N1TEVI = NOTETR( 5, N1TEVI )

               IF( NBTECR .GE. MXTE1S ) THEN
C                 SATURATION DU TABLEAU NOTE1S
                PRINT*,'chtrcf2s: SATURATION DU TABLEAU NOTE1S MXTE1S=',
     %                  MXTE1S
                  IERR = 3
                  GOTO 9999
               ENDIF
C              LE TETRAEDRE TRIANGLE N + NST4 EST CREE
               NBTECR = NBTECR + 1
               NOTE1S( NBTECR ) = NTE

               IF( NTRCF .GT. 0 ) THEN
                  NOTETR( 1, NTE ) = LEFACO( 1, NTRCF )
                  NOTETR( 2, NTE ) = LEFACO( 2, NTRCF )
                  NOTETR( 3, NTE ) = LEFACO( 3, NTRCF )
               ELSE
                  L = -NTRCF
                  NOTETR( 1, NTE ) = NO0FAR( 1, L )
                  NOTETR( 2, NTE ) = NO0FAR( 2, L )
                  NOTETR( 3, NTE ) = NO0FAR( 3, L )
               ENDIF
               NOTETR( 4, NTE ) = NST4

C              VOLUME SIGNE DU TETRAEDRE NTE
               V = VOLTET( PTXYZD( 1, NOTETR(1,NTE) ),
     %                     PTXYZD( 1, NOTETR(2,NTE) ),
     %                     PTXYZD( 1, NOTETR(3,NTE) ),
     %                     PTXYZD( 1, NOTETR(4,NTE) ) )

               IF( V .LT. 0D0 ) THEN
C                 LE VOLUME EST SUPPOSE POSITIF
C                 PERMUTATION DES SOMMETS 2 et 3 POUR L'OBTENIR
                  I                = NOTETR( 2, NTE )
                  NOTETR( 2, NTE ) = NOTETR( 3, NTE )
                  NOTETR( 3, NTE ) = I
               ENDIF

C              LES TETRAEDRES OPPOSES SONT INCONNUS
               NOTETR( 5, NTE ) = -1
               NOTETR( 6, NTE ) = -1
               NOTETR( 7, NTE ) = -1
               NOTETR( 8, NTE ) = -1

C              MISE A JOUR DE N1TETS
               DO I = 1, 4
                  N1TETS( NOTETR(I,NTE) ) = NTE
               ENDDO

            ENDDO

         ENDDO

C        RECENSEMENT DES TETRAEDRES OPPOSES AUX TETRAEDRES NOTECF
C        NON EUX MEMES DANS NOTECF
         NBTEOP = NBTECR
         DO N = 1, NBTECF
            NTE = NOTECF( N )
            DO 50 I=1,4
               NTEOP = NOTETR( 4+I, NTE )
               IF( NTEOP .GT. 0 ) THEN
                  DO M=1,NBTECF
                     IF( NTEOP .EQ. NOTECF(M) ) GOTO 50
                  ENDDO
C                 NTEOP N'EST PAS DANS NOTECF. IL EST AJOUTE
                  IF( NBTEOP .GE. MXTE1S ) THEN
C                    SATURATION DU TABLEAU NOTE1S
                PRINT*,'chtrcf2s: SATURATION DU TABLEAU NOTE1S MXTE1S=',
     %                  MXTE1S
                     IERR = 4
                     GOTO 9999
                  ENDIF
                  NBTEOP = NBTEOP + 1
                  NOTE1S( NBTEOP ) = NTEOP
               ENDIF
 50         ENDDO
         ENDDO

C        SUPPRESSION DES TETRAEDRES INITIAUX AVEC 3 SOMMETS CF
         DO N = 1, NBTECF
            NTE = NOTECF( N )
C           DESTRUCTION DU TETRAEDRE NTE DU TABLEAU LEFACO
            CALL SUTELEFA( NTE, NOTETR, 1, MXFACO, LEFACO )
C           DESTRUCTION DU TETRAEDRE NTE DU TABLEAU NOTETR
            DO M=1,8
               NOTETR( M, NTE ) = 0
            ENDDO
            NOTETR( 5, NTE ) = N1TEVI
            N1TEVI = NTE
         ENDDO

C        MISE A JOUR DES TETRAEDRES OPPOSES PAR LES FACES
         CALL MJOPTE( NBTEOP, NOTE1S, N1TETS, NOTETR, NUDTETR,
     %                N1TEVI, PTXYZD, NBFANR )

         IF(  NBFANR .GT. 0 ) THEN
            PRINT*,'chtrcf2s:',NBFANR,
     %             ' FACES DE TETRAEDRE OPPOSE INCONNU'
            IERR = 5
            GOTO 9999
         ENDIF

C        LE NO DES TETRAEDRES CREES REMPLACENT
C        LE NO DES TETRAEDRES INITIAUX DANS LE TABLEAU NOTECF DE L'ETOILE
         CALL TRTATA( NOTE1S, NOTECF, NBTECR )
         NBTECF = NBTECR

C        AJOUTER EVENTUELLEMENT LES 4 FACES DES TETRAEDRES CREES
C        DANS LE TABLEAU LEFACO
         DO N = 1, NBTECF
            NTE = NOTECF( N )
            CALL AJTELEFA( NTE, NOTETR, 1, MXFACO, LEFACO )
         ENDDO

         MODIF = 1

C        TRACE DES TETRAEDRES CREES ET DES TETRAEDRES OPPOSES
      KTITRE='chtrcf2s:         TETRAEDRES AVEC 3 SOMMETS CF et 2 POINTS
     % DESSUS DESSOUS'
         WRITE( KTITRE(11:17), '(I7)' ) NBTECF
         CALL SANSDBL( KTITRE, L )
         CALL TRFETO4( KTITRE(1:L), PTXYZD, 0, NFETOI,
     %                 NBTRCF, NOTRCF, LEFACO, NO0FAR,
     %                 NBTECF, NOTECF, NOTETR )
ccc         CALL TRFETO13( KTITRE, PTXYZD, NBTECF, NOTECF, NOTETR )

      ENDIF

 9999 RETURN
      END
