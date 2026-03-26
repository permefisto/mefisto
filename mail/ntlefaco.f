      SUBROUTINE NTLEFACO( MXSOMM, NBSOMM, PTXYZD,
     %                     MXTETR, NOTETR, N1TETS,
     %                     N1FASC, MXFACO, LEFACO,
     %                     MXFAPE, NBFAPE, NOFAPE, N1TEVI, NUDTETR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  MISE A JOUR DU NO NOTETR DE TETRAEDRE CONTENANT UNE FACE LEFACO
C -----  POUR DETERMINER LES FACES LEFACO FRONTIERE PERDUES DANS LA
C        TETRAEDRISATION NOTETR ACTUELLE

C ENTREES:
C --------
C MXSOMM : INDICE MAXIMAL DU TABLEAU N1FASC et PTXYZD et N1TETS
C MXTETR : NOMBRE MAXIMAL DE TETRAEDRES DU TABLEAU NOTETR
C NOTETR : LISTE DES TETRAEDRES
C          SOMMET1,    SOMMET2,    SOMMET3,    SOMMET4,
C          TETRAEDRE1, TETRAEDRE2, TETRAEDRE3, TETRAEDRE4
C          DE L'AUTRE COTE DE LA FACE
C          1: 123      2: 234      3: 341      4: 412
C          NOTETR(5,N) = TETRAEDRE VIDE SUIVANT SI NOTETR(1,NT)=0
C MXFACO : NOMBRE MAXIMAL DE FACES DU TABLEAU LEFACO
C MXSOMM : INDICE MAXIMAL DU TABLEAU N1FASC
C PTXYZD : TABLEAU DES COORDONNEES DES POINTS
C          PAR POINT : X  Y  Z DISTANCE_SOUHAITEE
C MXFAPE : NOMBRE MAXIMAL DE FACES PERDUES DE LEFACO DU TABLEAU NOFAPE

C MODIFIE:
C --------
C N1TETS : N1TETS(I) NUMERO NOTETR D'UN TETRAEDRE AYANT POUR SOMMET I
C LEFACO : FACE DU CONTOUR OU INTERFACES ENTRE VOLUMES
C          IL CONTIENT DANS CET ORDRE
C          1:   =0 POUR UNE FACE VIDE
C          123: NO (DANS PTXYZD) DU SOMMET 1 < SOMMET 2 < SOMMET 3
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
C              =0  SINON C-A-D AUCUN TETRAEDRE CONTIENT CETTE FACE
C              <0  SINON C-A-D AUCUN TETRAEDRE CONTIENT CETTE FACE
C                  MAIS EN PLUS ELLE EST MARQUEE TRAITEE
CCCC       LEFACO(12,*) NO DE FACE OC

C SORTIES:
C --------
C NBSOMM : PLUS GRAND NUMERO PTXYZD D'UN SOMMET DE NOTETR
C N1FASC : N1FASC(I) NUMERO LEFACO D'UN TRIANGLE DE SOMMET I
C NBFAPE : NOMBRE DE FACES LEFACO APPARTENANT A AUCUN TETRAEDRE NOTETR
C NOFAPE : >0  NUMERO LEFACO DES FACES DANS AUCUN TETRAEDRE NOTETR
C          <0 -NUMERO LEFACO DES FACES DANS AUCUN TETRAEDRE NOTETR
C              MAIS EN PLUS ELLE EST MARQUEE ABANDONNEE MOMENTANEMENT
C N1TEVI : NUMERO DU PREMIER TETRAEDRE VIDE
C NUDTETR: NUMERO du DERNIER TETRAEDRE OCCUPE dans NOTETR
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC et St PIERRE du PERRAY    Mars 2016
C MODIFS : ALAIN PERRONNET              St PIERRE du PERRAY Fevrier 2018
C MODIFS : ALAIN PERRONNET              St PIERRE du PERRAY Juillet 2020
C2345X7..............................................................012
      include"./incl/langue.inc"
      COMMON / TRTETR / STOPTE, TRACTE
      LOGICAL           STOPTE, TRACTE, TRACTE0
      DOUBLE PRECISION  PTXYZD(4,MXSOMM)
      INTEGER           NOTETR(8,MXTETR), N1TETS(MXSOMM),
     %                  LEFACO(11,0:MXFACO),
     %                  N1FASC(MXSOMM), NOFAPE(MXFAPE)
      INTEGER           NOSOTR(3), NOSOFATE(1:3,1:4)
C                       NO LOCAL DES 3 SOMMETS DES 4 FACES DU TETRAEDRE
      DATA              NOSOFATE/ 1,3,2, 2,3,4, 3,1,4, 4,1,2 /
C     => LA NORMALE A LA FACE EST ORIENTEE VERS L'EXTERIEUR DU TETRAEDRE

      TRACTE0 = TRACTE

C     REMISE A ZERO DU NUMERO NOTETR DE TETRAEDRE DES FACES LEFACO
      DO K = 1, MXFACO
         LEFACO( 11, K ) = 0
      ENDDO

C     REMISE A ZERO DE N1FASC et N1TETS
      NBSOM0 = NBSOMM
      NBSOMM = 0
      DO K = 1, MXSOMM
         N1FASC( K ) = 0
         N1TETS( K ) = 0
      ENDDO

C     BOUCLE SUR TOUS LES TETRAEDRES ACTUELS DE NOTETR
C     ------------------------------------------------
      NBPBOP = 0
      DO 10 NTE = 1, MXTETR

C        VERIFICATION DE LA COHERENCE DES TETRAEDRES OPPOSES AUX 4 FACES DE NTE
         IF( NOTETR(1,NTE) .LE. 0 .AND. NOTETR(8,NTE) .NE. 0 ) THEN
C           MISE A ZERO COMPLETE DU TETRAEDRE NTE
C           LE NO DU TETRAEDRE VIDE SUIVANT EST MIS A JOUR PLUS LOIN
            DO K=1,8
               NOTETR(K,NTE) = 0
            ENDDO
            GOTO 10
         ENDIF

         DO NFTE=1,4

C           LE TETRAEDRE OPPOSE A LA FACE NFTE DU TETRAEDRE NTE
            NTEOP = NOTETR( 4+NFTE, NTE )
            IF( NTEOP .GT. 0 .AND. NOTETR(1,NTEOP) .GT. 0 ) THEN

C              LES 3 SOMMETS DE LA FACE NFTE DE NTE
               NOSOTR(1) = NOTETR( NOSOFATE(1,NFTE), NTE )
               NOSOTR(2) = NOTETR( NOSOFATE(2,NFTE), NTE )
               NOSOTR(3) = NOTETR( NOSOFATE(3,NFTE), NTE )
               CALL TRI3NO( NOSOTR, NOSOTR )

C              LA FACE NFTE DE NTE EST ELLE VRAIMENT UNE FACE DE NTEOP?
               CALL NO1F1T( NOSOTR, NOTETR(1,NTEOP), NFOP )
               IF( NFOP .LE. 0 ) THEN
C                 NON: LA FACE NFTE DE NTE N'EST PAS UNE FACE DE NTEOP
                  NBPBOP = NBPBOP + 1
                  PRINT*
                  PRINT*,'ntlefaco: Probleme: la FACE',NFTE,
     %                   ' du TETRAEDRE',NTE,
     %                   ' NON COMMUNE au TETRAEDRE OPPOSE',NTEOP
                  PRINT*,'ntlefaco: NOTETR(',NTE,')=',
     %                   (NOTETR(K,NTE),K=1,8)
                  PRINT*,'ntlefaco: NOTETR(',NTEOP,')=',
     %                   (NOTETR(K,NTEOP),K=1,8)

C                 LE TETRAEDRE OPPOSE A LA FACE NFTE DE NTE EST DONC INCONNU
                  NOTETR( 4+NFTE, NTE ) = -1

                  PRINT*
               ELSE
C                 OUI: LE TETRAEDRE OPPOSE A LA FACE NFOP DE NTEOP
C                      DOIT ETRE NTE POUR ETRE CORRECT
                  NTEOP2 = NOTETR(4+NFOP,NTEOP)
                  IF( NTEOP2 .NE. NTE ) THEN

C                    OPPOSITION INCORRECTE
                     PRINT*
                     PRINT*,'ntlefaco: Face',NOSOTR,' du TETRAEDRE',NTE,
     %                      ' EST OPPOSE au TETRAEDRE',NTEOP,
     %                      ' OPPOSE au TETRAEDRE',NTEOP2,
     %                      ' QUI DEVRAIT ETRE',NTE
                     PRINT*,'ntlefaco: NOTETR(',NTE,')=',
     %                      (NOTETR(K,NTE),K=1,8)
                     PRINT*,'ntlefaco: NOTETR(',NTEOP,')=',
     %                      (NOTETR(K,NTEOP),K=1,8)
                     PRINT*,'ntlefaco: NOTETR(',NTEOP2,')=',
     %                      (NOTETR(K,NTEOP2),K=1,8)

C                    LA FACE NFTE DE NTE APPARTIENT ELLE A 3 TETRAEDRES?
                     IF( NTEOP2.GT.0 .AND. NOTETR(1,NTEOP2).GT.0 ) THEN
                        CALL NO1F1T( NOSOTR, NOTETR(1,NTEOP2), NFOP2 )
                        IF( NFOP2 .GT. 0 ) THEN

C                          LA FACE NFTE DE NTE EST UNE FACE DE NTEOP2 ET
C                          ELLE APPARTIENT AUX 3 TETRAEDRES NTE NTEOP NTEOP2
                           PRINT*,'ntlefaco: PB: la FACE',NOSOTR,
     %                            ' EST FACE des 3 TETRAEDRES',
     %                             NTE, NTEOP, NTEOP2
                           PRINT*,'ntlefaco: NOTETR(',NTE,')=',
     %                            (NOTETR(K,NTE),K=1,8)
                           PRINT*,'ntlefaco: NOTETR(',NTEOP,')=',
     %                            (NOTETR(K,NTEOP),K=1,8)
                           PRINT*,'ntlefaco: NOTETR(',NTEOP2,')=',
     %                            (NOTETR(K,NTEOP2),K=1,8)
                           PRINT*

C                          le tetraedre nte a t il ete oublie
C                          d'etre detruit?
                           call retetrop( nte, mxtetr, notetr, nbteop )
                           if( nbteop .eq. 0 ) then
C                             le tetraedre nte est detruit
                              do k=1,8
                                 notetr(k,nte) = 0
                              enddo
ccc                              notetr(5,nte) = n1tevi
ccc                              n1tevi = nte
                              goto 10
                           endif

ccc                        ELSE

cccC                          LA FACE NFTE DE NTE N'EST PAS UNE FACE DE NTEOP2
cccC                          CORRECTION: LE TETRAEDRE OPPOSE DE NTEOP EST NTE
ccc                           NOTETR(4+NFOP,NTEOP) = NTE

                        ENDIF
                     ENDIF
                  ENDIF

               ENDIF

            ENDIF
         ENDDO


         IF( NOTETR(1,NTE) .GT. 0 ) THEN

C           LE TETRAEDRE NTE N'EST PAS VIDE
            NUDTETR = NTE

C           MISE A JOUR DU NO DE TETRAEDRE DE SES 4 FACES DANS LEFACO
            DO NFTE=1,4

C              NO D'UN TETRAEDRE DE SOMMET NS
               NS = NOTETR( NFTE, NTE )

C              LE TETRAEDRE NTE A POUR SOMMET NS
               N1TETS( NS ) = NTE

C              LE PLUS GRAND NO DE SOMMET D'UN TETRAEDRE EST IL NS?
               IF( NS .GT. NBSOMM ) THEN
                  NBSOMM = NS
               ENDIF

C              LA FACE NFTE DU TETRAEDRE NTE EST ELLE DANS LEFACO?
               CALL NULEFT( NFTE, NTE, NOTETR, MXFACO, LEFACO, NFLEFA )
               IF( NFLEFA .GT. 0 ) THEN

C                 LE TETRAEDRE NTE CONTIENT LA FACE NFLEFA DE LEFACO
                  LEFACO( 11, NFLEFA ) = NTE

C                 UN NUMERO DE FACE LEFACO AUX 3 SOMMETS
                  DO K=1,3
                     N1FASC( LEFACO(K,NFLEFA) ) = NFLEFA
                  ENDDO

               ENDIF

            ENDDO

         ELSE

            IF( NOTETR(1,NTE) .LT. 0 ) THEN
C              VALEUR INCORRECTE DU NO DU PREMIER SOMMET DE NTE
               PRINT*,'ntlefaco: Pb avec NOTETR(',NTE,')=',
     %                (NOTETR(k,NTE),k=1,8)
            ENDIF

         ENDIF

 10   ENDDO

C     CHAINAGE DES TETRAEDRES VIDES
      N1TEVI = 0
      DO NTE = MXTETR, 1, -1
         IF( NOTETR(1,NTE) .LE. 0 ) THEN
C           MARQUEURS DE NTE TETRAEDRE VIDE
            NOTETR(1,NTE) = 0
            NOTETR(5,NTE) = N1TEVI
            N1TEVI = NTE
         ENDIF
      ENDDO

C     BILAN
      NBFALE = 0
      NBFA1T = 0
      NBFAPE = 0

      DO NFLE = 1, MXFACO
         IF( LEFACO( 1, NFLE ) .NE. 0 ) THEN

C           FACE LEFACO UTILISEE
            NBFALE = NBFALE + 1

C           1 TETRAEDRE LA CONTIENT ELLE?
            NTEFA = LEFACO( 11, NFLE ) 
            IF( NTEFA .GT. 0 ) THEN
C              LA FACE NFLE APPARTIENT A UN TETRAEDRE
               NBFA1T = NBFA1T + 1
            ELSE
C              LA FACE NFLE APPARTIENT A AUCUN TETRAEDRE
               NBFAPE = NBFAPE + 1
               IF( NTEFA .LT. 0 ) THEN
C                 LA FACE NFLE DE LEFACO EST MARQUEE ABANDONNEE MOMENTANEMENT
                  NOFAPE( NBFAPE ) = -NFLE
               ELSE
                  NOFAPE( NBFAPE ) = NFLE
               ENDIF
ccc               print*,'ntlefaco: Face PERDUE LEFACO(',NFLE,'):',
ccc     %                (LEFACO(k,NFLE),k=1,11)
            ENDIF

         ENDIF
      ENDDO

cccC     TRI SELON LES SURFACES DECROISSANTES DES FACES PERDUES NOFAPE
ccc      CALL TRIFAP( NBFAPE, NOFAPE, LEFACO, PTXYZD )

C     AFFICHAGES des MODIFICATIONS EFFECTUEES
      PRINT*
      IF( LANGAG .EQ. 0 ) THEN

      PRINT*,'ntlefaco:',NBSOM0,' INITIAL PLUS GRAND No des SOMMETS des 
     %TETRAEDRES'
      PRINT*,'ntlefaco:',NBSOMM,' FINAL   PLUS GRAND No des SOMMETS des 
     %TETRAEDRES'
      PRINT*,'ntlefaco:',NUDTETR,' No du DERNIER TETRAEDRE OCCUPE dans N
     %OTETR PARMI les',MXTETR,' TETRAEDRES UTILISABLES'
      PRINT*,'ntlefaco:',N1TEVI,' No du PREMIER TETRAEDRE VIDE'
      PRINT*,'ntlefaco:',NBFALE,' FACES FRONTIERE UTILISEES sur',MXFACO,
     %                          ' UTILISABLES'
      PRINT*,'ntlefaco:',NBFA1T,' FACES FRONTIERE APPARTENANT a 1 TETRAE
     %DRE au moins donc RETROUVEES'
      PRINT*,'ntlefaco:',NBFAPE,' FACES FRONTIERE APPARTENANT a 0 TETRAE
     %DRE donc PERDUES'
      PRINT*,'ntlefaco:',NBPBOP,' FACES NON COMMUNES A 2 TETRAEDRES OPPO
     %SES'

      ELSE

      PRINT*,'ntlefaco:',NBSOM0,' INITIAL TETRAHEDRA VERTEX GREATEST NUM
     %BER'
      PRINT*,'ntlefaco:',NBSOMM,' FINAL   TETRAHEDRA VERTEX GREATEST NUM
     %BER'
      PRINT*,'ntlefaco:',NUDTETR,' OCCUPIED TETRAHEDRON LAST NUMBER amon
     %g',MXTETR,' USABLE TETRAHEDRA'
      PRINT*,'ntlefaco:',N1TEVI,' FIRST EMPTY TETRAHEDRON NUMBER'
      PRINT*,'ntlefaco:',NBFALE,' BOUNDARY FACES are USED over',
     %                   MXFACO,' USABLES'
      PRINT*,'ntlefaco:',NBFA1T,' BOUNDARY FACES into  1 TETRAHEDRON at
     %least are RECOVERED'
      PRINT*,'ntlefaco:',NBFAPE,' BOUNDARY FACES into NO TETRAHEDRON are
     % LOST'
      PRINT*,'ntlefaco:',NBPBOP,' NOT COMMON FACES BETWEEN TWO OPPOSED T
     %ETRAHEDRA'

      ENDIF
      PRINT*

C     TRACE EVENTUEL DES NBFAPE FACES PERDUES DANS LEFACO
ccc      tracte  = .true.
      CALL TRFAPE( NBFAPE, NOFAPE, MXFACO, LEFACO, PTXYZD )
      TRACTE = TRACTE0

      RETURN
      END
