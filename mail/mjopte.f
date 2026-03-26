      SUBROUTINE MJOPTE( NBTETRA, NOTETRA, N1TETS, NOTETR, NUDTETR,
     %                   N1TEVI,  PTXYZD,  NBFANR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  MISE A JOUR DES TETRAEDRES OPPOSES DANS UNE LISTE DE TETRAEDRES
C -----  ET SINON DANS TOUS LES TETRAEDRES ACTUELS de NOTETR

C ENTREES:
C --------
C NBTETRA: NOMBRE DE TETRAEDRES DE LA LISTE NOTETRA
C NOTETRA: NUMERO NOTETR DES TETRAEDRES DE LA LISTE
C N1TETS : N1TETS(I) NUMERO D'UN TETRAEDRE AYANT POUR SOMMET I
C NUDTETR: PLUS GRAND NUMERO DANS NOTETR DE TETRAEDRE OCCUPE
C PTXYZD : PAR POINT : X  Y  Z  DISTANCE_SOUHAITEE

C MODIFIES:
C ---------
C NOTETR : LISTE DES TETRAEDRES
C          SOMMET1,    SOMMET2,    SOMMET3,    SOMMET4,
C          TETRAEDRE1, TETRAEDRE2, TETRAEDRE3, TETRAEDRE4
C          DE L'AUTRE COTE DE LA FACE (0:FRONTIERE -1:INCONNU)
C          1: 123      2: 234      3: 341      4: 412
C N1TEVI : NO DU PREMIER TETRAEDRE VIDE DANS NOTETR

C SORTIE :
C --------
C NBFANR : NOMBRE DE FACES DE TETRAEDRE OPPOSE NON RETROUVE
C          LE NUMERO DE TETRAEDRE OPPOSE INCONNU EST MIS A -1
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY Decembre 2016
C2345X7..............................................................012
      PARAMETER        (MXFANR=2048, MXTAUX=2048)
      COMMON / TRTETR / STOPTE, TRACTE
      LOGICAL           STOPTE, TRACTE
      LOGICAL           TRACTE0
      CHARACTER*80      KTITRE
      DOUBLE PRECISION  PTXYZD(4,*)

      INTEGER  NOTETRA(NBTETRA), NOTETR(8,*), N1TETS(*)
      INTEGER  NOFANR(2,MXFANR), NOTAUX(MXTAUX), NOSOTR(3)
      INTEGER  NOSOFATE(3,4)
      DATA     NOSOFATE  / 1,3,2,  2,3,4,  3,1,4,  4,1,2 /
      
      TRACTE0 = TRACTE
 5    NBTAUX  = 0

C     NOMBRE DE FACES DE TETRAEDRE OPPOSE NON RETROUVE
      NBFANR = 0
      DO NT1 = 1, NBTETRA

C        NO NOTETR DU TETRAEDRE NT1
         NTE1 = NOTETRA( NT1 )
         IF( NTE1 .GT. 0 .AND. NOTETR(1,NTE1) .GT. 0 ) THEN

            DO L=1,4
C              MISE A JOUR PEUT ETRE INUTILE...
C              NUMERO D'UN TETRAEDRE DE SOMMET NOTETR(L,NTE1)
               N1TETS( NOTETR(L,NTE1) ) = NTE1
            ENDDO

            DO 20 NF1=1,4

C              NO NOTETR DU TETRAEDRE OPPOSE A LA FACE NF1 DE NTE1
               NTE2 = NOTETR( 4+NF1, NTE1 )

               IF( NTE2 .EQ. 0 ) THEN
C                 LA FACE NF1 DE NTE1 EST FRONTIERE
                  GOTO 20
               ENDIF

C              LES 3 SOMMETS DE LA FACE NF1 DE NTE1
               DO L=1,3
C                 NUMERO DU SOMMET L DE LA FACE NF1 DE NTE1
                  NOSOTR(L) = NOTETR( NOSOFATE(L,NF1), NTE1 )
               ENDDO

               IF( NTE2 .GT. 0 .AND. NOTETR(1,NTE2) .GT. 0 ) THEN

C                 VERIFICATION CONFIRMATION
C                 NTE2 EST IL BIEN OPPOSE A LA FACE NF1 DE NTE1? <=>
C                 LA FACE NF1 DE NTE1 NOSOTR EST ELLE UNE FACE NF2 DE NTE2?
C                 ---------------------------------------------------------
                  CALL NUFATRTE( NOSOTR, NOTETR(1,NTE2), NF2 )

                  IF( NF2 .GT. 0 ) THEN

C                    OUI: LA FACE NF1 DE NTE1 EST LA FACE NF2 DE NTE2
C                    ------------------------------------------------
                     NOTETR( 4+NF1, NTE1 ) = NTE2
                     NOTETR( 4+NF2, NTE2 ) = NTE1
                     GOTO 20

                  ELSE

C                    NON: NTE2 N'EST PAS OPPOSE A LA FACE NF1 DE NTE1
C                    ---- -------------------------------------------
ccc                     PRINT*
ccc                     PRINT*,'mjopte: NOSOFR=',NOSOTR,' NF1=',NF1
ccc                     PRINT*,'mjopte: NTE1=',NTE1,
ccc     %               ' NOTETR(',NTE1,')=',(NOTETR(K,NTE1),K=1,8)
ccc                     PRINT*,'mjopte: NTE2=',NTE2,
ccc     %               ' NOTETR(',NTE2,')=',(NOTETR(K,NTE2),K=1,8)

C                    LA FACE NF1 DE NTE1 EST MARQUEE DOUBLEMENT INCONNUE
                     NOTETR( 4+NF1, NTE1 ) = -2

                  ENDIF

               ENDIF


C              ICI NF2=0 : Le TETRAEDRE NTE2 OPPOSE A LA FACE NF1 est INCONNU
C              RECHERCHE DE LA FACE NF1 DE NTE1 PARMI TOUTES
C              LES FACES DES TETRAEDRES DE NOTETRA SAUF NTE1
C              --------------------------------------------------------------
               DO NT2 = 1, NBTETRA

C                 NO NOTETR DU TETRAEDRE NT2 SUPPOSE >0
                  NTE2 = NOTETRA( NT2 )
                  IF( NTE2 .NE. NTE1 ) THEN
                     IF( NTE2 .GT. 0 .AND. NOTETR(1,NTE2) .GT. 0 ) THEN

C                       LA  FACE NF1 DE NTE1 NOSOTR EST ELLE
C                       UNE FACE NF2 DE NTE2?
                        CALL NUFATRTE( NOSOTR, NOTETR(1,NTE2), NF2 )
                        IF( NF2 .GT. 0 ) THEN

C                          OUI: LA FACE NF1 DE NTE1 EST
C                               LA FACE NF2 DE NTE2
C                          ----------------------------
                           NOTETR( 4+NF1, NTE1 ) = NTE2
                           NOTETR( 4+NF2, NTE2 ) = NTE1
                           GOTO 20

                        ENDIF

                     ENDIF
                  ENDIF

               ENDDO


cccC              LA FACE NOSOTR EST ELLE FRONTIERE?
cccC              i.e. UN DE SES SOMMETS EST IL UN SOMMET (No de 1 A 6)
cccC              DE L'OCTAEDRE INITIAL?
cccC              -----------------------------------------------------
ccc               DO L=1,3
ccc                  IF( 1 .LE. NOSOTR(L) .AND. NOSOTR(L) .LE. 6 ) THEN
cccC                    LA FACE EST FRONTIERE 
ccc                     NOTETR( 4+NF1, NTE1 ) = 0
ccc                     GOTO 20
ccc                  ENDIF
ccc               ENDDO


C              FACE NON RETROUVEE DANS LES FACES DES NBTETRA TETRAEDRES
C              RECHERCHE DE LA FACE NF1 DE NTE1 PARMI TOUTES
C              LES FACES DES TETRAEDRES DE LA TETRAEDRISATION
C              A PARTIR DES TETRAEDRES DU PREMIER SOMMET DE NOSOTR
C              --------------------------------------------------------
               CALL TETR1F( NOSOTR, N1TETS, NOTETR,
     %                      NBTE1F, MXTAUX, NOTAUX, IERR )

               IF( IERR .NE. 0 ) GOTO 10

               IF( NBTE1F .EQ. 1 ) THEN

C                 LA FACE NF1 DE NTE1 EST UNIQUE => SUR LA FRONTIERE
C                 --------------------------------------------------
                  NOTETR( 4+NF1, NTE1 ) = 0
                  GOTO 20

               ENDIF


               IF( NBTE1F .LE. 0 .OR. NBTE1F .GT. 2 ) THEN
                  PRINT*,'mjopte: ANOMALIE LA FACE',NF1,
     %                  ' de NOTETR(',NTE1,')=',(NOTETR(KK,NTE1),KK=1,8)
                  PRINT*,'mjopte: FACE NOSOTR=',NOSOTR,' DE',
     %                    NBTE1F,' TETRAEDRES'

                  IF( NBTE1F .GT. 2 ) THEN
                     DO NT2 = 1, NBTE1F
C                       NO NOTETR DU TETRAEDRE NT2
                        NTE2 = NOTAUX( NT2 )
                        IF( NTE2.GT.0 .AND. NOTETR(1,NTE2).GT.0 ) THEN
                           PRINT*,' NOTETR(',NTE2,')=',
     %                             (NOTETR(KK,NTE2),KK=1,8)
                        ENDIF
                     ENDDO

C                    TRACE DES FACES NOFANR NON RETROUVEES
                     TRACTE = .TRUE.
                KTITRE='mjopte:       TETRAEDRES AYANT UNE FACE COMMUNE'
                    WRITE( KTITRE(8:12),'(I5)') NBTE1F
                    CALL TRACEETOILE( KTITRE, PTXYZD, NOTETR, NOSOTR(1),
     %                                NBTE1F, NOTAUX )
                  ENDIF
               ENDIF


C              LA FACE NF1 DE NTE1 APPARTIENT A NBTE1F>1 TETRAEDRES
C              ----------------------------------------------------
               DO NT2 = 1, NBTE1F

C                 NO NOTETR DU TETRAEDRE NT2
                  NTE2 = NOTAUX( NT2 )
                  IF( NTE2 .GT. 0 .AND. NOTETR(1,NTE2) .GT. 0 ) THEN
                     IF( NTE2 .NE. NTE1 ) THEN

C                       LA  FACE NF1 DE NTE1 NOSOTR EST ELLE
C                       UNE FACE NF2 DE NTE2?
                        CALL NUFATRTE( NOSOTR, NOTETR(1,NTE2), NF2 )
                        IF( NF2 .GT. 0 ) THEN
C                          OUI: LA FACE NF1 DE NTE1 EST
C                               LA FACE NF2 DE NTE2
C                          ----------------------------
                           NOTETR( 4+NF1, NTE1 ) = NTE2
                           NOTETR( 4+NF2, NTE2 ) = NTE1
                           GOTO 20
                        ENDIF

                     ENDIF
                  ENDIF

               ENDDO


C              RECHERCHE DE LA FACE NF1 DE NTE1 PARMI TOUTES
C              LES FACES DES NUDTETR TETRAEDRES ACTUELS DE NOTETR!
C              ---------------------------------------------------
 10            PRINT*
               PRINT*,'mjopte: BIZARRE TETRAEDRE(',NTE1,'):',
     %                (NOTETR(M,NTE1),M=1,8),' de FACE NOSOTR',NOSOTR
               PRINT*,'mjopte: FACE RECHERCHEE dans les',NUDTETR,
     %                ' TETRAEDRES ACTUELS'

               DO NTE2 = 1, NUDTETR

                  IF( NTE2 .NE. NTE1 .AND. NOTETR(1,NTE2) .GT. 0 ) THEN

C                    LA  FACE NF1 DE NTE1 NOSOTR EST ELLE
C                    UNE FACE NF2 DE NTE2?
                     CALL NUFATRTE( NOSOTR, NOTETR(1,NTE2), NF2 )
                     IF( NF2 .GT. 0 ) THEN

C                       OUI: LA FACE NF1 DE NTE1 EST
C                            LA FACE NF2 DE NTE2
C                       ----------------------------
                        NOTETR( 4+NF1, NTE1 ) = NTE2
                        NOTETR( 4+NF2, NTE2 ) = NTE1

                        PRINT*,'mjopte: RECHERCHEE PARMI les',NUDTETR,
     %                         ' TETRAEDRES la FACE NOSOTR',NOSOTR,
     %                         ' est RETROUVEE dans les 2 TETRAEDRES:'
                        PRINT*,'mjopte: TETRAEDRE(',NTE1,'):',
     %                         (NOTETR(M,NTE1),M=1,8)
                        PRINT*,'mjopte: TETRAEDRE(',NTE2,'):',
     %                         (NOTETR(M,NTE2),M=1,8)

                        GOTO 20
                     ENDIF

                  ENDIF

               ENDDO


C              LA FACE NF1 DE NTE1 N'A PAS DE TETRAEDRE OPPOSE
C              ===============================================
               NTE2 = NOTETR( 4+NF1, NTE1 )

               IF( NTE2 .LT. 0 ) THEN

                  PRINT*
                 PRINT*,'mjopte:',NTE2,' TETRAEDRE OPPOSE A LA FACE',NF1
     %              ,' de NTE1=',NTE1,' PARMI LES',NUDTETR,' TETRAEDRES'
                  PRINT*,'mjopte: AVANT NOTETR(',NTE1,')=',
     %                   (NOTETR(L,NTE1),L=1,8)
                  DO K=1,3
                     NS = NOSOTR( K )
                 PRINT*,'mjopte: PTXYZD(',NS,')=',(PTXYZD(KK,NS),KK=1,4)
                  ENDDO

                  IF( NBFANR .LT. MXFANR ) THEN
                     NBFANR = NBFANR + 1
C                    NUMERO NOTETR DU TETRAEDRE DE FACE SANS OPPOSE
                     NOFANR( 1, NBFANR ) = NTE1
C                    NUMERO DE LA FACE DANS LE TETRAEDRE
                     NOFANR( 2, NBFANR ) = NF1
                  ENDIF

C                 LE TETRAEDRE OPPOSE EST DONC TRIPLEMENT INCONNU
                  NOTETR( 4+NF1, NTE1 ) = -3

cccC                 LE TETRAEDRE OPPOSE EST DONC FRONTIERE
ccc                  NOTETR( 4+NF1, NTE1 ) = 0

                  PRINT*,'mjopte: FIN  La FACE',NF1,' de NOTETR(',NTE1,
     %                   ')=',(NOTETR(L,NTE1),L=1,8),' est INTROUVABLE'

               ENDIF

 20         ENDDO

         ENDIF

      ENDDO

C     EXISTE T IL DES TETRAEDRES AYANT 2 TETRAEDRES OPPOSES DOUBLES?
C     SI OUI: TRAITEMENT PAR SUPPRESSION DES 2 TETRAEDRES IDENTIQUES ET
C     RESOLUTION DES TETRAEDRES OPPOSES RESTANTS
C     -----------------------------------------------------------------
      DO 50 NT1 = 1, NBTETRA

C        NO NOTETR DU TETRAEDRE NT1
         NTE1 = NOTETRA( NT1 )

         IF( NTE1 .GT. 0 .AND. NOTETR(1,NTE1) .GT. 0 ) THEN

C           UN TETRAEDRE OPPOSE AU TETRAEDRE NTE1 EST IL DOUBLE OPPOSE?
            CALL TETOPDOU( NTE1, NOTETR, NF1, NF2 )
            IF( NF1 .GT. 0 ) THEN
C              OUI:
               NTE2 = NOTETR( 4+NF1, NTE1 )

C              NTE1 est il IDENTIQUE a NTE2?
               DO 30 K=1,4
                  NS1 = NOTETR( K, NTE1 )
                  DO L=1,4
                     NS2 = NOTETR( L, NTE2 )
                     IF( NS1 .EQ. NS2 ) THEN
                        GOTO 30
                     ENDIF
                  ENDDO
C                 NON: NTE1 EST DIFFERENT DE NTE2
                  PRINT*
             PRINT*,'ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ'
                  PRINT*,'mjopte: Probleme: le TETRAEDRE',NTE1,
     %                   ' est DIFFERENT du TETRAEDRE',NTE2
                  PRINT*,'mjopte: NOTETR(',NTE1,')=',
     %                   (NOTETR(kk,NTE1),kk=1,8)
                  PRINT*,'mjopte: NOTETR(',NTE2,')=',
     %                   (NOTETR(kk,NTE2),kk=1,8)
             PRINT*,'ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ'
                  PRINT*
                  GOTO 50
 30            ENDDO

C              OUI: NTE1 = NTE2. SUPPRESSION DES 2 TETRAEDRES
C              S'ILS SONT TETRAEDRES OPPOSES ILS DEVIENNENT INCONNUS
               DO K=1,4
                  NTOP = NOTETR( 4+K, NTE1 )
                  IF( NTOP .NE. NTE2 ) THEN
                     DO L=1,4
                        NT = NOTETR( 4+L, NTOP )
                        IF( NT .EQ. NTE1 .OR. NT .EQ. NTE2 ) THEN
                           NOTETR( 4+L, NTOP ) = -4
                        ENDIF
                     ENDDO
                  ENDIF
               ENDDO
               DO K=1,4
                  NTOP = NOTETR( 4+K, NTE2 )
                  IF( NTOP .NE. NTE1 ) THEN
                     DO L=1,4
                        NT = NOTETR( 4+L, NTOP )
                        IF( NT .EQ. NTE1 .OR. NT .EQ. NTE2 ) THEN
                           NOTETR( 4+L, NTOP ) = -4
                        ENDIF
                     ENDDO
                  ENDIF
               ENDDO

             PRINT*
             PRINT*,'ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ'
             PRINT*,'mjopte: RETRAIT des TETRAEDRES DOUBLES',NTE1, NTE2
             PRINT*,'mjopte: NOTETR(',NTE1,')=',(NOTETR(kk,NTE1),kk=1,8)
             PRINT*,'mjopte: NOTETR(',NTE2,')=',(NOTETR(kk,NTE2),kk=1,8)
             PRINT*,'ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ'
             PRINT*

C              RETRAIT DE NTE1 et NTE2 DU TABLEAU NOTETRA
               NOTETRA( NT1 ) = 0
               DO K=1,NBTETRA
                  IF( NOTETRA(K) .EQ. NTE2 ) THEN
                     NOTETRA(K) = 0
                     GOTO 33
                  ENDIF
               ENDDO

C              COMPRESSION DU TABLEAU NOTETRA
 33            L = 0
               DO K=1,NBTETRA
                  NT = NOTETRA( K )
                  IF( NT .NE. 0 ) THEN
                     L = L + 1
                     NOTETRA( L ) = NT
                  ENDIF
               ENDDO
               NBTETRA = L

C              DESTRUCTION DU TETRAEDRE NTE1 DANS NOTETR
               DO K=1,8
                  NOTETR( K, NTE1 ) = 0
               ENDDO
C              MISE A JOUR DU 1-ER TETRAEDRE VIDE
               NOTETR( 5, NTE1 ) = N1TEVI
               N1TEVI = NTE1

C              DESTRUCTION DU TETRAEDRE NTE2 DANS NOTETR
               DO K=1,8
                  NOTETR( K, NTE2 ) = 0
               ENDDO
C              MISE A JOUR DU 1-ER TETRAEDRE VIDE
               NOTETR( 5, NTE2 ) = N1TEVI
               N1TEVI = NTE2

               GOTO 5

            ENDIF
         ENDIF

 50   ENDDO


C     AFFICHAGE DU NOMBRE DES FACES DE TETRAEDRES OPPOSES NON RETROUVES
C     -----------------------------------------------------------------
      IF( NBFANR .GT. 0 ) THEN
         PRINT*,'mjopte:',NBFANR,' FACES de TETRAEDRE OPPOSE INCONNU et 
     %NON RETROUVEES dans les TETRAEDRES ACTUELS'

C        TRACE DES FACES NOFANR NON RETROUVEES
         TRACTE = .TRUE.
         KTITRE='mjopte:       FACES de TETRAEDRE OPPOSE NON RETROUVEES'
         WRITE( KTITRE(8:12),'(I5)') NBFANR
         CALL TRFETO12( KTITRE,  PTXYZD,
     %                  NBTETRA, NOTETRA, NOTETR,
     %                  NBFANR,  NOFANR )
      ENDIF


C     VERIFIER L'ABSENCE DE TETRAEDRE OPPOSE NEGATIF (INCONNU)
C     VERIFIER L'ABSENCE DE TETRAEDRES OPPOSES DOUBLES
C     VERIFIER L'ABSENCE DE FACE COMMUNE A AU MOINS 3 TETRAEDRES
C     VERIFIER LA TOTALE OPPOSITION DES NBTETRA TETRAEDRES 
C     ----------------------------------------------------------
      CALL VEOPTE( NBTETRA, NOTETRA, NOTETR, PTXYZD, NBFANR )

      IF( NBFANR .GT. 0 ) THEN
C        AFFICHAGE DU NOMBRE DES FACES NON RETROUVEES
C        NOFANR( 1, 1:NBFANR ) = NO NOTETR DU TETRAEDRE DE FACE SANS OPPOSE
C        NOFANR( 2, 1:NBFANR ) = NO DE LA FACE DANS LE TETRAEDRE
         PRINT*,'mjopte:',NBFANR,' FACES FRONTIERE ou NON OPPOSEE NON RE
     %TROUVEES dans les',NBTETRA,' TETRAEDRES'
      ENDIF

      TRACTE = TRACTE0
      RETURN
      END
