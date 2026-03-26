      SUBROUTINE EC2TMTNT( NPt,    PTXYZD, INFACO,  MXFACO, LEFACO,
     %                     IVOLTE, NVOLTE, MXTEET,  NBTEET, NOTEET,
     %                     NOTETR, N1TEVI, NUDTETR, N1TETS,
     %                     NBEC2T3T, NBECMTNT )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   AMELIORER LA QUALITE DES TETRAEDRES DE L'ETOILE DU POINT NPt
C -----   DEFINIE PAR LES NBTEET TETRAEDRES NOTEET DU TABLEAU NOTETR
C         PAR ECHANGE 2T->3T ou mT->2m-4T des TETRAEDRES

C ENTREES:
C --------
C NPt    : NUMERO DU POINT A TETRAEDRISER DE COORDONNEES DANS PTXYZD
C PTXYZD : X Y Z DISTANCE SOUHAITEE DES POINTS
C INFACO : = 0 PAS DE TABLEAU LEFACO NI DE CONSERVATION DES
C              FACES FRONTIERE CAR LE NO DE VOLUME CONNU PAR NVOLTE
C          = 1 EXISTENCE DU TABLEAU LEFACO ET CONSERVATION DES
C              FACES DE LA FRONTIERE
C MXFACO : NOMBRE MAXIMAL DECLARABLE DE FACES DU CONTOUR
C LEFACO : FACE DU CONTOUR OU INTERFACES ENTRE VOLUMES
C          IL CONTIENT DANS CET ORDRE
C          NUMERO (DANS PTXYZD) DU SOMMET 1, SOMMET 2, SOMMET 3
C          NUMERO (DANS NUVOPA 0 SINON) DU VOLUME1 , VOLUME2 DE LA FACE
C          NUMERO (DANS LEFACO) DE LA FACE ADJACENTE PAR L'ARETE 1 2 3
C          ATTENTION: UNE ARETE PEUT APPARTENIR A PLUS DE 2 FACES
C          => CHAINAGE CIRCULAIRE DE CES FACES DANS LEFACO
C          LEFACO(9,*)  -> FACE SUIVANTE (*=0:VIDE, *<>0:NON VIDE)
C          HACHAGE AVEC LA SOMME DES 3 SOMMETS MODULO MXFACO
C          LF = MOD( NOSOFA(1)+NOSOFA(2)+NOSOFA(3) , MXFACO ) + 1
C          NF = LEFACO( 10, LF ) LE NUMERO DE LA 1-ERE FACE DANS LEFACO
C          SI LA FACE NE CONVIENT PAS. PASSAGE A LA SUIVANTE
C          NF  = LEFACO( 9, NF )  ...
C IVOLTE : 0 PAS       DE TABLEAU NVOLTE A L'APPEL de chmtnt
C          1 EXISTENCE DU TABLEAU NVOLTE A L'APPEL de chmtnt
C NVOLTE : NUMERO DU VOLUME (1 A NBVOPA) DE CHAQUE TETRAEDRE
C          -1 SI VOLUME INCONNU
C MXTEET : NOMBRE MAXIMAL DE NUMEROS DANS NOTEET
C NBTEET : NOMBRE DE TETRAEDRES DE L'ETOILE DU POINT NPt EN ENTREE

C MODIFIES:
C ---------
C NOTEET : NUMERO DANS NOTETR DES NBTEET TETRAEDRES DE L'ETOILE DE NPt
C NOTETR : LISTE DES TETRAEDRES
C          SOMMET1,    SOMMET2,    SOMMET3,    SOMMET4,
C          TETRAEDRE1, TETRAEDRE2, TETRAEDRE3, TETRAEDRE4
C          DE L'AUTRE COTE DE LA FACE
C          1: 123      2: 234      3: 341      4: 412
C N1TEVI : NO DU PREMIER TETRAEDRE VIDE DANS NOTETR
C NUDTETR: PLUS GRAND NUMERO DANS NOTETR DE TETRAEDRE OCCUPE
C N1TETS : NUMERO DANS NOTETR D'UN TETRAEDRE DE CHAQUE SOMMET

C SORTIES:
C --------
C NBEC2T3T: NOMBRE D'ECHANGES DES TETRAEDRES 2T->3T
C NBECMTNT: NOMBRE D'ECHANGES DES TETRAEDRES 3T->2T
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR:ALAIN PERRONNET LJLL UPMC PARIS & St PIERRE du PERRAY  Mai 2017
C23456...............................................................012
      PARAMETER        (MXTOLD=4096, MXTNEW=4096)
      INTEGER           NOTETR(8,*), N1TETS(*),NOTEET(MXTEET),
     %                  NOTOLD(MXTOLD), NOTNEW(MXTNEW), NVOLTE(*)
      INTEGER           LEFACO(1:11,0:MXFACO)
      DOUBLE PRECISION  PTXYZD(4,*)
      COMMON / TRTETR / STOPTE, TRACTE
      LOGICAL           STOPTE, TRACTE, TRACTE0
      CHARACTER*128     KTITRE

C     NUMERO LOCAL DES 2 SOMMETS DES 6 ARETES DU TETRAEDRE
      INTEGER           NOSOARTE(2,6)
      DATA              NOSOARTE/ 1,2, 2,3, 3,1, 1,4, 2,4, 3,4 /

      TRACTE0 = TRACTE

C     ESSAI D'AMELIORER LES TETRAEDRES DE L'ETOILE PAR 2T->3T ou mT->2m-4T
C     ====================================================================
      NBEC2T3T = 0
      NBECMTNT = 0
      LHPILE = NBTEET

C     TRAITEMENT SOUS FORME DE PILE DU TABLEAU NOTEET
 1    IF( LHPILE .GT. 0 ) THEN

         NTE    = NOTEET( LHPILE )
         LHPILE = LHPILE - 1

         IF( NTE .GT. 0 .AND. NOTETR(1,NTE) .GT. 0 ) THEN

            DO NF=1,4

C              TENTATIVE DE CHANGER 2T->3T SI MIN QUALITE EST MAXIMISE
C              -------------------------------------------------------
               CALL CH2T3T( INFACO, MXFACO, LEFACO, 1, IVOLTE, NVOLTE,
     %                      NTE,    NF,     NTOP,
     %                      PTXYZD, N1TETS, NOTETR, N1TEVI, NUDTETR,
     %                      NOTNEW, IER )

               IF( IER .EQ. 0 ) THEN

C                 REUSSITE: 2T->3T
                  NBEC2T3T = NBEC2T3T + 1

ccc               PRINT*,'ec2tmtnt: NPt=',NPt,' 2T=>3T AMELIORE l''ETOILE',
ccc     %         NOTNEW,' Qualites=',QUATE(1),QUATE(2),' =>',
ccc     %         QUATE(3),QUATE(4),QUATE(5)
              
C                 SUPPRESSION DE NTOP DE NOTEET
                  DO K = 1, LHPILE
                     IF( NOTEET(K) .EQ. NTOP ) THEN
                        DO L = K+1, LHPILE
                           NOTEET( L-1 ) = NOTEET( L )
                        ENDDO
                        LHPILE = LHPILE - 1
                        GOTO 20
                     ENDIF
                  ENDDO

C                 LES 3 TETRAEDRES NOUVEAUX NOTNEW(1:3) SONT EMPILES
 20               IF( LHPILE+3 .GT. MXTEET ) THEN
                     GOTO 9900
                  ENDIF
                  DO K=3,1,-1
                     LHPILE = LHPILE + 1
                     NOTEET( LHPILE ) = NOTNEW( K )
                  ENDDO

C                 TRACE DES LHPILE TETRAEDRES DE LA PILE
                  IF( LHPILE .GT. 512 ) THEN
                  tracte = .true.
          KTITRE='ec2tmtnt: NPt=          PILE de            TETRAEDRES'
                  WRITE(KTITRE(15:22),'(I8)') NPt
                  WRITE(KTITRE(33:39),'(I7)') LHPILE
                  CALL SANSDBL( KTITRE, NBC )
                  CALL TRACEETOILE( KTITRE(1:NBC), PTXYZD, NOTETR, NPt,
     %                              LHPILE, NOTEET )
                  ENDIF


C                 RETOUR EN HAUT DE PILE
                  GOTO 1

               ENDIF

            ENDDO


C           TENTATIVE DE CHANGER mT->2m-4T SI MIN QUALITE EST MAXIMISE
C           ----------------------------------------------------------
C           PARCOURS DES ARETES DE NTE POUR TROUVER L'ARETE
C           AYANT LE MOINS DE TETRAEDRES S'ENROULANT AUTOUR
            NARMIN = 0
            NBTMIN = NUDTETR
            DO 30 NAR = 1, 6

C              NUMERO PTXYZD DES 2 SOMMETS DE L'ARETE NAR
               NS1 = NOTETR( NOSOARTE( 1, NAR ), NTE )
               NS2 = NOTETR( NOSOARTE( 2, NAR ), NTE )
               CALL TETR1A( NS1,    NS2,    N1TETS, NOTETR,
     %                      NBTOLD, MXTOLD, NOTOLD, IERR )
               IF( IERR .NE. 0 ) THEN
                  GOTO 30
               ENDIF
               IF( NBTOLD .LT. NBTMIN ) THEN
                  NBTMIN = NBTOLD
                  NARMIN = NAR
               ENDIF

 30         ENDDO

            IF( NARMIN .GT. 0 .AND. NBTMIN .LE. 10 ) THEN

               CALL CHMTNT( INFACO, MXFACO, LEFACO, IVOLTE, NVOLTE,
     %                      PTXYZD, NTE,    NARMIN,
     %                      N1TETS, NOTETR, N1TEVI, NUDTETR,
     %                      MXTOLD, NBTOLD, NOTOLD,
     %                      MXTNEW, NBTNEW, NOTNEW, IER )
               IF( IER .EQ. 0 ) THEN

C                 REUSSITE: mT->2m-4T
                  NBECMTNT = NBECMTNT + 1

C                 RETRAIT DES TETRAEDRES NOTOLD(1:NBTOLD) DE LA PILE
                  DO 40 M = 1, NBTOLD

C                    NO NOTETR DU TETRAEDRE A RETIRER DE LA PILE
                     NTEOL = NOTOLD( M )

                     DO K = 1, LHPILE
                        IF( NOTEET( K ) .EQ. NTEOL ) THEN
                           DO L = K+1, LHPILE
                              NOTEET( L-1 ) = NOTEET( L )
                           ENDDO
                           LHPILE = LHPILE - 1
                           GOTO 40
                        ENDIF
                     ENDDO

 40               ENDDO

C                 EMPILEMENT DES NOTNEW(1:NBTNEW)
                  DO 50 M = NBTNEW, 1, -1
                     IF( LHPILE .GE. MXTEET ) THEN
                        GOTO 9900
                     ENDIF
                     LHPILE = LHPILE + 1
                     NOTEET( LHPILE ) = NOTNEW( M )
 50               ENDDO

ccc                  PRINT*,'ec2tmtnt: NPt=',NPt,NBTOLD,'T->',
ccc     %            2*NBTOLD-4,'T AMELIORENT l''ETOILE'
ccc                  PRINT*,'de TETRAEDRES:',(NOTEET(K),K=1,LHPILE)
ccc                  PRINT*,'ec2tmtnt: Qualites Old=',
ccc     %                   (QUATE(K),K=1,NBTOLD),
ccc     %                   '  Qualites New=',
ccc     %                   (QUATE(K),K=NBTOLD+1,NBTOLD+NBTNEW)


C                 TRACE DES LHPILE TETRAEDRES DE LA PILE
                  IF( LHPILE .GT. 512 ) THEN
                  tracte = .true.
          KTITRE='EC2tmtnt: NPt=          PILE de            TETRAEDRES'
                  WRITE(KTITRE(15:22),'(I8)') NPt
                  WRITE(KTITRE(33:39),'(I7)') LHPILE
                  CALL SANSDBL( KTITRE, NBC )
                  CALL TRACEETOILE( KTITRE(1:NBC), PTXYZD, NOTETR, NPt,
     %                              LHPILE, NOTEET )
                  ENDIF


C                 RETOUR EN HAUT DE PILE
                  GOTO 1

               ENDIF

            ENDIF

         ENDIF

C        RETOUR EN HAUT DE PILE
         GOTO 1

      ENDIF

      TRACTE = TRACTE0
      RETURN

 9900 PRINT*,'ec2tmtnt: NPt=',NPt,
     %       ' PILE NOTEET SATUREE MXTEET=',MXTEET
      TRACTE = .TRUE.

      RETURN
      END
