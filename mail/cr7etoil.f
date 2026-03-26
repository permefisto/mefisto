      SUBROUTINE CR7ETOIL( QUAMINEX, KTITRE, PTXYZD, NFLPER,
     %                     MXFACO, LEFACO, NO0FAR, NBTRCF, NOTRCF,
     %                     NBSTIS, NOSTIS, MXSTCF, NBSTCF, NOSTCF, 
     %                     NBCF,   MXARCF, N1ARCF, NOARCF, 
     %                     N1TEVI, NOTETR, N1TETS, NUDTETR,
     %                     MXPTIN, NBPTIN, PTINTERS,
     %                     MXFETO, N1FEOC, N1FEVI, NFETOI,
     %                     MXTECF, NBTECF, NOTECF, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CONSTRUCTION DE L'ETOILE DES NBTRCF FACES PERDUES
C -----    AJOUT DES TETRAEDRES OPPOSES AUX TETRAEDRES DE MAUVAISE QUALITE
C          RETRAIT DES TETRAEDRES DE L'ETOILE AYANT
C          AU MOINS 3 FACES SIMPLES DE L'ETOILE ET AU PLUS 1 SOMMET DU CF
C          ou 2 FACES SIMPLES DE L'ETOILE ET 0 SOMMET DU CF

C ENTREES:
C --------
C QUAMINEX:QUALITE MINIMALE EXIGEE POUR CREER UN TETRAEDRE
C KTITRE : TITRE D'UN TRACE
C PTXYZD : TABLEAU DES COORDONNEES DES POINTS
C          PAR POINT : X  Y  Z DISTANCE_SOUHAITEE

C NFLPER : NUMERO LEFACO DE LA FACE PERDUE A TRAITER
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
cccC          12: = NO FACEOC DE 1 A NBFACES D'OC
C NO0FAR : NUMERO DES 3 SOMMETS DES FACES AJOUTEES AU CF

C NBTRCF : NOMBRE DE FACES DE NOTRCF
C NOTRCF : >0 NUMERO DANS LEFACO DES TRIANGLES PERDUS  DU CF
C          <0 NUMERO DANS NO0FAR DES TRIANGLES AJOUTES AU CF
C NBSTIS : NOMBRE DE SOMMETS ISOLES DU CF
C NOSTIS : NUMERO PTXYZD DES NBSTIS SOMMETS ISOLES
C NBSTCF : NOMBRE DE SOMMETS DES ARETES PERIPHERIQUES DU CF
C NOSTCF : NUMERO PTXYZD DES NBSTCF SOMMETS DES ARETES PERIPHERIQUES

C NBCF   : NOMBRE DE LIGNES FERMEES PERIPHERIQUES DES FACES PERDUES
C MXARCF : MAXIMUM D'ARETES DECLARABLES DANS N1ARCF et NOARCF
C N1ARCF : POINTE SUR LE DEBUT DES ARETES DE CHAQUE LIGNE FERMEE DU CF
C          0 POUR LES PLACES VIDES
C NOARCF : 1:NUMERO DU SOMMET DE L'ARETE DE LA LIGNE DU CONTOUR FERME
C          2:NUMERO DANS NOARCF DE L'ARETE SUIVANTE DU CF
C          3:NUMERO DANS LEFACO DU TRIANGLE ADJACENT OPPOSE A L'ARETE

C N1TEVI : NO DU PREMIER TETRAEDRE VIDE DANS NOTETR
C NOTETR : LISTE DES TETRAEDRES
C          SOMMET1,    SOMMET2,    SOMMET3,    SOMMET4,
C          TETRAEDRE1, TETRAEDRE2, TETRAEDRE3, TETRAEDRE4
C          DE L'AUTRE COTE DE LA FACE
C N1TETS : N1TETS(NS) NUMERO D'UN TETRAEDRE AYANT POUR SOMMET NS
C NUDTETR: NU NOTETR DU DERNIER TETRAEDRE ACTIF

C MXPTIN : NOMBRE MAXIMAL DE POINTS D'INTERSECTION ARETE TRIANGLE
C NBPTIN : NOMBRE DE POINTS D'INTERSECTION
C PTINTERS : XYZ DES NBPTIN POINTS D'INTERSECTION
C MXFETO : NOMBRE MAXIMAL DE FACES DECLARABLES DANS LE CHAINAGE NFETOI

C SORTIES:
C --------
C N1FEOC : NUMERO NFETOI DE LA PREMIERE FACE SIMPLE DES TETRAEDRES
C N1FEVI : NUMERO NFETOI DE LA PREMIERE FACE VIDE DE NFETOI
C NFETOI : AU DEBUT VERSION1  FACES DE L'ETOILE
C          1: NUMERO DU TETRAEDRE DANS NOTETR AYANT CETTE FACE
C          2: NUMERO LOCAL AU TETRAEDRE DE LA FACE DE L'ETOILE
C             UN SIGNE NEGATIF INDIQUE UN TRAITEMENT EFFECTUE
C          3: NON UTILISE ICI
C          4: NUMERO DE CETTE FACE DANS LEFACO, 0 SI PAS DANS LEFACO
C          5: CHAINAGE SUIVANT DES FACES OCCUPEES ET VIDES
C          ENSUITE VERSION2
C          1: NUMERO DU TETRAEDRE DANS NOTETR OPPOSE A CETTE FACE
C          2: NUMERO PTXYZD DU SOMMET 1 DE LA FACE
C          3: NUMERO PTXYZD DU SOMMET 2 DE LA FACE
C          4: NUMERO PTXYZD DU SOMMET 3 DE LA FACE
C             S1S2xS1S3 EST DIRIGE VERS L'INTERIEUR DE L'ETOILE
C          5: CHAINAGE SUIVANT DES FACES OCCUPEES ET VIDES
C NBTECF : NOMBRE DE TETRAEDRES DE L'ETOILE
C NOTECF : NUMERO NOTETR DES NBTECF TETRAEDRES DE L'ETOILE
C IERR   : = 0 SI PAS D'ERREUR DETECTEE
C          =-1 SI LA FACE NFLPER A ETE RETROUVEE DANS UN TETRAEDRE
C          = 5 SATURATION du TABLEAU NOTECF
C          > 0 ERREUR DETECTEE DANS L'UN DES APPELS de SOUS-PROGRAMMES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  St PIERRE du PERRAY             Janvier 2018
C2345X7..............................................................012
      PARAMETER        (MXSTSU=256)
      CHARACTER*(*)     KTITRE

      DOUBLE PRECISION  PTXYZD(4,*)
      DOUBLE PRECISION  PTINTERS(3,MXPTIN)
      DOUBLE PRECISION  ARMIN, ARMAX, SURFTR(4), V
      REAL              Q5TE(5)

      INTEGER           NOTETR(8,*), N1TETS(*),
     %                  LEFACO(11,0:MXFACO), NO0FAR(3,*),
     %                  NOTRCF(NBTRCF), NOSTCF(MXSTCF), NOSTIS(NBSTIS),
     %                  N1ARCF(0:*), NOARCF(1:3,1:*), NFETOI(5,MXFETO),
     %                  NOTECF(MXTECF),
     %                  NOSOTR(3),  NTEOLD(3), NTENEW(2), NVOLTE(1)

      INTEGER           NOSOARTE(2,6), NOSOFATE(3,4)
      DATA              NOSOARTE / 1,2,  2,3,  3,1,  1,4,  2,4,  3,4 /
      DATA              NOSOFATE / 1,3,2,  2,3,4,  3,1,4,  4,1,2 /
C     => LA NORMALE A LA FACE EST ORIENTEE VERS L'EXTERIEUR DU TETRAEDRE

      IERR   = 0
ccc      PRINT*,'cr7etoil: Recherche FACE PERDUE(',NFLPER,')=',NFLPER

      NBTECF0 = NBTECF

      IF( NBTECF .LE. 0 ) THEN

C     INACTIF...........................................................
C     ------------------------------------------------------------------
C     6') SI UNE ARETE DU CF N'A PAS DE TETRAEDRE ADJACENT DANS L'ETOILE
C        ALORS AJOUT DU TETRAEDRE D'ANGLE MINIMAL S'ENROULANT AUTOUR.
C        DE MEME POUR TOUTE ARETE JOIGNANT 2 SOMMETS DU CF
C     ------------------------------------------------------------------
      NBARPCF = 0

      IF( NBTRCF .EQ. 1 ) THEN

C        1 SEUL TRIANGLE LEFACO FORME LE CF: NBSTCF=3
C        --------------------------------------------
         DO N1 = 1, NBSTCF

            NBARPCF = NBARPCF + 1

C           LES 2 SOMMETS DE L'ARETE N DU TRIANGLE CF
            NS1 = NOSTCF( N1 )

            IF( N1 .LT. NBSTCF ) THEN
               N2 = N1 + 1
            ELSE
               N2 = 1
            ENDIF
            NS2 = NOSTCF( N2 )

C           LE 3-EME SOMMET DU TRIANGLE CF
            IF( N2 .LT. NBSTCF ) THEN
               N3 = N2 + 1
            ELSE
               N3 = 1
            ENDIF
            NS3 = NOSTCF( N3 )

C           AJOUT DES TETRAEDRES D'ARETE NS1 NS2 ET D'ANGLE DE COSINUS
C           MINIMUM AVEC LE TRIANGLE NOSTCF NS1 NS2 NS3
            NBT0 = NBTECF
            CALL TETR1ACMI( NS1,NS2, NS3,  N1TETS, NOTETR, PTXYZD,
     %                      NBSTCF, NOSTCF,
     %                      NBTECF, MXTECF, NOTECF, NTEMX, IERR )

            IF( NBT0 .LT. NBTECF ) THEN

C              UNE DES FACES DES TETRAEDRES AJOUTES EST ELLE LA FACE
C              LEFACO NFLPER PERDUE?
               NOSOTR(1) = LEFACO(1,NFLPER)
               NOSOTR(2) = LEFACO(2,NFLPER)
               NOSOTR(3) = LEFACO(3,NFLPER)
               DO NN = NBT0+1, NBTECF
                  NTE =  NOTECF( NN )
                  CALL NUFATRTE( NOSOTR, NOTETR(1,NTE), NF )
                  IF( NF .GT. 0 ) THEN

C                    LA FACE NF DU TETRAEDRE NTE EST LA FACE NFLPER DE LEFACO
                     LEFACO( 11, NFLPER ) = NTE
                     PRINT*,'cr7etoil: FACE RETROUVEE LEFACO(',NFLPER,
     %                      '): St',NOSOTR,' est LA FACE',NF,
     %                      ' du TETRAEDRE NOTETR(', NTE,')=',
     %                       (NOTETR(MM,NTE),MM=1,8)
                     IERR = -1
                     GOTO 9999

                  ENDIF
               ENDDO
            ENDIF

         ENDDO

      ELSE IF( NBTRCF .EQ. 2 ) THEN

C        2 TRIANGLES LEFACO FORMENT LE CF: NBSTCF=4
C        ------------------------------------------
         DO N = 1, NBTRCF

C           NUMERO LEFACO DE LA FACE PERDUE N DU CF
            NOTR = NOTRCF( N )

            DO N1 = 1,3

               NBARPCF = NBARPCF + 1

C              LES 2 SOMMETS DE L'ARETE N1 DU TRIANGLE CF
               NS1 = LEFACO( N1, NOTR )

               IF( N1 .LT. 3 ) THEN
                  N2 = N1 + 1
               ELSE
                  N2 = 1
               ENDIF
               NS2 = LEFACO( N2, NOTR )

C              LE 3-EME SOMMET DU TRIANGLE CF
               IF( N2 .LT. 3 ) THEN
                  N3 = N2 + 1
               ELSE
                  N3 = 1
               ENDIF
               NS3 = LEFACO( N3, NOTR )

C              AJOUT DES TETRAEDRES D'ARETE NS1 NS2 ET D'ANGLE DE
C              COSINUS MINIMUM AVEC LE TRIANGLE NOSTCF NS1 NS2 NS3
               NBT0 = NBTECF
               CALL TETR1ACMI( NS1, NS2, NS3,  N1TETS, NOTETR, PTXYZD,
     %                         NBSTCF, NOSTCF,
     %                         NBTECF, MXTECF, NOTECF, NTEMX, IERR )

               IF( NBT0 .LT. NBTECF ) THEN

C                 UNE DES FACES DES TETRAEDRES AJOUTES EST ELLE LA FACE
C                 LEFACO NFLPER PERDUE?
                  NOSOTR(1) = LEFACO(1,NOTR)
                  NOSOTR(2) = LEFACO(2,NOTR)
                  NOSOTR(3) = LEFACO(3,NOTR)
                  DO NN = NBT0+1, NBTECF
                     NTE =  NOTECF( NN )
                     CALL NUFATRTE( NOSOTR, NOTETR(1,NTE), NF )
                     IF( NF .GT. 0 ) THEN

C                       LA FACE NF DU TETRAEDRE NTE EST LA FACE NFLPER DE LEFACO
ccc                        PRINT*,'cr7etoil: NBTRCF=2 la FACE',NF,
ccc     %                         ' DU TETRAEDRE NOTETR(', NTE,')=',
ccc     %                         (NOTETR(MM,NTE),MM=1,8),
ccc     %                         ' EST la FACE',NFLPER,' RETROUVEE'

                        LEFACO( 11, NFLPER ) = NTE
C                       LA FACE PERDUE NFLPER DE LEFACO EST RETROUVEE
                        PRINT *,'cr7etoil: FACE RETROUVEE LEFACO',NFLPER
     %                         ,' NBTRCF=',NBTRCF
                        IERR = -1
                        GOTO 9999
  
                     ENDIF
                  ENDDO

               ENDIF

            ENDDO

         ENDDO

C        AJOUT DES EVENTUELS TETRAEDRES D'ARETES LA DIAGONALE S1-S3 DU CF
         NS1 = NOSTCF( 1 )
         NS2 = NOSTCF( 3 )
C        AJOUT DE TOUS LES TETRAEDRES D'ARETE S1 S3
         CALL TETR1A( NS1,   NS2,   N1TETS, NOTETR,
     %                NBTEA, MXTECF-NBTECF, NOTECF(NBTECF+1), IERR )
         NBTECF = NBTECF + NBTEA

C        UNICITE DES NO DE TETRAEDRES DU TABLEAU NOTECF
         CALL UNITABL( NOTECF, NBTECF )

C        AJOUT DES EVENTUELS TETRAEDRES D'ARETES LA DIAGONALE S2-S4 DU CF
         NS1 = NOSTCF( 2 )
         NS2 = NOSTCF( 4 )
C        AJOUT DE TOUS LES TETRAEDRES D'ARETE S2 S4 DU CF
         CALL TETR1A( NS1,   NS2,   N1TETS, NOTETR,
     %                NBTEA, MXTECF-NBTECF, NOTECF(NBTECF+1), IERR )
         NBTECF = NBTECF + NBTEA

C        UNICITE DES NO DE TETRAEDRES DU TABLEAU NOTECF
         CALL UNITABL( NOTECF, NBTECF )

      ELSE

C        PLUS DE 2 TRIANGLES LEFACO FORMENT LE CF
C        ----------------------------------------
         DO N1 = 1, NBSTCF

C           TOUTES LES ARETES POSSIBLES DU CF ISSUES DE NS1 SONT ENVISAGEES
            NS1 = NOSTCF( N1 )

            DO 8 N2 = N1+1, NBSTCF

               NBARPCF = NBARPCF + 1

C              LE SOMMET NS2 DE L'ARETE ENVISAGEE
               NS2 = NOSTCF( N2 )

               IF(N2 .EQ. N1+1 .OR. (N1.EQ.1 .AND. N2.EQ.NBSTCF))THEN

C                 NS1-NS2 EST UNE ARETE DU CF
C                 RECHERCHE DE LA FACE NOTRCF AYANT CETTE ARETE
                  DO N = 1, NBTRCF
                     NOTR = NOTRCF( N )
                     DO K1=1,3
                        IF( LEFACO(K1,NOTR) .EQ. NS1 ) THEN
                           DO K2 = 1, 3
                              IF( LEFACO(K2,NOTR).EQ. NS2 ) THEN
C                                NOTR CONTIENT L'ARETE NS1-NS2
C                                RECHERCHE DE SON 3-EME SOMMET
                                 DO K3 = 1, 3
                                    IF( K3.NE.K1 .AND. K3.NE.K2 )THEN
                                       NS3 = LEFACO(K3,NOTR)
                                       GOTO 3
                                    ENDIF
                                 ENDDO
                              ENDIF
                           ENDDO
                        ENDIF
                     ENDDO
                  ENDDO

C                 RECHERCHE DU TETRAEDRE D'ARETE NS1 NS2 D'ANGLE DE 
C                 COSINUS MINIMAL AVEC LA FACE PERDUE NOTR de
C                 SOMMETS NS1 NS2 NS3
 3                CALL TETR1ACMI( NS1,NS2,NS3, N1TETS,NOTETR, PTXYZD,
     %                            NBSTCF, NOSTCF,
     %                            NBTECF, MXTECF, NOTECF, NTEMX,IERR)

               ELSE

cccC        ESSAI 2:    L'ARETE NS1-NS2 EST ELLE EXTERIEURE AU CF?
cccC                    LE SOMMET CF PRECEDANT NS2
ccc                     IF( N2 .EQ. 1 ) THEN
ccc                        N3 = NBSTCF
ccc                     ELSE
ccc                        N3 = N2-1
ccc                     ENDIF
ccc                     NS3 = NOSTCF( N3 )
ccc
cccC                    RECHERCHE DU TETRAEDRE D'ARETE NS1 NS2 D'ANGLE DE
cccC                    COSINS MINIMAL AVEC LA FACE de SOMMETS NS2 NS1 NS3
ccc                     CALL TETR1ACMI( NS2,NS1,NS3, N1TETS,NOTETR, PTXYZD,
ccc     %                               NBSTCF, NOSTCF,
ccc     %                               NBTECF, MXTECF, NOTECF, NTEMX, IERR )

ccc         ESSAI 1:    N'EVITE PAS LES TETRAEDRES DES FACES ENCOCHES....
ccc                     NOSOTR(1) = NS1
ccc                     NOSOTR(2) = NS3
ccc                     NOSOTR(3) = NS2
cccC                    CETTE FACE EST ELLE DANS LEFACO?
ccc                     CALL NULETR( NOSOTR, MXFACO, LEFACO,  N )
ccc                     IF( N .GT. 0 ) THEN
cccC                       OUI: EST ELLE UNE FACE DU CF?
ccc                        DO NN = 1, NBTRCF
ccc                           IF( N .EQ. NOTRCF(NN) ) GOTO 6
ccc                        ENDDO
cccC                       NON: ELLE EST EVITEE POUR NE PAS REMPLIR UNE ENCOCHE
ccc                        GOTO 8
ccc                     ENDIF
cccC                    NON: AJOUT DE TOUS LES TETRAEDRES D'ARETE NS1 NS2
ccc 6                   CALL TETR1A( NS1,   NS2,   N1TETS, NOTETR,
ccc     %                            NBTEA, MXTECF-NBTECF, NOTECF(NBTECF+1), IERR )
ccc                     NBTECF = NBTECF + NBTEA

C                 AJOUT DE TOUS LES TETRAEDRES D'ARETE NS1 NS2 NON PERIPHERIQUE
C                 AVEC LA VERSION VDR1CF REMPLISSANT LES ENCOCHES DU CF
                  CALL TETR1A( NS1,   NS2,   N1TETS, NOTETR,
     %                         NBTEA, MXTECF-NBTECF, NOTECF(NBTECF+1),
     %                         IERR )
                  NBTECF = NBTECF + NBTEA
C                 NBTECF:EN ENTREE NOMBRE TETRAEDRES DEJA RANGES DANS NOTECF
C                        EN SORTIE AJOUT DES TETRAEDRES DE SOMMETS NS1 NS2
C                        = 0 SI SATURATION DU TABLEAU NOTECF
C                        =-1 SI UN SOMMET NS N'EST PAS UN SOMMET DE N1TETS(NS)

C                 UNICITE DES NO DE TETRAEDRES DU TABLEAU NOTECF
                  CALL UNITABL( NOTECF, NBTECF )

               ENDIF

 8          ENDDO

         ENDDO
      ENDIF

C     INACTIF...........................................................
      ENDIF



C     ---------------------------------------------------------------------
C     7) TENTATIVE DE SUPPRIMER LES TETRAEDRES DE QUALITE MEDIOCRE ou SINON
C        AJOUT DES TETRAEDRES OPPOSES AUX 4 FACES DES TETRAEDRES MEDIOCRES
C        POUR LES ENCAPSULER ET LES SUPPRIMER LORS DE LA TETRAEDRISATION
C     ---------------------------------------------------------------------
      NBTEC0 = NBTECF
      NB3T2T = 0
      DO 25 N=1,NBTEC0

         NTE = NOTECF( N )
         IF( NTE .LE. 0 ) GOTO 25
         IF( NOTETR(1,NTE) .LE. 0 ) GOTO 25

C        VOLUME ET QUALITE DE NTE
         CALL QUATETD( PTXYZD(1,NOTETR(1,NTE)),
     %                 PTXYZD(1,NOTETR(2,NTE)),
     %                 PTXYZD(1,NOTETR(3,NTE)),
     %                 PTXYZD(1,NOTETR(4,NTE)),
     %                 ARMIN, ARMAX, SURFTR, V, Q )

         IF( Q .LE. QUAMINEX ) THEN

C           TENTATIVE DE SUPPRIMER CE TETRAEDRE DE MAUVAISE QUALITE
C           A PARTIR D'UN COUPLE DE TETRAEDRES OPPOSES ET ADJACENTS
C           -------------------------------------------------------
ccc            PRINT*
ccc            PRINT*,'cr7etoil: 7) SUPPRIMER? MEDIOCRE TETRA(',
ccc     %              NTE,')=',(NOTETR(L,NTE),L=1,8),' V=',V,' Q=',Q

            DO 20 NAR = 1, 6

C              NUMERO PTXYZD DES 2 SOMMETS DE L'ARETE NAR
               NS1 = NOTETR( NOSOARTE( 1, NAR ), NTE )
               NS2 = NOTETR( NOSOARTE( 2, NAR ), NTE )

C              NS1-NS2 ARETE SIMPLE DES NBTRCF FACES DU CF?
               CALL ARARCF( NS1,   NS2,    NBCF, N1ARCF, NOARCF,
     %                      NLFCF, NAVANT, NA1 )
C              NA1: >0  NO NOARCF DE L'ARETE NS1->NS2 DANS NOARCF
C                   <0 -NO NOARCF DE L'ARETE NS2->NS1 DANS NOARCF
C                   =0  ARETE NS1-NS2 NON ARETE SIMPLE DES FACES DU CF
               IF( NA1 .NE. 0 ) GOTO 20

C              NON: TENTATIVE D'ECHANGER L'ARETE COMMUNE NAR
C                   DE 3 TETRAEDRES PAR UNE FACE A 2 TETRAEDRES
               CALL CH3T2T( 1,      MXFACO, LEFACO, 0, 0, NVOLTE,
     %                      NTE,    NAR,    PTXYZD,
     %                      N1TETS, NOTETR, N1TEVI, NUDTETR,
     %                      NTEOLD, NTENEW, Q5TE,   IERR )
C              IERR=0 SI PAS D'ERREUR ET 3T -> 2T REALISE

               IF( IERR .EQ. 0 ) THEN

                  PRINT*,'cr7etoil: 7) ECHANGE des TETRAEDRES',NTEOLD,
     %                   ' en', NTENEW,' POUR SUPPRIMER le TETRA(',
     %                   NTE,')  V=',V,' Q=',Q

                  NB3T2T = NB3T2T + 1

C                 SUPPRESSION DE NOTECF DES 3 TETRAEDRES NTEOLD
                  DO 15 M=1,3
                     NT = NTEOLD(M)
                     DO L=1,NBTECF
                        IF( NOTECF(L) .EQ. NT ) THEN
                           NOTECF(L) = -NT
                           GOTO 15
                        ENDIF
                     ENDDO
 15               ENDDO

C                 COMPRESSION DU TABLEAU NOTECF
                  M      = NBTECF
                  NBTECF = 0
                  DO L=1,M
                     NT = NOTECF(L)
                     IF( NT .GT. 0 .AND. NOTETR(1,NT) .GT. 0 ) THEN
                        NBTECF = NBTECF + 1
                        NOTECF(NBTECF) = NT
                     ENDIF
                  ENDDO

C                 AJOUT DANS NOTECF DES 2 TETRAEDRES NTENEW ET
C                 MISE A JOUR DE LEFACO POUR LES FACES DES NOUVEAUX TETRAEDRES
                  DO M=1,2

C                    LE TETRAEDRE M AJOUTE
                     NT = NTENEW(M)

C                    AJOUT DANS NOTECF DU TETRAEDRES NTENEW(M)
                     NBTECF = NBTECF + 1
                     NOTECF(NBTECF) = NT

                     DO NF=1,4
C                       LA FACE NF DU TETRAEDRE NT EST ELLE DANS LEFACO?
                        CALL NULEFT( NF, NT, NOTETR, MXFACO, LEFACO,
     %                               NOFACO )
                        IF( NOFACO .GT. 0 ) THEN

C                          LE TETRAEDRE NT CONTIENT LA FACE NOFACO DE LEFACO
                           LEFACO( 11, NOFACO ) = NT

                           IF( NOFACO .EQ. NFLPER ) THEN
C                             LA FACE PERDUE INITIALE EST RETROUVEE
                              PRINT*,'cr7etoil: FACE PERDUE LEFACO',
     %                               NFLPER,' RETROUVEE'
ccc                              GOTO 9999
                           ENDIF
                        ENDIF
                     ENDDO
                  ENDDO

                  GOTO 25

               ENDIF
 20         ENDDO

         ENDIF

 25   ENDDO

      NBTEC0 = NBTECF
      NBTECF = 0
      DO 30 N=1,NBTEC0
         NTE = NOTECF( N )
         IF( NTE .LE. 0 ) GOTO 30
         IF( NOTETR(1,NTE) .EQ. 0 ) GOTO 30
         NBTECF = NBTECF + 1
         NOTECF( NBTECF ) = NTE
 30   ENDDO
      NBTEC0 = NBTECF

      IF( NB3T2T .GT. 0 ) THEN
      KTITRE='cr7etoil: 7)       TETRAEDRES apres       ECHANGES 3T->2T'
         WRITE(KTITRE(14:18),'(I5)') NBTECF
         WRITE(KTITRE(37:41),'(I5)') NB3T2T
         CALL SANSDBL( KTITRE, L )
         PRINT*,KTITRE(1:L)
         CALL TRCFFAPE( KTITRE, PTXYZD, MXFACO, LEFACO, NO0FAR, NOTETR,
     %                  NBTECF, NOTECF, NBTRCF, NOTRCF,
     %                  MXARCF, NBCF,   N1ARCF, NOARCF,
     %                  NBSTIS, NOSTIS, NBPTIN, PTINTERS )
      ENDIF


C     ------------------------------------------------------------------
C     8) AJOUT DES TETRAEDRES OPPOSES AUX TETRAEDRES DE MAUVAISE QUALITE
C     ------------------------------------------------------------------
 80   NBTEAJ = 0
      DO 85 N=1,NBTECF

         NTE = NOTECF( N )
         IF( NTE .LE. 0 ) GOTO 85
         IF( NOTETR(1,NTE) .LE. 0 ) GOTO 85

C        VOLUME ET QUALITE DE NTE
         CALL QUATETD( PTXYZD(1,NOTETR(1,NTE)),
     %                 PTXYZD(1,NOTETR(2,NTE)),
     %                 PTXYZD(1,NOTETR(3,NTE)),
     %                 PTXYZD(1,NOTETR(4,NTE)),
     %                 ARMIN, ARMAX, SURFTR, V, Q )

CCC         IF( Q .LE. QUAMINEX ) THEN
         IF( Q .LE. 1E-3 ) THEN

C           NTE EST UN TETRAEDRE DE MAUVAISE QUALITE
C           AJOUT DES TETRAEDRES OPPOSES AUX FACES
            DO 82 NF=1,4
               NTEOP = NOTETR( 4+NF, NTE )
               IF( NTEOP .GT. 0 .AND. NOTETR(1,NTEOP) .GT. 0 ) THEN
C                 NTEOP EST IL DEJA DANS NOTECF?
                  DO K=1,NBTECF
                     IF( NOTECF(K) .EQ. NTEOP ) GOTO 82
                  ENDDO
C                 NON: IL EST AJOUTE
                  IF( NBTECF .GE. MXTECF ) THEN
C                    ABANDON
                     GOTO 9900
                  ENDIF
                  NBTEAJ = NBTEAJ + 1
                  NBTECF = NBTECF + 1
                  NOTECF( NBTECF ) = NTEOP
                  PRINT*,'cr7etoil: 8) AJOUt TETRA OPPOSE(',NTEOP,
     %                   ')=',(NOTETR(L,NTEOP),L=1,8),' V=',V,' Q=',Q
               ENDIF
 82         ENDDO

         ENDIF

 85   ENDDO

      IF( NBTEAJ .GT. 0 ) THEN
         KTITRE='cr7etoil: 8) AJOUT          TETRAEDRES OPPOSES aux TETR
     %AEDRES de QUALITE MEDIOCRE'
         WRITE(KTITRE(20:24),'(I5)') NBTEAJ
         CALL SANSDBL( KTITRE, L )
         PRINT*,KTITRE(1:L)
         CALL TRFETO4( KTITRE, PTXYZD, N1FEOC, NFETOI,
     %                 NBTRCF, NOTRCF, LEFACO, NO0FAR,
     %                 NBTECF, NOTECF, NOTETR )
         GOTO 80
      ENDIF


C     ------------------------------------------------------------------
C     9) RETRAIT DES TETRAEDRES DE L'ETOILE AYANT
C        AU MOINS 3 FACES SIMPLES DE L'ETOILE ET AU PLUS 1 SOMMETS DU CF
C        ou 2 FACES SIMPLES DE L'ETOILE ET 0 SOMMET DU CF
C     ------------------------------------------------------------------
C     CONSTRUCTION NFETOI VERSION 1 DES TRIANGLES FACES SIMPLES
C     DES TETRAEDRES. LES FACES VUES 2 FOIS SONT ELIMINEES
      NBTEC0 = NBTECF
 90   CALL CRFETOI1( NBTECF, NOTECF, NOTETR,
     %               MXFETO, N1FEOC, N1FEVI, NFETOI )

cccC     TRACE DES FACES SIMPLES DE L'ETOILE DU CF NFETOI VERSION 1
ccc      KTITRE='9)          TETRAEDRES avant RETRAIT TETRA 3F+2S ou 2F+0S'
ccc      CALL TRFETO4( KTITRE, PTXYZD, N1FEOC, NFETOI,
ccc     %              NBTRCF, NOTRCF, LEFACO, NO0FAR,
ccc     %              NBTECF, NOTECF, NOTETR )

      M = 0
      DO 98 N = 1, NBTECF

C        LE TETRAEDRE N DE L'ETOILE
         NTE = NOTECF( N )

C        PARCOURS DES FACES SIMPLES DE L'ETOILE EN VERSION 1
         NBFASI = 0
         NF1 = N1FEOC
 96      IF( NF1 .GT. 0 ) THEN
C           NUMERO NOTETR DU TETRAEDRE DE LA FACE SIMPLE
            IF( NTE .EQ. ABS( NFETOI(1,NF1) ) ) NBFASI = NBFASI + 1
C           LA FACE SUIVANTE DE L'ETOILE
            NF1 = NFETOI(5,NF1)
            GOTO 96
         ENDIF

C        NOMBRE DE SOMMETS DE NTE SUR LE CF OU ISOLES
         NBSTFR=0
         DO 97 K=1,4

            NS = NOTETR(K,NTE)
            DO L=1,NBSTCF
               IF( NS .EQ. NOSTCF(L) ) THEN
                  NBSTFR = NBSTFR + 1
                  GOTO 97
               ENDIF
            ENDDO

            DO L=1,NBSTIS
               IF( NS .EQ. NOSTIS(L) ) THEN
                  NBSTFR = NBSTFR + 1
                  GOTO 97
               ENDIF
            ENDDO

 97      ENDDO

ccc         IF( (NBFASI .GE. 3 .AND. NBSTFR .LE. 2)  .OR.

         IF( (NBFASI .GE. 3 .AND. NBSTFR .LE. 1)  .OR.
     %       (NBFASI .EQ. 2 .AND. NBSTFR .LE. 0) ) THEN
C           LE TETRAEDRE NTE EST SUPPRIME DE L'ETOILE
            PRINT*,'cr7etoil: FACE PERDUE LEFACO',NFLPER,
     %             ' SUPPRESSION TETRAEDRE',NTE,
     %             ' NBFASI=',NBFASI,' NBSTFR=',NBSTFR,
     %             ' St:', (NOTETR(L,NTE),L=1,8)
            GOTO 98
         ENDIF

C        TETRAEDRE CONSERVE
         M = M+1
         NOTECF(M) = NOTECF(N)

 98   ENDDO

      IF( M .LT. NBTECF ) THEN
         NBTECF = M
         GOTO 90
      ENDIF

      IERR = 0

      IF( NBTEC0 .NE. NBTECF ) THEN
      KTITRE='cr7etoil: 9)          TETRAEDRES apres RETRAIT TETRA 3F+1S
     % ou 2F+0S'
         WRITE(KTITRE(14:18),'(I5)') NBTECF
         CALL SANSDBL( KTITRE, L )
         PRINT*,KTITRE(1:L)
         CALL TRFETO4( KTITRE, PTXYZD, N1FEOC, NFETOI,
     %                 NBTRCF, NOTRCF, LEFACO, NO0FAR,
     %                 NBTECF, NOTECF, NOTETR )
      ENDIF


C     ----------------------------------------------------------------
C     10) SI UN SOMMET DU CF N'EST PAS UN SOMMET DES NBTECF TETRAEDRES
C         ALORS TOUS LES TETRAEDRES DE CE SOMMET SONT AJOUTES A NOTECF
C     ----------------------------------------------------------------
C     CONSTRUCTION DU TABLEAU DES NBSTCF NUMEROS DES SOMMETS DU CF
C     A PARTIR DE LA LISTE NOTRCF DES TRIANGLES DU CF
      CALL CRSTCF2( NBTRCF, NOTRCF, LEFACO, NO0FAR,
     %              MXSTCF, NBSTCF, NOSTCF, IERR )

      NBTEC0 = NBTECF
      DO 100 N = 1, NBSTCF

         NS = NOSTCF( N )

         DO K = 1, NBTECF
            NTE = NOTECF( K )
            IF( NTE .GT. 0 .AND. NOTETR(1,NTE) .GT. 0 ) THEN
               DO L = 1, 4
                  IF( NS .EQ. NOTETR(L,NTE) ) GOTO 100
               ENDDO
            ENDIF
         ENDDO

C        NS N'EST PAS UN SOMMET DES NBTECF TETRAEDRES
C        AJOUT DES TETRAEDRES DE CE SOMMET
         CALL TETR1S( NS,    N1TETS,        NOTETR,
     %                NBTNS, MXTECF-NBTECF, NOTECF(NBTECF+1), IERR )
         NBTECF = NBTECF + NBTNS
         IF( IERR .NE. 0 ) THEN
            PRINT*,'cr7etoil: SATURATION du TABLEAU NOTECF NBTECF=',
     %              NBTECF
         ENDIF

C        UNICITE DES NO DE TETRAEDRES DU TABLEAU NOTECF
         CALL UNITABL( NOTECF, NBTECF )

 100  ENDDO

      IF( NBTECF .NE. NBTEC0 ) THEN

C        CONSTRUCTION NFETOI VERSION 1 DES TRIANGLES FACES SIMPLES DES
C        TETRAEDRES ENVELOPPANTS. LES FACES VUES 2 FOIS SONT ELIMINEES
         CALL CRFETOI1( NBTECF, NOTECF, NOTETR,
     %                  MXFETO, N1FEOC, N1FEVI, NFETOI )

C        TRACE DES FACES SIMPLES DE L'ETOILE DU CF NFETOI VERSION 1
         KTITRE = 'cr7etoil: 10)          TETRAEDRES apres AJOUT des TET
     %RAEDRES des SOMMETS MANQUANTS du CF'
         WRITE(KTITRE(15:19),'(I5)') NBTECF
         CALL TRFETO4( KTITRE, PTXYZD, N1FEOC, NFETOI,
     %                 NBTRCF, NOTRCF, LEFACO, NO0FAR,
     %                 NBTECF, NOTECF, NOTETR )

      ENDIF

      GOTO 9999


 9900 PRINT*,'cr7etoil: SATURATION du TABLEAU NOTECF NBTECF=',
     %        NBTECF,' AUGMENTER MXTECF=',MXTECF
      IERR = 5

 9999 IF( NBTECF0 .NE. NBTECF ) THEN
         PRINT*,'cr7etoil: les',NBTECF0,' TETRAEDRES INITIAUX DONNENT',
     %           NBTECF,' TETRAEDRES de l''ETOILE'
      ENDIF

      RETURN
      END
