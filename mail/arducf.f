      SUBROUTINE ARDUCF( PTXYZD, NBFAPE, NOFAPE, MXFACO, LEFACO, NO0FAR,
     %                   MXETOI, NASIDUCF, MXTRCF, NBTRCF, NOTRCF,
     %                   MXARCF, NBCF,     N1ARCF, NOARCF, NBARCF,
     %                   MXSTCF, NBSTCF, NOSTCF,
     %                   MXSTIS, NBSTIS, NOSTIS, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CONSTRUCTION DES SOMMETS DE LA FRONTIERE ET INTERNES AU CF
C -----    DES NBTRCF FACES PERDUES ADJACENTES (LEFACO ou NO0FAR) DU CF
C
C ENTREES:
C --------
C PTXYZD : TABLEAU DES COORDONNEES DES POINTS
C          PAR POINT : X  Y  DISTANCE_SOUHAITEE
C NBFAPE : NOMBRE DE  FACES PERDUES DE   LEFACO
C NOFAPE : NUMERO DES FACES PERDUES DANS LEFACO
C MXFACO : NOMBRE MAXIMAL DE TRIANGLES PERMIS POUR LA TRIANGULATION
C NO0FAR : NUMERO DES 3 SOMMETS DE LA FACE AJOUTEE AU CF
C MXTRCF : NOMBRE MAXIMAL DE TRIANGLES DANS NOTRCF
C MXSTCF : NOMBRE MAXIMAL DE SOMMETS DU CF DECLARABLES DANS NOSTCF
C MXSTIS : NOMBRE MAXIMAL DE SOMMETS ISOLES DU CF DECLARABLES DANS NOSTIS
C
C ENTREES ET SORTIES:
C -------------------
C LEFACO : ENSEMBLE DES FACES=TRIANGLES DE LA PEAU OU INTERFACES ENTRE VOLUMES
C          IL CONTIENT DANS CET ORDRE
C          123: NUMERO (DANS PTXYZD) DU SOMMET 1 < SOMMET 2 < SOMMET 3
C          45 : NUMERO (DANS NUVOPA 0 SINON) DU VOLUME1, VOLUME2 DE LA FACE
C          678: NUMERO (DANS LEFACO) DE LA FACE ADJACENTE PAR L'ARETE 1 2 3
C               ATTENTION: UNE ARETE PEUT APPARTENIR A PLUS DE 2 FACES
C                          => CHAINAGE CIRCULAIRE DE CES FACES DANS LEFACO
C          9  : LEFACO(9,*)  -> FACE SUIVANTE (*=0:VIDE, *<>0:NON VIDE)
C               HACHAGE AVEC LA SOMME DES 3 SOMMETS MODULO MXFACO
C               LF = MOD( NOSOFA(1)+NOSOFA(2)+NOSOFA(3) , MXFACO ) + 1
C               NF = LEFACO( 10, LF ) LE NUMERO DE LA 1-ERE FACE DANS LEFACO
C               SI LA FACE NE CONVIENT PAS. PASSAGE A LA SUIVANTE
C               NF  = LEFACO( 9, NF )  ...
C          10 : LEFACO(10,*) PREMIERE FACE DANS LE HACHAGE
C          11 : LEFACO(11,.) = NO NOTETR D'UN TETRAEDRE AYANT CETTE FACE, 0 SINON
cccC          12 : LEFACO(12,.) = NO FACEOC DE 1 A NBFACES D'OC

C MXETOI : NOMBRE MAXIMAL D'ARETES SIMPLES DES FACES PERDUES DU CF
C NASIDUCF:LISTE DES ARETES SIMPLES DES FACES PERDUES DU CF
C          POINTEE POUR LES ARETES OCCUPEES PAR N1ASDUCF
C          INUTILISABLE EN SORTIE CAR TABLEAU AUXILIAIRE ECRASE

C SORTIES :
C ---------
C NBCF   : NOMBRE DE LIGNES FERMEES PERIPHERIQUES DES FACES PERDUES
C NBTRCF : NOMBRE DE TRIANGLES DU TABLEAU NOTRCF SUR LEFACO
C          C-A-D DU POLYGONE ENCORE DIT ENSUITE ETOILE
C NOTRCF : TABLEAU DU NUMERO DANS LEFACO DES TRIANGLES DE L'ETOILE
C N1ARCF : NUMERO DU 1-ER SOMMET OU 1-ERE ARETE DU CONTOUR FERME
C NOARCF : 1:NUMERO DU SOMMET DE L'ARETE DE LA LIGNE DU CONTOUR FERME
C          2:NUMERO DANS NOARCF DE L'ARETE SUIVANTE DU CF
C          3:NUMERO DANS LEFACO DU TRIANGLE ADJACENT OPPOSE A L'ARETE
C            ou =0 SI ARETE D'UNE FACE NO0FAR
C NBARCF : NOMBRE D'ARETES   DU CF
C NBSTCF : NOMBRE DE SOMMETS DU CF
C NOSTCF : NUMERO PTXYZD DES NBSTCF SOMMETS PERIPHERIQUES DU CF
C NBSTIS : NOMBRE DE SOMMETS ISOLES DANS L'ETOILE
C NOSTIS : NUMERO DES SOMMETS ISOLES N'APPARTENANT PAS AU CONTOUR
C IERR   : 0 SI PAS D'ERREUR
C          1 ANOMALIE 2 TETRAEDRES NON REELLEMENT OPPOSES PAR UNE FACE
C          3 NBTRCF>MXTRCF => TROP DE TRIANGLES PERDUS DANS NOTRCF
C          4 NBSTIS>MXSTIS => TROP DE SOMMETS ISOLES
C          5 NBSTCF>MXSTCF => TROP DE SOMMETS DU CF
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET St PIERRE DU PERRAY & VEULETTES  Fevrier 2018
C2345X7..............................................................012
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)
      COMMON / TRTETR / STOPTE, TRACTE
      LOGICAL           STOPTE, TRACTE, TRACTE0
      DOUBLE PRECISION  PTXYZD(4,*)
      INTEGER           LEFACO(11,0:MXFACO),
     %                  NO0FAR(3,*),
     %                  NASIDUCF(4,MXETOI),
     %                  NOTRCF(MXTRCF),
     %                  NOSTCF(MXSTCF),
     %                  N1ARCF(0:MXETOI),
     %                  NOARCF(3,MXARCF),
     %                  NOSTIS(MXSTIS),
     %                  NOFAPE(NBFAPE),
     %                  NOSOTR(3)
      CHARACTER*128     KTITRE

      TRACTE0  = TRACTE
      NBTRCF00 = NBTRCF

      KTITRE='         FACES PERDUES du CF Recherche des POINTS ISOLES'
      WRITE(KTITRE(1:6),'(I6)') NBTRCF
      CALL TRDUCF( KTITRE, PTXYZD, NBTRCF, NOTRCF, 0, NOSTIS,
     %             LEFACO, NO0FAR )

C     PROTECTION DE NBFAPE INITIAL
 1    NBFAP0 = NBFAPE

C     NOMBRE DE CIRCUITS FERMES PERIPHERIQUES
      NBCF   = 0

C     NOMBRE D'ARETES DU CF
      NBARCF = 0

C     LE NOMBRE DE SOMMETS ISOLES DE L'ETOILE
      NBSTIS = 0

C     LE NOMBRE DES ARETES SIMPLES DES TRIANGLES DU CF
      NBASIDUCF = 0

C     CONSTRUCTION DES ARETES SIMPLES DES TRIANGLES DU CF
C     ---------------------------------------------------
      IF( NBTRCF .EQ. 1 ) THEN

C        UN SEUL TRIANGLE DANS NOTRCF
C        ============================
         NTRCF = NOTRCF( 1 )

C        CONSTRUCTION DIRECTE DES 3 ARETES SIMPLES
         DO N1 = 1, 3

C           FACE LEFACO
C           NO DES 2 SOMMETS DE L'ARETE N1 DE LA FACE NTRCF
            NS1 = LEFACO( N1, NTRCF )
            IF( N1 .EQ. 3 ) THEN
               N2 = 1
            ELSE
               N2 = N1+1
            ENDIF
            NS2 = LEFACO( N2, NTRCF )

C           LE TRIANGLE OPPOSE A NTRCF PAR SON ARETE N1
            NTROP = LEFACO( 5+N1, NTRCF )

C           LE SOMMET NBARCF D'UNE ARETE PERIPHERIQUE
            NBARCF = NBARCF + 1
            NOARCF( 1, NBARCF ) = NS1
C           L'ARETE SUIVANTE
            NOARCF( 2, NBARCF ) = NBARCF + 1
C           LE NUMERO DU TRIANGLE OPPOSE AU TRIANGLE DU CF AU DELA DE CETTE ARETE
            NOARCF( 3, NBARCF ) = NTROP

C           LE NUMERO PTXYZD DU SOMMET N1
            NOSTCF( N1 ) = NS1

         ENDDO
C        CHAINAGE EN BOUCLE
         NOARCF( 2, NBARCF ) = 1
         N1ARCF( 1 ) = 1
         NBCF   = 1
         NBSTCF = 3
         NBSTIS = 0
         GOTO 500

      ENDIF

C     PLUS D'UN TRIANGLE DANS NOTRCF
C     ==============================
      DO 20 K = 1, NBTRCF

         NTRCF = NOTRCF( K )
         NTRCN = -NTRCF

C        TRAITEMENT DES 3 ARETES DU TRIANGLE NTRCF DE LEFACO ou NO0FAR
         DO 10 N1 = 1, 3

            IF( NTRCF .GT. 0 ) THEN

C              FACE LEFACO
C              NO DES 2 SOMMETS DE L'ARETE N1 DE LA FACE NTRCF
               NS1 = LEFACO( N1, NTRCF )
               IF( N1 .EQ. 3 ) THEN
                  N2 = 1
               ELSE
                  N2 = N1+1
               ENDIF
               NS2 = LEFACO( N2, NTRCF )

C              LE TRIANGLE OPPOSE A NTRCF PAR SON ARETE N1
               NTROP = LEFACO( 5+N1, NTRCF )

            ELSE

C              FACE NO0FAR
C              NO DES 2 SOMMETS DE L'ARETE N1 DE LA FACE NTRCN
               NS1 = NO0FAR( N1, NTRCN )
               IF( N1 .EQ. 3 ) THEN
                  N2 = 1
               ELSE
                  N2 = N1+1
               ENDIF
               NS2 = NO0FAR( N2, NTRCN )

C              LE TRIANGLE OPPOSE A NTRCF PAR SON ARETE N1
               NTROP = 0

            ENDIF

C           CETTE ARETE NS1-NS2 EST ELLE DEJA UNE ARETE DES TRIANGLES 
C           QUI PRECEDENT DU CF?
            DO 5 KK = 1, NBTRCF

               IF( KK .EQ. K ) GOTO 5

               NTRCF2 = NOTRCF( KK )
               NTRCN2 = -NTRCF2

               DO NN1=1,3

                  IF( NTRCF2 .GT. 0 ) THEN

C                    FACE LEFACO
C                    NO DES 2 SOMMETS DE L'ARETE NN1 DE LA FACE NTRCF2
                     NSA1 = LEFACO( NN1, NTRCF2 )
                     IF( NN1 .EQ. 3 ) THEN
                        NN2 = 1
                     ELSE
                        NN2 = NN1+1
                     ENDIF
                     NSA2 = LEFACO( NN2, NTRCF2 )

                  ELSE

C                    FACE NO0FAR
                     NTRCN2 = -NTRCF2
C                    NO DES 2 SOMMETS DE L'ARETE NN1 DE LA FACE NTRCN2
                     NSA1 = NO0FAR( NN1, NTRCN2 )
                     IF( NN1 .EQ. 3 ) THEN
                        NN2 = 1
                     ELSE
                        NN2 = NN1+1
                     ENDIF
                     NSA2 = NO0FAR( NN2, NTRCN2 )

                  ENDIF

                  IF( NSA1 .EQ. NS2 .AND. NSA2 .EQ. NS1 .OR.
     %                NSA1 .EQ. NS1 .AND. NSA2 .EQ. NS2 ) THEN

C                    L'ARETE NS1-NS2 EST UNE ARETE DOUBLE DES TRIANGLES DU CF
C                    PASSAGE A L'ARETE SUIVANTE
                     GOTO 10


ccc  2020/01/02
ccc  ===========================================================================
ccc  MIS EN COMMENTAIRE CAR LEFACO(1,N) < LEFACO(2,N) < LEFACO(3,N) INCOMPATIBLE
ccc  AVEC LE VECTOR NORMAL DONNANT UNE INFORMATION INTERNE OU EXTERNE
ccc  DE PLUS UNE FACE D'UN INTERFACE LAQUELLE DES 2 NORMALES STOCKER?
ccc  ===========================================================================
ccc                  ELSE IF( NSA1 .EQ. NS1 .AND. NSA2 .EQ. NS2 ) THEN
cccC                    L'ARETE NS1-NS2 EST UNE ARETE DOUBLE DES TRIANGLES DU CF
cccC                    MAIS LE SENS EST A CHANGER POUR QUE LES NORMALES
cccC                    SOIENT DIRIGEES DANS LE MEME SENS
ccc                     IF( NTRCF2 .GT. 0 ) THEN
ccc                        M                = LEFACO(2,NTRCF2)
ccc                        LEFACO(2,NTRCF2) = LEFACO(3,NTRCF2) 
ccc                        LEFACO(3,NTRCF2) = M
ccc                        M                = LEFACO(6,NTRCF2)
ccc                        LEFACO(6,NTRCF2) = LEFACO(8,NTRCF2) 
ccc                        LEFACO(8,NTRCF2) = M
ccc                     ELSE
ccc                        M                = NO0FAR(2,NTRCN2)
ccc                        NO0FAR(2,NTRCN2) = NO0FAR(3,NTRCN2) 
ccc                        NO0FAR(3,NTRCN2) = M
ccc                     ENDIF
cccC                    PASSAGE A L'ARETE SUIVANTE
ccc                     GOTO 10


                  ENDIF

               ENDDO

 5          ENDDO

C           L'ARETE NS1-NS2 EST UNE ARETE SIMPLE DES TRIANGLES DU CF
C           STOCKAGE DANS NASIDUCF
            NBASIDUCF = NBASIDUCF + 1
C           LE TRIANGLE OPPOSE AU CF
            NASIDUCF( 1, NBASIDUCF ) = NTROP
C           LE SOMMET 1 DE L'ARETE
            NASIDUCF( 2, NBASIDUCF ) = NS1
C           LE SOMMET 2 DE L'ARETE
            NASIDUCF( 3, NBASIDUCF ) = NS2
C           CHAINAGES SUR LE SUIVANT
            NASIDUCF( 4, NBASIDUCF ) = NBASIDUCF+1

 10      ENDDO

 20   ENDDO

C     FIN DU CHAINAGE SUR LE SUIVANT
      NASIDUCF( 4, NBASIDUCF ) = 0

C     LA PREMIERE ARETE SIMPLE DU TABLEAU NASIDUCF
      N1ASDUCF = 1


C     LE CF FAIT IL 2 BOUCLES? 
C     I.E. UN SOMMET D'ARETE PERIPHERIQUE EST IL VU 4 FOIS?
C     =====================================================
      DO N = 1, NBASIDUCF
         DO K = 2, 3

C           LE NO DU SOMMET K DE L'ARETE N SIMPLE DU CF
            NS = NASIDUCF( K, N )
            NBVUNS = 0

            DO NN = 1, NBASIDUCF
               DO KK = 2, 3
                  IF( NASIDUCF( KK, NN ) .EQ. NS ) THEN
                     NBVUNS = NBVUNS + 1
                  ENDIF
               ENDDO
            ENDDO

            IF( NBVUNS .GE. 4 ) THEN

C              NS EST UN SOMMET COMMUN A AU MOINS 2 CIRCUITS
               print*,'arducf: le SOMMET',NS,
     %                ' APPARTIENT A 2 CIRCUITS FERMES'

               IF( K .EQ. 3 ) THEN
C                 NS DEVIENT LE PREMIER SOMMET DE L'ARETE N
                  NASIDUCF( 3, N ) = NASIDUCF( 2, N )
                  NASIDUCF( 2, N ) = NS
               ENDIF

C              NS DEVIENT LE PREMIER SOMMET DE LA PREMIERE ARETE
C              DES ARETES SIMPLES DU CF SAUF SI C'EST DEJA LE CAS
               IF( N .EQ. 1 .AND. K .EQ. 2 ) GOTO 30
C              LA PREMIERE ARETE SIMPLE DU TABLEAU NASIDUCF
               N1ASDUCF = N
C              CHAINAGE SUR LE SUIVANT DU DERNIER
               NASIDUCF( 4, NBASIDUCF ) = 1
C              L'ARETE PRECEDANT N DEVIENT LA DERNIERE DU CHAINAGE
               NASIDUCF( 4, N-1 ) = 0
               GOTO 30

            ENDIF

         ENDDO
      ENDDO

C     CHAINAGE DES SOMMETS DE ARCF VIDES
 30   N1ARCF(0) = 1
      DO I = 1, MXARCF
         NOARCF(2,I) = I+1
      ENDDO
      NOARCF(2,MXARCF) = 0

C     DECOUPAGE EN NBCF CONTOURS FERMES
      NBCF   = 0
      NBARCF = 0

C     LES ARETES SONT REORDONNEES CONSECUTIVEMENT POUR FORMER UNE LIGNE FERMEE
C     ========================================================================
C     UN CF PERIPHERIQUE DES TRIANGLES
 120  NBCF = NBCF + 1

C     LE NOMBRE DE SOMMETS DU CONTOUR FERME NBCF DE L'ETOILE
      NBARCF = NBARCF + 1

C     LE 1-ER SOMMET OU ARETE DU CONTOUR FERME NBCF
      N1ARCF( NBCF ) = NBARCF

C     LA PREMIERE ARETE A POUR PREMIER SOMMET NS0
      NA1 = N1ASDUCF
      NS0 = NASIDUCF(2,NA1)
      NS1 = NASIDUCF(3,NA1)

C     LE SOMMET NBARCF D'UNE ARETE PERIPHERIQUE
      NOARCF( 1, NBARCF ) = NS0
C     LE SOMMET SUIVANT
      NOARCF( 2, NBARCF ) = NBARCF + 1
C     LE NUMERO DU TRIANGLE OPPOSE AU TRIANGLE DU CF AU DELA DE CETTE ARETE
      NOARCF( 3, NBARCF ) = NASIDUCF(1,NA1)

C     DESTRUCTION DE L'ARETE NA1 DE NASIDUCF
      NA1 = NASIDUCF(4,NA1)
      N1ASDUCF = NA1


C     RECHERCHE DANS NASIDUCF DE L'ARETE SUIVANTE DE SOMMET DE DEPART NS1
 150  IF( NS1 .NE. NS0 ) THEN

C        RECHERCHE DANS NASIDUCF DE L'ARETE A PARTIR DE N1ASDUCF
         NA0 = N1ASDUCF
         NA1 = NA0
 160     IF( NA1 .GT. 0 ) THEN

C           LE NUMERO DU PREMIER SOMMET DE L'ARETE
            IF ( NS1 .NE. NASIDUCF(2,NA1) .AND.
     %           NS1 .NE. NASIDUCF(3,NA1) ) THEN
C              PASSAGE A L'ARETE SUIVANTE
               NA0 = NA1
               NA1 = NASIDUCF(4,NA1)
               GOTO 160
            ENDIF

C           ARETE PERIPHERIQUE RETROUVEE
            IF( NBARCF .GE. MXARCF ) THEN
               GOTO 9990
            ENDIF
            NBARCF = NBARCF + 1

C           NUMERO DU 1-ER SOMMET DE L'ARETE PERIPHERIQUE
            NOARCF( 1, NBARCF ) = NS1

C           LE NUMERO DU SECOND SOMMET DE L'ARETE
            IF( NS1 .EQ. NASIDUCF(2,NA1) ) THEN
               NS1 = NASIDUCF(3,NA1)
            ELSE
               NS1 = NASIDUCF(2,NA1)
            ENDIF

C           L'ARETE SUIVANTE
C           LE SOMMET SUIVANT DU DERNIER SOMMET EST LE PREMIER
C           CHAINAGE CIRCULAIRE
            NOARCF( 2, NBARCF ) = NBARCF+1

C           LE NUMERO DU TRIANGLE DE L'AUTRE COTE DE CELUI DU CF
            NOARCF( 3, NBARCF ) = NASIDUCF(1,NA1)
C
C           SUPPRESSION DE L'ARETE NASIDUCF DE SOMMET NS1
            IF( N1ASDUCF .EQ. NA1 ) THEN
                N1ASDUCF = NASIDUCF(4,NA1)
            ELSE
                NASIDUCF(4,NA0) = NASIDUCF(4,NA1)
            ENDIF

C           PASSAGE A L'ARETE SUIVANTE PERIPHERIQUE
            GOTO 150

         ENDIF

C        PROBLEME ARETE NON RETROUVEE : L'ETOILE NE SE REFERME PAS
         PRINT*,'arducf: ANOMALIE arducf A CORRIGER'
         PRINT*,'arducf: CF',NBCF,' NON REFERME Sommet',NS1,
     %                   ' NON RETROUVE'
         DO I=1,NBARCF
            PRINT*, 'NOARCF(',I,')=',(NOARCF(J,I),J=1,3)
         ENDDO
         PRINT*,NBTRCF,' TRIANGLES DE L''ETOILE'
         DO I=1,NBTRCF
            NTRCF = NOTRCF(I)
            IF( NTRCF .GT. 0 ) THEN
               PRINT*,'LEFACO(',NTRCF,')=',(LEFACO(J,NTRCF),J=1,11)
            ELSE
               PRINT*,'NO0FAR(',-NTRCF,')=',(NO0FAR(J,-NTRCF),J=1,3)
            ENDIF
         ENDDO
         IERR = 1
         GOTO 9999

      ENDIF

C     CHAINAGE EN BOUCLE
      NOARCF( 2, NBARCF ) = N1ARCF( NBCF )

C     CHAINAGE DES SOMMETS DE ARCF VIDES
      N1ARCF(0) = NBARCF+1

C     LE CF A T IL TOUTES SES ARETES PERIPHERIQUES ?
C     ==============================================
      IF( N1ASDUCF .NE. 0 .OR. NBASIDUCF .NE. NBARCF ) THEN
         GOTO 120

cccC        IL MANQUE DES ARETES => AJOUT DES TRIANGLES OPPOSES AUX FACES PERDUES
ccc 192     IF( N1ASDUCF .GT. 0 ) THEN
ccc            NTRCF = NASIDUCF( 1, N1ASDUCF )
cccC           TRACE DE LA FACE
ccc            CALL TRFATR( NCCYAN, NCROUG, LEFACO(1,NTRCF), PTXYZD )
ccc            DO I=1,NBFAPE
ccc               IF( NOFAPE(I) .EQ. NTRCF ) GOTO 198
ccc            ENDDO
cccC           FACE AJOUTEE COMME ETANT PERDUE
ccc            NBFAPE = NBFAPE + 1
ccc            NOFAPE( NBFAPE ) = NTRCF
cccC           SUPPRESSION DE L'ARETE SIMPLE
ccc 198        N1ASDUCF = NASIDUCF(4,N1ASDUCF )
ccc            GOTO 192
ccc         ENDIF
cccC        REDEPART
ccc         GOTO 1

      ENDIF

      IF( NBTRCF .GE. 57 ) THEN
C        AFFICHAGE DE NOARCF
         CALL AFNOARCF( NBCF, N1ARCF, NOARCF )
      ENDIF


      IF( NBTRCF .GT. 0 ) GOTO 300

CCCC    CE QUI EST DESSOUS jusqu'a 300 EST INUTILE ...  EN TEST...

C     LE CF A T IL DES TRIANGLES EXTERIEURS ENCOCHES?
C     ===============================================
C     RECHERCHE D'UN TRIANGLE AVEC 2 ARETES SUCCESSIVES DU CF
C     ET APPARTENANT A LEFACO SANS ETRE PERDUE
      NBTRCF0 = NBTRCF
      DO 220 NCF = 1, NBCF

C        LA PREMIERE ARETE DU CF
         NA0 = N1ARCF( NCF )
         NA1 = NA0

C        LES 3 ARETES CONSECUTIVES
 210     NA2 = NOARCF( 2, NA1 )
         NA3 = NOARCF( 2, NA2 )

C        LEURS 3 SOMMETS
         NOSOTR(1) = NOARCF(1,NA1)
         NOSOTR(2) = NOARCF(1,NA2)
         NOSOTR(3) = NOARCF(1,NA3)

C        LE TRIANGLE DE SOMMETS NOSOTR EST IL UN TRIANGLE DES
C        NBTRCF FACES NOTRCF DE LEFACO ou NO0FAR?
         CALL NTNOTRCF( NOSOTR, NBTRCF, NOTRCF, LEFACO, NO0FAR,
     %                  NFLEFA )

         IF( NFLEFA .GT. 0 ) THEN

C           OUI: LA FACE EST RETROUVEE DANS LEFACO
C           EST ELLE PERDUE C-A-D ELLE N'APPARTIENT PAS A UN TETRAEDRE

            IF( LEFACO(11,NFLEFA) .GT. 0 ) THEN
C              ELLE APPARTIENT A 1 TETRAEDRE. ELLE N'EST PAS PERDUE

C              ELLE RESSEMBLE A UNE ENCOCHE DANS LE CF

               DO K=1,NBTRCF
                  IF( ABS( NOTRCF(K) ) .EQ. NFLEFA ) GOTO 215
               ENDDO
C              AJOUT DE LA FACE NFLEFA A LA QUEUE DES TRIANGLES DU CF
               IF( NBTRCF .GE. MXTRCF ) GOTO 9990
               NBTRCF = NBTRCF + 1
               NOTRCF( NBTRCF ) = NFLEFA
               print*,'arducf: Ajout au CF de la Face LEFACO(',NFLEFA,
     %                '):',(LEFACO(kkk,NFLEFA),kkk=1,11)
               IF( NFLEFA.GT.0 .AND. LEFACO(11,NFLEFA).GT.0 ) THEN
C                 LA FACE NON PERDUE EST MARQUEE PERDUE
C                 POUR QUE NOTRCF SOIT LA LISTE DES FACES PERDUES DU CF
                  LEFACO(11,NFLEFA) = 0
               ENDIF

C              SUPPRESSION DE NA2 DES ARETES DU CF
C              NASIDUCF N'EST PLUS A JOUR MAIS N'EST PLUS UTILISE
 215           NA2 = NA3
C              CHAINAGE DE L'ARETE NA1 SUR NA3
               NOARCF(2,NA1) = NA3
            ENDIF

         ENDIF

C        PASSAGE A L'ARETE SUIVANTE
         IF( NA3 .NE. NA0 ) THEN
            NA1 = NA2
            GOTO 210
         ENDIF

 220  ENDDO

      IF( NBTRCF0 .LT. NBTRCF ) THEN
C        AU MOINS UNE FACE LEFACO A ETE AJOUTEE AU CF
         GOTO 1
      ENDIF


C     LE CF A T IL UN CIRCUIT INTERNE ?
C     =================================
C     RECHERCHE D'UN SOMMET 2 FOIS DANS UN CF
 300  DO 320 NCF = 1, NBCF

         NA0 = N1ARCF( NCF )
         NA1 = NA0
C
 301     NA2 = NOARCF( 2, NA1 )
C
 302     IF( NA2 .NE. NA0 ) THEN
            IF( NOARCF(1,NA1) .EQ. NOARCF(1,NA2) ) THEN

C              CROISEMENT DANS LES ARETES DU CF
               PRINT*,'arducf: NA1=',NA1,':',(NOARCF(kkk,NA1),kkk=1,3)
               PRINT*,'arducf: NA2=',NA2,':',(NOARCF(kkk,NA2),kkk=1,3)
               PRINT*,'arducf: CROISEMENT des ARETES du CF'

C              AFFICHAGE DE NOARCF
               CALL AFNOARCF( NBCF, N1ARCF, NOARCF )

C              TOUS LES TRIANGLES OPPOSES ENTRE NA1 ET LE PRECEDENT DE NA2
C              SONT CONSIDERES COMME PERDUS
C              ATTENTION: LEQUEL DES 2 PARCOURS POSSIBLES EST IL LE BON?
C              CHOIX POUVANT ETRE INCORRECT => LE PLUS COURT EN NOMBRE D'ARETES
               NBAR = 0
               NAZ  = NA1
 303           NAZ  = NOARCF(2,NAZ)
               NBAR = NBAR + 1
               IF( NAZ .NE. NA2 ) GOTO 303

               IF( NBAR .GT. NBARCF-NBAR ) THEN
C                 L'AUTRE PARCOURS DE NA2 A NA1 EST CHOISI
                  NAZ = NA2
                  NA2 = NA1
                  NA1 = NAZ
               ENDIF

C              LA FACE LEFACO ADJACENTE A L'ARETE NA1
 306           NTRCF = NOARCF( 3, NA1 )
               IF( NTRCF .LE. 0 ) THEN
                  print*,'arducf: TRIANGLE LEFACO ADJACENT=',NTRCF,' ??'
                  GOTO 315
               ENDIF

               DO I=1,NBFAPE
                  IF( NOFAPE(I) .EQ. NTRCF ) GOTO 310
               ENDDO

C              FACE AJOUTEE COMME ETANT PERDUE
               NBFAPE = NBFAPE + 1
               NOFAPE( NBFAPE ) = NTRCF

 310           DO K=1,NBTRCF
                  IF( ABS( NOTRCF(K) ) .EQ. NTRCF ) GOTO 315
               ENDDO

C              AJOUT DE LA FACE NTRCF AUX FACES PERDUES NOTRCF
               IF( NBTRCF .GE. MXTRCF ) GOTO 9990
               NBTRCF = NBTRCF + 1
               NOTRCF( NBTRCF ) = NTRCF
               print*,'arducf: AJOUT au CF de la Face LEFACO(',NTRCF,
     %                '):',(LEFACO(kkk,NTRCF),kkk=1,11)
               IF( NTRCF.GT.0 .AND. LEFACO(11,NTRCF).GT.0 ) THEN
C                 LA FACE NON PERDUE EST MARQUEE PERDUE
C                 POUR QUE NOTRCF SOIT LA LISTE DES FACES PERDUES DU CF
                  LEFACO(11,NTRCF) = 0
               ENDIF

 315           NA1 = NOARCF( 2, NA1 )
               IF( NA1 .NE. NA2 ) GOTO 306

C              REDEPART
               GOTO 1

            ELSE

C              PASSAGE A L'ARETE SUIVANTE
               NA2 = NOARCF( 2, NA2 )
               GOTO 302

            ENDIF
         ENDIF

C        PASSAGE A L'ARETE SUIVANTE
         NA1 = NOARCF( 2, NA1 )
         IF( NA1 .NE. NA0 ) GOTO 301

 320  ENDDO

C     ICI LE CF N'A PAS DE CROISEMENT et NBCF CONTOURS FERMES


C     CONSTRUCTION DU TABLEAU DES NBSTCF NUMEROS DES SOMMETS DU CF
C     ============================================================
      CALL CRSTCF( NBCF,   N1ARCF, NOARCF, NBARCF,
     %             MXSTCF, NBSTCF, NOSTCF, IERR )
      IF( IERR .NE. 0 ) GOTO 9999


C     CONSTRUCTION DU TABLEAU DES NBSTIS SOMMETS ISOLES ou INTERNES
C     DU CF C-A-D QU'ILS APPARTIENNENT AUX TRIANGLES DU CF MAIS 
C     N'APPARTIENNENT PAS AUX SOMMETS PERIPHERIQUES DU CONTOUR FERME
C     ==============================================================
      NBSTIS = 0
      DO J=1,NBTRCF
         NTRCF = NOTRCF(J)
         DO 400 I=1,3

C           NS SOMMET I DU TRIANGLE NTRCF
            IF( NTRCF .GT. 0 ) THEN
               NS = LEFACO( I, NTRCF )
            ELSE
               NS = NO0FAR( I, -NTRCF )
            ENDIF

            DO M=1,NBSTCF
               IF( NS .EQ. NOSTCF(M) ) GOTO 400
            ENDDO

C           LE SOMMET NS EST IL DEJA RETROUVE COMME SOMMET ISOLE?
            DO M=1,NBSTIS
               IF( NS .EQ. NOSTIS(M) ) GOTO 400
            ENDDO

C           NS EST DECOUVERT SOMMET ISOLE. IL EST AJOUTE A NOSTIS
            IF( NBSTIS .GE. MXSTIS ) THEN
               PRINT*,'arducf: TABLEAU NOSTIS SATURE. AUGMENTER MXSTIS='
     %               ,MXSTIS
               IERR = 4
               GOTO 9992
            ENDIF
ccc            TRACTE = .TRUE.
            NBSTIS = NBSTIS + 1
            NOSTIS( NBSTIS ) = NS
            PRINT*,'arducf:',NS,' EST UN SOMMET ISOLE'

 400     ENDDO
      ENDDO

cccC     AFFICHAGE A SUPPRIMER ENSUITE
ccc      print*,'fin arducf: NBTRCF=',NBTRCF
ccc      DO L=1,NBTRCF
ccc         print*,'fin arducf: NOTRCF(',L,')=',(LEFACO(K,NOTRCF(L)),K=1,3)
ccc      ENDDO
ccc      print*,'fin arducf: NOTRCF=',(NOTRCF(K),K=1,NBTRCF)
ccc      print*,'fin arducf: NBARCF=',NBARCF,' NBCF=',NBCF,
ccc     %       ' N1ARCF=',(N1ARCF(K),K=0,NBCF)
ccc      DO L=1,NBARCF
ccc         print*,'fin arducf: NOARCF(',L,')=',(NOARCF(K,L),K=1,3)
ccc      ENDDO

ccc 500  PRINT*,'arducf: NOMBRE de CONTOURS FERMES NBCF=',NBCF,' des',
ccc     %        NBTRCF,' FACES PERDUES avec',NBSTIS,' SOMMETS ISOLES'

ccc      if( nbcf .gt. 1 .OR. nbstis .gt. 0 ) then
ccc         tracte = .true.
ccc      endif

cccC     TRACE DES FACES TRIANGULAIRES DU CF, ARETES et SOMMETS
ccc      KTITRE='         FACES PERDUES du CF,        SOMMETS ISOLES,      
ccc     %  SOMMETS CF,        ARETES'
ccc      WRITE(KTITRE(1:6),  '(I6)') NBTRCF
ccc      WRITE(KTITRE(31:36),'(I6)') NBSTIS
ccc      WRITE(KTITRE(54:59),'(I6)') NBSTCF
ccc      WRITE(KTITRE(73:78),'(I6)') NBARCF
ccc      CALL  TRARTRCF( KTITRE, PTXYZD, NBTRCF, NOTRCF, 
ccc     %                NBCF,   N1ARCF, NOARCF, NBSTIS, NOSTIS,
ccc     %                NBSTCF, NOSTCF, LEFACO, NO0FAR )

 500  tracte = tracte0
      GOTO 9999


C     TABLEAU NOTRCF DE TAILLE INSUFFISANTE
 9990 PRINT*,'arducf: TROP DE FACES DANS LE CF. AUGMENTER MXTRCF=',
     %        MXTRCF
      IERR = 3
      GOTO 9999


C     TABLEAU NOSTIS DE TAILLE INSUFFISANTE
 9992 PRINT*,'arducf: TROP DE SOMMETS ISOLES DANS LE CF. AUGMENTER MXSTI
     %IS=',MXSTIS
      IERR = 4


 9999 IF( NBTRCF00 .LT. NBTRCF ) THEN
         PRINT*,'arducf: AUGMENTATION de',NBTRCF-NBTRCF00,
     %          ' FACES PERDUES NOTRCF du CF'
      ENDIF

      TRACTE = TRACTE0
      RETURN
      END
