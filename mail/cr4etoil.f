      SUBROUTINE CR4ETOIL( QTEAME, KTITRE, PTXYZD, NOTETR,
     %                     MXFACO, LEFACO, NO0FAR, NBTRCF, NOTRCF,
     %                     NBCF,   MXARCF, N1ARCF, NOARCF, 
     %                     MXFETO, N1FEOC, N1FEVI, NFETOI,
     %                     MXTECF, NBTECF, NOTECF )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CONSTRUCTION 4 DE L'ETOILE DES NBTRCF FACES PERDUES
C -----    DE LA FRONTIERE PAR AJOUT DES TETRAEDRES OPPOSES COMMUNS
C          A 2 FACES SIMPLES DE L'ETOILE SI L'ARETE COMMUNE N'APPARTIENT
C          PAS AU CF ET SI L'ANGLE DES 2 FACES EST SUPERIEUR A 200 DEGRES

C ENTREES:
C --------
C QTEAME : QUALITE MINIMALE EXIGEE POUR CREER UN TETRAEDRE
C KTITRE : TITRE D'UN TRACE
C PTXYZD : TABLEAU DES COORDONNEES DES POINTS
C          PAR POINT : X  Y  Z DISTANCE_SOUHAITEE
C NOTETR : LISTE DES TETRAEDRES
C          SOMMET1,    SOMMET2,    SOMMET3,    SOMMET4,
C          TETRAEDRE1, TETRAEDRE2, TETRAEDRE3, TETRAEDRE4
C          DE L'AUTRE COTE DE LA FACE

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

C NBCF   : NOMBRE DE LIGNES FERMEES PERIPHERIQUES DES FACES PERDUES
C MXARCF : MAXIMUM D'ARETES DECLARABLES DANS N1ARCF et NOARCF
C N1ARCF : POINTE SUR LE DEBUT DES ARETES DE CHAQUE LIGNE FERMEE DU CF
C          0 POUR LES PLACES VIDES
C NOARCF : 1:NUMERO DU SOMMET DE L'ARETE DE LA LIGNE DU CONTOUR FERME
C          2:NUMERO DANS NOARCF DE L'ARETE SUIVANTE DU CF
C          3:NUMERO DANS LEFACO DU TRIANGLE ADJACENT OPPOSE A L'ARETE

C MXFETO : NOMBRE MAXIMAL DE FACES DECLARABLES DANS LE CHAINAGE NFETOI

C SORTIES:
C --------
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
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  Saint PIERRE du PERRAY             Juin 2018
C2345X7..............................................................012
      CHARACTER*(*)     KTITRE
      DOUBLE PRECISION  PTXYZD(4,*)
      DOUBLE PRECISION  ARMIN, ARMAX, SURFTR(4), VOLUTE,
     %                  ANGLE, PIS4, DATAN

      INTEGER           NOTETR(8,*), LEFACO(11,0:MXFACO), NO0FAR(3,*),
     %                  NOTRCF(NBTRCF),
     %                  N1ARCF(0:MXARCF), NOARCF(3,MXARCF),
     %                  NFETOI(5,MXFETO), NOTECF(MXTECF)

      INTEGER           NOSOFATE(3,4)
      DATA              NOSOFATE / 1,2,3,  2,4,3,  3,4,1,  4,2,1 /
C     => LA NORMALE A LA FACE EST ORIENTEE VERS L'INTERIEUR DU TETRAEDRE

      IF( NBTECF .LE. 0 ) RETURN

C     CONSTANTE de CONVERSION Pi/4 RADIANS <-> 45 DEGRES
      PIS4 = DATAN( 1D0 )

C     ---------------------------------------------------------
C     CONSTRUCTION DE L'ETOILE DES FACES SIMPLES DES TETRAEDRES
C     DANS NFETOI EN VERSION 1: TETRAEDRE OPPOSE ET NO FACE
C     AJOUT DES TETRAEDRES OPPOSES COMMUNS A 2 FACES SIMPLES DE
C     L'ETOILE SI AUCUNE DE SES ARETES EST ARETE SIMPLE DU CF
C     ---------------------------------------------------------
C     SUPPRESSION DES TETRAEDRES DESACTIVES DANS NOTECF
      NBTECF0 = NBTECF
      NBTECF  = 0
      DO 30 K = 1, NBTECF0
         NTE = NOTECF( K )
         IF( NTE .GT. 0 .AND. NOTETR(1,NTE) .EQ. 0 ) GOTO 30
         NBTECF = NBTECF + 1
         NOTECF( NBTECF ) = NTE
 30   ENDDO
      NBTECF0 = NBTECF

C     CALCUL DE LA QUALITE MINIMALE DES NBTECF TETRAEDRES DE NOTETR
C     POUR MISE AU POINT. A SUPPRIMER ENSUITE
      CALL QUVOTETD( QTEAME, PTXYZD, NBTECF, NOTECF, NOTETR,
     %               QUAMIN, QUAMOY, NBTEQM )

C     CONSTRUCTION NFETOI DES TRIANGLES FACES SIMPLES DES TETRAEDRES
C     ENVELOPPANTS. LES FACES VUES 2 FOIS SONT ELIMINEES
      CALL CRFETOI1( NBTECF, NOTECF, NOTETR,
     %               MXFETO, N1FEOC, N1FEVI, NFETOI )

C     TRACE DES FACES SIMPLES DE L'ETOILE DU CF
      KTITRE ='cr4etoil:           TETRAEDRES et FACES SIMPLES'
      WRITE( KTITRE(11:15), '(I5)' ) NBTECF
      CALL TRFETO4( KTITRE, PTXYZD, N1FEOC, NFETOI,
     %              NBTRCF, NOTRCF, LEFACO, NO0FAR,
     %              NBTECF, NOTECF, NOTETR )

 41   NF1 = N1FEOC
 42   IF( NF1 .GT. 0 ) THEN

C           NUMERO NOTETR DU TETRAEDRE DE LA FACE SIMPLE
            NTE = NFETOI(1,NF1)
C           NUMERO DE LA FACE SIMPLE DANS LE TETRAEDRE NTE
            I = ABS( NFETOI(2,NF1) )

C           SON VOLUME ET SA QUALITE
            CALL QUATETD( PTXYZD(1,NOTETR(1,NTE)),
     %                    PTXYZD(1,NOTETR(2,NTE)),
     %                    PTXYZD(1,NOTETR(3,NTE)),
     %                    PTXYZD(1,NOTETR(4,NTE)),
     %                    ARMIN, ARMAX, SURFTR, VOLUTE, QUALTE )

ccc            IF( QUALTE .LE. QTEAME ) THEN
ccc               PRINT*,'cr4etoil: TETRA(',NTE,')=',
ccc     %                 (NOTETR(L,NTE),L=1,8),' V=',VOLUTE,' Q=',QUALTE
ccc            ENDIF

C           TETRAEDRE OPPOSE PAR LA FACE I AU TETRAEDRE NTE
            NTEOP = NOTETR(4+I,NTE)
            IF( NTEOP .LE. 0 ) GOTO 48
            IF( NOTETR(1,NTEOP) .EQ. 0 ) GOTO 48

C           EXISTE T IL UNE AUTRE FACE DE L'ETOILE AYANT CE MEME
C           TETRAEDRE OPPOSE?
            NF2 = NFETOI(5,NF1)
 45         IF( NF2 .GT. 0 ) THEN

C              FACE I2 DU TETRAEDRE NTE2 DE LA FACE SIMPLE NF2
               NTE2 = NFETOI(1,NF2)
               I2   = ABS( NFETOI(2,NF2) )

C              TETRAEDRE OPPOSE PAR LA FACE I2 AU TETRAEDRE NTE2
               NTEOP2 = NOTETR(4+I2,NTE2)

               IF( NTEOP2 .EQ. NTEOP .AND. NTE2 .NE. NTE ) THEN

C                 RECHERCHE DE L'ARETE COMMUNE AUX TETRAEDRES NTE ET NTE2
                  DO J1=1,3
C                    LES 2 SOMMETS DE L'ARETE J DE LA FACE I DE NTE
                     IF( J1 .EQ. 3 ) THEN
                        J2=1
                     ELSE
                        J2 = J1+1
                     ENDIF
                     NS1 = NOTETR( NOSOFATE(J1,I), NTE )
                     NS2 = NOTETR( NOSOFATE(J2,I), NTE )

                     DO JJ1=1,3
C                       LES 2 SOMMETS DE L'ARETE JJ1 DE LA FACE I2 DE NTE2
                        IF( JJ1 .EQ. 3 ) THEN
                           JJ2=1
                        ELSE
                           JJ2 = JJ1+1
                        ENDIF
                        NSA1 = NOTETR( NOSOFATE(JJ1,I2), NTE2 )
                        NSA2 = NOTETR( NOSOFATE(JJ2,I2), NTE2 )

                        IF( (NS1.EQ.NSA1 .AND. NS2.EQ.NSA2) .OR.
     %                      (NS1.EQ.NSA2 .AND. NS2.EQ.NSA1) ) THEN

C                          NS1-NS2 EST ELLE UNE ARETE SIMPLE DU CF?
                           CALL ARARCF( NS1, NS2, NBCF, N1ARCF, NOARCF,
     %                                  NFLCF, NAVANT, NA1 )
                           IF( NA1 .NE. 0 ) THEN
C                             OUI: NS1-NS2 ARETE DU CF
C                              => TETRAEDRE NTEOP NON AJOUTABLE
                              GOTO 47
                           ENDIF

C                          RECHERCHE DU 3-EME SOMMET DE LA FACE DE NTE
                           IF( J2 .EQ. 3 ) THEN
                              J3=1
                           ELSE
                              J3=J2+1
                           ENDIF
                           NS3 = NOTETR( NOSOFATE(J3,I), NTE )

C                          RECHERCHE DU 3-EME SOMMET DE LA FACE DE NTE2
                           IF( JJ2 .EQ. 3 ) THEN
                              JJ3=1
                           ELSE
                              JJ3=JJ2+1
                           ENDIF
                           NS4 = NOTETR( NOSOFATE(JJ3,I2), NTE2 )

C                          CALCUL DE L'ANGLE AUTOUR DE L'ARETE NS1-NS2
C                          ANGLE: ANGLE ENTRE LES PLANS DES 2 TRIANGLES
C                                 S1S2S3 et S2S1S4
C                                 DANS L'INTERVALLE [0, 2Pi] RADIANS
                           CALL ANG2TR3D( PTXYZD(1,NS1),
     %                                    PTXYZD(1,NS2),
     %                                    PTXYZD(1,NS3),
     %                                    PTXYZD(1,NS4),
     %                                    ANGLE, IERR )
                           IF( IERR .NE. 0 ) THEN
                              ANGLE = 0D0
                           ELSE
                              ANGLE = ANGLE * 45D0 / PIS4
                           ENDIF

ccc                           IF( ANGLE .GT. 520D0 ) THEN
ccc                           IF( ANGLE .GT. 89D0 ) THEN
                           IF( ANGLE .GT. 200D0 ) THEN

C                             AJOUT DU TETRAEDRE NTEOP
                              PRINT*,'cr4etoil: ANGLE=',ANGLE,' degres',
     %                               ' face1:',NS2,NS1,NS3,
     %                               ' face2:',NS2,NS1,NS4
                              GOTO 46

                           ENDIF

                        ENDIF

                     ENDDO
                  ENDDO
                  GOTO 47

C                 AJOUT DU TETRAEDRE NTEOP COMMUN AUX 2 FACES A L'ETOILE
 46               IF( NBTECF .GE. MXTECF ) THEN
C                    ABANDON
                     GOTO 9000
                  ENDIF
                  NBTECF = NBTECF + 1
                  NOTECF( NBTECF ) = NTEOP
                  PRINT*,'cr4etoil: ajout ',NBTECF,
     %                   ' du tetraedre NOTETR(',nteop,')=',
     %                   (notetr(kk,nteop),kk=1,8)

                  DO J=1,4
C                    LA FACE J DU TETRAEDRE NTEOP EST
C                    SOIT RETIREE, SOIT AJOUTEE EN FIN DE CHAINAGE
                     CALL AJFAET1( NTEOP, J, NOTETR,
     %                             N1FEOC, N1FEVI, NFETOI, NFS )
                     IF( N1FEVI .EQ. -1 ) GOTO 9999
                  ENDDO
                  GOTO 41

               ENDIF

C              LA FACE SUIVANTE DE L'ETOILE DE MEME TETRAEDRE A RETROUVER
 47            NF2 = NFETOI(5,NF2)
               GOTO 45
            ENDIF

C           LA FACE SUIVANTE DE L'ETOILE A ANALYSER
 48         NF1 = NFETOI(5,NF1)
            GOTO 42

      ENDIF

C     UNE SEULE FOIS UN TETRAEDRE DANS NOTECF
      CALL UNITABL( NOTECF, NBTECF )

      IF( NBTECF0 .NE. NBTECF ) THEN
C        TRACE DES FACES SIMPLES DE L'ETOILE DU CF NFETOI VERSION 1
         KTITRE = 'cr4etoil:           TETRAEDRES ACTUELS'
         WRITE( KTITRE(11:15), '(I5)' ) NBTECF
         CALL SANSDBL( KTITRE, L )
         PRINT*,KTITRE(1:L)
         CALL TRFETO4( KTITRE, PTXYZD, N1FEOC, NFETOI,
     %                 NBTRCF, NOTRCF, LEFACO, NO0FAR,
     %                 NBTECF, NOTECF, NOTETR )
      ELSE
         PRINT*,'cr4etoil: PAS DE TETRAEDRE AJOUTE'
      ENDIF

cccC        ------------------------------------------------------
cccC        11') RETRAIT DES TETRAEDRES DE L'ETOILE
cccC             DONT LA PROJECTION DU BARYCENTRE SUR LE PLAN DE
cccC             CHAQUE TRIANGLE DU CF N'EST INTERNE A AUCUN D'EUX
cccC        ------------------------------------------------------
ccc         CALL SUTEBACF( KTITRE, PTXYZD, LEFACO, NO0FAR,
ccc     %                  NBTRCF, NOTRCF, NBSTCF, NOSTCF, NOTETR,
ccc     %                  NBTECF, NOTECF )

      GOTO 9999

C     PROBLEME A TRAITER
 9000 PRINT*,'cr4etoil: AUGMENTER MXTECF=',MXTECF

 9999 RETURN
      END
