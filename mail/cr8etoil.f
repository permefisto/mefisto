      SUBROUTINE CR8ETOIL( QUAMINEX, KTITRE, PTXYZD, NOTETR, N1TETS,
     %                     NFP,    NBFAPE, NOFAPE, 
     %                     MXFACO, LEFACO, NO0FAR,
     %                     NBTRCF, NOTRCF, NBSTCF, NOSTCF,
     %                     NBCF,   MXARCF, N1ARCF, NOARCF, 
     %                     MXFETO, N1FEOC, N1FEVI, NFETOI,
     %                     MXTECF, NBTECF, NOTECF )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CONSTRUCTION 8 DE L'ETOILE DES NBTRCF FACES PERDUES
C -----    DE LA FRONTIERE PAR AJOUT DES TETRAEDRES D'ARETE COMMUNE
C          N'APPARTENANT PAS AU CF, AYANT AU MOINS UN SOMMET DU CF ET
C          D'ANGLE DES 2 FACES SUPERIEUR A 210 DEGRES

C ENTREES:
C --------
C QUAMINEX:QUALITE MINIMALE EXIGEE POUR CREER UN TETRAEDRE
C KTITRE : TITRE D'UN TRACE
C PTXYZD : TABLEAU DES COORDONNEES DES POINTS
C          PAR POINT : X  Y  Z DISTANCE_SOUHAITEE
C NOTETR : LISTE DES TETRAEDRES
C          SOMMET1,    SOMMET2,    SOMMET3,    SOMMET4,
C          TETRAEDRE1, TETRAEDRE2, TETRAEDRE3, TETRAEDRE4
C          DE L'AUTRE COTE DE LA FACE
C N1TETS : N1TETS(NS) NUMERO D'UN TETRAEDRE AYANT POUR SOMMET NS

C NFP    : NUMERO NOFAPE DE LA FACE PERDUE LEFACO A TRAITER
C NBFAPE : NOMBRE DE FACES PERDUES DE LEFACO
C NOFAPE : NUMERO DES NBFAPE FACES PERDUES

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
C NBSTCF : NOMBRE DE SOMMETS DU CF STOCKES DANS NOSTCF
C NOSTCF : NUMERO DES NBSTCF SOMMETS DU CF DES FACES PERDUES

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
      COMMON / TRTETR / STOPTE, TRACTE
      LOGICAL           STOPTE, TRACTE, TRACTE0

      CHARACTER*(*)     KTITRE
      DOUBLE PRECISION  PTXYZD(4,*)
      DOUBLE PRECISION  ANGLE, PIS4, DATAN

      INTEGER           NOTETR(8,*), NOFAPE(NBFAPE), N1TETS(*),
     %                  LEFACO(11,0:MXFACO), NO0FAR(3,*),
     %                  NOTRCF(NBTRCF), NOSTCF(NBSTCF),
     %                  N1ARCF(0:MXARCF), NOARCF(3,MXARCF),
     %                  NFETOI(5,MXFETO), NOTECF(MXTECF)

      INTEGER           NOSOFATE(3,4)
      DATA              NOSOFATE / 1,2,3,  2,4,3,  3,4,1,  4,2,1 /
C     => LA NORMALE A LA FACE EST ORIENTEE VERS L'INTERIEUR DU TETRAEDRE

      NFLPER  = NOFAPE( NFP )
      TRACTE0 = TRACTE

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
      DO 5 K = 1, NBTECF0
         NTE = NOTECF( K )
         IF( NTE .LE. 0 ) GOTO 5
         IF( NTE .GT. 0 .AND. NOTETR(1,NTE) .EQ. 0 ) GOTO 5
         NBTECF = NBTECF + 1
         NOTECF( NBTECF ) = NTE
 5    ENDDO
      NBTECF0 = NBTECF

C     CALCUL DE LA QUALITE MINIMALE DES NBTECF TETRAEDRES DE NOTETR
C     POUR MISE AU POINT. a supprimer ensuite
 10   CALL QUVOTETD( QUAMINEX, PTXYZD, NBTECF, NOTECF, NOTETR,
     %               QUAMIN,   QUAMOY, NBTEQM )

C     CONSTRUCTION NFETOI DES TRIANGLES FACES SIMPLES DES TETRAEDRES
C     ENVELOPPANTS. LES FACES VUES 2 FOIS SONT ELIMINEES
      CALL CRFETOI1( NBTECF, NOTECF, NOTETR,
     %               MXFETO, N1FEOC, N1FEVI, NFETOI )

cccC     TRACE DES FACES SIMPLES DE L'ETOILE DU CF
ccc      KTITRE='cr8etoil: Debut          TETRAEDRES et FACES SIMPLES'
ccc      WRITE( KTITRE(17:21), '(I5)' ) NBTECF
ccc      CALL TRFETO4( KTITRE, PTXYZD, N1FEOC, NFETOI,
ccc     %              NBTRCF, NOTRCF, LEFACO, NO0FAR,
ccc     %              NBTECF, NOTECF, NOTETR )


C     PARCOURS DES ARETES DOUBLES DES FACES SIMPLES DE L'ETOILE
C     NON ARETE DU CF
      NF1 = N1FEOC
 30   IF( NF1 .GT. 0 ) THEN

C        NUMERO NOTETR DU TETRAEDRE DE LA FACE SIMPLE NF1
         NTE1 = NFETOI(1,NF1)
C        NUMERO DE LA FACE SIMPLE DANS LE TETRAEDRE NTE1
         I1 = ABS( NFETOI(2,NF1) )

C        PARCOURS DES 3 ARETES DE LA FACE I1 DE NTE1
         DO 60 J1=1,3

C           LES 2 SOMMETS DE L'ARETE J1 DE LA FACE I1 DE NTE1
            NS1 = NOTETR( NOSOFATE(J1,I1), NTE1 )
            IF( J1 .EQ. 3 ) THEN
               J2=1
            ELSE
               J2 = J1+1
            ENDIF
            NS2 = NOTETR( NOSOFATE(J2,I1), NTE1 )

C           NS1-NS2 EST ELLE UNE ARETE SIMPLE DU CF?
            CALL ARARCF( NS1, NS2, NBCF, N1ARCF, NOARCF,
     %                   NFLCF, NAVANT, NACF )
            IF( NACF .NE. 0 ) THEN
C              OUI: NS1-NS2 ARETE DU CF. A NE PAS REMPLIR DE TETRAEDRES
               GOTO 60
            ENDIF

            NBS = 0
            DO N=1,NBSTCF
               NSCF = NOSTCF( N )
               IF( NSCF .EQ. NS1 ) THEN
                  NBS = NBS + 1
               ENDIF
               IF( NSCF .EQ. NS2 ) THEN
                  NBS = NBS + 1
               ENDIF
            ENDDO

            IF( NBS .EQ. 0 ) GOTO 60

C           ARETE AYANT AU MOINS UN SOMMET DU CF MAIS NON ARETE SIMPLE DU CF

C           RECHERCHE DE LA FACE ADJACENTE NF2 PAR CETTE ARETE NS1-NS2
            NF2 = N1FEOC
 50         IF( NF2 .GT. 0 ) THEN
               IF( NF2 .NE. NF1 ) THEN

C                 LE NO DU TETRAEDRE NTE INTERNE A L'ETOILE
                  NTE2 = NFETOI(1,NF2)
C                 LE NO NOFA LOCAL DE LA FACE
                  I2 = ABS( NFETOI(2,NF2) )

C                 PARCOURS DES 3 ARETES DU TRIANGLE NOSOTR
                  DO K1=1,3

C                    NO DES 2 SOMMETS DE L'ARETE K1
                     NSA1 = NOTETR( NOSOFATE(K1,I2), NTE2 )
                     IF( K1 .EQ. 3 ) THEN
                        K2 = 1
                     ELSE
                        K2 = K1 + 1
                     ENDIF
                     NSA2 = NOTETR( NOSOFATE(K2,I2), NTE2 )

                     IF( ( NSA1 .EQ. NS2 .AND. NSA2 .EQ. NS1 ) .OR.
     %                   ( NSA1 .EQ. NS1 .AND. NSA2 .EQ. NS2 ) ) THEN

C                       L'ARETE NS1-NS2 RETROUVEE DANS LA FACE NF2
C                       EST COMMUNE AUX FACES NF1 et NF2
C                       CALCUL DE L'ANGLE ENTRE CES 2 FACES
C                       RECHERCHE DU 3-EME SOMMET DE LA FACE DE NTE1
                        IF( J2 .EQ. 3 ) THEN
                           J3=1
                        ELSE
                           J3=J2+1
                        ENDIF
                        NS3 = NOTETR( NOSOFATE(J3,I1), NTE1 )

C                       RECHERCHE DU 3-EME SOMMET DE LA FACE DE NTE2
                        IF( K2 .EQ. 3 ) THEN
                           K3=1
                        ELSE
                           K3=K2+1
                        ENDIF
                        NS4 = NOTETR( NOSOFATE(K3,I2), NTE2 )

C                       CALCUL DE L'ANGLE AUTOUR DE L'ARETE NS1-NS2
C                       ANGLE: ANGLE ENTRE LES PLANS DES 2 TRIANGLES
C                              S1S2S3 et S2S1S4
C                              DANS L'INTERVALLE [0, 2Pi] RADIANS
                        CALL ANG2TR3D( PTXYZD(1,NS1),
     %                                 PTXYZD(1,NS2),
     %                                 PTXYZD(1,NS3),
     %                                 PTXYZD(1,NS4),
     %                                 ANGLE, IERR )
                        IF( IERR .NE. 0 ) THEN
                           ANGLE = 0D0
                        ELSE
                           ANGLE = ANGLE * 45D0 / PIS4
                        ENDIF
                        IF( ANGLE .GT. 210D0 ) THEN

                           PRINT*,'cr8etoil: ANGLE=',ANGLE,' degres',
     %                            ' face1:',NS1,NS2,NS3,
     %                            ' face2:',NS2,NS1,NS4


C        TRACE DES FACES SIMPLES DE L'ETOILE DU CF NFETOI VERSION 1
C        AVANT AJOUT DES TETRAEDRES D'ARETE NS1-NS2
ccc         TRACTE = .TRUE.
         KTITRE='cr8etoil: avant          TETRAEDRES NS1=        NS2=   
     %       '
         WRITE( KTITRE(17:21), '(I5)' ) NBTECF
         WRITE( KTITRE(41:47), '(I7)' ) NS1
         WRITE( KTITRE(53:59), '(I7)' ) NS2
         CALL SANSDBL( KTITRE, L )
         PRINT*,KTITRE(1:L)
         CALL TRFETO4( KTITRE, PTXYZD, N1FEOC, NFETOI,
     %                 NBTRCF, NOTRCF, LEFACO, NO0FAR,
     %                 NBTECF, NOTECF, NOTETR )


C                          OUI: AJOUT DE TOUS LES TETRAEDRES D'ARETE NS1-NS2
                           NBTEC0 = NBTECF
                           CALL TETR1A( NS1,   NS2,   N1TETS, NOTETR,
     %                                  NBTEA, MXTECF-NBTECF,
     %                                  NOTECF(NBTECF+1), IERR )

                           IF( NBTEA .LE. 0 .OR. IERR .GT. 0 ) THEN
C                             AUCUN TETRAEDRE AJOUTE ou ERREUR
                              GOTO 60
                           ENDIF

                           NBTECF = NBTECF + NBTEA

C                          UNICITE DES NO DE TETRAEDRES DU TABLEAU NOTECF
                           CALL UNITABL( NOTECF, NBTECF )


C        TRACE DES FACES SIMPLES DE L'ETOILE DU CF NFETOI VERSION 1
C        APRES AJOUT DES TETRAEDRES D'ARETE NS1-NS2
                           tracte = .true.
         KTITRE='cr8etoil: apres          TETRAEDRES NS1=        NS2=   
     %      '
         WRITE( KTITRE(17:21), '(I5)' ) NBTECF
         WRITE( KTITRE(41:47), '(I7)' ) NS1
         WRITE( KTITRE(53:59), '(I7)' ) NS2
         CALL SANSDBL( KTITRE, L )
         PRINT*,KTITRE(1:L)
         CALL TRFETO4( KTITRE, PTXYZD, N1FEOC, NFETOI,
     %                 NBTRCF, NOTRCF, LEFACO, NO0FAR,
     %                 NBTECF, NOTECF, NOTETR )



                           IF( NBTEC0 .NE. NBTECF ) THEN
C                             DES TETRAEDRES ONT ETE AJOUTES
C                             RETOUR AU DEPART DU CALCUL
                              GOTO 10
                           ENDIF

                        ENDIF

                     ENDIF

                  ENDDO
               ENDIF

C              PASSAGE A LA FACE SUIVANTE
               NF2 = NFETOI(5,NF2)
               GOTO 50

            ENDIF

 60      ENDDO

C        LA FACE SUIVANTE DE L'ETOILE A ANALYSER
         NF1 = NFETOI(5,NF1)
         GOTO 30

      ENDIF

      IF( NBTECF0 .NE. NBTECF ) THEN
C        TRACE DES FACES SIMPLES DE L'ETOILE DU CF NFETOI VERSION 1
         KTITRE = 'cr8etoil: Fin avec          TETRAEDRES'
         WRITE( KTITRE(20:24), '(I5)' ) NBTECF
         CALL SANSDBL( KTITRE, L )
         PRINT*,KTITRE(1:L)
         CALL TRFETO4( KTITRE, PTXYZD, N1FEOC, NFETOI,
     %                 NBTRCF, NOTRCF, LEFACO, NO0FAR,
     %                 NBTECF, NOTECF, NOTETR )
      ELSE
         PRINT*,'cr8etoil: Fin AUCUN TETRAEDRE AJOUTE a NOTECF'
      ENDIF

      TRACTE = TRACTE0
      RETURN
      END
