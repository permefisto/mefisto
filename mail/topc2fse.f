      SUBROUTINE TOPC2FSE( PTXYZD, NBCF,   N1ARCF, NOARCF,
     %                     MXFETO, N1FEOC, N1FEVI, NFETOI,
     %                     MXTECF, NBTECF, NOTECF, NOTETR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    RECHERCHER LES COUPLES DE FACES SIMPLES DE L'ETOILE NFETOI
C -----    D'ARETE COMMUNE NON DU CF ET AYANT UN TETRAEDRE OPPOSE
C          IDENTIQUE ET L'AJOUTER AU TABLEAU NOTECF
cccC          SI AU MOINS UN DES TETRAEDRES EST DE QUALITE MEDIOCRE

C ENTREES:
C --------
C PTXYZD : TABLEAU DES COORDONNEES DES POINTS
C          PAR POINT : X  Y  Z DISTANCE_SOUHAITEE
C NBCF   : NOMBRE DE LIGNES FERMEES PERIPHERIQUES DES FACES PERDUES
C N1ARCF : NUMERO DU 1-ER SOMMET OU 1-ERE ARETE DES LIGNES DU CONTOUR FERME
C          N1ARCF(0) POINTEUR SUR LA PREMIERE ARETE VIDE
C NOARCF : 1:NUMERO DU SOMMET DE L'ARETE DE LA LIGNE DU CONTOUR FERME
C          2:NUMERO DANS NOARCF DE L'ARETE SUIVANTE DU CF
C          3:NUMERO DANS LEFACO DU TRIANGLE ADJACENT OPPOSE A L'ARETE
C MXFETO : NOMBRE MAXIMAL DE FACES SIMPLES DECLARABLES DANS NFETOI

C MXTECF : NOMBRE MAXIMAL DE NUMERO NOTETR DE TETRAEDRES DECLARABLES
C          DANS NOTECF
C NBTECF : NOMBRE DE NUMEROS DE TETRAEDRES DE NOTECF
C NOTETR : LISTE DES TETRAEDRES
C          SOMMET1,    SOMMET2,    SOMMET3,    SOMMET4,
C          TETRAEDRE1, TETRAEDRE2, TETRAEDRE3, TETRAEDRE4
C          DE L'AUTRE COTE DE LA FACE
C          1: 123      2: 234      3: 341      4: 412

C SORTIES:
C --------
C NBTECF : NOMBRE DE NUMEROS DE TETRAEDRES DE NOTECF PEUT ETRE MODIFIE
C NOTECF : LES NBTECF NUMERO NOTETR DES TETRAEDRES DE L'ETOILE
C N1FEOC : POINTEUR SUR LA PREMIERE FACE SIMPLE DE L'ETOILE
C          CHAINAGE SUIVANT DANS NFETOI(5,*)
C N1FEVI : POINTEUR SUR LA PREMIERE FACE VIDE DE L'ETOILE
C          CHAINAGE SUIVANT DANS NFETOI(5,*)
C NFETOI : EN SORTIE VERSION 2
C          1: > 0 NUMERO NOTETR DU TETRAEDRE AU DELA DE LA FACE
C             = 0 PAS DE TETRAEDRE AU DELA DE LA FACE (FRONTIERE)
C             =-1 TETRAEDRE OPPOSE INCONNU A RETROUVER (NON FRONTIERE)
C          2: NUMERO PTXYZD DU 1-ER SOMMET DE LA FACE
C          3: NUMERO PTXYZD DU 2-ME SOMMET DE LA FACE
C          4: NUMERO PTXYZD DU 3-ME SOMMET DE LA FACE
C             S1S2xS1S3 EST LA NORMALE DIRIGEE VERS L'INTERIEUR DE
C             L'ETOILE
C          5: CHAINAGE SUIVANT DES FACES OCCUPEES ET VIDES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET  St PIERRE du PERRAY                 Juin 2018
C23456...............................................................012
      DOUBLE PRECISION  PTXYZD(4,*)
      INTEGER           N1ARCF(0:NBCF), NOARCF(1:3,1:*),
     %                  NFETOI(5,MXFETO), NOTETR(8,*), NOTECF(MXTECF),
     %                  NFA(16), NAF(16), NOSOTR(3)
      DOUBLE PRECISION  ARMIN, ARMAX, SURFTR(4), V

      IF( NBTECF .LE. 0 ) GOTO 9999

C     CONSTRUCTION NFETOI VERSION 1 DES TRIANGLES FACES SIMPLES DES
C     TETRAEDRES ENVELOPPANTS. LES FACES VUES 2 FOIS SONT ELIMINEES
 5    CALL CRFETOI1( NBTECF, NOTECF, NOTETR,
     %               MXFETO, N1FEOC, N1FEVI, NFETOI )

C     MODIFIER LES VALEURS DU TABLEAU NFETOI DE L'ETOILE de la
C     Version1 en la Version2
      CALL V12NFETOI( N1FEOC, NFETOI, NOTETR, NBFETO )

C     PARCOURS DES ARETES DES FACES SIMPLES DE L'ETOILE NOTECF
      NF = N1FEOC

 10   IF( NF .GT. 0 ) THEN

         DO 30 NA = 1, 3

            NS1 = NFETOI( 1+NA, NF)
            IF( NA .EQ. 3 ) THEN
               NA2 = 1
            ELSE
               NA2 = NA + 1
            ENDIF
            NS2 = NFETOI( 1+NA2, NF)

C           NS1-NS2 EST ELLE UNE ARETE DU CF?
            CALL ARARCF( NS1, NS2, NBCF, N1ARCF, NOARCF,
     %                   NLFCF, NAVANT, NA1 )
            IF( NA1 .NE. 0 ) THEN
C              OUI: NS1-NS2 EST L'ARETE NA1 DE NOARCF
               GOTO 30
            ENDIF

C           RECHERCHE DES FACES SIMPLES D'ARETE COMMUNE NON ARETE DU CF
C           AYANT UN MEME TETRAEDRE OPPOSE
            CALL ARFNSIET( NS1,    NS2, N1FEOC, NFETOI, NOTETR,
     %                     NBFA12, NFA, NAF )

            IF( NBFA12 .GE. 2 ) THEN

               DO 20 K1=1,NBFA12
                  NF1 = NFA( K1 )
                  NTEOP1 = NFETOI( 1, NF1 )

                  DO K2=K1+1,NBFA12
                     NF2 = NFA( K2 )
                     NTEOP2 = NFETOI( 1, NF2 )

                     IF( NTEOP1 .EQ. NTEOP2 ) THEN

C                       QUALITE DU TETRAEDRE NTE1 OPPOSE A NTEOP1?
                        DO M=1,3
                           NOSOTR(M) = NFETOI(1+M,NF1)
                        ENDDO
                        CALL NUFATRTE( NOSOTR, NOTETR(1,NTEOP1), NUFA1 )
                        NTE1 = NOTETR( 1+NUFA1, NTEOP1 )

C                       VOLUME ET QUALITE DE NTE1
                        CALL QUATETD( PTXYZD(1,NOTETR(1,NTE1)),
     %                                PTXYZD(1,NOTETR(2,NTE1)),
     %                                PTXYZD(1,NOTETR(3,NTE1)),
     %                                PTXYZD(1,NOTETR(4,NTE1)),
     %                                ARMIN, ARMAX, SURFTR, V, Q1 )

C                       QUALITE DU TETRAEDRE NTE2 OPPOSE A NTEOP2?
                        DO M=1,3
                           NOSOTR(M) = NFETOI(1+M,NF2)
                        ENDDO
                        CALL NUFATRTE( NOSOTR, NOTETR(1,NTEOP2), NUFA2 )
                        NTE2 = NOTETR( 1+NUFA2, NTEOP2 )

C                       VOLUME ET QUALITE DE NTE2
                        CALL QUATETD( PTXYZD(1,NOTETR(1,NTE2)),
     %                                PTXYZD(1,NOTETR(2,NTE2)),
     %                                PTXYZD(1,NOTETR(3,NTE2)),
     %                                PTXYZD(1,NOTETR(4,NTE2)),
     %                                ARMIN, ARMAX, SURFTR, V, Q2 )

ccc                        IF( Q1 .LT. 1E-3 .OR. Q2 .LT. 1E-3 ) THEN
ccc                        IF( Q1 .LT. 0.02 .OR. Q2 .LT. 0.02 ) THEN

C                       L'UN AU MOINS DES TETRAEDRES NTE1 NTE2 EST
C                       DE MEDIOCRE QUALITE

C                       CE TETRAEDRE OPPOSE UNIQUE EST AJOUTE A NOTECF
C                       S'IL N'EST PAS DEJA DANS NOTECF
                        DO M=1,NBTECF
                           IF( NOTECF(M) .EQ. NTEOP1 ) GOTO 20
                        ENDDO
                        IF( NBTECF .GE. MXTECF ) THEN
                           PRINT*,'topc2fse: AUGMENTER MXTECF=',MXTECF
                           GOTO 9999
                        ENDIF
                        NBTECF = NBTECF + 1
                        NOTECF( NBTECF ) = NTEOP1
                         PRINT*,'topc2fse: AJOUT a NOTECF DU TETRAEDRE',
     %                           NTEOP1,' St:',(NOTETR(M,NTEOP1),M=1,8),
     %                           ' Q1=',Q1,' Q2=',Q2

C                        RECONSTRUCTION DES FACES SIMPLES DE L'ETOILE
                         GOTO 5

ccc                         ENDIF
                     ENDIF

                  ENDDO
 20            ENDDO

            ENDIF

C           PASSAGE A L'ARETE SUIVANTE
 30      ENDDO

C        PASSAGE A LA FACE SUIVANTE
         NF = NFETOI(5,NF)
         GOTO 10

      ENDIF

 9999 RETURN
      END
