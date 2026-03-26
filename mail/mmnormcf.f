      SUBROUTINE MMNORMCF( NBTRCF, NOTRCF, LEFACO )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    PERMUTER LES ARETES DES FACES NOTRCF POUR AVOIR UNE SEULE
C -----    DIRECTION DE LA SURFACE (VERS L'INTERIEUR D'UNE ETOILE)

C ENTREES:
C --------
C NBTRCF : NOMBRE DE TRIANGLES LEFACO
C NOTRCF : NUMERO DANS LEFACO DES NBTRCF TRIANGLES

C ENTREE et SORTIE:
C -----------------
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
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  Saint PIERRE du PERRAY          Janvier 2018
C2345X7..............................................................012
      INTEGER    LEFACO(11,0:*), NOTRCF(NBTRCF)

      DO K = 1, NBTRCF

C        LE TRIANGLE K DE LEFACO
         NT = NOTRCF(K)

         IF( NT .GT. 0 ) THEN

            DO 10 L=1,3

C              L'ARETE L DE NT
               NS1 = LEFACO(L,NT)
               IF( L .EQ. 3 ) THEN
                  LL = 1
               ELSE
                  LL=L+1
               ENDIF
               NS2 = LEFACO(LL,NT)

C              LE TRIANGLE OPPOSE A L'ARETE L
               NTOP = LEFACO(5+L,NT)
               IF( NTOP .GT. 0 ) THEN

C                 BOUCLE SUR LES ARETES DE NTOP
                  DO M=1,3

C                    L'ARETE M DE NTOP
                     NSS1 = LEFACO(M,NTOP)
                     IF( M .EQ. 3 ) THEN
                        MM = 1
                     ELSE
                        MM=M+1
                     ENDIF
                     NSS2 = LEFACO(MM,NTOP)

                     IF( NS1 .EQ. NSS2 .AND. NS2 .EQ. NSS1 ) THEN

C                       PAS DE CHANGEMENT DE SENS DE LEFACO(NTOP)
                        GOTO 10

                     ELSE IF( NS1 .EQ. NSS1 .AND. NS2 .EQ. NSS2 ) THEN

C                       CHANGEMENT DE SENS DE LEFACO(NTOP)
                        N              = LEFACO(2,NTOP)
                        LEFACO(2,NTOP) = LEFACO(3,NTOP)
                        LEFACO(3,NTOP) = N
                        N              = LEFACO(6,NTOP)
                        LEFACO(6,NTOP) = LEFACO(8,NTOP)
                        LEFACO(8,NTOP) = N
                        GOTO 10

                     ENDIF

                  ENDDO
               ENDIF

 10         ENDDO

         ENDIF

      ENDDO

      RETURN
      END
