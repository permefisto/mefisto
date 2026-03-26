      SUBROUTINE ARARCF1( NS1, NS2, NBTRCF, NOTRCF, LEFACO, NO0FAR,
     %                    NFLEFA, NARETE )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    L'ARETE DE SOMMET NS1-NS2 EST ELLE UNE ARETE DES NBTRCF
C -----    FACES NOTRCF DE LEFACO ou NO0FAR
C ENTREES:
C --------
C NS1 NS2: NUMERO DES 2 SOMMETS DE L'ARETE A TESTER
C NBTRCF : NOMBRE DE FACES DE NOTRCF
C NOTRCF : NUMERO DANS LEFACO DES NBTRCF FACES
C          >0 FACE DE LEFACO
C          <0 FACE DE NO0FAR
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
C NO0FAR : NUMERO DES 3 SOMMETS DE LA FACE AJOUTEE AU CF

C SORTIES:
C --------
C NFLEFA : >0 NUMERO LEFACO DE LA PREMIERE FACE TROUVEE D'ARETE NS1-NS2
C          <0 NUMERO NO0FAR DE LA PREMIERE FACE TROUVEE D'ARETE NS1-NS2
C          0 SI NS1-NS2 N'EST PAS UNE ARETE DE NOTCF DE LEFACO
C NARETE : NUMERO DE 1 A 3 DE L'ARETE NS1-NS2 DE NFLEFA DE LEFACO
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  VEULETTES sur MER               Janvier 2018
C2345X7..............................................................012
      INTEGER    LEFACO(11,0:*), NOTRCF(NBTRCF), NO0FAR(3,*)
      DO K = 1, NBTRCF

C        LE TRIANGLE K DE LEFACO
         NFLEFA = NOTRCF( K )

         DO NARETE=1,3

C              L'ARETE NARETE DE NFLEFA
               IF( NFLEFA .GT. 0 ) THEN

C                 FACE DE LEFACO
                  NSL1 = LEFACO(NARETE,NFLEFA)
                  IF( NARETE .EQ. 3 ) THEN
                     NAR = 1
                  ELSE
                     NAR=NARETE+1
                  ENDIF
                  NSL2 = LEFACO(NAR,NFLEFA)

               ELSE

C                 FACE NO0FAR
                  NF = -NFLEFA
                  NSL1 = NO0FAR(NARETE,NF)
                  IF( NARETE .EQ. 3 ) THEN
                     NAR = 1
                  ELSE
                     NAR=NARETE+1
                  ENDIF
                  NSL2 = NO0FAR(NAR,NF)

               ENDIF

               IF( ( NS1 .EQ. NSL2 .AND. NS2 .EQ. NSL1 ) .OR.
     %             ( NS1 .EQ. NSL1 .AND. NS2 .EQ. NSL2 ) ) THEN

C                 ARETE RETROUVEE
                  GOTO 9999

               ENDIF

         ENDDO

      ENDDO

C     ARETE NS1-NS2 NON RETROUVEE
      NFLEFA = 0
      NARETE = 0

 9999 RETURN
      END
