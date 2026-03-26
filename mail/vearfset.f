      SUBROUTINE VEARFSET( NOTETR, N1FEOC, NFETOI, NBARNODO )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  RECENSER LES ARETES NON DOUBLES DES FACES SIMPLES D'UNE ETOILE
C -----
C ENTREES:
C --------
C NOTETR : LISTE DES TETRAEDRES
C          SOMMET1,    SOMMET2,    SOMMET3,    SOMMET4,
C          TETRAEDRE1, TETRAEDRE2, TETRAEDRE3, TETRAEDRE4
C          DE L'AUTRE COTE DE LA FACE
C          1: 123      2: 234      3: 341      4: 412
C N1FEOC : POINTEUR SUR LA PREMIERE FACE DE L'ETOILE
C          CHAINAGE SUIVANT DANS NFETOI(5,*)
C NFETOI : VERSION 1
C          1: NUMERO DU TETRAEDRE DANS NOTETR
C          2: NUMERO LOCAL AU TETRAEDRE DE LA FACE DE L'ETOILE
C             UN SIGNE NEGATIF INDIQUE UN TRAITEMENT EFFECTUE
C          3  NON UTILISE ICI
C          4: NUMERO DE CETTE FACE DANS LEFACO, 0 SI PAS DANS LEFACO
C          5: CHAINAGE SUIVANT DES FACES OCCUPEES ET VIDES
C SORTIE :
C --------
C NBARNODO: NOMBRE D'ARETES NON DOUBLE DES FACES SIMPLES DE L'ETOILE NFETOI
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Saint Pierre du Perray           Fevrier 2019
C23456...............................................................012
      PARAMETER (MXFAAR=8)
      INTEGER    NOTETR(8,*), NFETOI(5,*), NOSOAR1(3), NOSOAR2(3),
     %           NOFAARND(MXFAAR)

      INTEGER    NOSOFATE(3,4)
C     NO DES SOMMETS DES FACES POUR QUE VU DU SOMMET MANQUANT
C     LES SOMMETS SOIENT VUS DANS LE SENS DIRECT
      DATA              NOSOFATE/ 1,2,3, 2,4,3, 3,4,1, 4,2,1 /

C     NOMBRE D'ARETES NON DOUBLES DES FACES SIMPLES DE L'ETOILE
      NBARNODO = 0

      NF1 = N1FEOC

C     BOUCLE SUR LES FACES SIMPLES DE L'ETOILE
 10   IF( NF1 .GT. 0 ) THEN

C        NFETOI VERSION 1
C        LE NO DU TETRAEDRE NT1 INTERNE A L'ETOILE ET NO LOCAL DE LA FACE
         NT1   = NFETOI( 1, NF1 )
         NFNT1 = NFETOI( 2, NF1 )

C        NUMERO DES 3 SOMMETS DE LA FACE NFNT1 DU TETRAEDRE NT1
         NOSOAR1( 1 ) = NOTETR( NOSOFATE(1,NFNT1), NT1 )
         NOSOAR1( 2 ) = NOTETR( NOSOFATE(2,NFNT1), NT1 )
         NOSOAR1( 3 ) = NOTETR( NOSOFATE(3,NFNT1), NT1 )

         NS1 = NOSOAR1( 3 )
         DO I1=1,3

C           L'ARETE NS1-NS2 DE LA FACE NFNT1 DU TETRAEDRE NT1
            NS2 = NOSOAR1( I1 )

C           L'ARETE NS1-NS2 APPARTIENT A COMBIEN DE FACES SIMPLES?
            NBS1S2 = 0

            NF2 = N1FEOC

C           BOUCLE SUR LES ARETES DES FACES SIMPLES DE L'ETOILE
C           ---------------------------------------------------
 20         IF( NF2 .GT. 0 ) THEN

C              LE NO DU TETRAEDRE NT2 INTERNE A L'ETOILE ET NO LOCAL DE LA FACE
               NT2   = NFETOI( 1, NF2 )
               NFNT2 = NFETOI( 2, NF2 )

C              NUMERO DES 3 SOMMETS DE LA FACE NFNT2 DU TETRAEDRE NT2
               NOSOAR2( 1 ) = NOTETR( NOSOFATE(1,NFNT2), NT2 )
               NOSOAR2( 2 ) = NOTETR( NOSOFATE(2,NFNT2), NT2 )
               NOSOAR2( 3 ) = NOTETR( NOSOFATE(3,NFNT2), NT2 )

               NS3 = NOSOAR2( 3 )
               DO I2=1,3

C                 L'ARETE NS1-NS2 DE LA FACE NFNT2 DU TETRAEDRE NT2
                  NS4 = NOSOAR2( I2 )

                  IF( (NS1 .EQ. NS3 .AND. NS2 .EQ. NS4) .OR.
     %                (NS1 .EQ. NS4 .AND. NS2 .EQ. NS3)  ) THEN
                     NBS1S2 = NBS1S2 + 1
C                    NO DE LA FACE SIMPLE
                     NOFAARND( NBS1S2 ) = NF2
                     GOTO 30
                  ENDIF

                  NS3 = NS4
               ENDDO

C              PASSAGE A LA FACE SIMPLE NF2 SUIVANTE
 30            NF2 = NFETOI( 5, NF2 )
               GOTO 20

            ENDIF

C           BILAN DE L'ARETE NS1-NS2 NON DOUBLE DES FACES SIMPLES?
C           ------------------------------------------------------
            IF( NBS1S2 .NE. 2 ) THEN
C              UNE ARETE NON DOUBLE DES FACES SIMPLES
               NBARNODO = NBARNODO + 1
               PRINT*
               PRINT*,'vearfset: l''ARETE',NS1,'->',NS2,
     %                ' APPARTIENT AUX',NBS1S2,' FACES SIMPLES:',
     %                (NOFAARND(K),K=1,NBS1S2)
               DO K=1,NBS1S2
C                 AFFICHAGE DES TETRAEDRES AVEC CETTE ARETE NS1-NS2
                  NF2   = NOFAARND(K)
                  NT2   = NFETOI( 1, NF2 )
                  NFNT2 = NFETOI( 2, NF2 )
                  PRINT*,'vearfset: FACE',NFNT2,' NOTETR(',NT2,')=',
     %                   (NOTETR(kk,NT2),kk=1,8)
               ENDDO
            ENDIF

C           PASSAGE A L'ARETE SUIVANTE DE LA FACE SIMPLE NF1
            NS1 = NS2
         ENDDO

C        PASSAGE A LA FACE SIMPLE SUIVANTE
         NF1 = NFETOI( 5, NF1 )
         GOTO 10

      ENDIF

      IF( NBARNODO .GT. 0 ) THEN
         PRINT*,'vearfset:',NBARNODO,' ARETES NON DOUBLE des FACES SIMPL
     %ES de l''ETOILE NFETOI'
      ENDIF

      RETURN
      END
