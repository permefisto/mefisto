      SUBROUTINE TETR1ANMI( NS1,NS2, NS3,  N1TETS, NOTETR, PTXYZD,
     %                      NBTE1A, MXTE1A, NOTE1A, NTEMIN, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   TROUVER LE TETRAEDRE DE SOMMETS NS1 et NS2 ET DONT LES
C -----   2 FACES DE CETTE ARETE FORMENT UN ANGLE MINIMAL AVEC LE
C         TRIANGLE NS1 NS2 NS3 ET LE CHEVAUCHE
C
C ENTREES:
C --------
C NS1,NS2: NUMERO DES 2 SOMMETS DE L'ARETE COMMUNE A TOUS LES TETRAEDRES
C NS3    : NUMERO DU 3-EME SOMMET DEFINISSANT LE TRIANGLE
C N1TETS : N1TETS(I) NUMERO D'UN TETRAEDRE AYANT POUR SOMMET I
C NOTETR : LISTE DES TETRAEDRES
C          SOMMET1,    SOMMET2,    SOMMET3,    SOMMET4,
C          TETRAEDRE1, TETRAEDRE2, TETRAEDRE3, TETRAEDRE4
C          DE L'AUTRE COTE DE LA FACE DE NUMEROTATION
C          1: 123      2: 234      3: 341      4: 412
C PTXYZD : TABLEAU DES COORDONNEES DES POINTS
C          PAR POINT : X  Y  Z DISTANCE_SOUHAITEE
C MXTE1A : NOMBRE D'ENTIERS DU TABLEAU NOTE1A
C
C SORTIES:
C --------
C NBTE1A : >0 NOMBRE DE TETRAEDRES D'ARETE NS1-NS2
C          =0 PAS DE TETRAEDRES RETROUVES
C NOTE1A : NOTE1A(I) = NUMERO DU TETRAEDRE I DE SOMMETS NS1 et NS2
C NTEMIN : NUMERO NOTETR DU TETRAEDRE DE SOMMETS NS1 et NS2 ET DONT LES
C          2 FACES FORMENT UN ANGLE MINIMAL AVEC LE TRIANGLE NS1 NS2 NS3
C          = 0 SI NON TROUVE
C IERR   : = 0 PAS D'ERREUR
C          = 1 SI L'ARETE NS1-NS2 N'A PAS DE TETRAEDRES ENROULANTS
C          = 2 SI SATURATION DU TABLEAU NOTE1A
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC et St PIERRE DU PERRAY Novembre 2017
C....................................................................012
      INTEGER           NOTETR(8,*),
     %                  N1TETS(*),
     %                  NOTE1A(MXTE1A)
      DOUBLE PRECISION  PTXYZD(4,*)
      DOUBLE PRECISION  COSMAX, COS2FA(2), COSFMX(2), VOLTET, V1, V2
      INTEGER           NFA(2), NFAMX(2), NS4(2), NOSOFA(3,4)
      DATA              NOSOFA  / 1,3,2,  2,3,4,  3,1,4,  4,1,2 /
C     => LA NORMALE A LA FACE EST ORIENTEE VERS L'EXTERIEUR DU TETRAEDRE

      NTEMIN = 0

      CALL TETR1A( NS1,    NS2,    N1TETS, NOTETR,
     %             NBTE1A, MXTE1A, NOTE1A, IERR )
C     NBTE1A: NOMBRE DE TETRAEDRES AYANT NS1 NS2 COMME SOMMETS
      IF( IERR .NE. 0 ) GOTO 9999

C     RECHERCHE DES TETRAEDRES D'ARETE NS1 NS2 ET FAISANT UN ANGLE MIN
C     AVEC LE TRIANGLE NS1 NS2 NS3

C     RECHERCHE DE L'ANGLE DES FACES D'ARETE NS1-NS2 DES TETRAEDRES
      COSMAX = -2D0

      DO NT = 1, NBTE1A

C        NUMERO NOTETR DU TETRAEDRE NT D'ARETE NS1 NS2
         NTE = NOTE1A( NT )

         DO N1=1,4
C           RECHERCHE DU NUMERO DE SOMMET DE NS1 DANS NTE
            IF( NS1 .EQ. NOTETR(N1,NTE) ) GOTO 2
         ENDDO

 2       DO N2=1,4
C           RECHERCHE DU NUMERO DE SOMMET DE NS2 DANS NTE
            IF( NS2 .EQ. NOTETR(N2,NTE) ) GOTO 3
         ENDDO

C        NUMERO DANS NTE DES 2 FACES D'ARETE NS1-NS2
 3       IF( N1 .GT. N2 ) THEN
C           PERMUTATION DES 2 NUMEROS POUR QUE N1<N2
            MM = N1
            N1 = N2
            N2 = MM
         ENDIF

C        NUMERO DES 2 FACES DU TETRAEDRE NTE D'ARETE N1<N2
         IF( N2 .NE. 4 ) THEN

C           ICI N2=2 ou 3 CAR N1<N2
            IF( N1 .EQ. 1 ) THEN
C              ARETE N1 => NUMERO DES 2 FACES DE NTE
               IF( N2 .EQ. 2 ) THEN
C                 ARETE N1=1 N2=2
                  NFA(1) = 1
                  NFA(2) = 4
               ELSE
C                 ARETE N1=1 N2=3
                  NFA(1) = 1
                  NFA(2) = 3
               ENDIF
            ELSE
C              ARETE N1=2 N2=3
               NFA(1) = 1
               NFA(2) = 2
            ENDIF
            GOTO 50

         ELSE

C           ICI N2=4 et N1=1 ou 2 ou 3
            GOTO( 44, 45, 46 ),N1
C           ARETE 3+N1 => NUMERO DES 2 FACES DE NTE

C           ARETE N1=1 N2=4
 44         NFA(1) = 3
            NFA(2) = 4
            GOTO 50

C           ARETE N1=2 N2=4
 45         NFA(1) = 2
            NFA(2) = 4
            GOTO 50

C           ARETE N1=3 N2=4
 46         NFA(1) = 2
            NFA(2) = 3
            GOTO 50

         ENDIF

C        PARCOURS DES 2 FACES DE NTE D'ARETE NS1-NS2
 50      DO MM=1,2

C           NS4(MM) NUMERO DU SOMMET DE LA FACE NFA(MM)) DE NTE
C                   DIFFERENT DE NS1 ET NS2
            DO NN=1,3
               NS = NOTETR( NOSOFA(NN,NFA(MM)), NTE )
               IF( NS .NE. NS1 .AND. NS .NE. NS2 ) THEN
                  GOTO 60
               ENDIF
            ENDDO

C           COSINUS DE L'ANGLE ENTRE LES NORMALES AUX 2 FACES
C           ORIENTEES DANS LE MEME SENS
 60         NS4( MM ) = NS
            CALL COS2TD( PTXYZD(1,NS1), PTXYZD(1,NS2), PTXYZD(1,NS3),
     %                   PTXYZD(1,NS1), PTXYZD(1,NS2), PTXYZD(1,NS),
     %                   COS2FA( MM ), IERR1, IERR2 )

            IF( IERR1 .NE. 0 .OR. IERR2 .NE. 0 ) THEN
C              UNE FACE est REDUITE a UNE ARETE ou un SOMMET
C              => UNE NORMALE N'EST PAS CALCULABLE
               COS2FA( MM ) = -9D0
            ENDIF

         ENDDO

         DO MM = 1, 2
            IF( COS2FA(MM) .GT. COSMAX ) THEN

C              LES SOMMETS NS4 SONT ILS DE CHAQUE COTE DU TRIANGLE NS1-NS2-NS3?
               V1 = VOLTET( PTXYZD(1,NS1), PTXYZD(1,NS2), PTXYZD(1,NS3),
     %                      PTXYZD(1,NS4(1)) )
               V2 = VOLTET( PTXYZD(1,NS1), PTXYZD(1,NS2), PTXYZD(1,NS3),
     %                      PTXYZD(1,NS4(2)) )
               IF( V1*V2 .LE. 0D0 ) THEN
                  COSMAX   = COS2FA(MM)
                  NTEMIN   = NTE
                  NFAMX(1) = NFA(1)
                  NFAMX(2) = NFA(2)
                  COSFMX(1) = COS2FA(1)
                  COSFMX(2) = COS2FA(2)
               ENDIF

            ENDIF
         ENDDO

      ENDDO

      IF( NTEMIN .GT. 0 ) THEN
         PRINT *,'tetr1anmi: Arete',NS1,NS2,' NS3=',NS3,
     %           ' NOTETR(',NTEMIN,')=',(NOTETR(NN,NTEMIN),NN=1,8),
     %           ' avec ANGLE MINIMUM de COSMAX=', COSMAX
      ENDIF

 9999 RETURN
      END
