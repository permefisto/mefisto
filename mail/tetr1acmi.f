      SUBROUTINE TETR1ACMI( NS1, NS2, NS3,  N1TETS, NOTETR, PTXYZD,
     %                      NBSTCF, NOSTCF,
     %                      NBTE1A, MXTE1A, NOTE1A, NTECMIN, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   AJOUTER LE TETRAEDRE NTECMIN DE SOMMETS NS1 et NS2 ET DONT LES
C -----   2 FACES FORMENT UN ANGLE DE COSINUS MINIMAL AVEC LE
C         TRIANGLE NS1 NS2 NS3
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
C MODIFIES:
C ---------
C ATTENTION: LE PREMIER APPEL A FAIRE AVEC NBTE1A=0
C NBTE1A : EN ENTREE NOMBRE DE TETRAEDRES RANGES DANS NOTE1A
C          EN SORTIE AJOUT DES TETRAEDRES AYANT NS1 et NS2 COMME SOMMETS
C          = VALEUR EN ENTREE SI LE SOMMET NS2 APPARTIENT A AUCUN
C            DES TETRAEDRES DE SOMMET NS1
C            I.E. L'ARETE NS1-NS2 N'EST PAS DANS LA TETRAEDRISATION
C          > 0 SINON
C NOTE1A : NOTE1A(I) = NUMERO DU TETRAEDRE I DE SOMMETS NS1 et NS2

C SORTIES:
C --------
C NTECMIN: = 0 SI PAS DE TETRAEDRE
C          > 0 NUMERO NOTETR DU TETRAEDRE DE SOMMETS NS1 et NS2 ET DONT LES
C              2 FACES FORMENT UN ANGLE MINIMAL AVEC LE TRIANGLE NS1 NS2 NS3
C IERR   : = 0 PAS D'ERREUR
C          = 1 SI L'ARETE NS1-NS2 N'A PAS DE TETRAEDRES ENROULANTS
C          = 2 SI SATURATION DU TABLEAU NOTE1A
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC Saint PIERRE DU PERRAY Mars    2017
C MODIFS : ALAIN PERRONNET           Saint PIERRE DU PERRAY Janvier 2018
C....................................................................012
      INTEGER           NOTETR(8,*), N1TETS(*), NOTE1A(MXTE1A),
     %                  NOSTCF(NBSTCF)
      DOUBLE PRECISION  PTXYZD(4,*)
      DOUBLE PRECISION  COSMIN, COS2FA(2), COSFMI(2), ANGLE, ANGMIN,
     %                  PIS4, DATAN
      INTEGER           NFA(2), NFAMI(2), NOSOFA(3,4)
      DATA              NOSOFA  / 1,3,2,  2,3,4,  3,1,4,  4,1,2 /
C     => LA NORMALE A LA FACE EST ORIENTEE VERS L'EXTERIEUR DU TETRAEDRE

      NBTE1A0 = NBTE1A
      NTECMIN = 0
      PIS4 = DATAN( 1D0 )

C     AJOUT DES TETRAEDRES NON DEJA DANS NOTE1A D'ARETE NS1 NS2
      CALL TETR1A( NS1,   NS2,   N1TETS, NOTETR,
     %             NBTEA, MXTE1A-NBTE1A, NOTE1A(NBTE1A+1), IERR )
      NBTE1A = NBTE1A + NBTEA
C     NBTE1A: EN ENTREE NOMBRE DE TETRAEDRES DEJA RANGES DANS NOTE1A
C             EN SORTIE AJOUT DES TETRAEDRES AYANT NS1 NS2 COMME SOMMETS
      IF( IERR .NE. 0 ) GOTO 9999

      IF( NBTE1A .EQ. NBTE1A0 ) RETURN

C     UNICITE DES NO DE TETRAEDRES DU TABLEAU NOTE1A
      CALL UNITABL( NOTE1A, NBTE1A )

C     RECHERCHE DES TETRAEDRES D'ARETE NS1 NS2 ET FAISANT UN ANGLE MIN
C     AVEC LE TRIANGLE NS1 NS2 NS3

C     RECHERCHE DE L'ANGLE DES FACES D'ARETE NS1-NS2 DES TETRAEDRES
      COSMIN = 2D0

      DO NT = NBTE1A0+1, NBTE1A

C        NUMERO NOTETR DU TETRAEDRE NT D'ARETE NS1 NS2
         NTE = NOTE1A( NT )

C        L'UN DES SOMMETS DE NTE, AUTRE QUE NS1 NS2, EST IL SOMMET DU CF?
         DO N1=1,4
            NS = NOTETR(N1,NTE)
            IF( NS .NE. NS1 .AND. NS .NE. NS2 ) THEN
               DO N2=1,NBSTCF
                  IF( NS .EQ. NOSTCF(N2) ) THEN
C                    OUI: NTE TETRAEDRE AVEC 3 SOMMETS DU CF EST AJOUTE
                     NTECMIN = NTE
ccc                     print*,'tetr1acmi: NTECMIN=',NTECMIN,
ccc     %                      ' avec 3 SOMMETS du CF', NS1, NS2, NS,
ccc     %                      (NOTETR(NN,NTECMIN),NN=1,8)
                     GOTO 1000
                  ENDIF
               ENDDO
            ENDIF
         ENDDO

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

C           NS4 NUMERO DU SOMMET DE LA FACE NFA(MM)) DE NTE
C               DIFFERENT DE NS1 ET NS2
            DO NN=1,3
               NS4 = NOTETR( NOSOFA(NN,NFA(MM)), NTE )
               IF( NS4 .NE. NS1 .AND. NS4 .NE. NS2 ) THEN
                  GOTO 60
               ENDIF
            ENDDO

C           COSINUS DE L'ANGLE ENTRE LES NORMALES AUX 2 FACES
C           ORIENTEES DANS LE MEME SENS
 60         CALL COS2TD( PTXYZD(1,NS1), PTXYZD(1,NS2), PTXYZD(1,NS3),
     %                   PTXYZD(1,NS1), PTXYZD(1,NS2), PTXYZD(1,NS4),
     %                   COS2FA( MM ), IERR1, IERR2 )

            IF( IERR1 .NE. 0 .OR. IERR2 .NE. 0 ) THEN
C              UNE FACE est REDUITE a UNE ARETE ou un SOMMET
C              => UNE NORMALE N'EST PAS CALCULABLE
               COS2FA( MM ) = 9D0
            ENDIF

            CALL ANG2TR3D( PTXYZD(1,NS1), PTXYZD(1,NS2), PTXYZD(1,NS3),
     %                     PTXYZD(1,NS4),  ANGLE, IERR2 )

            print*,'tetr1acmi: NTE=',NTE,' NS1=',NS1,' NS2=',NS2,
     %             ' NS3=',NS3,' NS4=',NS4,' COS2FA=',COS2FA( MM ),
     %             ' ANGLE=',ANGLE*45D0/PIS4

         ENDDO

         DO MM = 1, 2
            IF( COS2FA(MM) .LT. COSMIN ) THEN
               COSMIN   = COS2FA(MM)
               NTECMIN  = NTE
               NFAMI(1) = NFA(1)
               NFAMI(2) = NFA(2)
               COSFMI(1) = COS2FA(1)
               COSFMI(2) = COS2FA(2)
               ANGMIN = ANGLE
            ENDIF

         ENDDO

      ENDDO

      print*,'tetr1acmi: AJOUT NTECMIN=',NTECMIN,':',
     %       (NOTETR(NN,NTECMIN),NN=1,8),
     %       ' No Faces Min=',NFAMI,' COSMIN',COSFMI,' ANGMIN=',ANGMIN


C     LE TETRAEDRE NTECMIN EST AJOUTE ( ECRASE LES ANCIENNES VALEURS )
 1000 NBTE1A = NBTE1A0 + 1
      NOTE1A( NBTE1A ) = NTECMIN

cccC     DES TETRAEDRES D'ARETE NS1-NS2 SEULS LES 2 TETRAEDRES
cccC     DE FACES OFFRANT UN ANGLE MIN AVEC NS1 NS2 NS3 SONT GARDES
ccc      IF( COSFMI(1) .GE. COSFMI(2) ) THEN
ccc         MM = 1
ccc      ELSE
ccc         MM = 2
ccc      ENDIF

cccC     LE TETRAEDRE OPPOSE A L'AUTRE FACE DE NTECMIN EST IL A AJOUTER?
ccc      NTEOP = NOTETR( 4+NFAMI(MM), NTECMIN )
ccc      IF( NTEOP .GT. 0 ) THEN

cccC        NTEOP EST IL DEJA DANS NOTE1A?
ccc         DO NT = 1, NBTE1A
ccc            IF( NOTE1A( NT ) .EQ. NTEOP ) GOTO 9999
ccc         ENDDO

cccC        NON: LE TETRAEDRE OPPOSE A L'AUTRE FACE DE NTECMIN EST AJOUTE
ccc         NBTE1A = NBTE1A + 1
ccc         NOTE1A( NBTE1A ) = NTEOP

ccc      ENDIF


ccc      PRINT *,'tetr1acmi: +',NBTE1A-NBTE1A0,' TETRAEDRES DE SOMMETS',
ccc     %                       NS1,NS2,' COSFMI=',COSFMI
ccc      DO NT = NBTE1A0+1, NBTE1A
ccc         NTE = NOTE1A(NT)
ccc         PRINT *,'TETRAEDRE',NT,' NOTETR(',NTE,')=',
ccc     %           (NOTETR(NN,NTE),NN=1,8)
ccc      ENDDO

 9999 RETURN
      END
