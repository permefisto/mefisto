      SUBROUTINE TETR1RCMI( NS1, NS2, NS3,  NOTETR, PTXYZD,
     %                      NBSTCF, NOSTCF, NBTECF, NOTECF, NTCOSMIN )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   RETIRER LE TETRAEDRE NTCOSMIN DE SOMMETS NS1 et NS2 ET DONT
C -----   UNE DE SES 2 FACES FORMENT UN ANGLE DE COSINUS MINIMAL AVEC
C         LE TRIANGLE NS1 NS2 NS3

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
C NBTECF : INDICE DANS NOTECF DU DERNIER TETRAEDRE D'ARETE NS1-NS2
C MXTECF : NOMBRE D'ENTIERS DU TABLEAU NOTECF

C MODIFIE:
C --------
C NOTECF : NOTECF(I) = NUMERO DU TETRAEDRE I
C          NOTECF( NTCOSMIN ) = -NTECMIN RENDU <0 POUR LE TETRAEDRE Cos MIN

C SORTIE :
C --------
C NTCOSMIN: = 0 SI PAS DE TEL TETRAEDRE DANS NOTECF
C           > 0 NUMERO NOTECF DU TETRAEDRE DE SOMMETS NS1 et NS2 ET DONT
C               UNE DES 2 FACES FORMENT UN ANGLE DE COSINUS MINIMAL AVEC
C               LE TRIANGLE NS1 NS2 NS3
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET           Saint PIERRE DU PERRAY     mai 2018
C....................................................................012
      INTEGER           NOTETR(8,*), NOTECF(NBTECF), NOSTCF(NBSTCF)
      DOUBLE PRECISION  PTXYZD(4,*)
      DOUBLE PRECISION  COSMIN, COS2FA(2), COSFMI(2), PIS4, DATAN
      INTEGER           NFA(2), NFAMI(2), NOSOFA(3,4)
      DATA              NOSOFA  / 1,3,2,  2,3,4,  3,1,4,  4,1,2 /
C     => LA NORMALE A LA FACE EST ORIENTEE VERS L'EXTERIEUR DU TETRAEDRE

      NTCOSMIN = 0
      NTECMIN  = 0
      PIS4 = DATAN( 1D0 )

C     RECHERCHE DU TETRAEDRE D'ARETE NS1 NS2 ET FAISANT UN ANGLE
C     DE COSINUS MINIMAL AVEC LE TRIANGLE NS1 NS2 NS3

C     RECHERCHE DE L'ANGLE DES FACES D'ARETE NS1-NS2 DES TETRAEDRES
      COSMIN = 2D0

      DO 80 NT = 1, NBTECF

C        NUMERO NOTETR DU TETRAEDRE NT D'ARETE NS1 NS2
         NTE = NOTECF( NT )
         IF( NTE .LE. 0 ) GOTO 80
         IF( NOTETR(1,NTE) .LE. 0 ) GOTO 80

C        NTE TETRAEDRE DE SOMMETS NS1 NS2?
         DO N1=1,4
            NS = NOTETR(N1,NTE)
            IF( NS .EQ. NS1 ) GOTO 10
         ENDDO
         GOTO 80

 10      DO N2=1,4
            NS = NOTETR(N2,NTE)
            IF( NS .EQ. NS2 ) GOTO 20
         ENDDO
         GOTO 80

C        L'UN DES SOMMETS DE NTE, AUTRE QUE NS1 NS2, EST IL SOMMET DU CF?
C        SI OUI, LE TETRAEDRE NE DOIT PAS ETRE RETIRE ET N'EST PAS CONSIDERE
 20      DO M1=1,4
            NS = NOTETR(M1,NTE)
            IF( NS .NE. NS1 .AND. NS .NE. NS2 ) THEN
               DO M2=1,NBSTCF
                  IF( NS .EQ. NOSTCF(M2) ) THEN
C                    TETRAEDRE AVEC 3 SOMMETS DU CF A NE PAS PRENDRE EN COMPTE
                     GOTO 80
                  ENDIF
               ENDDO
            ENDIF
         ENDDO

C        NUMERO DANS NTE DES 2 FACES D'ARETE NS1-NS2
         IF( N1 .GT. N2 ) THEN
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
            NOFACLOC = NFA(MM)

            DO NN=1,3
               NS4 = NOTETR( NOSOFA(NN,NOFACLOC), NTE )
               IF( NS4 .NE. NS1 .AND. NS4 .NE. NS2 ) THEN
                  GOTO 60
               ENDIF
            ENDDO

 60         IF( NS4 .EQ. 0 ) THEN
               print*,'tetr1rcmi: PROBLEME NTE=',NTE,' NS1=',NS1,
     %                ' NS2=',NS2,' NS3=',NS3,' NS4=',NS4,'?...'
               NOTECF( NT ) = -NTE
               GOTO 80
            ENDIF

C           COSINUS DE L'ANGLE ENTRE LES NORMALES AUX 2 FACES
C           ORIENTEES DANS LE MEME SENS
            CALL COS2TD( PTXYZD(1,NS1), PTXYZD(1,NS2), PTXYZD(1,NS3),
     %                   PTXYZD(1,NS1), PTXYZD(1,NS2), PTXYZD(1,NS4),
     %                   COS2FA( MM ), IERR1, IERR2 )
            IF( IERR1 .NE. 0 .OR. IERR2 .NE. 0 ) THEN
               COS2FA( MM ) = 9D0
            ENDIF

ccc            LA SECONDE NORMALE EST DIRIGEE DANS LE SENS OPPOSE / cos2td !
ccc            CALL ANG2TR3D( PTXYZD(1,NS1), PTXYZD(1,NS2), PTXYZD(1,NS3),
ccc     %                     PTXYZD(1,NS4),  ANGLE, IERR2 )

ccc            print*,'tetr1rcmi: NTE=',NTE,' NS1=',NS1,' NS2=',NS2,
ccc     %             ' NS3=',NS3,' NS4=',NS4,' COS2FA=',COS2FA( MM )
cccccc     %            ,' ANGLE=',ANGLE*45D0/PIS4

         ENDDO

         DO MM = 1, 2
            IF( COS2FA(MM) .LT. COSMIN ) THEN
               COSMIN   = COS2FA(MM)
               NTCOSMIN = NT
               NTECMIN  = NTE
               NFAMI(1) = NFA(1)
               NFAMI(2) = NFA(2)
               COSFMI(1) = COS2FA(1)
               COSFMI(2) = COS2FA(2)
            ENDIF
         ENDDO

 80   ENDDO

      IF( NTECMIN .GT. 0 ) THEN

C        SUPPRESSION DU TETRAEDRE DE COS MIN et 2 TETRAEDRES OPPOSES
C        SONT RETIRES DE L'ETOILE DES TETRAEDRES
C        -----------------------------------------------------------
ccc         print*,'tetr1rcmi: Retrait NTeCosMIN=',NTECMIN,':',
ccc     %          (NOTETR(NN,NTECMIN),NN=1,8),
ccc     %          ' NO Faces Min=',NFAMI,' CosMIN',COSFMI
C        LE TETRAEDRE NTECMIN EST RETIRE EN RENDANT SON NUMERO NOTECF NEGATIF
         NOTECF( NTCOSMIN ) = -NTECMIN

         DO 85 MM=1,2

C           LA FACE NFA(MM)) DE NTECMIN D'ARETE NS1 ET NS2
            NOFACLOC = NFA(MM)

C           LE TETRAEDRE OPPOSE A LA FACE NFAMI(MM)
            NTEOP = NOTETR( 4+NOFACLOC, NTECMIN )

            IF( NTEOP .GT. 0 ) THEN
C              RECHERCHE DE NTEOP DANS NOTECF
               DO K=1,NBTECF
                  IF( NOTECF(K) .EQ. NTEOP ) THEN
C                    LE TETRAEDRE NTEOP EST RETIRE EN RENDANT
C                    SON NUMERO NOTECF NEGATIF
                     NOTECF( K ) = -NTEOP
                     print*,'tetr1rcmi: Retrait NTEOP=',NTEOP,':',
     %                      (NOTETR(NN,NTEOP),NN=1,8)
                     GOTO 85
                  ENDIF
               ENDDO
            ENDIF

 85      ENDDO

      ENDIF

      RETURN
      END
