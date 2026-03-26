      SUBROUTINE VETEVIDE( N1TEVI, MXTETR, NOTETR, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT: VERIFIER SI TOUT TETRAEDRE VIDE A SON 1-ER NUMERO NOTETR(1,*) NUL
C ---- VERIFIER LE CHAINAGE DES TETRAEDRES VIDES (Boucle infinie...)
C      EN CAS D'ERREUR, CORRECTION du CHAINAGE des TETRAEDRES VIDES

C ENTREES:
C --------
C N1TEVI : NO DU PREMIER TETRAEDRE VIDE DANS NOTETR
C MXTETR : NOMBRE DE TETRAEDRES DECLARABLES DANS NOTETR
C NOTETR : LISTE DES TETRAEDRES
C          SOMMET1,    SOMMET2,    SOMMET3,    SOMMET4,
C          TETRAEDRE1, TETRAEDRE2, TETRAEDRE3, TETRAEDRE4
C          DE L'AUTRE COTE DE LA FACE
C          1: 123      2: 234      3: 341      4: 412
C SORTIES:
C --------
C IERR   : 0 SI PAS D'ERREUR
C          1 SI UN NOTETR(1,NTVI)>0
C          2 SI BOUCLE INFINIE des TETRAEDRES VIDES
C          -> CORRECTION de N1TEVI et du CHAINAGE des TETRAEDRES VIDES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  Saint Pierre du Perray         Novembre 2019
C2345X7..............................................................012
      INTEGER  NOTETR(8,MXTETR)

C     LE PREMIER TETRAEDRE VIDE
      NTVI   = N1TEVI
      NBTEVI = 0

 10   IF( NTVI .GT. 0 ) THEN

         IF( NOTETR(1,NTVI) .NE. 0 ) THEN

C           NTVI TETRAEDRE VIDE QUI N'EST PAS VIDE!
C           ---------------------------------------
            print*
            print*,'probleme dans vetevide:  @@@@@@@@@@@@@@@@@@@@@@@@@@'
            print*,'vetevide: notetr(1,',n1tevi,') devrait etre =0!!!!'
            print*,'vetevide: n1tevi=',n1tevi,' notetr(',n1tevi,')=',
     %                   (notetr(kk,n1tevi),kk=1,8)
            print*,'probleme dans vetevide:  @@@@@@@@@@@@@@@@@@@@@@@@@@'
            ierr = 1
            goto 50
         ENDIF

C        CONTROLE DE BOUCLE INFINIE du CHAINAGE des TETRAEDRES VIDES
         NBTEVI = NBTEVI + 1
         IF( NBTEVI .LE. MXTETR ) THEN
C           TETRAEDRE VIDE SUIVANT
            NTVI = NOTETR( 5, NTVI )
            GOTO 10
         ELSE
            GOTO 30
         ENDIF

      ENDIF

C     PAS DE TETRAEDRE VIDE INCORRECT
C     -------------------------------
      ierr = 0
      GOTO 9999

 30   ierr = 2
      print*
      print*,'vetevide: BOUCLE INFINIE des TETRAEDRES VIDES'

C     CORRECTION de N1TEVI et du CHAINAGE des TETRAEDRES VIDES
C     --------------------------------------------------------
 50   print*,'vetevide: CORRECTION de N1TEVI et du CHAINAGE des TETRAEDR
     %ES VIDES'
      N1TEVI = 0
      DO NTVI=MXTETR,1,-1
         IF( NOTETR(1,NTVI) .EQ. 0 ) THEN
C           NTVI est UN TETRAEDRE VIDE
            NOTETR(5,NTVI) = N1TEVI
            N1TEVI = NTVI
         ENDIF
      ENDDO
      print*

 9999 RETURN
      END
