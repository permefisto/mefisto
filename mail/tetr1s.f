      SUBROUTINE TETR1S( NOST,   N1TETS, NOTETR,
     %                   NBTE1S, MXTE1S, NOTE1S, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    LISTER DANS NOTE1S LE NUMERO NOTETR DES NBTE1S TETRAEDRES
C -----    DE SOMMET NOST

C ENTREES:
C --------
C NOST   : NUMERO DU SOMMET COMMUN A TOUS LES TETRAEDRES
C N1TETS : N1TETS(I) NUMERO D'UN TETRAEDRE AYANT POUR SOMMET I
C NOTETR : LISTE DES TETRAEDRES
C          SOMMET1,    SOMMET2,    SOMMET3,    SOMMET4,
C          TETRAEDRE1, TETRAEDRE2, TETRAEDRE3, TETRAEDRE4
C          DE L'AUTRE COTE DE LA FACE DE NUMEROTATION
C          1: 123      2: 234      3: 341      4: 412
C MXTE1S : NOMBRE D'ENTIERS DU TABLEAU NOTE1S

C SORTIES:
C --------
C NBTE1S : NOMBRE DE TETRAEDRES DE SOMMET NOST SI IERR=0
C NOTE1S : NOTE1S(I) >0 NUMERO NOTETR DU TETRAEDRE I DE SOMMET NOST
C IERR   : = 0 PAS D'ERREUR DETECTEE
C          = 1 SI LE SOMMET NOST N'EST PAS UN SOMMET ACTIF
C              CAS D'UN TETRAEDRE N1TETS(NOST)=<0 DESACTIVE
C          = 2 SI LE TETRAEDRE N1TETS(NOST) N'EST PAS ACTIF
C          = 3 SI SATURATION DU TABLEAU NOTE1S
C          = 4 SI UN TETRAEDRE DE LA PILE EST VIDE
C          = 5 SI UN TETRAEDRE OPPOSE EST VIDE
C          = 6 SI LE SOMMET NOST N'EST PAS DANS LE TETRAEDRE N1TETS(NOST)
C                 CELA PEUT ETRE LE CAS d''un SOMMET DESACTIVE'
C          = 7 SI UN TETRAEDRE OPPOSE A UNE FACE N'A PAS UN SOMMET DE
C                 LA FACE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS UPMC    Decembre 2005
C MODIFS : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY  Octobre 2014
C MODIFS : ALAIN PERRONNET             St PIERRE du PERRAY  Janvier 2019
C2345X7..............................................................012
      PARAMETER        (MXPILE=16384)
      INTEGER           NOTETR(8,*),
     %                  N1TETS(*),
     %                  NOTE1S( MXTE1S ) 
      INTEGER           LAPILE( MXPILE )
      INTEGER           NOFA1STE( 3, 4 )
      DATA              NOFA1STE / 1,3,4, 2,4,1, 3,1,2, 4,2,3 /  

      IERR   = 0
      NBTE1S = 0

C     NTE0 UN PREMIER TETRAEDRE DE SOMMET NOST
C     ----------------------------------------
 10   NTE0 = N1TETS( NOST )

      IF( NTE0 .LE. 0 ) THEN
         IERR = 1
         PRINT*,'tetr1s: Le Point',NOST,' APPARTIENT a AUCUN TETRAEDRE. 
     % N1TETS(',NOST,')=',NTE0,' IERR=',IERR
         GOTO 9999
      ENDIF

      IF( NOTETR(1,NTE0) .LE. 0 ) THEN

C        ERREUR: NTE0 EST UN TETRAEDRE VIDE
         PRINT*,'ERREUR tetr1s: St=',NOST,' NOTETR(',NTE0,')=',
     %          (NOTETR(I,NTE0),I=1,8)

C        RECHERCHE D'UN TETRAEDRE DE SOMMET NOST DANS LA
C        TETRAEDRISATION ACTUELLE
C        -----------------------------------------------
         MAXTETR = 0
         DO I=5,8
            NTE = NOTETR(I,NTE0)
            IF( NTE .GT. MAXTETR ) THEN
               MAXTETR = NTE
            ENDIF
         ENDDO

C        MISE A JOUR DE N1TETS ET MAXTETR JUSQU'A TROUVER NOST
         DO NTE=1,MAXTETR
            IF( NOTETR(1,NTE) .GT. 0 ) THEN
               DO N=1,4

                  NTEOP = NOTETR(4+N,NTE)
                  IF( NTEOP .GT. MAXTETR ) THEN
                     MAXTETR = NTEOP
                  ENDIF

                  NS = NOTETR(N,NTE)
                  N1TETS( NS ) = NTE

                  IF( NS .EQ. NOST ) THEN
                     GOTO 10
                  ENDIF

               ENDDO
            ENDIF
         ENDDO

C        NOST N'EST PAS UN SOMMET D'UN TETRAEDRE ACTIF DE LA TETRAEDRISATION
         IERR = 2
         PRINT*
         PRINT*,'Attention tetr1s: No St',NOST,' N''EST PAS UN SOMMET D'
     %'UN TETRAEDRE ACTIF IERR=',IERR
         PRINT*
         GOTO 9999

      ENDIF

C     RECHERCHE DU NUMERO LOCAL I DU SOMMET NOST DANS LE TETRAEDRE NTE0
      DO I=1,4
         IF( NOST .EQ. NOTETR(I,NTE0) ) GOTO 20
      ENDDO

      IERR = 6
      PRINT*
      PRINT*,'tetr1s: No St',NOST,
     %       ' NON SOMMET du TETRAEDRE INITIAL',NTE0,' de SOMMETS',
     %        (NOTETR(J,NTE0),J=1,8)
      PRINT*,'tetr1s: REVOIR la VALEUR de N1TETS(',NOST,')=',NTE0,
     %       ' et/ou NPSOFR(',NOST,'). IERR=',IERR
      PRINT*
      GOTO 9999

C     NOST EST LE SOMMET I DU TETRAEDRE NTE0
C     NTE0 EST LE PREMIER TETRAEDRE EMPILE DE SOMMET NOST
 20   LHPILE = 1
      LAPILE( LHPILE ) = NTE0


C     TRAITEMENT DU TETRAEDRE DE SOMMET NOST EN HAUT DE PILE
C     ------------------------------------------------------
 30   IF( LHPILE .GT. 0 ) THEN

C        LE TETRAEDRE NT CONTIENT NOST
         NT = LAPILE( LHPILE )

C        LE HAUT DE PILE EST TRAITE
         LHPILE = LHPILE - 1

         IF( NOTETR(1,NT) .LE. 0 ) THEN
            IERR = 4
            PRINT*,'ERREUR tetr1s: St',NOST,' NOTETR(',NT,')=',
     %      (NOTETR(J,NT),J=1,8),' de PREMIER SOMMET',NOTETR(1,NT),
     %       ' ? IERR=',IERR

            GOTO 30
         ENDIF

C        RECHERCHE DU NUMERO LOCAL I DU SOMMET NOST DANS LE TETRAEDRE NT
         DO I=1,4
            IF( NOST .EQ. NOTETR(I,NT) ) GOTO 40
         ENDDO

         PRINT*
         PRINT*,'tetr1s: No St',NOST,
     %          ' NON SOMMET du TETRAEDRE EMPILE',NT,' de SOMMETS',
     %          (NOTETR(kk,NT),kk=1,8)
         PRINT*,'tetr1s: LES TETRAEDRES DE LA PILE DE HAUT EN BAS'
         DO M = LHPILE, 1, -1
            NTE = LAPILE( M )
            PRINT*,'tetr1s: NOTETR(',NTE,')=',(NOTETR(kk,NTE),kk=1,8)
         ENDDO
         PRINT*
         GOTO 30

C        NT EST UN TETRAEDRE DE SOMMET NOST DONC IL EST A RANGER
C        DANS NOTE1S SI CE N'EST DEJA FAIT
 40      DO N = 1, NBTE1S
            IF( NT .EQ. NOTE1S(N) ) GOTO 30
         ENDDO
         IF( NBTE1S .GE. MXTE1S ) THEN
            GOTO 9000
         ENDIF
         NBTE1S = NBTE1S + 1
         NOTE1S( NBTE1S ) = NT

C        TRAITEMENT DES 3 TETRAEDRES ADJACENTS PAR LES 3 FACES DE NT
C        DE SOMMET NOST, I DANS NT
         DO 50 M=1,3

C           NO DE LA FACE M CONTENANT LE SOMMET I DE NT
            NF = NOFA1STE( M, I )

C           LE TETRAEDRE NTOP OPPOSE A LA FACE NF DE NT EXISTE T IL?
            NTOP = NOTETR( 4+NF, NT )

            IF( NTOP .LE. 0 ) GOTO 50
            IF( NOTETR(1,NTOP) .LE. 0 ) THEN
               IERR = 5
               PRINT*
               PRINT*,'tetr1s: Probleme du TETRAEDRE  ',NT,':',
     %                (NOTETR(kk,NT),kk=1,8)
               PRINT*,'tetr1s: a pour TETRAEDRE OPPOSE',NTOP,':',
     %               (NOTETR(kk,NTOP),kk=1,8),' qui est VIDE IERR=',IERR
               PRINT*

               GOTO 50
            ENDIF

C           VERIFICATION NOST SOMMET DE NTOP?
            DO N = 1, 4
               IF( NOST .EQ. NOTETR(N,NTOP) ) GOTO 45
            ENDDO

            IERR = 7
            PRINT*
            PRINT*,'tetr1s: Recherche des TETRAEDRES de SOMMET',NOST
            PRINT*,'tetr1s: Probleme du TETRAEDRE      ',NT,':',
     %             (NOTETR(kk,NT),kk=1,8),' FACE',NF
            PRINT*,'tetr1s: qui a pour TETRAEDRE OPPOSE',NTOP,':',
     %             (NOTETR(kk,NTOP),kk=1,8),' SANS SOMMET',NOST,
     %             ' IERR=',IERR
            PRINT*
            GOTO 9999

C           NTOP EST EMPILE SI CE N'EST DEJA FAIT
 45         DO N = LHPILE, 1, -1
               IF( LAPILE(N) .EQ. NTOP ) GOTO 50
            ENDDO

C           NTOP EST EMPILE
            IF( LHPILE .GE. MXPILE ) THEN
               GOTO 8000
            ENDIF
            LHPILE = LHPILE + 1
            LAPILE( LHPILE ) = NTOP

 50      ENDDO

C        PASSAGE AU TETRAEDRE DE SOMMET NOST SUIVANT A TRAITER
         GOTO 30

      ENDIF

C     NBTE1S NOMBRE TOTAL DE TETRAEDRES de SOMMET NOST STOCKES DANS NOTE1S
      GOTO 9999


C     SATURATION DE LA PILE
 8000 PRINT*,'tetr1s: St ',NOST,' PILE SATUREE. AUGMENTER MXPILE=',
     %MXPILE,' => RETOUR AVEC NBTE1S=',NBTE1S,' IERR=',IERR
      GOTO 9999


C     AJOUT IMPOSSIBLE CAR LE TABLEAU NOTE1S EST SATURE
 9000 IERR = 3
      PRINT*,'tetr1s: SATURATION tableau NOTE1S pour le St ',NOST,
     %       ' AUGMENTER MXTE1S=',MXTE1S,' => RETOUR AVEC NBTE1S=',
     %        NBTE1S,' IERR=',IERR

ccc      PRINT*,'tetr1s: les 10 DERNIERS TETRAEDRES de NOTE1S'
ccc      DO M = NBTE1S, MAX(1,NBTE1S-10), -1
ccc         NT = NOTE1S( M )
ccc         PRINT*,'tetr1s: M=',M,' NOTETR(',NT,')=',(NOTETR(kk,NT),kk=1,8)
ccc      ENDDO


 9999 RETURN
      END
