      SUBROUTINE AJTESTCF( KTITRE, PTXYZD, NOTETR, N1TETS,
     %                     NFLPER, MXFACO, LEFACO, NO0FAR,
     %                     NBTRCF, NOTRCF, NBSTCF, NOSTCF,
     %                     MXTECF, NBTECF, NOTECF, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    AJOUTER AUX TETRAEDRES DE L'ETOILE LES TETRAEDRES D'UN SOMMET
C -----    DU CF NON PRESENT PARMI LES SOMMETS DE L'ETOILE ET 2 SOMMETS
C          DE L'ETOILE

C ENTREES:
C --------
C KTITRE : TITRE D'UN TRACE
C PTXYZD : TABLEAU DES COORDONNEES DES POINTS
C          PAR POINT : X  Y  Z DISTANCE_SOUHAITEE
C NOTETR : LISTE DES TETRAEDRES
C          SOMMET1,    SOMMET2,    SOMMET3,    SOMMET4,
C          TETRAEDRE1, TETRAEDRE2, TETRAEDRE3, TETRAEDRE4
C          DE L'AUTRE COTE DE LA FACE
C N1TETS : N1TETS(NS) NUMERO D'UN TETRAEDRE AYANT POUR SOMMET NS

C NFLPER : NUMERO LEFACO DE LA FACE PERDUE A TRAITER
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
cccC          12: = NO FACEOC DE 1 A NBFACES D'OpenCascade
C NO0FAR : NUMERO DES 3 SOMMETS DES FACES AJOUTEES AU CF

C NBTRCF : NOMBRE DE FACES DE NOTRCF
C NOTRCF : >0 NUMERO DANS LEFACO DES TRIANGLES PERDUS  DU CF
C          <0 NUMERO DANS NO0FAR DES TRIANGLES AJOUTES AU CF
C NBSTCF : NOMBRE DE SOMMETS NOSTCF DU CF
C NOSTCF : NUMERO DANS PTXYZD DES NBSTCF SOMMETS DU CF

C SORTIES:
C --------
C NBTECF : NOMBRE DE TETRAEDRES DE L'ETOILE
C NOTECF : NUMERO NOTETR DES NBTECF TETRAEDRES DE L'ETOILE
C IERR   : =0  SI PAS D'ERREUR DETECTEE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Veulettes sur mer               Decembre 2019
C2345X7..............................................................012
      PARAMETER        (MXSTET=2048, MXTNSCF=2048)
      CHARACTER*(*)     KTITRE
      DOUBLE PRECISION  PTXYZD(4,*)
      INTEGER           NOTETR(8,*), N1TETS(*), NFETOI(3),
     %                  LEFACO(11,0:MXFACO), NO0FAR(3,*),
     %                  NOTRCF(NBTRCF), NOSTCF(NBSTCF), NOTECF(MXTECF),
     %                  NOSTET(MXSTET), NOTNSCF(MXTNSCF)

      IERR = 0
      NBTECF00 = NBTECF

C     -------------------------------------------------------------------
C     CONSTRUCTION Du TABLEAU DES NUMEROS PTXYZD DES SOMMETS DE L'ETOILE
C     -------------------------------------------------------------------
      NBTECF0 = NBTECF
      NBTECF  = 0
      NBSTET  = 0
      DO N = 1, NBTECF0
         NTE = NOTECF( N )
         IF( NTE .GT. 0 .AND. NOTETR(1,NTE) .GT. 0 ) THEN
            NBTECF = NBTECF + 1
            NOTECF( NBTECF ) = NTE
            DO 10 K=1,4
               NS = NOTETR( K, NTE )
C              SOMMET NS DEJA STOCKE?
               DO M=1,NBSTET
                  IF( NS .EQ. NOSTET(M) ) GOTO 10
               ENDDO
C              NON STOCKE: NS EST UN SOMMET DU CF DE PLUS
               IF( NBSTET .GE. MXSTET ) THEN
                  PRINT*,'ajtestcf: AUGMENTER MXSTET=',MXSTET
                  IERR = 1
                  GOTO 9999
               ENDIF
               NBSTET = NBSTET + 1
               NOSTET( NBSTET ) = NS
 10         ENDDO
         ENDIF
      ENDDO

      print*,'ajtestcf: Face PERDUE',NFLPER,' Tableau nostet des',NBSTET
     %      ,' sommets du CF :'
ccc      print*,'ajtestcf:',(nostet(n),n=1,nbstet)

C     RECHERCHE D'UN SOMMET DU CF NON SOMMET DE L'ETOILE
      DO 40 N = 1, NBSTCF

C        LE SOMMET N DU CF
         NSCF = NOSTCF( N )

         DO 20 M = 1, NBSTET
            NSET = NOSTET( M )
            IF( NSCF .EQ. NSET ) GOTO 40
 20      ENDDO

C        LE SOMMET NSCF N'EST PAS UN SOMMET DU CF
C        ----------------------------------------
C        RECHERCHE DE TOUS LES TETRAEDRES DE SOMMET NSCF
         CALL TETR1S( NSCF,    N1TETS,  NOTETR,
     %                NBTNSCF, MXTNSCF, NOTNSCF, IERR )

C        RECHERCHE DES TETRAEDRES DE SOMMET NSCF + 2 SOMMETS DE L'ETOILE
         NBTECF0 = NBTECF
         DO 35 K = 1, NBTNSCF
            NTE = NOTNSCF( K )

C           NOMBRE DE SOMMETS DE L'ETOILE DU TETRAEDRE NTE
            NBSET = 0
            DO 30 I = 1, 4
               NS = NOTETR( I, NTE )
               DO M = 1, NBSTET
                  IF( NS .EQ. NOSTET( M ) ) THEN
                     NBSET = NBSET + 1
                     GOTO 30
                  ENDIF
               ENDDO
 30         ENDDO

            IF( NBSET .GE. 2 ) THEN
C              TETRAEDRE 1St du CF + 2 St ETOILE:
C              IL EST AJOUTE A L'ETOILE SI CE N'EST DEJA FAIT
               DO  M = 1, NBTECF
                  IF( NOTECF(M) .EQ. NTE ) GOTO 35
               ENDDO
               IF( NBTECF .GE. MXTECF ) THEN
                 PRINT*,'ajtestcf: Tableau NOTECF SATURE MXTECF=',MXTECF
                  IERR = 2
                  GOTO 9999
               ENDIF
               NBTECF = NBTECF + 1
               NOTECF( NBTECF ) = NTE
               print*,'ajtestcf: Face PERDUE',NFLPER,' SOMMET',NSCF,
     %                ' NON SUR l''ETOILE: Ajout du TETRAEDRE',
     %                 NTE,':',(NOTETR(kk,NTE),kk=1,8)
            ENDIF

 35      ENDDO

         IF( NBTECF0 .LT. NBTECF ) THEN
      KTITRE='ajtestcf:        TETRAEDRES AJOUTES POUR ATTRAPER le SOMME
     %T          '
            WRITE(KTITRE(11:15),'(I5)') NBTECF - NBTECF0
            WRITE(KTITRE(61:67),'(I7)') NBSTET
            CALL SANSDBL( KTITRE, L )
            PRINT*, KTITRE(1:L)
            CALL TRFETO4( KTITRE(1:L), PTXYZD, 0, NFETOI,
     %                    NBTRCF, NOTRCF, LEFACO, NO0FAR,
     %                    NBTECF, NOTECF, NOTETR )
         ENDIF

 40   ENDDO

 9999 print*,'ajtestcf: Face PERDUE',NFLPER,': Ajout de',NBTECF-NBTECF00
     %      ,' TETRAEDRES'
      
      RETURN
      END
