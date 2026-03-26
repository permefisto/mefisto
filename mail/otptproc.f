      SUBROUTINE OTPTPROC( NUPRPTSF, NBSOMM, MXSOMM, PTXYZD, NPSOFR,
     %                     MXARBO, LARBRO, MXARBT, LARBRT, NUOTPT,
     %                     MXFACO, LEFACO,
     %                     MXPILE, MNPILE, NBSTSU )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    SUPPRIMER DE LA LISTE DES SOMMETS A TETRAEDRISER LES SOMMETS
C -----    D'OT TROP PROCHES D'UN POINT OU D'UNE FACE LEFACO
C          AJUSTAGE DE LA DISTANCE SOUHAITEE DES SOMMETS D'OT

C ENTREES:
C --------
C NUPRPTSF: NUMERO DU PREMIER POINT INTERNE OU FRONTALIER DANS PTXYZD
C NBSOMM : NOMBRE ACTUEL DE POINTS DECLARES DANS PTXYZD
C MXSOMM : NOMBRE MAXIMAL DE SOMMETS UTILISABLES DANS PTXYZD
C PTXYZD : PAR POINT : X  Y  Z  DISTANCE_SOUHAITEE
C MXARBO : NOMBRE MAXIMAL D'OCTAEDRES DANS LARBRO
C LARBRO : ARBRE-14 DES OCTAEDRES ( FOND DE LA TETRAEDRISATION )
C      LARBRO(0,0) : NO DU 1-ER OCTAEDRE VIDE DANS LARBRO
C      LARBRO(1,0) : MAXIMUM DU 1-ER INDICE DE LARBRO (ICI -1:20)
C      LARBRO(2,0) : MAXIMUM DECLARE DU 2-EME INDICE DE LARBRO
C                    (ICI = MXARBO)
C
C      LARBRO(-1:20,1) : RACINE DE L'ARBRE (OCTAEDRE SANS PERE)
C
C      LARBRO(-1,J) : NO DU PERE DE L'OCTAEDRE J DANS UN DES 2 ARBRES
C                     >0 => DANS LARBRO
C                     <0 => DANS LARBRT
C      LARBRO(0,J)  : 1 A 14 NO DE FILS DE L'OCTAEDRE J POUR SON PERE
C                     + 100 * NO TYPE DE L'OT J
C                     NO TYPE DE L'OT : 0 SI OCTAEDRE
C                                       1 SI TETRAEDRE T RONDE (T1)
C                                       2 SI TETRAEDRE T       (T2)
C   SI LARBRO(0,J)>0 ALORS J EST UN OCTAEDRE OCCUPE
C      SI LARBRO(1,.)>0 ALORS
C         LARBRO(1:14,J): NO (>0) LARBRO DES 14 SOUS-OCTA-TETRAEDRES
C      SINON
C         LARBRO(1:14,J):-NO PTXYZD DES 0 A 14 POINTS INTERNES DE L'OCTA J
C                         0  SI PAS DE POINT
C                        ( J EST ALORS UNE FEUILLE DE L'ARBRE )
C
C      LARBRO(15:20,J) : NO PTXYZD DES 6 SOMMETS DE L'OCTAEDRE J
C   SINON
C      LARBRO(0,J): -ADRESSE DANS LARBRO DE L'OCTAEDRE VIDE SUIVANT
C MXARBT : NOMBRE MAXIMAL DE TETRAEDRES DANS LARBRT
C LARBRT : ARBRE-5 DES TETRAEDRES ( FOND DE LA TETRAEDRISATION )
C      LARBRT(0,0) : NO DU 1-ER TETRAEDRE VIDE DANS LARBRT
C      LARBRT(1,0) : MAXIMUM DU 1-ER INDICE DE LARBRT (ICI -1:9)
C      LARBRT(2,0) : MAXIMUM DECLARE DU 2-EME INDICE DE LARBRT
C                     (ICI = MXARBT)
C
C      LARBRT(-1,J) : NO DU PERE DU TETRAEDRE J DANS UN DES 2 ARBRES
C                     >0 => DANS LARBRO
C                     <0 => DANS LARBRT
C      LARBRT(0,J) : 0 A 4 NO DE FILS DU TETRAEDRE J POUR SON PERE
C                    + 100 * NO TYPE DE L'OT J
C                     NO TYPE DE L'OT : 0 SI OCTAEDRE
C                                       1 SI TETRAEDRE T RONDE (T1)
C                                       2 SI TETRAEDRE T       (T2)
C
C   SI LARBRT(0,J)>0 ALORS J EST UN TETRAEDRE OCCUPE
C      SI LARBRT(1,J)>0 ALORS
C         LARBRT(1:5,J): NO (>0) LARBRT DES 5 SOUS-OCTA-TETRAEDRES
C      SINON
C         LARBRT(1:5,J):-NO PTXYZD DES 0 A 5 POINTS INTERNES AU TETRA J
C                         0  SI PAS DE POINT
C                        ( J EST ALORS UNE FEUILLE DU ARBRE )
C
C      LARBRT(6:9,J) : NO PTXYZD DES 4 SOMMETS DU TETRAEDRE J
C   SINON
C      LARBRT(0,J): ADRESSE DANS LARBRT DU TETRAEDRE VIDE SUIVANT
C
C NUOTPT : NUMERO D'OT (>0 SI LARBRO, <0 SI LARBRT) DE CHAQUE POINT PTXYZD
C          SI SOMMET D'OT, C'EST LE NUMERO D'OT DU PERE
C          SI POINT INTERNE OU FRONTALIER, C'EST LE NUMERO D'OT LE CONTENANT
C
C NBFACO : NOMBRE DE FACES TRIANGULAIRES RECENSEES SUR LE CONTOUR
C          ET LES INTERFACES ENTRE VOLUMES
C MXFACO : NOMBRE MAXIMAL DECLARABLE DANS LEFACO DE FACES DU CONTOUR
C LEFACO : FACE DU CONTOUR OU INTERFACES ENTRE VOLUMES
C          IL CONTIENT DANS CET ORDRE
C          NUMERO (DANS PTXYZD) DU SOMMET 1, SOMMET 2, SOMMET 3
C          NUMERO (DANS NUVOPA 0 SINON) DU VOLUME1 , VOLUME2 DE LA FACE
C          NUMERO (DANS LEFACO) DE LA FACE ADJACENTE PAR L'ARETE 1 2 3
C
C          ATTENTION: UNE ARETE PEUT APPARTENIR A PLUS DE 2 FACES
C          => CHAINAGE CIRCULAIRE DE CES FACES DANS LEFACO
C          LEFACO(9,*)  FACE SUIVANTE (*=0:VIDE, *<>0:NON VIDE)
C          HACHAGE AVEC LA SOMME DES 3 SOMMETS MODULO MXFACO
C          LF = MOD( NOSOFA(1)+NOSOFA(2)+NOSOFA(3) , MXFACO ) + 1
C          LE NUMERO DE LA 1-ERE FACE DANS LEFACO
C          NF = LEFACO( 10, LF )
C          SI LA FACE NE CONVIENT PAS. PASSAGE A LA SUIVANTE
C          NF  = LEFACO( 9, NF )  ...
C          LEFACO(10,*) PREMIERE FACE DANS LE HACHAGE
C N1FASC : N1FASC(NS)=NUMERO (DANS LEFACO) D'UNE FACE DE SOMMET NS
C
C ENTREES ET SORTIES:
C -------------------
C NPSOFR :  (1 000 000) * (NO SURFACE + 1) + NUMERO DU SOMMET DANS
C          LA NUMEROTATION DES SOMMETS DE LA SURFACE
C          -4 SI SOMMET DE LA GRILLE DES TETRAEDRES
C          -1 SI SOMMET DE LA GRILLE DES TETRAEDRES ELIMINE
C MXPILE : NOMBRE D'ENTIERS DECLARES POUR LE TABLEAU PILE
C MNPILE : ADRESSE MCN DU TABLEAU PILE
C          CE TABLEAU DOIT ETRE DECLARE EN ENTREE

C SORTIE :
C --------
C NBSTSU : NOMBRE DE SOMMETS SUPPRIMES MAIS PAS RETIRES de NBSOMM
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 1992
C2345X7..............................................................012
      include"./incl/langue.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)
      DOUBLE PRECISION  PTXYZD(4,MXSOMM),BAR(3), D2
      INTEGER           LEFACO(1:11,0:MXFACO),
     &                  LARBRO(-1:20,0:MXARBO),
     &                  LARBRT(-1:9,0:MXARBT),
     &                  NUOTPT(1:MXSOMM),
     &                  NPSOFR(1:MXSOMM)
      INTEGER           NUOTVO(256)

C     SUPPRESSION DES SOMMETS D'OT TROP PROCHES D'UN POINT INTERNE
C     OU D'UN SOMMET D'UNE FACE DE LEFACO ( SURFACE FERMEE TRIANGULEE )
C     LE SOMMET D'OT EST SUPPRIME SI SA DISTANCE AU POINT EST INFERIEURE
C     A LA MOITIE DE LA DISTANCE SOUHAITEE AUTOUR DU POINT
C     ==================================================================

C     LA BOUCLE SUR LES POINTS
      NBSTSU = 0
      DO 100 N = NUPRPTSF, NBSOMM
         IF( NUOTPT(N) .EQ. 0 ) GOTO 100
         IF( NPSOFR(N) .LE. 0 ) GOTO 100

C        N EST UN POINT INTERNE OU SOMMET D'UNE FACE DES SURFACES FERMEES
C        LA DISTANCE SOUHAITEE MAXIMALE A UN SOMMET D'OT
ccc         DISTST = REAL( PTXYZD(4,N) * 0.5D0 )
         DISTST = REAL( PTXYZD(4,N) * 0.333D0 )

C        RETROUVER LES OT VOISINS DU POINT N DONT UN SOMMET
C        EST A DISTANCE INFERIEURE A SQRT(DISTST)
         CALL OTVADP( N,      DISTST, PTXYZD, LARBRO, LARBRT, NUOTPT,
     %                NBOTVO, NUOTVO )

C        TOUT POINT DE CES OT VOISINS EST EVALUE EN DISTANCE A N
C        -------------------------------------------------------
         DO 80 K=1,NBOTVO

C           INITIALISATION DE LA PILE
            LHPILE = 1
            MCN( MNPILE ) = NUOTVO( K )
C
C           TANT QUE LA PILE N'EST PAS VIDE
 10         IF( LHPILE .GT. 0 ) THEN
C
C              L'OT A TRAITER EST DEPILE
               LHPILE = LHPILE - 1
               NUOT   = MCN( MNPILE + LHPILE )
               IF( NUOT .GT. 0 ) THEN
C
C                 OCTAEDRE
                  IF( LARBRO(1,NUOT) .GT. 0 ) THEN
C                    NUOT SANS POINT : LES 14 FILS SONT EMPILES
                     IF( LHPILE + 14 .GT. MXPILE ) THEN
C                       SATURATION DE LA PILE => ELLE EST AGRANDIE
                        CALL TNMCAU( 'ENTIER', MXPILE, MXPILE+2048,
     %                                MXPILE, MNPILE )
                        MXPILE = MXPILE + 2048
                     ENDIF
                     DO 20 I=1,14
                        MCN( MNPILE-1+LHPILE+I ) = LARBRO(I,NUOT)
 20                  CONTINUE
                     LHPILE = LHPILE + 14
                  ENDIF
                  DO 30 I=15,20
C                    LE NUMERO PTXYZD DU SOMMET I-14 DE L'OT
                     NS = LARBRO(I,NUOT)
                     IF( NPSOFR(NS) .EQ. -4 ) THEN
C                       CALCUL DE LA DISTANCE DE NS AU POINT N
                        D2 = ( PTXYZD(1,N) - PTXYZD(1,NS) ) ** 2
     %                     + ( PTXYZD(2,N) - PTXYZD(2,NS) ) ** 2
     %                     + ( PTXYZD(3,N) - PTXYZD(3,NS) ) ** 2
                        IF( D2 .LT. DISTST**2 ) THEN
C                          LE SOMMET NS EST TROP PRES DU POINT N
C                          IL EST SUPPRIME DES POINTS A TETRAEDRISER
ccc                           PRINT*,'otptproc: SUPPRESSION du SOMMET',NS,
ccc     %                            ' TROP PROCHE DU SOMMET',N,
ccc     %                            ' a la DISTANCE=',SQRT(D2)
ccc                           PRINT*,'otptproc: PTXYZD(',NS,')=',
ccc     %                             (PTXYZD(kk,NS),kk=1,4)
ccc                           PRINT*,'otptproc: PTXYZD(',N,')=',
ccc     %                             (PTXYZD(kk,N),kk=1,4)
                           NPSOFR( NS ) = -1
                           NBSTSU = NBSTSU + 1
                        ENDIF
                     ENDIF
 30               CONTINUE

               ELSE

C                 TETRAEDRE
                  NUOTT = -NUOT
                  IF( LARBRT(1,NUOTT) .GT. 0 ) THEN
C                    NUOTT SANS POINT : LES 5 FILS SONT EMPILES
                     IF( LHPILE + 5 .GT. MXPILE ) THEN
C                       SATURATION DE LA PILE => ELLE EST AGRANDIE
                        CALL TNMCAU( 'ENTIER', MXPILE, MXPILE+2048,
     %                                MXPILE, MNPILE )
                        MXPILE = MXPILE + 2048
                     ENDIF
                     DO 50 I=1,5
                        MCN( MNPILE-1+LHPILE+I ) = LARBRT(I,NUOTT)
 50                  CONTINUE
                     LHPILE = LHPILE + 5
                  ENDIF
                  DO 60 I=6,9
C                    LE NUMERO PTXYZD DU SOMMET I-5 DE L'OT
                     NS = LARBRT(I,NUOTT)
                     IF( NPSOFR(NS) .EQ. -4 ) THEN
C                       CALCUL DE LA DISTANCE DE NS AU POINT N
                        D2 = ( PTXYZD(1,N) - PTXYZD(1,NS) ) ** 2
     %                     + ( PTXYZD(2,N) - PTXYZD(2,NS) ) ** 2
     %                     + ( PTXYZD(3,N) - PTXYZD(3,NS) ) ** 2
                        IF( D2 .LT. DISTST**2 ) THEN
C                          LE SOMMET NS EST TROP PRES DU POINT N
C                          IL EST SUPPRIME DES POINTS A TETRAEDRISER
                           NPSOFR( NS ) = -1
                           NBSTSU = NBSTSU + 1
                        ENDIF
                     ENDIF
 60               CONTINUE
               ENDIF
               GOTO 10
            ENDIF
 80      CONTINUE
 100   CONTINUE

C     SUPPRESSION DES SOMMETS D'OT TROP PROCHES D'UNE FACE DE LEFACO
C     LE SOMMET D'OT EST SUPPRIME SI SA DISTANCE AU BARYCENTRE DE LA FACE
C     EST INFERIEURE A LA MOITIE DE LA PLUS LONGUE ARETE DU TRIANGLE
C     ===================================================================
      DO 200 N=1,MXFACO
         IF( LEFACO(1,N) .EQ. 0 ) GOTO 200

C        FACE OCCUPEE : CALCUL DU BARYCENTRE ET MAX DES DISTANCES
C        SOUHAITEES AUTOUR DES 3 SOMMETS DE LA FACE
         DISTST = 0
         BAR(1) = 0
         BAR(2) = 0
         BAR(3) = 0

         DO 120 I=1,3
            NS     = LEFACO(I,N)
            DO 115 J=1,3
C              CALCUL DU BARYCENTRE
               BAR(J) = BAR(J) + PTXYZD(J,NS)
 115        CONTINUE
C           LA DISTANCE MAXIMALE SOUHAITEE POUR LES 3 SOMMETS
            DISTST = REAL( MAX( DBLE(DISTST), PTXYZD(4,NS) ) )
 120     CONTINUE

C        LES 3 COORDONNEES DU BARYCENTRE DE LA FACE N DE LEFACO
         BAR(1) = BAR(1) / 3
         BAR(2) = BAR(2) / 3
         BAR(3) = BAR(3) / 3

C        LE CARRE DE LA LONGUEUR/2 DE LA PLUS LONGUE ARETE DE LA FACE N
         DISTST = DISTST * 0.5

C        RETROUVER LES OT VOISINS D'UN POINT PT DONT UN SOMMET
C        EST A DISTANCE INFERIEURE A DISTST
         NBOTVO = 0
         DO 125 I=1, 3
            NS = LEFACO(I,N)
            CALL OTVADP( NS,    DISTST, PTXYZD, LARBRO, LARBRT, NUOTPT,
     %                   NBOTV, NUOTVO(NBOTVO+1) )
C           AJOUT SANS DUPLICATION DES NBOTV OT VOISINS
            CALL AJPESR( NBOTV, NUOTVO(NBOTVO+1), NBOTVO, NUOTVO )
 125     CONTINUE

C        BOUCLE SUR LES SOMMETS DES OT VOISINS DE LA FACE N DE LEFACO
         DO 180 K=1,NBOTVO

C           INITIALISATION DE LA PILE
            LHPILE = 1
            MCN(MNPILE) = NUOTVO( K )

C           TANT QUE LA PILE N'EST PAS VIDE
 130        IF( LHPILE .GT. 0 ) THEN

C              L'OT A TRAITER EST DEPILE
               LHPILE = LHPILE - 1
               NUOT   = MCN( MNPILE+LHPILE )
               IF( NUOT .GT. 0 ) THEN
c
C                 OCTAEDRE
                  IF( LARBRO(1,NUOT) .GT. 0 ) THEN
C                    NUOT SANS POINT : LES 14 FILS SONT EMPILES
                     IF( LHPILE + 14 .GT. MXPILE ) THEN
C                       SATURATION DE LA PILE => ELLE EST AGRANDIE
                        CALL TNMCAU( 'ENTIER', MXPILE, MXPILE+2048,
     %                                MXPILE, MNPILE )
                        MXPILE = MXPILE + 2048
                     ENDIF
                     DO 135 I=1,14
                        MCN( MNPILE-1+LHPILE+I ) = LARBRO(I,NUOT)
 135                 CONTINUE
                     LHPILE = LHPILE + 14
                  ENDIF
                  DO 140 I=15,20
C                    LE NUMERO PTXYZD DU SOMMET I-14 DE L'OT
                     NS = LARBRO(I,NUOT)
                     IF( NPSOFR(NS) .EQ. -4 ) THEN
C                       CALCUL DE LA DISTANCE DE NS AU POINT N
                        D2 = ( BAR(1) - PTXYZD(1,NS) ) ** 2
     %                     + ( BAR(2) - PTXYZD(2,NS) ) ** 2
     %                     + ( BAR(3) - PTXYZD(3,NS) ) ** 2
                        IF( D2 .LT. DISTST**2 ) THEN
C                          LE SOMMET NS EST TROP PRES DU POINT N
C                          IL EST SUPPRIME DES POINTS A TETRAEDRISER
                           NPSOFR( NS ) = -1
                           NBSTSU = NBSTSU + 1
                        ENDIF
                     ENDIF
 140              CONTINUE

               ELSE

C                 TETRAEDRE
                  NUOTT = -NUOT
                  IF( LARBRT(1,NUOTT) .GT. 0 ) THEN
C                    NUOTT SANS POINT : LES 5 FILS SONT EMPILES
                     IF( LHPILE + 5 .GT. MXPILE ) THEN
C                       SATURATION DE LA PILE => ELLE EST AGRANDIE
                        CALL TNMCAU( 'ENTIER', MXPILE, MXPILE+2048,
     %                                MXPILE, MNPILE )
                        MXPILE = MXPILE + 2048
                     ENDIF
                     DO 150 I=1,5
                        MCN( MNPILE-1+LHPILE+I ) = LARBRT(I,NUOTT)
 150                 CONTINUE
                     LHPILE = LHPILE + 5
                  ENDIF
                  DO 160 I=6,9
C                    LE NUMERO PTXYZD DU SOMMET I-5 DE L'OT
                     NS = LARBRT(I,NUOTT)
                     IF( NPSOFR(NS) .EQ. -4 ) THEN
C                       CALCUL DE LA DISTANCE DE NS AU POINT N
                        D2 = ( BAR(1) - PTXYZD(1,NS) ) ** 2
     %                     + ( BAR(2) - PTXYZD(2,NS) ) ** 2
     %                     + ( BAR(3) - PTXYZD(3,NS) ) ** 2
                        IF( D2 .LT. DISTST**2 ) THEN
C                          LE SOMMET NS EST TROP PRES DU POINT N
C                          IL EST SUPPRIME DES POINTS A TETRAEDRISER
                           NPSOFR( NS ) = -1
                           NBSTSU = NBSTSU + 1
                        ENDIF
                     ENDIF
 160              CONTINUE
               ENDIF
               GOTO 130
C
            ENDIF
 180     CONTINUE
 200  CONTINUE

C     AJUSTAGE DE LA DISTANCE SOUHAITEE DES SOMMETS D'OT
      DO K=1,NBSOMM
         IF( NPSOFR(K) .EQ. -4 .OR. NPSOFR(K) .EQ. -1 ) THEN
C           RECHERCHE DE L'OT DE CE SOMMET
            NUOT = NUOTPT( K )
C           CARRE DE LA LONGUEUR DE SON ARETE
            CALL LOAROT( NUOT, LARBRO, LARBRT, PTXYZD, TAILAR )
C           LA DISTANCE SOUHAITEE AUTOUR DU SOMMET
            PTXYZD( 4, K ) = SQRT( TAILAR )
         ENDIF
      ENDDO

      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,*)'otptproc: NOMBRE SOMMETS D''OT TROP PROCHES SUP
     %PRIMES =',NBSTSU
      ELSE
         WRITE(IMPRIM,*)'otptproc: SUPPRESSED TOO NEAR OT VERTICES Numbe
     %r=',NBSTSU
      ENDIF

      RETURN
      END
