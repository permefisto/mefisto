      SUBROUTINE OTARDI( NOFOTI, ARETGR, ARETOCG, HEXA,
     %                   MXARBO, LARBRO, MXARBT, LARBRT, NUOTPT,
     %                   NBSOMM, MXSOMM, PTXYZD,
ccc     %                NPSOFR,
     %                   LAVIDE, MXARET, LARETE, MXPILE, LAPILE,
     %                   IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CREATION DES OCTAEDRES ET TETRAEDRES POUR QUE TOUT POINT
C -----    APPARTIENNE A UN OT DONT LES ARETES SOIENT INFERIEURES
C          A SA DISTANCE SOUHAITEE AUTOUR AVEC PRISE EN COMPTE
C          DE LA FONCTION TAILLE_IDEALE(x,y,z) SI ELLE EXISTE

C ENTREES:
C --------
C NOFOTI : NUMERO DE LA FONCTION TAILLE_IDEALE(x,y,z) DES ARETES
C ARETGR : LONGUEUR SOUHAITEE DE L'ARETE DE LA GRILLE REGULIERE
C ARETOCG: DISTANCE SOUHAITEE des SOMMETS DE L'OCTAEDRE ENGLOBANT
C HEXA   : MIN ET MAX DES 3 COORDONNEES D'UN HEXAEDRE ENGLOBANT LES POINTS
C MXARBO : NOMBRE MAXIMAL D'OCTAEDRES DANS LARBRO
C MXARBT : NOMBRE MAXIMAL DE TETRAEDRES DANS LARBRT
C MXSOMM : NOMBRE MAXIMAL DE SOMMETS PERMIS POUR LA TETRAEDRISATION
C MXARET : NOMBRE MAXIMAL D'ARETES DES OT

C AUXILIAIRES:
C ------------
C MXPILE : MAXIMUM DE OT DANS LA PILE
C LAPILE : LES OT DE LA PILE

C ENTREES ET SORTIES:
C -------------------
C NBSOMM : NOMBRE ACTUEL DE POINTS DECLARES DANS PTXYZD
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

C NUOTPT : NUMERO D'OT (>0 SI LARBRO, <0 SI LARBRT) DE CHAQUE POINT PTXYZD
C PTXYZD : X Y Z DISTANCE SOUHAITEE DES POINTS

cccC NPSOFR : (I) =  0 SI LE POINT EST AJOUTE NON SUR LA SURFACE DE L'OBJET
cccC              =  1 SI LE POINT EST AJOUTE SUR UNE ARETE DE LA SURFACE
cccC              =  2 SI LE POINT EST AJOUTE COMME BARYCENTRE SUR LA SURFACE
cccC              =  1 000 000  + NUMERO DU  POINT INTERNE UTILISATEUR
cccC              = (1 000 000) * (NO SURFACE + 1) + NUMERO DU SOMMET DANS
cccC                   LA NUMEROTATION DES SOMMETS DE LA SURFACE
cccC              = -4 SI LE POINT D'OT N'EST PAS SUR LA FRONTIERE
cccC              = -NPSOFR(I) SI POINT DEPLACE SUR LA SURFACE

C LAVIDE : PREMIERE ADRESSE LARETE A EXPLORER POUR TROUVER UNE ARETE VIDE

C SORTIES:
C --------
C LARETE : TABLEAU DE HACHAGE DES ARETES DES OT
C          LARETE(1,J) = SOMMET 1 DANS PTXYZD
C          LARETE(2,J) = SOMMET 2 DANS PTXYZD
C          LARETE(3,J) = NUMERO DANS LARETE DE L'ARETE SUIVANTE
C          LARETE(4,J) = MILIEU DE L'ARETE DANS PTXYZD
C          HACHAGE(NS1,NS2) = (NS1 + NS2) MODULO MXARET + 1
C IERR   : 0  SI PAS D'ERREUR,  >0 SI ERREUR
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC       JUIN 1997
C MODIF  : ALAIN PERRONNET Saint Pierre du Perray         Septembre 2019
C2345X7..............................................................012
      include"./incl/gsmenu.inc"
      INTEGER           LARBRO(-1:20,0:MXARBO),
     %                  LARBRT(-1:9 ,0:MXARBT),
     %                  NUOTPT(1:MXSOMM),
ccc     %                  NPSOFR(1:MXSOMM),
     %                  LARETE(1:4 ,1:MXARET)
      DOUBLE PRECISION  PTXYZD(4,MXSOMM)
      INTEGER           LAPILE(1:MXPILE)
      REAL              HEXA(6,2), ARETOCG, ARETGR

      IF( NOFOTI .GT. 0 ) THEN

C        LA DISTANCE SOUHAITEE EST CELLE FOURNIE PAR TAILLE_IDEALE(x,y,z)
         DO I=7,NBSOMM
            CALL OTTAID( NOFOTI, ARETOCG, HEXA, PTXYZD(1,I),
     %                   PTXYZD(4,I), IERR )
            IF( IERR .NE. 0 ) RETURN
         ENDDO
         ARETG2 = 0

      ELSE

C        DESCENTE DES ARBRES JUSQU'A ARETE < ARETGR * SQRT(2)
         ARETG2 = ( ARETGR ** 2 ) * 2

      ENDIF

C     SI L'OCTAEDRE ENGLOBANT N'EST PAS ENCORE SUBDIVISE
C     ALORS SA SUBDIVISION EST IMPOSEE
C     ==================================================
      NUMOT = 1
      IF( LARBRO(1,1) .LE. 0 ) THEN

C        DECOMPOSITION DE L'OCTAEDRE ENGLOBANT
         NBSOM0 = NBSOMM
         CALL SBOCTA( 0,      HEXA,   NUMOT,
     %                MXARBO, LARBRO, MXARBT, LARBRT, NUOTPT,
     %                NBSOMM, MXSOMM, PTXYZD,
     %                LAVIDE, MXARET, LARETE, IERR )
         IF( IERR .NE. 0 ) RETURN

         IF( NOFOTI .GT. 0 ) THEN
C           MISE A JOUR DE TAILLE_IDEALE DES NOUVEAUX SOMMETS
            DO I=NBSOM0+1,NBSOMM
               CALL OTTAID( NOFOTI, ARETOCG, HEXA, PTXYZD(1,I),
     %                      PTXYZD(4,I), IERR )
               IF( IERR .NE. 0 ) RETURN
            ENDDO
         ENDIF

      ENDIF


C     PARCOURS DES ARBRES LARBRO et LARBRT
C     ====================================
C     INITIALISATION DE LA PILE AVEC L'OCTAEDRE ENGLOBANT
      LHPILE = 1
      LAPILE(1) = 1

C     TANT QUE L'ARETE > DISTANCE_SOUHAITEE d'1 SOMMET : DECOMPOSER L'OT
C     .................................................  ...............
 10   IF( LHPILE .GT. 0 ) THEN

C        L'OT DU HAUT DE LA PILE EST DEPILE
         NUMOT  = LAPILE( LHPILE )
         LHPILE = LHPILE - 1
C
         IF( NUMOT .GT. 0 ) THEN
C
C           OCTAEDRE
            IF( LARBRO(1,NUMOT) .GT. 0 ) THEN
C              OCTAEDRE NON FEUILLE. SES FILS SONT EMPILES
               IF( LHPILE+14 .GT. MXPILE ) THEN
C                 SATURATION DE LA PILE
                  NBLGRC(NRERR) = 1
                  KERR(1) ='otardi: LAPILE SATUREE MXPILE=             '
                  WRITE(KERR(1)(31:37),'(I7)') MXPILE
                  CALL LEREUR
                  IERR = 8
                  RETURN
               ENDIF
C              LES 14 OT SONT EMPILES
               DO I=14,1,-1
C                 LEUR NUMERO DANS LARBR
                  LAPILE(LHPILE+I) = LARBRO(I,NUMOT)
               ENDDO
               LHPILE = LHPILE + 14
               GOTO 10
            ENDIF

         ELSE

C           TETRAEDRE
            NUMOTT = - NUMOT
            IF( LARBRT(1,NUMOTT) .GT. 0 ) THEN
C              TETRADRE NON FEUILLE. SES FILS SONT EMPILES
               IF( LHPILE+5 .GT. MXPILE ) THEN
C                 SATURATION DE LA PILE
                  NBLGRC(NRERR) = 1
                  KERR(1) ='otardi: LAPILE SATUREE MXPILE=             '
                  WRITE(KERR(1)(31:37),'(I7)') MXPILE
                  CALL LEREUR
                  IERR = 8
                  RETURN
               ENDIF
C              LES 5 OT SONT EMPILES
               DO I=5,1,-1
C                 LEUR NUMERO DANS LARBRT OU LARBRO
                  LAPILE(LHPILE+I) = LARBRT(I,NUMOTT)
               ENDDO
               LHPILE = LHPILE + 5
               GOTO 10
            ENDIF

         ENDIF

C        ICI, L'OT NUMOT EST UN OT FEUILLE
C        .................................
C        ARETE2 = CARRE DE LA LONGUEUR DE L'ARETE DE L'OT NUMOT
         CALL LOAROT( NUMOT, LARBRO, LARBRT, PTXYZD, ARETE2 )
C
         IF( NOFOTI .GT. 0 ) THEN
cccC           LE MINIMUM DES DISTANCES SOUHAITEES AUX SOMMETS DE L'OT
C           LE MAXIMUM DES DISTANCES SOUHAITEES AUX SOMMETS DE L'OT
            IF( NUMOT .GT. 0 ) THEN
C              OCTAEDRE => 6 SOMMETS
ccc               ARETG2 = REAL( MIN( PTXYZD(4,LARBRO(15,NUMOT)),
               ARETG2 = REAL( MAX( PTXYZD(4,LARBRO(15,NUMOT)),
     %                             PTXYZD(4,LARBRO(16,NUMOT)),
     %                             PTXYZD(4,LARBRO(17,NUMOT)),
     %                             PTXYZD(4,LARBRO(18,NUMOT)),
     %                             PTXYZD(4,LARBRO(19,NUMOT)),
     %                             PTXYZD(4,LARBRO(20,NUMOT))  ) )
            ELSE
C              TETRAEDRE => 4 SOMMETS
               NUMOTT = - NUMOT
ccc               ARETG2 = REAL( MIN( PTXYZD(4,LARBRT(6,NUMOTT)),
               ARETG2 = REAL( MAX( PTXYZD(4,LARBRT(6,NUMOTT)),
     %                             PTXYZD(4,LARBRT(7,NUMOTT)),
     %                             PTXYZD(4,LARBRT(8,NUMOTT)),
     %                             PTXYZD(4,LARBRT(9,NUMOTT)) ) )
            ENDIF

cccC           LE CARRE POUR LA COMPARAISON ET LA RACINE CARREE(2)
ccc            ARETG2 = ( ARETG2 ** 2 ) * 2

C           LE CARRE POUR LA COMPARAISON SQRT(1.21)=1.1
            ARETG2 = ( ARETG2 ** 2 ) * 1.21
         ENDIF
C
C        SI MIN DES DISTANCES SOUHAITEES DES SOMMETS DE L'OT
C        EST SUPERIEURE A LA TAILLE DE L'ARETE DE L'OT
C        ALORS L'OT N'EST PAS DECOMPOSE
         IF( ARETG2 .GE. ARETE2 ) GOTO 10
C
         IF( NUMOT .GT. 0 ) THEN
C
C           OCTAEDRE : A DECOMPOSER SI INTERNE A L'HEXAEDRE
C           ===============================================
C           L'OCTAEDRE EST IL EXTERIEUR A L'HEXAEDRE ENGLOBANT DES POINTS?
            CALL EXTHEX( 6, LARBRO(15,NUMOT), PTXYZD, HEXA, NONOUI )
            IF( NONOUI .GT. 0 ) GOTO 10
C
C           DECOMPOSITION DE L'OCTAEDRE INTERNE A L'HEXAEDRE
            NBSOM0 = NBSOMM
            CALL SBOCTA( 0,      HEXA,   NUMOT,
     %                   MXARBO, LARBRO, MXARBT, LARBRT, NUOTPT,
     %                   NBSOMM, MXSOMM, PTXYZD,
     %                   LAVIDE, MXARET, LARETE, IERR )
C
         ELSE
C
C           TETRAEDRE : A DECOMPOSER SI INTERNE A L'HEXAEDRE
C           ================================================
            NUMOTT = - NUMOT
C           SON TYPE
            NTYOT  = LARBRT(0,NUMOTT) / 100
C           LE TETRAEDRE EST IL EXTERIEUR A L'HEXAEDRE ENGLOBANT DES POINTS?
            CALL EXTHEX( 4, LARBRT(6,NUMOTT), PTXYZD, HEXA, NONOUI )
            IF( NONOUI .GT. 0 ) GOTO 10
C
C           DECOMPOSITION DU TETRAEDRE INTERNE A L'HEXAEDRE
            NBSOM0 = NBSOMM
            CALL SBTETR( 0,      HEXA,   NTYOT,  NUMOT,
     %                   MXARBO, LARBRO, MXARBT, LARBRT, NUOTPT,
     %                   NBSOMM, MXSOMM, PTXYZD,
     %                   LAVIDE, MXARET, LARETE, IERR )
C
         ENDIF
         IF( IERR .NE. 0 ) RETURN
C
         IF( NOFOTI .GT. 0 ) THEN
C           MISE A JOUR DE TAILLE_IDEALE DES NOUVEAUX SOMMETS
            DO I=NBSOM0+1,NBSOMM
               CALL OTTAID( NOFOTI, ARETOCG, HEXA, PTXYZD(1,I),
     %                      PTXYZD(4,I), IERR )
               IF( IERR .NE. 0 ) RETURN
            ENDDO
         ENDIF
C
cccC        L'OT PERE QUI VIENT D'ETRE DECOMPOSE EST A NOUVEAU EMPILE
cccC        AFIN QUE SES FILS SOIENT EUX AUSSI EVENTUELLEMENT DECOMPOSES
ccc         LHPILE = LHPILE + 1
         GOTO 10
      ENDIF


C     TOUT OT FEUILLE D'UN POINT EST DECOMPOSE SI
C     SON ARETE EST SUPERIEURE A LA DISTANCE SOUHAITEE AUTOUR DU POINT
C     ================================================================
C     (PAS DE DO A CAUSE DE LA MISE A JOUR DE NBSOMM)
      NP = 7
C
 30   IF( NP .LE. NBSOMM ) THEN

C        ACTUELLEMENT TOUS LES POINTS INTERNES ET FRONTALIERS ONT ETE
C        AJOUTES ET LEURS NUMEROS SONT INFERIEURS A NBSOMM
ccc         IF( NPSOFR(NP) .LE. 0 ) GOTO 100

C        RECHERCHE DE L'OT FEUILLE DU POINT NP
 40      CALL NOTFPT( PTXYZD(1,NP), PTXYZD, NUOTPT(NP),
     %                LARBRO, LARBRT, NUOT )
         IF( NUOT .NE. 0 ) THEN

C           LA TAILLE DES ARETES DE L'OT FEUILLE NUOT
C           ARETE2 = CARRE DE LA LONGUEUR DE L'ARETE DE L'OT NUOT
            CALL LOAROT( NUOT, LARBRO, LARBRT, PTXYZD, ARETE2 )
            IF( ARETE2 .GT. (PTXYZD(4,NP)**2)*2 ) THEN

C              OT TROP GRAND A DECOMPOSER
               NBSOM0 = NBSOMM
               CALL SBOT( 0,      HEXA,   NUOT,
     %                    MXARBO, LARBRO, MXARBT, LARBRT, NUOTPT,
     %                    NBSOMM, MXSOMM, PTXYZD,
     %                    LAVIDE, MXARET, LARETE,
     %                    IERR )
               IF( IERR .NE. 0 ) RETURN

               IF( NOFOTI .GT. 0 .AND. NBSOM0 .LT. NBSOMM ) THEN
C                 MISE A JOUR DE TAILLE_IDEALE DES NOUVEAUX SOMMETS
                  DO I=NBSOM0+1,NBSOMM
                     CALL OTTAID( NOFOTI, ARETOCG, HEXA, PTXYZD(1,I),
     %                            PTXYZD(4,I), IERR )
                     IF( IERR .NE. 0 ) RETURN
                  ENDDO
               ENDIF
               GOTO 40

            ENDIF

         ENDIF

C        LE POINT SUIVANT A TRAITER
ccc100
         NP = NP + 1
         GOTO 30

      ENDIF

      RETURN
      END
