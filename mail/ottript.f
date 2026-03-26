      SUBROUTINE OTTRIPT( NU1RPT, NBSOMM, NPSOFR, LARBRO, LARBRT,
     %                    MXQUEU, MNQUEU, NBPOIN, NUPOIN, NUTEST )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRI DES POINTS ET SOMMETS D'OT POUR LES TETRAEDRISER
C -----    LES SOMMETS D'OT SONT PLACES EN PREMIER SELON LA TAILLE DES OT
C          DES PLUS GRANDS AUX PLUS PETITS

C ENTREES:
C --------
C NU1RPT : NUMERO DU PREMIER POINT STOCKE DANS PTXYZD
C NBSOMM : NOMBRE ACTUEL DE POINTS ACTIFS DANS PTXYZD
C NPSOFR : (I) =  0  SI POINT AJOUTE NON SUR LA SURFACE DE L'OBJET
C              =  1 SI LE POINT EST AJOUTE SUR UNE ARETE DE LA SURFACE
C              =  2 SI LE POINT EST AJOUTE COMME BARYCENTRE SUR LA SURFACE
C              =  1 000 000  + NUMERO DU  POINT INTERNE UTILISATEUR
C              = (1 000 000) * (NO SURFACE + 1) + NUMERO DU SOMMET DANS
C                   LA NUMEROTATION DES SOMMETS DE LA SURFACE
C              = -4 SI LE POINT EST EXTERIEUR A L'OBJET
C              = -NPSOFR(I) SI POINT DEPLACE SUR LA SURFACE
C LARBRO : ARBRE-14 DES OCTAEDRES ( FOND DE LA TETRAEDRISATION )
C      LARBRO(0,0) : NO DU 1-ER OCTAEDRE VIDE DANS LARBRO
C      LARBRO(1,0) : MAXIMUM DU 1-ER INDICE DE LARBRO (ICI -1:20)
C      LARBRO(2,0) : MAXIMUM DECLARE DU 2-EME INDICE DE LARBRO
C                    (ICI = MXARBO)

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

C      LARBRO(15:20,J) : NO PTXYZD DES 6 SOMMETS DE L'OCTAEDRE J
C   SINON
C      LARBRO(0,J): -ADRESSE DANS LARBRO DE L'OCTAEDRE VIDE SUIVANT

C LARBRT : ARBRE-5 DES TETRAEDRES ( FOND DE LA TETRAEDRISATION )
C      LARBRT(0,0) : NO DU 1-ER TETRAEDRE VIDE DANS LARBRT
C      LARBRT(1,0) : MAXIMUM DU 1-ER INDICE DE LARBRT (ICI -1:9)
C      LARBRT(2,0) : MAXIMUM DECLARE DU 2-EME INDICE DE LARBRT
C                     (ICI = MXARBT)

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

C      LARBRT(6:9,J) : NO PTXYZD DES 4 SOMMETS DU TETRAEDRE J
C   SINON
C      LARBRT(0,J): ADRESSE DANS LARBRT DU TETRAEDRE VIDE SUIVANT

C MXQUEU : NOMBRE D'ENTIERS DE LA QUEUE . CE NOMBRE DOIT ETRE PAIR
C MNQUEU : ADRESSE MCN DU TABLEAU DE LA QUEUE

C SORTIES:
C --------
C NBPOIN : NOMBRE DE POINTS A TETRAEDRISER
C NUPOIN : NUMERO DANS PTXYZD DES POINTS A TETRAEDRISER(0:MXSOMM AU PLUS)

C AUXILIAIRE:
C -----------
C NUTEST : TABLEAU AUXILIAIRE DE NBSOMM ENTIERS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC   NOVEMBRE 1992
C2345X7..............................................................012
      include"./incl/langue.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)
      INTEGER           LARBRO(-1:20,0:*),
     %                  LARBRT(-1:9 ,0:*),
     %                  NPSOFR(1:*),
     %                  NUPOIN(0:*),
     %                  NUTEST(1:NBSOMM)

C     INITIALISATION A ZERO DU TABLEAU NUTEST
      DO NS=1,NBSOMM
         NUTEST( NS ) = 0
      ENDDO
      NUPOIN(0) = 0

C     LES SOMMETS D'OT SONT PLACES EN PREMIER SELON LA TAILLE DES OT
C     DES PLUS GRANDS AUX PLUS PETITS
C     ==============================================================
      NBPOIN = 0
      LHQUEU = 1
C     L'OT RACINE EST MISE DANS LA QUEUE
      MCN( MNQUEU ) = 1
C     L'ADRESSE MCN DE LA FIN DE QUEUE
      MNQ    = MNQUEU
C     LA FIN DE L'ANNEAU DANS MCN POUR LA QUEUE CIRCULAIRE
      MNQFIN = MNQUEU + MXQUEU - 1

C     TANT QUE LA QUEUE N'EST PAS VIDE FAIRE
C     ======================================
 10   IF( LHQUEU .GT. 0 ) THEN

C        L'OT A TRAITER
         NUOT   = MCN( MNQ )

C        L'OT EST RETIRE DE LA QUEUE
         LHQUEU = LHQUEU - 1
         MNQ    = MNQ + 1
         IF( MNQ .GT. MNQFIN ) THEN
C           FIN DE L'ANNEAU
            MNQ = MNQUEU
         ENDIF

C        LES FILS DE L'OT SONT MIS DANS LA QUEUE
         IF( NUOT .GT. 0 ) THEN

C           OCTAEDRE
            IF( LARBRO(1,NUOT) .GT. 0 ) THEN
C              LES 14 FILS SONT MIS DANS LA QUEUE
               IF( LHQUEU + 14 .GT. MXQUEU ) THEN
C                 SATURATION DE LA QUEUE => ELLE EST AGRANDIE
                  CALL TNMCDC( 'ENTIER', MXQUEU+2048, MNQQ )
C                 LA QUEUE REDEMARRE AU DEBUT DU TABLEAU
                  DO 15 I=0,LHQUEU-1
                     MN = MNQ + I
                     IF( MN .GT. MNQFIN ) MN = MN - MXQUEU
                     MCN(MNQQ+I) = MCN(MN)
 15               CONTINUE
                  CALL TNMCDS( 'ENTIER', MXQUEU, MNQUEU )
                  MNQUEU = MNQQ
                  MXQUEU = MXQUEU + 2048
                  MNQFIN = MNQUEU + MXQUEU - 1
C                 LE POINTEUR SUR LA FIN DE QUEUE
                  MNQ    = MNQUEU
               ENDIF
               MN = MNQ + LHQUEU - 1
               DO 20 I=1,14
                  MN = MN + 1
                  IF( MN .GT. MNQFIN ) MN = MN - MXQUEU
                  MCN( MN ) = LARBRO(I,NUOT)
 20            CONTINUE
               LHQUEU = LHQUEU + 14
            ENDIF

C           TRAITEMENT DES 6 SOMMETS DE L'OT
            DO 22 I=15,20
               NS = LARBRO(I,NUOT)
               IF( NPSOFR( NS ) .EQ. -1 ) GOTO 22
               IF( NUTEST( NS ) .EQ. 0 ) THEN
C                 AJOUT SANS REPETITION DU SOMMET NS DANS LE TABLEAU NUPOIN
                  NBPOIN = NBPOIN + 1
                  NUPOIN( NBPOIN ) = NS
                  NUTEST(   NS   ) = 1
               ENDIF
 22         CONTINUE

         ELSE

C           TETRAEDRE
            NUOTT  = -NUOT
            IF( LARBRT(1,NUOTT) .GT. 0 ) THEN
C              LES 5 FILS SONT MIS DANS LA QUEUE
               IF( LHQUEU + 5 .GT. MXQUEU ) THEN
C                 SATURATION DE LA QUEUE => ELLE EST AGRANDIE
                  CALL TNMCDC( 'ENTIER', MXQUEU+2048, MNQQ )
                  DO 25 I=0,LHQUEU-1
                     MN = MNQ + I
                     IF( MN .GT. MNQFIN ) MN = MN - MXQUEU
                     MCN(MNQQ+I) = MCN(MN)
 25               CONTINUE
                  CALL TNMCDS( 'ENTIER', MXQUEU, MNQUEU )
                  MNQUEU = MNQQ
                  MXQUEU = MXQUEU + 2048
                  MNQFIN = MNQUEU + MXQUEU - 1
C                 LE POINTEUR SUR LA FIN DE QUEUE
                  MNQ    = MNQUEU
               ENDIF
               MN = MNQ + LHQUEU - 1
               DO 30 I=1,5
                  MN = MN + 1
                  IF( MN .GT. MNQFIN ) MN = MN - MXQUEU
                  MCN( MN ) = LARBRT(I,NUOTT)
 30            CONTINUE
               LHQUEU = LHQUEU + 5
            ENDIF

C           TRAITEMENT DES 4 SOMMETS DE L'OT
            DO 32 I=6,9
               NS = LARBRT(I,NUOTT)
               IF( NPSOFR( NS ) .EQ. -1 ) GOTO 32
               IF( NUTEST( NS ) .EQ. 0 ) THEN
C                 AJOUT SANS REPETITION DU SOMMET NS DANS LE TABLEAU NUPOIN
                  NBPOIN = NBPOIN + 1
                  NUPOIN( NBPOIN ) = NS
                  NUTEST(   NS   ) = 1
               ENDIF
 32         CONTINUE
         ENDIF

C        RETOUR EN TETE DE QUEUE
         GOTO 10
      ENDIF
      WRITE(IMPRIM,*)
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,*) 'DERNIER SOMMET D''OT       =',NBPOIN
      ELSE
         WRITE(IMPRIM,*) 'LAST VERTEX of OT          =',NBPOIN
      ENDIF

C     STOCKAGE DES POINTS INTERNES ET DES SOMMETS DE FACES DES SURFACES FERMEES
      DO 190 NS=NU1RPT,NBSOMM
         IF( NPSOFR(NS) .GT. 0 ) THEN
C           AJOUT SANS REPETITION DU POINT NS DANS LE TABLEAU NUPOIN
C           EN EFFET UN POINT PEUT ETRE UN SOMMET D'OT
            IF( NUTEST( NS ) .EQ. 0 ) THEN
C              AJOUT SANS REPETITION DU SOMMET NS DANS LE TABLEAU NUPOIN
               NBPOIN = NBPOIN + 1
               NUPOIN( NBPOIN ) = NS
               NUTEST(   NS   ) = 1
            ENDIF
         ENDIF
 190  CONTINUE
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,*) 'DERNIER POINT SUR LA PEAU =',NBPOIN
      ELSE
         WRITE(IMPRIM,*) 'LAST POINT on the BOUNDARY =',NBPOIN
      ENDIF

CCCC     ESSAI DES POINTS SUPPRIMES A DETRUIRE SI LA TETRAEDRISATION
CCCC     ALENTOUR EST INCORRECTE
CCC      DO 300 NS=1, NBSOMM
CCC         IF( NPSOFR(NS) .EQ. -1 ) THEN
CCCC           POINT SUPPRIME QUI EST A NOUVEAU AJOUTE
CCC            IF( NUTEST( NS ) .EQ. 0 ) THEN
CCCC              AJOUT SANS REPETITION DU SOMMET NS DANS LE TABLEAU NUPOIN
CCC               NBPOIN = NBPOIN + 1
CCC               NUPOIN( NBPOIN ) = NS
CCC               NUTEST(   NS   ) = 1
CCC            ENDIF
CCC         ENDIF
CCC 300  CONTINUE
CCC      WRITE(IMPRIM,*) 'DERNIER POINT TROP PROCHE =',NBPOIN

      RETURN
      END
