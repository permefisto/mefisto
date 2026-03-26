      SUBROUTINE TTAJSF( ARETGR, NOFOTI, NBVOPA, NUVOPA,
     %                   NBSOMM, MXSOMM, PTXYZD, NPSOFR,
     %                   NBFACO, MXFACO, LEFACO, N1FASC,
     %                   MXAREF, LARETF,
     %                   HEXAPAVE, NBIPAV, ECHPAV, N1SPAVE, NOPTSUIV,
     %                   NBST1SF, MNST1SF, MNTR1SF, MNIDST1SF,
     %                   QTRMIN, QTRMOY, AREMIN, AREMAX, COSMIN, COSMAX,
     %                   IERR  )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    AJOUTER LES SOMMETS DES TRIANGLES DES SURFACES CONTOURS
C -----    METTRE A JOUR NPSOFR NUMERO DE SURFACE DES SOMMETS DES FACES
C          SI UN SOMMET DES TRIANGLES LEFACO EST TRES PROCHE D'UN
C          SOMMET DES TETRAEDRES DE L'HEXAEDRE ENGLOBANT
C          ALORS CE SOMMET DEVIENT LE SOMMET LEFACO

C ENTREES:
C --------
C ARETGR : TAILLE MAXIMALE DES ARETES
C NOFOTI : NUMERO DE LA FONCTION 'TAILLE_IDEALE' DES ARETES SINON 0
C NBVOPA : NOMBRE DE  VOLUMES DU VOLUME PARTITIONNE
C NUVOPA : NUMERO DES VOLUMES DU VOLUME PARTITIONNE
C MXSOMM : NOMBRE MAXIMAL DE SOMMETS UTILISABLES DANS PTXYZD
C MXARBO : NOMBRE MAXIMAL D'OCTAEDRES DANS LARBRO
C MXARBT : NOMBRE MAXIMAL DE TETRAEDRES DANS LARBRT
C MXARET : NOMBRE MAXIMAL D'ARETES DES OT
C MXFACO : NOMBRE MAXIMAL DECLARABLE DE FACES DU CONTOUR
C MXAREF : NOMBRE MAXIMAL D'ARETES DES FACES LEFACO
C HEXAPAVE:MIN ET MAX DES COORDONNEES DU PAVAGE
C NBIPAV : NOMBRE D'ARETES DANS LA DIRECTION I DU PAVAGE
C ECHPAV : ECHELLE DANS LA DIRECTION I DU PAVAGE

C ENTREES ET SORTIES:
C -------------------
C NBSOMM : NOMBRE ACTUEL DE POINTS DECLARES DANS PTXYZD
C PTXYZD : PAR POINT : X  Y  Z  DISTANCE_SOUHAITEE
C NPSOFR :  (1 000 000) * (NO SURFACE + 1) + NUMERO DU SOMMET DANS
C          LA NUMEROTATION DES SOMMETS DE LA SURFACE
C          -4 SI SOMMET DE LA GRILLE DES TETRAEDRES
C          -1 SI SOMMET DE LA GRILLE DES TETRAEDRES ELIMINE
C
C NBFACO : NOMBRE DE FACES TRIANGULAIRES RECENSEES SUR LE CONTOUR
C          ET LES INTERFACES ENTRE VOLUMES
C LEFACO : FACE DU CONTOUR OU INTERFACES ENTRE VOLUMES
C          IL CONTIENT DANS CET ORDRE
C          NUMERO (DANS PTXYZD) DU SOMMET 1, SOMMET 2, SOMMET 3
C          NUMERO (DANS NUVOPA 0 SINON) DU VOLUME1 , VOLUME2 DE LA FACE
C          NUMERO (DANS LEFACO) DE LA FACE ADJACENTE PAR L'ARETE 1 2 3
C
C          ATTENTION: UNE ARETE PEUT APPARTENIR A PLUS DE 2 FACES
C          => CHAINAGE CIRCULAIRE DE CES FACES DANS LEFACO
C          LEFACO(9,*)  -> FACE SUIVANTE (*=0:VIDE, *<>0:NON VIDE)
C          HACHAGE AVEC LA SOMME DES 3 SOMMETS MODULO MXFACO
C          LF = MOD( NOSOFA(1)+NOSOFA(2)+NOSOFA(3) , MXFACO ) + 1
C          NF = LEFACO( 10, LF ) LE NUMERO DE LA 1-ERE FACE DANS LEFACO
C          SI LA FACE NE CONVIENT PAS => PASSAGE A LA SUIVANTE
C          NF = LEFACO( 9, NF )  ...
C          LEFACO(11,.) = NO NOTETR D'UN TETRAEDRE AYANT CETTE FACE, 0 SINON
C N1FASC : N1FASC(NS)=NUMERO (DANS LEFACO) D'UNE FACE DE SOMMET NS
C
C LARETF : NS1 NS2 ARETE SUIVANTE PREMIERE FACE LEFACO AYANT CETTE ARETE
C          HACHAGE(NS1,NS2)=(NS1 + NS2) MODULO MXAREF + 1
C N1SPAVE: NO DANS PTXYZD DU 1-ER SOMMET DE CHAQUE PAVE
C NOPTSUIV: NO DU POINT SUIVANT DANS LE CHAINAGE DES POINTS DES PAVES
C
C SORTIES:
C --------
C NBST1SF: NOMBRE DE SOMMETS IDENTIFIES DE LA PREMIERE SURFACE
C          CONVEXE ENGLOBANTE DU PREMIER VOLUME
C MNST1SF: ADRESSE MCN DU TABLEAU XYZSOMMET DE LA 1-ERE SURFACE
C MNTR1SF: ADRESSE MCN DU TABLEAU NSEF DE LA 1-ERE SURFACE
C MNIDST1SF:ADRESSE MCN DU TABLEAU MNIDST ANCIEN NO => NO POINT IDENTIFIE
C          DANS PTXYZD DE LA SURFACE FERMEE ENGLOBANTE

C QTRMIN : QUALITE MINIMALE DES TRIANGLES DES SURFACES FERMEES
C QTRMOY : QUALITE MOYENNE  DES TRIANGLES DES SURFACES FERMEES
C AREMIN : LONGUEUR DE LA PLUS PETITE ARETE DE LA TRIANGULATION INITIALE
C AREMAX : LONGUEUR DE LA PLUS GRANDE ARETE DE LA TRIANGULATION INITIALE
C COSMIN : COSINUS MINIMAL DES ANGLES DES TRIANGLES
C COSMAX : COSINUS MAXIMAL DES ANGLES DES TRIANGLES

C IERR   : 0 SI PAS D'ERREUR
C          1 SURFACE DE TYPE TRIANGLE STRUCTURE OU EXISTENCE D'UN QUADRANGLE
C          2 SATURATION DES SOMMETS
C          3 SI TRIANGULATION AVEC UNE ARETE DE LONGUEUR NULLE
C          MULTIPLE DE 10   SI TRIANGLE DEGENERE EN ARETE
C          MULTIPLE DE 100  SI TRIANGLE APPARTIENT A 3 VOLUMES
C          MULTIPLE DE 1000 SI ARETE APPARTENANT A UNE SEULE FACE
C                                    SURFACE OUVERTE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC ET ST PIERRE DU PERRAY   Avril 2008
C MODIFS : ALAIN PERRONNET LJLL UPMC ET ST PIERRE DU PERRAY   Aout  2014
C MODIFS : ALAIN PERRONNET LJLL UPMC et St PIERRE DU PERRAY   Mars  2016
C2345X7..............................................................012
      include"./incl/langue.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/gsmenu.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/a_surface__definition.inc"
      include"./incl/a_volume__definition.inc"
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
C
      DOUBLE PRECISION  PTXYZD(4,MXSOMM), DTAILL, XYZD(4)
      INTEGER           NUVOPA(1:NBVOPA),
     &                  LEFACO(1:11,0:MXFACO),
     &                  N1FASC(1:MXSOMM),
     &                  NPSOFR(1:MXSOMM),
     &                  LARETF(1:4,1:MXAREF)
      DOUBLE PRECISION  HEXAPAVE(3,2), ECHPAV(3)
      INTEGER           NBIPAV(3), N1SPAVE(*), NOPTSUIV(MXSOMM)
C
      REAL              XYZ(1:3,1:3),DIST(3),VECT(3,2)
      INTEGER           NOSOFA(1:3),NSS(2)
      EQUIVALENCE      (NSS(1),NS1),(NSS(2),NS2)
C
C     INITIALISATION DE LEFACO
C     ------------------------
C     MISE A ZERO DU PREMIER SOMMET DES FACES DU CONTOUR
C     VALEUR TEST DE NON UTILISATION DE LA FACE
C     LEFACO(9,*) EST LE CHAINAGE DES FACES VIDES
      NBFACO = 0
      DO N=0,MXFACO
         LEFACO( 1,N) = 0
         LEFACO( 6,N) = 0
         LEFACO( 7,N) = 0
         LEFACO( 8,N) = 0
         LEFACO( 9,N) = N+1
         LEFACO(10,N) = 0
         LEFACO(11,N) = 0
      ENDDO
C     CHAINAGE DES FACES VIDES DANS LEFACO
      LEFACO( 9, MXFACO) = 0
C     LA PREMIERE FACE VIDE DANS LEFACO
      LEFACO(10, 0     ) = 1
C
C     CONVERSION RADIANS DEGRES
      RADEGR = 45. / ATAN( 1. )
C
C     LE COSINUS MINIMAL ET MAXIMAL
      COSMX = COS(   8.0 / RADEGR )
      COSMN = COS( 160.0 / RADEGR )
C
      COSMAX = -2.
      COSMIN =  2.
      NBANMI =  0
      NBANMA =  0
      AREMAX = 0.
      AREMIN = RINFO( 'GRAND' )
      QTRMIN = 2
      QTRMOY = 0
C
C     BOUCLE SUR LES VOLUMES DE LA PARTITION
C     ======================================
      NBTRIA = 0
      DO 500 N = 1, NBVOPA
C
C        LE NUMERO DU VOLUME N
         NOVO = NUVOPA( N )
C        LE TABLEAU 'DEFINITION' DE CE VOLUME
         CALL LXNLOU( NTVOLU , NOVO , NTLXVO , MNDFVO )
         CALL LXTSOU( NTLXVO , 'DEFINITION' , NTDFVO , MNDFVO )
         NBSF1V = MCN( MNDFVO + WBSF1V )
         MNSV   = MNDFVO + WUSF1V - 1
C
C        BOUCLE SUR LES SURFACES FERMEES DU VOLUME N
C        ===========================================
         NBSOEF = 0
         NBFASU = 0
         MNSOF  = 0
         DO 490 IS = 1, NBSF1V

C           LE NUMERO DE LA SURFACE FERMEE IS DU CONTOUR DU VOLUME N
            NOSU = MCN( MNSV + IS )
C           LE LEXIQUE DE LA SURFACE
            CALL LXNLOU( NTSURF, NOSU, NTLXSU, MNFASU )
C           LE TABLEAU 'NSEF' DE CETTE SURFACE
            CALL LXTSOU( NTLXSU, 'NSEF', NTFASU, MNFASU )

C           LE TYPE STRUCTURE OU NON DE LA SURFACE
            NUTYMA = MCN( MNFASU + WUTYMA )
C           LE NOMBRE DE FACES DE LA SURFACE
            IF( NUTYMA .EQ. 0 ) THEN
C              SURFACE NON STRUCTUREE
               NBFASU = MCN( MNFASU + WBEFOB )
               NBSOEF = MCN( MNFASU + WBSOEF )
               MNSOF  = MNFASU + WUSOEF - 1 - NBSOEF
            ELSE IF( NUTYMA .EQ. 3 ) THEN
C              TRIANGLE STRUCTURE
               NBLGRC(NRERR) = 1
C              POUR LEVER CETTE RESTRICTION SI BESOIN EST
C              UTILISER ~/util/NSEFPA.f et NSEFNS.f
C              EN FAIT UN TRIANGLE STRUCTURE NE PEUT ETRE UNE SURFACE FERMEE
               KERR(1) ='ttajsf: TRIANGLE STRUCTURE NON PROGRAMME'
               CALL LEREUR
               IERR = 1
               GOTO 9999
            ENDIF

C           LE TABLEAU 'XYZSOMMET' DE CETTE SURFACE
            CALL LXTSOU( NTLXSU, 'XYZSOMMET', NTSOSU, MNSOSU )
C           DECLARATION DU TABLEAU D'IDENTIFICATION DES SOMMETS DANS PTXYZD
            MNIDST = 0
            NBSOSU = MCN( MNSOSU + WNBSOM )
            CALL TNMCDC( 'ENTIER', NBSOSU+1, MNIDST )
            NSTMAX = 0
C
C           IDENTIFICATION DES SOMMETS
C           --------------------------
            DO 100 NS1=1,NBSOSU
C
C              LES 3 COORDONNEES DE CE SOMMET
               MNS = MNSOSU + WYZSOM + 3 * NS1 - 3

C              LES 3 COORDONNEES DU SOMMET NS1
C              RECHERCHE DU PLUS PROCHE POINT PTXYZD DE NS1
               XYZD(1) = RMCN( MNS )
               XYZD(2) = RMCN( MNS + 1 )
               XYZD(3) = RMCN( MNS + 2 )
C
C              DISTANCE SOUHAITEE AUTOUR DE CE SOMMET NS1
               IF( NOFOTI .GT. 0 ) THEN
C                 LA FONCTION TAILLE_IDEALE(X,Y,Z) EXISTE
C                 CALCUL DE XYZD(4)
                  CALL FONVAL( NOFOTI ,3, XYZD(1), NCODEV, DTAILL )
                  IF( NCODEV .LE. 0 ) THEN
                     IERR   = 1
                     DTAILL = ARETGR
ccc                     RETURN
                  ENDIF
                  IF( DTAILL .LT. 0D0 ) THEN
                     WRITE(IMPRIM,10000) DPARAF
10000                FORMAT('ATTENTION: TAILLE_IDEALE(',
     %               G14.6,',',G14.6,',',G14.6,')<=0! => RENDUE >0' )
                     DTAILL = -DTAILL
                  ENDIF
                  IF( DTAILL .EQ. 0D0 ) THEN
                     WRITE(IMPRIM,10001) DPARAF
10001                FORMAT('ERREUR TAILLE_IDEALE(',
     %               G14.6,',',G14.6,',',G14.6,')=0!' )
                     IERR = 2
                     DTAILL = ARETGR
                  ENDIF
                  XYZD(4) = DTAILL
               ELSE
                  XYZD(4) = ARETGR
               ENDIF

C              AJOUT DU POINT INTERNE IMPOSE I A PTXYZD ET LA GRILLE
               CALL REAJSTPV(HEXAPAVE, NBIPAV, ECHPAV, N1SPAVE,NOPTSUIV,
     %                       XYZD(1), XYZD(4), NBSOMM, MXSOMM, PTXYZD,
     %                       NS )
               IF( NS .EQ. 0 ) THEN
                  NBLGRC(NRERR) = 1
                  IF( LANGAG .EQ. 0 ) THEN
                     KERR(1)='ttajsf: AUGMENTER LE NOMBRE DE SOMMETS'
                  ELSE
                     KERR(1)='ttajsf: AUGMENT the VERTICES NUMBER'
                  ENDIF
                  CALL LEREUR
                  IERR = 2
                  GOTO 450
               ENDIF

C              LE NUMERO NS DANS PTXYZD DU SOMMET NS1 DE LA SURFACE
               MCN( MNIDST + NS1 ) = NS

C              LE NUMERO DE LA SURFACE ET DU SOMMET LOCAL
               NPSOFR( NS ) = 1 000 000 * (NOSU+1) + NS1

C              RECHERCHE DU NUMERO DE SOMMET PTXYZD MAXIMAL PAR SURFACE
               IF( NS .GT. NSTMAX ) NSTMAX = NS

 100        ENDDO
C
C           BOUCLE SUR LES TRIANGLES OU FACES DE LA SURFACE FERMEE IS
C           =========================================================
            NBTRIA = NBTRIA + NBFASU
            DO 400 I=1,NBFASU

C              LE NUMERO DES SOMMETS DE LA FACE I
               MNSOF = MNSOF + NBSOEF
C
               IF( MCN(MNSOF+4) .GT. 0 ) THEN
C                 QUADRANGLE INTERDIT!
                  NBLGRC(NRERR) = 2
                  IF( LANGAG .EQ. 0 ) THEN
                KERR(1)='ttajsf: LA PEAU DU VOLUME DOIT ETRE TRIANGULEE'
                     KERR(2)='PRESENCE D''UN QUADRANGLE A SUPPRIMER'
                  ELSE
                     KERR(1)='ttajsf: The BOUNDARY MUST BE TRIANGULATED'
                     KERR(2)='PRESENCE of a QUADRANGLE to DELETE'
                  ENDIF
                  CALL LEREUR
                  IERR = 1
                  GOTO 9999
               ENDIF
C
               DO J=1,3
                  NS1 = MCN( MNSOF + J )
                  NOSOFA( J ) = NS1
C                 LES 3 COORDONNEES DE CE SOMMET
                  MNS = MNSOSU + WYZSOM + 3 * NS1 - 3
C                 LES COORDONNEES DES 3 SOMMETS
                  XYZ(1,J) = RMCN( MNS )
                  XYZ(2,J) = RMCN( MNS + 1 )
                  XYZ(3,J) = RMCN( MNS + 2 )
               ENDDO
C
C              QUALITE DU TRIANGLE
               CALL QUALEF( 3, NOSOFA, MCN(MNSOSU+WNBSOM),
     %                      RMCN(MNSOSU+WYZSOM), SURFVOLU, Q, J )
               QTRMIN = MIN( QTRMIN, Q )
               QTRMOY = QTRMOY + Q
C
C              LA LONGUEUR DES 3 COTES DE LA FACE I
               DO J=1,3
                  IF( J .EQ. 3 ) THEN
                     K=1
                  ELSE
                     K = J+1
                  ENDIF
                  DIST(J) = DIST2P( XYZ(1,J) , XYZ(1,K) )
                  IF( DIST(J) .GT. AREMAX ) AREMAX = DIST(J)
                  IF( DIST(J) .LT. AREMIN ) AREMIN = DIST(J)
               ENDDO
C
C              EXISTE-T-IL UN ANGLE > 180-10=170 DEGRES ?
               COSMA = -2.
               COSMI =  2.
               DO J=1,3
                  K  = J - 1
                  IF( K .LE. 0 ) K = K + 3
                  NS = J + 1
                  IF( NS .GT. 3 ) NS = NS - 3
C                 LE COSINUS DE L'ANGLE J
                  DO KL=1,3
                     VECT(KL,1) = XYZ(KL,K)  - XYZ(KL,J)
                     VECT(KL,2) = XYZ(KL,NS) - XYZ(KL,J)
                  ENDDO
                  COSM = PROSCR( VECT(1,1) , VECT(1,2) , 3 )
                  COSM = COSM / ( DIST(K) * DIST(J) )
                  IF( COSM .LT. COSMI ) THEN
                     COSMI = COSM
                  ENDIF
                  IF( COSM .GT. COSMA ) THEN
                     COSMA = COSM
                  ENDIF
               ENDDO
C
               IF( COSMI .LT. COSMN ) THEN
C                 => L'ANGLE EST PLUS GRAND QUE 170. DEGRES
                  NBANMA = NBANMA + 1
               ELSE IF( COSMA .GT. COSMX ) THEN
C                 => L'ANGLE EST PLUS PETIT QUE 5 DEGRES
                  NBANMI = NBANMI + 1
               ENDIF
C
               IF( Q .LT. 0.1 .OR.
     %            (COSMA .GT. COSMX .OR. COSMI .LT. COSMN) ) THEN

C                 TRAITEMENT DES ERREURS D'ARRONDIS
                  IF( COSMA .GT. 1. ) THEN
                     COSMA = 1.
                  ENDIF
                  IF( COSMI .LT. -1. ) THEN
                     COSMI = -1.
                  ENDIF

                  ACOSMI = ACOS( COSMI )
                  ACOSMA = ACOS( COSMA )

                  IF( LANGAG .EQ. 0 ) THEN
                     WRITE(IMPRIM,10309) N,IS,I,Q,
     %                     ACOSMI*RADEGR, ACOSMA*RADEGR,
     %                     ((XYZ(K,L),K=1,3),L=1,3)
                  ELSE
                     WRITE(IMPRIM,20309) N,IS,I,Q,
     %                     ACOSMI*RADEGR, ACOSMA*RADEGR,
     %                     ((XYZ(K,L),K=1,3),L=1,3)
                  ENDIF
               ENDIF
C
10309    FORMAT(' VOLUME ',I3,' SURFACE FERMEE ',I3,
     %    ' TRIANGLE',I8,' QUALITE = ',F8.5,
     %    ' ANGLE MAX = ',F7.2,' ANGLE MIN = ',F7.2,' DEGRES'/
     %    3(' X=',G15.6,' Y=',G15.6,' Z=',G15.6/))
20309    FORMAT(' VOLUME ',I3,' CLOSED SURFACE ',I3,
     %    ' TRIANGLE',I8,' QUALITY = ',F8.5,
     %    ' ANGLE MAX = ',F7.2,' ANGLE MIN = ',F7.2,' DEGREES'/
     %    3(' X=',G15.6,' Y=',G15.6,' Z=',G15.6/))
C
               COSMAX = MAX( COSMA , COSMAX )
               COSMIN = MIN( COSMI , COSMIN )

ccc               DMAX   = MAX( DIST(1), DIST(2), DIST(3) )
ccc               DMOY   = ( DIST(1) + DIST(2) + DIST(3) ) / 3

C              DISTANCE SOUHAITEE DU J-EME SOMMET DU TRIANGLE I DE LA SURFACE
C              ==============================================================
               DO  J=1,3
C                 NS1 NUMERO DU SOMMET J DE LA FACE FRONTIERE I
                  NS1 = NOSOFA(J)
C
C                 LE NUMERO DE TETRAEDRISATION DU SOMMET
                  NS = MCN( MNIDST + NS1 )
                  NOSOFA( J ) = NS

CCCC                 LA DISTANCE SOUHAITEE EST LA PLUS PETITE DES
CCCC                 PLUS GRANDES ARETES D'UNE FACE DU SOMMET NS
CCC                  IF( PTXYZD(4,NS) .LE. 0D0 ) THEN
CCC                     PTXYZD(4,NS) = DMAX
CCC                  ELSE
CCC                     PTXYZD(4,NS) = MIN( PTXYZD(4,NS), DBLE(DMAX) )
CCC                  ENDIF

cccC                 LA DISTANCE SOUHAITEE EST LA PLUS
cccC                 GRANDE DES ARETES D'UNE FACE DU SOMMET NS
ccc                  PTXYZD(4,NS) = MAX( PTXYZD(4,NS), DBLE(DMAX) )

C                 LA DISTANCE SOUHAITEE EST LE MINIMUM DE ARETGR ET DE LA
C                 MOYENNE DES LONGUEURS D'ARETES ISSUES DE CE SOMMET
                  IF( J .EQ. 1 ) THEN
                     K = 3
                  ELSE
                     K = J-1
                  ENDIF
                  DMOY = ( DIST(J) + DIST(K) ) / 2

                  IF( PTXYZD(4,NS) .LE. 0D0 ) THEN
                     PTXYZD(4,NS) = DBLE(DMOY)
                  ELSE
                     PTXYZD(4,NS) = MIN( PTXYZD(4,NS), DBLE(DMOY) )
                  ENDIF

               ENDDO
C
C              TRAITEMENT DE LA FACE I
C              =======================
               IF( NOSOFA(1) .EQ. NOSOFA(2) .OR.
     %             NOSOFA(2) .EQ. NOSOFA(3) .OR.
     %             NOSOFA(3) .EQ. NOSOFA(1) ) THEN
                   NBLGRC(NRERR) = 3
                   WRITE(KERR(MXLGER)(1:4),'(I4)') NOSU
                   WRITE(KERR(MXLGER)(5:10),'(I6)') I
                   WRITE(KERR(MXLGER-1)(1:24),'(3I8)') NOSOFA
                   IF( LANGAG .EQ. 0 ) THEN
                     KERR(1) = 'ttajsf: LE TRIANGLE'//KERR(MXLGER)(5:10)
     %                    // ' DE LA SURFACE' // KERR(MXLGER)(1:4)
                      KERR(2) =
     %               'A 2 SOMMETS IDENTIQUES apres IDENTIFICATION'
                      KERR(3) = 'SOMMETS:' // KERR(MXLGER-1)(1:24)
                   ELSE
                    KERR(1) = 'ttajsf: The TRIANGLE'//KERR(MXLGER)(5:10)
     %                    // ' of SURFACE' // KERR(MXLGER)(1:4)
                      KERR(2) =
     %               'HAS 2 IDENTICAL VERTICES after IDENTIFICATION'
                      KERR(3) = 'VERTICES:' // KERR(MXLGER-1)(1:24)
                   ENDIF
                   CALL LEREUR
                   IERR = IERR + 10
                   GOTO 400
               ENDIF
C
C              TRI CROISSANT DES 3 NUMEROS DES SOMMETS DE LA FACE
               CALL TRI3NO( NOSOFA(1), NOSOFA(1) )
C
C              LA FACE EXISTE-T-ELLE DEJA DANS LEFACO?
C              HACHAGE AVEC LA SOMME DES 3 SOMMETS MODULO MXFACO
               LF  = MOD( NOSOFA(1)+NOSOFA(2)+NOSOFA(3) , MXFACO ) + 1
C              LE NUMERO DE LA 1-ERE FACE DANS LEFACO
               NF  = LEFACO( 10, LF )
 330           IF( NF .GT. 0 ) THEN
                  IF( NOSOFA(1) .EQ. LEFACO(1,NF) ) THEN
                     IF( NOSOFA(2) .EQ. LEFACO(2,NF) ) THEN
                        IF( NOSOFA(3) .EQ. LEFACO(3,NF) ) GOTO 350
                     ENDIF
                  ENDIF
C                 LA FACE N'EST PAS RETROUVEE. PASSAGE A LA SUIVANTE
                  NF  = LEFACO(9,NF)
                  GOTO 330
               ENDIF
C
C              LA FACE N'A PAS ETE RETROUVEE. ELLE EST AJOUTEE
C              RECHERCHE D'UNE FACE VIDE DANS LEFACO
               IF( LEFACO( 10, 0 ) .LE. 0 ) THEN
                  NBLGRC(NRERR) = 1
                  IF( LANGAG .EQ. 0 ) THEN
                     KERR(1) ='SATURATION DES FACES LEFACO'
                  ELSE
                     KERR(1) ='ARRAY LEFACO SATURATED'
                  ENDIF
                  CALL LEREUR
                  IERR = 5
                  GOTO 450
               ENDIF
C
C              IL EXISTE UNE FACE LEFACO NF VIDE
               NF = LEFACO( 10, 0 )
C              MISE A JOUR DE LA PREMIERE FACE VIDE
               LEFACO( 10, 0 ) = LEFACO( 9, NF )
C
C              LA FACE NF CREEE DANS LEFACO DEVIENT LA PREMIERE DU HACHAGE
               LEFACO(  9, NF ) = LEFACO( 10, LF )
               LEFACO( 10, LF ) = NF
C
C              LA FACE EST AJOUTEE DANS LE TABLEAU LEFACO
C              LES 3 SOMMETS
               LEFACO(1,NF) = NOSOFA(1)
               LEFACO(2,NF) = NOSOFA(2)
               LEFACO(3,NF) = NOSOFA(3)
C              LE NUMERO DU PREMIER VOLUME
               LEFACO(4,NF) = N
               LEFACO(5,NF) = 0
C              MISE A JOUR DU TABLEAU NO D'UNE FACE AYANT CES SOMMETS
               N1FASC( NOSOFA(1) ) = NF
               N1FASC( NOSOFA(2) ) = NF
               N1FASC( NOSOFA(3) ) = NF
C              UNE FACE DU CONTOUR DE PLUS
               NBFACO = NBFACO + 1
               GOTO 400
C
C              FACE RETROUVEE DANS LEFACO
C              LE SECOND NUMERO DE VOLUME EST MIS A JOUR
 350           IF( LEFACO(5,NF) .GT. 0 ) THEN
                  NBLGRC(NRERR) = 2
                  WRITE(KERR(MXLGER)(1:10),'(I10)') NF
                  IF( LANGAG .EQ. 0 ) THEN
                     KERR(1) ='ttajsf:LA FACE' // KERR(MXLGER)(1:10)
                     KERR(2) ='APPARTIENT A 3 VOLUMES => INTERDIT'
                  ELSE
                     KERR(1) ='ttajsf: The TRIANGLE'//KERR(MXLGER)(1:10)
                     KERR(2) ='BELONGS to 3 VOLUMES => FORBIDDEN'
                  ENDIF
                  CALL LEREUR
                  IERR = IERR + 100
                  GOTO 400
               ENDIF
C              LE NUMERO DE 1 A NBVOPA DU 2-EME VOLUME DE CETTE FACE
               LEFACO(5,NF) =  N
C
 400        ENDDO

 450        IF( N .EQ. 1 .AND. IS .EQ. 1 ) THEN
C              NOMBRE DE SOMMETS IDENTIFIES DE LA PREMIERE SURFACE
C              DU PREMIER VOLUME (PEUT ETRE CONVEXE ENGLOBANTE)
               NBST1SF = NSTMAX
               MNST1SF = MNSOSU
C              NOMBRE DE TRIANGLES-FACES DE LA PREMIERE SURFACE
C              DU PREMIER VOLUME (PEUT ETRE CONVEXE ENGLOBANTE)
ccc               NBTR1SF = NBFASU
               MNTR1SF = MNFASU
C              TABLEAU MNIDST ANCIEN NO => NO POINT IDENTIFIE DANS PTXYZD
               MNIDST1SF = MNIDST
            ENDIF

            IF( N .NE. 1 .OR. IS .NE. 1 ) THEN
C              PROTECTION DU TABLEAU MNIDST1SF
               IF(MNIDST .GT. 0) CALL TNMCDS( 'ENTIER',NBSOSU+1,MNIDST )
            ENDIF

            IF( IERR .NE. 0 ) GOTO 9999

 490     ENDDO
 500  ENDDO

C     MISE A JOUR DU CHAINAGE DES FACES ADJACENTES PAR LES ARETES
C     ===========================================================
C     INITIALISATION DU HACHAGE DES ARETES DES FACES LEFACO
      CALL AZEROI( 4*MXAREF, LARETF )
C     LA 1-ERE ARETE LIBRE DERRIERE CELLES ADRESSEES PAR LE MINIMUM
      LAVIDF = MXAREF

      DO 590 NF=1,MXFACO
         IF( LEFACO(1,NF) .LE. 0 ) GOTO 590

C        LA FACE EXISTE . BOUCLE SUR SES 3 ARETES
         DO J=1,3
            IF( J .EQ. 3 ) THEN
               NS1 = 1
               NS2 = 3
            ELSE
               NS1 = J
               NS2 = J + 1
            ENDIF
            NS1 = LEFACO( NS1, NF )
            NS2 = LEFACO( NS2, NF )
            IF( NS1 .GT. NS2 ) THEN
               NOAR = NS1
               NS1  = NS2
               NS2  = NOAR
            ENDIF

C           HACHAGE DE L'ARETE NS1 NS2 POUR LA RETROUVER
            CALL HACHAG( 2, NSS, 4, MXAREF, LARETF, 3,  LAVIDF, NOAR )

            IF( NOAR .EQ. 0 )  THEN
C              SATURATION DES ARETES DE LARETF
               NBLGRC(NRERR) = 1
               IF( LANGAG .EQ. 0 ) THEN
                  KERR(1) ='ttajsf: TABLEAU LARETF SATURE'
               ELSE
                  KERR(1) ='ttajsf: ARRAY LARETF SATURATED'
               ENDIF
               CALL LEREUR
               IERR = 3
               GOTO 9999
            ENDIF
C
            IF( NOAR .LT. 0 ) THEN
C              L'ARETE EST CREEE . AJOUT DU NUMERO DE FACE LEFACO
               LARETF(4,-NOAR) = NF
            ELSE
C              L'ARETE EXISTE : CHAINAGE CIRCULAIRE DES FACES LEFACO
C              LA PREMIERE FACE CHAINEE
               NF0 = LARETF(4,NOAR)
C              RECHERCHE DE L'ARETE I DE NF0  DE SOMMETS NS1 NS2
               DO I=1,3
                  IF( I .EQ. 3 ) THEN
                     II = 1
                  ELSE
                     II = I + 1
                  ENDIF
                  NSS1 = LEFACO(I,NF0)
                  NSS2 = LEFACO(II,NF0)
                  IF( ( NS1 .EQ. NSS1 .AND. NS2 .EQ. NSS2 )   .OR.
     %                ( NS2 .EQ. NSS1 .AND. NS1 .EQ. NSS2 ) ) GOTO 520
               ENDDO
C
 520           NFA = LEFACO( 5+I, NF0 )
               IF( NFA .EQ. 0 ) THEN
C                 PREMIERE FACE ADJACENTE
                  LEFACO(5+I,NF0) = NF
                  LEFACO(5+J,NF ) = NF0
               ELSE
C                 DEUXIEME FACE ADJACENTE
                  LEFACO(5+I,NF0) = NF
                  LEFACO(5+J,NF ) = NFA
               ENDIF
            ENDIF
         ENDDO
 590  ENDDO

C     VERIFICATION S'IL EXISTE DES ARETES LEFACO APPARTENANT
C     A UN SEUL TRIANGLE
C     ======================================================
      DO NF=1,MXFACO
         IF( LEFACO(1,NF) .GT. 0 ) THEN
C           LA FACE EXISTE . BOUCLE SUR SES 3 ARETES
            DO J=1,3
               IF( LEFACO(5+J,NF) .EQ. 0 ) THEN
C                  UNE SEULE FACE DE LEFACO POUR CETTE ARETE J : SURFACE OUVERTE
                   IF( J .EQ. 3 ) THEN
                      NS1 = 1
                      NS2 = 3
                   ELSE
                      NS1 = J
                      NS2 = J + 1
                   ENDIF
                   NBLGRC(NRERR) = 3
                   WRITE(KERR(MXLGER)( 1:10),'(I10)') LEFACO(NS1,NF)
                   WRITE(KERR(MXLGER)(11:20),'(I10)') LEFACO(NS2,NF)
                   IF( LANGAG .EQ. 0 ) THEN
                      KERR(1) = 'ttajsf: L''ARETE' // KERR(MXLGER)(1:20)
                      KERR(2) = 'APPARTIENT A UN SEUL TRIANGLE'
                      KERR(3) = 'SURFACE OUVERTE A REFERMER'
                   ELSE
                      KERR(1) = 'ttajsf: The EDGE' // KERR(MXLGER)(1:20)
                      KERR(2) = 'BELONGS TO ONE AND ONLY ONE TRIANGLE'
                      KERR(3) = 'OPENED SURFACE TO CLOSE'
                   ENDIF
                   CALL LEREUR
                   WRITE(IMPRIM,*) 'FACE OUVERTE ',NF,' :',
     %                             (LEFACO(K,NF),K=1,9)
                   IERR = IERR + 1000
               ENDIF
            ENDDO
         ENDIF
      ENDDO
C
C     BILAN SUR LA QUALITE DE LA TRIANGULATION DES SURFACES FERMEES
C     =============================================================
      QTRMOY = QTRMOY / NBFACO
      IF( NBANMI+NBANMA .GT. 0 ) THEN
         NBLGRC(NRERR) = 7
         WRITE(KERR(MXLGER)(1:10),'(F10.6)') QTRMIN
         WRITE(KERR(MXLGER-1)(1:10),'(F10.6)') QTRMOY
         WRITE(KERR(MXLGER-2)(1:10),'(I10)') NBANMI
         WRITE(KERR(MXLGER-3)(1:10),'(I10)') NBANMA
         IF( LANGAG .EQ. 0 ) THEN

         KERR(1)=KERR(MXLGER)(1:10)//' QUALITE MINIMALE DES TRIANGLES'
         KERR(2)=KERR(MXLGER-1)(1:10)//' QUALITE MOYENNE  DES TRIANGLES'
         KERR(3)=KERR(MXLGER-2)(1:10)//' ANGLES TROP PETITS'
         KERR(4)=KERR(MXLGER-3)(1:10)//' ANGLES TROP GRANDS'
         KERR(5)='OPTION 30 DES SURFACES A UTILISER'
         KERR(6)='AFIN DE FAVORISER LA TETRAEDRISATION'
         KERR(7)='ATTENTION: TETRAEDRISATION NON GARANTIE'

         ELSE

         KERR(1)=KERR(MXLGER)(1:10)//' TRIANGLES QUALITY MINIMUM'
         KERR(2)=KERR(MXLGER-1)(1:10)//' TRIANGLES QUALITY AVERAGE'
         KERR(3)=KERR(MXLGER-2)(1:10)//' ANGLES TOO SMALL'
         KERR(4)=KERR(MXLGER-3)(1:10)//' ANGLES TOO GREAT'
         KERR(5)='=> USE the OPTION 30 of SURFACES'
         KERR(6)='TO FAVORIZE a GOOD TETRAHEDRIZATION'
         KERR(7)='ATTENTION: TETRAHEDRIZATION NOT GUARANTED'

         ENDIF
         CALL LERESU
      ENDIF

      IF( LANGAG .EQ. 0 ) THEN
      WRITE(IMPRIM,10510) NBTRIA,QTRMIN,QTRMOY,
     %                    ACOS(COSMAX)*RADEGR,NBANMI,ACOS(COSMX)*RADEGR,
     %                    ACOS(COSMIN)*RADEGR,NBANMA,ACOS(COSMN)*RADEGR,
     %                    AREMIN, AREMAX
      ELSE
      WRITE(IMPRIM,20510) NBTRIA,QTRMIN,QTRMOY,
     %                    ACOS(COSMAX)*RADEGR,NBANMI,ACOS(COSMX)*RADEGR,
     %                    ACOS(COSMIN)*RADEGR,NBANMA,ACOS(COSMN)*RADEGR,
     %                    AREMIN, AREMAX
      ENDIF
10510 FORMAT(/' NOMBRE DE TRIANGLES DE LA PEAU =',I7/
     %        ' QUALITE DE LA TRIANGULATION DE LA PEAU DE L''OBJET'/
     %        ' QUALITE MINIMALE DES TRIANGLES =',G15.6/
     %        ' QUALITE MOYENNE  DES TRIANGLES =',G15.6/
     %        ' ANGLE MINIMAL    DES TRIANGLES =',F7.2,' DEGRES  et',
     %        I6,' ANGLES < ',F7.2,' DEGRES' /
     %        ' ANGLE MAXIMAL    DES TRIANGLES =',F7.2,' DEGRES  et',
     %        I6,' ANGLES > ',F7.2,' DEGRES' /
     %        ' LONGUEUR DE LA PLUS PETITE ARETE =',G15.6/
     %        ' LONGUEUR DE LA PLUS GRANDE ARETE =',G15.6)
20510 FORMAT(/' BOUNDARY TRIANGLES NUMBER =',I7/
     %        ' BOUNDARY TRIANGLES QUALITY MINIMUM',G15.6/
     %        ' BOUNDARY TRIANGLES QUALITY AVERAGE',G15.6/
     %        ' TRIANGLES ANGLE MINIMUM =',F7.2,' DEGREES  and',
     %        I6,' ANGLES < ',F7.2,' DEGREES' /
     %        ' TRIANGLES ANGLE MAXIMUM =',F7.2,' DEGREES  and',
     %        I6,' ANGLES > ',F7.2,' DEGREES' /
     %        ' LENGTH of the SMALLEST EDGE =',G15.6/
     %        ' LENGTH of the GREATEST EDGE =',G15.6)

      IF( AREMIN .LE. 0 ) IERR = 3

 9999 RETURN
      END
