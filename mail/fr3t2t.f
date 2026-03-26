      SUBROUTINE FR3T2T( PTXYZD, NBFAPE, NOFAPE,
     %                   MXFACO, LEFACO, N1FASC,
     %                   N1TETS, N1TEVI, NOTETR, NUDTETR,
     %                   MXTETRA1S, NOTETRA1S,
     %                   MXTECF, NOTECF,
     %                   NBFRE,  IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : ECHANGE DES 3TETRAEDRES D'ARETE COMMUNE EN 2 TETRAEDRES DE FACE COMMUNE
C -----
C ENTREES:
C --------
C PTXYZD : TABLEAU DES COORDONNEES DES POINTS
C          PAR POINT : X  Y  Z DISTANCE_SOUHAITEE
C NBFAPE : NOMBRE DE FACES PERDUES DANS LEFACO
C NOFAPE : NO DES FACES LEFACO PERDUES
C MXFACO : MAX DE FACES DECLARABLES DANS LEFACO
C LEFACO : FACE=TRIANGLE DE LA PEAU OU DES INTERFACES ENTRE VOLUMES
C          IL CONTIENT DANS CET ORDRE
C          NUMERO (DANS PTXYZD) DU SOMMET 1, SOMMET 2, SOMMET 3
C          NUMERO (DANS NUVOPA 0 SINON) DU VOLUME1, VOLUME2 DE LA FACE
C          NUMERO (DANS LEFACO) DE LA FACE ADJACENTE PAR L'ARETE 1 2 3
C
C          ATTENTION: UNE ARETE PEUT APPARTENIR A PLUS DE 2 FACES
C          => CHAINAGE CIRCULAIRE DE CES FACES DANS LEFACO
C          LEFACO(9,*)  -> FACE SUIVANTE (*=0:VIDE, *<>0:NON VIDE)
C          HACHAGE AVEC LA SOMME DES 3 SOMMETS MODULO MXFACO
C          LF = MOD( NOSOFA(1)+NOSOFA(2)+NOSOFA(3) , MXFACO ) + 1
C          NF = LEFACO( 10, LF ) LE NUMERO DE LA 1-ERE FACE DANS LEFACO
C          SI LA FACE NE CONVIENT PAS. PASSAGE A LA SUIVANTE
C          NF  = LEFACO( 9, NF )  ...
C          LEFACO(10,*) PREMIERE FACE DANS LE HACHAGE
C          LEFACO(11,.) = NO NOTETR D'UN TETRAEDRE AYANT CETTE FACE, 0 SINON
CCCC       LEFACO(12,.) = NO FACEOC DE 1 A NBFACES D'OC
C N1FASC : N1FASC(I) NUMERO D'UN TRIANGLE DE LEFACO DE SOMMET I
C N1TETS : N1TETS(I) NUMERO D'UN TETRAEDRE AYANT POUR SOMMET I
C N1TEVI : NO DU PREMIER TETRAEDRE VIDE DANS NOTETR
C NOTETR : LISTE DES TETRAEDRES
C          SOMMET1,    SOMMET2,    SOMMET3,    SOMMET4,
C          TETRAEDRE1, TETRAEDRE2, TETRAEDRE3, TETRAEDRE4
C          DE L'AUTRE COTE DE LA FACE
C          1: 123      2: 234      3: 341      4: 412
C
C MXTETRA1S: NOMBRE MAX DE TETRAEDRES D'UN SOMMET DE LA TETRAEDRISATION
C NOTETRA1S: NO DES TETRAEDRES D'UN SOMMET DE LA TETRAEDRISATION
C MXTECF : NOMBRE MAX DE TETRAEDRES     D'UN CF LEFACO
C NOTECF : NUMERO NOTETR DES TETRAEDRES D'UN CF LEFACO
C
C SORTIES :
C ---------
C NBFRE  : NOMBRE DE FACES PERDUES RECUPEREES
C IERR   : 0 SI PAS D'ERREUR , 1 SATURATION D'UN TABLEAU
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LABORATOIRE JLL PARIS UPMC          MARS 2006
C2345X7..............................................................012
      include"./incl/langue.inc"
      PARAMETER        (MXETOI=1000, MX1ACF=100, MXARCF=1000, MXNFNF=16)
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      DOUBLE PRECISION  PTXYZD(1:4,1:*)
      INTEGER           NOTETR(8,1:*), N1TETS(1:*)
      INTEGER           NOFAPE(NBFAPE),
     %                  LEFACO(11,0:MXFACO),
     %                  N1FASC(*)
      INTEGER           NOTETRA1S(MXTETRA1S)
      INTEGER           NOTECF(MXTECF)
      DOUBLE PRECISION  PT(3), COBARY(3),
     %                  A, X21, Y21, Z21
C
      INTEGER           NOSOAR(2,6)
      DATA              NOSOAR / 1,2, 2,3, 3,1, 1,4, 2,4, 3,4 /
C
C     NOMBRE DE FACES RETROUVEES
      NBFRE = 0
      NFP   = 0
C
 1    IERR = 0
      NFP  = NFP + 1
      IF( NFP .LE. NBFAPE ) THEN
C
C        LE NUMERO DE LA FACE PERDUE DANS LEFACO
         NF0 = NOFAPE( NFP )
         IF( NF0 .LE. 0 ) GOTO 1
C
CCC         PRINT *,'====================================================='
CCC         PRINT *,'fr3t2t: FACE PERDUE',NFP,' LEFACO(',NF0,')  PERFOREE?'
C
C        RECENSENSEMENT DES SOMMETS DES ARETES PERFORANT UNE FACE LEFACO
C        LISTE DES TETRAEDRES AYANT UNE ARETE COMMUNE PERFORANTE
         DO 90 NS=1,3
C           LE NO DE SOMMET NS DE LA FACE PERDUE
            NSF = LEFACO(NS,NF0)
C
C           RECHERCHE DES NBTE TETRAEDRES DE SOMMET NSF
            CALL TETR1S( NSF,  N1TETS,    NOTETR,
     %                   NBTE, MXTETRA1S, NOTETRA1S, IERR )
            IF( IERR .NE. 0 ) GOTO 1

            DO 85 NT=1,NBTE
C              UN TETRAEDRE DE SOMMET NSF
               NTE = NOTETRA1S( NT )
C
C              PARCOURS DES 6 ARETES DU TETRAEDRE NTE
               DO 80 NA = 1,6
C
C                 LES 2 SOMMETS DE L'ARETE NA DU TETRAEDRE NTE
                  NS1 = NOTETR( NOSOAR(1,NA), NTE )
                  IF( NS1 .EQ. NSF ) GOTO 80
                  NS2 = NOTETR( NOSOAR(2,NA), NTE )
                  IF( NS2 .EQ. NSF ) GOTO 80
C
C                 L'ARETE NS1-NS2 N'A PAS LE SOMMET NSF
C                 CALCUL DU POINT D'INTERSECTION DE NS1-NS2 ET DU TRIANGLE NF0
                  CALL INDRPL( PTXYZD(1,NS1), PTXYZD(1,NS2),
     %                         PTXYZD(1,LEFACO(1,NF0)),
     %                         PTXYZD(1,LEFACO(2,NF0)),
     %                         PTXYZD(1,LEFACO(3,NF0)),
     %                         PT, IERR )
                  IF( IERR .NE. 0 ) THEN
C                    S1=S2 OU S1-S2 PARALLELE A LA FACE FRONTIERE
                     IERR = 0
                     GOTO 80
                  ENDIF
C
C                 IL EXISTE UN POINT D'INTERSECTION S1-S2
C                 AVEC LA FACE FRONTIERE PERDUE NF0
C                 LE POINT D'INTERSECTION EST IL DANS LE TRIANGLE NF0?
                  CALL PTDSTR(PT,PTXYZD,LEFACO(1,NF0),COBARY,NONOUI)
                  IF( NONOUI .GT. 0 ) THEN
C                    OUI: PT EST DANS OU SUR LE TRIANGLE
C                         PT EST IL ENTRE NS1 ET NS2 ?
                     X21 = PTXYZD(1,NS2) - PTXYZD(1,NS1)
                     Y21 = PTXYZD(2,NS2) - PTXYZD(2,NS1)
                     Z21 = PTXYZD(3,NS2) - PTXYZD(3,NS1)
                     A = ( (PT(1)-PTXYZD(1,NS1)) * X21
     %                   + (PT(2)-PTXYZD(2,NS1)) * Y21
     %                   + (PT(3)-PTXYZD(3,NS1)) * Z21 )
     %                   / ( X21 * X21 + Y21 * Y21 + Z21 * Z21 )
                     IF( 1D-3 .LT. A .AND. A .LT. 0.999D0 ) THEN
C                       OUI: PT EST ENTRE NS1 ET NS2 ET DANS LE TRIANGLE FRONTIE
                        DO 55 MM=1,3
                           IF( 0.999D0 .LE. COBARY(MM) .AND.
     %                         COBARY(MM) .LE. 1.0001D0 ) THEN
C                             PT EST UN SOMMET DU TRIANGLE OU TRES PRES
ccc                              WRITE(IMPRIM,*)
ccc     %                    'fr3t2t: PT TROP PROCHE SOMMET DU TRIANGLE A='
ccc     %                    ,A,' COBARY=',COBARY
                              GOTO 1
                           ENDIF
 55                     CONTINUE
C
C                       OUI: PT EST ENTRE NS1 ET NS2 ET DANS LE TRIANGLE FRONTIE
C                            ET N'EST PAS TROP PRES D'UN SOMMET DU TRIANGLE
C                            => LE TETRAEDRE EST JUGE INTERSECTANT LE CF
C
C                       RECHERCHE DES TETRAEDRES D'ARETE NS1-NS2
                        CALL TETR1A( NS1,    NS2,    N1TETS, NOTETR,
     %                               NBTECF, MXTECF, NOTECF, IERR )
                        IF( NBTECF .NE. 3 ) GOTO 1
C
C                       UNE ARETE COMMUNE A 3TE => 2TE AVEC LA FACE COMMUNE NF0
                        CALL TE3TE2( NF0, LEFACO, N1FASC, NOTECF,
     %                              PTXYZD,N1TETS,NOTETR,N1TEVI,NUDTETR,
     %                              IERR )
                        IF( IERR .NE. 0 ) GOTO 1
C
C                       ECHANGE CORRECT UNE FACE RETROUVEE DE PLUS
                        NBFRE = NBFRE + 1
ccc        WRITE(IMPRIM,*)'fr3t2t: 3TE => 2TE FACE ',NF0,
ccc     %  ' POUR LES ST',NS1,NS2
C                       LA FACE NF0=NOTRCF(M) LEFACO TRAITEE N'EST PLUS PERDUE
                        NOFAPE( NFP ) = -NF0
                        GOTO 1
                     ENDIF
                  ENDIF
 80            CONTINUE
 85         CONTINUE
 90      CONTINUE
         GOTO 1
      ENDIF

      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10090) NBFRE
10090    FORMAT(' fr3t2t :',I7,' FACES PERDUES RETROUVEES' )
      ELSE
         WRITE(IMPRIM,20090) NBFRE
20090   FORMAT(' fr3t2t  :',I7,' LOST FACES HAVE BEEN RECOVERED')
      ENDIF

      RETURN
      END
