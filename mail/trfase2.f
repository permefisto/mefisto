      SUBROUTINE TRFASE2( NBSOM,  XYZSOM, MNNSEF, NBFATR, NOFATR,
     %                    NBFASE, NOSTFA, NUEFFA, DIOEIL,
     %                    MOFACE, LFACES,
     %                    NBFACT, XYZPLAN )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACER LES FACES FRONTALIERES ET LES FACES DES EF AYANT UNE
C -----    INTERSECTION AVEC LE PLAN NOAXE DE COORDONNEE COPLAN ET
C          LES 4 ARETES DU RECTANGLE DE VISUALISATION DU PLAN DE SECTION
C
C ENTREES:
C --------
C NBSOM  : NOMBRE DE SOMMETS DU MAILLAGE DU VOLUME
C XYZSOM : 3 COORDONNEES DES SOMMETS DU MAILLAGE
C MNNSEF : ADRESSE MCN DU TMS NSEF DU VOLUME
C NBFATR : NOMBRE DE FACES A TRACER (FRONTALIERES OU SECTIONNEES)
C NOFATR : NO>0 DE LA FACE DANS LFACES DES FACES FRONTALIERES
C            <0 DE LA FACE DANS NOSTFA DES FACES SECTIONNEES
C NBFASE : NOMBRE DE FACES DES EF SECTIONNES PAR LE PLAN
C NOSTFA : NO DES SOMMETS DANS XYZSOM  DES FACES DES EF SECTIONNES
C NUEFFA : NO DE L'EF 3D DES FACES SECTIONNEES
C DIOEIL : DISTANCE DU BARYCENTRE A LA FACE (FRONTALIERE OU SECTIONNEE)
C MOFACE : NOMBRE DE MOTS PAR FACE DU TABLEAU LFACES DES FACES
C LFACES : TABLEAU NUMERO DES FACES, LIEN, NO DES CUBES, ...
C          LFACES(1,I)= NO DU 1-ER  SOMMET DE LA FACE
C          LFACES(2,I)= NO DU 2-EME SOMMET > 1-ER  SOMMET
C          LFACES(3,I)= NO DU 3-EME SOMMET DE LA FACE
C          LFACES(4,I)= NO DU 4-EME SOMMET DE LA FACE
C                       0 SI TRIANGLE
C NBFACT : NOMBRE DE FACES FRONTALIERES + DES FACES DES EF SECTIONNES
C XYZPLAN: 3 COORDONNEES DES 4 SOMMETS DU RECTANGLE DU PLAN DE SECTION
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Laboratoire JL LIONS UPMC Paris Decembre 2006
C....................................................................012
      include"./incl/xyzext.inc"
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
C
      include"./incl/pp.inc"
      COMMON   MCN(MOTMCN)
      REAL     XYZSOM(3,NBSOM)
      INTEGER  NOFATR(NBFATR), NOSOEF(64)
      INTEGER  NOSTFA(4,NBFASE), NUEFFA(NBFASE)
      REAL     DIOEIL(NBFATR)
      INTEGER  LFACES(MOFACE,*)
      REAL     XYZPLAN(3,4)
      REAL     XYZ(3,4), XYZP(3), CNORFA(3)
      CHARACTER*10  NMSOMM
C
C     LES CARACTERISTIQUES DES EF DE CE MAILLAGE INITIAL
      CALL NSEFPA( MCN(MNNSEF),
     %             NUTYMA, NBSOEL, NBSOEF, NBTGEF,
     %             LDAPEF, LDNGEF, LDTGEF, NBEFOB,
     %             NX    , NY    , NZ    ,
     %             IERR   )
C     NUTYMA : 'NUMERO DE TYPE DU MAILLAGE'    ENTIER
C              0 : 'NON STRUCTURE'      , 2 : 'SEGMENT    STRUCTURE',
C              3 : 'TRIANGLE  STRUCTURE', 4 : 'QUADRANGLE STRUCTURE',
C              5 : 'TETRAEDRE STRUCTURE', 6 : 'PENTAEDRE  STRUCTURE',
C              7 : 'HEXAEDRE  STRUCTURE'
C     NBSOEL : NOMBRE DE SOMMETS DES EF
C              0 SI MAILLAGE NON STRUCTURE
C     NBSOEF : NOMBRE DE SOMMETS DE STOCKAGE DES EF
C              ( TETRAEDRE NBSOEL=4  NBSOEF=8 )
C     NBEFOB : NOMBRE DE EF DU MAILLAGE
C     NX, NY, NZ : LE NOMBRE D'ARETES DANS LES DIRECTION X Y Z
C                CF LE TMS ~td/d/a___nsef
C     PARCOURS DES FACES DES EF VOLUMIQUES DU VOLUME
C
      IF ( LCRITR .GT. 0 ) THEN
C        PALETTE ARC EN CIEL DES QUALITES DES EF
         CALL PALCDE( 12 )
      ELSE
C        PALETTE DES GRIS
         CALL PALCDE( 10 )
      ENDIF
C
C     LA DIRECTION DE VISEE PTV-OEIL
      DIREVI(1) = AXOEIL(1) - AXOPTV(1)
      DIREVI(2) = AXOEIL(2) - AXOPTV(2)
      DIREVI(3) = AXOEIL(3) - AXOPTV(3)
C
C     DISTANCES MIN ET MAX DES FACES A L'OEIL
      DMIN = DIOEIL(1)
      DMAX = DIOEIL(NBFATR)
C
C     CALCUL DES COULEURS DES FACES AVEC PRISE EN COMPTE
C     DE LA DIRECTION DE VISEE ET L'ELOIGNEMENT
C     POID : POIDS DE LA DIRECTION DE VISEE DANS CE CALCUL
C            ATTENTION: 0 < POID < 1
      EP1   = 0.05
      EP2   = SQRT(EP1)
      EP3   = SQRT(1+EP1) - EP2
      DELTA = DMAX - DMIN
      IF ( DELTA .EQ. 0 ) THEN
         DELTA = 1.
         POID  = 1.
      ELSE
         POID  = 0.6
      ENDIF
      POID1 = 1. - POID
C
C     REDUCTION DES FACES
      REDUCF = PREDUF * 0.01
      REDUC1 = 1.0 - REDUCF
C
C     NOMBRE DE COULEURS DISPONIBLES
      NBCOUL = NDCOUL - N1COUL + 1
C
      DO 500 NF = 1, NBFATR
C
C        LE NUMERO DE LA FACE LA PLUS ELOIGNEE NON TRACEE
         NOF = NOFATR( NF )
C
         IF( NOF .GE. 0 ) THEN
C
C           FACE FRONTALIERE A TRACER EN FIL DE FER
C           =======================================
C           LE NOMBRE DE SOMMETS DE LA FACE
            IF( LFACES(4,NOF) .GT. 0 ) THEN
               NBSF = 4
            ELSE
               NBSF = 3
            ENDIF
C
            DO 10 J=1,NBSF
C              LES COORDONNEES DU SOMMET J
               NSOM = LFACES(J,NOF)
               XYZ(1,J) = XYZSOM(1,NSOM)
               XYZ(2,J) = XYZSOM(2,NSOM)
               XYZ(3,J) = XYZSOM(3,NSOM)
 10         CONTINUE
C
C           TRACE DES ARETES DE LA FACE SANS REMPLISSAGE
            K = NBSF
            DO 30 J=1,NBSF
               CALL TRAIT3D( NCGRIM, XYZ(1,K), XYZ(1,J) )
               K = J
 30         CONTINUE
C
         ELSE IF( NOF .LT. 0 .AND. NOF .GE. -NBFACT ) THEN
C
C           FACE D'UN EF SECTIONNE PAR LE PLAN
C           ==================================
            NOF = -NOF
C
C           LE NOMBRE DE SOMMETS DE LA FACE
            IF( NOSTFA(4,NOF) .GT. 0 ) THEN
               NBSF = 4
            ELSE
               NBSF = 3
            ENDIF
C
            DO 110 J=1,NBSF
C              LES COORDONNEES DU SOMMET J
               NSOM = NOSTFA(J,NOF)
               XYZ(1,J) = XYZSOM(1,NSOM)
               XYZ(2,J) = XYZSOM(2,NSOM)
               XYZ(3,J) = XYZSOM(3,NSOM)
 110         CONTINUE
C
            IF( IAVNEF .NE. 0 .OR. REDUCF .GT. 0.0 ) THEN
C              REDUCTION DES FACES DEMANDEE
C              CALCUL DES COORDONNEES XYZP DU BARYCENTRE DE LA FACE
               CALL COBAPO( NBSF, XYZ, XYZP )
C              L'HOMOTHETIE DE CENTRE LE BARYCENTRE DE LA FACE
               DO 120 J=1,NBSF
                  XYZ(1,J) = XYZ(1,J) * REDUC1 + XYZP(1) * REDUCF
                  XYZ(2,J) = XYZ(2,J) * REDUC1 + XYZP(2) * REDUCF
                  XYZ(3,J) = XYZ(3,J) * REDUC1 + XYZP(3) * REDUCF
 120           CONTINUE
            ENDIF
C
C           COULEUR PAR DEFAUT DE LA FACE
            NCF = NCBLAN
C
            IF ( LCRITR .GT. 0 ) THEN
C
C              LA COULEUR SELON LA QUALITE DE L'EF 3D
C              --------------------------------------
C              LE NUMERO DE L'EF
               NEF = NUEFFA(NOF)
C
C              LE NUMERO DES SOMMETS DE L'EF NEF
               CALL NSEFNS( NEF, NUTYMA, NBSOEF, NBTGEF,
     %                      LDAPEF, LDNGEF, LDTGEF,
     %                      MNNSEF, NX, NY, NZ,
     %                      NCOGEL, NUGEEF, NUEFTG, NOSOEF, IERR )
C
C              CALCUL DE LA QUALITE DE L'ELEMENT FINI
               CALL QUALEF( NCOGEL,   NOSOEF, NBSOM, XYZSOM,
     %                      SURFVOLU, QUALIT, IERR )
               IF ( IERR .NE. 0 ) RETURN
C
C              LA COULEUR VISUALISE LA QUALITE
               NCF = NCOQUA( QUALIT )
C
            ELSE
C
C              COULEURS SELON DIRECTION DE VISEE ET ELOIGNEMENT
C              ------------------------------------------------
C              L'ELOIGNEMENT DE LA FACE NF
               R  = ( DIOEIL(NF)  - DMIN ) / DELTA
               R2 = ( SQRT(R+EP1) - EP2  ) / EP3
C
C              L'OTHOGONALITE A LA FACE
C              CNORFA LES COORDONNEES DE LA NORMALE A LA FACE (NORME=1)
               CALL NORF34( NBSF, NOSTFA(1,NOF), XYZSOM, CNORFA, IERR )
               IF( IERR .NE. 0 ) GOTO 130
C
C              LE PRODUIT SCALAIRE DIRECTION VISEE ET NORMALE A LA FACE
               R = PROSCR( DIREVI , CNORFA , 3 )
               R = 1. - ABS(R)
C
               IF( IAVELO .NE. 0 ) THEN
C                 LA COULEUR PONDEREE PAR LA DIRECTION ET L'ELOIGNEMENT
                  R =  R * POID + POID1 * (1. - R2)
               ENDIF
C
C              LA COULEUR DU TRACE
               NCF = NINT( NBCOUL * R + N1COUL )
            ENDIF
C
C           LE TRACE DE LA FACE EN ROUGE ET ARETES EN NOIR
 130        CALL FACE3D( NCF, NCNOIR, NBSF, XYZ )
C
C           TRACE EVENTUEL DU NO DES SOMMETS DE LA FACE
            IF( IAVNSO .NE. 0 ) THEN
               DO 140 J=1,NBSF
                  NSOM = NOSTFA(J,NOF)
                  WRITE( NMSOMM , '(I8)' ) NSOM
                  CALL SANSBL( NMSOMM, L )
                  CALL TEXTE3D( NCONSO, XYZSOM(1,NSOM), NMSOMM(1:L) )
 140           CONTINUE
            ENDIF
C
C           TRACE EVENTUEL DU NO DE L'EF
            IF( IAVNEF .NE. 0 ) THEN
               IF( REDUCF .EQ. 0.0 ) THEN
C                 CALCUL DES COORDONNEES XYZP DU BARYCENTRE DE LA FACE
                  CALL COBAPO( NBSF, XYZ, XYZP )
               ENDIF
               NEF = ABS( NUEFFA(NOF) )
               NMSOMM = 'EF '
               WRITE( NMSOMM(3:10) , '(I8)' ) NEF
               CALL SANSBL( NMSOMM, L )
               CALL TEXTE3D( NCONEF, XYZP, NMSOMM(1:L) )
            ENDIF
C
         ELSE
C
C           TRACE DE L'UNE DES  4 ARETES DU RECTANGLE
C           DE VISUALISATION DU PLAN
C           =========================================
            NS1 = -NOF - NBFACT
            IF( NS1 .GT. 1 ) THEN
               NS0 = NS1 - 1
            ELSE
               NS0 = 4
            ENDIF
            CALL XVEPAISSEUR( 3 )
            CALL TRAIT3D( NCROSE, XYZPLAN(1,NS0), XYZPLAN(1,NS1) )
            CALL XVEPAISSEUR( NEPARF )
C
         ENDIF
 500  CONTINUE
C
C     FINITION DU TRACE DES 3 PLANS PAR TRACE DES 3 ARETES VUES
      CALL T3PLAP
C
C     RETOUR AU TRACE CONTINU DES TRAITS
      CALL XVTYPETRAIT( LIGCON )
C
      RETURN
      END
