      SUBROUTINE TRFASEOB2( NBCOOR, NBPOI,  XYZPOI,
     %                      NBSO,   NSPSOM, NBAR,  NSLARE, NBFA, NSSFAC,
     %                      MOARFR, MXARFR, LAREFR,
     %                      NOBARY, DIBARY,
     %                      NUMXPO, NUMXLI, NUMXSU, NBPPLS,
     %                      NBFASE, NOSTFA, QUALIEF,
     %                      NBFACT, XYZPLAN )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACER LES ARETES FRONTALIERES, FACES DES SURFACES, ARETES
C -----    DES LIGNES ET POINTS D'UN OBJET
C
C ENTREES:
C --------
C NBCOOR : NOMBRE DE COORDONNEES D'UN POINT OU NOEUD (3 ou 6)
C NBPOI  : NOMBRE DE POINTS STOCKES DANS XYZPOI
C XYZPOI : XYZ DES POINTS DES EF DU MAILLAGE
C NBSO   : NOMBRE DE SOMMETS=POINTS UTILISATEUR
C NSPSOM : NUMERO XYZPOI DU SOMMET ET NO DU POINT
C NBAR   : NOMBRE D' ARETES SUR UNE LIGNE   UTILISATEUR
C NSLARE : NUMERO XYZPOI DES 2 SOMMETS ET NO DE LA LIGNE
C NBFA   : NOMBRE DE FACES  SUR UNE SURFACE UTILISATEUR
C NSSFAC : NUMERO XYZPOI DES 3 OU 4 SOMMETS ET NO DE LA SURFACE
C MOARFR : NOMBRE DE MOTS PAR ARETE FRONTALIERE DU TABLEAU LAREFR
C MXARFR : NOMBRE DE FACES DU TABLEAU LAREFR
C L1ARFR : NUMERO DANS LAREFR DE LA PREMIERE ARETE FRONTALIERE
C LAREFR : TABLEAU NUMERO DES 2 SOMMETS ET LIEN
C          LAREFR(1,I)= NO DU 1-ER  SOMMET DE L'ARETE FRONTALIERE
C          LAREFR(2,I)= NO DU 2-EME SOMMET > 1-ER  SOMMET
C          LAREFR(3,I)= 0 OU NUMERO DE L'ARETE FRONTALIERE SUIVANTE
C NBBARY : NOMBRE DE BARYCENTRES TRIES
C NOBARY : NUMERO DE L'ITEM A TRACER SELON SON ORDRE
C          LES PLUS ELOIGNES D'ABORD (ALGORITHME DU PEINTRE)
C DIBARY : COTE AXONOMETRIQUE DES BARYCENTRES
C NUMXPO : NUMERO MAXIMAL DES POINTS   DE L'OBJET
C NUMXLI : NUMERO MAXIMAL DES LIGNES   DE L'OBJET
C NUMXSU : NUMERO MAXIMAL DES SURFACES DE L'OBJET
C NBFASE : NOMBRE DE FACES DES EF SECTIONNES PAR LE PLAN
C NOSTFA : NO DES SOMMETS DANS XYZPOI  DES FACES DES EF SECTIONNES
C QUALIEF: QUALITE DE L'EF 3D DE CHAQUE FACE (POUR DEFINIR LA COULEUR)
C NBFACT : NOMBRE DE FACES FRONTALIERES + DES FACES DES EF SECTIONNES
C XYZPLAN: 3 COORDONNEES DES 4 SOMMETS DU RECTANGLE DU PLAN DE SECTION
C
C MODIFIE:
C --------
C NBPPLS : TEMOIN DE TRACE D'UN POINT OU LIGNE OU SURFACE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Laboratoire JL LIONS UPMC Paris Decembre 2006
C ...................................................................012
      PARAMETER     (LIGCON=0, LIGTIR=1)
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      REAL           XYZPOI(1:NBCOOR,1:NBPOI)
      INTEGER        NSPSOM(2,NBSO)
      INTEGER        NSLARE(3,NBAR)
      INTEGER        NSSFAC(5,NBFA)
      INTEGER        LAREFR(1:MOARFR,1:MXARFR)
      INTEGER        NOBARY(1:NBFACT)
      REAL           DIBARY(1:NBFACT)
      INTEGER        NBPPLS(1:NUMXPO+NUMXLI+NUMXSU)
      INTEGER        NOSTFA(4,NBFASE)
      REAL           QUALIEF(NBFASE)
      REAL           XYZPLAN(3,4)
      REAL           XYZ(3,4), XYZP(3)
      CHARACTER*24   NOMOBJ
      CHARACTER*10   NMSOMM
C
C     LE NOMBRE-1 DE COULEURS DISPONIBLES
      NBCOUL = NDCOUL - N1COUL
C
C     LES ARETES SONT TRACEES EN CONTINU
      CALL XVTYPETRAIT( LIGCON )
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
      DMIN = DIBARY(1)
      DMAX = DIBARY(NBFACT)
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
      CALL T3PLAV
      CALL TRAXES
      CALL XVEPAISSEUR( 1 )
C     AUCUN ITEM SUR L'ECRAN
      CALL ITEMS0
C
      NBSAF = NBSO + NBAR + NBFA
      DO 200 K = 1, NBFACT
C
C        LE NUMERO DE L'ITEM LE PLUS ELOIGNE ET RESTANT A TRACER
         N = NOBARY( K )
C
         IF( N .GE. 0 ) THEN
            IF( N .LE. NBSO ) THEN
C
C              TRACE DU SOMMET=POINT
C              ----------------------
C              LE NUMERO DU POINT-UTILISATEUR
               NUOB = NSPSOM(2,N)
C              LE NOM DU POINT
               CALL NMOBNU( 'POINT',  NUOB, NOMOBJ )
               CALL ITEMP3( XYZPOI(1,NSPSOM(1,N)), NOMOBJ, NUOB )
C
            ELSE IF( N .LE. NBSO+NBAR ) THEN
C
C              TRACE DE L'ARETE D'UNE LIGNE-UTILISATEUR
C              ----------------------------------------
               CALL XVEPAISSEUR( 3 )
               N = N - NBSO
C              LE NUMERO DE LA LIGNE
               NUOB = NSLARE(3,N)
C              LA COULEUR DE LA LIGNE
CCC             NCOUL = MOD( NUOB+NBSO, NBCOUL ) + N1COUL
               NCOUL = NCNOIR
C              NO DES 2 SOMMETS EXTREMITES
               NSL1 = NSLARE(1,N)
               NSL2 = NSLARE(2,N)
               IF( PREDUA .GT. 0 ) THEN
C                 REDUCTION DE L'ARETE
                  REDUCF = PREDUA * 0.01
                  REDUC1 = 1.0 - REDUCF
                  XYZ(1,1)=REDUCF*XYZPOI(1,NSL1) + REDUC1*XYZPOI(1,NSL2)
                  XYZ(2,1)=REDUCF*XYZPOI(2,NSL1) + REDUC1*XYZPOI(2,NSL2)
                  XYZ(3,1)=REDUCF*XYZPOI(3,NSL1) + REDUC1*XYZPOI(3,NSL2)
                  XYZ(1,2)=REDUC1*XYZPOI(1,NSL1) + REDUCF*XYZPOI(1,NSL2)
                  XYZ(2,2)=REDUC1*XYZPOI(2,NSL1) + REDUCF*XYZPOI(2,NSL2)
                  XYZ(3,2)=REDUC1*XYZPOI(3,NSL1) + REDUCF*XYZPOI(3,NSL2)
C                 LE TRACE DE L'ARETE REDUITE
                  CALL TRAIT3D( NCOUL, XYZ(1,1), XYZ(1,2) )
               ELSE
C                 LE TRACE DE L'ARETE
                  CALL TRAIT3D( NCOUL, XYZPOI(1,NSL1), XYZPOI(1,NSL2) )
               ENDIF
C
               IF( NBPPLS(NUMXPO+NUOB) .EQ. 0 ) THEN
C                 LE NOM DE LA LIGNE EST TRACE
                  CALL XVCOULEUR( NCGRIS )
                  CALL NMOBNU( 'LIGNE',  NUOB, NOMOBJ )
                  CALL ITEML3( XYZPOI(1,NSL2), NOMOBJ, NUOB )
                  NBPPLS(NUMXPO+NUOB) = 1
               ENDIF
               CALL XVEPAISSEUR( 1 )
C
            ELSE IF( N .LE. NBSAF ) THEN
C
C              TRACE DE LA FACE D'UNE SURFACE
C              ------------------------------
C              COULEUR SELON LE NUMERO DE LA SURFACE
               N = N - NBSO - NBAR
C              LE NUMERO DE LA SURFACE
               NUOB = NSSFAC(5,N)
C              LA COULEUR DE LA SURFACE
cccc              AVEC CALL PALCDE(1)
ccc               NCOUL = MOD( NUOB+NBSO+NBAR, NBCOUL ) + N1COUL
C              LA COULEUR DE LA SURFACE ENTRE 1 ET 8 =>
C              NCBLAN, NCROUG, NCJAUN, NCVERT, NCCYAN, NCBLEU, NCMAGE, NCGRIS
               NCOUL = MOD( NUOB+NBSO+NBAR,8 ) + 1
C
C              LE NOMBRE DE SOMMETS DE LA FACE A TRACER
               IF( NSSFAC(4,N) .EQ. 0 ) THEN
                  NBSF = 3
               ELSE
                  NBSF = 4
               ENDIF
C
C              LE TRACE DE LA FACE
C              LES XYZ DES SOMMETS DE LA FACE
               DO 180 I=1,NBSF
                  DO 170 J=1,3
                     XYZ(J,I) = XYZPOI(J,NSSFAC(I,N))
 170              CONTINUE
 180           CONTINUE
C
ccc               IF( PREDUF .GT. 0 ) THEN
C              REDUCTION DE LA FACE POUR LA DISSOCIER DES FACES DE SECTION
               REDUCF = 0.8
               REDUC1 = 1.0 - REDUCF
C              CALCUL DES COORDONNEES XYZP DU BARYCENTRE DE LA FACE
               CALL COBAPO( NBSF, XYZ, XYZP )
C              L'HOMOTHETIE DE CENTRE LE BARYCENTRE DE LA FACE
               DO J=1,NBSF
                  XYZ(1,J) = XYZ(1,J) * REDUC1 + XYZP(1) * REDUCF
                  XYZ(2,J) = XYZ(2,J) * REDUC1 + XYZP(2) * REDUCF
                  XYZ(3,J) = XYZ(3,J) * REDUC1 + XYZP(3) * REDUCF
               ENDDO
ccc               ENDIF
ccc               CALL FACE3D( NCOUL, NCNOIR, NBSF, XYZ )
               CALL FACE3D( NCOUL, NCBLAN, NBSF, XYZ )
C
               IF( NBPPLS(NUMXPO+NUMXLI+NUOB) .EQ. 0 ) THEN
C                 LE NOM DE LA SURFACE EST TRACE
                  CALL NMOBNU( 'SURFACE',  NUOB, NOMOBJ )
                  CALL XVCOULEUR( NCGRIS )
                  CALL ITEMS3( XYZPOI(1,NSSFAC(2,N)), NOMOBJ, NUOB )
                  NBPPLS(NUMXPO+NUMXLI+NUOB) = 1
               ENDIF
C
            ELSE IF( N .GT. NBSAF ) THEN
C
C              TRACE DE L' ARETE FRONTALIERE DES VOLUMES DE L'OBJET
C              ----------------------------------------------------
C              LE NUMERO DE L'ARETE DANS LE TABLEAU LAREFR
               N = N - NBSAF
               CALL TRAIT3D( NCGRIM, XYZPOI(1,LAREFR(1,N)),
     %                               XYZPOI(1,LAREFR(2,N)) )
            ENDIF
C
          ELSE
C
C           ICI N<0
            IF( N .GE. -NBFASE ) THEN
C
C              FACE D'UN EF SECTIONNE PAR LE PLAN
C              ----------------------------------
               N = -N
C
C              LE NOMBRE DE SOMMETS DE LA FACE
               IF( NOSTFA(4,N) .GT. 0 ) THEN
                  NBSF = 4
               ELSE
                  NBSF = 3
               ENDIF
C
               DO 110 J=1,NBSF
C                 LES COORDONNEES DU SOMMET J
                  NSOM = NOSTFA(J,N)
                  XYZ(1,J) = XYZPOI(1,NSOM)
                  XYZ(2,J) = XYZPOI(2,NSOM)
                  XYZ(3,J) = XYZPOI(3,NSOM)
 110           CONTINUE
C
               IF( PREDUF .GT. 0 ) THEN
C                 UNE REDUCTION DEMANDEE EST UTILISEE
                  REDUCF = PREDUF * 0.01
               ELSE
C                 UNE REDUCTION DE 10% EST IMPOSEE POUR VOIR DEDANS!
                  REDUCF = 0.1
               ENDIF
               IF( IAVNEF .NE. 0 .OR. REDUCF .GT. 0.0 ) THEN
C                 REDUCTION DES FACES DEMANDEE
C                 CALCUL DES COORDONNEES XYZP DU BARYCENTRE DE LA FACE
                  CALL COBAPO( NBSF, XYZ, XYZP )
C                 L'HOMOTHETIE DE CENTRE LE BARYCENTRE DE LA FACE
                  REDUC1 = 1.0 - REDUCF
                  DO 120 J=1,NBSF
                     XYZ(1,J) = XYZ(1,J) * REDUC1 + XYZP(1) * REDUCF
                     XYZ(2,J) = XYZ(2,J) * REDUC1 + XYZP(2) * REDUCF
                     XYZ(3,J) = XYZ(3,J) * REDUC1 + XYZP(3) * REDUCF
 120              CONTINUE
               ENDIF
C
C              COULEUR PAR DEFAUT DE LA FACE
               NCF = NCBLAN
C
               IF ( LCRITR .GT. 0 ) THEN
C
C                 LA COULEUR VISUALISE LA QUALITE DE L'EF 3D
                  NCF = NCOQUA( QUALIEF(N) )
C
               ELSE
C
C                 COULEURS SELON DIRECTION DE VISEE ET ELOIGNEMENT
C                 L'ELOIGNEMENT DE LA FACE K
                  R  = ( DIBARY(K)  - DMIN ) / DELTA
                  R2 = ( SQRT(R+EP1) - EP2  ) / EP3
C
C                 L'ORTHOGONALITE A LA FACE
C                 CNORFA LES COORDONNEES DE LA NORMALE A LA FACE (NORME=1)
                  CALL NORF34( NBSF, NOSTFA(1,N), XYZPOI, CNORFA, IERR )
                  IF( IERR .NE. 0 ) GOTO 130
C
C                 LE PRODUIT SCALAIRE DIRECTION VISEE ET NORMALE A LA FACE
                  R = PROSCR( DIREVI , CNORFA , 3 )
                  R = 1. - ABS(R)
C
                  IF( IAVELO .NE. 0 ) THEN
C                    LA COULEUR PONDEREE PAR LA DIRECTION ET L'ELOIGNEMENT
                     R =  R * POID + POID1 * (1. - R2)
                  ENDIF
C
C                 LA COULEUR DU TRACE
                  NCF = NINT( NBCOUL * R + N1COUL )
               ENDIF
C
C              LE TRACE DE LA FACE EN COULEUR NCF ET ARETES EN NOIR
 130           CALL FACE3D( NCF, NCNOIR, NBSF, XYZ )
C
C              TRACE EVENTUEL DU NO DES SOMMETS DE LA FACE
               IF( IAVNSO .NE. 0 ) THEN
                  DO 140 J=1,NBSF
                     NSOM = NOSTFA(J,N)
                     WRITE( NMSOMM , '(I8)' ) NSOM
                     CALL SANSBL( NMSOMM, L )
                     CALL TEXTE3D( NCONSO, XYZPOI(1,NSOM), NMSOMM(1:L) )
 140              CONTINUE
               ENDIF
C
            ELSE
C
C              TRACE DE L'UNE DES  4 ARETES DU RECTANGLE
C              DE VISUALISATION DU PLAN DE SECTION DES EF 3D
C              ---------------------------------------------
               NS1 = -N - NBFASE
               IF( NS1 .GT. 1 ) THEN
                  NS0 = NS1 - 1
               ELSE
                  NS0 = 4
               ENDIF
C              L'ARETE EST TRACEE EN CONTINU
               CALL XVEPAISSEUR( 3 )
               CALL TRAIT3D( NCROSE, XYZPLAN(1,NS0), XYZPLAN(1,NS1) )
               CALL XVEPAISSEUR( NEPARF )
C
            ENDIF
         ENDIF
 200  CONTINUE
C
C     FINITION DU TRACE DES 3 PLANS PAR TRACE DES 3 ARETES VUES
      CALL T3PLAP
C
C     RETOUR AU TRACE CONTINU DES TRAITS
      CALL XVTYPETRAIT( LIGCON )
      END
