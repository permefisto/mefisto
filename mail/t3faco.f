      SUBROUTINE T3FACO( NMOBJT , NUOBJT ,
     %                   NBMOFA , NBFACE , LFACES , NBSOM, XYZSOM ,
     %                   NBFAFR , NOFAFR , DIOEIL , MNNSEF )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACER EN COULEURS LES FACES DU CONTOUR D'UN VOLUME
C -----    UNE FACE EST TRACEE SI ELLE APPARTIENT A UN SEUL ELEMENT FINI
C
C ENTREES:
C --------
C NMOBJT : NOM DE L'OBJET A TRACER
C NUOBJT : NUMERO DU VOLUME
C NBMOFA : NOMBRE DE MOTS PAR FACE
C NBFACE : NOMBRE DE FACES
C LFACES : TABLEAU ENTIER DU NO DES SOMMETS DES FACES DU VOLUME
C NBSOM  : NOMBRE DE SOMMETS DU MAILLAGE
C XYZSOM : COORDONNEES DES SOMMETS
C NBFAFR : NOMBRE DE FACES FRONTALIERES
C NOFAFR : NUMERO DES FACES SELON LEUR DISTANCE CROISSANTE A L'OEIL
C DIOEIL : DISTANCE A L'OEIL DES FACES
C MNNSEF : ADRESSE MCN DU TABLEAU NSEF DU VOLUME
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS         MAI 1994
C2345X7..............................................................012
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      include"./incl/gsmenu.inc"
C
C     LE NUMERO D'UNITE DU CLAVIER, FENETRE DES AFFICHAGES ET
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
C
C     DECLARATION DU SUPER-TABLEAU NUMERIQUE MCN
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
      CHARACTER*(*)     NMOBJT
      CHARACTER*8       NMSOMM
C
      INTEGER           LFACES(1:NBMOFA,1:NBFACE), NOFAFR(1:NBFAFR)
      REAL              XYZSOM(1:3,1:NBSOM), DIOEIL(1:NBFAFR)
      REAL              CNORFA(1:3), XYZP(1:3),
     %                  XYZ(3,4)
      INTEGER           NOSOEL(1:64)
C
      IERR = 0
C
C     LA DIRECTION DE VISEE PTV-OEIL
      DIREVI(1) = AXOEIL(1) - AXOPTV(1)
      DIREVI(2) = AXOEIL(2) - AXOPTV(2)
      DIREVI(3) = AXOEIL(3) - AXOPTV(3)
C
C     LA DIRECTION DE VISEE NORMALISEE A 1.
      CALL NORMER( 3, DIREVI, IERR )
      IF( IERR .NE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'OEIL ET POINT VU IDENTIQUES'
         CALL LEREUR
         RETURN
      ENDIF
C
C     CHOIX DE LA COULEUR SELON LE TERMINAL COULEUR OU NON
      NBCOUL = NDCOUL - N1COUL
C
C     PREPARATION EN CAS DE TRACE DE TRACE EN FONCTION DE LA QUALITE
C     LES PARAMETRES DES NO SOMMET DU MAILLAGE
      IF ( LCRITR .GT. 0 ) THEN
         CALL NSEFPA( MCN(MNNSEF),
     %                NUTYMA, NBSOEL, NBSOEF, NBTGEF,
     %                LDAPEF, LDNGEF, LDTGEF, NBEFOB,
     %                NX    , NY    , NZ    ,
     %                IERR   )
      ENDIF
C
C     CALCUL DES COULEURS DES FACES AVEC PRISE EN COMPTE
C     DE LA DIRECTION DE VISEE ET L'ELOIGNEMENT
C     POID : POIDS DE LA DIRECTION DE VISEE DANS CE CALCUL
C            ATTENTION: 0 < POID < 1
C     ----------------------------------------------------
C     DISTANCES MIN ET MAX
      DMIN = DIOEIL(1)
      DMAX = DIOEIL(NBFAFR)
C
C     CALCUL DES COULEURS DES FACES AVEC PRISE EN COMPTE
C     DE LA DIRECTION DE VISEE ET L'ELOIGNEMENT
C     POID : POIDS DE LA DIRECTION DE VISEE DANS CE CALCUL
C            ATTENTION: 0 < POID < 1
C     ----------------------------------------------------
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
C     LE TRACE DES FACES EN COMMENCANT PAR LES PLUS ELOIGNEES
C     =======================================================
      DO 200 NF = 1, NBFAFR
C
C        LE NUMERO DE LA FACE LA PLUS ELOIGNEE NON TRACEE
         NF1 = NOFAFR( NF )
C
C        LE NOMBRE DE SOMMETS DE LA FACE
         IF( LFACES(4,NF1) .GT. 0 ) THEN
            NAF = 4
         ELSE
            NAF = 3
         ENDIF
C
         DO 120 J=1,NAF
C           LES COORDONNEES DU SOMMET J
            NSOM = LFACES(J,NF1)
            XYZ(1,J) = XYZSOM(1,NSOM)
            XYZ(2,J) = XYZSOM(2,NSOM)
            XYZ(3,J) = XYZSOM(3,NSOM)
 120     CONTINUE
C
         IF( IAVNEF .NE. 0 .OR. REDUCF .GT. 0.0 ) THEN
C           REDUCTION DES FACES
C           CALCUL DES COORDONNEES XYZP DU BARYCENTRE DE LA FACE
            CALL COBAPO( NAF, XYZ, XYZP )
C           L'HOMOTHETIE DE CENTRE LE BARYCENTRE DE LA FACE
            DO 160 J=1,NAF
               XYZ(1,J) = XYZ(1,J) * REDUC1 + XYZP(1) * REDUCF
               XYZ(2,J) = XYZ(2,J) * REDUC1 + XYZP(2) * REDUCF
               XYZ(3,J) = XYZ(3,J) * REDUC1 + XYZP(3) * REDUCF
 160        CONTINUE
         ENDIF
C
         IF( LCRITR .EQ. 0 ) THEN
C
C           COULEURS SELON DIRECTION DE VISEE ET ELOIGNEMENT
C           ================================================
C           L'ELOIGNEMENT DE LA FACE
            R  = ( DIOEIL(NF)  - DMIN ) / DELTA
            R2 = ( SQRT(R+EP1) - EP2  ) / EP3
C
C           L'OTHOGONALITE A LA FACE
C           CNORFA LES COORDONNEES DE LA NORMALE A LA FACE (NORME=1)
            CALL NORF34( NAF , LFACES(1,NF1) , XYZSOM , CNORFA, IERR )
            IF( IERR .NE. 0 ) GOTO 200
C
C           LE PRODUIT SCALAIRE DIRECTION VISEE ET NORMALE A LA FACE
            R = PROSCR( DIREVI , CNORFA , 3 )
            R = 1. - ABS(R)
C
            IF( IAVELO .NE. 0 ) THEN
C              LA COULEUR PONDEREE PAR LA DIRECTION ET L'ELOIGNEMENT
               R =  R * POID + POID1 * (1. - R2)
            ENDIF
C
C           LA COULEUR DU TRACE
            J = NINT( NBCOUL * R + N1COUL )
C
         ELSEIF ( LCRITR .GT. 0 ) THEN
C
C           COULEURS SELON LES QUALITES DES ELEMENTS FINIS
C           ==============================================
C           ELEMENT AUQUEL APPARTIENT LA FACE
            NELT = ABS( LFACES(6,NF1) )
C
C           LE NUMERO DES NBSOEF SOMMETS DU CUBE NELT
            CALL NSEFNS( NELT   , NUTYMA , NBSOEF , NBTGEF,
     %                   LDAPEF , LDNGEF , LDTGEF ,
     %                   MNNSEF , NX , NY , NZ ,
     %                   NCOGEL , NUGEEF , NUEFTG, NOSOEL , IERR )
C           CALCUL DE LA QUALITE DE L'ELEMENT FINI
            CALL QUALEF( NCOGEL ,  NOSOEL , NBSOM, XYZSOM ,
     %                   SURFVOLU, QUALIT , IERR )
            IF ( IERR .NE. 0 ) RETURN
C           LA COULEUR VISUALISE LA QUALITE
            J = NCOQUA( QUALIT )
         ENDIF
C
C        LE TRACE DE LA FACE
         CALL FACE3D( J, NCOUAF, NAF, XYZ )
C
C        TRACE EVENTUEL DU NO DE L'EF
         IF( IAVNEF .NE. 0 ) THEN
            IF( REDUCF .EQ. 0.0 ) THEN
C              CALCUL DES COORDONNEES XYZP DU BARYCENTRE DE LA FACE
               CALL COBAPO( NAF, XYZ, XYZP )
            ENDIF
            NE = ABS( LFACES(6,NF1) )
            WRITE( NMSOMM , '(I8)' ) NE
            CALL SANSBL( NMSOMM, L )
            CALL TEXTE3D( NCONEF, XYZP, NMSOMM(1:L) )
         ENDIF
C
C        TRACE EVENTUEL DU NO DES SOMMETS
         IF( IAVNSO .NE. 0 ) THEN
            DO 180 J=1,NAF
               NSOM = LFACES(J,NF1)
               WRITE( NMSOMM , '(I8)' ) NSOM
               CALL SANSBL( NMSOMM, L )
               CALL TEXTE3D( NCONSO, XYZSOM(1,NSOM), NMSOMM(1:L) )
 180        CONTINUE
         ENDIF
 200  CONTINUE
C
C     LE TRACE DE LA POIGNEE DU VOLUME
C     CALCUL DES COORDONNEES XYZP DU BARYCENTRE DE LA DERNIERE FACE
      CALL COBAPO( NAF, XYZ, XYZP )
      CALL ITEMV3( XYZP, NMOBJT, NUOBJT )

      RETURN
      END
