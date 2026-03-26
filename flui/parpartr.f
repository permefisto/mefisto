      SUBROUTINE PARPARTR( KNOMOB, MNNPEF, TIMES,  NCAS0, NCAS1,
     %                     NBSFTR, NUSFTR, MNXYZSF, MNSEFSF,
     %                     NBPART, XVRVIPART, RAYPARMI,RAYPARMX,
     %                     NBXYZPART, xyzvpart, VITMXPAR,
     %                     NBPARINA, NUPARINA,
     %                     NBITEMS, NOITEMS, DISTEMS, NODRITEM )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  TRACER 3D DU PARCOURS DES PARTICULES DANS UN FLUIDE 3D MAILLE EN
C -----  TETRAEDRES DURANT L'INTERVALLE DE TEMPS OU LA VITESSE EST CONNUE
C         sur des INTERPOLATIONS de TAYLOR-HOOD ou BREZZI-FORTIN

C ENTREES:
C --------
C KNOMOB : NOM DE L'OBJET
C MNNPEF : ADRESSE MCN DES TABLEAUX ELEMENTS DE CETTE TOPOLOGIE
C TIMES  : (NCAST0:NCAST1) TEMPS DU CALCUL DES NBVPFILE VECTEURS
C NCAS0  : NUMERO DU PREMIER VECTEUR VITESSE+PRESSION A TRAITER
C NCAS1  : NUMERO DU DERNIER VECTEUR VITESSE+PRESSION A TRAITER
C          STOCKE SUR DES FICHIERS DANS LE REPERTOIRE DU PROJET

C NBSFTR  : NOMBRE DE SURFACES TRIANGULEES A TRACER
C NUSFTR  : NUMERO DES NBSFTR SURFACES A TRACER
C MNXYZSF : ADRESSE MCN DU TMS XYZSOMMET DES NBSFTR SURFACES A TRACER
C MNSEFSF : ADRESSE MCN DU TMS NSEF      DES NBSFTR SURFACES A TRACER

C NBPART  : NOMBRE DE PARTICULES DE PARCOURS A TRACER
C XVRVIPART:XYZ + VITESSE INITIALE + RAYON de la BOULE des NBPART PARTICULES
C RAYPARMI: RAYON MINIMAL de l'UNE des BOULES-PARTICULES
C RAYPARMX: RAYON MAXIMAL de l'UNE des BOULES-PARTICULES

C NBXYZPART:NOMBRE DE POINTS et V STOCKES DU PARCOURS DE CHAQUE PARTICULE
C           <0 MARQUAGE DE PARCOURS INACHEVE DE LA PARTICULE
c xyzvpart: XYZ et VXYZ EN CHAQUE POINT D'INTEGRATION DU PARCOURS
C           DES PARTICULES xyzvpart( 1:6, 1:ABS(NBXYZPART(K)), 1:NBPART )
C VITMXPAR: MODULE MAXIMAL DE LA VITESSE DES PARTICULES LORS DES PARCOURS

C NBPARINA: NOMBRE DE PARTICULES DE CALCUL DU PARCOURS INACHEVE
C NUPARINA: NUMERO DES NBPARINA PARTICULES DE PARCOURS INACHEVE

C NBITEMS : NOMBRE DES ITEMS (SEGMENT ou TRIANGLE) A TRACER
C NOITEMS : NUMERO DES ITEMS APRES LE TRI SELON LA DISTANCE A L'OEIL
C DISTEMS : DISTANCE A L'OEIL DES ITEMS
C NBPARSF : NOMBRE SOMME DU NOMBRE DE PARTICULES ET SURFACES A TRACER
C NODRITEM: NO DU DERNIER ITEM DES SEGMENTS DES PARTICULES PUIS
C           DES TRIANGLES DES SURFACES A TRACER
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  Saint PIERRE du PERRAY         Decembre 2020
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/a_objet__topologie.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/a___npef.inc"
      include"./incl/a___vecteur.inc"
      include"./incl/a___xyzpoint.inc"
      include"./incl/a___xyznoeud.inc"
      include"./incl/a___face.inc"
      include"./incl/ponoel.inc"
      include"./incl/donele.inc"
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      include"./incl/ctemps.inc"
      include"./incl/xyzext.inc"
      include"./incl/nctyef.inc"

      include"./incl/pp.inc"
      COMMON          MCN(MOTMCN)
      REAL            RMCN(1)
      EQUIVALENCE    (MCN(1),RMCN(1))

      CHARACTER*(*)    KNOMOB
      CHARACTER*24     KNOMFGIF, KNMSURF
      CHARACTER*8      KNOPART

      INTEGER          NUSFTR(NBSFTR), MNXYZSF(NBSFTR), MNSEFSF(NBSFTR),
     %                 NBXYZPART(NBPART), NUPARINA(NBPARINA),
     %                 NOITEMS(NBITEMS), MNXYZS(3),
     %                 NODRITEM(0:NBPART+NBSFTR), NOSOTR(3)
      REAL             TIMES(NCAS0:NCAS1),
     %                 XVRVIPART(8,NBPART), RAYPARMI, RAYPARMX,
     %                 XYZPT0(3), XYZPT1(3), VPART0(3), VPART1(3),
     %                 DISTEMS(NBITEMS), BARYCENT(3), CNORFA(3),
     %                 XYZTRIA(3,3), COUL(3), R, R2, EP1, EP2, EP3,
     %                 DISTMIN, DISTMAX, DISTMIMX, POID, POID1

      type typ_r4tab2i
         REAL, dimension(:,:), pointer :: r4tab2i
      end type typ_r4tab2i
      type( typ_r4tab2i ), dimension(1:NBPART) :: xyzvpart

      DOUBLE PRECISION DTXYZ(4), DOUI
      INTRINSIC        SQRT, REAL, INT

      PRINT*
      IF( LANGAG .EQ. 0 ) THEN
         PRINT*,'parpartr: A partir du temps',TIMES(NCAS0),
     %          ' trace du PARCOURS des',NBPART,' PARTICULES dont',
     %           NBPARINA,' PARCOURS de CALCUL INACHEVE'
         IF( NBPARINA .GT. 0 ) THEN
           PRINT*,'parpartr: Liste des Particules de PARCOURS INACHEVE:'
           PRINT*,(NUPARINA(K),K=1,NBPARINA)
         ENDIF
      ELSE
         PRINT*,'parpartr: From the time',TIMES(NCAS0),
     %          ' TRIP of',NBPART,' PARTICLES with',
     %           NBPARINA,' COMPUTED UNFINISHED RUNS'
         IF( NBPARINA .GT. 0 ) THEN
          PRINT*,'parpartr: List of UNFINISHED RUN PARTICLE NUMBER:'
          PRINT*,(NUPARINA(K),K=1,NBPARINA)
         ENDIF
      ENDIF

      IF( LANGAG .EQ. 0 ) THEN
         PRINT*,'parpartr:',NBITEMS,' ITEMS a TRACER'
      ELSE
         PRINT*,'parpartr:',NBITEMS,' ITEMS to DRAW'
      ENDIF

      DO K = 1, NBSFTR
C        LE NOM DE LA SURFACE K A TRACER EST AFFICHE
         CALL NMOBNU( 'SURFACE', NUSFTR(K), KNMSURF )
         IF( LANGAG .EQ. 0 ) THEN
            PRINT*,'parpartr: TRACE de la SURFACE: ',KNMSURF
         ELSE
            PRINT*,'parpartr: DRAWING of SURFACE: ',KNMSURF
         ENDIF
      ENDDO

      IF( NCAS0 .GE. NCAS1 .OR. NBPART .LE. 0 )RETURN

C     NOM DU FICHIER VIDEO // 'path'
      CALL VIDEONM( 5, 'path', KNOMFGIF )

C     LES DONNEES SUR LES TETRAEDRES DU MAILLAGE DE L'OBJET
C     -----------------------------------------------------
C     L'ADRESSE MCN DU DEBUT DU TABLEAU NPEF DE CE TYPE D'EF
      MNELE = MCN( MNNPEF )

C     LE NUMERO DU TYPE DES ELEMENTS FINIS
      NUTYEL = MCN( MNELE + WUTYEL )

C     LE NOMBRE DE TELS ELEMENTS FINIS ici des TETRAEDRES
      NBELEM = MCN( MNELE + WBELEM )

C     LES CARACTERISTIQUES DE L'ELEMENT FINI TETRAEDRE
      CALL ELTYCA( NUTYEL )

C     NO D'INTERPOLATION DES COMPOSANTES DE LA VITESSE
C     NOMBRE DE NOEUDS VITESSE DE L'EF
      IF( NUTYEL .EQ. 19 ) THEN
C        TETRAEDRE BREZZI-FORTIN
         NOINTE = 3003
         NBNOE  = 5
         NBNOED = 4
      ELSE IF( NUTYEL .EQ. 20 ) THEN
C        TETRAEDRE TAYLOR-HOOD
         NOINTE = 3002
         NBNOE  = 10
         NBNOED = NBNOE
      ELSE
C        EF NON TRAITE ICI
         IERR = 1
         GOTO 9999
      ENDIF

C     NOMBRE DE COULEURS DISPONIBLES AU DELA DES COULEURS PRIMAIRES
      NBCOUL  = NDCOUL - NDCORE
      NBCOUL2 = NBCOUL / 2
      NDCOUL2 = NDCORE + NBCOUL2
      N1COULGRIS = NDCOUL2 + 1
C     NOMBRE DE GRIS - 1 CAR SOMME A N1COULGRIS
      NBCOULGRIS = NDCOUL - N1COULGRIS

C     LA PALETTE 17: ARC EN CIEL POUR LA VITESSE + GRISE POUR LES SURFACES
      CALL PALCDE(17)

C     NOMBRE MINIMAL DE SEGMENTS POUR TRACER LE PARCOURS D'UNE PARTICULE
      MINSEGATR = 1

C     CALCUL DES COULEURS DES FACES AVEC PRISE EN COMPTE
C     DE LA DIRECTION DE VISEE ET L'ELOIGNEMENT
      EP1   = 0.05
      EP2   = SQRT(EP1)
      EP3   = SQRT(1+EP1) - EP2

ccc      LCRITR = 0
      LCRITR = -1
      AXOAVA = 0
      AXOARR = 0

      TEMPS = TIMES( NCAS0 )

C     OPTIONS DE LA VISEE 3D POUR VOIR LE PARCOURS DES PARTICULES
C     ===========================================================
 10   CALL VISE3D( NMTCL )
      IF( NMTCL .LT. 0 ) GOTO 9999

C     INITIALISATION DE TRANSLATION ORBITE ZOOM
      IF( LORBITE .NE. 0 ) THEN
         CALL ORBITE0( NOTYEV )
         IF( NOTYEV .EQ. 0 ) GOTO 10
      ENDIF

C     L'ECRAN PIXELS EST EFFACE (PAS DE SCINTILLEMENT)
      CALL EFFACEMEMPX

C     la FONCTION 'REGION(t,x,y,z)' a t elle ete DONNEE
C     PAR l'UTILISATEUR?
C     REGION(t,x,y,z)=1 SI LE PARCOURS ARRIVANT en (t,x,y,z)
C                          EST DEMANDE
C                    =0 SINON
      CALL LXNMNO( NTFONC, 'REGION', NOFORE, I )
C     NOFORE>0 SI CETTE FONCTION REGION(t,x,y,z) EXISTE

C     INITIALISATION DE LA DIRECTION DE VISEE PTV-OEIL NORMALISEE A 1
C     ---------------------------------------------------------------
C     DISTANCE DE L'OEIL AU POINT VISE
 30   DPTVOEIL = 0
      DO K=1,3
         DIREVI(K) = AXOEIL(K) - AXOPTV(K)
         DPTVOEIL  = DPTVOEIL + DIREVI(K)**2
      ENDDO
      DPTVOEIL = SQRT( DPTVOEIL )

C     AXONOMETRIE DE L'OEIL POUR CALCULER LA DISTANCE A L'OEIL
      CALL XYZAXO( AXOEIL, AXOOEIL )
ccc      PRINT*,'parpartr: AXOOEIL=',AXOOEIL,'=DPTVOEIL=',DPTVOEIL

      IF( DPTVOEIL .NE. 0 ) THEN
         DO K=1,3
            DIREVI(K) = DIREVI(K) / DPTVOEIL
         ENDDO
      ELSE
         NBLGRC(NRERR) = 1
         KERR(1) = 'OEIL et POINT VU sont IDENTIQUES. VISEE IMPOSSIBLE'
         CALL LEREUR
         GOTO 10
      ENDIF

C     CONSTRUCTION DE LA DISTANCE A L'OEIL DES ITEMS APRES AXONOMETRIE
C     ================================================================
      NODRITEM(0) = 0
      NOITEM = 0
      DO K = 1, NBPART

C        CALCUL DU BARYCENTRE DES NBSTSEG-1 SEGMENTS
C        ENTRE 2 POINTS DU PARCOURS DE LA PARTICULE K
C        --------------------------------------------
         NBSTSEG = ABS( NBXYZPART( K ) )
         IF( NBSTSEG .LE. MINSEGATR ) THEN
C           PAS DE SEGMENT A TRACER POUR LE PARCOURS DE LA PARTICULE K
            GOTO 40
         ENDIF

         DO N = 2, NBSTSEG

C           XYZ DU BARYCENTRE DU SEGMENT N-1 DU PARCOURS DE LA PARTICULE K
C           --------------------------------------------------------------
            DO M=1,3
               BARYCENT(M) = ( xyzvpart(k)%r4tab2i(m,n-1)
     %                       + xyzvpart(k)%r4tab2i(m,n  ) ) / 2
            ENDDO

C           ELOIGNEMENT DU POINT A L'OEIL PAR PASSAGE EN AXONOMETRIE
            CALL XYZAXO( BARYCENT, BARYCENT )

C           STOCKAGE POUR LE TRI
            NOITEM = NOITEM + 1
            NOITEMS( NOITEM ) = NOITEM
            DISTEMS( NOITEM ) = -BARYCENT(3)

         ENDDO

C        NO DU DERNIER ITEM DU PARCOURS DE LA PARTICULE K
 40      NODRITEM(K) = NOITEM

      ENDDO

C     LES TRIANGULATIONS DES NBSFTR SURFACES A TRACER
C     -----------------------------------------------
      DISFMIN = 1E28
      DISFMAX =-1E28
      ZMIN =  1E28
      ZMAX = -1E28

      DO K = 1, NBSFTR

C        NOMBRE DE TRIANGLES DE LA SURFACE K
         NBTRSF = MCN( MNSEFSF(K) + WBEFOB )

         MNEF = MNSEFSF(K) + WUSOEF -5
         DO NTR = 1, NBTRSF

C           XYZ DU BARYCENTRE DU TRIANGLE N DE LA SURFACE K
C           -----------------------------------------------
            DO M=1,3
               BARYCENT(M) = 0
            ENDDO
            DO I=1,3
C              NO DU SOMMET I DU TRIANGLE NTR DE LA SURFACE K
               NS = MCN( MNEF + 4*NTR + I )
C              ADRESSE MCN-1 DE XYZ DU SOMMET I DU TRIANGLE NTR DE LA SURFACE K
               MNXYZS(I) = MNXYZSF(K) + WYZSOM-4 + 3*NS
            ENDDO
            DO M=1,3
               BARYCENT(M) = ( RMCN( MNXYZS(1)+M ) +
     %                         RMCN( MNXYZS(2)+M ) +
     %                         RMCN( MNXYZS(3)+M ) ) / 3
            ENDDO

            ZN = BARYCENT(3)
            IF( ZN .LT. ZMIN ) ZMIN = ZN
            IF( ZN .GT. ZMAX ) ZMAX = ZN

C           ELOIGNEMENT DU POINT A L'OEIL PAR PASSAGE EN AXONOMETRIE
            CALL XYZAXO( BARYCENT, BARYCENT )

            IF( BARYCENT(3) .GT. DPTVOEIL ) THEN
C              LE POINT EST DERRIERE L'OEIL => DISTANCE<0
               R = -R
            ENDIF

C           STOCKAGE DE L'ITEM NOITEM POUR LE TRI
            NOITEM = NOITEM + 1
            NOITEMS( NOITEM ) = NOITEM
            R = -BARYCENT(3)
            DISTEMS( NOITEM ) = R
            IF( R .LT. DISFMIN ) DISFMIN = R
            IF( R .GT. DISFMAX ) DISFMAX = R

         ENDDO

C        NO DU DERNIER ITEM DES TRIANGLES DE LA SURFACE K A TRACER
         NODRITEM(NBPART+K) = NOITEM

      ENDDO
      DISFMIMX = DISFMAX - DISFMIN
      ZMINMAX  = ZMAX - ZMIN

C     LE TRI CROISSANT DES DISTANCES (COTE -Z AXONOMETRIQUE) A L'OEIL
C     ===============================================================
      CALL TRITRP( NOITEM, DISTEMS, NOITEMS )
C     NOITEMS(N): NUMERO INITIAL DANS LES ITEMS APRES LE TRI
 
C     DISTANCES MIN ET MAX des BARYCENTRES des ITEMS a l'OEIL
      DISTMIN =  DISTEMS( 1 )
      ZZZMIN  =  DISTMIN
      DISTMAX =  DISTEMS( NOITEM )
      ZZZMAX  =  DISTMAX
      DISTMIMX = DISTMAX - DISTMIN

C     POID : POIDS DE LA DIRECTION DE VISEE DANS CE CALCUL
C            ATTENTION: 0 < POID < 1
C     ----------------------------------------------------
      IF( DISTMIMX .EQ. 0 ) THEN
         DISTMIMX = 1.0
         POID     = 1.0
      ELSE
         POID = 0.6
      ENDIF
      POID1 = 1. - POID

C     TRACE DES AXES 3D
C     -----------------
      CALL TRAXE3

C     TRACE DES 12 ARETES DE L'HEXAEDRE ENGLOBANT
C     -------------------------------------------
      CALL TRHEXAGL

C     TRACE DES ITEMS SELON L'ALGORITHME DU PEINTRE
C     =============================================
      DO 500 NIT = NOITEM, 1, -1

C        LE NUMERO DE L'ITEM AVANT LE TRI
         NITEM = NOITEMS( NIT )

         IF( NITEM .LE. NODRITEM(NBPART) ) THEN

C           L'ITEM EST UN SEGMENT DE L'UN DES NBPART PARCOURS DES PARTICULES
C           ----------------------------------------------------------------
C           QUELLE PARTICULE?
            DO K=1,NBPART

               IF( NITEM .LE. NODRITEM(K) ) THEN

C                 LE SEGMENT NITEM APPARTIENT AU PARCOURS DE LA PARTICULE K
C                 NO DU SEGMENT DU PARCOURS DE LA PARTICULE K
                  NOSEGMK = NITEM - NODRITEM(K-1)

C                 TRACE DU SEGMENT NOSEGMK SI DANS LA REGION
                  IF( NOFORE .GT. 0 ) THEN
C                    la FONCTION REGION(t,x,y,z) EXISTE
C                    LE POINT FINAL DU PARCOURS DE LA PARTICULE K
                     NBSTSEG = NODRITEM(K) - NODRITEM(K-1)
                     DTXYZ(1) = TEMPS
                     DO M=1,3
                        DTXYZ(1+M) = xyzvpart(K)%r4tab2i( M, NBSTSEG+1 )
                     ENDDO
C                    REGION(t,x,y,z)=1 SI LE PARCOURS d'une PARTICULE
C                                      ARRIVANT en (t,x,y,z) DOIT ETRE TRACE
C                                   =0 SINON
                     CALL FONVAL( NOFORE, 4, DTXYZ, NCODEV, DOUI )
                     IF( NINT( DOUI ) .NE. 1 ) GOTO 500
                  ENDIF

C                 XYZ et VXYZ DEBUT DU SEGMENT NOSEGMK DE LA PARTICULE K
C                 XYZ et VXYZ FIN   DU SEGMENT NOSEGMK DE LA PARTICULE K
                  DO M=1,3
                     XYZPT0(M) = xyzvpart(K)%r4tab2i( M,   NOSEGMK   )
                     XYZPT1(M) = xyzvpart(K)%r4tab2i( M,   NOSEGMK+1 )
                     VPART0(M) = xyzvpart(K)%r4tab2i( 3+M, NOSEGMK   )
                     VPART1(M) = xyzvpart(K)%r4tab2i( 3+M, NOSEGMK+1 )
                  ENDDO

C                 TRACE DE LA PARTICULE K AU POINT DE DEPART ET AFFICHAGE
                  IF( NOSEGMK .EQ. 1 ) THEN
                     KNOPART(1:2) = '+D'
                     WRITE(KNOPART(3:8),'(I6)') K
                     CALL SANSBL( KNOPART, NBK )
C                    '+Depart' + NO PARTICULE
                     CALL SYMBOLE3D( NCNOIR, XYZPT0, KNOPART(1:NBK) )
                  ENDIF

C                 NORME DE LA VITESSE VPART0 DE LA PARTICULE EN XYZPT0
                  VPNORM = SQRT(VPART0(1)**2+VPART0(2)**2+VPART0(3)**2)

C                 TRACE DU SEGMENT XYZPT0-XYZPT1 AVEC UNE COULEUR
C                 SELON LA VITESSE VPART0 DE LA PARTICULE EN XYZPT0
                  NCOUL = NINT( N1COUL + VPNORM / VITMXPAR * NBCOUL2 )
                  IF( NCOUL .LT. N1COUL ) THEN
                     NCOUL = N1COUL
                  ENDIF
                  IF( NCOUL .GT. NDCOUL2 ) THEN
                     NCOUL = NDCOUL2
                  ENDIF

C                 EPAISSEUR DU TRAIT DU SEGMENT EN FONCTION DU RAYON
C                 DE LA PARTICULE
C                 RAYPARMI: RAYON MINIMAL de l'UNE des BOULES-PARTICULES
C                 RAYPARMX: RAYON MAXIMAL de l'UNE des BOULES-PARTICULES
C                 RAYPARTK: RAYON DE LA BOULE PARTICULE K
                  RAYPARTK = XVRVIPART(7,K)
                  NBEPAIS  = NINT( 3 * (RAYPARTK-RAYPARMI)
     %                               / (RAYPARMX-RAYPARMI) )

C                 LE SEGMENT NOSEGMK DU PARCOURS DE LA PARTICULE K
C                 EST TRACE AVEC NBEPAIS EPAISSEURS
                  CALL XVEPAISSEUR( NBEPAIS )

                  CALL TRAIT3D( NCOUL, XYZPT0, XYZPT1 )

                  IF( NITEM .EQ. NODRITEM(K) ) THEN

C                    DERNIER SEGMENT DU PARCOURS AVEC IMPACT
C                    IMPACT DU PARCOURS DE LA PARTICULE K
C                    TRACE d'UNE BOULE AU POINT D'IMPACT XYZPT1
C                    LE RAYON DE LA BOULE PARTICULE K SUIT LA VALEUR
C                    DU RAYON DE LA PARTICULE
                     RAYPAR = XVRVIPART( 7, K )
                     IF( RAYPARMX.GT.RAYPARMI ) THEN
                        RAYOSP = 0.006 + 0.01 *
     %                          (RAYPAR-RAYPARMI)/(RAYPARMX-RAYPARMI)
                     ELSE
                        RAYOSP = 0.016
                     ENDIF

C      CHOIX DE LA COULEUR DE TRACE DE LA SURFACE et des ARETES DE LA BOULE
C      COMMON / TRVAR1 / NOPACL, N1CORE, NDCORE, N1COEL, NDCOEL,
C     %                  N1COUL, NDCOUL,
C     %                  NCNOIR=0,  NCROUG=1,  NCVERT=2,  NCBLEU=3,
C     %                  NCCYAN=4,  NCJAUN=5,  NCMAGE=6,  NCBLAN=7,
C     %                  NCGRIS=8,  NCGRIM=9,  NCGRIC=10, NCBEIG=11,
C     %                  NCORAN=12, NCSAUM=13, NCROSE=14, NCTURQ=15
                     NCSFBOUL = 1 + NINT( 5 *
     %                         (RAYPAR-RAYPARMI)/(RAYPARMX-RAYPARMI) )

                     IF( NBXYZPART( K ) .GT. 0 ) THEN
C                       FIN DE PARCOURS de la PARTICULE K
C                       TRACE D'UNE BOULE A PARTIR D'UN ICOSAEDRE
                        CALL TRABOULE( RAYOSP, XYZPT1, NCSFBOUL, NCNOIR)
                     ELSE
C                       PARCOURS INACHEVE FAUTE DE TEMPS DE LA PARTICULE K
C                       TRACE D'UN TETRAEDRE EQUILATERAL d'ARETE RAYOSP
                        CALL TRATEFIN( RAYOSP*1.8, XYZPT1, NCSFBOUL )
                     ENDIF

C                    TRACE DU NUMERO K DE LA PARTICULE AU CENTRE DE LA BOULE
                     KNOPART(1:2) = '  '
                     CALL SANSBL( KNOPART, NBK )
                     CALL SYMBOLE3D( NCGRIS, XYZPT1, KNOPART(1:NBK) )

                  ENDIF

C                 FIN DU TRACE DU SEGMENT NITEM
                  GOTO 500

               ENDIF

            ENDDO

         ELSE

C           L'ITEM EST UN TRIANGLE D'UNE DES NBSFTR SURFACES
C           ------------------------------------------------
C           QUELLE SURFACE?
            DO KSF = NBPART+1, NBPART+NBSFTR

               IF( NITEM .LE. NODRITEM(KSF) ) THEN

C                 LE TRIANGLE NITEM APPARTIENT A LA SURFACE K
                  K = KSF - NBPART

C                 NO DU TRIANGLE DANS LA TRIAGULATION
                  NTR = NITEM - NODRITEM(KSF-1)
                  IF( NTR .LE. 0 ) THEN
                     PRINT*,'parpartr: ERREUR NTR=',NTR,' NITEM=',NITEM
                     GOTO 500
                  ENDIF

                  MNEF = MNSEFSF(K) + WUSOEF -5
                  MNXY = MNXYZSF(K) + WYZSOM -4 
                  DO I=1,3
C                    NO DU SOMMET I DU TRIANGLE NTR DE LA SURFACE K
                     NS = MCN( MNEF + 4*NTR + I )
                     NOSOTR(I) = NS
C                    XYZ DU SOMMET I DU TRIANGLE NTR
                     MNX = MNXY + 3*NS
                     DO M=1,3
                        XYZTRIA(M,I) = RMCN( MNX + M )
                     ENDDO
                  ENDDO

C                 CALCUL DE LA COULEUR DE LA FACE AVEC PRISE EN COMPTE
C                 DE LA DIRECTION DE VISEE ET L'ELOIGNEMENT
C                 ----------------------------------------------------
C                 LA COULEUR DE LA FACE EST UNIQUE
                  R  = ( DISTEMS(NIT) - DISTMIN ) / DISTMIMX
                  R2 = ( SQRT(R+EP1) - EP2 ) / EP3

C                 CNORFA LES COORDONNEES DE LA NORMALE A LA FACE (NORME=1)
                  MNSTS = MNXYZSF(K) + WYZSOM
                  CALL NORF34( 3, NOSOTR, RMCN(MNSTS),  CNORFA, IERR)
                  IF( IERR .NE. 0 ) GOTO 500
C
C                 LE PRODUIT SCALAIRE DIRECTION VISEE ET NORMALE A LA FACE
                  R = PROSCR( DIREVI, CNORFA, 3 )
                  R = 1. - ABS(R)

C                 LA COULEUR DE LA FACE PONDEREE PAR LE COSINUS
C                 DE L'ANGLE (NORMALE,DIRECTION DE VISEE) et ELOIGNEMENT
                  RCOUL = R * POID + POID1 * R2

                  IF( LCRITR .EQ. 0 ) THEN

C                    UNE COULEUR PAR TRIANGLE TRACE
                     NCOUFA = NINT( N1COULGRIS + NBCOULGRIS * RCOUL )

C                    TRACE DE LA FACE P1 SANS SES ARETES
                     CALL FACE3D( NCOUFA, -1, 3, XYZTRIA )

                  ELSE IF( LCRITR .LT. 0 ) THEN

C                    LA COULEUR EST INTERPOLEE SELON LA COULEUR AUX 3 SOMMETS
                     DO I=1,3
C                       LA COULEUR DU SOMMET I DU TRIANGLE NTR SELON
C                       LA COTE (SOLEIL A LA VERTICALE) DU SOMMET I
C                       PONDEREE PAR LA DIRECTION ET L'ELOIGNEMENT
                        COUL(I)=(2*RCOUL+(XYZTRIA(3,I)-ZMIN)/ZMINMAX)/3
                        IF( COUL(I) .LT. 0.0 ) COUL(I)=0.0
                        IF( COUL(I) .GT. 1.0 ) COUL(I)=1.0
                        COUL(I) = N1COULGRIS + NBCOULGRIS * COUL(I)

C                       PROTECTION a SUPPRIMER SI JAMAIS RENCONTRE...
                        IF( COUL(I) .LT. N1COULGRIS ) THEN
                           PRINT *,'parpartr: ZZZ COUL(',I,')=',COUL(I),
     %                             ' < ',N1COULGRIS
                           COUL(I) = N1COULGRIS
                        ENDIF
                        IF( COUL(I) .GT. NDCOUL ) THEN
                           PRINT *,'parpartr: ZZZ COUL(',I,')=',COUL(I),
     %                             ' >',NDCOUL
                           COUL(I) = NDCOUL
                        ENDIF
                     ENDDO

C                    LE TRACE EFFECTIF DU TRIANGLE DE SOMMETS 123
                     CALL TRIACOUL3D( XYZTRIA, COUL )

                  ENDIF

C                 FIN DU TRACE DE L'ITEM TRIANGLE NITEM
                  GOTO 500

               ENDIF

            ENDDO

         ENDIF

 500  ENDDO


C     LE TRACE FINAL DES PARCOURS DES PARTICULES AVEC LE TITRE
C     --------------------------------------------------------
      CALL LEGPARTI( KNOMOB, NBPART, NCAS0, 0., VITMXPAR )

C     MISE SUR FICHIER KnomfgifBoImage.xwd puis KnomfgifNoImage.jpg
C     DE LA PIXMAP de la FENETRE X11 ACTUELLE
      CALL VIDEO1( KNOMFGIF, NCAS0 )

C     ATTENDRE POUR LIRE LE TRACE
      CALL ATTENDSEC( TEMP2TRAC )

C     CONSTRUIRE le FICHIER VIDEO Nomfic.gif A PARTIR DES FICHIERS
C     CONSTRUITS de NOMS KnomfgifNoImag.jpg
      CALL VIDEOFIN( KNOMFGIF )


C     RETOUR POUR UNE NOUVELLE VISEE
C     ------------------------------
      IF( LORBITE .NE. 0 ) THEN
C        ORBITE BOUTON ENFONCE et DEPLACE?
         CALL ORBITE1( NOTYEV )
         IF( NOTYEV .NE. 0 ) GOTO 30
      ENDIF
      GOTO 10

C     SORTIE DU TRACE 3D DU PARCOURS DES PARTICULES
 9999 RETURN
      END
