SUBROUTINE BGVITBF( tn, tn1, Rho, CoGrPr, NDIM, NBSOM, NBNOVI, &
                    NDDLNO, XYZSOM, &
                    NUTYEL, NBNOEF, NBELEM, NUNOEF, NTDLHB, &
                    NBVVEF, NBSFEF, NBLAEF, &
                    NUVVEF, NUSFEF, NULAEF, &
                    NUMILI, NUMALI, LTDELI, &
                    NUMISU, NUMASU, LTDESU, &
                    NUMIVO, NUMAVO, LTDEVO, &
                    MOARET, MXARET, MNLARE, &
                    MOFACE, MXFACE, LFACES, &
                    VXYZPNtn,VITMXtn, &
                    PRESStn, NOBARY, AGGAUSS, &
                    BG,      NBCHTL, IERR )
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! BUT: CONSTRUIRE ET ASSEMBLER LE SECOND MEMBRE GLOBAL BG DU PROBLEME
! ---- VG U(tn+1) = (Rho - dt Mhu Laplacien ) U(tn+1) =
!      = Rho U(tn,X(tn+1,x) - dt CoGrPr GRAD P(tn) + dt Force(tn+1)
!      CALCUL DU SECOND MEMBRE DU SYSTEME LINEAIRE:
!      Integrale Rho tP1B Ui(tn,X(tn;tn+1,bl) dX
!    - Integrale tP1B dt CoGrPr GRAD P(tn) dX + dt Integrale tP1B Force(tn) dX

! ENTREES:
! --------
! tn, tn1: INTERVALLE DE TEMPS tn tn+1  => DT PAS du TEMPS = tn1-tn
! DT     : PAS CONSTANT DU TEMPS
! Rho    : DENSITE DE MASSE  Integrale Rho P1B P1B dX
! CoGrPr : COEFFICIENT DU GRADIENT DE PRESSION DANS LES EDP
! NDIM   : DIMENSION DE L'ESPACE DES COORDONNEES 2 ou 3
! NBSOM  : NOMBRE DE SOMMETS DU MAILLAGE
! NBNOVI : NOMBRE DE NOEUDS VITESSE SOMMETS et BARYCENTRE DES EF
! XYZSOM : XYZSOM(3,NBNOVI)  3 COORDONNEES DES NOEUDS DU MAILLAGE
! NDDLNO : NUMERO DU DERNIER DL DE CHAQUE NOEUD VITESSE (0:NBNOVI)
!          A UTILISER POUR RETROUVER LES DL PRESSION DES VECTEURS VXYZPNtn
! NUTYEL : NUMERO DU TYPE D'EF ( 13 TRIANGLE, 19 TETRAEDRE BREZZI-FORTIN
! NBNOEF : NOMBRE DE SOMMETS D'UN EF DE CE TYPE
!          ( 3 POUR LE TRIANGLE, 4 POUR LE TETRAEDRE)
! NBELEM : NOMBRE D'EF DU MAILLAGE
! NUNOEF : NUNOEF(NBELEM,NBNOEF) NO DES NOEUDS DES NBELEM EF

! NBVVEF : NOMBRE DE VOLUME  POUR 1 EF DE CE TYPE
! NBSFEF : NOMBRE DE FACES   POUR 1 EF DE CE TYPE
! NBLAEF : NOMBRE D'ARETES   POUR 1 EF DE CE TYPE
! NUVVEF : TABLEAU DU NUMERO DE VOLUME DE CHAQUE ELEMENT FINI
! NUSFEF : TABLEAU DU NUMERO DE SURFACE DES FACES   DE CHAQUE EF
! NULAEF : TABLEAU DU NUMERO DE LIGNE   DES ARETES  DE CHAQUE EF

! NUMILI : NUMERO MINIMAL DES OBJETS LIGNES
! NUMALI : NUMERO MAXIMAL DES OBJETS LIGNES
! LTDELI : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES DU FLUIDE
!          DES OBJETS LIGNES

! NUMISU : NUMERO MINIMAL DES OBJETS SURFACES
! NUMASU : NUMERO MAXIMAL DES OBJETS SURFACES
! LTDESU : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES DU FLUIDE
!          DES OBJETS SURFACES

! NUMIVO : NUMERO MINIMAL DES OBJETS VOLUMES
! NUMAVO : NUMERO MAXIMAL DES OBJETS VOLUMES
! LTDEVO : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES DU FLUIDE
!          DES OBJETS VOLUMES

! MOARET : NOMBRE DE MOTS POUR LES INFO D'UNE ARETE
! MXARET : NOMBRE MAXIMAL D'ARETES DECLAREES
! MNLARE : EN 2D SEULEMENT ADRESSE MCN DU 1-ER MOT DU TABLEAU LARETE

! MOFACE : NOMBRE DE MOTS POUR LES INFO D'UNE FACE
! MXFACE : NOMBRE MAXIMAL D'FACES DECLAREES
! LFACES : TABLEAU DES FACES DU MAILLAGE  (cf hachag.f)
!          LFACES(1,I)= NO DU 1-ER  SOMMET DE LA FACE
!          LFACES(2,I)= NO DU 2-EME SOMMET > 1-ER  SOMMET
!          LFACES(3,I)= NO DU 3-EME SOMMET DE LA FACE
!          LFACES(4,I)= NO DU 4-EME SOMMET DE LA FACE
!                       0 SI TRIANGLE
!          LFACES(5,I)= CHAINAGE HACHAGE SUR FACE SUIVANTE
!                       1 CUBE EST UN TETRA ou PENTA ou HEXAEDRE
!          LFACES(6,I)= NUMERO DU 1-ER  CUBE CONTENANT CETTE FACE
!                       0 SI PAS DE 1-ER  CUBE
!          LFACES(7,I)= NUMERO DU 2-EME CUBE CONTENANT CETTE FACE
!                       ou CHAINAGE SUR LA FACE FRONTALIERE SUIVANTE
!          LFACES(8,I)= NUMERO DE LA FACE A TANGENTE DANS LE TABLEAU NUTGFA
!          LFACES(9,I)= NUMERO DU TYPE DU DERNIER EF DE CETTE FACE

! VXYZPNtn: VECTEUR GLOBAL des DL VITESSES-PRESSIONS (1:NTDL) au TEMPS tn
! VITMXtn : NORME DE LA VITESSE MAXIMALE AU TEMPS PRE!EDANT
! PRESStn : DL DE LA PRESSION AUX SOMMETS DE LA TETRAEDRISATION
! NOBARY : NBNOEF+1 NO DU DL BARYCENTRE DANS L'EF
! AGGAUSS: LIGNE DE AE DU DEGRE DE LIBERTE AU BARYCENTRE
! NTDLVP : NOMBRE TOTAL DE DL VITESSE (NDIM*NBNOVI) + PRESSION NBSOM

! SORTIES:
! --------
! BG     : VECTEUR GLOBAL SECOND MEMBRE ASSEMBLE
! NBCHTL : NOMBRE DE REMONTEES DE LA CARACTERISTIQUE TROP LONGUES
! IERR   : 0 PAS DE PROBLEME, >0 PROBLEME RENCONTRE
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! AUTEUR : ALAIN PERRONNET LJLL UPMC & St Pierre du Perray  Octobre 2012
!23456---------------------------------------------------------------012
!$ USE OMP_LIB
      IMPLICIT NONE
      include"./incl/langue.inc"
      include"./incl/threads.inc"
      include"./incl/donflu.inc95"
      REAL              DT, tn1, tn, XYZSOM(3,NBSOM)
      DOUBLE PRECISION  DELTAT, Rho, CoGrPr, &
                        BG(NBNOVI*NDIM), AGGAUSS(NOBARY,NBELEM)

      INTEGER           NUNOEF(NBELEM,NBNOEF), NDDLNO(0:NBNOVI)
      INTEGER           NUMILI,NUMALI,LTDELI(1:MXDOFL,NUMILI:NUMALI)
      INTEGER           NUMISU,NUMASU,LTDESU(1:MXDOFL,NUMISU:NUMASU)
      INTEGER           NUMIVO,NUMAVO,LTDEVO(1:MXDOFL,NUMIVO:NUMAVO)
      INTEGER           NUVVEF(NBVVEF,NBELEM)
      INTEGER           NUSFEF(NBSFEF,NBELEM)
      INTEGER           NULAEF(NBLAEF,NBELEM)
      INTEGER           LFACES(MOFACE,MXFACE)
      INTEGER           NDIM,   NBSOM, NOBARY, NTDLHB, &
                        NBVVEF, NBSFEF, NBLAEF, &
                        NUTYEL, NBELEM, NBNOEF, &
                        MOARET, MXARET, MNLARE, &
                        MOFACE, MXFACE, NBNOVI, NBCHTL, IERR, IER
      INTEGER           NUELEM, M, K, MEK, NSK, NONOEF(5)

      DOUBLE PRECISION  BE1(19), BE2(15), D
      DOUBLE PRECISION  PRESStn(NBSOM), VXYZPNtn(1:*), VITMXtn
      DOUBLE PRECISION  t_cpu_0, t_cpu_1

      call cpu_time(t_cpu_0)
!!!      print*,'ENTREE de bgvitbf.f95 avec 1 ATOMIC sur BG'

!     LE PAS DE TEMPS EN DOUBLE PRECISION
      DT     = tn1 - tn
      DELTAT = DT

!     MISE A ZERO DE BG
      CALL AZEROD( NBNOVI*NDIM, BG )

!///////////////////////////////////////////////////////////////////////
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(NUELEM,M,K,MEK,NSK,NONOEF,BE1,BE2,D,IER)
!     BOUCLE SUR LES EF de BREZZI-FORTIN
!$OMP DO SCHEDULE( STATIC, NBELEM/NBTHREADS )
      DO NUELEM = 1, NBELEM

!        NO DES NOEUDS DE L'ELEMENT FINI NUELEM
         DO M = 1, NBNOEF
            NONOEF(M) = NUNOEF(NUELEM,M)
         ENDDO
!        NO DU BARYCENTRE
         NONOEF(NOBARY) = NBSOM + NUELEM
!
         IF( NUTYEL .EQ. 13 ) THEN

!           TRIANGLE BREZZI_FORTIN
!           1) Integrale Rho tP1B ui(tn,X(tn;tn+1,x) dx
!           SECOND MEMBRE ELEMENTAIRE du TRANSPORT pour NAVIER-STOKES
!           avec 4+4 COEFFICIENTS
            CALL F2SNVP1BP1( NUELEM, 1,       NBELEM, NUNOEF, &
                             MOARET, MXARET,  MNLARE, &
                             DT,     Rho,     XYZSOM, NDDLNO, &
                             NTDLHB, VXYZPNtn, VITMXtn, &
                             BE1,    IER )
            IF( IER .EQ. 4 ) THEN
!$OMP ATOMIC
               NBCHTL = NBCHTL + 1
            ENDIF

!           ASSEMBLAGE DES DEGRES DE LIBERTE 4 ET 8 DE BE1 DANS BG
!           COMPRESSION DE 8 A 6 DEGRES DE LIBERTE DE BE1
            CALL GAS2P1BP1( NBSOM, NUELEM, BE1, NBNOVI, BG )

!           2) -Integrale tP1B dt CoGrPr GRAD P(tn) + dt tP1B Force(tn) dX
            CALL F2EX2P1BP1( XYZSOM, NONOEF, &
                             NULAEF(1,NUELEM), NUMILI, NUMALI, LTDELI, &
                             NUSFEF(1,NUELEM), NUMISU, NUMASU, LTDESU, &
                             DELTAT, CoGrPr,   NBSOM,  PRESStn, &
                             BE2 )

!           ASSEMBLAGE DES DEGRES DE LIBERTE 4 ET 8 DE BE2 DANS BG
!           COMPRESSION DE 8 A 6 DEGRES DE LIBERTE DE BE2
            CALL GAS2P1BP1( NBSOM, NUELEM, BE2, NBNOVI, BG )

         ELSE

!           TETRAEDRE BREZZI-FORTIN
!           1) Integrale Rho tP1B ui(tn,X(tn;tn+1,x) dx
!           SECOND MEMBRE ELEMENTAIRE du TRANSPORT pour NAVIER-STOKES
!           avec 5+5+5 COEFFICIENTS
            CALL F3SNVP1BP1( NUELEM, NBELEM,  NUNOEF, &
                             MOFACE, MXFACE,  LFACES, &
                             DT,     Rho,     XYZSOM, NDDLNO, &
                             NTDLHB, VXYZPNtn, VITMXtn, &
                             BE1,    IER )
            IF( IER .EQ. 4 ) THEN
!              NOMBRE DE TETRAEDRES PARCOURUS TROP GRAND
!$OMP ATOMIC
               NBCHTL = NBCHTL + 1
            ENDIF

!           ASSEMBLAGE  DES DEGRES DE LIBERTE 5 10 15 DE BE1 DANS BG
!           COMPRESSION DES DEGRES DE LIBERTE 5 10 15 DE BE1
            CALL GAS3P1BP1( NBSOM, NUELEM, BE1, NBNOVI, BG )

!           2) -Integrale tP1B dt CoGrPr GRAD P(tn) + dt tP1B Force(tn) dX
            CALL F3EX2P1BP1( XYZSOM, NONOEF, &
                             NUSFEF(1,NUELEM), NUMISU, NUMASU, LTDESU, &
                             NUVVEF(1,NUELEM), NUMIVO, NUMAVO, LTDEVO, &
                             DELTAT, CoGrPr,   NBSOM,  PRESStn, &
                             BE2 )

!           ASSEMBLAGE  DES DEGRES DE LIBERTE 5 10 15 DE BE2 DANS BG
!           COMPRESSION DES DEGRES DE LIBERTE 5 10 15 DE BE2
            CALL GAS3P1BP1( NBSOM, NUELEM, BE2, NBNOVI, BG )

         ENDIF

!!!         IF( NUELEM .EQ. 5 ) THEN
!!!            call affvect( 'NEF=5 BE1 X(tn;tn+1,x)', NDIM*NOBARY, BE1 )
!!!            call affvect( 'NEF=5 BE2', NDIM*NOBARY, BE2 )
!!!         ENDIF

!        ASSEMBLAGE DE BE1 + BE2 dans BG = VITX1, VITY1, VITZ1
!        Rho U(tn,X(tn+1,x)) - dt CoGrPr GRAD P(tn) + dt FOmega(tn+1)
         DO M=1,NBNOEF

!           NUMERO DU M-EME NOEUD DE L'ELEMENT FINI
            NSK = NONOEF( M )
            MEK = M

!           ASSEMBLAGE DES 2 COEFFICIENTS ELEMENTAIRES DANS LE COEFFICIENT GLOBAL
            DO K=1,NDIM
               D = BE1( MEK ) + BE2( MEK )
!$OMP ATOMIC
               BG( NSK ) = BG( NSK ) + D

               NSK = NSK + NBNOVI
               MEK = MEK + NBNOEF
            ENDDO

         ENDDO

      ENDDO     ! FIN DE LA BOUCLE DES EF POUR ASSEMBLER BG
!$OMP END DO
!$OMP BARRIER
!$OMP END PARALLEL
!///////////////////////////////////////////////////////////////////////

!///////////////////////////////////////////////////////////////////////
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( NUELEM, M, NONOEF )
!     LA BOUCLE SUR LES ELEMENTS FINIS BREZZI-FORTIN POUR FAIRE GAUSS
!     SUR LE SECOND MEMBRE DU SYSTEME LINEAIRE A RESOUDRE
!     ---------------------------------------------------------------
!$OMP DO SCHEDULE( STATIC, NBELEM/NBTHREADS )
      DO NUELEM = 1, NBELEM

!        NO DES NOEUDS DE L'ELEMENT FINI NUELEM
         DO M = 1, NBNOEF
            NONOEF(M) = NUNOEF(NUELEM,M)
         ENDDO
!!!      NO DU BARYCENTRE
!!!      NONOEF(NOBARY) = NBSOM + NUELEM

         IF( NUTYEL .EQ. 13 ) THEN

!           TRIANGLE BREZZI-FORTIN P1BULLEP3 EN VITESSE et P1 EN PRESSION
!           -------------------------------------------------------------
!           METHODE DE GAUSS SUR LES 2 DL DU BARYCENTRE SUR LE SECOND
!           MEMBRE DU SYSTEME A RESOUDRE
            CALL GAR2P1BP1( NBSOM, NUELEM, NONOEF, AGGAUSS(1,NUELEM), &
                            NBNOVI, BG )

         ELSE

!           TETRAEDRE BREZZI-FORTIN P1BULLEP4 EN VITESSE et P1 EN PRESSION
!           --------------------------------------------------------------
!           METHODE DE GAUSS SUR LES 3 DL DU BARYCENTRE SUR LE SECOND
!           MEMBRE DU SYSTEME A RESOUDRE
            CALL GAR3P1BP1( NBSOM, NUELEM, NONOEF, AGGAUSS(1,NUELEM), &
                            NBNOVI, BG )

         ENDIF

      ENDDO     ! FIN DE LA BOUCLE DES EF POUR ASSEMBLER BG
!$OMP END DO
!$OMP BARRIER
!$OMP END PARALLEL
!///////////////////////////////////////////////////////////////////////

!!!      call affvect( 'BGVITBF BG=', 20,          BG )
!!!      call afl1ve(  'BGVITBF BG=', NDIM*NBNOVI, BG )


!     Le TEMPS CPU donne par le systeme en fin d'execution de bgvitth.f
      call cpu_time(t_cpu_1)
 
      IF( LANGAG .EQ. 0 ) THEN
         print*,'SORTIE de bgvitbf.f95 avec 1 ATOMIC sur BG en',&
                 t_cpu_1-t_cpu_0,' CPU secondes'
      ELSE
         print*,'EXIT bgvitbf.f95 with 1 ATOMIC on BG. CPU seconds=',&
                 t_cpu_1-t_cpu_0
      ENDIF

      IERR = 0
      RETURN
END SUBROUTINE BGVITBF
