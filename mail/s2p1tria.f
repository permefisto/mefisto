       SUBROUTINE S2P1TRIA( NDIM, MOXYZS, MNXYZS, MONSEF, MNNSEF, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   CLIQUER 2 SOMMETS + 1 POINT EXTERIEUR ET EN FAIRE UN TRIANGLE
C -----

C ENTREES:
C --------
C NDIM   : DIMENSION (2 ou 3) DE L'ESPACE DE LA SURFACE
C MOXYZS : NOMBRE DE MOTS DECLARES DU TMS XYZSOMMET
C MNXYZS : ADRESSE DU TMS XYZSOMMET DE LA SURFACE
C MONSEF : NOMBRE DE MOTS DECLARES DU TMS NSEF
C MNNSEF : ADRESSE DU TMS NSEF DE LA SURFACE

C SORTIES:
C --------
C IERR   : 0 SI PAS D'ERREUR
C          1 SI POINT EN DEHORS DE LA TRIANGULATION
C          2 SI ECHANGE REFUSE POUR QUADRANGLE NON CONVEXE OU DEGENERE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : Alain PERRONNET Labo J-L. LIONS  UPMC   PARIS  SEPTEMBRE 2007
C MODIFS : Alain PERRONNET LJLL UPMC & StPIERRE du PERRAY SEPTEMBRE 2010
C2345X7..............................................................012
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/pp.inc"
      COMMON         MCN(MOTMCN)
      REAL          RMCN(1)
      EQUIVALENCE  (RMCN(1),MCN(1))
      INTEGER       NUMSTO(3)
      REAL          XYZPT(3), XYZST(3)

      IERR = 0
C     NOMBRE DE TRIANGLES ou QUADRANGLES du MAILLAGE
      NBTQ  = MCN( MNNSEF + WBEFOB )

C     NOMBRE DE SOMMETS du MAILLAGE
      NBSOM = MCN( MNXYZS + WNBSOM )

C     CLIC DES 2 SOMMETS DU MAILLAGE
C     ==============================
      DO K=2,3

C        RECHERCHE DU SOMMET LE PLUS PROCHE NUMSTO(K) DU POINT CLIQUE
C        ============================================================
         CALL SESTCLIC( RMCN(MNXYZS+WYZSOM), ITSTCLIC, NSCLIC )
         IF( NSCLIC .LE. 0 ) THEN
C           NSCLIC=0 SI LE POINT CLIQUE EST HORS MAILLAGE ou ABANDON DEMANDE
            IERR = 1
            GOTO 9999
         ENDIF

         NUMSTO(K) = NSCLIC
         MN = MNXYZS + WYZSOM + 3*NSCLIC -3
         PRINT*,'s2p1tria: SOMMET',K,' CLIQUE: No=',NSCLIC,
     %          ' X=',RMCN(MN),' Y=',RMCN(MN+1),' Z=',RMCN(MN+2)

      ENDDO

      IF( NDIM .EQ. 2 ) THEN

C        CLIC DU POINT EXTERIEUR AU MAILLAGE 2D
C        ======================================
         CALL SAIPTC( NOTYEV, NX, NY, NOCHAR )
         IF( NOTYEV .LE. 0 ) THEN
C           -1 => CARACTERE TAPE AU CLAVIER AVEC NOCHAR
C            0 => ABANDON
            IERR = 1
            GOTO 9999
         ENDIF

C        LES COORDONNEES DU POINT CLIQUE
         XYZPT(1) = XOB2PX( NX )
         XYZPT(2) = YOB2PX( NY )
         XYZPT(3) = 0
C
C        RECHERCHE DU TRIANGLE OU QUADRANGLE CONTENANT LE POINT
         MNSOEL = MNNSEF + WUSOEF
         CALL TQPTCLIC( NX,NY,NDIM,MCN(MNXYZS+WYZSOM), NBTQ,MCN(MNSOEL),
     %                NUMTQ )

C        DEJA DANS UN TRIANGLE ou QUADRANGLE?
         IF( NUMTQ .GT. 0 ) THEN
            NBLGRC(NRERR)=2
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1)='POINT INTERNE AU MAILLAGE'
               KERR(2)='CHOISIR LE POINT A L''EXTERIEUR DU MAILLAGE'
            ELSE
               KERR(1)='POINT INSIDE THE MESH'
               KERR(2)='CHOOSE IT AT THE EXTERIOR OF THE MESH'
            ENDIF
            CALL LEREUR
            IERR=1
            GOTO 9999
         ENDIF

      ELSE

C        CLIC DU POINT EXTERIEUR AU MAILLAGE 3D
C        ======================================
C        EN 3D: LECTURE DES 3 COORDONNEES XYZ DU POINT EXTERNE AU MAILLAGE
         CALL INVITE( 142 )
         NCVALS = 0
         CALL LIRXYZ( NCVALS, XYZPT )
         IF( NCVALS .LE. 0 ) GOTO 9999

C        LE POINT EST IL DANS UN TRIANGLE OU QUADRANGLE VISIBLE?
         CALL XYZAXO( XYZPT, XYZST )
         NX = NUPXEX( XYZST(1) )
         NY = NUPXEY( XYZST(2) )
         CALL NOSTCLIC( NX, NY, NDIM, NBSOM, RMCN(MNXYZS+WYZSOM),
     %                  NBTQ,  MCN(MNNSEF+WUSOEF),
     %                  NUMTQ, XYZST, NUMSTO(1) )
         IF( NUMSTO(1) .LE. 0 ) GOTO 9999
         IF( NUMTQ .GT. 0 ) THEN
            NBLGRC(NRERR)=2
            IF( LANGAG .EQ. 0 ) THEN
           KERR(1)='ATTENTION: Le POINT PROJETE est INTERNE AU MAILLAGE'
           KERR(2)='IL PEUT ETRE CORRECT mais A VERIFIER'
            ELSE
             KERR(1)='ATTENTION: The PROJECTED POINT is INSIDE THE MESH'
             KERR(2)='THIS MAY BE WELL but VERIFY IT'
            ENDIF
            CALL LERESU
         ENDIF

      ENDIF

C     LE POINT EXTERNE EST STOCKE AU DELA DES SOMMETS ACTUELS
C     MISE A JOUR DU TMS XYZSOMMET
C     =======================================================
      L = WYZSOM + 3*NBSOM
C     ON AJUSTE LA TAILLE DU TABLEAU
      IF ( MOXYZS .LT. L+3 ) THEN
         CALL TNMCAU( 'ENTIER', MOXYZS, L+3, MOXYZS, MNXYZS )
         IF( MNXYZS .LE. 0 ) THEN
C           PAS ASSEZ DE MCN
            IERR = 2
            GOTO 9999
         ENDIF
         MOXYZS = L+3
         NBSOM  = MCN(MNXYZS+WNBSOM)
      ENDIF
C     STOCKAGE DE XYZ DU NOUVEAU POINT CLIQUE EXTERIEUR AU MAILLAGE
      RMCN(MNXYZS+L  ) = XYZPT(1)
      RMCN(MNXYZS+L+1) = XYZPT(2)
      RMCN(MNXYZS+L+2) = XYZPT(3)
      NBSOM = NBSOM + 1
      MCN(MNXYZS+WNBSOM) = NBSOM
C     NBSOM NUMERO DE SOMMET DU POINT CLIQUE
      NUMSTO(1) = NBSOM

      CALL ECDATE(MCN(MNXYZS))
      MCN(MNXYZS+MOTVAR(6))=NONMTD('~>>>XYZSOMMET')

C     LE TRIANGLE PT 2 SOMMETS EST IL DE SURFACE POSITIVE?
      IF( NDIM .EQ. 2 ) THEN
         MN = MNXYZS + WYZSOM - 3
         SURF = SURTR2( RMCN(MN+3*NUMSTO(1)),
     %                  RMCN(MN+3*NUMSTO(2)),
     %                  RMCN(MN+3*NUMSTO(3)) )
         IF( SURF .LE. 0 ) THEN
C           PAR PERMUTATION ST 2-3 LA SURFACE DEVIENT POSITIVE
            K         = NUMSTO(2)
            NUMSTO(2) = NUMSTO(3)
            NUMSTO(3) = K
         ENDIF
      ENDIF

C     AJOUT DU TRIANGLE NUMSTO(1:3) AU TMS NSEF
C     =========================================
      L = WUSOEF + 4 * NBTQ
      IF ( MONSEF .LT. L+4 ) THEN
         CALL TNMCAU( 'ENTIER', L, L+4, L, MNNSEF )
         MONSEF = L+4
      ENDIF
      MN = MNNSEF + WUSOEF + 4 * NBTQ - 1
      DO 40 K=1,3
         MCN(MN+K) = NUMSTO(K)
 40   CONTINUE
      MCN(MN+4) = 0

C     UN TRIANGLE DE PLUS
      NBTQ = NBTQ + 1
      MCN(MNNSEF+WBEFOB) = NBTQ

      CALL ECDATE( MCN(MNNSEF) )
      MCN(MNNSEF+MOTVAR(6))=NONMTD( '~>>>NSEF' )

 9999 RETURN
      END
