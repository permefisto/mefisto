      SUBROUTINE TRAFAC( NMSURF, NUSURF, MNNSEF, MNXYZS )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACER LES FACES (TRIANGLES ET/OU QUADRANGLES) D'UNE SURFACE
C -----
C ENTREES:
C --------
C NMSURF : NOM DE LA SURFACE A TRACER
C NUSURF : NUMERO DE SURFACE DANS LE LEXIQUE DES SURFACES
C MNNSEF : ADRESSE MCN DU TABLEAU 'NSEF'      DE LA SURFACE A TRACER
C MNXYZS : ADRESSE MCN DU TABLEAU 'XYZSOMMET' DE LA SURFACE A TRACER
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR:PERRONNET ALAIN ANALYSE NUMERIQUE LJLL UPMC PARIS     MARS 1991
C MODIF :PERRONNET ALAIN ANALYSE NUMERIQUE LJLL UPMC PARIS NOVEMBRE 2003
C2345X...............................................................012
      include"./incl/gsmenu.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      include"./incl/traaxe.inc"
      include"./incl/xyzext.inc"
      include"./incl/a___xyzsommet.inc"
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
      CHARACTER*(*)     NMSURF
      CHARACTER*32      KNSURF

C     EPAISSEUR DES TRAITS
      CALL XVEPAISSEUR( NEPARF )
C     NUMERO DE LA COULEUR DE TRACE DES ARETES DES FACES
      NCOUAF = NCGRIS

C     L'OPTION FINALE DE TRACE DES FACES
      IF( NDIMLI .EQ. 2 ) THEN

C        TRACE DES FACES 2D SELON OU NON LEUR QUALITE
C        ============================================
         IF( LORBITE .LE. 0 ) GOTO 20

C        INITIALISATION DU ZOOM DEPLACEMENT
         CALL ZOOM2D0( NOTYEV )
         IF( NOTYEV .EQ. 0 ) GOTO 9000
         GOTO 20
C
C        ZOOM OU TRANSLATION ACTIFS
 10      CALL ZOOM2D1( NOTYEV )
         IF( NOTYEV .EQ. 0 ) GOTO 9000
C
C        TRACE EFFECTIF
 20      NETAXE = 0
         CALL TRAXE2
         CALL T21FAC( NMSURF, NUSURF, MNNSEF, MNXYZS )
         L = NUDCNB( NMSURF )
         KNSURF = 'SURFACE ' // NMSURF(1:L)
         L = NUDCNB( KNSURF )
         CALL TRFINS( KNSURF(1:L) )
C
C        REPRISE DE TRANSLATION ZOOM
         IF( LORBITE .GT. 0 ) GOTO 10
C
      ELSE
C
C        TRACE DES FACES 3D SELON OU NON LEUR QUALITE
C        ============================================
         IF( LORBITE .LE. 0 ) GOTO 40
C
C        INITIALISATION DE L'ORBITE ZOOM DEPLACEMENT
         CALL ORBITE0( NOTYEV )
         GOTO 40
C
C        ORBITE OU ZOOM OU TRANSLATION ACTIFS
 30      CALL ORBITE1( NOTYEV )
         IF( NOTYEV .EQ. 0 ) GOTO 9000
         NETAXE = 0
         CALL TRAXE3
C
C        LE TRACE EFFECTIF
 40      IF( IAVFAC .EQ. 0 ) THEN
C
C           TRACE FIL DE FER DES ARETES DES FACES FRONTIERES
            CALL T31FDF( NMSURF, NUSURF, MNNSEF, MNXYZS )
C
         ELSE
C
C           TRACE DES FACES AVEC ARETES OU NON, REDUITES OU NON, ...
            CALL T31FCO( NMSURF, NUSURF, MNNSEF, MNXYZS )
C
         ENDIF

C        TRACE EFFECTIF DU MAILLAGE SUR L'ECRAN
C        --------------------------------------
         L = NUDCNB( NMSURF )
         KNSURF = 'SURFACE ' // NMSURF(1:L)
         L = NUDCNB( KNSURF )
         CALL TRFINS( KNSURF(1:L) )
C
C        REPRISE DE L'ORBITE
C        ===================
         IF( LORBITE .GT. 0 ) GOTO 30
      ENDIF
C
 9000 RETURN
      END
