      SUBROUTINE RECTTR( NUR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : TRACE DU RECTANGLE DE NUMERO NURECT
C ----- VERSION xvue
C
C ENTREE :
C --------
C NUR : NUMERO DU RECTANGLE A TRACER
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS    JANVIER 1997
C23456---------------------------------------------------------------012
      include"./incl/epombr.inc"
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/gsmenu.inc"
      include"./incl/trvari.inc"
      INTEGER   DY
C
      IF( INTERA .LT. 3 ) RETURN
C
      IF( DYRECT(NUR) .LE. 0. ) RETURN
      IF( DXRECT(NUR) .LE. 0. ) RETURN
C
C     TRACE DU REMPLISSAGE DU RECTANGLE
C     =================================
      CALL XVCOULEUR( NFRECT( NUR ) )
      CALL XVFRECTANGLE( XRECT(NUR),  YRECT(NUR),
     %                   DXRECT(NUR), DYRECT(NUR) )
C
      IF( NUR .EQ. NRMENU ) THEN
C        COULEUR SPECIALE POUR CERTAINES LIGNES DU MENU
C        NOMBRE DE LIGNES DU MENU
         NBL = NBLGRC( NRMENU )
C        HAUTEUR D'UNE LIGNE
         DY  = DYRECT( NRMENU ) / NBL
C        LA PREMIERE LIGNE DE MENU: TITRE DU MENU
         CALL XVCOULEUR( NCSAUM )
         CALL XVFRECTANGLE( XRECT(NRMENU),  YRECT(NRMENU),
     %                      DXRECT(NRMENU), DY )
C        LA LIGNE ABANDON
         CALL XVCOULEUR( NCVERT )
         CALL XVFRECTANGLE( XRECT(NRMENU),  YRECT(NRMENU)+(NBL-2)*DY,
     %                      DXRECT(NRMENU), DY )
C        LA LIGNE DOCUMENTATION
         CALL XVCOULEUR( NCBLEU )
         CALL XVFRECTANGLE( XRECT(NRMENU),  YRECT(NRMENU)+(NBL-1)*DY,
     %                      DXRECT(NRMENU), DY )
      ENDIF
C
C     TRACE DE L'OMBRAGE DU RECTANGLE
C     ===============================
      CALL XVEPAISSEUR( NBEPOM )
      CALL XVCOULEUR( NCNOIR )
C
C     LES COORDONNEES DES 2 SOMMETS EXTREMAUX DU RECTANGLE
      NPXMI = XRECT( NUR ) + NBEPOM
      NPYMI = YRECT( NUR ) + DYRECT( NUR ) + NBEPO2 + 2
      NPXMA = XRECT( NUR ) + DXRECT( NUR ) + NBEPO2 + 1
      NPYMA = YRECT( NUR ) + NBEPOM
      CALL XVFTRAIT( NPXMI, NPYMI, NPXMA+NBEPO2, NPYMI )
      CALL XVFTRAIT( NPXMA, NPYMI, NPXMA,        NPYMA )
C
C     TRACE DE LA LIGNE DU BORD DU RECTANGLE
C     ======================================
      CALL XVCOULEUR( NBRECT( NUR ) )
      CALL XVEPAISSEUR( NBEPO2 )
      CALL XVFBORDRECTANGLE( XRECT(NUR),  YRECT(NUR),
     %                       DXRECT(NUR), DYRECT(NUR) )
C
C     EFFET DE CADRE ECLAIRE EN HAUT A GAUCHE
C     =======================================
      CALL XVEPAISSEUR( 2 )
C
C     LES TRAITS BLANCS
      CALL XVCOULEUR( NCBLAN )
C     GAUCHE DE LA GAUCHE DE BAS EN HAUT
      CALL XVFTRAIT( XRECT(NUR)-NBEPO4,
     %               YRECT(NUR)+DYRECT(NUR)+NBEPO4,
     %               XRECT(NUR)-NBEPO4,
     %               YRECT(NUR)-NBEPO4 )
C     HAUT DU HAUT DE GAUCHE A DROITE
      CALL XVFTRAIT( XRECT(NUR)-NBEPO4,
     %               YRECT(NUR)-NBEPO4,
     %               XRECT(NUR)+DXRECT(NUR)+NBEPO4,
     %               YRECT(NUR)-NBEPO4 )
C     HAUT DU BAS DE GAUCHE A DROITE
      CALL XVFTRAIT( XRECT(NUR)+NBEPO4,
     %               YRECT(NUR)+DYRECT(NUR)-NBEPO4,
     %               XRECT(NUR)+DXRECT(NUR)-NBEPO4,
     %               YRECT(NUR)+DYRECT(NUR)-NBEPO4 )
C     GAUCHE DE LA DROITE DE BAS EN HAUT
      CALL XVFTRAIT( XRECT(NUR)+DXRECT(NUR)-NBEPO4,
     %               YRECT(NUR)+DYRECT(NUR)-NBEPO4,
     %               XRECT(NUR)+DXRECT(NUR)-NBEPO4,
     %               YRECT(NUR)+NBEPO4 )
C
C     LES TRAITS NOIRS
      CALL XVCOULEUR( NCNOIR )
C     BAS DU HAUT DE GAUCHE A DROITE
      CALL XVFTRAIT( XRECT(NUR)+NBEPO4,
     %               YRECT(NUR)+NBEPO4,
     %               XRECT(NUR)+DXRECT(NUR)-NBEPO4,
     %               YRECT(NUR)+NBEPO4 )
C     DROITE DE LA GAUCHE DU HAUT EN BAS
      CALL XVFTRAIT( XRECT(NUR)+NBEPO4,
     %               YRECT(NUR)+NBEPO4,
     %               XRECT(NUR)+NBEPO4,
     %               YRECT(NUR)+DYRECT(NUR)-NBEPO4 )
C     BAS DU BAS DE GAUCHE A DROITE
      CALL XVFTRAIT( XRECT(NUR)-NBEPO4,
     %               YRECT(NUR)+DYRECT(NUR)+NBEPO4,
     %               XRECT(NUR)+DXRECT(NUR)+NBEPO4,
     %               YRECT(NUR)+DYRECT(NUR)+NBEPO4 )
C     DROITE DE LA DROITE DE BAS EN HAUT
      CALL XVFTRAIT( XRECT(NUR)+DXRECT(NUR)+NBEPO4,
     %               YRECT(NUR)+DYRECT(NUR)+NBEPO4,
     %               XRECT(NUR)+DXRECT(NUR)+NBEPO4,
     %               YRECT(NUR)-NBEPO4 )
C
C     LE TRAIT NOIR DU BAS DE L'OMBRAGE
      CALL XVEPAISSEUR( 0 )
      NPXMA = NPXMA + 3
      NPYMI = NPYMI + 3
      CALL XVFTRAIT( NPXMI, NPYMI, NPXMA, NPYMI )
      CALL XVFTRAIT( NPXMA, NPYMI, NPXMA, NPYMA )
      END
