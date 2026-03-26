      SUBROUTINE LIGLUE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : TRACER LES LIGNES DEJA LUES DANS KLG
C -----
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS        MARS 1990
C23456---------------------------------------------------------------012
      include"./incl/gsmenu.inc"
      include"./incl/trvari.inc"
      include"./incl/lu.inc"
      INTEGER        NA(1)
      INTRINSIC      INT
C
C     LE NUMERO MAXIMAL DU DERNIER CARACTERE NON BLANC DE KLG
      MDRECT(NRLGLU) = MXCANB( LHKLG, KLG )
C
C     LE COIN SUPERIEUR GAUCHE ET LA LARGEUR ET HAUTEUR
      DXRECT(NRLGLU) = MAX( DXRECT(NRMENU),
     %              2*ECARLR(NRLGLU)+MDRECT(NRLGLU)*(DXKARC(NRLGLU)+1) )
      DXRECT(NRLGLU) = MAX( DXRECT(NRMENU), DXRECT(NRLGLU) )
C     UNE LARGEUR MINIMALE
      DXRECT(NRLGLU) = MAX( 200, DXRECT(NRLGLU) )
      XRECT (NRLGLU) = LAPXFE - 2*ECARRC - DXRECT(NRLGLU)
C
      NBLGRC(NRLGLU) = MAX( 0, LHKLG-1 )
      DYRECT(NRLGLU) = DYLGRC(NRLGLU) * NBLGRC(NRLGLU)
      IF( DYRECT(NRLGLU) .GT. 0 ) THEN
         YRECT(NRLGLU) =INT(YRECT(NRMENU) + DYRECT(NRMENU) + 1.5*ECARRC)
      ELSE IF( DYRECT(NRERR) .GT. 0 ) THEN
         YRECT(NRLGLU) =INT(YRECT(NRERR) + DYRECT(NRERR) + 1.5*ECARRC)
      ELSE
         YRECT(NRLGLU) = 0
      ENDIF
C
      IF( MDRECT(NRLGLU) .LE. 0 ) RETURN
C
C     LE TRACE DU TEXTE
      NA(1) = 0
      CALL RECTTX( NRLGLU, KLG, 0, NA )
C
      RETURN
      END
