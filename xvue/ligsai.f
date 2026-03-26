      SUBROUTINE LIGSAI
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : TRACER LA LIGNE DE SAISIE LUE DANS KLGSA
C -----
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS        MARS 1996
C23456---------------------------------------------------------------012
      include"./incl/gsmenu.inc"
      include"./incl/trvari.inc"
      include"./incl/langue.inc"
      include"./incl/lu.inc"
      INTEGER     NA(1)
      INTRINSIC   INT
C
C     AJOUT DU TEXTE INITIAL
      IF( LANGAG .EQ. 0 ) THEN
         KLGSA(1) = 'Clavier:' // KLG(LHKLG)
      ELSE
         KLGSA(1) = 'Keyboard:' // KLG(LHKLG)
      ENDIF
      NBLGRC(NRLGSA) = 1
C
C     LE NUMERO MAXIMAL DU DERNIER CARACTERE NON BLANC DE LA LIGNE
      MDRECT(NRLGSA) = MAX( 20, MXCANB( 1, KLGSA ) )
C
C     LE NOMBRE DE PIXELS DU TEXTE DE KLGSA
      IF( MDRECT(NRLGSA) .EQ. 20 ) THEN
         LAPX = DXKARC(NRLGSA) * MDRECT(NRLGSA)
      ELSE
         LAPX = LAMXPXTXT(1,KLGSA) + DXKARC(NRLGSA)
      ENDIF
C
C     LA LARGEUR EN PIXELS ET LA POSITION EN X
      DXRECT(NRLGSA) = ECARLR(NRLGSA) * 2 + LAPX
      DXRECT(NRLGSA) = MAX( DXRECT(NRMENU),
     %                      DXRECT(NRLGLU),
     %                      DXRECT(NRINVI),
     %                      DXRECT(NRLGSA) )
      XRECT (NRLGSA) = LAPXFE - 2*ECARRC - DXRECT(NRLGSA)
C
C     LA HAUTEUR EN PIXELS ET LA POSITION EN Y
      NBLGRC(NRLGSA) = 1
       YRECT(NRLGSA) = 0
      DYRECT(NRLGSA) = DYLGRC(NRLGSA) * NBLGRC(NRLGSA)
      IF( DYRECT(NRLGSA) .GT. 0 ) THEN
         IF( DYRECT(NRINVI) .GT. 0 ) THEN
            YRECT(NRLGSA) =INT(YRECT(NRINVI)+DYRECT(NRINVI)+1.5*ECARRC)
         ELSE IF( DYRECT(NRLGLU) .GT. 0 ) THEN
            YRECT(NRLGSA) =INT(YRECT(NRLGLU)+DYRECT(NRLGLU)+1.5*ECARRC)
         ELSE IF( DYRECT(NRMENU) .GT. 0 ) THEN
            YRECT(NRLGSA) =INT(YRECT(NRMENU)+DYRECT(NRMENU)+1.5*ECARRC)
         ELSE IF( DYRECT(NRERR) .GT. 0 ) THEN
            YRECT(NRLGSA) =INT(YRECT(NRERR) +DYRECT(NRERR) +1.5*ECARRC)
         ENDIF
      ENDIF
C
C     LE TRACE DU TEXTE
      NA(1) = 0
      CALL RECTTX( NRLGSA, KLGSA, 0, NA )
C
      RETURN
      END
