      SUBROUTINE T2FLEC( NOCOUL, X, Y, DISTCM, DIRXCM, DIRYCM )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : TRACE UNE FLECHE PARTANT DE (X,Y) DANS LA DIRECTION (DIRX,DIRY)
C ----- DE LONGUEUR DISTCM DE TRIANGLE EQUILATERAL DANS LE SENS USUEL
C       SI DISTCM > 0 , DANS L AUTRE SENS SI DISTCM < 0
C
C ENTREES :
C ---------
C NOCOUL : NUMERO DE LA COULEUR DE LA FLECHE A TRACER
C X, Y   : COORDONNEES OBJET DU POINT DE DEPART DE LA FLECHE SUR
C          LA SURFACE DE VISUALISATION
C DISTCM : LONGUEUR DE LA FLECHE
C          SON SIGNE DETERMINE LE SENS DU TRIANGLE DE LA FLECHE
C DIRXCM : ABCISSE  EN CM DE LA DIRECTION DE LA FLECHE
C DIRYCM : ORDONNEE EN CM DE LA DIRECTION DE LA FLECHE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS         MAI 1999
C23456---------------------------------------------------------------012
      include "./incl/mecoit.inc"
      include "./incl/trvari.inc"
      REAL      XY(2,4)
C
C     SI LA TAILLE DE LA FLECHE EST < 0.1 CM ELLE N EST PAS TRACEE
      D = ABS( DISTCM )
      IF( D .LT. 0.1 ) RETURN
C
C     LE NOMBRE DE CM DE LA LARGEUR DE LA FENETRE
      R = LAPXFE * 0.1 / CXMMPX
C
C     1 CM DE LARGEUR SUR L'ECRAN REPRESENTE R UNITES OBJETS
      R = ( XOBMAX - XOBMIN ) / R
C
C     MISE A L'ECHELLE DES 2 COMPOSANTES
      D  = D * R
      XD = DIRXCM * R
      YD = DIRYCM * R
C
C     VECTEUR NORMAL A LA FLECHE
      XN = - YD / D
      YN =   XD / D
C
C     FIN DE LA FLECHE
      XY(1,1) = X + XD
      XY(2,1) = Y + YD
C
C     TRACE DE LA FLECHE SANS LE TRIANGLE TERMINAL
      CALL XVEPAISSEUR( NEPFLE )
      CALL TRAIT2D( NOCOUL, X, Y, XY(1,1), XY(2,1) )
C
C     DEMI-LONGUEUR DU COTE DU TRIANGLE EQUILATERAL
      D = D * 0.06
      IF( NEPFLE .GT. 1 ) D = D * NEPFLE * 0.75
C
C     POINT INTERSECTION FLECHE ET TRIANGLE
      XY(1,4) = X + XD * 0.75
      XY(2,4) = Y + YD * 0.75
C
C     LE TRIANGLE TERMINAL EST AIGUISE
      D = D * 0.75

      IF( DISTCM .GT. 0 ) THEN
C
C        TRACTION ---|> EN NOIR
         XY(1,2) = XY(1,4) + XN * D
         XY(2,2) = XY(2,4) + YN * D
C
         XY(1,3) = XY(1,4) - XN * D
         XY(2,3) = XY(2,4) - YN * D
C
C        TRACE DU TRIANGLE TERMINAL
ccc         CALL TRAIT2D( NCNOIR, XY(1,1), XY(2,1), XY(1,2), XY(2,2) )
ccc         CALL TRAIT2D( NCNOIR, XY(1,2), XY(2,2), XY(1,3), XY(2,3) )
ccc         CALL TRAIT2D( NCNOIR, XY(1,3), XY(2,3), XY(1,1), XY(2,1) )
         CALL TRIA2D( XY, NCNOIR )
C
      ELSE
C
C        COMPRESSION ---<| EN ROUGE
         XY(1,2) = XY(1,1) + XN * D
         XY(2,2) = XY(2,1) + YN * D
C
         XY(1,3) = XY(1,1) - XN * D
         XY(2,3) = XY(2,1) - YN * D
C
         XY(1,1) = XY(1,4)
         XY(2,1) = XY(2,4)
C
C        TRACE DU TRIANGLE TERMINAL
ccc         CALL TRAIT2D( NCROUG, XY(1,2), XY(2,2), XY(1,3), XY(2,3) )
ccc         CALL TRAIT2D( NCROUG, XY(1,3), XY(2,3), XY(1,4), XY(2,4) )
ccc         CALL TRAIT2D( NCROUG, XY(1,4), XY(2,4), XY(1,2), XY(2,2) )
         CALL TRIA2D( XY, NCROUG )
C
      ENDIF
C
      RETURN
      END
