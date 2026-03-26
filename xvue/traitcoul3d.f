      SUBROUTINE TRAITCOUL3D( XYZ, COUL )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACE LE TRAIT ENTRE 2 SOMMETS DE COULEURS COUL
C -----    SELON LES COULEURS INTERMEDIAIRES  (PALETTE 11 RECOMMANDEE)
C          ATTENTION XYZ EST EN COORDONNEES OBJETS 3D
C
C ENTREES:
C --------
C XYZ    : 3 COORDONNEES DES 2 SOMMETS DU TRAIT
C COUL   : NUMERO REEL DE LA COULEUR AUX 2 SOMMETS DU TRAIT
C          REEL DE VALEURS COMPRISES ENTRE N1COUL ET NDCOUL
C          (cf ~/incl/trvari.inc)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR:ALAIN PERRONNET LJLL UPMC PARIS & St Pierre du Perray JUIN 2009
C23456---------------------------------------------------------------012
      include"./incl/trvari.inc"
      include"./incl/minint.inc"
      REAL           XYZ(3,2), COUL(2)
      REAL           AXYZ(3), XYDIF(2), XY1(2), XY2(2)
      INTEGER        XYPX(2,2)
C
ccc      print *
ccc      print *,'XYZ=',XYZ,'  COUL=',COUL
C
C     TRANSFORMATION EN PIXELS DANS LA FENETRE XV
      DO 10 I=1,2
C        TRANSFORMATION EN COORDONNEES AXONOMETRIQUES
         CALL XYZAXO( XYZ(1,I), AXYZ )
C        SI UN POINT EST EXTERIEUR AUX 2 PLANS LA FACE N'EST PAS TRACEE
         IF( AXOARR .NE. 0 .OR. AXOAVA .NE. 0 ) THEN
C           AXOARR ET AXOAVA SONT ACTIFS
            IF( AXYZ(3) .LT. AXOARR ) RETURN
            IF( AXYZ(3) .GT. AXOAVA ) RETURN
         ENDIF
C        TRANSFORMATION EN PIXELS DANS LA FENETRE XV
         NX = NUPXEX( AXYZ(1) )
         NY = NUPXEY( AXYZ(2) )
C        SI LE NUMERO PIXEL EST INCORRECT ABANDON DU TRACE DE LA FACE
         IF( NX .EQ. MININT .OR. NY .EQ. MININT ) RETURN
         XYPX(1,I) = NX
         XYPX(2,I) = NY
 10   CONTINUE
C
C     LES 2 COULEURS EXTREMITES
      NC1 = INT( COUL(1) )
      NC2 = INT( COUL(2) )
ccc      print *,'XYPX=',XYPX,'  COUL=',NC1,NC2
      IF( NC1 .EQ. NC2 ) THEN
C
C        TRACE EFFECTIF DU TRAIT AVEC LA COULEUR COMMUNE
         CALL XVCOULEUR( NC1 )
         CALL XVTRAIT( XYPX(1,1), XYPX(2,1), XYPX(1,2), XYPX(2,2) )
C
      ELSE
C
C        TRACE PAR TRAIT DES POINTS AYANT MEME COULEUR
         XYDIF(1) = (XYPX(1,2)-XYPX(1,1)) / ( COUL(2)-COUL(1) )
         XYDIF(2) = (XYPX(2,2)-XYPX(2,1)) / ( COUL(2)-COUL(1) )
C
         IF( NC1 .LE. NC2 ) THEN
            INCR = 1
         ELSE
            INCR = -1
         ENDIF
C
C        TRACE EFFECTIF DU PREMIER TRAIT DE COULEUR NC1
         CALL XVCOULEUR( NC1 )
         XY1(1) = XYPX(1,1) + XYDIF(1) * ( NC1-COUL(1) )
         XY1(2) = XYPX(2,1) + XYDIF(2) * ( NC1-COUL(1) )
         CALL XVTRAIT( XYPX(1,1),    XYPX(2,1),
     %                 NINT(XY1(1)), NINT(XY1(2)) )
ccc      print *,'trait',nc1,':',(XYPX(K,1),k=1,2),' ->',(XY1(K),k=1,2)
C
C        TRACE DES TRAITS INTERMEDIAIRES DE COULEUR NC
         DO 20 NC = NC1+INCR, NC2-INCR, INCR
C
C           POINT DE FIN DE LA COULEUR NC
            XY2(1) = XYPX(1,1) + XYDIF(1) * ( NC-COUL(1) )
            XY2(2) = XYPX(2,1) + XYDIF(2) * ( NC-COUL(1) )
C           TRACE EFFECTIF
            CALL XVCOULEUR( NC )
            CALL XVTRAIT( NINT(XY1(1)), NINT(XY1(2)),
     %                    NINT(XY2(1)), NINT(XY2(2)) )
ccc      print *,'trait',nc,':',(XY1(k),k=1,2),' ->',(XY2(K),k=1,2)
            XY1(1) = XY2(1)
            XY1(2) = XY2(2)
C
 20      CONTINUE
C
C        TRACE EFFECTIF DU DERNIER TRAIT DE COULEUR NC2
         CALL XVCOULEUR( NC2 )
         CALL XVTRAIT( NINT(XY1(1)), NINT(XY1(2)),
     %                 XYPX(1,2),    XYPX(2,2) )
ccc      print *,'trait',nc2,':',(XY1(k),k=1,2),' ->',(XYPX(K,2),k=1,2)
C
      ENDIF
C
      RETURN
      END
