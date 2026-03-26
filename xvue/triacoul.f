      SUBROUTINE TRIACOUL( XYPX, COUL )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  TRACE LE TRIANGLE DE SOMMETS PIXELS XYPX ET DE COULEURS COUL
C -----  SELON LES COULEURS INTERMEDIAIRES  (PALETTE 11 RECOMMANDEE)
C        ATTENTION (X,Y) EN COORDONNEES PIXELS INTEGER*2
C                  COUL  EST UN TABLEAU DE REELS DE VALEURS COMPRISES
C                        ENTRE N1COUL ET NDCOUL (cf ~/incl/trvari.inc)
C
C ENTREES:
C --------
C XYPX   : 2 COORDONNEES DES 3 SOMMETS
C COUL   : NUMERO DES 3 COULEURS AUX 3 SOMMETS DU TRIANGLE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC  PARIS    OCTOBRE 1994
C2345X7..............................................................012
      INTEGER*2  XYPX(1:2,1:3)
      REAL       COUL(1:3)
      INTRINSIC  INT2
C
      INTEGER*2  XY(2,3), XYF(2,5)
C
      IF(COUL(1) .LT. 0 .OR. COUL(2) .LT. 0 .OR. COUL(3) .LT. 0) RETURN
C
C     RECHERCHE DE LA COULEUR MINIMALE C1 ET MAXIMALE C3
      I1 = 1
      I3 = 1
      C1 = COUL(1)
      C3 = COUL(1)
      IF( COUL(2) .LT. C1 ) THEN
         C1 = COUL(2)
         I1 = 2
      ELSE IF( COUL(2) .GT. C3 ) THEN
         C3 = COUL(2)
         I3 = 2
      ENDIF
      IF( COUL(3) .LT. C1 ) THEN
         C1 = COUL(3)
         I1 = 3
      ELSE IF( COUL(3) .GT. C3 ) THEN
         C3 = COUL(3)
         I3 = 3
      ENDIF
C
C     PARTIE ENTIERE DE LA COULEUR
      NC1 = INT( C1 )
      NC3 = INT( C3 )
C
      IF( NC1 .EQ. NC3 ) THEN
C
C        LES 3 COULEURS A TRACER SONT EGALES -> TRIANGLE MONO-COULEUR
         CALL XVCOULEUR( NC1 )
         CALL XVFACE( 3, XYPX )
         RETURN
      ENDIF
C
C     RECHERCHE DE LA COULEUR INTERMEDIAIRE
      DO 10 I2=1,3
        IF( I2 .NE. I1 .AND. I2 .NE. I3 ) GOTO 15
 10   CONTINUE
 15   C2  = COUL(I2)
      NC2 = INT( C2 )
C
C     LES 3 SOMMETS SONT REORDONNES SELON C1<=C2<=C3
      DO 20 NC=1,2
         XYF(NC,1) = XYPX(NC,I1)
         XY (NC,1) = XYPX(NC,I1)
         XY (NC,2) = XYPX(NC,I2)
         XY (NC,3) = XYPX(NC,I3)
 20   CONTINUE
C
C     RECHERCHE DE L'ARETE P13-P12 OU LA COULEUR NC1 DEVIENT NC1+1
      NCP1 = NC1 + 1
      IF( NC1 .NE. NC2 ) THEN
         C21 = 1.0 / (C2-C1)
         XYF(1,2) = INT2( NINT( (XY(1,1)*(C2-NCP1)
     %                          +XY(1,2)*(NCP1-C1))*C21 ) )
         XYF(2,2) = INT2( NINT( (XY(2,1)*(C2-NCP1)
     %                          +XY(2,2)*(NCP1-C1))*C21 ) )
         C31 = 1.0 / (C3-C1)
         XYF(1,3) = INT2( NINT( (XY(1,1)*(C3-NCP1)
     %                          +XY(1,3)*(NCP1-C1))*C31 ) )
         XYF(2,3) = INT2( NINT( (XY(2,1)*(C3-NCP1)
     %                          +XY(2,3)*(NCP1-C1))*C31 ) )
C        TRACE DU TRIANGLE SOMMET 1 - P12 - P13
         CALL XVCOULEUR( NC1 )
         CALL XVFACE( 3, XYF )
C
C        LES QUADRANGLES SUIVANTS
         XYF(1,1) = XYF(1,3)
         XYF(2,1) = XYF(2,3)
         DO 30 NC=NC1+1, NC2-1
            NCP1 = NC + 1
            XYF(1,3) = INT2( NINT( (XY(1,1)*(C2-NCP1)
     %                             +XY(1,2)*(NCP1-C1))*C21 ) )
            XYF(2,3) = INT2( NINT( (XY(2,1)*(C2-NCP1)
     %                             +XY(2,2)*(NCP1-C1))*C21 ) )
            XYF(1,4) = INT2( NINT( (XY(1,1)*(C3-NCP1)
     %                             +XY(1,3)*(NCP1-C1))*C31 ) )
            XYF(2,4) = INT2( NINT( (XY(2,1)*(C3-NCP1)
     %                             +XY(2,3)*(NCP1-C1))*C31 ) )
C           TRACE DU QUADRANGLE P13-P12-P12-P13
            CALL XVCOULEUR( NC )
            CALL XVFACE( 4, XYF )
C           3->2 ET 4->1
            XYF(1,2) = XYF(1,3)
            XYF(2,2) = XYF(2,3)
            XYF(1,1) = XYF(1,4)
            XYF(2,1) = XYF(2,4)
 30      CONTINUE
C
         IF( NC2 .EQ. NC3 ) THEN
C
C           IL RESTE LE QUADRANGLE P13-P12-S2-S3
            XYF(1,3) = XY(1,2)
            XYF(2,3) = XY(2,2)
            XYF(1,4) = XY(1,3)
            XYF(2,4) = XY(2,3)
C           TRACE DU QUADRANGLE P13-P12-S2-S3
            CALL XVCOULEUR( NC2 )
            CALL XVFACE( 4, XYF )
            RETURN
         ENDIF
C
C        CAS NC1<NC2<NC3
C        LE PENTAGONE P13-P12-S2-P23-P13
         XYF(1,3) = XY(1,2)
         XYF(2,3) = XY(2,2)
         NCP1 = NC2 + 1
         C32  = 1.0 / (C3-C2)
         XYF(1,4) = INT2( NINT( (XY(1,2)*(C3-NCP1)
     %                          +XY(1,3)*(NCP1-C2))*C32 ) )
         XYF(2,4) = INT2( NINT( (XY(2,2)*(C3-NCP1)
     %                          +XY(2,3)*(NCP1-C2))*C32 ) )
         XYF(1,5) = INT2( NINT( (XY(1,1)*(C3-NCP1)
     %                          +XY(1,3)*(NCP1-C1))*C31 ) )
         XYF(2,5) = INT2( NINT( (XY(2,1)*(C3-NCP1)
     %                          +XY(2,3)*(NCP1-C1))*C31 ) )
C        TRACE DU PENTAGONE P13-P12-SOMMET2-P23-P13
         CALL XVCOULEUR( NC2 )
         CALL XVFACE( 5, XYF )
C        4->2 ET 5->1
         XYF(1,2) = XYF(1,4)
         XYF(2,2) = XYF(2,4)
         XYF(1,1) = XYF(1,5)
         XYF(2,1) = XYF(2,5)
         GOTO 50
      ENDIF
C
C     CAS : NC1=NC2<NC3
C     LE QUADRANGLE SOMMET1-SOMMET2-P23-P13
C     XYF(1,1) = XY(1,1)   DEJA FAIT
C     XYF(2,1) = XY(2,1)
      XYF(1,2) = XY(1,2)
      XYF(2,2) = XY(2,2)
      NCP1 = NC2 + 1
      C31 = 1.0 / (C3-C1)
      C32 = 1.0 / (C3-C2)
      XYF(1,3) = INT2( NINT( (XY(1,2)*(C3-NCP1)
     %                       +XY(1,3)*(NCP1-C2))*C32 ) )
      XYF(2,3) = INT2( NINT( (XY(2,2)*(C3-NCP1)
     %                       +XY(2,3)*(NCP1-C2))*C32 ) )
      XYF(1,4) = INT2( NINT( (XY(1,1)*(C3-NCP1)
     %                       +XY(1,3)*(NCP1-C1))*C31 ) )
      XYF(2,4) = INT2( NINT( (XY(2,1)*(C3-NCP1)
     %                       +XY(2,3)*(NCP1-C1))*C31 ) )
C     TRACE DU QUADRANGLE SOMMET1-SOMMET2-P23-P13
      CALL XVCOULEUR( NC2 )
      CALL XVFACE( 4, XYF )
C     3->2 ET 4->1
      XYF(1,2) = XYF(1,3)
      XYF(2,2) = XYF(2,3)
      XYF(1,1) = XYF(1,4)
      XYF(2,1) = XYF(2,4)
C
C     LES QUADRANGLES P13-P23-P23-P13
 50   DO 60 NC=NC2+1, NC3-1
         NCP1 = NC + 1
         XYF(1,3) = INT2( NINT( (XY(1,2)*(C3-NCP1)
     %                          +XY(1,3)*(NCP1-C2))*C32 ) )
         XYF(2,3) = INT2( NINT( (XY(2,2)*(C3-NCP1)
     %                          +XY(2,3)*(NCP1-C2))*C32 ) )
         XYF(1,4) = INT2( NINT( (XY(1,1)*(C3-NCP1)
     %                          +XY(1,3)*(NCP1-C1))*C31 ) )
         XYF(2,4) = INT2( NINT( (XY(2,1)*(C3-NCP1)
     %                          +XY(2,3)*(NCP1-C1))*C31 ) )
C        TRACE DU QUADRANGLE P13-P23-P23-P13
         CALL XVCOULEUR( NC )
         CALL XVFACE( 4, XYF )
C        3->2 ET 4->1
         XYF(1,2) = XYF(1,3)
         XYF(2,2) = XYF(2,3)
         XYF(1,1) = XYF(1,4)
         XYF(2,1) = XYF(2,4)
 60   CONTINUE
C
C     LE TRIANGLE P13-P12-SOMMET 3
      XYF(1,3) = XY(1,3)
      XYF(2,3) = XY(2,3)
C     TRACE DU TRIANGLE P13-P23-SOMMET 3
      CALL XVCOULEUR( NC3 )
      CALL XVFACE( 3, XYF )
      RETURN
      END
