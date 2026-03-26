      SUBROUTINE BSPLAC( NBCOMP ,
     %                   DEGREX , LRX , LUX , RX , UX , NBAXQB,
     %                   DEGREY , LRY , LUY , RY , UY , NBAYQB ,
     %                   SPLINE , XYZSOM,
     %                   NBTGS  , XYZTGS , XYZTST )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LES 3 COORDONNEES DES (NBAXQB+1)*(NBAYQB+1) SOMMETS
C -----    DES QUADRANGLES D'UNE SURFACE B-SPLINE D'INTERPOLATION
C
C          RAISON GEOMETRIQUE = 1 EN X ET Y
C          NOMBRE D'ARETES MULTIPLE DU NOMBRE D'INTERVALLES DES POINTS
C          D'INTERPOLATION EN X ET Y
C
C          SOUS CES CONDITIONS RESPECT DES POINTS D'INTERPOLATION COMME
C          SOMMETS DU MAILLAGE
C
C ENTREES:
C --------
C NBCOMP : NOMBRE DE COMPOSANTES DE SPLINE ( 3 POLYNOME, 4 RATIONNELLE )
C DEGREX : DEGRE DES POLYNOMES DE LA B-SPLINE
C LRX    : NOMBRE DE POINTS-1 DE CONTROLE (R) DE LA LIGNE
C LUX    : NOMBRE DE NOEUDS-1 DES POINTS D'INTERPOLATION EN X
C RX     : LES EXTREMITES DES INTERVALLES PARAMETRE EN X
C UX     : LES VALEURS DES NOEUDS EN X
C NBAXQB : NOMBRE D'ARETES EN X
C DEGREY : DEGRE DES POLYNOMES DE LA B-SPLINE
C LRY    : NOMBRE DE POINTS-1 DE CONTROLE (R) DE LA LIGNE
C LUY    : NOMBRE DE NOEUDS-1 DES POINTS D'INTERPOLATION EN Y
C RY     : LES EXTREMITES DES INTERVALLES PARAMETRE EN Y
C UY     : LES VALEURS DES NOEUDS EN Y
C NBAYQB : NOMBRE D'ARETES EN Y
C SPLINE : LES COEFFICIENTS DES POLYNOMES SUR CHAQUE INTERVALLE
C NBTGS  : NOMBRE DES TANGENTES A CALCULER
C
C SORTIES:
C --------
C XYZSOM : 3 COORDONNEES DES SOMMETS DES QUADRANGLES
C XYZTGS : 3 COORDONNEES DES 8 TANGENTES DES QUADRANGLES
C
C TABLEAU AUXILIAIRE:
C -------------------
C XYZTST : 3 COORDONNEES DES 2 TANGENTES EN CHACUN DES SOMMETS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS         MAI 1996
C2345X7..............................................................012
      INTEGER  DEGREX,DEGREY
      REAL     RX(0:LRX),RY(0:LRY),UX(0:LUX),UY(0:LUY),
     %         SPLINE(0:DEGREX,0:DEGREY,0:LRX-1,0:LRY-1,1:NBCOMP),
     %         XYZSOM(1:3,0:NBAXQB,0:NBAYQB)
      REAL     XYZTGS(1:3,1:8,1:NBAXQB,1:NBAYQB),
     %         XYZTST(1:3,1:2,0:NBAXQB,0:NBAYQB)
C
C     LE NOMBRE DE SOUS ARETES ENTRE POINTS D'INTERPOLATION
      NBAX = NBAXQB / (LUX-1)
      NBAY = NBAYQB / (LUY-1)
C
C     PASSAGE AUX 3 COORDONNEES DES SOMMETS DE LA LIGNE DANS R**3
C     ===========================================================
      NUY  = 0
      IY   = 0
      Y    = UY(0)
      DO 200 NY=0,NBAYQB
C
C        RECHERCHE DU [RY(IY),RY(IY-1)[ CONTENANT Y
 30      IF( Y .GE. RY(IY+1) ) THEN
C           PASSAGE A L'INTERVALLE SUIVANT DE RY
            IY = IY + 1
            IF( IY .LT. LRY ) GOTO 30
C           LE DERNIER SOMMET EST CALCULE PAR PROLONGEMENT
            IY = LRY - 1
         ENDIF
         RRY = Y - RY(IY)
C
         NUX  = 0
         IX   = 0
         X    = UX(0)
         DO 190 NX=0,NBAXQB
C
C           RECHERCHE DU [RX(IX),RX(IX-1)[ CONTENANT X
 50         IF( X .GE. RX(IX+1) ) THEN
C              PASSAGE A L'INTERVALLE SUIVANT DE RX
               IX = IX + 1
               IF( IX .LT. LRX ) GOTO 50
C              LE DERNIER SOMMET EST CALCULE PAR PROLONGEMENT
               IX = LRX - 1
            ENDIF
            RRX = X - RX(IX)
C
            DO 80 J=1,3
               AA = 0
               DO 70 M=DEGREY,0,-1
                  A = SPLINE(DEGREX,M,IX,IY,J)
                  DO 60 L=DEGREX-1,0,-1
C                    DOUBLE METHODE DE HORNER
                     A = A * RRX + SPLINE(L,M,IX,IY,J)
 60               CONTINUE
                  AA = AA * RRY + A
 70            CONTINUE
               XYZSOM(J,NX,NY) = AA
 80         CONTINUE
C
            IF( NBTGS .GT. 0 ) THEN
C
C              CALCUL DE LA DERIVEE/PARAMETRE1 EN CE SOMMET
               DO 130 J=1,3
                  AA = 0
                  DO 120 M=DEGREY,0,-1
                     A = DEGREX * SPLINE(DEGREX,M,IX,IY,J)
                     DO 110 L=DEGREX-1,1,-1
C                       DOUBLE METHODE DE HORNER
                        A = A * RRX + L * SPLINE(L,M,IX,IY,J)
 110                 CONTINUE
                     AA = AA * RRY + A
 120              CONTINUE
                  XYZTST(J,1,NX,NY) = AA
 130           CONTINUE
C
C              CALCUL DE LA DERIVEE/PARAMETRE2 EN CE SOMMET
               DO 160 J=1,3
                  AA = 0
                  DO 150 L=DEGREX,0,-1
                     A = DEGREY * SPLINE(L,DEGREY,IX,IY,J)
                     DO 140 M=DEGREY-1,1,-1
C                       DOUBLE METHODE DE HORNER
                        A = A * RRY + M * SPLINE(L,M,IX,IY,J)
 140                 CONTINUE
                     AA = AA * RRX + A
 150              CONTINUE
                  XYZTST(J,2,NX,NY) = AA
 160           CONTINUE
            ENDIF
C
C           LA VALEUR DU PARAMETRE SUIVANT EN X
            IF( MOD(NX,NBAX) .EQ. 0 .AND. NUX .LT. LUX ) NUX = NUX + 1
            X = X + ( UX(NUX)-UX(NUX-1) ) / NBAX
 190     CONTINUE
C
C        LA VALEUR DU PARAMETRE SUIVANT EN Y
         IF( MOD(NY,NBAY) .EQ. 0 .AND. NUY .LT. LUY ) NUY = NUY + 1
         Y = Y + ( UY(NUY)-UY(NUY-1) ) / NBAY
 200  CONTINUE
C
C     REPARTITION DES TANGENTES DANS CHAQUE QUADRANGLE
C     ================================================
      IF( NBTGS .LE. 0 ) RETURN
      NUY  = 1
      Y    = UY(0)
      PASY = ( UY(1)-UY(0) ) / NBAY
C
      DO 500 NY=1,NBAYQB
C
C        LE PAS EN Y EST POUR L'ARETE NY EST PASY
         NUX  = 1
         X    = UX(0)
         PASX = ( UX(1) - UX(0) ) / NBAX
C
         DO 490 NX=1,NBAXQB
C           LE PAS EN X EST POUR L'ARETE NX EST PASX
C           LES 8 TANGENTES DU QUADRANGLE (NX,NY)
            DO 450 J=1,3
C              LA COORDONNEE J DES 8 TANGENTES DU QUADRANGLE (NX,NY)
C              LA TANGENTE SELON LE PARAMETRE X AU SOMMET 1 DU QUADRANGLE
               XYZTGS(J,1,NX,NY) =  XYZTST(J,1,NX-1,NY-1) * PASX
C              LA TANGENTE SELON LE PARAMETRE Y AU SOMMET 1 DU QUADRANGLE
               XYZTGS(J,2,NX,NY) =  XYZTST(J,2,NX-1,NY-1) * PASY
C
C              LA TANGENTE SELON LE PARAMETRE X AU SOMMET 2 DU QUADRANGLE
               XYZTGS(J,3,NX,NY) =  XYZTST(J,2,NX,NY-1) * PASY
C              LA TANGENTE SELON LE PARAMETRE Y AU SOMMET 2 DU QUADRANGLE
               XYZTGS(J,4,NX,NY) = -XYZTST(J,1,NX,NY-1) * PASX
C
C              LA TANGENTE SELON LE PARAMETRE X AU SOMMET 3 DU QUADRANGLE
               XYZTGS(J,5,NX,NY) = -XYZTST(J,1,NX,NY) * PASX
C              LA TANGENTE SELON LE PARAMETRE Y AU SOMMET 3 DU QUADRANGLE
               XYZTGS(J,6,NX,NY) = -XYZTST(J,2,NX,NY) * PASY
C
C              LA TANGENTE SELON LE PARAMETRE X AU SOMMET 4 DU QUADRANGLE
               XYZTGS(J,7,NX,NY) = -XYZTST(J,2,NX-1,NY) * PASY
C              LA TANGENTE SELON LE PARAMETRE Y AU SOMMET 4 DU QUADRANGLE
               XYZTGS(J,8,NX,NY) =  XYZTST(J,1,NX-1,NY) * PASX
 450        CONTINUE
C
C           LA VALEUR DU PARAMETRE SUIVANT EN X
            IF( MOD(NX,NBAX) .EQ. 0 .AND. NUX .LT. LUX ) THEN
               NUX  = NUX + 1
               PASX = ( UX(NUX)-UX(NUX-1) ) / NBAX
            ENDIF
            X = X + PASX
 490     CONTINUE
C
C        LA VALEUR DU PARAMETRE SUIVANT EN Y
         IF( MOD(NY,NBAY) .EQ. 0 .AND. NUY .LT. LUY ) THEN
            NUY  = NUY + 1
            PASY = ( UY(NUY)-UY(NUY-1) ) / NBAY
         ENDIF
         Y = Y + PASY
 500  CONTINUE
      END
