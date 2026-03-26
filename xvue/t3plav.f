      SUBROUTINE T3PLAV
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  TRACER DES 3 FACES DU TRIEDRE LE PLUS ELOIGNE DE L'OEIL
C -----  POUR TERMINER, IL FAUT CALL T3PLAP POUR TRACER LES 3 ARETES VUES
C        LES PARAMETRES DE TRACE SONT DANS $MEFISTO/incl/mecoit.inc
C        ET LE COMMON / T3PLAN /
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS         DECEMBRE 1997
C2345X7..............................................................012
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      include"./incl/xyzext.inc"
      COMMON / T3PLAN / XYZTRI(3,2), H(3), NDIV(3), NOFA3P(3), NOET3P
C
      REAL              XYZ1(3), XYZ2(3), XYZQUA(3,4)
C
C     PAS DE TRACE SI NTY3PL EST NEGATIF OU NUL
      NOET3P = 0
      IF( NTY3PL .LE. 0 ) RETURN
C
C     LONGUEUR MAXIMALE D'UN COTE DE L'HEXAEDRE ENGLOBANT
      DISMAX = MAX( COOEXT(1,2) - COOEXT(1,1),
     %              COOEXT(2,2) - COOEXT(2,1),
     %              COOEXT(3,2) - COOEXT(3,1) ) / 10.0
C
C     NOMBRE DE SUDIVISIONS DES 3 COTES
      DO 5 K=1,3
         X12     = COOEXT(K,2) - COOEXT(K,1)
         NDIV(K) = NINT( X12 / DISMAX ) + 2 * MAR3PL
         IF( NDIV(K) .LE. 2 ) NDIV(K) = 3
         H(K)    = X12 / ( NDIV(K) - 2 * MAR3PL )
 5    CONTINUE
C
C     L'HEXAEDRE DE TRACE
      XYZTRI(1,1) = COOEXT(1,1) - H(1) * MAR3PL
      XYZTRI(2,1) = COOEXT(2,1) - H(2) * MAR3PL
      XYZTRI(3,1) = COOEXT(3,1) - H(3) * MAR3PL
      XYZTRI(1,2) = COOEXT(1,2) + H(1) * MAR3PL
      XYZTRI(2,2) = COOEXT(2,2) + H(2) * MAR3PL
      XYZTRI(3,2) = COOEXT(3,2) + H(3) * MAR3PL
C
C     COORDONNEES DES BARYCENTRES DES FACES DE L'HEXAEDRE
      X12 = ( XYZTRI(1,1) + XYZTRI(1,2) ) * 0.5
      Y12 = ( XYZTRI(2,1) + XYZTRI(2,2) ) * 0.5
      Z12 = ( XYZTRI(3,1) + XYZTRI(3,2) ) * 0.5
C
C     TRACE DE LA FACE LA PLUS ELOIGNEE DE L'OEIL PARMI LES 2 FACES OPPOSEES
C     ======================================================================
C     L'EPAISSEUR DES TRAITS
      IF( NEP3PL .GE. 0 ) CALL XVEPAISSEUR( NEP3PL )
C
C     FACES ARRIERE ET DEVANT
C     -----------------------
C     LE BARYCENTRE DE LA FACE ARRIERE DE L'HEXAEDRE
      XYZ1(1) = XYZTRI(1,1)
      XYZ1(2) = Y12
      XYZ1(3) = Z12
C     LE BARYCENTRE DE LA FACE DEVANT DE L'HEXAEDRE
      XYZ2(1) = XYZTRI(1,2)
      XYZ2(2) = Y12
      XYZ2(3) = Z12
C
C     LES POINTS DANS L'AXONOMETRIE
      CALL XYZAXO( XYZ1, XYZ1 )
      CALL XYZAXO( XYZ2, XYZ2 )
      IF( XYZ1(3) .LE. XYZ2(3) ) THEN
C        TRACE DE LA FACE MIN
         L = 1
      ELSE
C        TRACE DE LA FACE MAX
         L = 2
      ENDIF
C     LE NUMERO DE LA FACE ARRIERE OU DEVANT TRACEE
      NOFA3P(1) = L
C     L'ABSCISSE FIXE
      X = XYZTRI(1,L)
C
C     TRACE SELON LE TYPE DEMANDE
      IF( NTY3PL .EQ. 1 ) THEN
C
C        8 ARETES ENGLOBANTES => LES 4 ARETES DE LA FACE A TRACER
         XYZ1(1) = X
         XYZ1(2) = XYZTRI(2,1)
         XYZ1(3) = XYZTRI(3,1)
         XYZ2(1) = X
         XYZ2(2) = XYZTRI(2,2)
         XYZ2(3) = XYZTRI(3,1)
         CALL TRAIT3D( NC13PL, XYZ1, XYZ2 )
C
         XYZ1(1) = X
         XYZ1(2) = XYZTRI(2,2)
         XYZ1(3) = XYZTRI(3,2)
         CALL TRAIT3D( NC13PL, XYZ2, XYZ1 )
C
         XYZ2(1) = X
         XYZ2(2) = XYZTRI(2,1)
         XYZ2(3) = XYZTRI(3,2)
         CALL TRAIT3D( NC13PL, XYZ1, XYZ2 )
C
         XYZ1(1) = X
         XYZ1(2) = XYZTRI(2,1)
         XYZ1(3) = XYZTRI(3,1)
         CALL TRAIT3D( NC13PL, XYZ2, XYZ1 )
C
      ELSE IF( NTY3PL .EQ. 2 ) THEN
C
C        TRACE SIMPLE DE LA GRILLE SANS REMPLISSAGE EN COULEUR NC13PL
         DO 10 K=0,NDIV(2)
C           LE POINT INITIAL
            XYZ1(1) = X
            XYZ1(2) = XYZTRI(2,1) + K * H(2)
            XYZ1(3) = XYZTRI(3,1)
C           LE POINT FINAL
            XYZ2(1) = X
            XYZ2(2) = XYZ1(2)
            XYZ2(3) = XYZTRI(3,2)
C           LE TRACE DU TRAIT
            CALL TRAIT3D( NC13PL, XYZ1, XYZ2 )
 10      CONTINUE
C
         DO 12 K=0,NDIV(3)
C           LE POINT INITIAL
            XYZ1(1) = X
            XYZ1(2) = XYZTRI(2,1)
            XYZ1(3) = XYZTRI(3,1) + K * H(3)
C           LE POINT FINAL
            XYZ2(1) = X
            XYZ2(2) = XYZTRI(2,2)
            XYZ2(3) = XYZ1(3)
C           LE TRACE DU TRAIT
            CALL TRAIT3D( NC13PL, XYZ1, XYZ2 )
 12      CONTINUE
C
      ELSE IF( NTY3PL .EQ. 3 .OR. NTY3PL .EQ. 4 ) THEN
C
C        TRACE D'UN DAMIER A 2 COULEURS NC23PL ET NC33PL ET DE GRILLE NC13PL
C        LE FOND
         XYZQUA(1,1) = X
         XYZQUA(2,1) = XYZTRI(2,1)
         XYZQUA(3,1) = XYZTRI(3,1)
C
         XYZQUA(1,2) = X
         XYZQUA(2,2) = XYZTRI(2,2)
         XYZQUA(3,2) = XYZTRI(3,1)
C
         XYZQUA(1,3) = X
         XYZQUA(2,3) = XYZTRI(2,2)
         XYZQUA(3,3) = XYZTRI(3,2)
C
         XYZQUA(1,4) = X
         XYZQUA(2,4) = XYZTRI(2,1)
         XYZQUA(3,4) = XYZTRI(3,2)
         CALL FACE3D( NC23PL, NC13PL, 4, XYZQUA )
C
         IF( NTY3PL .EQ. 4 ) THEN
C
C           TRACE DES CARRES DE COULEUR NC33PL DU DAMIER
            DO 16 L=0,NDIV(3)-1
               DO 14 K=MOD(L,2),NDIV(2)-1,2
C                 LE POINT INITIAL DU CARRE K,L
                  Y = XYZTRI(2,1) + K * H(2)
                  Z = XYZTRI(3,1) + L * H(3)
C
                  XYZQUA(1,1) = X
                  XYZQUA(2,1) = Y
                  XYZQUA(3,1) = Z
C
                  XYZQUA(1,2) = X
                  XYZQUA(2,2) = Y + H(2)
                  XYZQUA(3,2) = Z
C
                  XYZQUA(1,3) = X
                  XYZQUA(2,3) = Y + H(2)
                  XYZQUA(3,3) = Z + H(3)
C
                  XYZQUA(1,4) = X
                  XYZQUA(2,4) = Y
                  XYZQUA(3,4) = Z + H(3)
                  CALL FACE3D( NC33PL, NC13PL, 4, XYZQUA )
 14            CONTINUE
 16         CONTINUE
         ENDIF
      ENDIF
C
C     FACES GAUCHE ET DROITE
C     ----------------------
C     LE BARYCENTRE DE LA FACE GAUCHE DE L'HEXAEDRE
      XYZ1(1) = X12
      XYZ1(2) = XYZTRI(2,1)
      XYZ1(3) = Z12
C     LE BARYCENTRE DE LA FACE DROITE DE L'HEXAEDRE
      XYZ2(1) = X12
      XYZ2(2) = XYZTRI(2,2)
      XYZ2(3) = Z12
C
C     LES POINTS DANS L'AXONOMETRIE
      CALL XYZAXO( XYZ1, XYZ1 )
      CALL XYZAXO( XYZ2, XYZ2 )
      IF( XYZ1(3) .LE. XYZ2(3) ) THEN
C        TRACE DE LA FACE MIN
         L = 1
      ELSE
C        TRACE DE LA FACE MAX
         L = 2
      ENDIF
C     LE NUMERO DE LA FACE GAUCHE OU DROITE TRACEE
      NOFA3P(2) = L
C     L'ORDONNEE FIXE
      Y = XYZTRI(2,L)
C
C     TRACE SELON LE TYPE DEMANDE
      IF( NTY3PL .EQ. 1 ) THEN
C
C        8 ARETES ENGLOBANTES => LES 4 ARETES DE LA FACE A TRACER
         XYZ1(1) = XYZTRI(1,1)
         XYZ1(2) = Y
         XYZ1(3) = XYZTRI(3,1)
         XYZ2(1) = XYZTRI(1,2)
         XYZ2(2) = Y
         XYZ2(3) = XYZTRI(3,1)
         CALL TRAIT3D( NC13PL, XYZ1, XYZ2 )
C
         XYZ1(1) = XYZTRI(1,2)
         XYZ1(2) = Y
         XYZ1(3) = XYZTRI(3,2)
         CALL TRAIT3D( NC13PL, XYZ2, XYZ1 )
C
         XYZ2(1) = XYZTRI(1,1)
         XYZ2(2) = Y
         XYZ2(3) = XYZTRI(3,2)
         CALL TRAIT3D( NC13PL, XYZ1, XYZ2 )
C
         XYZ1(1) = XYZTRI(1,1)
         XYZ1(2) = Y
         XYZ1(3) = XYZTRI(3,1)
         CALL TRAIT3D( NC13PL, XYZ2, XYZ1 )
C
      ELSE IF( NTY3PL .EQ. 2 ) THEN
C
C        TRACE SIMPLE DE LA GRILLE SANS REMPLISSAGE EN COULEUR NC13PL
         DO 20 K=0,NDIV(1)
C           LE POINT INITIAL
            XYZ1(1) = XYZTRI(1,1) + K * H(1)
            XYZ1(2) = Y
            XYZ1(3) = XYZTRI(3,1)
C           LE POINT FINAL
            XYZ2(1) = XYZ1(1)
            XYZ2(2) = Y
            XYZ2(3) = XYZTRI(3,2)
C           LE TRACE DU TRAIT
            CALL TRAIT3D( NC13PL, XYZ1, XYZ2 )
 20      CONTINUE
C
         DO 22 K=0,NDIV(3)
C           LE POINT INITIAL
            XYZ1(1) = XYZTRI(1,1)
            XYZ1(2) = Y
            XYZ1(3) = XYZTRI(3,1) + K * H(3)
C           LE POINT FINAL
            XYZ2(1) = XYZTRI(1,2)
            XYZ2(2) = Y
            XYZ2(3) = XYZ1(3)
C           LE TRACE DU TRAIT
            CALL TRAIT3D( NC13PL, XYZ1, XYZ2 )
 22      CONTINUE
C
      ELSE IF( NTY3PL .EQ. 3 .OR. NTY3PL .EQ. 4 ) THEN
C
C        TRACE D'UN DAMIER A 2 COULEURS NC23PL ET NC33PL ET DE GRILLE NC13PL
C        LE FOND
         XYZQUA(1,1) = XYZTRI(1,1)
         XYZQUA(2,1) = Y
         XYZQUA(3,1) = XYZTRI(3,1)
C
         XYZQUA(1,2) = XYZTRI(1,2)
         XYZQUA(2,2) = Y
         XYZQUA(3,2) = XYZTRI(3,1)
C
         XYZQUA(1,3) = XYZTRI(1,2)
         XYZQUA(2,3) = Y
         XYZQUA(3,3) = XYZTRI(3,2)
C
         XYZQUA(1,4) = XYZTRI(1,1)
         XYZQUA(2,4) = Y
         XYZQUA(3,4) = XYZTRI(3,2)
         CALL FACE3D( NC23PL, NC13PL, 4, XYZQUA )
C
         IF( NTY3PL .EQ. 4 ) THEN
C
C           TRACE DES CARRES DE COULEUR NC33PL DU DAMIER
            DO 26 L=0,NDIV(3)-1
               DO 24 K=MOD(L,2),NDIV(1)-1,2
C                 LE POINT INITIAL DU CARRE K,L
                  X = XYZTRI(1,1) + K * H(1)
                  Z = XYZTRI(3,1) + L * H(3)
C
                  XYZQUA(1,1) = X
                  XYZQUA(2,1) = Y
                  XYZQUA(3,1) = Z
C
                  XYZQUA(1,2) = X + H(1)
                  XYZQUA(2,2) = Y
                  XYZQUA(3,2) = Z
C
                  XYZQUA(1,3) = X + H(1)
                  XYZQUA(2,3) = Y
                  XYZQUA(3,3) = Z + H(3)
C
                  XYZQUA(1,4) = X
                  XYZQUA(2,4) = Y
                  XYZQUA(3,4) = Z + H(3)
                  CALL FACE3D( NC33PL, NC13PL, 4, XYZQUA )
 24            CONTINUE
 26         CONTINUE
         ENDIF
      ENDIF
C
C     FACES BASSE ET HAUTE
C     --------------------
C     LE BARYCENTRE DE LA FACE BASSE DE L'HEXAEDRE
      XYZ1(1) = X12
      XYZ1(2) = Y12
      XYZ1(3) = XYZTRI(3,1)
C     LE BARYCENTRE DE LA FACE HAUTE DE L'HEXAEDRE
      XYZ2(1) = X12
      XYZ2(2) = Y12
      XYZ2(3) = XYZTRI(3,2)
C
C     LES POINTS DANS L'AXONOMETRIE
      CALL XYZAXO( XYZ1, XYZ1 )
      CALL XYZAXO( XYZ2, XYZ2 )
      IF( XYZ1(3) .LE. XYZ2(3) ) THEN
C        TRACE DE LA FACE MIN
         L = 1
      ELSE
C        TRACE DE LA FACE MAX
         L = 2
      ENDIF
C     LE NUMERO DE LA FACE BASSE OU HAUTE TRACEE
      NOFA3P(3) = L
C     LA COTE FIXE
      Z = XYZTRI(3,L)
C
C     TRACE SELON LE TYPE DEMANDE
      IF( NTY3PL .EQ. 1 ) THEN
C
C        8 ARETES ENGLOBANTES => LES 4 ARETES DE LA FACE A TRACER
         XYZ1(1) = XYZTRI(1,1)
         XYZ1(2) = XYZTRI(2,1)
         XYZ1(3) = Z
         XYZ2(1) = XYZTRI(1,2)
         XYZ2(2) = XYZTRI(2,1)
         XYZ2(3) = Z
         CALL TRAIT3D( NC13PL, XYZ1, XYZ2 )
C
         XYZ1(1) = XYZTRI(1,2)
         XYZ1(2) = XYZTRI(2,2)
         XYZ1(3) = Z
         CALL TRAIT3D( NC13PL, XYZ2, XYZ1 )
C
         XYZ2(1) = XYZTRI(1,1)
         XYZ2(2) = XYZTRI(2,2)
         XYZ2(3) = Z
         CALL TRAIT3D( NC13PL, XYZ1, XYZ2 )
C
         XYZ1(1) = XYZTRI(1,1)
         XYZ1(2) = XYZTRI(2,1)
         XYZ1(3) = Z
         CALL TRAIT3D( NC13PL, XYZ2, XYZ1 )
C
      ELSE IF( NTY3PL .EQ. 2 ) THEN
C
C        TRACE SIMPLE DE LA GRILLE SANS REMPLISSAGE EN COULEUR NC13PL
         DO 30 K=0,NDIV(1)
C           LE POINT INITIAL
            XYZ1(1) = XYZTRI(1,1) + K * H(1)
            XYZ1(2) = XYZTRI(2,1)
            XYZ1(3) = Z
C           LE POINT FINAL
            XYZ2(1) = XYZ1(1)
            XYZ2(2) = XYZTRI(2,2)
            XYZ2(3) = Z
C           LE TRACE DU TRAIT
            CALL TRAIT3D( NC13PL, XYZ1, XYZ2 )
 30      CONTINUE
C
         DO 32 K=0,NDIV(2)
C           LE POINT INITIAL
            XYZ1(1) = XYZTRI(1,1)
            XYZ1(2) = XYZTRI(2,1) + K * H(2)
            XYZ1(3) = Z
C           LE POINT FINAL
            XYZ2(1) = XYZTRI(1,2)
            XYZ2(2) = XYZ1(2)
            XYZ2(3) = Z
C           LE TRACE DU TRAIT
            CALL TRAIT3D( NC13PL, XYZ1, XYZ2 )
 32      CONTINUE
C
      ELSE IF( NTY3PL .EQ. 3 .OR. NTY3PL .EQ. 4 ) THEN
C
C        TRACE D'UN DAMIER A 2 COULEURS NC23PL ET NC33PL ET DE GRILLE NC13PL
C        LE FOND
         XYZQUA(1,1) = XYZTRI(1,1)
         XYZQUA(2,1) = XYZTRI(2,1)
         XYZQUA(3,1) = Z
C
         XYZQUA(1,2) = XYZTRI(1,2)
         XYZQUA(2,2) = XYZTRI(2,1)
         XYZQUA(3,2) = Z
C
         XYZQUA(1,3) = XYZTRI(1,2)
         XYZQUA(2,3) = XYZTRI(2,2)
         XYZQUA(3,3) = Z
C
         XYZQUA(1,4) = XYZTRI(1,1)
         XYZQUA(2,4) = XYZTRI(2,2)
         XYZQUA(3,4) = Z
         CALL FACE3D( NC23PL, NC13PL, 4, XYZQUA )
C
         IF( NTY3PL .EQ. 4 ) THEN
C
C           TRACE DES CARRES DE COULEUR NC33PL DU DAMIER
            DO 36 L=0,NDIV(2)-1
               DO 34 K=MOD(L,2),NDIV(1)-1,2
C                 LE POINT INITIAL DU CARRE K,L
                  X = XYZTRI(1,1) + K * H(1)
                  Y = XYZTRI(2,1) + L * H(2)
C
                  XYZQUA(1,1) = X
                  XYZQUA(2,1) = Y
                  XYZQUA(3,1) = Z
C
                  XYZQUA(1,2) = X + H(1)
                  XYZQUA(2,2) = Y
                  XYZQUA(3,2) = Z
C
                  XYZQUA(1,3) = X + H(1)
                  XYZQUA(2,3) = Y + H(2)
                  XYZQUA(3,3) = Z
C
                  XYZQUA(1,4) = X
                  XYZQUA(2,4) = Y + H(2)
                  XYZQUA(3,4) = Z
                  CALL FACE3D( NC33PL, NC13PL, 4, XYZQUA )
 34            CONTINUE
 36         CONTINUE
         ENDIF
      ENDIF
C
C     RETOUR A L'EPAISSEUR STANDARD DES TRAITS
      CALL XVEPAISSEUR( 0 )
C
C     ETAT DE SORTIE CORRECT LES TABLEAUX DU COMMON/T3PLAN/ SONT INITIALISES
      NOET3P = 1
      RETURN
      END
