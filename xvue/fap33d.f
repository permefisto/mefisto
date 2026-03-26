      SUBROUTINE FAP33D( NCF, NCA, POREDF, NBS, X, Y, Z, NOSTEF,ROSTEF )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  TRACE LA FACE DEGRE 3 3D DE SOMMETS ET TANGENTES RANGEES DANS
C -----  X,Y,Z  (SOMMETS 1 A NBS, TANGENTES NBS+1 A 3*NBS) ET
C        SELON LES COULEURS NCF POUR LA FACE ET NCA POUR LES ARETES
C
C        ATTENTION (X,Y,Z) EN COORDONNEES OBJET 3D
C        LA TRANSFORMATION EN PIXELS EST ASSUREE DANS CE SP
C        LES AUTRES CARACTERISTIQUES DU TRACE SONT CELLES ACTUELLES
C
C ENTREES:
C --------
C NCF    : NUMERO DE LA COULEUR DE REMPLISSAGE DE LA FACE A TRACER
C          <0 PAS DE TRACE DE L'INTERIEUR DE LA FACE
C NCA    : NUMERO DE LA COULEUR DU CONTOUR     DE LA FACE A TRACER
C          <0 PAS DE TRACE DES ARETES DU CONTOUR DE LA FACE
C POREDF : POURCENTAGE DE REDUCTION DE LA FACE ( 0<=POREDF<100 )
C NBS    : NOMBRE DE SOMMETS ( 3 ou 4 ) DE LA FACE A TRACER
C X      : ABSCISSE OBJET DES NBS SOMMETS ET TANGENTES DE LA FACE
C Y      : ORDONNEE OBJET DES NBS SOMMETS ET TANGENTES DE LA FACE
C Z      : COTE     OBJET DES NBS SOMMETS ET TANGENTES DE LA FACE
C NOSTEF : NUMERO DES SOMMETS ET DES TANGENTES DE LA FACE
C          RANGEES SUIVANT L'ORDRE
C          TRIANGLE  : NO  SOMMET1, NS2 , NS3,
C                      NO TANGENTE1(S1S2), NT2(S1S3),
C                      NO TANGENTE3(S2S3), NT4(S2S1),
C                      NO TANGENTE5(S3S1), NT6(S3S2)
C          QUADRANGLE: NO  SOMMET1, NS2 , NS3, NS4,
C                      NO TANGENTE1(S1S2), NT2(S1S3),   NT3(S2S3), NT4(S2S1),
C                      NO TANGENTE5(S3S4), NT6(S3S2),   NT7(S4S1), NT8(S4S3)
C          CE CHOIX PERMET UNE BOUCLE SUR LES TANGENTES PAR LES SOMMETS
C ROSTEF : IDEM NOSTEF A L'APPELet EN SORTIE (en EQUIVALENCE)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC  PARIS    OCTOBRE 1996
C2345X7..............................................................012
      include"./incl/mxsuaf.inc"
      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)
      REAL              X(12),Y(12),Z(12), XYZ(3), AXYZ(3),
     %                  CNORFA(3), COUL(3)
      REAL              RCOULT(MXSUBT), RCOULQ(0:MXSUAF,0:MXSUAF)
      EQUIVALENCE      (RCOULT,RCOULQ)
      INTEGER           NOSTEF(1:12)
      REAL              ROSTEF(1:12)
      INTEGER*2         XPX(0:MXSUAF,0:MXSUAF),
     %                  YPX(0:MXSUAF,0:MXSUAF)
      INTEGER*2         XTPX(1:MXSUBT),
     %                  YTPX(1:MXSUBT)
      EQUIVALENCE      (XPX,XTPX), (YPX,YTPX)
      INTEGER*2         XYPX(2,5)
      REAL              PROFON(4)
      INTRINSIC         INT2
      include"./incl/trvari.inc"
      include"./incl/minint.inc"
      include"./incl/mecoit.inc"
C
C     FONCTION INTERNE
      NUSTRI(I,J) = I+1 + J*(J+1)/2
C
C     SI COULEUR NEGATIVE => RETOUR
      IF( NCF .LT. 0 .AND. NCA .LT. 0 ) RETURN
C
C     VERIFICATION
      IF( NBS .LE. 2 .OR. NBS. GT. 4 ) RETURN
C
C     LA REDUCTION DE LA FACE ENTRE 0 ET 1
C     0.005 = 0.01/3 CAR TRACE ICI DANS L'ELEMENT DE REFERENCE
C     ET NON PAR RAPPROCHEMENT DU BARYCENTRE COMME DANS LE CAS D'EF P1
      REDUCF = POREDF * 0.00333
      REDUC1 = 1.0 - REDUCF
C
C     RECHERCHE DU NOMBRE DE SUBDIVISIONS DE L'EF REDUIT DE REFERENCE
C     SELON LE NOMBRE DE PIXELS OBTENUS APRES MISE A L'ECHELLE PIXELS
C     RECHERCHE DES 2 SOMMETS DE LA FACE LES PLUS ELOIGNES DE L'OEIL
C     ===============================================================
      NBPX = 0
C
      DO 5 I=1,NBS
C        TRANSFORMATION EN COORDONNEES AXONOMETRIQUES
         XYZ(1) = X(I)
         XYZ(2) = Y(I)
         XYZ(3) = Z(I)
         CALL XYZAXO( XYZ, AXYZ )
C        SI UN POINT EST EXTERIEUR AUX 2 PLANS LA FACE N'EST PAS TRACEE
         IF( AXOARR .NE. 0 .OR. AXOAVA .NE. 0 ) THEN
C           AXOARR ET AXOAVA SONT ACTIFS
            IF( AXYZ(3) .LT. AXOARR ) RETURN
            IF( AXYZ(3) .GT. AXOAVA ) RETURN
         ENDIF
C        LA DISTANCE A L'OEIL DU SOMMET I
         PROFON(I) = AXYZ(3)
C        TRANSFORMATION EN PIXELS DANS LA FENETRE XV
         NX = NUPXEX( AXYZ(1) )
         IF( NX .LE. MININT ) RETURN
C        SI LE NUMERO PIXEL EST INCORRECT ABANDON DU TRACE DE LA FACE
         NY = NUPXEY( AXYZ(2) )
         IF( NY .LE. MININT ) RETURN
         XYPX(1,I) = INT2( NX )
         XYPX(2,I) = INT2( NY )
         IF( I .GT. 1 ) THEN
            I1   = ABS( XYPX(1,I)-XYPX(1,I-1) )
            I2   = ABS( XYPX(2,I)-XYPX(2,I-1) )
            NBPX = MAX( NBPX, I1, I2 )
         ENDIF
C
C        LES COORDONNEES PIXELS DU SOMMET I + PREMIERE TANGENTE
         I1  = NBS + 2 * I - 1
         XYZ(1) = X(I) + X(I1)
         XYZ(2) = Y(I) + Y(I1)
         XYZ(3) = Z(I) + Z(I1)
         CALL XYZAXO( XYZ, AXYZ )
         NX1 = NUPXEX( AXYZ(1) )
         NY1 = NUPXEY( AXYZ(2) )
         NBPX = MAX( NBPX, ABS(NX1-NX) , ABS(NY1-NY) )
 5    CONTINUE
      IF( NBPX .EQ. 0 ) RETURN
      NBSUBA = MIN( MAX(2,NBPX/8), NBSUAF )
C
C     CALCUL DES IMAGES P3 DES SOMMETS DU MAILLAGE DE LA FACE REDUITE
C     TRACE DES SOUS-EF, IMAGE P3 DES SOUS-EF DE L'EF DE REFERENCE
C     ===============================================================
      H = ( REDUC1 - REDUCF ) / NBSUBA
C
      IF( NBS .EQ. 3 ) THEN
C
C        TRIANGLE HCT-REDUIT
C        ===================
C        RECHERCHE DU SOMMET DU TRIANGLE LE PLUS ELOIGNE DE L'OEIL
C        POUR TRAITER L'ALGORITHME DU PEINTRE
         P = MIN( PROFON(1), PROFON(2), PROFON(3) )
         DO 6 K=1,3
            IF( P .EQ. PROFON(K) ) GOTO 7
 6       CONTINUE
 7       IF( K .EQ. 2 ) THEN
C           PERMUTATION DES SOMMETS POUR AMENER CE SOMMET 2 EN POSITION 1
            CALL PECI33( ROSTEF )
            CALL PECI33( X )
            CALL PECI33( Y )
            CALL PECI33( Z )
         ELSE IF( K .EQ. 3 ) THEN
C           PERMUTATION DES SOMMETS POUR AMENER CE SOMMET 3 EN POSITION 1
            CALL PECM33( ROSTEF )
            CALL PECM33( X )
            CALL PECM33( Y )
            CALL PECM33( Z )
         ENDIF
C
C        LES COORDONNEES X Y Z DES SOMMETS DES SOUS TRIANGLES DU TRIANGLE UNITE
         V = REDUCF
         K = 0
         DO 20 J=0,NBSUBA
            U = REDUCF + J * H
            V = REDUCF
            DO 10 I=0,J
               K  = K + 1
C              LES 3 COORDONNES HCT AU POINT (U,V)
               CALL XYZHCT( U, V, X, Y, Z,
     %                      XYZ(1), XYZ(2), XYZ(3) )
C              LES COORDONNEES EN PIXELS SUR L'ECRAN
               CALL XYZAXO( XYZ(1), AXYZ )
               XTPX(K) = INT2( NUPXEX( AXYZ(1) ) )
               YTPX(K) = INT2( NUPXEY( AXYZ(2) ) )
C
               IF( LCRITR .EQ. 0 ) THEN
C                 LES 3 COORDONNES HCT DE LA NORMALE AU POINT (U,V)
                  CALL NORHCT( U, V, X, Y, Z,
     %                         CNORFA, IERR )
C                 CALCUL DE LA COULEUR EN CE POINT A PARTIR
C                 DE LA NORMALE ET DES ECLAIRAGES
                  CALL ECLAIR( CNORFA, RCOULT(K) )
               ENDIF
C
C              PASSAGE AU POINT SUIVANT
               U = U - H
               V = V + H
 10         CONTINUE
 20      CONTINUE
C
C        TRACE EVENTUEL DES 2 TANGENTES DU SOMMET 1 LE PLUS ELOIGNE
         IF( IAVTGF .NE. 0 ) THEN
C
C           LA COULEUR DE TRACE POINTILLE DES TANGENTES
            CALL XVCOULEUR( NCOTGF )
C           TRACE EN  CONTINU(0) ou POINTILLE (1) ou DOUBLE POINTILLE(2)
            CALL XVTYPETRAIT( NTRTGF )
            CALL XVEPAISSEUR( 0 )
C
C           AU SOMMET 1
            XYZ(1) = X(1)
            XYZ(2) = Y(1)
            XYZ(3) = Z(1)
C           LES COORDONNEES EN PIXELS SUR L'ECRAN
            CALL XYZAXO( XYZ, AXYZ )
            NX1 = NUPXEX( AXYZ(1) )
            NY1 = NUPXEY( AXYZ(2) )
C
            DO 21 J=1,2
               IF( NOSTEF(3+J) .NE. 0 ) THEN
C                 LA J-EME TANGENTE AU SOMMET 1
                  XYZ(1) = X(1) + X(3+J)
                  XYZ(2) = Y(1) + Y(3+J)
                  XYZ(3) = Z(1) + Z(3+J)
C                 LES COORDONNEES EN PIXELS SUR L'ECRAN DE LA POINTE DE LA TG
                  CALL XYZAXO( XYZ, AXYZ )
                  NX2 = NUPXEX( AXYZ(1) )
                  NY2 = NUPXEY( AXYZ(2) )
C                 LE TRACE DE LA TANGENTE
                  CALL XVTRAIT( NX1, NY1,  NX2, NY2 )
               ENDIF
 21         CONTINUE
            NUPRTG = 2
         ENDIF
C
C        TRACE DES SOUS TRIANGLES ET DES SOUS ARETES DES 3 COTES
         DO 40 J=1,NBSUBA
            DO 30 I=1,J
C              LES COORDONNEES PIXELS DES 3 SOMMETS DU SOUS-TRIANGLE
               NS1 = NUSTRI(I-1,J-1)
               NS2 = NUSTRI(I-1,J  )
               NS3 = NUSTRI(I  ,J  )
C
               XYPX(1,1) = XTPX( NS1 )
               XYPX(1,2) = XTPX( NS2 )
               XYPX(1,3) = XTPX( NS3 )
               XYPX(1,4) = XTPX( NS1 )
C
               XYPX(2,1) = YTPX( NS1 )
               XYPX(2,2) = YTPX( NS2 )
               XYPX(2,3) = YTPX( NS3 )
               XYPX(2,4) = YTPX( NS1 )
C
C              TRACE DE LA SOUS-FACE TRIANGULAIRE
               IF( NCF .GE. 0 ) THEN
                  IF( LCRITR .EQ. 0 ) THEN
C                    TRACE SELON LA COULEUR AUX 3 SOMMETS DU TRIANGLE
                     COUL(1) = RCOULT( NS1 )
                     COUL(2) = RCOULT( NS2 )
                     COUL(3) = RCOULT( NS3 )
                     CALL TRIACOUL( XYPX, COUL )
                  ELSE
C                    TRACE DU TRIANGLE SUIVANT LA COULEUR NCF
                     CALL XVCOULEUR( NCF )
                     CALL XVFACE( 4, XYPX )
                  ENDIF
               ENDIF
C
C              TRACE EVENTUEL DE L'ARETE DE LA SOUS FACE
               IF( NCA .GE. 0 ) THEN
                  IF( I .EQ. 1 ) THEN
C                    L'ARETE 1 DU TRIANGLE TRACE PLEIN AVEC NEPARF EPAISSEURS
                     CALL XVCOULEUR( NCA )
C                    TRACE EN  CONTINU(0) ou POINTILLE (1) ou DOUBLE POINTILLE(2
                     CALL XVTYPETRAIT( NTRTGF )
                     CALL XVEPAISSEUR( NEPARF )
                     NX1 = XYPX(1,1)
                     NY1 = XYPX(2,1)
                     NX2 = XYPX(1,2)
                     NY2 = XYPX(2,2)
                     CALL XVTRAIT( NX1, NY1, NX2, NY2 )
                  ENDIF
                  IF( J .EQ. NBSUBA ) THEN
C                    L'ARETE 2 DU TRIANGLE TRACE PLEIN AVEC NEPARF EPAISSEURS
                     CALL XVCOULEUR( NCA )
                     CALL XVTYPETRAIT( 0 )
                     CALL XVEPAISSEUR( NEPARF )
                     NX1 = XYPX(1,2)
                     NY1 = XYPX(2,2)
                     NX2 = XYPX(1,3)
                     NY2 = XYPX(2,3)
                     CALL XVTRAIT( NX1, NY1, NX2, NY2 )
                  ENDIF
                  IF( I .EQ. J ) THEN
C                    L'ARETE 2 DU TRIANGLE TRACE PLEIN AVEC NEPARF EPAISSEURS
                     CALL XVCOULEUR( NCA )
                     CALL XVTYPETRAIT( 0 )
                     CALL XVEPAISSEUR( NEPARF )
                     NX1 = XYPX(1,3)
                     NY1 = XYPX(2,3)
                     NX2 = XYPX(1,1)
                     NY2 = XYPX(2,1)
                     CALL XVTRAIT( NX1, NY1, NX2, NY2 )
                  ENDIF
               ENDIF
C
C              LES COORDONNEES PIXELS DES 3 SOMMETS DU SOUS-TRIANGLE
               IF( I .NE. J ) THEN
                  NS2 = NUSTRI(I,J-1)
                  XYPX(1,1) = XTPX( NS1 )
                  XYPX(1,2) = XTPX( NS2 )
                  XYPX(1,3) = XTPX( NS3 )
                  XYPX(1,4) = XTPX( NS1 )
C
                  XYPX(2,1) = YTPX( NS1 )
                  XYPX(2,2) = YTPX( NS2 )
                  XYPX(2,3) = YTPX( NS3 )
                  XYPX(2,4) = YTPX( NS1 )
C
C                 TRACE DE LA SOUS-FACE TRIANGULAIRE
                  IF( NCF .GE. 0 ) THEN
                     IF( LCRITR .EQ. 0 ) THEN
C                       TRACE SELON LA COULEUR AUX 3 SOMMETS DU TRIANGLE
                        COUL(1) = RCOULT( NS1 )
                        COUL(2) = RCOULT( NS2 )
                        COUL(3) = RCOULT( NS3 )
                        CALL TRIACOUL( XYPX, COUL )
                     ELSE
C                       TRACE DU TRIANGLE SUIVANT LA COULEUR NCF
                        CALL XVCOULEUR( NCF )
                        CALL XVFACE( 4, XYPX )
                     ENDIF
                  ENDIF
               ENDIF
 30         CONTINUE
 40      CONTINUE
CCCC
CCCC        TRACE EVENTUEL DES ARETES DU TRIANGLE HCTR
CCC         IF( NCA .LT. 0 ) RETURN
CCCC        LA COULEUR DE L'ARETE
CCC         CALL XVCOULEUR( NCA )
CCC         CALL XVEPAISSEUR( NEPARF )
CCCC
CCCC        LA PREMIERE ARETE
CCC         NS1 = 1
CCC         NX1 = XTPX( NS1 )
CCC         NY1 = YTPX( NS1 )
CCC         DO 50 J=1,NBSUBA
CCC            NS2 = NUSTRI( 0, J )
CCC            NX2 = XTPX( NS2 )
CCC            NY2 = YTPX( NS2 )
CCC            CALL XVTRAIT( NX1, NY1,  NX2, NY2 )
CCC            NX1 = NX2
CCC            NY1 = NY2
CCC 50      CONTINUE
CCCC
CCCC        LA SECONDE ARETE
CCC         DO 60 J=1,NBSUBA
CCC            NS2 = NUSTRI( J, NBSUBA )
CCC            NX2 = XTPX( NS2 )
CCC            NY2 = YTPX( NS2 )
CCC            CALL XVTRAIT( NX1, NY1,  NX2, NY2 )
CCC            NX1 = NX2
CCC            NY1 = NY2
CCC 60      CONTINUE
CCCC
CCCC        LA TROISIEME ARETE
CCC         DO 70 J=NBSUBA-1,0,-1
CCC            NS2 = NUSTRI( J, J )
CCC            NX2 = XTPX( NS2 )
CCC            NY2 = YTPX( NS2 )
CCC            CALL XVTRAIT( NX1, NY1,  NX2, NY2 )
CCC            NX1 = NX2
CCC            NY1 = NY2
CCC 70      CONTINUE
C
      ELSE
C
C        QUADRANGLE dVS-REDUIT
C        =====================
C        RECHERCHE DU SOMMET DU TRIANGLE LE PLUS ELOIGNE DE L'OEIL
C        POUR TRAITER L'ALGORITHME DU PEINTRE
         P = MIN( PROFON(1), PROFON(2), PROFON(3) , PROFON(4) )
         DO 106 K=1,4
            IF( P .EQ. PROFON(K) ) GOTO 107
 106     CONTINUE
 107     IF( K .EQ. 2 ) THEN
C           PERMUTATION DES SOMMETS POUR AMENER CE SOMMET 2 EN POSITION 1
            CALL PECI44( ROSTEF )
            CALL PECI44( X )
            CALL PECI44( Y )
            CALL PECI44( Z )
            CALL PECI4R( PROFON )
         ELSE IF( K .EQ. 4 ) THEN
C           PERMUTATION DES SOMMETS POUR AMENER CE SOMMET 4 EN POSITION 1
            CALL PECM44( ROSTEF )
            CALL PECM44( X )
            CALL PECM44( Y )
            CALL PECM44( Z )
            CALL PECM4R( PROFON )
         ELSE IF( K .EQ. 3 ) THEN
C           PERMUTATION DES SOMMETS POUR AMENER CE SOMMET 3 EN POSITION 1
            CALL PEC244( ROSTEF )
            CALL PEC244( X )
            CALL PEC244( Y )
            CALL PEC244( Z )
            CALL PEC24R( PROFON )
         ENDIF
         IF( PROFON(4) .LT. PROFON(2) ) THEN
C           INVERSION DU SENS DES 4 SOMMETS
C           AINSI L'ARETE LA PLUS ELOIGNEE EST L'ARETE 1 DU QUADRANGLE
C           CE QUI CONDUIT A TRACER SES SOUS QUADRANGLES ADJACENTS D'ABORD
            CALL PE44IN( ROSTEF )
            CALL PE44IN( X )
            CALL PE44IN( Y )
            CALL PE44IN( Z )
            CALL PE4INV( PROFON )
         ENDIF
C
C        LES COORDONNEES X Y DES SOMMETS DES SOUS CARRES DU CARRE UNITE
         V = REDUCF
         DO 120 J=0,NBSUBA
            U = REDUCF
            DO 110 I=0,NBSUBA
               CALL XYZDVS( U, V, X, Y, Z,
     %                      XYZ(1), XYZ(2), XYZ(3) )
C              LES COORDONNEES EN PIXELS SUR L'ECRAN
               CALL XYZAXO( XYZ, AXYZ )
               XPX(I,J) = INT2( NUPXEX( AXYZ(1) ) )
               YPX(I,J) = INT2( NUPXEY( AXYZ(2) ) )
C
               IF( LCRITR .EQ. 0 ) THEN
C                 LES 3 COORDONNES DE LA NORMALE AU POINT (U,V)
                  CALL NORDVS( U, V, X, Y, Z,
     %                         CNORFA, IERR )
C                 CALCUL DE LA COULEUR EN CE POINT A PARTIR
C                 DE LA NORMALE ET DES ECLAIRAGES
                  CALL ECLAIR( CNORFA, RCOULQ(I,J) )
               ENDIF
C
C              PASSAGE AU SOMMET SUIVANT
               U = U + H
 110        CONTINUE
            V = V + H
 120     CONTINUE
C
C        TRACE EVENTUEL DES TANGENTES AUX 2 SOMMETS DE L'ARETE LA PLUS ELOIGNEE
         IF( IAVTGF .NE. 0 ) THEN
C
C           LA COULEUR DE TRACE POINTILLE DES TANGENTES
            CALL XVCOULEUR( NCOTGF )
C           TRACE EN  CONTINU(0) ou POINTILLE (1) ou DOUBLE POINTILLE(2)
            CALL XVTYPETRAIT( NTRTGF )
            CALL XVEPAISSEUR( 0 )
            K = NBS
            DO 122 I=1,2
C              LE SOMMET I
               XYZ(1) = X(I)
               XYZ(2) = Y(I)
               XYZ(3) = Z(I)
C              LES COORDONNEES EN PIXELS SUR L'ECRAN
               CALL XYZAXO( XYZ, AXYZ )
               NX1 = NUPXEX( AXYZ(1) )
               NY1 = NUPXEY( AXYZ(2) )
C
               DO 121 J=1,2
C                 LA J-EME TANGENTE AU SOMMET I
                  K = K + 1
                  IF( NOSTEF(K) .NE. 0 ) THEN
                     XYZ(1) = X(I) + X(K)
                     XYZ(2) = Y(I) + Y(K)
                     XYZ(3) = Z(I) + Z(K)
C                    LES COORDONNEES EN PIXELS SUR L'ECRAN DE LA POINTE DE LA TG
                     CALL XYZAXO( XYZ, AXYZ )
                     NX2 = NUPXEX( AXYZ(1) )
                     NY2 = NUPXEY( AXYZ(2) )
C                    LE TRACE DE LA TANGENTE
                     CALL XVTRAIT( NX1, NY1,  NX2, NY2 )
                  ENDIF
 121           CONTINUE
 122        CONTINUE
C           LES 2 TANGENTES DU SECOND POINT ELOIGNE SERA RETRACE
            NUPRTG = 2
         ENDIF
C
         DO 140 J=1,NBSUBA
            DO 130 I=1,NBSUBA
C              LES COORDONNEES PIXELS DES 4 SOMMETS DE LA SOUS-FACE
               XYPX(1,1) = XPX( I-1, J-1 )
               XYPX(1,2) = XPX( I  , J-1 )
               XYPX(1,3) = XPX( I  , J   )
               XYPX(1,4) = XPX( I-1, J   )
               XYPX(1,5) = XYPX(1,1)
C
               XYPX(2,1) = YPX( I-1, J-1 )
               XYPX(2,2) = YPX( I  , J-1 )
               XYPX(2,3) = YPX( I  , J   )
               XYPX(2,4) = YPX( I-1, J   )
               XYPX(2,5) = XYPX(2,1)
C
C              TRACE DE LA SOUS-FACE QUADRANGULAIRE
               IF( NCF .GE. 0 ) THEN
                  IF( LCRITR .EQ. 0 ) THEN
C                    TRACE SELON LA COULEUR AUX 3 SOMMETS DU TRIANGLE 123
                     COUL(1) = RCOULQ( I-1, J-1 )
                     COUL(2) = RCOULQ( I  , J-1 )
                     COUL(3) = RCOULQ( I  , J   )
                     CALL TRIACOUL( XYPX, COUL )
C                    TRACE SELON LA COULEUR AUX 3 SOMMETS DU TRIANGLE 341
                     COUL(1) = RCOULQ( I  , J   )
                     COUL(2) = RCOULQ( I-1, J   )
                     COUL(3) = RCOULQ( I-1, J-1 )
                     CALL TRIACOUL( XYPX(1,3), COUL )
                  ELSE
C                    TRACE DU QUADRANGLE SUIVANT LA COULEUR NCF
                     CALL XVCOULEUR( NCF )
                     CALL XVFACE( 5, XYPX )
                  ENDIF
               ENDIF
C
C              TRACE EVENTUEL DE L'ARETE DE LA SOUS FACE
               IF( NCA .GE. 0 ) THEN
                  IF( J .EQ. 1 ) THEN
C                    L'ARETE 1 DU QUADRANGLE TRACE PLEIN AVEC NEPARF EPAISSEURS
                     CALL XVCOULEUR( NCA )
                     CALL XVTYPETRAIT( 0 )
                     CALL XVEPAISSEUR( NEPARF )
                     NX1 = XPX(I-1,0)
                     NY1 = YPX(I-1,0)
                     NX2 = XPX(I  ,0)
                     NY2 = YPX(I  ,0)
                     CALL XVTRAIT( NX1, NY1, NX2, NY2 )
                  ENDIF
                  IF( I .EQ. NBSUBA ) THEN
C                    L'ARETE 2 DU QUADRANGLE TRACE PLEIN AVEC NEPARF EPAISSEURS
                     CALL XVCOULEUR( NCA )
                     CALL XVTYPETRAIT( 0 )
                     CALL XVEPAISSEUR( NEPARF )
                     NX1 = XPX(NBSUBA,J-1)
                     NY1 = YPX(NBSUBA,J-1)
                     NX2 = XPX(NBSUBA,J  )
                     NY2 = YPX(NBSUBA,J  )
                     CALL XVTRAIT( NX1, NY1, NX2, NY2 )
                  ENDIF
                  IF( J .EQ. NBSUBA ) THEN
C                    L'ARETE 3 DU QUADRANGLE TRACE PLEIN AVEC NEPARF EPAISSEURS
                     CALL XVCOULEUR( NCA )
                     CALL XVTYPETRAIT( 0 )
                     CALL XVEPAISSEUR( NEPARF )
                     NX1 = XPX(I-1,NBSUBA)
                     NY1 = YPX(I-1,NBSUBA)
                     NX2 = XPX(I  ,NBSUBA)
                     NY2 = YPX(I  ,NBSUBA)
                     CALL XVTRAIT( NX1, NY1, NX2, NY2 )
                  ENDIF
                  IF( I .EQ. 1 ) THEN
C                    L'ARETE 4 DU QUADRANGLE TRACE PLEIN AVEC NEPARF EPAISSEURS
                     CALL XVCOULEUR( NCA )
                     CALL XVTYPETRAIT( 0 )
                     CALL XVEPAISSEUR( NEPARF )
                     NX1 = XPX(0,J-1)
                     NY1 = YPX(0,J-1)
                     NX2 = XPX(0,J  )
                     NY2 = YPX(0,J  )
                     CALL XVTRAIT( NX1, NY1, NX2, NY2 )
                  ENDIF
               ENDIF
 130        CONTINUE
 140     CONTINUE
CCCC
CCCC        TRACE EVENTUEL DES ARETES DU QUADRANGLE dVSR
CCC         IF( NCA .LT. 0 ) RETURN
CCCC        LA COULEUR DE L'ARETE
CCC         CALL XVCOULEUR( NCA )
CCCCCCC        TRACE DES TRAITS EN CONTINU
CCCCCC         CALL XVTYPETRAIT( 0 )
CCC         CALL XVEPAISSEUR( NEPARF )
CCCC
CCCC        LA PREMIERE ARETE
CCC         NX1 = XPX( 0, 0 )
CCC         NY1 = YPX( 0, 0 )
CCC         DO 150 I=1,NBSUBA
CCC            NX2 = XPX( I, 0 )
CCC            NY2 = YPX( I, 0 )
CCC            CALL XVTRAIT( NX1, NY1,  NX2, NY2 )
CCC            NX1 = NX2
CCC            NY1 = NY2
CCC 150     CONTINUE
CCCC
CCCC        LA SECONDE ARETE
CCC         DO 160 J=1,NBSUBA
CCC            NX2 = XPX( NBSUBA, J )
CCC            NY2 = YPX( NBSUBA, J )
CCC            CALL XVTRAIT( NX1, NY1,  NX2, NY2 )
CCC            NX1 = NX2
CCC            NY1 = NY2
CCC 160     CONTINUE
CCCC
CCCC        LA TROISIEME ARETE
CCC         DO 170 I=NBSUBA-1,0,-1
CCC            NX2 = XPX( I, NBSUBA )
CCC            NY2 = YPX( I, NBSUBA )
CCC            CALL XVTRAIT( NX1, NY1,  NX2, NY2 )
CCC            NX1 = NX2
CCC            NY1 = NY2
CCC 170     CONTINUE
CCCC
CCCC        LA QUATRIEME ARETE
CCC         DO 180 J=NBSUBA-1,0,-1
CCC            NX2 = XPX( 0, J )
CCC            NY2 = YPX( 0, J )
CCC            CALL XVTRAIT( NX1, NY1,  NX2, NY2 )
CCC            NX1 = NX2
CCC            NY1 = NY2
CCC 180     CONTINUE
CCC
      ENDIF
C
C     TRACE EVENTUEL DES TANGENTES AUX SOMMETS EN BLANC
      IF( IAVTGF .NE. 0 ) THEN
C
C        LA COULEUR DE TRACE POINTILLE DES TANGENTES
         CALL XVCOULEUR( NCOTGF )
C        TRACE EN  CONTINU(0) ou POINTILLE (1) ou DOUBLE POINTILLE(2)
         CALL XVTYPETRAIT( NTRTGF )
         CALL XVEPAISSEUR( 0 )
C
C        NUPRTG = NUMERO DU PROCHAIN SOMMET DE TGS A TRACER
         K = NBS + 2 * NUPRTG - 2
         DO 220 I=NUPRTG,NBS
C
C           LE SOMMET I
            XYZ(1) = X(I)
            XYZ(2) = Y(I)
            XYZ(3) = Z(I)
C           LES COORDONNEES EN PIXELS SUR L'ECRAN
            CALL XYZAXO( XYZ, AXYZ )
            NX1 = NUPXEX( AXYZ(1) )
            NY1 = NUPXEY( AXYZ(2) )
C
            DO 210 J=1,2
C
C              LA J-EME TANGENTE AU SOMMET I
               K = K + 1
               IF( NOSTEF(K) .NE. 0 ) THEN
                  XYZ(1) = X(I) + X(K)
                  XYZ(2) = Y(I) + Y(K)
                  XYZ(3) = Z(I) + Z(K)
C                 LES COORDONNEES EN PIXELS SUR L'ECRAN DE LA POINTE DE LA TG
                  CALL XYZAXO( XYZ, AXYZ )
                  NX2 = NUPXEX( AXYZ(1) )
                  NY2 = NUPXEY( AXYZ(2) )
C
C                 LE TRACE DE LA TANGENTE
                  CALL XVTRAIT( NX1, NY1,  NX2, NY2 )
               ENDIF
C
 210        CONTINUE
C
 220     CONTINUE
      ENDIF
C
C     AU DELA LES TRAITS SONT CONTINUS ET FINS
      CALL XVTYPETRAIT( 0 )
      CALL XVEPAISSEUR( 0 )
      END
