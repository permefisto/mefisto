      SUBROUTINE FAP32D( NCF, NCA, POREDF, NBS, X, Y, NOSOEL )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  TRACE LA FACE DEGRE 3 2D DE SOMMETS ET TANGENTES RANGEES DANS
C -----  X ET Y  (SOMMETS 1 A NBS, TANGENTES NBS+1 A 3*NBS) ET
C        SELON LES COULEURS NCF POUR LA FACE ET NCA POUR LES ARETES
C
C        ATTENTION (X,Y) EN COORDONNEES OBJET 2D
C        LA TRANSFORMATION EN PIXELS EST ASSUREE DANS CE SP
C        LES AUTRES CARACTERISTIQUES DU TRACE SONT CELLES ACTUELLES
C
C ENTREE :
C --------
C NCF    : NUMERO DE LA COULEUR DE REMPLISSAGE DE LA FACE A TRACER
C          <0 PAS DE TRACE DE L'INTERIEUR DE LA FACE
C NCA    : NUMERO DE LA COULEUR DU CONTOUR     DE LA FACE A TRACER
C          <0 PAS DE TRACE DES ARETES DU CONTOUR DE LA FACE
C POREDF : POURCENTAGE DE REDUCTION DE LA FACE ( 0<=POREDF<100 )
C NBS    : NOMBRE DE SOMMETS ( 3 ou 4 ) DE LA FACE A TRACER
C X      : ABSCISSE OBJET DES NBS SOMMETS ET TANGENTES DE LA FACE
C Y      : ORDONNEE OBJET DES NBS SOMMETS ET TANGENTES DE LA FACE
C NOSOEL : NUMERO DES SOMMETS ET DES EVENTUELLES TANGENTES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC  PARIS        MAI 1996
C2345X7..............................................................012
      include"./incl/mxsuaf.inc"
      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)
      REAL              X(12),Y(12)
      INTEGER           NOSOEL(12)
      INTRINSIC         INT2
C
      INTEGER*2         XPX(0:MXSUAF,0:MXSUAF),
     %                  YPX(0:MXSUAF,0:MXSUAF)
      INTEGER*2         XTPX(1:MXSUBT),
     %                  YTPX(1:MXSUBT)
      EQUIVALENCE      (XPX,XTPX), (YPX,YTPX)
      INTEGER*2         XYPX(2,5)
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
C     ===============================================================
      NBPX = 0
      NX1  = NUPXEX( X(NBS) )
      NY1  = NUPXEY( Y(NBS) )
      IF( NX1 .EQ. MININT .OR. NY1 .EQ. MININT ) RETURN
      DO 5 I=1,NBS
C        LES COORDONNEES PIXELS DU SOMMET I
         NX2 = NUPXEX( X(I) )
         NY2 = NUPXEY( Y(I) )
         IF( NX2 .EQ. MININT .OR. NY2 .EQ. MININT ) RETURN
         NBPX = MAX( NBPX, ABS(NX2-NX1) , ABS(NY2-NY1) )
C        LES COORDONNEES PIXELS DU SOMMET I + SECONDE TANGENTE
         I1  = NBS + 2*I
         NX3 = NUPXEX( X(I) + X(I1) )
         NY3 = NUPXEY( Y(I) + Y(I1) )
         NBPX = MAX( NBPX, ABS(NX3-NX2) , ABS(NY3-NY2) )
C        LES COORDONNEES PIXELS DU SOMMET I + PREMIERE TANGENTE
         I1  = I1-1
         NX3 = NUPXEX( X(I) + X(I1) )
         NY3 = NUPXEY( Y(I) + Y(I1) )
         NBPX = MAX( NBPX, ABS(NX3-NX2) , ABS(NY3-NY2) )
C        PASSAGE AU SUIVANT
         NX1 = NX2
         NY1 = NY2
 5    CONTINUE
      IF( NBPX .EQ. 0 ) RETURN
C     LE NOMBRE DE SUBDIVISIONS DE L'ARETE DU PARAMETRE
      NBSUBA = MIN( MAX(2,NBPX/8), NBSUAF )
C
C     CALCUL DES IMAGES P3 DES SOMMETS DU MAILLAGE DE LA FACE REDUITE
C     TRACE DES SOUS-EF, IMAGE P3 DES SOUS-EF DE L'EF DE REFERENCE
C     ===============================================================
      H = ( REDUC1 - 2*REDUCF ) / NBSUBA
C
      IF( NBS .EQ. 3 ) THEN
C
C        TRIANGLE HCT-REDUIT
C        ===================
C        LES COORDONNEES X Y DES SOMMETS DES SOUS TRIANGLES DU TRIANGLE UNITE
         V = REDUCF
         K = 0
         DO 20 J=0,NBSUBA
            U = REDUCF + J * H
            V = REDUCF
            DO 10 I=0,J
               K  = K + 1
C              LES 2 COORDONNES HCT AU POINT (U,V)
               CALL XYHCT( U, V, X, Y, XX, YY )
C              LES COORDONNEES EN PIXELS SUR L'ECRAN
               XTPX(K) = INT2( NUPXEX( XX ) )
               YTPX(K) = INT2( NUPXEY( YY ) )
C              PASSAGE AU POINT SUIVANT
               U = U - H
               V = V + H
 10         CONTINUE
 20      CONTINUE
C
         DO 40 J=1,NBSUBA
            DO 30 I=1,J
C              LES COORDONNEES PIXELS DES 3 SOMMETS DU SOUS-TRIANGLE
               NS1 = NUSTRI(I-1,J-1)
               NS2 = NUSTRI(I-1,J  )
               NS3 = NUSTRI(I  ,J  )
               XYPX(1,1) = XTPX( NS1 )
               XYPX(1,2) = XTPX( NS2 )
               XYPX(1,3) = XTPX( NS3 )
               XYPX(1,4) = XTPX( NS1 )
C
               XYPX(2,1) = YTPX( NS1 )
               XYPX(2,2) = YTPX( NS2 )
               XYPX(2,3) = YTPX( NS3 )
               XYPX(2,4) = YTPX( NS1 )
C              TRACE DES ARETES ET DE LA SOUS-FACE TRIANGULAIRE
               IF( NCF .GE. 0 ) THEN
                  CALL XVCOULEUR( NCF )
                  CALL XVFACE( 4, XYPX )
CCC               CALL XVFACETRAITS( NCF, NCBLAN, 4, XYPX )
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
C                 TRACE DES ARETES ET DE LA SOUS-FACE TRIANGULAIRE
                  IF( NCF .GE. 0 ) THEN
                     CALL XVCOULEUR( NCF )
                     CALL XVFACE( 4, XYPX )
CCC                  CALL XVFACETRAITS( NCF, NCBLAN, 4, XYPX )
                  ENDIF
               ENDIF
 30         CONTINUE
 40      CONTINUE
C
C        TRACE EVENTUEL DES ARETES DU TRIANGLE HCTR
         IF( NCA .LT. 0 ) RETURN
C        LA COULEUR DE L'ARETE
         CALL XVCOULEUR( NCA )
C        TRACE DES TRAITS EN CONTINU
         CALL XVTYPETRAIT( 0 )
         CALL XVEPAISSEUR( NEPARF )
C
C        LA PREMIERE ARETE
         NS1 = 1
         NX1 = XTPX( NS1 )
         NY1 = YTPX( NS1 )
         DO 50 J=1,NBSUBA
            NS2 = NUSTRI( 0, J )
            NX2 = XTPX( NS2 )
            NY2 = YTPX( NS2 )
            CALL XVTRAIT( NX1, NY1,  NX2, NY2 )
            NX1 = NX2
            NY1 = NY2
 50      CONTINUE
C
C        LA SECONDE ARETE
         DO 60 J=1,NBSUBA
            NS2 = NUSTRI( J, NBSUBA )
            NX2 = XTPX( NS2 )
            NY2 = YTPX( NS2 )
            CALL XVTRAIT( NX1, NY1,  NX2, NY2 )
            NX1 = NX2
            NY1 = NY2
 60      CONTINUE
C
C        LA TROISIEME ARETE
         DO 70 J=NBSUBA-1,0,-1
            NS2 = NUSTRI( J, J )
            NX2 = XTPX( NS2 )
            NY2 = YTPX( NS2 )
            CALL XVTRAIT( NX1, NY1,  NX2, NY2 )
            NX1 = NX2
            NY1 = NY2
 70      CONTINUE
C
      ELSE
C
C        QUADRANGLE dVS-REDUIT
C        =====================
C        LES COORDONNEES X Y DES SOMMETS DES SOUS CARRES DU CARRE UNITE
         V = REDUCF
         DO 120 J=0,NBSUBA
            U = REDUCF
            DO 110 I=0,NBSUBA
               CALL XYDVS( U, V, X, Y,  XX, YY )
C              LES COORDONNEES PIXELS SUR L'ECRAN
               XPX(I,J) = INT2( NUPXEX( XX ) )
               YPX(I,J) = INT2( NUPXEY( YY ) )
C              PASSAGE AU SOMMET SUIVANT
               U = U + H
 110        CONTINUE
            V = V + H
 120     CONTINUE
C
         DO 140 J=1,NBSUBA
            DO 130 I=1,NBSUBA
C              LES COORDONNEES PIXELS DES 4 SOMMETS DE LA SOUS-FACE
               XYPX(1,1) = XPX(I-1,J-1)
               XYPX(1,2) = XPX(I  ,J-1)
               XYPX(1,3) = XPX(I  ,J  )
               XYPX(1,4) = XPX(I-1,J  )
               XYPX(1,5) = XYPX(1,1)
C
               XYPX(2,1) = YPX(I-1,J-1)
               XYPX(2,2) = YPX(I  ,J-1)
               XYPX(2,3) = YPX(I  ,J  )
               XYPX(2,4) = YPX(I-1,J  )
               XYPX(2,5) = XYPX(2,1)
C
C              TRACE DES ARETES ET DE LA SOUS-FACE QUADRANGULAIRE
               IF( NCF .GE. 0 ) THEN
                  CALL XVCOULEUR( NCF )
                  CALL XVFACE( 4, XYPX )
CCC               CALL XVFACETRAITS( NCF, NCBLAN, 5, XYPX )
               ENDIF
 130        CONTINUE
 140     CONTINUE
C
C        TRACE EVENTUEL DES ARETES DU QUADRANGLE dVSR
         IF( NCA .LT. 0 ) RETURN
C        LA COULEUR DE L'ARETE
         CALL XVCOULEUR( NCA )
C        TRACE DES TRAITS EN CONTINU
         CALL XVTYPETRAIT( 0 )
         CALL XVEPAISSEUR( NEPARF )
C
C        LA PREMIERE ARETE
         NX1 = XPX( 0, 0 )
         NY1 = YPX( 0, 0 )
         DO 150 I=1,NBSUBA
            NX2 = XPX( I, 0 )
            NY2 = YPX( I, 0 )
            CALL XVTRAIT( NX1, NY1,  NX2, NY2 )
            NX1 = NX2
            NY1 = NY2
 150     CONTINUE
C
C        LA SECONDE ARETE
         DO 160 J=1,NBSUBA
            NX2 = XPX( NBSUBA, J )
            NY2 = YPX( NBSUBA, J )
            CALL XVTRAIT( NX1, NY1,  NX2, NY2 )
            NX1 = NX2
            NY1 = NY2
 160     CONTINUE
C
C        LA TROISIEME ARETE
         DO 170 I=NBSUBA-1,0,-1
            NX2 = XPX( I, NBSUBA )
            NY2 = YPX( I, NBSUBA )
            CALL XVTRAIT( NX1, NY1,  NX2, NY2 )
            NX1 = NX2
            NY1 = NY2
 170     CONTINUE
C
C        LA QUATRIEME ARETE
         DO 180 J=NBSUBA-1,0,-1
            NX2 = XPX( 0, J )
            NY2 = YPX( 0, J )
            CALL XVTRAIT( NX1, NY1,  NX2, NY2 )
            NX1 = NX2
            NY1 = NY2
 180     CONTINUE
      ENDIF
C
C     TRACE EVENTUEL DES TANGENTES AUX SOMMETS EN BLANC
      IF( IAVTGF .NE. 0 ) THEN
C
C        LA COULEUR DE TRACE DES TANGENTES
         CALL XVCOULEUR( NCOTGF )
C        TRACE FIN EN POINTILLE
C        TRACE EN  CONTINU(0) ou POINTILLE (1) ou DOUBLE POINTILLE(2)
         CALL XVTYPETRAIT( NTRTGF )
         CALL XVEPAISSEUR( 0 )
C
         IF( NBS .EQ. 3 ) THEN
C           TRIANGLE: DECALAGE DE 1 POUR ATTEINDRE LES TGS
            LEDECA = 1
         ELSE
C           QUADRANGLE: PAS DE DECALAGE POUR ATTEINDRE LES TGS
            LEDECA = 0
         ENDIF
         K = NBS
         DO 220 I=1,NBS
C
C           LE SOMMET I: LES COORDONNEES EN PIXELS SUR L'ECRAN
            NX1 = NUPXEX( X(I) )
            NY1 = NUPXEY( Y(I) )
C
            DO 210 J=1,2
C
C              LA J-EME TANGENTE AU SOMMET I
               K = K + 1
               IF( NOSOEL(K+LEDECA) .EQ. 0 ) GOTO 210
C              IL EXISTE UNE VRAIE TANGENTE
C              LES COORDONNEES EN PIXELS SUR L'ECRAN DE LA POINTE DE LA TG
               NX2 = NUPXEX( X(I) + X(K) )
               NY2 = NUPXEY( Y(I) + Y(K) )
C
C              LE TRACE DE LA TANGENTE
               CALL XVTRAIT( NX1, NY1,  NX2, NY2 )
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
