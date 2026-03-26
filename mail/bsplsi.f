      SUBROUTINE BSPLSI( NUTYSU,
     %                   DEGREX, LUX, LTX, LRX, UX, TX, RX,
     %                   DEGREY, LUY, LTY, LRY, UY, TY, RY,
     %                   NBPOIN, POINTS, MNXYZL, AX, BX, AY, BY,
     %                   FACM  , MATRIX,
     %                   NORXTX, NORYTY, SX,
     %                   SPLINE, IERR )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    GENERER LES COEFFICIENTS DES POLYNOMES SUR CHAQUE QUADRANGLE
C -----    D'UNE SURFACE B-SPLINE D'INTERPOLATION

C ENTREES:
C --------
C NUTYSU : NUMERO DU TYPE DE LA SURFACE BSPLINE D'INTERPOLATION
C          3:DEFINIE PAR DES POINTS DE L'UTILISATEUR
C          4:DEFINIE PAR LES SOMMETS DE LIGNES (MEME NOMBRE)
C DEGREX : DEGRE  EN X DES POLYNOMES DE LA B-SPLINE
C DEGREY : DEGRE  EN Y DES POLYNOMES DE LA B-SPLINE

C LUX    : NOMBRE-1 DE NOEUDS D'INTERPOLATION EN X
C LUY    : NOMBRE-1 DE NOEUDS D'INTERPOLATION EN Y
C LTX    : NOMBRE-1 DE NOEUDS EN X POUR LES BX(J,M)
C LTY    : NOMBRE-1 DE NOEUDS EN Y POUR LES BY(J,M)
C LRX    : NOMBRE-1 DE NOEUDS TX IDENTIFIES EN X
C LRY    : NOMBRE-1 DE NOEUDS TY IDENTIFIES EN Y

C UX     : LES VALEURS DES NOEUDS D'INTERPOLATION EN X
C UY     : LES VALEURS DES NOEUDS D'INTERPOLATION EN Y
C TX     : LES VALEURS DES NOEUDS DE BX(J,M) EN X
C TY     : LES VALEURS DES NOEUDS DE BY(J,M) EN Y
C RX     : LES VALEURS DES NOEUDS IDENTIFIES DE TX
C RY     : LES VALEURS DES NOEUDS IDENTIFIES DE TY

C BX     : LES VALEURS BX(J,M)
C BY     : LES VALEURS BY(J,M)

C NBPOIN : LE NOMBRE TOTAL DE POINTS D'INTERPOLATION
C POINTS : LES 3 COORDONNEES DES POINTS D'INTERPOLATION (NUTYSU=3)
C MNXYZL : ADRESSE MCN DU TABLEAU XYZSOMMET DES LIGNES  (NUTYSU=4)
C AX,AY  : VALEURS INTERMEDIAIRES
C BX,BY  : VALEURS INTERMEDIAIRES

C FACM   : 1/M!
C MATRIX : LA MATRICE A INVERSER

C NORXTX : NUMERO DU DERNIER NOEUD TX DE CHAQUE POINT RX
C NORYTY : NUMERO DU DERNIER NOEUD TY DE CHAQUE POINT RY
C SX     : TABLEAU AUXILIAIRE

C SORTIES:
C --------
C SPLINE : LES COEFFICIENTS DES POLYNOMES SUR CHAQUE INTERVALLE
C IERR   : =0 SI PAS D'ERREUR
C          >0 SI ERREUR RENCONTREE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS       JUIN  1990
C2345X7..............................................................012
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/pp.inc"
      COMMON       MCN(MOTMCN)
      REAL         RMCN(1)
      EQUIVALENCE (MCN(1),RMCN(1))

      INTEGER      NORXTX(0:LRX), NORYTY(0:LRY), DEGREX, DEGREY,
     %             MNXYZL(LUY)

      REAL         UX(0:LUX), UY(0:LUY),
     %             RX(0:LRX), RY(0:LRY),
     %             TX(0:LTX+DEGREX), TY(0:LTY+DEGREY),
     %             BX(-DEGREX:1,0:DEGREX),
     %             AX(-DEGREX:0,0:DEGREX),
     %             BY(-DEGREY:1,0:DEGREY),
     %             AY(-DEGREY:0,0:DEGREY),
     %             FACM(0:*),
     %             POINTS(1:LUX,1:LUY,1:3),
     %             SX(0:DEGREX,1:3,-DEGREY:0),
     %             SPLINE(0:DEGREX,0:DEGREY,1:LRX,1:LRY,1:3)

      DOUBLE PRECISION SS, MATRIX(NBPOIN,NBPOIN+3)

C     LE TABLEAU DES 1 / M!
C     ---------------------
      FACM(0) = 1.0
      FACM(1) = 1.0
      SS      = 1.D0
      DO M=2,MAX(DEGREX,DEGREY)
        SS = SS * M
        FACM(M) = REAL( 1.D0 / SS )
      ENDDO

C     CALCUL DES PARAMETRES UX ET UY
C     ------------------------------
      SS = LUX
      DO IX=0,LUX
         UX(IX) = REAL( IX / SS )
      ENDDO

      SS = LUY
      DO IY=0,LUY
         UY(IY) = REAL( IY / SS )
      ENDDO

C     LE CALCUL DES NOEUDS TX RX NORXTX A PARTIR DES UX
C     -------------------------------------------------
      CALL BSPLS1( DEGREX, LUX, LTX, LRX, UX,
     %             TX, RX, NORXTX )

C     LE CALCUL DES NOEUDS TY RY NORXTY A PARTIR DES UY
C     -------------------------------------------------
      CALL BSPLS1( DEGREY, LUY, LTY, LRY, UY,
     %             TY, RY, NORYTY )

C     GENERATION DE LA MATRICE DES BX(IX,DX)(UX) * BY(IY,DY)(UY)
C     ==========================================================
C     MISE A ZERO DE LA MATRICE
      CALL AZEROD( NBPOIN * NBPOIN, MATRIX )

C     LES LIGNES DUES AUX POINTS INTERNES D'INTERPOLATION
C     ---------------------------------------------------
      DO IY=1,LUY-2

C        CALCUL DE BY(-DEGREY:0)(UY(I))
         CALL BSPLS2( UY(IY), DEGREY, LTY, LRY, TY, RY, NORYTY,
     %                BY, NRY )
C        LE NUMERO DU DERNIER NOEUD
         NORTY = NORYTY( NRY )

         DO IX=1,LUX-2

C           CALCUL DE BX(-DEGREX:0)(UX(I))
            CALL BSPLS2( UX(IX), DEGREX, LTX, LRX, TX, RX, NORXTX,
     %                   BX, NRX )
            NORTX = NORXTX( NRX )

C           LE NUMERO DE LA LIGNE DE LA MATRICE
            M = 1 + IX + LUX * IY
            DO KY=-DEGREY,0
C              LE VERITABLE NUMERO DE 0 A LUY - 1 DU NOEUD EN Y
               JY = NORTY + KY

               DO KX=-DEGREX,0
C                 LE VERITABLE NUMERO DE 0 A LUX - 1 DU NOEUD EN X
                  JX = NORTX + KX

C                 RANGEMENT DE BX(JX,DX)(UX(IX) * BY(JY,DY)(UY(IY)
C                 DANS LA COLONNE N
                  N = 1 + JX + LUX * JY
                  MATRIX(M,N) = BX(KX,0) * BY(KY,0)
               ENDDO
            ENDDO

         ENDDO
      ENDDO

C     LES 2 COTES HORIZONTAUX
C     -----------------------
      DO IX=1,LUX-2
C        CALCUL DE BX(-DEGREX:0)(UX(I))
         CALL BSPLS2( UX(IX), DEGREX, LTX, LRX, TX, RX, NORXTX,
     %                BX, NRX )
         NORTX = NORXTX( NRX )
         DO KX=-DEGREX,0
C           LE VERITABLE NUMERO DE 0 A LUX - 1 DU NOEUD EN X
            JX = NORTX + KX
C           RANGEMENT DE BX(JX,DX)(UX(IX) * 1
            M = 1 + IX
            N = 1 + JX
            MATRIX(M,N) = BX(KX,0)
            M = 1 + IX + LUX * (LUY-1)
            N = 1 + JX + LUX * NORYTY(LRY-1)
            MATRIX(M,N) = BX(KX,0)
         ENDDO
      ENDDO

C     LES 2 COTES VERTICAUX
C     ---------------------
      DO IY=1,LUY-2
C        CALCUL DE BY(-DEGREY:0)(UY(I))
         CALL BSPLS2( UY(IY), DEGREY, LTY, LRY, TY, RY, NORYTY,
     %                BY, NRY )
C        LE NUMERO DU DERNIER NOEUD
         NORTY = NORYTY( NRY )
C        LE NUMERO DE LA LIGNE DE LA MATRICE
         DO KY=-DEGREY,0
C           LE VERITABLE NUMERO DE 0 A LUY - 1 DU NOEUD EN Y
            JY = NORTY + KY
C           RANGEMENT DE BY(JY,DY)(UY(IY) * 1
            M = 1 + LUX * IY
            N = 1 + LUX * JY
            MATRIX(M,N) = BY(KY,0)
            M = LUX + LUX * IY
            N = 1 + NORXTX(LRX-1) + LUX * JY
            MATRIX(M,N) = BY(KY,0)
         ENDDO
      ENDDO

C     LES LIGNES REDUITES A L'IDENTITE CORRESPONDANT AUX 4 COINS DU CONTOUR
C     ---------------------------------------------------------------------
C     LE COIN(1,1)
      MATRIX(1,1) = 1D0
C     LE COIN(LUX,1)
      MATRIX(LUX,LUX) = 1D0
C     LE COIN(LUX,LUY)
      IY = LUX * LUY
      MATRIX(IY,IY) = 1D0
C     LE COIN(1,LUY)
      IY = IY - LUX + 1
      MATRIX(IY,IY) = 1D0

C     GENERATION DES 3 SECONDS MEMBRES DU SYSTEME =
C          3 COORDONNEES DES POINTS D'INTERPOLATION
C     =============================================
      IF( NUTYSU .EQ. 3 ) THEN

C        DEFINITION PAR LES POINTS D'INTERPOLATION DE L'UTILISATEUR
         DO IY = 1, LUY
            DO IX = 1, LUX
               DO NBC = 1, 3
                  MATRIX(IX+LUX*(IY-1),NBPOIN+NBC) = POINTS(IX,IY,NBC)
               ENDDO
            ENDDO
         ENDDO
C
      ELSE IF( NUTYSU .EQ. 4 ) THEN

C        DEFINITION DES POINTS D'INTERPOLATION PAR LES SOMMETS DES LIGNES EN X
         DO IY=1,LUY
C           LES POINTS DE LA LIGNE IY SONT RECHERCHES
C           L'ADRESSE MCN DU TABLEAU XYZ DES SOMMETS DE LA LIGNE
            MN  = MNXYZL( IY ) + WYZSOM -1
            DO IX = 1, LUX
               DO NBC = 1, 3
                  MN = MN + 1
                  MATRIX(IX+LUX*(IY-1),NBPOIN+NBC) = RMCN(MN)
               ENDDO
            ENDDO
         ENDDO

      ENDIF

C     INVERSION DU SYSTEME LINEAIRE
C     =============================
      CALL GAUSPT( NBPOIN, 3, MATRIX, IERR )
      IF( IERR .NE. 0 ) THEN
         NBLGRC(NRERR) = 3
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'SURFACE B-SPLINE: MAUVAIS CHOIX DES PARAMETRES'
            KERR(2) = 'ASSOCIES AUX POINTS INTERPOLATION'
            KERR(3) = 'SYSTEME NON INVERSIBLE'
         ELSE
            KERR(1) = 'B-SPLINE SURFACE: BAD CHOICE OF PARAMETERS'
            KERR(2) = 'ASSOCIATED AT INTERPOLATION POINTS'
            KERR(3) = 'The SYSTEM IS NOT INVERSIBLE'
         ENDIF
         CALL LEREUR
         RETURN
      ENDIF

C     GENERATION DES POLYNOMES DE LA B-SPLINE
C     =======================================
      DO JY=0,LRY-1

C        L'INTERVALLE JY+1 EN Y :  LA VALEUR DU PARAMETRE Y EN JY
C        ------------------------
         RYJY   = RY(JY)
C        LE NUMERO DU DERNIER NOEUD EN RY(JY)
         NORTY  = NORYTY( JY )

C        CALCUL DES BY IY,M (RY(JY))
C        ---------------------------
         CALL AZEROR( (DEGREY+2)*(DEGREY+1), BY )
         BY(0,0) = 1.0
         DO M=1,DEGREY
            DO J=-M,0
C              LE VRAI INDICE DE B
               JB = NORTY + J
C              EVALUATION DES FRACTIONS DE LA FORMULE
               IF( TY(JB+M) .NE. TY(JB) ) THEN
                  U1 = ( RYJY - TY(JB) ) / ( TY(JB+M) - TY(JB) )
               ELSE
                  U1 = 0
               ENDIF
               JBM1 = JB + M + 1
               IF( TY(JBM1) .NE. TY(JB+1) ) THEN
                  U2 = ( TY(JBM1) - RYJY     )
     %               / ( TY(JBM1) - TY(JB+1) )
               ELSE
                  U2 = 0
               ENDIF
               BY(J,M) = U1 * BY(J,M-1) + U2 * BY(J+1,M-1)
            ENDDO
         ENDDO

         DO JX=0,LRX-1

C           GENERATION DU POLYNOME EN X SUR L'INTERVALLE RX(JX)
C           DE TOUTES LES LIGNES NORTY-DEGREY A NORTY
C           ---------------------------------------------------
            CALL BSPLS4( JX, DEGREX, LRX, RX , LTX, TX,
     %                   NORTY , DEGREY, LUX, LUY,
     %                   MATRIX(1,NBPOIN+1),
     %                   BX, AX, FACM, NORXTX, SX )

            DO KX=0,DEGREX

C              LA KX DERIVATION EN X
C              ---------------------
               DO NBC = 1, 3

C                 CALCUL DES AY(JX,KX,IY,0,NBC)
                  DO J=-DEGREY,0
C                    LA MULTIPLICATION PAR KX!
                     AY(J,0) = SX( KX, NBC, J ) / FACM(KX)
                  ENDDO

C                 LE CALCUL DES AY IY,M
C                 ---------------------
                  DO M=1,DEGREY
                     KM = DEGREY - M + 1
                     DO J=-DEGREY+M,0
C                       LE VRAI INDICE DE A
                        JB = NORTY + J
                        IF( TY(JB+KM) .NE. TY(JB) ) THEN
                           AY(J,M) = KM * (AY(J,M-1)-AY(J-1,M-1))
     %                                  / (TY(JB+KM)-TY(JB))
                        ELSE
                           AY(J,M) = 0
                        ENDIF
                     ENDDO
                  ENDDO

C                 LE CALCUL PROPREMENT DIT DE DDS/DKX DKY S(RX(JX),RY(JY))
C                 --------------------------------------------------------
                  DO M=0,DEGREY
                     SS = 0
                     DO J=-DEGREY+M,0
                        SS = SS + AY(J,M) * BY(J,DEGREY-M)
                     ENDDO
                     SPLINE(KX,M,JX+1,JY+1,NBC) =
     %                                REAL( FACM(M) * FACM(KX) * SS )
                  ENDDO

               ENDDO

            ENDDO

         ENDDO

      ENDDO

      RETURN
      END
