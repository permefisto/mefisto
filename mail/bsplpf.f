      SUBROUTINE BSPLPF( KDEGRE , LR , NBINBS ,
     %                   NBSOLI , RAGEBL , R , S ,
     %                   XYZSOM , U , XYZTGS , NBTGS )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LES 3 COORDONNEES DES NBSOLI-1 SOMMETS DE LA LIGNE
C -----    D'UNE B-SPLINE POLYNOMIALE UNIFORME FERMEE
C          TELS QUE LES LONGUEURS DES ARETES VERIFIENT LA RAISON DE LA
C          PROGRESSION GEOMETRIQUE
C
C ENTREES:
C --------
C KDEGRE : DEGRE DES POLYNOMES DE LA B-SPLINE
C LR     : NOMBRE DE PARAMETRES DIFFERENTS DE T
C NBINBS : NOMBRE D'INTERVALLES DE CALCUL DE LA B-SPLINE
C NBSOLI : NOMBRE DE SOMMETS DE LA LIGNE
C RAGEBL : RAISON DE LA PROGRESSION GEOMETRIQUE DES SOMMETS
C R      : LES ABSCISSES PARAMETRE AYANT POUR IMAGE LES POINTS CONTROLE
C S      : LES COEFFICIENTS DES POLYNOMES SUR CHAQUE INTERVALLE
C NBPCBL : NOMBRE DE POINTS DE CONTROLE
C XYZPC  : LES 3 COORDONNEES DES NBPCBL POINTS DE CONTROLE
C NBTGS  : NOMBRE DE TANGENTES DE LA LIGNE
C
C SORTIES:
C --------
C XYZSOM : 3 COORDONNEES DES POINTS DE CONTROLE
C XYZTGS : 3 COORDONNEES DES TANGENTES
C
C TRAVAIL:
C --------
C U      : ABSCISSES CURVILIGNES SEGMENT DE REFERENCE ET SUR B-SPLINE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
CAUTEURS:A.PERRONNET C.DOURSAT ANALYSE NUMERIQUE UPMC PARIS JUILLET 1996
C2345X7..............................................................012
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      REAL     R(0:LR),
     %         S(0:KDEGRE,0:NBINBS-1,1:3),
     %         XYZSOM(1:3,1:NBSOLI),
     %         XYZTGS(1:3,1:NBTGS) ,
     %         U(1:4,1:NBSOLI+1),
     %         XYZ(3)
C
      U(1,1) = R(0)
      U(2,1) = R(0)
      U(3,1) = R(0)
      U(4,1) = R(0)
C
      U(1,NBSOLI+1) = R(LR)
      U(2,NBSOLI+1) = R(LR)
      U(3,NBSOLI+1) = R(LR)
      U(4,NBSOLI+1) = R(LR)
C
C     PROGRESSION GEOMETRIQUE OU NON
      IF( RAGEBL .EQ. 1.0 ) THEN
C
C        SOMMETS EQUIDISTANTS DANS L'ESPACE DU PARAMETRE
C        ===============================================
         PAS = ( R(LR) - R(0) ) / NBSOLI
         DO 20 I=2,NBSOLI+1
            U(1,I) = R(0) + (I-1) * PAS
            U(4,I) = U(1,I)
 20      CONTINUE
C
      ELSE
C
C        SOMMETS EN PROGRESSION GEOMETRIQUE DANS L'ESPACE DU PARAMETRE
C        =============================================================
         PAS = (R(LR)-R(0)) * (1.0-RAGEBL) / (1.0-RAGEBL**NBSOLI)
         DO 30 I=2,NBSOLI+1
            U(1,I) = U(1,I-1) + PAS
            U(4,I) = U(1,I)
            PAS    = PAS * RAGEBL
 30      CONTINUE
      ENDIF
C
C     LES 3 COORDONNEES DU PREMIER POINT DE LA LIGNE FERMEE
      DO 36 J=1,3
         XYZSOM(J,1) = S(0,0,J)
 36   CONTINUE
C
      DO 150 ITER=1,256
        I      = 0
        IDENTT = 1
        IDENTQ = 0
        DO 80 N=2,NBSOLI
C          LA VALEUR DU PARAMETRE
           R0 = U(1,N)
 40        IF( R0 .GE. R(I+1) ) THEN
C
C             PASSAGE A L'INTERVALLE SUIVANT DE R
              I = I + 1
              GOTO 40
           ELSE
C
C             R0 EST DANS L'INTERVALLE R(I) R(I+1)
              RR = R0 - R(I)
C             LES 3 COMPOSANTES DU NUMERATEUR
              DO 70 J=1,3
                 A = S(KDEGRE,I,J)
                 DO 60 M=KDEGRE-1,0,-1
                    A = A * RR + S(M,I,J)
 60              CONTINUE
                 XYZ(J) = A
70           CONTINUE
              IF (ITER.NE.1) CALL XYZIDE(XYZ,XYZSOM(1,N),IDENTQ)
              IDENTT = MIN(IDENTT,IDENTQ)
              XYZSOM(1,N) = XYZ(1)
              XYZSOM(2,N) = XYZ(2)
              XYZSOM(3,N) = XYZ(3)
           ENDIF
  80    CONTINUE
C
C       CONVERGENCE DES ITERATIONS ?
        IF( IDENTT .EQ. 1 ) GOTO 200
C
C       CALCUL DE L ABSCISSE CURVILIGNE SUR B-SPLINE
        DO 100 N=2,NBSOLI+1
           U(2,N) = 0.0
           IF( N .LE. NBSOLI ) THEN
              NN = N
           ELSE
C             LE SOMMET NBSOLI+1 EST EN FAIT LE SOMMET 1
              NN = 1
           ENDIF
           DO 90 J=1,3
              U(2,N) = U(2,N) + ( XYZSOM(J,NN) - XYZSOM(J,N-1) )**2
  90       CONTINUE
           U(2,N) = U(2,N-1) + SQRT(U(2,N))
 100    CONTINUE
        RAPP = (U(1,NBSOLI+1)-U(1,1)) / (U(2,NBSOLI+1)-U(2,1))
        DO 110 N=2,NBSOLI+1
           U(2,N) = RAPP * (U(2,N)-U(2,1))
 110    CONTINUE
C
C       INVERSION
        I = 0
        DO 130 N=2,NBSOLI
C          LA VALEUR DU PARAMETRE
           U0 = U(4,N)
 120       IF( U0 .GE. U(2,I+1) ) THEN
C             PASSAGE A L'INTERVALLE SUIVANT DE U
              I = I + 1
              GOTO 120
           ELSE
C             U0 EST DANS L'INTERVALLE U(2,I) U(2,I+1)
              UU = (U0 - U(2,I)) / (U(2,I+1) - U(2,I))
              U(3,N) = (1.0E0-UU) * U(1,I) + UU * U(1,I+1)
           ENDIF
 130    CONTINUE
C       PREPARATION DE L'ITERATION SUIVANTE
        UU = U(1,NBSOLI+1)
        DO 140 N=2,NBSOLI+1
           U(1,N) = U(3,N)
 140    CONTINUE
        U(1,NBSOLI+1) = UU
 150  CONTINUE
C
      NBLGRC(NRERR) = 1
      KERR(1) = 'BSPLPF: NON CONVERGENCE DE L''ABSCISSE CURVILIGNE'
      CALL LEREUR
C
C     CALCUL DES TANGENTES AUX SOMMETS DE LA BSPLINE
C     ==============================================
 200  IF( NBTGS .LE. 0 ) RETURN
C     LES 2 TANGENTES DU PREMIER SOMMET
      PAS1 = U(1,2) - U(1,1)
      PAS0 = U(1,NBSOLI+1) - U(1,NBSOLI)
      DO 205 J=1,3
         XYZTGS(J,    1) =  S(1,0,J) * PAS1
         XYZTGS(J,NBTGS) = -S(1,0,J) * PAS0
 205  CONTINUE
C
      I = 0
      DO 290 N=2,NBSOLI
         PAS = U(1,N)
 210     IF ( PAS .GE. R(I+1) ) THEN
C           PASSAGE A L'INTERVALLE SUIVANT DE R
            I = I + 1
            GOTO 210
         ELSE
C           PAS EST DANS L'INTERVALLE R(I) R(I+1)
            RR   = PAS - R(I)
            PAS0 = PAS1
            PAS1 = U(1,N+1) - U(1,N)
C           CALCUL DES COORDONNEES DES TANGENTES SI LE DEGRE EST PLUS GRAND OU E
            DO 280 J=1,3
               A = KDEGRE * S(KDEGRE,I,J)
               DO 260 M=KDEGRE-1,1,-1
                  A = A * RR + M * S(M,I,J)
 260           CONTINUE
               XYZTGS(J,2*N-2) = -A * PAS0
               XYZTGS(J,2*N-1) =  A * PAS1
 280        CONTINUE
         ENDIF
 290  CONTINUE
      END
