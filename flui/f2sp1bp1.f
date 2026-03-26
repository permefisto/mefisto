      SUBROUTINE F2SP1BP1( X,
     %                     NOOBLA, NUMILI, NUMALI, LTDELI,
     %                     NOOBSF, NUMISU, NUMASU, LTDESU,
     %                     BE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DU SECOND MEMBRE DU TRIANGLE BREZZI FORTIN
C -----    P1+BULLE CONTINU POUR LES 2 COMPOSANTES DE LA VITESSE
C          P1       CONTINU POUR LA PRESSION
C
C ENTREES:
C --------
C X      : 2 COORDONNEES DES 3 SOMMETS DU TRIANGLE
C
C NOOBLA : NUMERO DES OBJETS LIGNES DES ARETES DE L'ELEMENT FINI
C NUMILI : NUMERO MINIMAL DES OBJETS LIGNES
C NUMALI : NUMERO MAXIMAL DES OBJETS LIGNES
C LTDELI : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES DU FLUIDE
C          DES OBJETS LIGNES
C
C NOOBSF : NUMERO DE SURFACE DE CET ELEMENT FINI
C NUMISU : NUMERO MINIMAL DES OBJETS SURFACES
C NUMASU : NUMERO MAXIMAL DES OBJETS SURFACES
C LTDESU : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES DU FLUIDE
C          DES OBJETS SURFACES
C
C SORTIES:
C --------
C BE     : BE(11) LE SECOND MEMBRE ELEMENTAIRE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET Laboratoire J-L.LIONS UPMC Paris Novembre 2008
C23456---------------------------------------------------------------012
      include"./incl/donflu.inc"
      COMMON /UNITES/ LECTEU,IMPRIM,NUNITE(30)
      INTRINSIC         SQRT
C
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE      (MCN(1), RMCN(1), DMCN(1))
C
      REAL              X(3,2)
      INTEGER           NOOBLA(3)
      INTEGER           LTDELI(1:MXDOFL,NUMILI:NUMALI)
      INTEGER           LTDESU(1:MXDOFL,NUMISU:NUMASU)
      DOUBLE PRECISION  BE(11)
C
      DOUBLE PRECISION  FGAMMA(2), FSURF(2,3)
      DOUBLE PRECISION  DELTA, D, XD, YD, VN(2)
      DOUBLE PRECISION  X21, Y21, X31, Y31
C
      DOUBLE PRECISION  TPLa(4,3)
C     TPLa(i,j) = Integrale sur e chapeau Pi Lambdaj dx dy
      DATA  TPLa/
     %  0.58333333333333333D-01,    0.16666666666666667D-01,
     %  0.16666666666666667D-01,    0.75D-01,
     %  0.16666666666666667D-01,    0.58333333333333333D-01,
     %  0.16666666666666667D-01,    0.75D-01,
     %  0.16666666666666667D-01,    0.16666666666666667D-01,
     %  0.58333333333333333D-01,    0.75D-01 /
C
C     PHASE PRELIMINAIRE
C     -------------------
C     INITIALISATION DES TABLEAUX UTILISES
      DO J=1,11
        BE(J) = 0.D0
      ENDDO
C
C     1°) CONTRIBUTION DES EFFORTS SURFACIQUES INTERPOLES P1
C     ------------------------------------------------------
      IF( LTDESU(LPFORC,NOOBSF) .GT. 0 ) THEN
C
        X21 = X(2,1) - X(1,1)
        X31 = X(3,1) - X(1,1)
C
        Y21 = X(2,2) - X(1,2)
        Y31 = X(3,2) - X(1,2)
C
C       LE DETERMINANT DE LA JACOBIENNE
        DELTA = ABS( X21*Y31 - X31*Y21 )
C
C       VALEUR DES EFFORTS SURFACIQUES AUX 3 SOMMETS DU TRIANGLE
        DO J=1,3
           XD = X(J,1)
           YD = X(J,2)
           CALL REFORC(3,NOOBSF,2, XD,YD,0D0, 0D0,0D0,0D0,
     %                 LTDESU(LPFORC,NOOBSF), FSURF(1,J) )
           FSURF(1,J) = FSURF(1,J) * DELTA
           FSURF(2,J) = FSURF(2,J) * DELTA
        ENDDO
C
C       COEFFICIENTS DU SECOND MEMBRE LIES AUX EFFORTS SURFACIQUES
        DO I=1,4
           DO J=1,3
              BE(I)   = BE(I)   + TPLa(I,J) * FSURF(1,J)
              BE(I+4) = BE(I+4) + TPLa(I,J) * FSURF(2,J)
           ENDDO
        ENDDO
      ENDIF
C
C     2°) CONTRIBUTION DES EFFORTS LINEIQUES
C     --------------------------------------
      DO 50 K=1,3
C
C        NO DE LA LIGNE DE L'ARETE K
         NL = NOOBLA(K)
         IF( NL .GT. 0 ) THEN
            IF( LTDELI(LPFORC,NL) .GT. 0 ) THEN
C
C              IL EXISTE UNE FORCE DONNEE SUR CETTE LIGNE
C              LE NUMERO DES SOMMETS DU COTE K DU TRIANGLE
               IF( K .NE. 3 ) THEN
                  KK = K+1
               ELSE
                  KK = 1
               ENDIF
C
C              LE VECTEUR ORTHOGONAL AU VECTEUR TANGENT
               VN(1) = X(KK,2) - X(K ,2)
               VN(2) = X(K ,1) - X(KK,1)
C
C              LE JACOBIEN
               D = SQRT( VN(1)**2 + VN(2)**2 )
C
C              LE VECTEUR NORMAL UNITAIRE
               VN(1) = VN(1) / D
               VN(2) = VN(2) / D
C
C              LE PRODUIT DU JACOBIEN PAR LE POIDS
               D = D * 0.5D0
C
C              CALCUL DES EFFORTS LIES AU SOMMET K
               XD = X(K,1)
               YD = X(K,2)
               CALL REFORC( 2, NL, 2,
     %                      XD,YD,0D0, VN(1),VN(2),0D0,
     %                      LTDELI(LPFORC,NL), FGAMMA)
               BE(K)   = BE(K)   + D * FGAMMA(1)
               BE(K+4) = BE(K+4) + D * FGAMMA(2)
C
C              CALCUL DES EFFORTS LIES AU SOMMET K+1
               XD = X(KK,1)
               YD = X(KK,2)
               CALL REFORC( 2, NL, 2,
     %                      XD,YD,0D0, VN(1),VN(2),0D0,
     %                      LTDELI(LPFORC,NL), FGAMMA)
               BE(KK)   = BE(KK)   + D * FGAMMA(1)
               BE(KK+4) = BE(KK+4) + D * FGAMMA(2)
            ENDIF
         ENDIF
C
 50   CONTINUE
C
      RETURN
      END
