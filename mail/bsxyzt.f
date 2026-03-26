      SUBROUTINE BSXYZT( KDEGRE, LR, R, NBINBS, S, NBPCBL, XYZPC,
     %                   NBSOLI,
     %                   MNPAST, MNXYZS, MNXYZT )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LES 3 COORDONNEES DES NBSOLI SOMMETS DE LA LIGNE
C -----    ET SES TANGENTES D'UNE B-SPLINE POLYNOMIALE UNIFORME OUVERTE
C          CALCULER LE TABLEAU NSEF ET CELA DE FACON QUE
C          LES LONGUEURS DES ARETES VERIFIENT LA FONCTION TAILLE_IDEALE(X,Y,Z)
C
C ENTREES:
C --------
C KDEGRE : DEGRE DES POLYNOMES DE LA B-SPLINE
C LR     : NOMBRE DE PARAMETRES DIFFERENTS DE T
C NBINBS : NOMBRE D'INTERVALLES DE CALCUL DE LA B-SPLINE
C R      : LES ABSCISSES PARAMETRE AYANT POUR IMAGE LES POINTS CONTROLE
C S      : LES COEFFICIENTS DES POLYNOMES SUR CHAQUE INTERVALLE
C NBPCBL : NOMBRE DE POINTS DE CONTROLE
C XYZPC  : LES 3 COORDONNEES DES NBPCBL POINTS DE CONTROLE
C NBSOLI : NOMBRE DE SOMMETS DE LA LIGNE
C          SI KDEGRE>1 ALORS (NBSOLI-1)*2 TANGENTES SONT CALCULEES
C MNPAST : ADRESSE MCN DU DEBUT DES PARAMETRES DES SOMMETS DE LA LIGNE
C MNXYZS : ADRESSE MCN DU DEBUT DES COORDONNEES DES NBSOLI SOMMETS DE LA LIGNE
C MNXYZT : ADRESSE MCN DU DEBUT DES COORDONNEES DES TANGENTES DE LA LIGNE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
CAUTEUR : A.PERRONNET  ANALYSE NUMERIQUE UPMC PARIS            AOUT 1997
C2345X7..............................................................012
      include"./incl/pp.inc"
      COMMON            RMCN(MOTMCN)
      REAL              R(0:LR),
     %                  S(0:KDEGRE,0:NBINBS-1,1:3),
     %                  XYZPC(1:3,1:NBPCBL),
     %                  TG(3)
C
C     LES 3 COORDONNEES DU PREMIER SOMMET DE LA LIGNE OUVERTE
C     C'EST LE PREMIER POINT DE CONTROLE  DE LA LIGNE OUVERTE
      RMCN(MNXYZS  ) = XYZPC(1,1)
      RMCN(MNXYZS+1) = XYZPC(2,1)
      RMCN(MNXYZS+2) = XYZPC(3,1)
C
C     LA TANGENTE INITIALE DU PREMIER ARC
      TG(1) = S(1,0,1)
      TG(2) = S(1,0,2)
      TG(3) = S(1,0,3)
C
      R1    = R(0)
      MNS   = MNXYZS - 1
      MNT   = MNXYZT - 4
C
      DO 100 K=2,NBSOLI-1
C
C        LE PARAMETRE DU SOMMET K
         R2 = RMCN(MNPAST-1+K)
C
C        RECHERCHE DE L'INTERVALLE I CONTENANT R2
         I = 0
 70      IF( R2 .GE. R(I+1) ) THEN
C           PASSAGE A L'INTERVALLE SUIVANT DE R
            I = I + 1
            GOTO 70
         ENDIF
C        R2 EST DANS L'INTERVALLE R(I) R(I+1)
         RR = R2 - R(I)
C
C        LES 3 COORDONNEES DU SOMMET K DE LA B-SPLINE
C        --------------------------------------------
         MNS = MNS + 3
         DO 76 J=1,3
            A = S(KDEGRE,I,J)
            DO 74 M=KDEGRE-1,0,-1
               A = A * RR + S(M,I,J)
 74         CONTINUE
            RMCN( MNS + J ) = A
 76      CONTINUE
C
C        CALCUL EFFECTIF DES 2 TANGENTES DE L'ARC K-1 DE LA B-SPLINE
C        -----------------------------------------------------------
         IF( KDEGRE .GT. 1 ) THEN
C
C           LES 2 TANGENTES DE CETTE ARETE K-1
C
C           1: LA TANGENTE INITIALE DE L'ARC K-1
            PAS = R2 - R1
            MNT = MNT + 3
            RMCN(MNT+1) = TG(1) * PAS
            RMCN(MNT+2) = TG(2) * PAS
            RMCN(MNT+3) = TG(3) * PAS
C
C           2: LA TANGENTE FINALE DE L'ARC K-1
            MNT = MNT + 3
            DO 90 J=1,3
               A = KDEGRE * S(KDEGRE,I,J)
               DO 80 M=KDEGRE-1,1,-1
                  A = A * RR + M * S(M,I,J)
 80            CONTINUE
               TG(J) = A
               RMCN( MNT + J ) = -A * PAS
 90         CONTINUE
         ENDIF
C
         R1 = R2
 100  CONTINUE
C
C     LES 3 COORDONNEES DU DERNIER SOMMET DE LA LIGNE OUVERTE
C     C'EST LE DERNIER POINT DE CONTROLE  DE LA LIGNE OUVERTE
      MNS = MNS + 3
      RMCN(MNS+1) = XYZPC(1,NBPCBL)
      RMCN(MNS+2) = XYZPC(2,NBPCBL)
      RMCN(MNS+3) = XYZPC(3,NBPCBL)
C
C     CALCUL DES COMPOSANTES DE LA TANGENTE INITIALE DU DERNIER ARC
      IF( KDEGRE .GT. 1 ) THEN
         PAS   = R(LR) - R1
         MNT   = MNT + 3
         RMCN(MNT+1) = TG(1) * PAS
         RMCN(MNT+2) = TG(2) * PAS
         RMCN(MNT+3) = TG(3) * PAS
C
C        CALCUL DES COMPOSANTES DE LA TANGENTE DU DERNIER SOMMET DU DERNIER ARC
         I   = LR - 1
         RR  = R(LR) - R(I)
         MNT = MNT + 3
         DO 120 J=1,3
            A = KDEGRE * S(KDEGRE,I,J)
            DO 110 M=KDEGRE-1,1,-1
               A = A * RR + M * S(M,I,J)
 110        CONTINUE
            RMCN(MNT+J) = -A * PAS
 120     CONTINUE
      ENDIF
      RETURN
      END
