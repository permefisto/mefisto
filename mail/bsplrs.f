      SUBROUTINE BSPLRS( LIGFER , NBCOMP , KDEGRE , LR , NBINBS ,
     %                   NBSOLI , RAGEBL , R , S , NBPCBL , XYZPC ,
     %                   XYZSOM , XYZTGS , NBTGS )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LES 3 COORDONNEES DES NBSOLI-1 SOMMETS DE LA LIGNE
C -----    D'UNE B-SPLINE RATIONNELLE OUVERTE OU FERMEE
C
C ENTREES:
C --------
C LIGFER : 0 SI LIGNE OUVERTE , 1 SI LIGNE FERMEE
C NBCOMP : 3 SI B-SPLINE POLYNOMIALE, 4 SI RATIONNELLE
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
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS       MAI   1990
C2345X7..............................................................012
      REAL     R(0:LR),
     %         S(0:KDEGRE,0:NBINBS-1,1:NBCOMP),
     %         XYZPC (1:3,1:NBPCBL),
     %         XYZTGS(1:3,1:NBTGS),
     %         XYZSOM(1:3,1:NBSOLI)
C
      XYZSOM(3,1)      = R(0)
      XYZSOM(3,NBSOLI) = R(LR)
C
C     PROGRESSION GEOMETRIQUE OU NON ?
      NBS1 = NBSOLI - 1
      IF( RAGEBL .EQ. 1.0 ) THEN
C
C        SOMMETS EQUIDISTANTS DANS L'ESPACE DU PARAMETRE
C        ===============================================
         PAS = ( R(LR) - R(0) ) / NBS1
         DO 20 I=2,NBSOLI
            XYZSOM(3,I) = R(0) + (I-1) * PAS
 20      CONTINUE
C
      ELSE
C
C        SOMMETS EN PROGRESSION GEOMETRIQUE DANS L'ESPACE DU PARAMETRE
C        =============================================================
         PAS = (R(LR)-R(0)) * (1.0-RAGEBL) / (1.0-RAGEBL**NBS1)
         DO 40 I=2,NBSOLI
            XYZSOM(3,I) = XYZSOM(3,I-1) + PAS
            PAS         = PAS * RAGEBL
 40      CONTINUE
      ENDIF
C
C     PASSAGE AUX NBCOMP COORDONNEES DES SOMMETS DE LA LIGNE
C     ======================================================
      PAS_INI = XYZSOM(3,2) - XYZSOM(3,1)
      PAS_FIN = XYZSOM(3,NBSOLI) - XYZSOM(3,NBS1)
      PAS1    = PAS_FIN
      I       = 0
      D       = 1.0
      DO 100 N=1,NBS1
C        LA VALEUR DU PARAMETRE ET INTERVALLE DU PARAMETRE
         DD   = 0
         R0   = XYZSOM(3,N)
         PAS0 = PAS1
         PAS1 = XYZSOM(3,N+1) - XYZSOM(3,N)
C
 50      IF( R0 .GE. R(I+1) ) THEN
C
C           PASSAGE A L'INTERVALLE SUIVANT DE R
            I = I + 1
            GOTO 50
C
         ELSE
C
C           R0 EST DANS L'INTERVALLE R(I) R(I+1)
            RR = R0 - R(I)
C
            IF( NBCOMP .EQ. 4 ) THEN
C              B-SPLINE RATIONNELLE : CALCUL DU DENOMINATEUR
               D = S(KDEGRE,I,4)
               DO 60 M=KDEGRE-1,0,-1
                  D = D * RR + S(M,I,4)
 60            CONTINUE
C              LA DERIVEE DU DENOMINATEUR
               DD = KDEGRE * S(KDEGRE,I,4)
               DO 62 M=KDEGRE-1,1,-1
                  DD = DD * RR + M * S(M,I,4)
 62            CONTINUE
            ENDIF
C
C           LES 3 COMPOSANTES DU NUMERATEUR
            DO 80 J=1,3
               A = S(KDEGRE,I,J)
               DO 70 M=KDEGRE-1,0,-1
                  A = A * RR + S(M,I,J)
 70            CONTINUE
C              LA J-EME COORDONNEE DU SOMMET N
               XYZSOM(J,N) = A / D
C
C              LA COMPOSANTE J DE LA DERIVEE DU NUMERATEUR
               IF( NBTGS .GT. 0 ) THEN
                  DN = KDEGRE * S(KDEGRE,I,J)
                  DO 75 M=KDEGRE-1,1,-1
                     DN = DN * RR + M * S(M,I,J)
 75               CONTINUE
                  IF( NBCOMP .EQ. 4 ) THEN
C                    B-SPLINE RATIONNELLE
                     DN = ( DN * D - A * DD ) / ( D * D )
                  ENDIF
                  IF ( N .EQ. 1 ) THEN
                     XYZTGS(J,1) = DN * PAS_INI
                     IF( LIGFER .EQ. 1 ) THEN
                        XYZTGS(J,NBTGS) = DN * PAS_FIN
                     ENDIF
                  ELSE
                     XYZTGS(J,2*N-2) =  DN * PAS0
                     XYZTGS(J,2*N-1) = -DN * PAS1
                  ENDIF
               ENDIF
 80         CONTINUE
         ENDIF
 100  CONTINUE
C
C     AJOUT DU DERNIER SOMMET
      IF( LIGFER .EQ. 0 ) THEN
C        C'EST LE DERNIER POINT DE CONTROLE
         XYZSOM(1,NBSOLI) = XYZPC(1,NBPCBL)
         XYZSOM(2,NBSOLI) = XYZPC(2,NBPCBL)
         XYZSOM(3,NBSOLI) = XYZPC(3,NBPCBL)
      ELSE
C        C'EST LE PREMIER POINT DE LA LIGNE
         XYZSOM(1,NBSOLI) = XYZSOM(1,1)
         XYZSOM(2,NBSOLI) = XYZSOM(2,1)
         XYZSOM(3,NBSOLI) = XYZSOM(3,1)
      ENDIF
C
      IF( LIGFER .EQ. 0 ) THEN
C        CALCUL DES COORDONNEES DE LA TANGENTE RELATIVE AU DERNIER SOMMET
         IF( NBTGS .GT. 0 ) THEN
            DD = 0
            RR = R(LR) - R(LR-1)
            I  = LR - 1
            IF( NBCOMP .EQ. 4 ) THEN
C              B-SPLINE RATIONNELLE : CALCUL DU DENOMINATEUR
               D = S(KDEGRE,I,4)
               DO 160 M=KDEGRE-1,0,-1
                  D = D * RR + S(M,I,4)
 160           CONTINUE
C              LA DERIVEE DU DENOMINATEUR
               DD = KDEGRE * S(KDEGRE,I,4)
               DO 165 M=KDEGRE-1,1,-1
                  DD = DD * RR + M * S(M,I,4)
 165           CONTINUE
            ENDIF
C
            DO 180 J=1,3
C              LE NUMERATEUR
               A = S(KDEGRE,I,J)
               DO 170 M=KDEGRE-1,0,-1
                  A = A * RR + S(M,I,J)
 170           CONTINUE
C              LA COMPOSANTE J DE LA DERIVEE DU NUMERATEUR
               DN = KDEGRE * S(KDEGRE,I,J)
               DO 175 M=KDEGRE-1,1,-1
                  DN = DN * RR + M * S(M,I,J)
 175           CONTINUE
               IF( NBCOMP .EQ. 4 ) THEN
C                 B-SPLINE RATIONNELLE
                  DN = ( DN * D - A * DD ) / ( D * D )
               ENDIF
               XYZTGS(J,NBTGS) = -DN * PAS_FIN
 180        CONTINUE
         ENDIF
      ENDIF
      END
