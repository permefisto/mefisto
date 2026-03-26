      SUBROUTINE PACTQT( NBCT, MNXY, NUCOTE, NBSOCT, RLONGC, MNCOCU )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CONSTRUIRE LE TABLEAU VALEURS DU PARAMETRE
C -----   (COORDONNEE CURVILIGNE HOMOGENE) DES SOMMETS DES NBCT COTES
C          D'UN QUADRANGLE (CF SP LISTQU)
C          D'UN TRIANGLE   (CF SP LISTTR)
C
C ENTREES:
C --------
C MNXY   : ADRESSE MCN DE LA PREMIERE COORDONNEE DU PREMIER SOMMET
C          DE LA LIGNE FORMANT LE CONTOUR DU QUADRANGLE OU TRIANGLE PLAN DONT LA
C          LONGUEUR DES COTES EST CELLE DES COTES DU QUADRANGLE OU TRIANGLE COUR
C          LES SOMMETS SONT NUMEROTES DANS LE SENS DIRECT A PARTIR DE S1
C
C NUCOTE : NUCOTE(I) NUMERO D'ENTREE DE 1 A NBCT DU COTE I DANS LE SENS DIRECT
C          SI NUCOTE(3)=-2 LE COTE 3 EST EN FAIT LE COTE 2 DONN'E A PARCOURIR
C                          EN SENS INVERSE DE SON RANGEMENT...
C
C QUADRANGLE: POUR LE SENS C1:S1S2 C2:S2S3 C3:S3S4 C4:S4S1
C
C                S4----<-------S3
C                |             |
C                |             |
C               \/             /\
C                |             |
C                |             |
C                S1---->-------S2
C
C TRIANGLE:   POUR LE SENS C1:S1S2 C2:S2S3 C3:S3S1
C
C                S3
C                |  \
C                |    \
C          C3   \/     /\ C2
C                |         \
C                |            \
C                S1---->-------S2
C                      C1
C
C
C NBSOCT : NOMBRE DE SOMMETS DES NBCT COTES DU QUADRANGLE
C RLONGC : LONGUEUR DES ARETES DROITES DE CHACUN DES NBCT COTES
C
C SORTIE :
C --------
C MNCOCU : ADRESSE MCN DU TABLEAU DES VALEURS DU PARAMETRE SUR
C          LES NBCT COTES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS         FEVRIER 1997
C234567..............................................................012
      INTEGER          NBSOCT(NBCT), NUCOTE(NBCT)
      REAL             RLONGC(NBCT)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL             RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
C
      NBST = NBSOCT(1)
      DO 10 J=2, NBCT
        NBST = NBST + NBSOCT(J)
 10   CONTINUE
C
C     RESERVATION DU TABLEAU SI CE N'EST DEJA FAIT
      CALL TNMCDC( 'REEL', NBST, MNCOCU )
      MNC = MNCOCU
C
C     LES NBCT-1 PREMIERS COTES DU POLYGONE
      MN  = MNXY
      DO 20 J=1,NBCT-1
         XL1 = 0
         RMCN(MNC) = 0
         DO 15 I=2,NBSOCT(ABS(NUCOTE(J)))
            XL1 = XL1 + DIST2P( RMCN(MN), RMCN(MN+3) ) / RLONGC(J)
            MN  = MN  + 3
            MNC = MNC + 1
            RMCN(MNC) = XL1
 15      CONTINUE
         RMCN(MNC) = 1.0
         MNC = MNC + 1
 20   CONTINUE
C
C     LE DERNIER COTE
      XL1 = 0
      RMCN(MNC) = 0
      DO 30 I=2,NBSOCT(ABS(NUCOTE(NBCT)))-1
C        BOUCLE SPECIALE CAR LE DERNIER POINT ALIAS LE PREMIER
C        POINT N'EST PAS STOCKE 2 FOIS DANS LA LIGNE FERMEE
         XL1 = XL1 + DIST2P( RMCN(MN), RMCN(MN+3) ) / RLONGC(NBCT)
         MN  = MN  + 3
         MNC = MNC + 1
         RMCN(MNC) = XL1
 30   CONTINUE
      MNC = MNC + 1
      RMCN(MNC) = 1.0
      END
