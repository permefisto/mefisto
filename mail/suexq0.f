      SUBROUTINE SUEXQ0( NBS1, NBS2, COSO,
ccc     %                   COSOBO, ALPHA,
     %                   ITRAV)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    MAILLER UN QUADRANGLE DANS UN PLAN DEFINI PAR SES BORDS
C -----    PAR LA METHODE ALGEBRIQUE CONTROLE ORTHOGONALITE
C
C ENTREES :
C ---------
C           NBS1   : NOMBRE DE POINTS SUR LIGNES "HORIZONTALES"
C           NBS2   : NOMBRE DE POINTS SUR LIGNES "VERTICALES"
C           COSO   : COORDONNEES DES SOMMETS DU BORD DU MAILLAGE
C
C SORTIES :
C ---------
C           COSO  : COORDONNEES DE TOUS LES SOMMETS DU MAILLAGE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: CHRISTOPHE DOURSAT ANALYSE NUMERIQUE UPMC JUILLET 1993
C23456---------------------------------------------------------------012
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU , IMPRIM , INTERA , NUNITE(29)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE       (MCN(1),RMCN(1))
      DIMENSION COSO(3,NBS1*NBS2)
ccc      DIMENSION COSOBO(3,2*(NBS1+NBS2)-3),ALPHA(2*(NBS1+NBS2)-3)
      DIMENSION ITRAV(NBS1*NBS2)
      INTEGER   SENSFR
C
C             REMPLISSAGE DE LA MATRICE COSOBO CONTENANT LES
C                   COORDONNEES DES NOEUDS DU BORD
C
C      DO 10 I=1,NBS1
C        COSOBO(1,I) = COSO(1,I)
C        COSOBO(2,I) = COSO(2,I)
C        COSOBO(3,I) = 0.E0
C   10 CONTINUE
C      DO 20 J=2,NBS2-1
C        NEUB = NBS1+J-1
C        NEU  = J*NBS1
C        COSOBO(1,NEUB) = COSO(1,NEU)
C        COSOBO(2,NEUB) = COSO(2,NEU)
C        COSOBO(3,NEUB) = 0.E0
C   20 CONTINUE
C      DO 30 I=NBS1,1,-1
C        NEUB = 2*NBS1+NBS2-I-1
C        NEU  = (NBS2-1)*NBS1+I
C        COSOBO(1,NEUB) = COSO(1,NEU)
C        COSOBO(2,NEUB) = COSO(2,NEU)
C        COSOBO(3,NEUB) = 0.E0
C   30 CONTINUE
C      DO 40 J=NBS2-1,1,-1
C        NEUB = 2*NBS1+2*NBS2-J-2
C        NEU  = (J-1)*NBS1+1
C        COSOBO(1,NEUB) = COSO(1,NEU)
C        COSOBO(2,NEUB) = COSO(2,NEU)
C        COSOBO(3,NEUB) = 0.E0
C   40 CONTINUE
C
C         CALCUL DES ANGLES ENTRE DEUX ARETES DU CONTOUR POUR
C
C      CALL CALANG(NEUB,COSOBO,1,SENSFR,ALPHA)
C      IF (SENSFR.EQ.0) THEN
C        NBLGRC(NRERR) = 2
C        KERR(1) = 'ERREUR SUEXQ0:'
C        KERR(2) = 'BORD DEGENERE'
C        CALL LEREUR
C        RETURN
C      ENDIF
C
C                          INTERPOLATION
C
      DO 110 J=2,NBS2-1
        DO 100 I=2,NBS1-1
          NEU  =    (J-1)*NBS1 +I
          NEUB =                I
          NEUH = (NBS2-1)*NBS1 +I
          NEUG =    (J-1)*NBS1 +1
          NEUD =    (J-1)*NBS1 +NBS1
          XI   = (I-1.)/(NBS1-1.)
          ET   = (J-1.)/(NBS2-1.)
          ALUN = 1 -XI
          ALDE =    XI
          BEUN = 1 -ET
          BEDE =    ET
C
C                     4 TERMES DE BORD ET 4 HEXSECS
C
          DO 90 K=1,3
            COSO (K,NEU) =   ALUN * COSO(K,NEUG)
     S                     + ALDE * COSO(K,NEUD)
     S                     + BEUN * COSO(K,NEUB)
     S                     + BEDE * COSO(K,NEUH)
     S                     - ALUN * BEUN * COSO(K,1)
     S                     - ALDE * BEUN * COSO(K,NBS1)
     S                     - ALUN * BEDE * COSO(K,(NBS2-1)*NBS1+1)
     S                     - ALDE * BEDE * COSO(K,NBS1*NBS2)
   90     CONTINUE
  100   CONTINUE
  110 CONTINUE
C
C     VERIFICATION DE LA BIJECTIVITE DU MAILLAGE
C
      NBELTR = -1
      CALL VERGEO(NBS1,NBS2,COSO,SENSFR,NBELTR,NBELTD,
     S            ITRAV)
      RETURN
      END
