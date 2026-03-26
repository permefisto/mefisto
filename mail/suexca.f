      SUBROUTINE SUEXCA(NBS1,NBS2,COSO,NCARRE,XLBAS,XLHAU,XLGAU,XLDRO,
     S                  XTRAV,NCARREMET)
C***********************************************************************
C BUT :    MAILLER UN QUADRANGLE DANS UN PLAN DEFINI PAR SES BORDS
C          PAR LA METHODE ALGEBRIQUE CONTROLE ORTHOGONALITE
C***********************************************************************
C
C ENTREE:
C           NCARRE : CHOIX DE LA METHODE 0 SEGMENTS
C                                        1 HERMITE
C                                        2 IMAGE DE SEGMENTS PAR TOPO
C           NBS1   : NOMBRE DE POINTS SUR LIGNES "HORIZONTALES"
C           NBS2   : NOMBRE DE POINTS SUR LIGNES "VERTICALES"
C           COSO   : COORDONNEES DES SOMMETS DU BORD DU MAILLAGE
C
C SORTIES:
C           COSO  : COORDONNEES DE TOUS LES SOMMETS DU MAILLAGE
C***********************************************************************
C AUTEUR: CHRISTOPHE DOURSAT ANALYSE NUMERIQUE UPMC JUILLET 1991
C****+7**************************************************************012
C***********************************************************************
C             DECLARATION DES VARIABLES ET DES TABLEAUX
C                  ET INITIALISATION DES VARIABLES
C***********************************************************************
      DIMENSION COSO(3,NBS1*NBS2), XTRAV(4,NBS1*NBS2)
C
      REAL M(2,3)
C
      TESTCV = 1.E-4
      XIMA   = 0
      DER    = 0
C***********************************************************************
C
C     CALCUL DU CARRE DE REFERENCE
C
C***********************************************************************
      DO 60 J=2,NBS2-1
        DO 50 I=2,NBS1-1
          NEU  =    (J-1)*NBS1 +I
          NEUB =                I
          NEUH = (NBS2-1)*NBS1 +I
          NEUG =    (J-1)*NBS1 +1
          NEUD =    (J-1)*NBS1 +NBS1
          XB = COSO(1,NEUB)
          XH = COSO(1,NEUH)
          YG = COSO(2,NEUG)
          YD = COSO(2,NEUD)
C-----------------------------------------------------------------------
C         CALCUL DU CARRE PAR INTERSECTION DE SEGMENTS = 0
C-----------------------------------------------------------------------
          X = XB+YG*(XH-XB)
          X = X / (1.E0-(XH-XB)*(YD-YG))
          Y = YG+XB*(YD-YG)
          Y = Y / (1.E0-(XH-XB)*(YD-YG))
C-----------------------------------------------------------------------
C         CALCUL DU CARRE PAR INTERSECTION D'HERMITES = 1
C-----------------------------------------------------------------------
          IF (NCARRE.EQ.1) THEN
            DO 10 K=1,200
              M(1,1) = 1
              M(1,2) = 6*(XH-XB)*(Y**2-Y)
              M(1,3) = -X-(XH-XB)*(2*Y**3-3*Y**2)+XB
              M(2,1) = 6*(YD-YG)*(X**2-X)
              M(2,2) = 1
              M(2,3) = -Y-(YD-YG)*(2*X**3-3*X**2)+YG
              DENOMI = M(1,1)*M(2,2)-M(2,1)*M(1,2)
              IF (DENOMI.NE.0) THEN
                DX = ( M(1,3)*M(2,2)-M(2,3)*M(1,2) ) / DENOMI
                DY = ( M(1,1)*M(2,3)-M(2,1)*M(1,3) ) / DENOMI
              ELSE
                DX = 0.E0
                DY = 0.E0
              ENDIF
              TEST = ABS(DX/X)
              TEST = MAX (TEST,ABS(DY/Y))
              X    = X+DX
              Y    = Y+DY
              IF (TEST.LT.TESTCV) THEN
                GOTO 20
              ENDIF
   10       CONTINUE
C            WRITE (*,*) " "
C            WRITE (*,*) "NON CONVERGENCE "
C            WRITE (*,*) " "
   20       CONTINUE
          ELSE IF (NCARRE.EQ.2) THEN
C-----------------------------------------------------------------------
C         CALCUL DU CARRE PAR INTERSECTION D'INVERSES D'HERMITES = 2
C-----------------------------------------------------------------------
            NMET = 1
            X0 = X
            Y0 = Y
   21       CONTINUE
            DO 30 K=1,200
             IF (NMET.EQ.1) THEN
              BET0   =  2.E0*Y**3-3.E0*Y**2+1.E0
              BET1   = -2.E0*Y**3+3.E0*Y**2
              DEN    =  BET0*XLBAS + BET1*XLHAU
              M(1,1) =  1
              M(1,2) = -6*(XH-XB)*(Y**2-Y)/DEN**2
              M(1,3) = -X+(BET0*XLBAS*XB + BET1*XLHAU*XH)/DEN
              BET0   =  2.E0*X**3-3.E0*X**2+1.E0
              BET1   = -2.E0*X**3+3.E0*X**2
              DEN    =  BET0*XLGAU + BET1*XLDRO
              M(2,1) = -6*(YD-YG)*(X**2-X)/DEN**2
              M(2,2) =  1
              M(2,3) = -Y+(BET0*XLGAU*YG + BET1*XLDRO*YD)/DEN
              DENOMI = M(1,1)*M(2,2)-M(2,1)*M(1,2)
              IF (DENOMI.NE.0) THEN
                DX = ( M(1,3)*M(2,2)-M(2,3)*M(1,2) ) / DENOMI
                DY = ( M(1,1)*M(2,3)-M(2,1)*M(1,3) ) / DENOMI
              ELSE
                DX = 0.E0
                DY = 0.E0
              ENDIF
              TEST = ABS(DX/X)
              TEST = MAX (TEST,ABS(DY/Y))
              X    = X+DX
              Y    = Y+DY
              IF ((X.LT.0.0).OR.(X.GT.1.0).OR.
     S            (Y.LT.0.0).OR.(Y.GT.1.0)) THEN
C                WRITE (*,*) " "
C                WRITE (*,*) "BASCULE NMET = 2"
C                WRITE (*,*) " "
                NMET = 2
                X = X0
                Y = Y0
                GOTO 21
              ENDIF
              IF (TEST.LT.TESTCV) THEN
                GOTO 40
              ENDIF
             ELSE
              BET0   =  2.E0*Y**3-3.E0*Y**2+1.E0
              BET1   = -2.E0*Y**3+3.E0*Y**2
              DEN    =  BET0*XLBAS + BET1*XLHAU
              X1 = ( BET0*XLBAS*XB + BET1*XLHAU*XH ) / DEN
              BET0   =  2.E0*X**3-3.E0*X**2+1.E0
              BET1   = -2.E0*X**3+3.E0*X**2
              DEN    =  BET0*XLGAU + BET1*XLDRO
              Y1 = ( BET0*XLGAU*YG + BET1*XLDRO*YD ) / DEN
              DX = ( X1 - X) / 2.0
              DY = ( Y1 - Y) / 2.0
              TEST = ABS(DX/X)
              TEST = MAX (TEST,ABS(DY/Y))
              X    = X+DX
              Y    = Y+DY
              IF (TEST.LT.TESTCV) THEN
                GOTO 40
              ENDIF
             ENDIF
   30       CONTINUE
C            WRITE (*,*) " "
C            WRITE (*,*) "NON CONVERGENCE PASSAGE : ",NMET
C            WRITE (*,*) " "
            IF (NMET.EQ.1) THEN
              NMET = 2
              X = X0
              Y = Y0
              GOTO 21
            ENDIF
   40       CONTINUE
          ELSE IF (NCARRE.EQ.3) THEN
C-----------------------------------------------------------------------
C         CALCUL DU CARRE PAR INTERSECTION D'INVERSES D'HERMITES = 2
C         AVEC PRISE EN COMPTE DE LA COURBURE
C-----------------------------------------------------------------------
            DO 130 K=1,200
              IF (NCARREMET.LT.4) THEN
              GAM0   =  XTRAV(1,I)
              GAM1   =  XTRAV(1,(NBS2-1)*NBS1+I)
              BET0   =   1.E0 - GAM0*Y
     S                 + (-3.E0+2.E0*GAM0)*Y**2 + (2.E0-GAM0)*Y**3
              BET1   =   ( 3.E0+     GAM1)*Y**2 +(-2.E0-GAM1)*Y**3
              DBET0  =        - GAM0
     S             +2.E0*(-3.E0+2.E0*GAM0)*Y+3.E0*(2.E0-GAM0)*Y**2
              DBET1  =2.E0*( 3.E0+     GAM1)*Y+3.E0*(-2.E0-GAM1)*Y**2
              ELSE
              GAM00   =  XTRAV(1,I)
              GAM10   =  XTRAV(2,I)
              GAM01   =  XTRAV(1,(NBS2-1)*NBS1+I)
              GAM11   =  XTRAV(2,(NBS2-1)*NBS1+I)
              BE0    =   1.E0 - 3.E0*Y**2 + 2.E0*Y**3
              BE1    =          3.E0*Y**2 - 2.E0*Y**3
              DBE0   =        - 6.E0*Y    + 6.E0*Y**2
              DBE1   =          6.E0*Y    - 6.E0*Y**2
              B10    =      Y - 2.E0*Y**2 +      Y**3
              B11    =        -      Y**2 +      Y**3
              DB10   =   1.E0 - 4.E0*Y    + 3.E0*Y**2
              DB11   =        - 2.E0*Y    + 3.E0*Y**2
              IF (NCARREMET.EQ.4) THEN
              CET0   = (1.E0-Y)*GAM00+Y*GAM10
              CET1   = (1.E0-Y)*GAM01+Y*GAM11
              DCET0  = GAM10-GAM00
              DCET1  = GAM11-GAM01
              BET0   = BE0 - CET0*B10
              DBET0  = DBE0 - DCET0*B10 - CET0*DB10
              BET1   = BE1 - CET1*B11
              DBET1  = DBE1 - DCET1*B11 - CET1*DB11
              ELSE IF (NCARREMET.EQ.5) THEN
              CET0   = (1.E0-Y**2)*GAM00+Y**2*GAM10
              CET1   = (1.E0-Y)**2*GAM01+(1.E0-(1.E0-Y)**2)*GAM11
              DCET0  = 2.E0*(GAM10-GAM00)
              DCET1  = 2.E0*(1.E0-Y)*(GAM11-GAM01)
              BET0   = BE0 - CET0*B10
              DBET0  = DBE0 - DCET0*B10 - CET0*DB10
              BET1   = BE1 - CET1*B11
              DBET1  = DBE1 - DCET1*B11 - CET1*DB11
              ELSE IF (NCARREMET.EQ.6) THEN
              CET0   = (1.E0-Y)**2*GAM00+(1.E0-(1.E0-Y)**2)*GAM10
              CET1   = (1.E0-Y**2)*GAM01+Y**2*GAM11
              DCET0  = 2.E0*(1.E0-Y)*(GAM10-GAM00)
              DCET1  = 2.E0*(GAM11-GAM01)
              BET0   = BE0 - CET0*B10
              DBET0  = DBE0 - DCET0*B10 - CET0*DB10
              BET1   = BE1 - CET1*B11
              DBET1  = DBE1 - DCET1*B11 - CET1*DB11
              ELSE IF (NCARREMET.EQ.7) THEN
              CET0   = 27.E0/4.E0*Y*(1.E0-Y)**2*GAM00
     S                +(1.E0-27.E0/4.E0*Y*(1.E0-Y)**2)*GAM10
              CET1   = (1.E0-27.E0/4.E0*Y**2*(1.E0-Y))*GAM01+
     S                 27.E0/4.E0*Y**2*(1.E0-Y)*GAM11
              DCET0  = -27.E0/4.E0*(1.E0-Y)*(1.E0-3.E0*Y)
     S                 *(GAM10-GAM00)
              DCET1  =  27.E0/4.E0*(2.E0*Y-3.E0*Y**2)
     S                 *(GAM11-GAM01)
              BET0   = BE0 - CET0*B10
              DBET0  = DBE0 - DCET0*B10 - CET0*DB10
              BET1   = BE1 - CET1*B11
              DBET1  = DBE1 - DCET1*B11 - CET1*DB11
              ELSE IF (NCARREMET.EQ.8) THEN
              CET0   = BE0*GAM00+BE1*GAM10
              CET1   = BE0*GAM01+BE1*GAM11
              DCET0  = DBE0*GAM00+DBE1*GAM10
              DCET1  = DBE0*GAM01+DBE1*GAM11
              BET0   = BE0 - CET0*B10
              DBET0  = DBE0 - DCET0*B10 - CET0*DB10
              BET1   = BE1 - CET1*B11
              DBET1  = DBE1 - DCET1*B11 - CET1*DB11
              ELSE IF(NCARREMET.EQ.9) THEN
                GAM0   =  XTRAV(1,I)
                GAM1   =  XTRAV(2,(NBS2-1)*NBS1+I)
                BET0   =   1.E0 - GAM0*Y
     S                   + (-3.E0+2.E0*GAM0)*Y**2 + (2.E0-GAM0)*Y**3
                BET1   =   ( 3.E0+     GAM1)*Y**2 +(-2.E0-GAM1)*Y**3
                DBET0  =        - GAM0
     S               +2.E0*(-3.E0+2.E0*GAM0)*Y+3.E0*(2.E0-GAM0)*Y**2
                DBET1  =2.E0*( 3.E0+     GAM1)*Y+3.E0*(-2.E0-GAM1)*Y**2
                DEN    =  BET0*XLBAS + BET1*XLHAU
                DER1   = -XLBAS*XLHAU*(XH-XB)
     S                    *(BET0*DBET1-BET1*DBET0)/DEN**2
                XIMA1  = (BET0*XLBAS*XB + BET1*XLHAU*XH)/DEN
                GAM0   =  XTRAV(2,I)
                GAM1   =  XTRAV(1,(NBS2-1)*NBS1+I)
                BET0   =   1.E0 - GAM0*Y
     S                   + (-3.E0+2.E0*GAM0)*Y**2 + (2.E0-GAM0)*Y**3
                BET1   =   ( 3.E0+     GAM1)*Y**2 +(-2.E0-GAM1)*Y**3
                DBET0  =        - GAM0
     S               +2.E0*(-3.E0+2.E0*GAM0)*Y+3.E0*(2.E0-GAM0)*Y**2
                DBET1  =2.E0*( 3.E0+     GAM1)*Y+3.E0*(-2.E0-GAM1)*Y**2
                DEN    =  BET0*XLBAS + BET1*XLHAU
                DER2   = -XLBAS*XLHAU*(XH-XB)
     S                    *(BET0*DBET1-BET1*DBET0)/DEN**2
                XIMA2  = (BET0*XLBAS*XB + BET1*XLHAU*XH)/DEN
                DEN2   =  BET0*XLBAS + BET1*XLHAU
                DER2   = -XLBAS*XLHAU*(XH-XB)
     S                    *(BET0*DBET1-BET1*DBET0)/DEN**2
                XIMA2  = (BET0*XLBAS*XB + BET1*XLHAU*XH)/DEN
                XIMA = (1.E0-2.E0*Y)**2*XIMA2+4.E0*Y*(1.E0-Y)*XIMA1
                DER = (1.E0-2.E0*Y)**2*DER2+4.E0*Y*(1.E0-Y)*DER1
     S               +4.E0*(1.E0-2.E0*Y)*(XIMA1-XIMA2)
              ELSE
              GAM00   =  XTRAV(3,I)
              GAM10   =  XTRAV(4,I)
              GAM01   =  XTRAV(3,(NBS2-1)*NBS1+I)
              GAM11   =  XTRAV(4,(NBS2-1)*NBS1+I)
              B00    =   1.E0 - 3.E0*Y**2 + 2.E0*Y**3
              B01    =          3.E0*Y**2 - 2.E0*Y**3
              DB00   =        - 6.E0*Y    + 6.E0*Y**2
              DB01   =          6.E0*Y    - 6.E0*Y**2
              BE0    = (1.E0-Y)*GAM00+Y*GAM10
              BE1    = (1.E0-Y)*GAM01+Y*GAM11
              DBE0   = GAM10-GAM00
              DBE1   = GAM11-GAM01
              GAM00   =  XTRAV(1,I)
              GAM10   =  XTRAV(2,I)
              GAM01   =  XTRAV(1,(NBS2-1)*NBS1+I)
              GAM11   =  XTRAV(2,(NBS2-1)*NBS1+I)
              B10    =      Y - 2.E0*Y**2 +      Y**3
              B11    =        -      Y**2 +      Y**3
              DB10   =   1.E0 - 4.E0*Y    + 3.E0*Y**2
              DB11   =        - 2.E0*Y    + 3.E0*Y**2
              CET0   = (1.E0-Y)*GAM00+Y*GAM10
              CET1   = (1.E0-Y)*GAM01+Y*GAM11
              DCET0  = GAM10-GAM00
              DCET1  = GAM11-GAM01
              BET0   = BE0*B00 - CET0*B10
              DBET0  = DBE0*B00+BE0*DB00 - DCET0*B10 - CET0*DB10
              BET1   = BE1*B01 - CET1*B11
              DBET1  = DBE1*B01+BE1*DB01 - DCET1*B11 - CET1*DB11
              ENDIF
              ENDIF
              IF (NCARREMET.NE.9) THEN
                DEN    =  BET0*XLBAS + BET1*XLHAU
                DER    = -XLBAS*XLHAU*(XH-XB)
     S                    *(BET0*DBET1-BET1*DBET0)/DEN**2
                XIMA   = (BET0*XLBAS*XB + BET1*XLHAU*XH)/DEN
              ENDIF
              M(1,1) =  1
              M(1,2) =  DER
              M(1,3) = -X+XIMA
              IF (NCARREMET.LT.4) THEN
              GAM0   =  XTRAV(1,(J-1)*NBS1+1)
              GAM1   =  XTRAV(1,J*NBS1)
              BET0   =   1.E0 - GAM0*X
     S                 + (-3.E0+2.E0*GAM0)*X**2 + (2.E0-GAM0)*X**3
              BET1   =   ( 3.E0+     GAM1)*X**2 +(-2.E0-GAM1)*X**3
              DBET0  =        - GAM0
     S             +2.E0*(-3.E0+2.E0*GAM0)*X+3.E0*(2.E0-GAM0)*X**2
              DBET1  =2.E0*( 3.E0+     GAM1)*X+3.E0*(-2.E0-GAM1)*X**2
              ELSE
              GAM00   =  XTRAV(1,(J-1)*NBS1+1)
              GAM10   =  XTRAV(2,(J-1)*NBS1+1)
              GAM01   =  XTRAV(1,J*NBS1)
              GAM11   =  XTRAV(2,J*NBS1)
              BE0    =   1.E0 - 3.E0*X**2 + 2.E0*X**3
              BE1    =          3.E0*X**2 - 2.E0*X**3
              DBE0   =        - 6.E0*X    + 6.E0*X**2
              DBE1   =          6.E0*X    - 6.E0*X**2
              B10    =      X - 2.E0*X**2 +      X**3
              B11    =        -      X**2 +      X**3
              DB10   =   1.E0 - 4.E0*X    + 3.E0*X**2
              DB11   =        - 2.E0*X    + 3.E0*X**2
              IF (NCARREMET.EQ.4) THEN
              CET0   = (1.E0-X)*GAM00+X*GAM10
              CET1   = (1.E0-X)*GAM01+X*GAM11
              DCET0  = GAM10-GAM00
              DCET1  = GAM11-GAM01
              BET0   = BE0 - CET0*B10
              DBET0  = DBE0 - DCET0*B10 - CET0*DB10
              BET1   = BE1 - CET1*B11
              DBET1  = DBE1 - DCET1*B11 - CET1*DB11
              ELSE IF (NCARREMET.EQ.5) THEN
              CET0   = (1.E0-X**2)*GAM00+X**2*GAM10
              CET1   = (1.E0-X)**2*GAM01+(1.E0-(1.E0-X)**2)*GAM11
              DCET0  = 2.E0*(GAM10-GAM00)
              DCET1  = 2.E0*(1.E0-X)*(GAM11-GAM01)
              BET0   = BE0 - CET0*B10
              DBET0  = DBE0 - DCET0*B10 - CET0*DB10
              BET1   = BE1 - CET1*B11
              DBET1  = DBE1 - DCET1*B11 - CET1*DB11
              ELSE IF (NCARREMET.EQ.6) THEN
              CET0   = (1.E0-X)**2*GAM00+(1.E0-(1.E0-X)**2)*GAM10
              CET1   = (1.E0-X**2)*GAM01+X**2*GAM11
              DCET0  = 2.E0*(1.E0-X)*(GAM10-GAM00)
              DCET1  = 2.E0*(GAM11-GAM01)
              BET0   = BE0 - CET0*B10
              DBET0  = DBE0 - DCET0*B10 - CET0*DB10
              BET1   = BE1 - CET1*B11
              DBET1  = DBE1 - DCET1*B11 - CET1*DB11
              ELSE IF (NCARREMET.EQ.7) THEN
              CET0   = 27.E0/4.E0*X*(1.E0-X)**2*GAM00
     S                +(1.E0-27.E0/4.E0*X*(1.E0-X)**2)*GAM10
              CET1   = (1.E0-27.E0/4.E0*X**2*(1.E0-X))*GAM01+
     S                 27.E0/4.E0*X**2*(1.E0-X)*GAM11
              DCET0  = -27.E0/4.E0*(1.E0-X)*(1.E0-3.E0*X)
     S                 *(GAM10-GAM00)
              DCET1  =  27.E0/4.E0*(2.E0*X-3.E0*X**2)
     S                 *(GAM11-GAM01)
              BET0   = BE0 - CET0*B10
              DBET0  = DBE0 - DCET0*B10 - CET0*DB10
              BET1   = BE1 - CET1*B11
              DBET1  = DBE1 - DCET1*B11 - CET1*DB11
              ELSE IF (NCARREMET.EQ.8) THEN
              CET0   = BE0*GAM00+BE1*GAM10
              CET1   = BE0*GAM01+BE1*GAM11
              DCET0  = DBE0*GAM00+DBE1*GAM10
              DCET1  = DBE0*GAM01+DBE1*GAM11
              BET0   = BE0 - CET0*B10
              DBET0  = DBE0 - DCET0*B10 - CET0*DB10
              BET1   = BE1 - CET1*B11
              DBET1  = DBE1 - DCET1*B11 - CET1*DB11
              ELSE IF(NCARREMET.EQ.9) THEN
                GAM0   =  XTRAV(1,(J-1)*NBS1+1)
                GAM1   =  XTRAV(2,J*NBS1)
                BET0   =   1.E0 - GAM0*X
     S                   + (-3.E0+2.E0*GAM0)*X**2 + (2.E0-GAM0)*X**3
                BET1   =   ( 3.E0+     GAM1)*X**2 +(-2.E0-GAM1)*X**3
                DBET0  =        - GAM0
     S               +2.E0*(-3.E0+2.E0*GAM0)*X+3.E0*(2.E0-GAM0)*X**2
                DBET1  =2.E0*( 3.E0+     GAM1)*X+3.E0*(-2.E0-GAM1)*X**2
                DEN    =  BET0*XLGAU + BET1*XLDRO
                DER1   = -XLGAU*XLDRO*(YD-YG)
     S                    *(BET0*DBET1-BET1*DBET0)/DEN**2
                XIMA1  = (BET0*XLGAU*YG + BET1*XLDRO*YD)/DEN
                GAM0   =  XTRAV(2,(J-1)*NBS1+1)
                GAM1   =  XTRAV(1,J*NBS1)
                BET0   =   1.E0 - GAM0*X
     S                   + (-3.E0+2.E0*GAM0)*X**2 + (2.E0-GAM0)*X**3
                BET1   =   ( 3.E0+     GAM1)*X**2 +(-2.E0-GAM1)*X**3
                DBET0  =        - GAM0
     S               +2.E0*(-3.E0+2.E0*GAM0)*X+3.E0*(2.E0-GAM0)*X**2
                DBET1  =2.E0*( 3.E0+     GAM1)*X+3.E0*(-2.E0-GAM1)*X**2
                DEN    =  BET0*XLGAU + BET1*XLDRO
                DER2   = -XLGAU*XLDRO*(YD-YG)
     S                    *(BET0*DBET1-BET1*DBET0)/DEN**2
                XIMA2  = (BET0*XLGAU*YG + BET1*XLDRO*YD)/DEN
                XIMA = (1.E0-2.E0*X)**2*XIMA2+4.E0*X*(1.E0-X)*XIMA1
                DER = (1.E0-2.E0*X)**2*DER2+4.E0*X*(1.E0-X)*DER1
     S               +4.E0*(1.E0-2.E0*X)*(XIMA1-XIMA2)
              ELSE
              GAM00   =  XTRAV(3,(J-1)*NBS1+1)
              GAM10   =  XTRAV(4,(J-1)*NBS1+1)
              GAM01   =  XTRAV(3,J*NBS1)
              GAM11   =  XTRAV(4,J*NBS1)
              B00    =   1.E0 - 3.E0*X**2 + 2.E0*X**3
              B01    =          3.E0*X**2 - 2.E0*X**3
              DB00   =        - 6.E0*X    + 6.E0*X**2
              DB01   =          6.E0*X    - 6.E0*X**2
              BE0    = (1.E0-X)*GAM00+X*GAM10
              BE1    = (1.E0-X)*GAM01+X*GAM11
              DBE0   = GAM10-GAM00
              DBE1   = GAM11-GAM01
              GAM00   =  XTRAV(1,(J-1)*NBS1+1)
              GAM10   =  XTRAV(2,(J-1)*NBS1+1)
              GAM01   =  XTRAV(1,J*NBS1)
              GAM11   =  XTRAV(2,J*NBS1)
              B10    =      X - 2.E0*X**2 +      X**3
              B11    =        -      X**2 +      X**3
              DB10   =   1.E0 - 4.E0*X    + 3.E0*X**2
              DB11   =        - 2.E0*X    + 3.E0*X**2
              CET0   = (1.E0-X)*GAM00+X*GAM10
              CET1   = (1.E0-X)*GAM01+X*GAM11
              DCET0  = GAM10-GAM00
              DCET1  = GAM11-GAM01
              BET0   = BE0*B00 - CET0*B10
              DBET0  = DBE0*B00+BE0*DB00 - DCET0*B10 - CET0*DB10
              BET1   = BE1*B01 - CET1*B11
              DBET1  = DBE1*B01+BE1*DB01 - DCET1*B11 - CET1*DB11
              ENDIF
              ENDIF
              IF (NCARREMET.NE.9) THEN
                DEN    =  BET0*XLGAU + BET1*XLDRO
                DER    = -XLGAU*XLDRO*(YD-YG)
     S                    *(BET0*DBET1-BET1*DBET0)/DEN**2
                XIMA   = (BET0*XLGAU*YG + BET1*XLDRO*YD)/DEN
              ENDIF
              M(2,1) =  DER
              M(2,2) =  1
              M(2,3) = -Y+XIMA
              DENOMI = M(1,1)*M(2,2)-M(2,1)*M(1,2)
              IF (DENOMI.NE.0) THEN
                DX = ( M(1,3)*M(2,2)-M(2,3)*M(1,2) ) / DENOMI
                DY = ( M(1,1)*M(2,3)-M(2,1)*M(1,3) ) / DENOMI
              ELSE
                DX = 0.E0
                DY = 0.E0
              ENDIF
              TEST = ABS(DX/X)
              TEST = MAX (TEST,ABS(DY/Y))
              X    = X+DX
              Y    = Y+DY
              IF (TEST.LT.TESTCV) THEN
                GOTO 140
              ENDIF
  130       CONTINUE
C            WRITE (*,*) " "
C            WRITE (*,*) "NON CONVERGENCE "
C            WRITE (*,*) " "
  140       CONTINUE
          ELSE
C-----------------------------------------------------------------------
C           CALCUL DU CARRE PAR INTERPOLATION DIRECTE HERMITE = 3
C-----------------------------------------------------------------------
            XI   = (I-1.)/(NBS1-1.)
            ET   = (J-1.)/(NBS2-1.)
            ALUN = 1 -XI
            ALDE =    XI
            BEUN = 1 -ET
            BEDE =    ET
C-----------------------------------------------------------------------
            X          =   ALUN * COSO(1,NEUG)
     S                   + ALDE * COSO(1,NEUD)
     S                   + BEUN * COSO(1,NEUB)
     S                   + BEDE * COSO(1,NEUH)
     S                   - ALUN * BEUN * COSO(1,1)
     S                   - ALDE * BEUN * COSO(1,NBS1)
     S                   - ALUN * BEDE * COSO(1,(NBS2-1)*NBS1+1)
     S                   - ALDE * BEDE * COSO(1,NBS1*NBS2)
C-----------------------------------------------------------------------
            Y          =   ALUN * COSO(2,NEUG)
     S                   + ALDE * COSO(2,NEUD)
     S                   + BEUN * COSO(2,NEUB)
     S                   + BEDE * COSO(2,NEUH)
     S                   - ALUN * BEUN * COSO(2,1)
     S                   - ALDE * BEUN * COSO(2,NBS1)
     S                   - ALUN * BEDE * COSO(2,(NBS2-1)*NBS1+1)
     S                   - ALDE * BEDE * COSO(2,NBS1*NBS2)
          ENDIF
          COSO(1,NEU) = X
          COSO(2,NEU) = Y
   50   CONTINUE
   60 CONTINUE
      RETURN
      END
