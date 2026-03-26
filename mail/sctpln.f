      SUBROUTINE SCTPLN( NOAXE,  VNORMAL,
     %                   PTPLAN, COPLAN0, COPLANMIN, COPLANMAX,
     %                   COPLAN, XYZPLAN )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    A PARTIR DU POINT PTPLAN RETROUVER LA COORDONNEE DU PLAN
C -----    SELON SA NORMALE ET LES COORDONNEES DES 4 SOMMETS DU CARRE
C          DE VISUALISATION DU PLAN
C ENTREES:
C --------
C NOAXE  : =1 AXE X
C          =2 AXE Y
C          =3 AXE Z
C          =4 AXE AVEC VNORMAL=NORMALE AU PLAN
C VNORMAL: VECTEUR NORMAL AU PLAN POUR NOAXE=4
C PTPLAN : PLAN CENTRAL DU CARRE ET PLAN DE SECTION ISSU D'UN POINTE SOURIS
C COPLAN0: COORDONNEE NORMALE AU POINT VU
C COPLANMIN:COORDONNEE NORMALE MINIMALE DE L'OBJET
C COPLANMAX:COORDONNEE NORMALE MAXIMALE DE L'OBJET
C
C SORTIES:
C --------
C COPLAN : COORDONNEE NORMALE ACTUELLE DU PLAN DE SECTION
C XYZPLAN: 3 COORDONNEES DES 4 SOMMETS DU CARRE VISUALISANT LE PLAN
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Laboratoire JL LIONS UPMC Paris Decembre 2006
C23456...............................................................012
      include"./incl/xyzext.inc"
C
      REAL  VNORMAL(3), PTPLAN(3)
      REAL  COPLAN0, COPLAN, COPLANMIN, COPLANMAX
      REAL  XYZPLAN(3,4)
C
      REAL  VNORMAL2(3), VNORMAL3(3)
C
C     ECART MAXIMAL MIN-MAX = DIAGONALE DE L'HEXAEDRE
      ECMX = SQRT( (COOEXT(1,2)-COOEXT(1,1))**2 +
     %             (COOEXT(2,2)-COOEXT(2,1))**2 +
     %             (COOEXT(3,2)-COOEXT(3,1))**2 )
C
C     DEPLACEMENT DE LA SOURIS SANS BOUTON ENFONCE
      IF ( NOAXE .LE. 3 ) THEN
C
C        COORDONNEE SELON L'AXE NOAXE
         COPLAN = PTPLAN(NOAXE)
C        PROJECTION EVENTUELLE SUR LES EXTREMES
         COPLAN = MIN( COPLAN, COPLANMAX )
         COPLAN = MAX( COPLAN, COPLANMIN )
         ECM    = ECMX * 0.1
C
C        LE RECTANGLE DANS LE PLAN DE SECTION
         IF( NOAXE .EQ. 1 ) THEN
C           RECTANGLE DANS LE PLAN YZ
            XYZPLAN(1,1) = COPLAN
            XYZPLAN(2,1) = COOEXT(2,1) - ECM
            XYZPLAN(3,1) = COOEXT(3,1) - ECM
C
            XYZPLAN(1,2) = COPLAN
            XYZPLAN(2,2) = COOEXT(2,2) + ECM
            XYZPLAN(3,2) = COOEXT(3,1) - ECM
C
            XYZPLAN(1,3) = COPLAN
            XYZPLAN(2,3) = COOEXT(2,2) + ECM
            XYZPLAN(3,3) = COOEXT(3,2) + ECM
C
            XYZPLAN(1,4) = COPLAN
            XYZPLAN(2,4) = COOEXT(2,1) - ECM
            XYZPLAN(3,4) = COOEXT(3,2) + ECM
C
         ELSE IF( NOAXE .EQ. 2 ) THEN
C           RECTANGLE DANS LE PLAN XZ
            XYZPLAN(1,1) = COOEXT(1,1) - ECM
            XYZPLAN(2,1) = COPLAN
            XYZPLAN(3,1) = COOEXT(3,1) - ECM
C
            XYZPLAN(1,2) = COOEXT(1,2) + ECM
            XYZPLAN(2,2) = COPLAN
            XYZPLAN(3,2) = COOEXT(3,1) - ECM
C
            XYZPLAN(1,3) = COOEXT(1,2) + ECM
            XYZPLAN(2,3) = COPLAN
            XYZPLAN(3,3) = COOEXT(3,2) + ECM
C
            XYZPLAN(1,4) = COOEXT(1,1) - ECM
            XYZPLAN(2,4) = COPLAN
            XYZPLAN(3,4) = COOEXT(3,2) + ECM
C
         ELSE
C           RECTANGLE DANS LE PLAN XY
            XYZPLAN(1,1) = COOEXT(1,1) - ECM
            XYZPLAN(2,1) = COOEXT(2,1) - ECM
            XYZPLAN(3,1) = COPLAN
C
            XYZPLAN(1,2) = COOEXT(1,2) + ECM
            XYZPLAN(2,2) = COOEXT(2,1) - ECM
            XYZPLAN(3,2) = COPLAN
C
            XYZPLAN(1,3) = COOEXT(1,2) + ECM
            XYZPLAN(2,3) = COOEXT(2,2) + ECM
            XYZPLAN(3,3) = COPLAN
C
            XYZPLAN(1,4) = COOEXT(1,1) - ECM
            XYZPLAN(2,4) = COOEXT(2,2) + ECM
            XYZPLAN(3,4) = COPLAN
C
         ENDIF
C
      ELSE
C
C        COORDONNEE SELON LA DIRECTION NORMALE AU PLAN
         COPLAN = PROSCR( PTPLAN, VNORMAL, 3 )
C        PROJECTION EVENTUELLE SUR LES EXTREMES
         COPLAN = MIN( COPLAN, COPLANMAX )
         COPLAN = MAX( COPLAN, COPLANMIN )
C
C        LARGEUR/2 DU CARRE A TRACER DE VISUALISATION DU PLAN
         ECM = ECMX * 0.5
C
C        DEFINITION D'UN AXE ORTHOGONAL A VNORMAL A PARTIR DE Z
C        VNORMAL2 = Z PV VNORMAL
         VNORMAL2(1) =-VNORMAL(2)
         VNORMAL2(2) = VNORMAL(1)
         VNORMAL2(3) = 0
C
         R = SQRT( VNORMAL2(1)**2 + VNORMAL2(2)**2 )
         IF( R .GE. 1E-3 ) THEN
C
C           Z N'EST PAS COLINEAIRE A VNORMAL
C           --------------------------------
            VNORMAL2(1) = VNORMAL2(1) / R
            VNORMAL2(2) = VNORMAL2(2) / R
C
C           VNORMAL3 = VNORMAL PV VNORMAL2
            VNORMAL3(1) = -VNORMAL2(2) * VNORMAL(3)
            VNORMAL3(2) =  VNORMAL2(1) * VNORMAL(3)
            VNORMAL3(3) =  VNORMAL2(2) * VNORMAL(1)
     %                      - VNORMAL2(1) * VNORMAL(2)
C
         ELSE IF( VNORMAL(3) .GE. 0 ) THEN
C
C           Z EST COLINEAIRE A VNORMAL
C           --------------------------
C           e1=z=VNORMAL   e2=x   e3=y
            VNORMAL2(1) = 1
            VNORMAL2(2) = 0
            VNORMAL2(3) = 0
            VNORMAL3(1) = 0
            VNORMAL3(2) = 1
            VNORMAL3(3) = 0
         ELSE
C
C           Z EST COLINEAIRE A -VNORMAL
C           ---------------------------
C           e1=-z   e2= x   e3=-y
            VNORMAL2(1) = 1
            VNORMAL2(2) = 0
            VNORMAL2(3) = 0
            VNORMAL3(1) = 0
            VNORMAL3(2) =-1
            VNORMAL3(3) = 0
         ENDIF
C
C        LE POINT CENTRAL DU PLAN DE SECTION SUR PTVAXO - VNORMAL
         CP = COPLAN - COPLAN0
         PTPLAN(1) = (COOEXT(1,1)+COOEXT(1,2))*0.5+CP*VNORMAL(1)
         PTPLAN(2) = (COOEXT(1,1)+COOEXT(1,2))*0.5+CP*VNORMAL(2)
         PTPLAN(3) = (COOEXT(1,1)+COOEXT(1,2))*0.5+CP*VNORMAL(3)
C
C        LES 4 COINS DU CARRE A TRACER
         DO 10 K=1,3
            XYZPLAN(K,1)=PTPLAN(K) - ECM*VNORMAL2(K) - ECM*VNORMAL3(K)
            XYZPLAN(K,2)=PTPLAN(K) + ECM*VNORMAL2(K) - ECM*VNORMAL3(K)
            XYZPLAN(K,3)=PTPLAN(K) + ECM*VNORMAL2(K) + ECM*VNORMAL3(K)
            XYZPLAN(K,4)=PTPLAN(K) - ECM*VNORMAL2(K) + ECM*VNORMAL3(K)
 10      CONTINUE
      ENDIF
C
      RETURN
      END
