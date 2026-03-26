C
C     FICHIER ~/incl/cm4vue.inc
C
C     PTAXE3 : LES 3 AXES SELON L'AXONOMETRIE
C     PV     : POINT VISE PAR L'AXONOMETRIE
C     OEIL   : POSITION DE L'OEIL
C     CADRSA : CADRE DE SAISIE MIN ET MAX DES COORDONNEES (.,1) ET (.,2)
C                              PAS DANS LA DIRECTION POUR (.,3)
      DOUBLE PRECISION  PTAXE3(3,3)
      REAL              PV(3),OEIL(3),CADRSA(3,3)
      COMMON / CM4VUE / PTAXE3,PV,OEIL,CADRSA,
     %                  YAXMIN,YAXMAX,ZAXMIN,ZAXMAX,
     %                  XZX1,XZX2,XZY1,XZY2,
     %                  YZX1,YZX2,YZY1,YZY2,
     %                  XYX1,XYX2,XYY1,XYY2,
     %                  AXX1,AXX2,AXY1,AXY2
