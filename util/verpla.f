      SUBROUTINE VERPLA( XYZ1 , D2D3 , NBS1 , NBS2 , COSO , IERR )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : VERIFICATION QUE LE TABLEAU DE COORDONNEES EST BIEN PLAN
C -----
C  (X       )    (          ) ( XX )    ( XX ) ( T         )( X      )
C  (Y  - P1 ) =  (  D2D3    ) ( YY ) ;  ( YY )=(   D2D3    )( Y - P1 )
C  (Z       )    (          ) ( ZZ )    ( ZZ ) (           )( Z      )
C
C          AVEC DANS LES APPLICATIONS ZZ = 0.
C     X  , Y  , Z  COORDONNEES DANS LE REPERE 3D
C     XX , YY , ZZ COORDONNEES DANS LE REPERE 3D DU PLAN P1 P2 P3 ET
C                                                DIRECTION ORTHOGONALE
C
C ENTREES :
C ---------
C XYZ1   : LE POINT P1 ORIGINE DU SYSTEME D'AXES DU PLAN
C D2D3   : LA MATRICE DE PASSAGE DEFINIE CI DESSUS
C           NBS1   : NOMBRE DE POINTS SUR LIGNES "HORIZONTALES"
C           NBS2   : NOMBRE DE POINTS SUR LIGNES "VERTICALES"
C           COSO   : COORDONNEES DES SOMMETS DU MAILLAGE
C
C SORTIE :
C --------
C IERR   : 0 SI LES POINTS SONT DANS UN MEME PLAN 1 SINON
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : CHRISTOPHE DOURSAT ANALYSE NUMERIQUE PARIS      MAI 1990
C........................................................................
      DOUBLE PRECISION D2D3(3,3)
      REAL             XYZ1(3),XYZ(3),XY3(3)
      REAL             COSO(3,NBS1*NBS2)
C
C     REVUE DE TOUS LES NOEUDS DE LA SURFACE STRUCTUREE
      XMIN = 1E20
      YMIN = 1E20
      XMAX =-1E20
      YMAX =-1E20
      ZMAX = 0.
      DO 20 J=1,NBS2
         DO 10 I=1,NBS1
            NEU = (J-1)*NBS1+I
            XYZ(1) = COSO(1,NEU)
            XYZ(2) = COSO(2,NEU)
            XYZ(3) = COSO(3,NEU)
            CALL CH3D3D( XYZ1 , D2D3 ,XYZ ,XY3 )
            XMIN = MIN(XMIN,XY3(1))
            XMAX = MAX(XMAX,XY3(1))
            YMIN = MIN(YMIN,XY3(2))
            YMAX = MAX(YMAX,XY3(2))
            ZMAX = MAX(ZMAX,ABS(XY3(3)))
 10      CONTINUE
 20   CONTINUE
      DREF = MAX(XMAX-XMIN,YMAX-YMIN)
      IF (ZMAX/DREF.LT.0.00001) THEN
         IERR = 0
      ELSE
         IERR = 1
      ENDIF
      END
