      SUBROUTINE LETRPT( XY, PXYD, MOSOAR, NOSOAR, MOARTR, NOARTR,
     %                   NOTRI0, NOTRI1, COBARY )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    RECHERCHER UN TRIANGLE NOTRI1 CONTENANT LE POINT XY
C -----    EN CHERCHANT LA MEILLEURE NORMALE DES COTES POUR SE RAPPROCHER
C
C ENTREES:
C --------
C XY     : LES 2 COORDONNEES DU POINT
C PXYD   : TABLEAU DES COORDONNEES 2D DES POINTS
C MOSOAR : NOMBRE MAXIMAL D'ENTIERS PAR ARETE DU TABLEAU NOSOAR
C NOSOAR : NUMERO DES 2 SOMMETS , NO LIGNE, 2 TRIANGLES DE L'ARETE,
C          CHAINAGE DES ARETES FRONTALIERES, CHAINAGE DU HACHAGE DES ARETES
C          HACHAGE DES ARETES = H(NS1,NS2)=MIN(NS1,NS2)
C MOARTR : NOMBRE MAXIMAL D'ENTIERS PAR ARETE DU TABLEAU NOARTR
C MXARTR : NOMBRE MAXIMAL DE TRIANGLES STOCKABLES DANS LE TABLEAU NOARTR
C NOARTR : LES 3 ARETES DES TRIANGLES +-ARETE1, +-ARETE2, +-ARETE3
C          ARETE1 = 0 SI TRIANGLE VIDE => ARETE2 = TRIANGLE VIDE SUIVANT
C NOTRI0 : NUMERO DANS NOTRIA DU 1-ER TRIANGLE SUSCEPTIBLE DE
C          CONTENIR LE POINT XN YN
C
C ENTREE ET SORTIE :
C ------------------
C NOTRI1 : > 0 NUMERO DANS NOARTR DU TRIANGLE CONTENANT XY
C          <=0 SI PAS DE TEL TRIANGLE DANS NOTRIA
C COBARY : 3 COORDONNEES BARYCENTRIQUES DU POINT XY DANS LE TRIANGLE NOTRI1
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & ST PIERRE DU PERRAY  FEVRIER 2008
C2345X7..............................................................012
      include"./incl/trvari.inc"
      INTEGER           NOSOAR(MOSOAR,*),
     %                  NOARTR(MOARTR,*)
      DOUBLE PRECISION  XY(2), PXYD(3,*)
      DOUBLE PRECISION  COBARY(3)
      DOUBLE PRECISION  XN, YN, X1, X2, X3, Y1, Y2, Y3, D, DMIN,
     %                  XMI, YMI, DMI, XNOR, YNOR, DNOR, COSA, COSMAX
      INTEGER           NOSOTR(3), NOSOT1(3), NOTROP(3)
C
C     ******************************************************************
C     LE TRIANGLE NOTRI1 CONTIENT IL LE POINT XN YN ?
C     ******************************************************************
      NOTRI1 = NOTRI0
      XN     = XY(1)
      YN     = XY(2)
C
 10   IF( NOTRI1 .LE. 0 ) GOTO 9999
      IF( NOARTR(1,NOTRI1) .EQ. 0 ) GOTO 9999
C
cccC     TRACE DU TRIANGLE NOTRI1
ccc      CALL MTTRTR( PXYD, NOTRI1, MOARTR, NOARTR, MOSOAR, NOSOAR,
ccc     %             NCTURQ, NCBEIG )
C
C     LE NO DES 3 SOMMETS DU TRIANGLE NOTRI1
      CALL NUSOTR( NOTRI1, MOSOAR, NOSOAR, MOARTR, NOARTR, NOSOTR )
      X1 = PXYD( 1, NOSOTR(1) )
      Y1 = PXYD( 2, NOSOTR(1) )
      X2 = PXYD( 1, NOSOTR(2) )
      Y2 = PXYD( 2, NOSOTR(2) )
      X3 = PXYD( 1, NOSOTR(3) )
      Y3 = PXYD( 2, NOSOTR(3) )
C
C     2 FOIS LA SURFACE DU TRIANGLE = DETERMINANT DE LA MATRICE
C     DE CALCUL DES COORDONNEES BARYCENTRIQUES DU POINT P
      D  = ( X2 - X1 ) * ( Y3 - Y1 ) - ( X3 - X1 ) * ( Y2 - Y1 )
C
      NCOTE = 0
      IF( D .LE. 0 ) THEN
C
C        TRIANGLE DEGENERE
C        -----------------
C        RECHERCHE DU TRIANGLE OPPOSE DE SOMMET LE PLUS PROCHE DE NP
         DMIN = 1D123
         DO 20 I=1,3
C
C           LE TRIANGLE NT1 OPPOSE A NOTRI1 PAR L'ARETE I
            NA = ABS( NOARTR(I,NOTRI1) )
            IF( NOSOAR(4,NA) .EQ. NOTRI1 ) THEN
               NT1 = NOSOAR(5,NOTRI1)
            ELSE
               NT1 = NOSOAR(4,NOTRI1)
            ENDIF
            NOTROP(I) = NT1
C
            IF( NT1 .GT. 0 ) THEN
C              LE SOMMET OPPOSE A L'ARETE I DANS LE TRIANGLE NT1
C              LE NO DES 3 SOMMETS DU TRIANGLE NT1
               CALL NUSOTR(NT1, MOSOAR, NOSOAR, MOARTR, NOARTR, NOSOT1)
               DO 11 K=1,3
                   NS1 = NOSOT1(K)
                   IF( NOSOAR(1,NA) .NE. NS1 .AND.
     %                 NOSOAR(2,NA) .NE. NS1 ) GOTO 12
 11            CONTINUE
 12            D = (PXYD(1,NS1)-XN)**2 + (PXYD(2,NS1)-YN)**2
               IF( D .LT. DMIN ) THEN
                  NCOTE = I
                  DMIN  = D
               ENDIF
            ENDIF
 20      CONTINUE
C
         IF( NCOTE .EQ. 0 ) GOTO 9999
         IF( NOTROP(NCOTE) .NE. NOTRI1 ) THEN
C           SUITE DU PARCOURS A TRAVERS L'ARETE NCOTE DE NOTRI1
            NOTRI1 = NOTROP(NCOTE)
            GOTO 10
         ELSE
C           BOUCLE INFINIE ENTRE 2 TRIANGLES => PAS DE TRIANGLE
            GOTO 9999
         ENDIF
      ENDIF
C
C     TRIANGLE NON DEGENERE:
C     ----------------------
C     LES 3 COORDONNEES BARYCENTRIQUES
      COBARY(1) = ( ( X2-XN ) * ( Y3-YN ) - ( X3-XN ) * ( Y2-YN ) ) / D
      COBARY(2) = ( ( X3-XN ) * ( Y1-YN ) - ( X1-XN ) * ( Y3-YN ) ) / D
      COBARY(3) = ( ( X1-XN ) * ( Y2-YN ) - ( X2-XN ) * ( Y1-YN ) ) / D
C
      IF( COBARY(1) .GE. 0D0 .AND. COBARY(1) .LE. 1D0 .AND.
     %    COBARY(2) .GE. 0D0 .AND. COBARY(2) .LE. 1D0 .AND.
     %    COBARY(3) .GE. 0D0 .AND. COBARY(3) .LE. 1D0 ) THEN
c
C        LE TRIANGLE NOTRI1 CONTIENT LE POINT XN YN
cccc
cccc        copie de mempx dans fenetre
ccc         call mempxfenetre
cccc        pour vider le buffer de x11
ccc         call xvvoir
cccc        saisie d'un point par clic de la souris
cccc        ou entree d'un caractere pour voir le trace
cccc         call saiptc( notyev, nx, ny, nochar )
         RETURN
C
      ELSE
C
C        RECHERCHE D'UN PROCHAIN TRIANGLE SUSCEPTIBLE DE CONTENIR LE POINT XY
C        RECHERCHE DE L'ARETE DE NORMALE POINTEE SUR LE POINT N A TRAITER
         COSMAX = -3D0
         DO 50 I=1,3
C
C           LES 2 SOMMETS DE L'ARETE NA DU TRIANGLE NOTRI1
            NA = NOARTR(I,NOTRI1)
            IF( NA .GT. 0 ) THEN
               NS1 = NOSOAR(1,NA)
               NS2 = NOSOAR(2,NA)
            ELSE
               NA = -NA
               NS2 = NOSOAR(1,NA)
               NS1 = NOSOAR(2,NA)
            ENDIF
C
C           LE MILIEU DE L'ARETE
            XMI = XN - ( PXYD(1,NS1) + PXYD(1,NS2) ) * 0.5D0
            YMI = YN - ( PXYD(2,NS1) + PXYD(2,NS2) ) * 0.5D0
            DMI = XMI * XMI + YMI * YMI
            IF( DMI .LE. 0.0 ) THEN
C              XN YN EST LE MILIEU DE L'ARETE NS1-NS2
               RETURN
            ENDIF
C
C           LA DIRECTION NORMALE A L'ARETE
            XNOR = PXYD(2,NS2) - PXYD(2,NS1)
            YNOR = PXYD(1,NS1) - PXYD(1,NS2)
            DNOR = XNOR * XNOR + YNOR * YNOR
C
C           LE COSINUS DE L'ANGLE : NORMALE-DIRECTION DU POINT XY
            COSA = ( XMI * XNOR + YMI * YNOR ) / SQRT( DNOR * DMI )
            IF( COSA .GT. COSMAX ) THEN
C              LE NOUVEAU MAXIMUM DES COSINUS
               COSMAX = COSA
               NCOTE  = NA
            ENDIF
 50      CONTINUE
C
C        LE NOUVEAU TRIANGLE A TESTER
         IF( NOSOAR(4,NCOTE) .EQ. NOTRI1 ) THEN
            NOTRI1 = NOSOAR(5,NCOTE)
         ELSE
            NOTRI1 = NOSOAR(4,NCOTE)
         ENDIF
         GOTO 10
      ENDIF
C
C     AUCUN TRIANGLE NE CONTIENT LE POINT XN YN
 9999 NOTRI1 = 0
C
      RETURN
      END
