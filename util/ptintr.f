      SUBROUTINE PTINTR( XP, YP, XYSTR, N1, N2, N3, NSIGNE )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    LE POINT EST IL INTERNE AU TRIANGLE DE SOMMETS NOSOTR
C -----    (CF SP PTDATR COMME SP DE DEPART!)
C
C ENTREES:
C --------
C XP, YP  : LES 2 COORDONNEES DU POINT
C XYSTR   : LES 2 COORDONNEES DES SOMMETS DE LA TRIANGULATION
C N1,N2,N3: LE NUMERO DES 3 SOMMETS DU TRIANGLE
C
C SORTIES:
C --------
C NSIGNE : >0 SI LE POINT EST DANS LE TRIANGLE OU SUR UNE DES 3 ARETES
C          =0 SI LE TRIANGLE EST DEGENERE OU INDIRECT OU NE CONTIENT PAS LE POIN
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC   NOVEMBRE 1997
C....................................................................012
      REAL          XYSTR(2,*)
      REAL          X1,X2,X3, Y1,Y2,Y3, D,DD, CB1,CB2,CB3
C
      X1 = XYSTR( 1 , N1 )
      Y1 = XYSTR( 2 , N1 )
C
      X2 = XYSTR( 1 , N2 )
      Y2 = XYSTR( 2 , N2 )
C
      X3 = XYSTR( 1 , N3 )
      Y3 = XYSTR( 2 , N3 )
C
C     2 FOIS LA SURFACE DU TRIANGLE = DETERMINANT DE LA MATRICE
C     DE CALCUL DES COORDONNEES BARYCENTRIQUES DU POINT P
      D  = ( X2 - X1 ) * ( Y3 - Y1 ) - ( X3 - X1 ) * ( Y2 - Y1 )
C
      IF( D .GT. 0 ) THEN
C
C        TRIANGLE NON DEGENERE
C        =====================
C        CALCUL DES 3 COORDONNEES BARYCENTRIQUES DU
C        POINT XP YP DANS LE TRIANGLE
         CB1 = ( ( X2-XP ) * ( Y3-YP ) - ( X3-XP ) * ( Y2-YP ) ) / D
         CB2 = ( ( X3-XP ) * ( Y1-YP ) - ( X1-XP ) * ( Y3-YP ) ) / D
         CB3 = ( ( X1-XP ) * ( Y2-YP ) - ( X2-XP ) * ( Y1-YP ) ) / D

CCC         IF( CB1 .GE. -0.00005 .AND. CB1 .LE. 1.00005 .AND.
ccc         IF( CB1 .GE. 0 .AND. CB1 .LE. 1 .AND.
ccc     %       CB2 .GE. 0 .AND. CB2 .LE. 1 .AND.
ccc     %       CB3 .GE. 0 .AND. CB3 .LE. 1 ) THEN

         IF( ABS(CB1) + ABS(CB2) + ABS(CB3) .LT. 1.0001 ) THEN 
C           LE TRIANGLE NOSOTR CONTIENT LE POINT
            NSIGNE = 1
         ELSE
            NSIGNE = 0
         ENDIF
C
      ELSE
C
C        TRIANGLE DEGENERE
C        =================
C        LE POINT EST IL DU MEME COTE QUE LE SOMMET OPPOSE DE CHAQUE ARETE?
         NSIGNE = 0
         DO 10 I=1,3
C           LE SINUS DE L'ANGLE P1 P2-P1 POINT
            X1  = XYSTR(1,N1)
            Y1  = XYSTR(2,N1)
            D   = ( XYSTR(1,N2) - X1 ) * ( YP - Y1 )
     %          - ( XYSTR(2,N2) - Y1 ) * ( XP - X1 )
            DD  = ( XYSTR(1,N2) - X1 ) * ( XYSTR(2,N3) - Y1 )
     %          - ( XYSTR(2,N2) - Y1 ) * ( XYSTR(1,N3) - X1 )
            CB1 = ( XYSTR(1,N2) - X1 ) ** 2
     %          + ( XYSTR(2,N2) - Y1 ) ** 2
            CB2 = ( XP - X1 ) ** 2
     %          + ( YP - Y1 ) ** 2
            CB3 = ( XYSTR(1,N3) - X1 ) ** 2
     %          + ( XYSTR(2,N3) - Y1 ) ** 2
            IF( ABS( DD ) .LE. 1E-4 * SQRT( CB1 * CB3 ) ) THEN
C              LE POINT 3 EST SUR L'ARETE 1-2
C              LE POINT DOIT Y ETRE AUSSI
               IF( ABS( D ) .LE. 1E-4 * SQRT( CB1 * CB2 ) ) THEN
C                 POINT SUR L'ARETE
                  NSIGNE = NSIGNE + 1
               ENDIF
            ELSE
C              LE POINT 3 N'EST PAS SUR L'ARETE . TEST DES SIGNES
               IF( D * DD .GE. 0 ) THEN
                  NSIGNE = NSIGNE + 1
               ENDIF
            ENDIF
C           PERMUTATION CIRCULAIRE DES 3 SOMMETS ET ARETES
            N  = N1
            N1 = N2
            N2 = N3
            N3 = N
 10      CONTINUE
         IF( NSIGNE .NE. 3 ) NSIGNE = 0
      ENDIF
      END
