      SUBROUTINE TRPTDT( NP, PXYD, NOTRIA, ZEMEPS, UNPEPS,
     %                   NOTRI1, CB1, CB2, CB3, NCOTE )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    RECHERCHER UN TRIANGLE NOTRI1 CONTENANT LE POINT PXYD(NP)
C -----
C
C ENTREES:
C --------
C NP     : NUMERO DU POINT DANS PXYD
C PXYD   : TABLEAU DES COORDONNEES 2D DES POINTS
C NOTRIA : LISTE DES TRIANGLES
C                 ------- ------- ------- -------- -------- --------
C  PAR TRIANGLE : SOMMET1 SOMMET2 SOMMET3 TR_VOIS1 TR_VOIS2 TR_VOIS3
C                 ------- ------- ------- -------- -------- --------
C                 SOMMET    EST LE NUMERO DU SOMMET
C                 TR_VOIS i EST LE NUMERO DANS NOTRIA DU TRIANGLE
C                               ADJACENT PAR L'ARETE i
C ZEMEPS :  -EPSILON POUR LA COORDONNEE BARYCENTRIQUE
C UNPEPS : 1+EPSILON POUR LA COORDONNEE BARYCENTRIQUE
C
C ENTREE ET SORTIE :
C ------------------
C NOTRI1 : EN ENTREE NUMERO DANS NOTRIA DU 1-ER TRIANGLE SUSCEPTIBLE DE
C                    CONTENIR LE POINT XN YN
C          EN SORTIE NUMERO DANS NOTRIA DU TRIANGLE CONTENANT NP
C                    <=0 SI PAS DE TEL TRIANGLE DANS NOTRIA
C                    - NUMERO DU DERNIER TRIANGLE DE RECHERCHE
C                      DANS LE PARCOURS VERS LE POINT NP
C CB1,CB2,CB3 : LES 3 COORDONNEES BARYCENTIQUES DU POINT NP DANS
C               LE TRIANGLE NOTRI1
C NCOTE  : NUMERO DE 1 A 3 DE L'ARETE DE NOTRI1 EN CAS D'ARETE FRONTIERE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC    JANVIER 1991
C....................................................................012
      INTEGER           NOTRIA(6,*)
      DOUBLE PRECISION  ZEMEPS, UNPEPS
      DOUBLE PRECISION  PXYD(3,*), CB1, CB2, CB3
      DOUBLE PRECISION  XN, YN, X1, X2, X3, Y1, Y2, Y3, D, DMIN,
     %                  XMI, YMI, DMI, XNOR, YNOR, DNOR, COSA, COSMAX
C
C     ******************************************************************
C     LE TRIANGLE NOTRI1 CONTIENT IL LE POINT XN YN ?
C     ******************************************************************
      NOTRI0 = NOTRI1
      XN     = PXYD(1,NP)
      YN     = PXYD(2,NP)
C
 10   IF( NOTRI1 .LE. 0 ) GOTO 9999
      IF( NOTRIA(1,NOTRI1) .LE. 0 ) GOTO 9999
C
      X1 = PXYD( 1 , NOTRIA(1,NOTRI1) )
      X2 = PXYD( 1 , NOTRIA(2,NOTRI1) )
      X3 = PXYD( 1 , NOTRIA(3,NOTRI1) )
C
      Y1 = PXYD( 2 , NOTRIA(1,NOTRI1) )
      Y2 = PXYD( 2 , NOTRIA(2,NOTRI1) )
      Y3 = PXYD( 2 , NOTRIA(3,NOTRI1) )
C
C     LES 3 COORDONNEES BARYCENTRIQUES DU POINT XN YN DANS LE TRIANGLE NOTRI1
      D   =  ( X3 - X2 ) * Y1 + ( X1 - X3 ) * Y2 + ( X2 - X1 ) * Y3
C
      IF( D .LE. 0 ) THEN
C
C        TRIANGLE DEGENERE
C        RECHERCHE DU TRIANGLE OPPOSE DE SOMMET LE PLUS PROCHE DE NP
         DMIN = 1D28
         DO 20 I=1,3
C           LE TRIANGLE OPPOSE
            NT1 = NOTRIA(3+I,NOTRI1)
            IF( NT1 .GT. 0 ) THEN
C              LE SOMMET DE NT NON SUR L'ARETE I
               IF( NOTRIA(4,NT1) .EQ. NOTRI1 ) THEN
                  IA = 3
               ELSE IF( NOTRIA(5,NT1) .EQ. NOTRI1 ) THEN
                  IA = 1
               ELSE
                  IA = 2
               ENDIF
               NS1 = NOTRIA(IA,NT1)
               D = (PXYD(1,NS1)-PXYD(1,NP))**2
     %           + (PXYD(2,NS1)-PXYD(2,NP))**2
               IF( D .LT. DMIN ) THEN
                  NCOTE = IA
                  DMIN  = D
               ENDIF
            ENDIF
 20      CONTINUE
C
         IF( NOTRIA( 3+NCOTE, NT1 ) .NE. NOTRI1 ) THEN
C
C           SUITE DU PARCOURS A TRAVERS L'ARETE NCOTE DE NT1
            NOTRI0 = NT1
            NOTRI1 = NOTRIA( 3+NCOTE, NT1 )
            GOTO 10
         ELSE
C
C           BOUCLE INFINIE ENTRE 2 TRIANGLES
            GOTO 9999
         ENDIF
      ENDIF
C
C     LES 3 COORDONNEES BARYCENTRIQUES
      CB1 = (( X3 - X2 ) * YN +
     %       ( XN - X3 ) * Y2 +
     %       ( X2 - XN ) * Y3 ) / D
      CB2 = (( X3 - XN ) * Y1 +
     %       ( X1 - X3 ) * YN +
     %       ( XN - X1 ) * Y3 ) / D
      CB3 = 1D0 - CB1 - CB2
C
      IF( CB1 .GT. ZEMEPS .AND. CB1 .LT. UNPEPS .AND.
     %    CB2 .GT. ZEMEPS .AND. CB2 .LT. UNPEPS .AND.
     %    CB3 .GT. ZEMEPS .AND. CB3 .LT. UNPEPS ) THEN
C
C        LE TRIANGLE NOTRI1 CONTIENT LE POINT XN YN
         RETURN
C
      ELSE
C
C        RECHERCHE D'UN PROCHAIN TRIANGLE SUSCEPTIBLE DE CONTENIR LE POINT N
C        RECHERCHE DE L'ARETE DE NORMALE POINTEE SUR LE POINT N A TRAITER
         COSMAX = -3D0
         NS1    = NOTRIA(1,NOTRI1)
         DO 50 I=1,3
C
C           LES 2 SOMMETS DE L'ARETE
            NS2 = I+1
            IF( NS2 .EQ. 4 ) NS2 = 1
            NS2 = NOTRIA(NS2,NOTRI1)
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
C           LE COSINUS DE L'ANGLE : NORMALE-DIRECTION DU POINT N
            COSA = ( XMI * XNOR + YMI * YNOR ) / SQRT( DNOR * DMI )
            IF( COSA .GT. COSMAX ) THEN
C              LE NOUVEAU MAXIMUM DES COSINUS
               COSMAX = COSA
               NCOTE  = I
            ENDIF
            NS1 = NS2
 50      CONTINUE
C
C        LE NOUVEAU TRIANGLE A TESTER AVEC PROTECTION DE L'ANCIEN
         NOTRI0 = NOTRI1
         NOTRI1 = NOTRIA(NCOTE+3,NOTRI1)
         GOTO 10
      ENDIF
C
C     AUCUN TRIANGLE NE CONTIENT LE POINT XN YN
 9999 NOTRI1 = - ABS( NOTRI0 )
      END
