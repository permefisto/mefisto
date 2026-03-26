         SUBROUTINE  REXYT1( XP, YP, M, XYSTR,
     %                       NTY, I, J, CB )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    RETROUVER DANS LE TRIANGLE RECTANGLE UNITE STRUCTUR'E
C -----    REGULIEREMENT, LE SOUS-TRIANGLE CONTENANT LE POINT
C          PUIS, CALCULER SES 3 COORDONNEES BARYCENTRIQUES
C          RETOURNER SON IDENTIFICATION
C
C ENTREES:
C --------
C XP,YP  : LES 2 COORDONNEES DU POINT DU TRIANGLE RECTANGLE UNITE A RETROUVER
C M      : NOMBRE DE SOMMETS PAR ARETE DU TRIANGLE
C XYSTR  : LES 2 COORDONNEES DES SOMMETS DES SOUS-TRIANGLES DU TRIANGLE UNITE
C
C SORTIES:
C --------
C NTY    : =1:TRIANGLE AVEC SOMMET 1 "VERS L'ORIGINE", -1 SINON
C I,J    : INDICES PERMETTANT DE CALCULER LE NUMERO DES 3 SOMMETS
C          DU SOUS-TRIANGLE CONTENANT LE POINT (CF NS1,NS2,NS3 PLUS BAS)
C CB     : LES 3 COORDONNEES BARYCENTRIQUES DU POINT (XP,YP) DANS LE
C          SOUS TRIANGLE (I,J) DE TYPE NTY DU TRIANGLE RECTANGLE UNITE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 1997
C23456---------------------------------------------------------------012
      REAL      XYSTR(2,*)
      REAL      CB(3)
C
C     LA FONCTION FORMULE DE CALCUL DU NUMERO DU SOMMET (IS,JS) DANS LE TRIANGLE
      NUSOTR( IS, JS ) = ( IS * IS - IS ) / 2 + JS
C
C     LE PREMIER TRIANGLE SUSCEPTIBLE DE CONTENIR LE POINT (XP,YP)
C     TRIANGLE AVEC SOMMET 1 "VERS L'ORIGINE"
      NTY = 1
      J   = M / 2
      I   = J + 1
      M1  = M-1
C
C     BOUCLE SUR LES SOUS-TRIANGLES SUSCEPTIBLES DE CONTENIR LE POINT
C     NUMERO DES 3 SOMMETS DE CE SOUS-TRIANGLE POSSIBLE
 5    IF( NTY .GT. 0 ) THEN
C        TRIANGLE DE SOMMET 1 "VERS L'ORIGINE"
         NS1 = NUSOTR(I  ,J  )
         NS2 = NUSOTR(I+1,J  )
         NS3 = NUSOTR(I+1,J+1)
      ELSE
C        TRIANGLE DE SOMMET 1 "NON VERS L'ORIGINE"
         NS1 = NUSOTR(I+1,J+1)
         NS2 = NUSOTR(I  ,J+1)
         NS3 = NUSOTR(I  ,J  )
      ENDIF
C
C     LE POINT (XP,YP) EST IL INTERNE AU SOUS-TRIANGLE DE
C     SOMMETS NS1, NS2, NS3?
      X1 = XYSTR( 1 , NS1 )
      Y1 = XYSTR( 2 , NS1 )
C
      X2 = XYSTR( 1 , NS2 )
      Y2 = XYSTR( 2 , NS2 )
C
      X3 = XYSTR( 1 , NS3 )
      Y3 = XYSTR( 2 , NS3 )
C
C     2 FOIS LA SURFACE DU TRIANGLE = DETERMINANT DE LA MATRICE
C     DE CALCUL DES COORDONNEES BARYCENTRIQUES DU POINT P
      D  = ( X2 - X1 ) * ( Y3 - Y1 ) - ( X3 - X1 ) * ( Y2 - Y1 )
C
      IF( D .LE. 0 ) THEN
C
C        TRIANGLE DEGENERE POURQUOI?
C        =================
         CALL XVPAUSE
         IERR = 1
         RETURN
      ENDIF
C
C     TRIANGLE NON DEGENERE
C     =====================
C     CALCUL DES 3 COORDONNEES BARYCENTRIQUES DU
C     POINT XP YP DANS LE SOUS-TRIANGLE
      CB(1) = ( ( X2-XP ) * ( Y3-YP ) - ( X3-XP ) * ( Y2-YP ) ) / D
      IF( ABS(CB(1)) .LE. 1E-6 ) CB(1) = 0
      CB(2) = ( ( X3-XP ) * ( Y1-YP ) - ( X1-XP ) * ( Y3-YP ) ) / D
      IF( ABS(CB(2)) .LE. 1E-6 ) CB(2) = 0
      CB(3) = 1 - CB(1) - CB(2)
C
C     CALCUL DU NOMBRE DE COORDONNEES BARYCENTRIQUES NEGATIVES
      NB = 0
      DO 10 K=1,3
         IF( ABS(CB(K)) .LT. 1E-6 ) CB(K)=0
         IF( CB(K) .LT. 0 ) NB = NB + 1
 10   CONTINUE
C
C     LE POINT EST IL DERRIERE UN SOMMET? C-A-D
C     EXISTE T IL 2 COORDONNEES BARYCENTIQUES<0?
      IF( NB .EQ. 2 ) THEN
C
C        LE POINT EST DERRIERE LE SOMMET SK
C        CAR SA SEULE COORDONNEE BARYCENTRIQUE >0 EST K
         IF( CB(1) .GT. 0 ) THEN
C
C           DERRIERE LE SOMMET 1
            IF( NTY .GT. 0 ) THEN
C                              3'____2'
C              TRIANGLES         \  /
C                                 1'
C                                 1
C                              2 /__\3
C              LES NOUVEAUX INDICES ET LE TYPE DU SOUS-TRIANGLE
               IF( I .EQ. 1 ) GOTO 50
               IF( J .EQ. 1 ) GOTO 50
               I   = I - 1
               J   = J - 1
               NTY = -NTY
               GOTO 5
C
            ELSE
C                               3____2
C              TRIANGLES         \  /
C                                 1
C                                 1'
C                              2'/__\3'
C
               IF( I .EQ. M1 ) GOTO 50
               I   = I + 1
               NTY = -NTY
               GOTO 5
            ENDIF
C
         ELSE IF( CB(2) .GT. 0 ) THEN
C
C           DERRIERE LE SOMMET 2
            IF( NTY .GT. 0 ) THEN
C
C              TRIANGLES             1
C                          2'___3'=2/__\3
C                            \ /
C                             1'
C              LES NOUVEAUX INDICES ET LE TYPE DU SOUS-TRIANGLE
               IF( I .EQ. M1 ) GOTO 50
               IF( J .EQ.  1 ) GOTO 50
               I   = I + 1
               J   = J - 1
               NTY = -NTY
               GOTO 5
            ELSE
C                                   1'
C                                  /  \
C                           3___2=2'___3'
C              TRIANGLES     \ /
C                             1
C              LES NOUVEAUX INDICES ET LE TYPE DU SOUS-TRIANGLE
               IF( I .EQ. 1  ) GOTO 50
               IF( J .EQ. M1 ) GOTO 50
               I   = I - 1
               J   = J + 1
               NTY = -NTY
               GOTO 5
            ENDIF
C
         ELSE
C
C           DERRIERE LE SOMMET 3
            IF( NTY .GT. 0 ) THEN
C
C              TRIANGLES    1
C                         2/__\3=3'___2'
C                                 \  /
C                                  1'
C              LES NOUVEAUX INDICES ET LE TYPE DU SOUS-TRIANGLE
               IF( I .EQ. M1 ) GOTO 50
               IF( J .EQ. M1 ) GOTO 50
               I = I + 1
               J = J + 1
               NTY = -NTY
               GOTO 5
C
            ELSE
C                            1'
C                           / \
C              TRIANGLES  2'___3'=3___2
C                                  \ /
C                                   1
               IF( I .EQ. 1 ) GOTO 50
               IF( J .EQ. 1 ) GOTO 50
               I   = I - 1
               J   = J - 1
               NTY = -NTY
               GOTO 5
            ENDIF
         ENDIF
      ENDIF
C
      IF( NB .EQ. 1 ) THEN
C
C        1 COORDONNEE BARYCENTRIQUE<0
         IF( CB(1) .LT. 0 ) THEN
C
            IF( NTY .GT. 0 ) THEN
C
C                            1
C                           / \
C              TRIANGLES   2___3
C                          2'__3'
C                           \ /
C                           1'
               IF( I .EQ. M1 ) GOTO 50
               I   = I + 1
               NTY = -NTY
               GOTO 5
C
            ELSE
C
C                           1'
C                          / \
C             TRIANGLES   2'__3'
C                         2___3
C                          \ /
C                           1
              IF( I .EQ. 1 ) GOTO 50
              I   = I - 1
              NTY = -NTY
              GOTO 5
C
           ENDIF
C
        ELSE IF( CB(2) .LT. 0 ) THEN
C
           IF( NTY .GT. 0 ) THEN
C
C                           1 3'__2'
C                          / \ \ /
C             TRIANGLES   2___3 1'
C
              NTY = -NTY
              GOTO 5
C
           ELSE
C
C                           1' 3___2
C                          / \  \ /
C             TRIANGLES   2'__3' 1
C
              NTY = -NTY
              GOTO 5
C
           ENDIF
C
        ELSE
C
C          CB(3)<0
           IF( NTY .GT. 0 ) THEN
C
C                          3'___2' 1
C                           \  /  / \
C              TRIANGLES      1' 2___3
C
               IF( J .EQ. 1 ) GOTO 50
               J   = J - 1
               NTY = -NTY
               GOTO 5
C
            ELSE
C
C                          3____2 1'
C                           \  / / \
C              TRIANGLES      1 2'__3'
C
               IF( J .EQ. M1 ) GOTO 50
               J   = J + 1
               NTY = -NTY
               GOTO 5
C
            ENDIF
         ENDIF
      ENDIF
C
C     PAS DE COORDONNEES NEGATIVES => POINT INTERNE OU SUR LES COTES DU TRIANGLE
C
C     LE POINT EST DANS LE SOUS-TRIANGLE DE SOMMETS NS1,NS2,NS3
C     AVEC LES COORDONNEES BARYCENTRIQUES CB(1), CB(2), CB(3)
C     CALCUL DE L'INTERPOLATION P1 DU POINT DANS LE SOUS-TRIANGLE
C
 50   RETURN
      END
