      SUBROUTINE NORARTRR( P1, P2, P3, U, IERR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CONSTRUIRE LE VECTOR NORMAL UNITAIRE A L'ARETE P1 P2
C -----    DIRIGE VERS P3

C ENTREES:
C --------
C P1     : XYZ DU SOMMET 1 DU TRIANGLE
C P2     : XYZ DU SOMMET 2 DU TRIANGLE
C P3     : XYZ DU SOMMET 3 DU TRIANGLE

C SORTIES:
C --------
C U      : VECTOR NORMAL UNITAIRE A L'ARETE P1 P2 DIRIGE VERS P3
C IERR   : =0 CALCUL CORRECT
C          =1 VECTEUR NORMAL AU TRIANGLE INCALCULABLE (TRIANGLE INCORRECT)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE DU PERRAY  Octobre 2017
C23456---------------------------------------------------------------012
      REAL              P1(3), P2(3), P3(3)
      DOUBLE PRECISION  U(3)
      DOUBLE PRECISION  P21(3), P31(3), VN(3), D, SQRT, DBLE, A, B, C, E

C     LE VECTEUR NORMAL AU TRIANGLE P1 P2 P3 : ( P2 - P1 ) x ( P3 - P1 )
      CALL VECNOR3( P1, P2, P3, VN )

C     CARRE DE LA LONGUEUR DU VECTEUR NORMAL
      D = VN(1)**2 + VN(2)**2 + VN(3)**2
      IF( D .LE. 0D0 ) THEN
C        LE TRIANGLE EST INCORRECT
         PRINT*,'norartrr: TRIANGLE DE VECTEUR NORMAL NUL'
         PRINT*,'norartrr: P1=',P1
         PRINT*,'norartrr: P2=',P2
         PRINT*,'norartrr: P3=',P3
         IERR = 1
         GOTO 9999
      ENDIF
      IERR = 0

C     LE VECTEUR NORMAL UNITAIRE (NORME 1) AU TRIANGLE P1 P2 P3
      D = SQRT( D )
      DO I=1,3
         VN(I) = VN(I) / D
      ENDDO

C     XYZ des VECTEURS P1P2 et P1P3
      DO I=1,3
         P21(I) = DBLE( P2(I) - P1(I) )
         P31(I) = DBLE( P3(I) - P1(I) )
      ENDDO

C     1-ERE EQUATION DE DEFINITION DE U : U EST NORMAL A VN
C     VN(1) * U(1) + VN(2) * U(2) + VN(3) * U(3) = 0
C     RECHERCHE DE LA PLUS GRANDE COMPOSANTE DE VN
      I1 = 0
      D  = 0D0
      DO I=1,3
         A = ABS( VN(I) )
         IF( A .GT. D ) THEN
            I1 = I
            D  = A
         ENDIF
      ENDDO

C     LES INDICES SUIVANTS I1 MODULO 3
      IF( I1 .EQ. 1 ) THEN
         I2 = 2
         I3 = 3
      ELSE IF( I1 .EQ. 2 ) THEN
         I2 = 3
         I3 = 1
      ELSE
         I2 = 1
         I3 = 2
      ENDIF

C     PAR CONSTRUCTION VN(I1) EST NON NUL =>
C     U(I1) = ( -VN(I2) * U(I2) -VN(I3) * U(I3) ) / VN(I1)

C     2-EME EQUATION DE DEFINITION DE U : U EST NORMAL A P1P2
C     P21(I1) * U(I1) + P21(I2) * U(I2) + P21(I3) * U(I3) = 0

      IF( P21(I1) .NE. 0D0 ) THEN

C        U(I1) = ( - P21(I2) * U(I2) - P21(I3) * U(I3) ) / P21(I1)
C        U(I1) = ( - VN(I2)  * U(I2) - VN(I3)  * U(I3) ) / VN(I1)
C        =>
C          ( - P21(I2)/P21(I1) + VN(I2)/VN(I1) ) * U(I2)
C        = (   P21(I3)/P21(I1) - VN(I3)/VN(I1) ) * U(I3)
C        i.e.
C             A * U(I2) = B * U(I3)
         A = - P21(I2)/P21(I1) + VN(I2)/VN(I1)
         B =   P21(I3)/P21(I1) - VN(I3)/VN(I1)

C        3-EME EQUATION DE DEFINITION DE U : U EST DE NORME 1
C        U(I1)**2 + U(I2)**2 + U(I3)**2 = 1
C        =>
C        [( -VN(I2) * B/A * U(I3) -VN(I3) * U(I3) ) / VN(I1)]**2
C      + [ B/A * U(I3)]**2 + [U(I3)]**2 = 1
C        =>
C        [( -VN(I2)*B/A -VN(I3) )/VN(I1)]**2 + [B/A]**2 + 1 ]  [U(I3)]**2 = 1
C        =>
C        [U(I3)]**2 = 1 / ([( -VN(I2)*B/A -VN(I3) )/VN(I1)]**2 + [B/A]**2 + 1 ])
         E = ( -VN(I2) * B/A - VN(I3) ) / VN(I1)
         C = E**2 + (B/A)**2 +  1

C        U(I3)**2 = 1 / C
         U(I3) = SQRT( 1D0 / C )

C        4-EME CONDITION DE DEFINITION DE U : SENS U = SENS DE P1P3
C        P31(I1) * U(I1) + P31(I2) * U(I2) + P31(I3) * U(I3) > 0
C        =>
C        P31(I1) * ( - VN(I2) * B/A * U(I3) - VN(I3)  * U(I3) ) / VN(I1) )
C      + P31(I2) * B/A * U(I3) + P31(I3) * U(I3) > 0
C        =>
C      ( P31(I1)*(-VN(I2)*B/A-VN(I3))/VN(I1) ) +P31(I2)*B/A +P31(I3) )*U(I3)> 0
         D = P31(I1)*E + P31(I2)*B/A + P31(I3)

         IF( D .LT. 0D0 ) THEN
C           LE SIGNE - DE SQRT EST LE BON SIGNE
            U(I3) = - U(I3)
         ENDIF

         U(I2) = B/A * U(I3)

         U(I1) = ( - VN(I2) * U(I2) - VN(I3) * U(I3) ) / VN(I1)
         GOTO 9999

      ELSE

C        VN(I1) =/ 0  et  P21(I1)=0D0
C        1-ERE EQUATION DE DEFINITION DE U : U EST NORMAL A VN
C        U(I1) = ( -VN(I2) * U(I2) -VN(I3) * U(I3) ) / VN(I1)

C        2-EME EQUATION DE DEFINITION DE U : U EST NORMAL A P1P2
C        P21(I2) * U(I2) + P21(I3) * U(I3) = 0

 20      IF( P21(I2) .NE. 0D0 ) THEN

C           VN(I1) =/ 0  et  P21(I1)=0D0 et P21(I2)=/0D0
C           U(I2) = -P21(I3)/P21(I2) * U(I3)
C           U(I1) = ((-VN(I2)*(-P21(I3))/P21(I2)-VN(I3))/VN(I1)) * U(I3)= 
            E =      (-VN(I2)*(-P21(I3))/P21(I2)-VN(I3))/VN(I1)
C           U(I1) = E * U(I3)

C           3-EME EQUATION DE DEFINITION DE U : U EST DE NORME 1
C           U(I1)**2 + U(I2)**2 + U(I3)**2 = 1

C        => [ E * U(I3)]**2 + [-P21(I3)/P21(I2) * U(I3)]**2 + [U(I3)]**2 = 1

C        => [E]**2 + [-P21(I3)/P21(I2)]**2 + 1 ) * [U(I3)]**2 = 1

C        => [U(I3)]**2 = 1 / ( [E]**2 + [-P21(I3)/P21(I2)]**2 + 1 )

            C = E**2 + (P21(I3)/P21(I2))**2 + 1D0

C           U(I3)**2 = 1 / C
            U(I3) = SQRT( 1D0 / C )

C           4-EME CONDITION DE DEFINITION DE U : SENS U = SENS DE P1P3
C           P31(I1) * U(I1) + P31(I2) * U(I2) + P31(I3) * U(I3) > 0
C           =>
C           P31(I1) * ( ((-VN(I2)*(-P21(I3))/P21(I2)-VN(I3))/VN(I1)) * U(I3) ) +
C           P31(I2) * ( -P21(I3)/P21(I2) * U(I3) ) + P31(I3) * U(I3) > 0
C           =>
C           ( P31(I1) * ( ((-VN(I2)*(-P21(I3))/P21(I2)-VN(I3))/VN(I1)) ) +
C           P31(I2) * ( -P21(I3)/P21(I2) ) + P31(I3) ) * U(I3) > 0

            D =  P31(I1)*( ((-VN(I2)*(-P21(I3))/P21(I2)-VN(I3))/VN(I1)))
     %         + P31(I2)*( -P21(I3)/P21(I2) ) + P31(I3)

            IF( D .LT. 0D0 ) THEN
C              LE SIGNE - DE SQRT EST LE BON SIGNE
               U(I3) = - U(I3)
            ENDIF

            U(I2) = -P21(I3)/P21(I2) * U(I3)

            U(I1) = ( - VN(I2) * U(I2) - VN(I3) * U(I3) ) / VN(I1)
            GOTO 9999

         ELSE IF( P21(I3) .NE. 0D0 ) THEN

C           VN(I1)=/0  et  P21(I1)=0D0 et P21(I2)=0D0 et P21(I3)=/0D0
C           PERMUTATION DE I2 et I3
            I  = I2
            I2 = I3
            I3 = I
            GOTO 20

         ELSE

C           VN(I1) =/ 0  et  P21(I1)=0D0 et P21(I2)=0D0 et P21(I3)=0D0
C           ARETE P1P2 REDUITE A UN POINT
            IERR = 1
            GOTO 9999

         ENDIF

      ENDIF

 9999 RETURN
      END
