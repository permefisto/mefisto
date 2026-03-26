      SUBROUTINE DIPTPL( Pt, PT1, PT2, PT3, DISPLA )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LA DISTANCE DISPLA D'UN POINT Pt AU PLAN DEFINI
C -----    PAR 3 POINTS PT1, PT2, PT3
C
C ENTREES:
C --------
C Pt     : X Y Z DU POINT
C PT1, PT2, PT3 : LES 3 POINTS DU PLAN
C
C SORTIES:
C --------
C DISPLA : >=0 DISTANCE DE Pt AU PLAN DEFINI PAR LES 3 POINTS
C          =-1 SI POINTS ALIGNES OU IDENTIQUES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 1992
C2345X7..............................................................012
      DOUBLE PRECISION  Pt(3), PT1(3), PT2(3), PT3(3), DISPLA,
     %                  A(4), COORPO(3,3)

      DO K=1,3
         COORPO(K,1) = PT1(K)
         COORPO(K,2) = PT2(K)
         COORPO(K,3) = PT3(K)
      ENDDO

C     LES 4 COEFFICIENTS DE L'EQUATION DU PLAN
      CALL EQPLAD( COORPO, A, IERR )

      IF( IERR .GT. 0 ) THEN

C        POINTS IDENTIQUES OU ALIGNES
         DISPLA = -1
         RETURN

      ELSE

C        DISTANCE
         DISPLA = ABS(A(1) * Pt(1) + A(2) * Pt(2) + A(3) * Pt(3) + A(4))
     %          / SQRT( A(1)**2 + A(2)**2 + A(3)**2 )

      ENDIF

      RETURN
      END
