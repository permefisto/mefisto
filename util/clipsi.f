      SUBROUTINE CLIPSI( R3 , S3 , R4 , S4 ,
     %                   R1 , R2 , TRACER  )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CLIPPER ENTRE R1 R2 LE SEGMENT (R3,S3) (R4,S4)
C -----
C
C ENTREES ET SORTIES :
C --------------------
C R3,S3 : COORDONNEES DU POINT INITIAL
C R4,S4 : COORDONNEES DU POINT FINAL
C
C ENTREES :
C ---------
C R1,R2 : MIN ET MAXIMUM DE L'INTERVALLE DE CLIPPAGE
C
C SORTIE :
C --------
C TRACER : VRAI SI LE SEGMENT EST A TRACER, FAUX SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS        MAI 1990
C23456---------------------------------------------------------------012
      LOGICAL TRACER
C
C     SI ABSCISSE DU POINT 1 > ABSCISSE DU POINT 2 ALORS PERMUTATION
      NUPERM = 0
      IF( R3 .GT. R4 ) THEN
         A  = R3
         R3 = R4
         R4 = A
         A  = S3
         S3 = S4
         S4 = A
         NUPERM = 1
      ENDIF
      IF( R4 .GE. R1 .AND. R3 .LE. R2 ) THEN
C        LE SEGMENT EST VISIBLE
         TRACER = .TRUE.
         IF( R3 .LT. R1 ) THEN
C           INTERSECTION DU SEGMENT ET R=R1
            S3 = ( (R4-R1) * S3 + (R1-R3) * S4 ) / (R4-R3)
            R3 = R1
         ENDIF
         IF( R4 .GT. R2 ) THEN
C           INTERSECTION DU SEGMENT ET R=R2
            S4 = ( (R4-R2) * S3 + (R2-R3) * S4 ) / (R4-R3)
            R3 = R2
         ENDIF
      ELSE
C        SEGMENT EN DEHORS DU SEGMENT R1-R2
         TRACER = .FALSE.
      ENDIF
C
      IF( NUPERM .NE. 0 ) THEN
C        REMISE DANS L'ORDRE INITIAL
         A  = R3
         R3 = R4
         R4 = A
         A  = S3
         S3 = S4
         S4 = A
      ENDIF
      END
