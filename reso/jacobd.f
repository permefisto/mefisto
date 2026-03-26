      SUBROUTINE JACOBD( NC,   RTOL, AR, BR, D,
     %                   EIGV, VEC, IERR )
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    EXTRACTION DES VALEURS ET VECTEURS PROPRES DU PROBLEME
C -----   ( AR - EIGV * BR ) * VEC = 0
C          AVEC AR,BR MATRICES CARREES PLEINES PAR LA METHODE DE JACOBI
C
C          CF LE RAPPORT SUR LES VALEURS PROPRES
C
C ENTREES:
C --------
C  NC    : ORDRE DE AR,BR,VEC,D
C  RTOL  : SEUIL RELATIF DE CONVERGENCE DE LA METHODE DE JACOBI
C  AR    : MATRICE CARREE D ORDRE NC
C  BR    : MATRICE CARREE D ORDRE NC
C  D     : TABLEAU AUXILIAIRE DE NC VARIABLES
C
C SORTIES:
C --------
C EIGV   : TABLEAU         DES NC VALEURS  PROPRES
C VEC    : TABLEAU (NC,NC) DES NC VECTEURS PROPRES
C IERR   : >0 SI NOMBRE D'ITERATIONS MAX ATTEINT SANS CONVERGENCE
C          =0 SI PAS D'ERREUR RENCONTREE
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEURS     : BENAZETH-GOURDIN SERVENIERE A.PERRONNET  JANVIER   1977
C MODIFICATION: A. PERRONNET ANALYSE NUMERIQUE ET INRIA  SEPTEMBRE 1979
C MODIFICATION: A. PERRONNET ANALYSE NUMERIQUE ET INRIA  NOVEMBRE  1989
C .....................................................................
      PARAMETER ( MXITER = 32 )
C     MXITER NOMBRE MAXIMAL D'ITERATIONS
      include"./incl/langue.inc"

      COMMON /UNITES/  LECTEU,IMPRIM,NUNITE(30)
      DOUBLE PRECISION AR(NC,NC),BR(NC,NC),VEC(NC,NC),EIGV(NC),D(NC)
      DOUBLE PRECISION XX,AIJ,BIJ,CX,
     +                 AJM,BIM,VECI,AK,AIM,BJJ,EPS,
     +                 ALFA,VECJ,AII,AJJ,BII,BJM,BETA

   45 FORMAT('NOMBRE D''ITERATIONS DE JACOBI = ',I6)
20045 FORMAT('NUMBER of JACOBI ITERATIONS = ',I6)
 1100 FORMAT('JACOBD: VALEURS PROPRES CALCULEES'/
     +      5(1X,I4,' : ',G15.7))
21100 FORMAT('JACOBD: COMPUTED EIGENVALUES'/
     +      5(1X,I4,' : ',G15.7))
 1110 FORMAT('PB JACOBD: NOMBRE MAXIMAL D ITERATIONS ATTEINT: ',I6)
21110 FORMAT('PB JACOBD: MAXIMUM of ITERATIONS REACHED: ',I6)

C     PREMIERES ESTIMATIONS DES VALEURS PROPRES ET DES VECTEURS PROPRES
      IERR = 0
      DO I = 1,NC
         D(I) = AR(I,I) / BR(I,I)
      ENDDO
      IF( NC .EQ. 1 ) GOTO 1000

      DO I = 1,NC
         DO J = 1,NC
            VEC(I,J) = 0
         ENDDO
         VEC(I,I) = 1
      ENDDO

C     K EST LE NOMBRE D'ITERATIONS
      K   = 1
      EPS = 1.D-3

C     LA K-EME ITERATION
 2500 IF( K .LE. 8 ) EPS = EPS * 0.01D0

      DO I = 1,NC-1
         DO 3 J = I+1,NC
C           CALCUL DES COUPLES (I,J). LES COEFFICIENTS SONT ILS PETITS?
            IF(    BR(I,J)**2 .LE. EPS*ABS(BR(I,I)*BR(J,J)) ) THEN
               IF( AR(I,J)**2 .LE. EPS*ABS(AR(I,I)*AR(J,J)) ) GOTO 3
            ENDIF
C
C           TRANSFORMATION A EFFECTUER: AR(I,J) OU BR(I,J) EST GRAND
            AII = AR(I,I) * BR(I,J) - BR(I,I) * AR(I,J)
            AJJ = AR(J,J) * BR(I,J) - BR(J,J) * AR(I,J)
            IF( ABS(AR(I,J)) .GE. ABS(BR(I,J)) ) THEN
C              TRAITEMENT AVEC AR(I,J)
               AK = ( AII * AR(J,J) - AJJ * AR(I,I) ) / AR(I,J)
            ELSE
C              TRAITEMENT AVEC BR(I,J)
               AK = ( AII * BR(J,J) - AJJ * BR(I,I) ) / BR(I,J)
            ENDIF
C
C           CALCUL DE XX
            CX = AK * AK+ 4.D0 * AII * AJJ
            IF( CX .LT. 0. ) THEN
CCCC               WRITE (IMPRIM,2001) K,I,J,CX,AR(I,I),AR(J,J),AR(I,J),
CCCC     +         BR(I,I),BR(J,J),BR(I,J),AII,AJJ,AK
CCCC               CX = 0
            GOTO 3
            ENDIF
CCCC 2001 FORMAT(' PB JACOBD: ITERATION ',I5,' INDICES I,J:',I6,',',I6,
CCCC     + ' CX=',G15.7,' AU LIEU DE >=0'/
CCCC     +       ' A(II), JJ, IJ : ',3G15.7/
CCCC     +       ' B(II), JJ, IJ : ',3G15.7/
CCCC     +       ' AII  ,AJJ, AK : ',3G15.7,' OUBLI DU COEFFICIENT I,J'/
C
            CX = SQRT(CX)
C
            IF( AK .GE. 0. ) THEN
               XX = ( AK + CX ) * 0.5D0
            ELSE
               XX = ( AK - CX ) * 0.5D0
            ENDIF
            IF(XX .NE. 0) THEN
C              XX = 0  = > AR ET BR SONT MULTIPLES SCALAIRES
C              CALCUL DE ALFA ET BETA
               ALFA =  AJJ / XX
               BETA = -AII / XX
            ELSE
               ALFA = 0
               BETA = -AR(I,J) / AR(J,J)
            ENDIF
C
C           TRANSFORMATION DE AR ET BR
            AIJ = AR(I,J)
            BIJ = BR(I,J)
            AJJ = AR(J,J)
            BJJ = BR(J,J)
            AII = AR(I,I)
            BII = BR(I,I)
            AR(I,J) = 0
            BR(I,J) = 0
            AR(J,J) = ALFA*ALFA*AII+2*ALFA*AIJ+AJJ
            BR(J,J) = ALFA*ALFA*BII+2*ALFA*BIJ+BJJ
            AR(I,I) = BETA*BETA*AJJ+2*BETA*AIJ+AII
            BR(I,I) = BETA*BETA*BJJ+2*BETA*BIJ+BII
            J1 = J+1
            II = I-1
            JJ = J-1
            DO M = 1,II
               AIM = AR(M,I)
               BIM = BR(M,I)
               AJM = AR(M,J)
               BJM = BR(M,J)
               AR(M,I) = AIM+BETA*AJM
               BR(M,I) = BIM+BETA*BJM
               AR(M,J) = ALFA*AIM+AJM
               BR(M,J) = ALFA*BIM+BJM
            ENDDO
            DO M = I+1,JJ
               AIM = AR(I,M)
               BIM = BR(I,M)
               AJM = AR(M,J)
               BJM = BR(M,J)
               AR(I,M) = AIM+BETA*AJM
               BR(I,M) = BIM+BETA*BJM
               AR(M,J) = ALFA*AIM+AJM
               BR(M,J) = ALFA*BIM+BJM
            ENDDO
            DO M = J1,NC
               AIM = AR(I,M)
               BIM = BR(I,M)
               AJM = AR(J,M)
               BJM = BR(J,M)
               AR(I,M) = AIM+BETA*AJM
               BR(I,M) = BIM+BETA*BJM
               AR(J,M) = ALFA*AIM+AJM
               BR(J,M) = ALFA*BIM+BJM
            ENDDO
            DO L = 1,NC
               VECI = VEC(L,I)
               VECJ = VEC(L,J)
               VEC(L,I) = VECI+BETA*VECJ
               VEC(L,J) = VECJ+ALFA*VECI
            ENDDO
 3       ENDDO
      ENDDO
C
C     CALCUL DES NOUVELLES ESTIMATIONS DES VALEURS PROPRES
C     COMPARAISON AVEC LES ESTIMATIONS PRECEDENTES
      J = 0
      DO I = 1,NC
         EIGV(I) = AR(I,I) / BR(I,I)
         IF( ABS(EIGV(I)-D(I)) .GT. ABS(RTOL*D(I)) ) J = 1
      ENDDO
      IF( J .GT. 0 ) GOTO 3000
C
C     VERIFICATION QUE LES COUPLES SONT INFERIEURS A LA TOLERANCE
      DO I = 1,NC-1
         DO J = I+1,NC
            IF( AR(I,J)**2 .GT. RTOL*ABS(AR(I,I)*AR(J,J)) ) GOTO 3000
            IF( BR(I,J)**2 .GT. RTOL*ABS(BR(I,I)*BR(J,J)) ) GOTO 3000
         ENDDO
      ENDDO

C     QUAND ON SORT DE CETTE BOUCLE LA CONVERGENCE EST ACCEPTEE !
C     ===========================================================
C     REMPLISSAGE DE L'AUTRE MOITIE DE AR ET DE BR
  117 DO I = 1,NC
         DO J = I+1,NC
            BR(J,I) = BR(I,J)
            AR(J,I) = AR(I,J)
         ENDDO
      ENDDO
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,45) K
      ELSE
         WRITE(IMPRIM,20045) K
      ENDIF
      RETURN
C
C     PREPARATION DE L'ITERATION SUIVANTE
C     -----------------------------------
 3000 IF( K .LE. MXITER ) THEN
         DO I = 1,NC
            D(I) = EIGV(I)
         ENDDO
         K = K + 1
         GOTO 2500
      ENDIF
C
C     TROP D'ITERATIONS SANS CONVERGENCE
      IERR = 1
      IF( LANGAG .EQ. 0 ) THEN
         WRITE (IMPRIM,1110) MXITER
         WRITE (IMPRIM,1100) (I,EIGV(I),I=1,NC)
      ELSE
         WRITE (IMPRIM,21110) MXITER
         WRITE (IMPRIM,21100) (I,EIGV(I),I=1,NC)
      ENDIF
      CALL IMPTAD( NC, NC, VEC )
      GOTO 117
C
C     MATRICES REDUITES A DES SCALAIRES NC=1
 1000 WRITE (IMPRIM,1100) (I,D(I),I = 1,NC)
      DO J = 1,NC
         VEC(J,1) = 0
      ENDDO
      VEC(1,1) = 1
      DO I = 1,NC
         EIGV(I) = D(I)
      ENDDO
C
      RETURN
      END
