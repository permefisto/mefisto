      SUBROUTINE BSPINO( KDEGRE, LU, LR, LT,
     %                   U, T, R, B, XYZ, A, S, FACM, MATRIX,
     %                   NORT, IERR )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    GENERER LES COEFFICIENTS DES POLYNOMES SUR CHAQUE INTERVALLE
C -----    D'UNE B-SPLINE OUVERTE D'INTERPOLATION
C
C ENTREES:
C --------
C KDEGRE : DEGRE DES POLYNOMES DE LA B-SPLINE
C LU     : NOMBRE DE POINTS D'INTERPOLATION
C LT     : NOMBRE DE POINTS-1 DE CONTROLE (T) DE LA LIGNE
C          OU COEFFICIENTS OBTENUS APRES INVERSION DU SYSTEME LINEAIRE
C NBNOET : NOMBRE DE NOEUDS T
C U      : LES ABSCISSES PARAMETRE AYANT POUR IMAGE LES POINTS
C          D'INTERPOLATION
C XYZ    : 3 COORDONNEES DES POINTS DE CONTROLE
C
C AUXILIAIRES :
C -------------
C T      : VALEURS DU PARAMETRE POUR CALCULER LES B(J,M)
C B      : VALEURS DES B(J,M)(TI)
C A      : VALEURS INTERMEDIAIRES
C FACM   : 1/M!
C MATRIX : LA MATRICE DES BJ,M(RI) A INVERSER
C NORT   : NUMERO DU DERNIER NOEUD T DE CHAQUE POINT R
C
C SORTIES:
C --------
C R      : LE PARAMETRE DES POINTS ASSOCIES A T
C S      : LES COEFFICIENTS DES POLYNOMES SUR CHAQUE INTERVALLE
C IERR   : =0 SI PAS D'ERREUR
C          >0 SI ERREUR RENCONTREE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS       JUIN  1990
C2345X7..............................................................012
      include"./incl/gsmenu.inc"
      INTEGER  NORT(0:LR)
      REAL     U(0:LU),R(0:LR),
     %         T(0:LT+KDEGRE),
     %         B(-KDEGRE:1,0:KDEGRE),
     %         A(-KDEGRE:0,0:KDEGRE),
     %         XYZ(1:3,0:LU-1),
     %         FACM(0:KDEGRE),
     %         S(0:KDEGRE,0:LR-1,1:3)
      DOUBLE PRECISION SS, MATRIX(LU,LU+3)
C
C     LE TABLEAU DES 1/M!
      FACM(0) = 1.0
      FACM(1) = 1.0
      SS      = 1.D0
      DO 10 M=2,KDEGRE
        SS = SS * M
        FACM(M) = REAL( 1.D0 / SS )
 10   CONTINUE
C
C     LE CALCUL DES NOEUDS T ET R NORT A PARTIR DES U
C     ===============================================
      CALL BSPLS1( KDEGRE, LU, LT, LR, U, T,  R, NORT )
C
C     GENERATION DE LA MATRICE DES BJ,M(RI)
C     =====================================
      LB = (KDEGRE+2)*(KDEGRE+1)
C
C     MISE A ZERO DE LA MATRICE
      CALL AZEROD( LU*LU, MATRIX )
C
      MATRIX(1,1) = 1.0D0
      DO 57 I=1,LU-2
C
C        QUEL EST L'INTERVALLE [RJ,RJ+1[ CONTENANT U(I) ?
         DO 20 J=0,LR-1
            IF( U(I) .GE. R(J) .AND. U(I) .LT. R(J+1) ) GOTO 30
 20      CONTINUE
         J = LR
C
C        LE NUMERO DU DERNIER NOEUD T AU POINT R(J)
 30      NOI = NORT( J )
C
C        CALCUL DES BJM(RI)
         CALL AZEROR( LB, B )
         B(0,0) = 1.0
C
         DO 50 M=1,KDEGRE
            DO 48 J=-M,0
C              LE VRAI INDICE DE B
               JB = NOI + J
C              EVALUATION DES FRACTIONS DE LA FORMULE
               IF( T(JB+M) .NE. T(JB) ) THEN
                  U1 = ( U(I) - T(JB) ) / ( T(JB+M) - T(JB) )
               ELSE
                  U1 = 0
               ENDIF
               IF( T(JB+M+1) .NE. T(JB+1) ) THEN
                  U2 = ( T(JB+M+1) - U(I) )
     %               / ( T(JB+M+1) - T(JB+1) )
               ELSE
                  U2 = 0
               ENDIF
               B(J,M) = U1 * B(J,M-1) + U2 * B(J+1,M-1)
 48         CONTINUE
 50      CONTINUE
C
C        RANGEMENT DE B(RI) DANS LA LIGNE I+1 DE LA MATRICE MATRIX
         DO 55 J=-KDEGRE,0
            MATRIX(1+I,1+NOI+J) = B(J,KDEGRE)
 55      CONTINUE
 57   CONTINUE
      MATRIX(LU,LU) = 1.0D0
C
C     GENERATION DES 3 SECONDS MEMBRES
C     ================================
      DO 60 I=1,LU
         DO 58 NC=1,3
            MATRIX(I,LU+NC) = XYZ(NC,I-1)
 58      CONTINUE
 60   CONTINUE
C
C     INVERSION DU SYSTEME LINEAIRE
C     =============================
      CALL GAUSPT( LU, 3, MATRIX, I )
      IF( I .NE. 0 ) THEN
         NBLGRC(NRERR) = 2
         KERR(1) = 'MAUVAIS CHOIX DES PARAMETRES'
         KERR(2) = 'ASSOCIES AUX POINTS'
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF
C
C     GENERATION DES POLYNOMES DE LA B-SPLINE
C     =======================================
      DO 1000 I=0,LR-1
C
C        LA VALEUR DU PARAMETRE
         SS  = R(I)
C        LE NUMERO DU DERNIER NOEUD DU POINT I
         NOI = NORT( I )
C
C        CALCUL DES BJM(RI)
         CALL AZEROR( LB, B )
         B(0,0) = 1.0
C
         DO 90 M=1,KDEGRE
            DO 88 J=-M,0
C              LE VRAI INDICE DE B
               JB = NOI + J
C              EVALUATION DES FRACTIONS DE LA FORMULE
               IF( T(JB+M) .NE. T(JB) ) THEN
                  U1 = REAL( ( SS - T(JB) ) / ( T(JB+M) - T(JB) ) )
               ELSE
                  U1 = 0
               ENDIF
               IF( T(JB+M+1) .NE. T(JB+1) ) THEN
                  U2 = REAL( ( T(JB+M+1) - SS )
     %                     / ( T(JB+M+1) - T(JB+1) ) )
               ELSE
                  U2 = 0
               ENDIF
               B(J,M) = U1 * B(J,M-1) + U2 * B(J+1,M-1)
 88         CONTINUE
 90      CONTINUE
C
C        CALCULS POUR CHAQUE COMPOSANTE DES POINTS
         DO 900 NC=1,3
C
C           INITIALISATION
            DO 120 J=-KDEGRE,0
               A(J,0) = REAL( MATRIX( 1 + NOI + J, LU + NC ) )
 120        CONTINUE
C
C           LE CALCUL PROPREMENT DIT
            DO 240 M=1,KDEGRE
               KM = KDEGRE - M + 1
               DO 230 J=-KDEGRE+M,0
C                 LE VRAI INDICE DE A
                  JA = NOI + J
                  IF( T(JA+KM) .NE. T(JA) ) THEN
                     A(J,M) = KM*(A(J,M-1)-A(J-1,M-1))/(T(JA+KM)-T(JA))
                  ELSE
                     A(J,M) = 0
                  ENDIF
 230           CONTINUE
 240        CONTINUE
C
C           LES COEFFICIENTS S(M,I,NC) DU POLYNOME
            DO 320 M=0,KDEGRE
               SS = 0
               DO 310 J=-KDEGRE+M,0
                  SS = SS + A(J,M) * B(J,KDEGRE-M)
 310           CONTINUE
               S(M,I,NC) = REAL( FACM(M) * SS )
 320        CONTINUE
C
 900     CONTINUE
 1000 CONTINUE
C
      RETURN
      END
