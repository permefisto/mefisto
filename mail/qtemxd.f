      SUBROUTINE QTEMXD( NS,     D,
     %                   XYZSOM, NBSOTE, NOSOTE, NO1TSO, NOTESO,
     %                   VOLUMT, QUALIT,
     %                   VOLUNS, QUALNS, NTQMIN, NBTENS )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++012
C BUT :    CALCULER LE COEFFICIENT MULTIPLICATEUR DE D LA DIRECTION
C -----    A PARTIR DU POINT P POUR OBTENIR LE POINT QUI MAXIMISE
C          LA QUALITE DES TETRAEDRES DE SOMMET NS

C ENTREES:
C --------
C NS     : NUMERO DU SOMMET DE QUALITE A AMELIORER
C D      : 3 COORDONNEES DE LA DIRECTION DE LA DROITE A PARTIR DE P
C NS     : NUMERO DU SOMMET DE QUALITE A AMELIORER
C XYZSOM : COORDONNEES X Y Z DES NBSOMM SOMMETS DES TETRAEDRES
C NBSOTE : VALEUR MAXIMALE DE DECLARATION DU PREMIER INDICE DE NOSOTE(>3)
C NOSOTE : NUMERO DES 4 SOMMETS DE CHAQUE TETRAEDRE
C NO1TSO : NUMERO DU 1-ER TETRAEDRE DANS NOTESO DE CHAQUE SOMMET
C NOTESO : NOTESO(1,*) NUMERO DANS NOSOTE DU TETRAEDRE
C          NOTESO(2,*) NUMERO DANS NOTESO DU TETRAEDRE SUIVANT
C                      0 SI C'EST LE DERNIER
C VOLUMT : VOLUME  DES TETRAEDRES DE LA TETRAEDRISATION
C QUALIT : QUALITE DES TETRAEDRES DE LA TETRAEDRISATION
C
C SORTIES:
C --------
C VOLUNS : VOLUME DES TETRAEDRES DE SOMMET NS APRES MAXIMISATION
C QUALNS : QUALITE DU SOMMET NS DE LA TETRAEDRISATION APRES MAXIMISATION
C NTQMIN : NUMERO NOSOTE DU TETRAEDRE DE SOMMET NS ET DE QUALITE MINIMALE
C NBTENS : NOMBRE DE TETRAEDRES DE SOMMET NS

C REMARQUE: AU PIRE EN SORTIE LES 3 COORDONNEES DE NS N'ONT PAS CHANGE
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++012
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC  PARIS   NOVEMBRE 1993
C2345X7..............................................................012
      PARAMETER        (EPS=0.001,ITEMAX=8)
      INTEGER           NOSOTE(NBSOTE,*),
     %                  NO1TSO(*),
     %                  NOTESO(2,*)
      REAL              XYZSOM(3,*),
     %                  VOLUMT(*),
     %                  QUALIT(*)
      REAL              P(3),
     %                  D(3)

C     PARAMETRES INITIAUX POUR OBTENIR   T1 < T2 < T3
C                                    QUALITE(T1)<QUALITE(T2)
C                                    QUALITE(T3)<QUALITE(T2)
      EPS2 = EPS * EPS
C     PROTECTION DES XYZ DU SOMMET INITIAL
      P(1) = XYZSOM(1,NS)
      P(2) = XYZSOM(2,NS)
      P(3) = XYZSOM(3,NS)

      T1 = -0.6
      T2 =  0.0
      T3 =  0.6

C     QUALITE DU POINT NS = POINT 2 CENSE PRES DU MAXIMUM
      CALL QUALST( NS, XYZSOM, NBSOTE, NOSOTE, NO1TSO, NOTESO,
     %             VOLUMT, QUALIT,
     %             VOLUNS, Q2, NTQMIN, NBTENS )
      IF( NBTENS .LE. 0 ) GOTO 9999
      Q20 = Q2 * EPS

      ITER = 0

 1    ITER = ITER + 1
      XYZSOM(1,NS) = P(1) + T1 * D(1)
      XYZSOM(2,NS) = P(2) + T1 * D(2)
      XYZSOM(3,NS) = P(3) + T1 * D(3)
      CALL QUALST( NS, XYZSOM, NBSOTE, NOSOTE, NO1TSO, NOTESO,
     %             VOLUMT, QUALIT,
     %             VOLUNS, Q1, NTQMIN, NBTENS )
      IF( NBTENS .LE. 0 ) GOTO 9999
      IF( Q1 .LE. 0 ) THEN
         IF( ITER .GT. 4 ) THEN
            T2 = 0
            GOTO 50
         ENDIF
         T1 = T1 * 0.5
         GOTO 1
      ENDIF

      ITER = 0

 3    ITER = ITER + 1
      XYZSOM(1,NS) = P(1) + T3 * D(1)
      XYZSOM(2,NS) = P(2) + T3 * D(2)
      XYZSOM(3,NS) = P(3) + T3 * D(3)
      CALL QUALST( NS, XYZSOM, NBSOTE, NOSOTE, NO1TSO, NOTESO,
     %             VOLUMT, QUALIT,
     %             VOLUNS, Q3, NTQMIN, NBTENS )
      IF( NBTENS .LE. 0 ) GOTO 9999
      IF( Q3 .LE. 0 ) THEN
         IF( ITER .GT. 4 ) THEN
            T2 = 0
            GOTO 50
         ENDIF
         T3 = T3 * 0.5
         GOTO 3
      ENDIF
C
      Q4   = 0.0
      T4   = 0
      ITER = 0
C
 10   IF( Q2 .LT. Q1 .OR. Q2 .LT. Q3 ) GOTO 20
      ITER = ITER + 1
      IF( ITER .LE. ITEMAX .AND.
     %    ABS(T2-T1) .GT. EPS .AND. ABS(T3-T2) .GT. EPS .AND.
     %    ABS(Q2-Q1) .GT. Q20 .AND. ABS(Q3-Q2) .GT. Q20 ) THEN

C        PARABOLE PASSANT PAR (T1,QUALITE(P+T1*D)),
C                             (T2,QUALITE(P+T2*D)),
C                             (T3,QUALITE(P+T3*D))
         D1 = (T1-T2) * (T1-T3)
         D2 = (T2-T3) * (T2-T1)
         D3 = (T3-T1) * (T3-T2)

         A  =    Q1/D1         +   Q2/D2       + Q3/D3
         IF( ABS(A) .LT. EPS2 ) GOTO 20
         B  = -( Q1*(T2+T3)/D1 + Q2*(T3+T1)/D2 + Q3*(T1+T2)/D3 )

         T4 = -B / ( 2 * A )
         XYZSOM(1,NS) = P(1) + T4 * D(1)
         XYZSOM(2,NS) = P(2) + T4 * D(2)
         XYZSOM(3,NS) = P(3) + T4 * D(3)
         CALL QUALST( NS,XYZSOM,NBSOTE,NOSOTE,NO1TSO,NOTESO,
     %                VOLUMT, QUALIT,
     %                VOLUNS, Q4, NTQMIN, NBTENS )
         IF( NBTENS .LE. 0 ) GOTO 9999

C        TEST DE MISE A JOUR
         IF( T2 .LT. T4 .AND. T4 .LT. T3 ) THEN
            T1 = T2
            T2 = T4
            Q1 = Q2
            Q2 = Q4
         ELSE IF( T1 .LT. T4 .AND. T4 .LT. T2 ) THEN
            T3 = T2
            T2 = T4
            Q3 = Q2
            Q2 = Q4
         ELSE IF( T4 .LE. T1 ) THEN
            IF( Q2 .GT. Q4 .AND. Q2 .GT. Q3 ) THEN
               T1 = T4
               Q1 = Q4
            ELSE
C              MAX EN DEHORS DE L'INTERVALLE
               GOTO 20
            ENDIF
         ELSE IF( T4 .GE. T3 ) THEN
            IF( Q2 .GT. Q4 .AND. Q2 .GT. Q1 ) THEN
               T3 = T4
               Q3 = Q4
            ELSE
C              MAX EN DEHORS DE L'INTERVALLE
               GOTO 20
            ENDIF
         ELSE
C           T4 = T2
            GOTO 20
         ENDIF
         GOTO 10
      ENDIF

C     CONVERGENCE RECHERCHE DU MAX Q2 DE LA QUALITE PARMI LES 4 POINTS
 20   IF( Q1 .GT. Q2 ) THEN
         Q2 = Q1
         T2 = T1
      ENDIF
      IF( Q3 .GT. Q2 ) THEN
         Q2 = Q3
         T2 = T3
      ENDIF
      IF( Q4 .GT. Q2 ) THEN
         T2 = T4
      ENDIF

C     LE POINT FINAL: VOLUME ET QUALITE
 50   XYZSOM(1,NS) = P(1) + T2 * D(1)
      XYZSOM(2,NS) = P(2) + T2 * D(2)
      XYZSOM(3,NS) = P(3) + T2 * D(3)
      CALL QUALST( NS, XYZSOM, NBSOTE, NOSOTE, NO1TSO, NOTESO,
     %             VOLUMT, QUALIT,
     %             VOLUNS, QUALNS, NTQMIN, NBTENS )

 9999 RETURN
      END
