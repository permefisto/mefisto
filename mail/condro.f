      SUBROUTINE CONDRO( P, Q,  RAYON1, RAYON2, HAUTEU,
     %                   NBRAC, RACINE )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DES 0, 1, 2, INFINITE DE  RACINES OU
C -----    POINTS D'INTERSECTION DE LA DROITE PQ AVEC LE CONE SUPPOSE INFINI
C          L'ORIGINE DU REPERE PROPRE AU TRONC DE CONE EST LE CENTRE
C          DU CERCLE DE RAYON1 => Z=0 DANS CE PLAN
C
C PROGRAMMATION DOUBLE PRECISION SANS SEUIL EPSILON ... FAUT VOIR!
C
C ENTREES:
C --------
C P  Q   : POINTS DE DEFINITION DE LA DROITE
C RAYON1 : RAYON DU CERCLE AU POINT AX1 (PLAN ORTHOGONAL A L'AXE DU CONE )
C RAYON2 : RAYON DU CERCLE AU POINT AX2 (PLAN ORTHOGONAL A L'AXE DU CONE )
C HAUTEU : HAUTEUR DU TRONC DE CONE
C
C SORTIES:
C --------
C NBRAC  : NOMBRE DE RACINES OU POINTS D'INTERSECTION
C          0, 1 SIMPLE, 2 SIMPLES, 3 UNE DOUBLE,
C          4 <=> UNE INFINITE PQ EST TANGENTE GENERATRICE DU CONE
C RACINE : COEFFICIENT L TEL QUE LE POINT D'INTERSECTION M VERIFIE
C          PM = L * PQ  EGALITE VECTORIELLE
C          NBRAC=0 => RACINE NON INITIALISE
C          NBRAC=1 => RACINE(1)=L1 RACINE SIMPLE
C          NBRAC=2 => RACINE(1)=L1 RACINE(2)=L2  LES 2 RACINES SIMPLES
C          NBRAC=3 => RACINE(1)=L1 RACINE(2)=L1  RACINE DOUBLE
C          NBRAC=4 => RACINE NON INITIALISE CAR PQ GENERATRICE DU CONE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS    DECEMBRE 1997
C2345X7..............................................................012
      DOUBLE PRECISION  P(3), Q(3), RAYON1, RAYON2, HAUTEU, RACINE(2)
      DOUBLE PRECISION  XPQ, YPQ, ZPQ, R, S, A, B, C
C
      XPQ = Q(1) - P(1)
      YPQ = Q(2) - P(2)
      ZPQ = Q(3) - P(3)
C
C     DETECTION DE PQ GENERATRICE DU CONE
      IF( P(3) .NE. Q(3) ) THEN
C
C        POINT DE LA DROITE PQ ET DU PLAN Z=0
         R = P(3) / ZPQ
         S = ( P(1) + R * XPQ ) ** 2
     %     + ( P(2) + R * YPQ ) ** 2
C
C        POINT DE LA DROITE PQ ET DU PLAN Z=HAUTEU
         R = ( P(3) - HAUTEU ) / ZPQ
         R = ( P(1) + R * XPQ ) ** 2 + ( P(2) + R * YPQ ) ** 2
C
         IF( S .EQ. RAYON1 ** 2 .AND. R .EQ. RAYON2 ** 2 ) THEN
C           LA DROITE PQ EST UNE GENERATRICE DU CONE
            NBRAC = 4
            RETURN
         ENDIF
      ENDIF
C
C     ICI PQ N'EST PAS UNE GENERATRICE DU CONE
C     CALCUL DES 3 COEFFICIENTS DU TRINOME
      R = ( RAYON2 - RAYON1 ) / HAUTEU
      S = RAYON1 + R * P(3)
C
      A = XPQ ** 2 + YPQ ** 2 - ( R * ZPQ ) ** 2
      B = P(1) * XPQ + P(2) * YPQ - S * R * ZPQ
      C = P(1) ** 2 + P(2) ** 2 - S ** 2
C
      IF( A .EQ. 0D0 ) THEN
C
C        UNE RACINE SIMPLE
         NBRAC = 1
         RACINE(1) = - C / ( 2 * B )
         RETURN
C
      ENDIF
C
C     IL Y A 0 OU 2 RACINES (SIMPLES OU DOUBLE)
      S = B ** 2 - A * C
      IF( S .LT. 0D0 ) THEN
C
C        0 RACINE
         NBRAC = 0
         RETURN
C
      ELSE IF( S .EQ. 0D0 ) THEN
C
C        1 RACINE DOUBLE
         NBRAC = 3
         RACINE(1) = - B / A
         RACINE(2) = RACINE(1)
C
      ELSE
C
C        2 RACINES SIMPLES
         NBRAC = 2
         S = SQRT( S )
C        CELLE DE PLUS GRAND MODULE EST EN POSITION 1
         IF( B .GE. 0D0 ) THEN
            RACINE(1) = - ( B + S ) / A
            RACINE(2) = - 2D0 * B / A - RACINE(1)
         ELSE
            RACINE(1) = ( S - B ) / A
            RACINE(2) = - 2D0 * B / A - RACINE(1)
         ENDIF
C
      ENDIF
      RETURN
      END
