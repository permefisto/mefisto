      SUBROUTINE PRL2P3( N, U, XY, RAYON, ANGTOT,
     %                   A, B, X, Y, AUX, TG )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    GENERER LES ARETES DE LA LIGNE ARC DE CERCLE
C -----    DEFINIE PAR 3 POINTS DE R**3 PAR PROJECTION L2 SUR
C          L'ESPACE DES POLYNOMES P3 HERMITE PAR MORCEAUX
C
C ENTREES:
C --------
C N      : NOMBRE DE POINTS DE L'ARC DE CERCLE
C U      : LES N VALEURS DU PARAMETRE U(1)=0 U(N)=1
C XY     : LES 2 COORDONNEES DES N SOMMETS DE L'ARC DE CERCLE EN 2D
C          CENTRE EN (0,0)
C RAYON  : RAYON DU CERCLE CENTRE EN (0,0)
C ANGTOT : ANGLE TOTAL DE L'ARC DE CERCLE A PARTIR DE L'AXE DES X
C
C AUXILIAIRES:
C ------------
C A      : SOUS DIAGONALE DE LA MATRICE QUI EST SYMETRIQUE
C B      : DIAGONALE DE LA MATRICE
C X      : PREMIER SECOND-MEMBRE (ABSCISSE)
C Y      : SECOND  SECOND-MEMBRE (ORDONNEE)
C AUX    : 2 VECTEURS DE N COMPOSANTES REELLES
C
C SORTIES:
C --------
C TG     : LES N COMPOSANTES ABSCISSE DES TANGENTES SUIVIES
C          DES N COMPOSANTES ORDONNEE DES TANGENTES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS    NOVEMBRE 1995
C2345X7..............................................................012
      REAL   U(N), XY(2,N), A(N), B(N), X(N), Y(N), AUX(N,2), TG(N,2)
C
C     CONSTRUCTION DE LA MATRICE SYMETRIQUE
C     LONGUEUR DE L'INTERVALLE 1 DU PARAMETRE
      HI = U(2) - U(1)
      AI = 0.0
      BI = ANGTOT * HI
C
C     LES INTEGRALES RIL DE 0 A 1 DE COS(A+Bt) * (t**L) dt
C     LES INTEGRALES RJL DE 0 A 1 DE SIN(A+Bt) * (t**L) dt
      COSA  = COS(AI)
      COSAB = COS(AI+BI)
      SINA  = SIN(AI)
      SINAB = SIN(AI+BI)
C
      RI0 = ( SINAB - SINA ) / BI
      RJ0 = ( COSA - COSAB ) / BI
C
      RI1 = ( SINAB - RJ0 ) / BI
      RJ1 = ( RI0 - COSAB ) / BI
C
      RI2 = ( SINAB - 2 * RJ1 ) / BI
      RJ2 = ( 2 * RI1 - COSAB ) / BI
C
      RI3 = ( SINAB - 3 * RJ2 ) / BI
      RJ3 = ( 3 * RI2 - COSAB ) / BI
C
C     PREMIER COEFFICIENT DIAGONAL
      B(1) = HI / 105.0
C
C     LES 2 SECONDS MEMBRES
      X(1) = HI * ( RAYON * ( RI1 - 2 * RI2 + RI3 )
     %              - 13.0/420.0 * XY(1,2) - 11.0/210.0 * XY(1,1) )
      Y(1) = HI * ( RAYON * ( RJ1 - 2 * RJ2 + RJ3 )
     %              - 13.0/420.0 * XY(2,2) - 11.0/210.0 * XY(2,1) )
C
C     ANGLE FINAL DE LA PREMIERE ARETE
      AI  = ANGTOT * HI
C
C     LES POINTS 2 A N-1
      DO 30 I=2,N-1
C
C        LONGUEUR DE L'INTERVALLE I DU PARAMETRE U
         HI1 = U(I+1) - U(I)
C
C        ANGLE DE L'ARETE I
         BI  = ANGTOT * HI1
C
C        LES INTEGRALES RIL DE 0 A 1 DE COS(A+Bt) * (t**L) dt
C        LES INTEGRALES RJL DE 0 A 1 DE SIN(A+Bt) * (t**L) dt
         COSA  = COS(AI)
         COSAB = COS(AI+BI)
         SINA  = SIN(AI)
         SINAB = SIN(AI+BI)
C
         RI10 = ( SINAB - SINA ) / BI
         RJ10 = ( COSA - COSAB ) / BI
C
         RI11 = ( SINAB - RJ10 ) / BI
         RJ11 = ( RI10 - COSAB ) / BI
C
         RI12 = ( SINAB - 2 * RJ11 ) / BI
         RJ12 = ( 2 * RI11 - COSAB ) / BI
C
         RI13 = ( SINAB - 3 * RJ12 ) / BI
         RJ13 = ( 3 * RI12 - COSAB ) / BI
C
         A(I) = - HI / 140.0
         B(I) = ( HI + HI1 ) / 105.0
C
C        LES 2 SECONDS MEMBRES
         X(I) = HI  * ( RAYON * ( RI3 - RI2 )
     %                 + 13.0/420.0 * XY(1,I-1) + 11.0/210.0 * XY(1,I) )
     %        + HI1 * ( RAYON * ( RI11 - 2 * RI12 + RI13 )
     %                 - 13.0/420.0 * XY(1,I+1) - 11.0/210.0 * XY(1,I) )
         Y(I) = HI  * ( RAYON * ( RJ3 - RJ2 )
     %                 + 13.0/420.0 * XY(2,I-1) + 11.0/210.0 * XY(2,I) )
     %        + HI1 * ( RAYON * ( RJ11 - 2 * RJ12 + RJ13 )
     %                 - 13.0/420.0 * XY(2,I+1) - 11.0/210.0 * XY(2,I) )
C
C        MISE A JOUR POUR LE SUIVANT
         HI  = HI1
         AI  = AI + BI
         RI2 = RI12
         RI3 = RI13
         RJ2 = RJ12
         RJ3 = RJ13
 30   CONTINUE
C
C     LA DERNIERE LIGNE DU SYSTEME
      A(N) = - HI / 140.0
      B(N) =   HI / 105.0
C
C     LES 2 SECONDS MEMBRES
      X(N) = HI  * ( RAYON * ( RI3 - RI2 )
     %                 + 13.0/420.0 * XY(1,N-1) + 11.0/210.0 * XY(1,N) )
      Y(N) = HI  * ( RAYON * ( RJ3 - RJ2 )
     %                 + 13.0/420.0 * XY(2,N-1) + 11.0/210.0 * XY(2,N) )
C
C     RESOLUTION DU SYSTEME TRIDIAGONAL SYMETRIQUE AVEC 2 SECONDS MEMBRES
      CALL TRDISY( N, A, B, X, AUX(1,1), AUX(1,2), TG(1,1) )
      CALL TRDISY( N, A, B, Y, AUX(1,1), AUX(1,2), TG(1,2) )
      END
