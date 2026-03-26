      SUBROUTINE INTARAR3( S1, S2, P1, P2, NBPTI, PTI, MUP )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     DETERMINER LE POINT D'INTERSECTION DES ARETES S1S2 ET P1P2
C -----     D'UN MEME PLAN DE R3 EN NE RETENANT LE POINT SI ET SEULEMENT
C           S'IL EST INTERNE AUX 2 ARETES (EXTREMITES COMPRISES)
C           PTI EST LE POINT LE PLUS PROCHE DES 2 DROITES CALCULE
C           PAR MOINDRES CARRES DE LA DISTANCE AVEC DROITES PARAMETREES
C           (Autre methode dans mail/mindis.f)
C
C ENTREES:
C --------
C S1,S2  : LES 2 POINTS QUI DEFINISSENT L'ARETE 1
C P1,P2  : LES 2 POINTS QUI DEFINISSENT L'ARETE 2
C
C SORTIES:
C --------
C NBPTI  : LE NOMBRE DE POINT D'INTERSECTION DES ARETES S1S2 ET P1P2
C          =0 SI LE POINT D'INTERSECTION N'EXISTE PAS OU EST EXTERNE
C             STRICTEMENT A L'UNE DES 2 ARETES
C          =1 SI LE POINT D'INTERSECTION EXISTE ET EST INTERNE AUX 2
C             ARETES (EXTREMITES INCLUSES)
C PTI    : LES 3 COORDONNEES DU POINT D'INTERSECTION INTERNE
C          (PEUT ETRE L'UN DES 4 SOMMETS DES 2 ARETES)
C MUP    : COORDONNEE BARYCENTRIQUE DE PTI SUR LE VECTEUR P1P2
C          MUP=0 PTI=P1;   MUP=1 PTI=P2;
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: PERRONNET Alain LJLL UPMC & St Pierre du Perray SEPTEMBRE 2011
C2345X7..............................................................012
      DOUBLE PRECISION  S1(3), S2(3), P1(3), P2(3), PTI(3)
      INTEGER           NBPTI
      DOUBLE PRECISION  A, B, C, D, XS21, YS21, ZS21,  XP21, YP21, ZP21,
     %                  XPS1, YPS1, ZPS1, MUS, MUS1, MUP, MUP1, DENOM,
     %                  NRS1S2, NRP1P2
      INTRINSIC         ABS, SQRT
C
      MUP = -1D10
C
C     VECTEUR S1P1
      XPS1 = P1(1) - S1(1)
      YPS1 = P1(2) - S1(2)
      ZPS1 = P1(3) - S1(3)
C
C     VECTEUR S1S2 ET SA NORME
      XS21 = S2(1) - S1(1)
      YS21 = S2(2) - S1(2)
      ZS21 = S2(3) - S1(3)
      A    = XS21**2   + YS21**2   + ZS21**2
C     NORME S1S2
      NRS1S2 = SQRT( A )
C
C     VECTEUR P1P2 ET SA NORME
      XP21 = P2(1) - P1(1)
      YP21 = P2(2) - P1(2)
      ZP21 = P2(3) - P1(3)
      D    = XP21 * XP21 + YP21 * YP21 + ZP21 * ZP21
C     NORME P1P2
      NRP1P2 = SQRT( D )
C
C     (S1S2,P1P2)
      B = XS21 * XP21 + YS21 * YP21 + ZS21 * ZP21
C
C     (S1S2,S1P1)
      C = XS21 * XPS1 + YS21 * YPS1 + ZS21 * ZPS1
C
C     LE DENOMINATEUR DU PARAMETRE MUP DE LA DROITE P1P2
      DENOM = B * B - A * D
C
      IF( ABS(DENOM) .LE. NRS1S2 * NRP1P2 * 1D-4 ) THEN
C        S1S2 // P1P2
         GOTO 9000
      ENDIF
C
C     LE PARAMETRE MUP DE LA DROITE P1P2
      MUP = ( XP21 * ( A * XPS1 - C * XS21 )
     %      + YP21 * ( A * YPS1 - C * YS21 )
     %      + ZP21 * ( A * ZPS1 - C * ZS21 ) ) / DENOM
C
C     LE POINT D'INTERSECTION EST IL INTERNE AU SEGMENT P1P2?
      IF( MUP .LT. 0D0 .OR. MUP .GT. 1D0 ) THEN
C        PTI EST EXTERIEUR AU SEGMENT P1P2
         GOTO 9000
      ENDIF
C
C     LE PARAMETRE MUS DE LA DROITE S1S2 DU POINT D'INTERSECTION
      MUS = ( MUP * B + C ) / A
      IF( MUS .LT. 0D0 .OR. MUS .GT. 1D0 ) THEN
C        PTI EST EXTERIEUR AU SEGMENT S1S2
         GOTO 9000
      ENDIF
C
C     LE POINT D'INTERSECTION PTI =( P1 + MUP * ( P2 - P1 )
C                                  + S1 + MUS * ( S2 - S1 ) ) / 2
C     COMME MILIEU DU SEGMENT DES 2 POINTS DES 2 DROITES
C     SI ELLES NE SONT PAS DANS LE MEME PLAN
      MUS1   = 1D0 - MUS
      MUP1   = 1D0 - MUP
      NBPTI  = 1
      PTI(1) = (MUS1 * S1(1) + MUS * S2(1) + MUP1 * P1(1) + MUP * P2(1))
     %         / 2D0
      PTI(2) = (MUS1 * S1(2) + MUS * S2(2) + MUP1 * P1(2) + MUP * P2(2))
     %         / 2D0
      PTI(3) = (MUS1 * S1(3) + MUS * S2(3) + MUP1 * P1(3) + MUP * P2(3))
     %         / 2D0
C
ccc      print *,'intarar3: DENOM=',DENOM,' MUP=',MUP,' MUS=',MUS,
ccc     %        ' voir si difference...'
ccc      print *,'intarar3: PTI =',(PTI(K),K=1,3)
ccc      print *,'intarar3: PTIS=',(MUS1 * S1(K) + MUS * S2(K),K=1,3)
ccc      print *,'intarar3: PTIP=',(MUP1 * P1(K) + MUP * P2(K),K=1,3)
ccc
      GOTO 9900
C
C     PAS DE POINT D'INTERSECTION PTI
 9000 NBPTI = 0
C
 9900 RETURN
      END
