      SUBROUTINE NLSE2P1D( DT,     DTALFA, DTBETA, NBSOM, XYZSOM,
     %                     NBTRIA, NUELEM, NOSOTR,
     %                     NOOBSF, NUMISU, NUMASU, LTDESU,
     %                     VG,     WG,
     %                     tPP,    tDPDP,  tPU2P, FOmega )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DES INTEGRALES EXACTES DES POLYNOMES DE BASE POUR
C -----    UN TRIANGLE 2P1D LAGRANGE DE DEGRE 1 EXCEPTEE POUR TPP
C          OU L'INTEGRATION EST AUX SOMMETS POUR PRODUIRE UNE MATRICE
C          DIAGONALE
C
C ENTREES:
C --------
C DT     : PAS DE TEMPS EN DOUBLE PRECISION
C DTALFA : COEFFICIENT du LAPLACIEN dt * Alfa
C DTBETA : COEFFICIENT de U**3      dt * Beta
C NBSOM  : NOMBRE DE SOMMETS DU MAILLAGE
C XYZSOM : COORDONNEES DES NBSOM SOMMETS
C NBTRIA : NOMBRE DE TRIANGLES
C NUELEM : NUMERO DU TRIANGLE A TRAITER
C NOSOTR : NUMERO DES 3 SOMMETS DES NBTRIA TRIANGLES
C
C NOOBSF : NUMERO DE L'OBJET SURFACE DE CET ELEMENT
C NUMISU : NUMERO MINIMAL DES OBJETS SURFACES
C NUMASU : NUMERO MAXIMAL DES OBJETS SURFACES
C LTDESU : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES SURFACIQUES DE L'ONDE
C VG, WG : ONDE UG=PARTIE REELLE VG + i PARTIE IMAGINAIRE WG
C
C SORTIES:
C --------
C tPP    : Integrale t[P][P] dX = TPP
C          AVEC INTEGRATION AUX SOMMETS POUR PRODUIRE UNE MATRICE DIAGONALE
C tDPDP  : dt ALFA Integrale t[DP][DP] dX (pour -LAPLACIEN)
C tPU2P  : dt BETA Integrale ( (P vn)**2 + (P wn)**2 )  t[P][P] dX
C FOmega : Integrale t[P][P] dX FeRI(t,St,V,W) AUX 3 SOMMETS ET 2 PARTIES
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Texas A & M University at QATAR  Fevrier 2011
C23456---------------------------------------------------------------012
      include"./incl/donthe.inc"
      include"./incl/a___source.inc"
      include"./incl/ctemps.inc"
      include"./incl/cthet.inc"
      include"./incl/cnonlin.inc"
C
      INTEGER           LTDESU(1:MXDOTH, NUMISU:NUMASU)
      DOUBLE PRECISION  DT, DTALFA, DTBETA, VG(NBSOM), WG(NBSOM)
      REAL              XYZSOM(3,NBSOM)
      INTEGER           NOSOTR(NBTRIA,3)
      DOUBLE PRECISION  X21, Y21, X31, Y31, X32, Y32, DELTA, S, DELTS6,
     %                  XYZPI(3) ,FOmegaST(2)
      DOUBLE PRECISION  tPP(3), TDPDP(3,3), tPU2P(3), FOmega(3,2)
C
C     NUMERO DES 3 SOMMETS DU TRIANGLE
C     --------------------------------
      NS1 = NOSOTR( NUELEM, 1 )
      NS2 = NOSOTR( NUELEM, 2 )
      NS3 = NOSOTR( NUELEM, 3 )
C
C     DETERMINANT DE LA TRANSFORMATION F: e ref --> e
C     -----------------------------------------------
      X21 = XYZSOM(1,NS2) - XYZSOM(1,NS1)
      X31 = XYZSOM(1,NS3) - XYZSOM(1,NS1)
      X32 = XYZSOM(1,NS3) - XYZSOM(1,NS2)
C
      Y21 = XYZSOM(2,NS2) - XYZSOM(2,NS1)
      Y31 = XYZSOM(2,NS3) - XYZSOM(2,NS1)
      Y32 = XYZSOM(2,NS3) - XYZSOM(2,NS2)
C
      DELTA = ABS( X21 * Y31 - X31 * Y21 )
C
C     Integrale t[DP][DP] dX = t[DP] t[DF]-1 [DF]-1 [DP] / (2 DELTA)
C     --------------------------------------------------------------
      S = DTALFA / ( 2D0 * DELTA )
      TDPDP(1,1) = (  X32*X32 + Y32*Y32 ) * S
      TDPDP(2,1) = ( -X31*X32 - Y31*Y32 ) * S
      TDPDP(3,1) = (  X21*X32 + Y21*Y32 ) * S
C
      TDPDP(1,2) = TDPDP(2,1)
      TDPDP(2,2) = (  X31*X31 + Y31*Y31 ) * S
      TDPDP(3,2) = ( -X21*X31 - Y21*Y31 ) * S
C
      TDPDP(1,3) = TDPDP(3,1)
      TDPDP(2,3) = TDPDP(3,2)
      TDPDP(3,3) = (  X21*X21 + Y21*Y21 ) * S
C
C     Integrale ( (P vn)**2 + (P wn)**2 )  t[P][P] dX
C     -----------------------------------------------
      DELTS6 = DELTA / 6D0
      S      = DTBETA * DELTS6
      DO K=1,3
         NS = NOSOTR( NUELEM, K )
         tPU2P(K) = ( VG(NS)**2 + WG(NS)**2 ) * S
      ENDDO
C
C     Integrale t[P][P] dX = TPP MATRICE DIAGONALE
C     --------------------------------------------
      tPP(1) = DELTS6
      tPP(2) = DELTS6
      tPP(3) = DELTS6
C
C     Integrale t[P][P] dX FeRI(t,St,V,W)
C     -----------------------------------
      DELTS6 = DELTS6 * DT
      DO K=1,3
C
C        NUMERO DU SOMMET K DU TRIANGLE
         NS = NOSOTR( NUELEM, K )
C
C        LES 3 COORDONNEES DU SOMMET K
         DO I=1,3
            XYZPI(I) = XYZSOM( I, NS )
         ENDDO
C
C        VALEUR DE L'ONDE EN CE SOMMET
         TEMPEL = VG( NS )
         ONDEPI = WG( NS )
C
C        LES 2 PARTIES REELLE et IMAGINAIRE de FOmega AU SOMMET K
C        LA VALEUR DES SOURCES SURFACIQUES EN CE SOMMET K
         CALL REWAVE( 3, NOOBSF, 3, XYZPI,
     %                   LTDESU(LPSOUR,NOOBSF), FOmegaST )
C
C        LES 2 VECTEURS ELEMENTAIRES
         FOmega(K,1) = DELTS6 * FOmegaST(1)
         FOmega(K,2) = DELTS6 * FOmegaST(2)
C
      ENDDO
C
      RETURN
      END
