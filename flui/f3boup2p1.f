      SUBROUTINE F3BOUP2P1( XYZNOE, NONOTE, TP2P2,
     %                      DELTAT, Rho, G, CoBOUS, TEMPER, TEMPER0,
     %                      VE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CONTRIBUTION AU SECOND MEMBRE DU TETRAEDRE TAYLOR HOOD
C -----    P2 CONTINU POUR LES 3 COMPOSANTES DE LA VITESSE
C          P2 CONTINU POUR LA TEMPERATURE
C         -Integrale tP2 P2 dx dt Rho G (-Z) CoBOUS (Temper(tn+1)-Temper(t0))
C          La pesanteur G est dirigee vers l'axe -Z

C ENTREES:
C --------
C XYZNOE : 3 COORDONNEES DES NOEUDS DES TETRAEDRES
C NONOTE : NONOTE(I) NO GLOBAL DU I-EME NOEUD DU TETRAEDRE I=1,...,10
C TP2P2  : Integrale tP2 P2 dX sur LE TETRAEDRE UNITE

C DELTAT : PAS DE TEMPS POUR L'INTEGRATION EN TEMPS
C Rho    : DENSITE DE MASSE  Integrale Rho P2 P2 dX
C G      : ACCELERATION DE LA PESANTEUR TERRESTRE SELON l'AXE -Z
C CoBOUS : COEFFICIENT dans l'APPROXIMATION DE BOUSSINESQ
C          Rho0 - Rho = Rho CoBOUS (Temper-Temper0)
C          CoBOUS est le COEFFICIENT DE DILATATION THERMIQUE du FLUIDE
C TEMPER : TEMPERATURE au TEMPS tn+1 DES NOEUDS du MAILLAGE
C TEMPER0: TEMPERATURE au TEMPS tn   DES NOEUDS du MAILLAGE

C SORTIE :
C --------
C VE     : DL ELEMENTAIRES DU SECOND MEMBRE EN VITESSE(1:10,1:3)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET  Saint Pierre du Perray               Mai 2022
C23456---------------------------------------------------------------012
      IMPLICIT  NONE
      include"./incl/langue.inc"
      include"./incl/donflu.inc"
      include"./incl/ctemps.inc"
      INTEGER         LECTEU, IMPRIM, NUNITE
      COMMON /UNITES/ LECTEU, IMPRIM, NUNITE(30)
C
      REAL              XYZNOE(3,*)
      DOUBLE PRECISION  TP2P2(10,10), TEMPER(*), TEMPER0(*)
      INTEGER           NONOTE(10), I, J, N
      DOUBLE PRECISION  DF(3,3), DELTAT, Rho, G, CoBOUS
      DOUBLE PRECISION  VE(10,3)
      DOUBLE PRECISION  DELTAe, X1, Y1, Z1, D, S
      INTRINSIC         ABS

      IF( CoBOUS .EQ. 0D0 ) GOTO 9999

C     CALCUL DE LA MATRICE JACOBIENNE ET DE SON DETERMINANT
      N = NONOTE(1)
      X1 = XYZNOE(1,N)
      Y1 = XYZNOE(2,N)
      Z1 = XYZNOE(3,N)
      DO I=1,3
         N = NONOTE(I+1)
         DF(I,1) = XYZNOE(1,N) - X1
         DF(I,2) = XYZNOE(2,N) - Y1
         DF(I,3) = XYZNOE(3,N) - Z1
      ENDDO

C     LE DETERMINANT DE DF
      DELTAe = ABS( DF(1,1) * ( DF(2,2) * DF(3,3) - DF(3,2) * DF(2,3) )
     %            + DF(2,1) * ( DF(3,2) * DF(1,3) - DF(1,2) * DF(3,3) )
     %            + DF(3,1) * ( DF(1,2) * DF(2,3) - DF(2,2) * DF(1,3) ))
      IF( DELTAe .LE. 0D0 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            PRINT *, 'f3ex2p2p1: ATTENTION EF de VOLUME*6=',
     %                DELTAe,' NON PRIS EN COMPTE'
         ELSE
            PRINT *, 'f3ex2p2p1: ATTENTION FE of VOLUME*6=',
     %                DELTAe,' is NOT TAKEN in ACCOUNT'
         ENDIF
         GOTO 9999
      ENDIF
C     ICI LE TETRAEDRE EST DE VOLUME NON NUL

C    -Integrale tP2 P2 dx  dt Rho G CoBOUS (Temper(tn+1)-Temper(t0))
C     g SELON L'AXE -Z   3-EME COMPOSANTE DE VE
      D = DELTAe * DELTAT * Rho * G
      DO I=1,10
         DO J=1,10
            N = NONOTE(J)
            S = CoBOUS * ( TEMPER(N) - TEMPER0(N) )
cccC           COMPOSANTE X DE LA VITESSE
ccc            VE(I,1) = VE(I,1) + TP2P2(I,J) * D * S
cccC           COMPOSANTE Y DE LA VITESSE
ccc            VE(I,2) = VE(I,2) + TP2P2(I,J) * D * S
C           COMPOSANTE Z DE LA VITESSE
            VE(I,3) = VE(I,3) + TP2P2(I,J) * D * S
         ENDDO
      ENDDO

 9999 RETURN
      END
