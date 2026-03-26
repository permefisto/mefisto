      SUBROUTINE F2BOUP2P1( XYZNOE, NONOTR, TP2P2,
     %                      DELTAT, Rho, G, CoBOUS, TEMPER, TEMPER0,
     %                      VE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CONTRIBUTION AU SECOND MEMBRE DU TRIANGLE TAYLOR HOOD
C -----    P2 CONTINU POUR LES 2 COMPOSANTES DE LA VITESSE
C          P2 CONTINU POUR LA TEMPERATURE
C         -Integrale tP2 P2 dx dt Rho (-Y) G CoBOUS (Temper(tn+1)-Temper(t0))
C          La pesanteur G est dirigee vers l'axe -Y

C ENTREES:
C --------
C XYZNOE : 3 COORDONNEES DES NOEUDS DES TRIANGLES
C NONOTR : NONOTR(I) NO GLOBAL DU I-EME NOEUD DU TRIANGLE I=1,...,6
C TP2P2  : Integrale tP2 P2 dX sur LE TRIANGLE UNITE

C DELTAT : PAS DE TEMPS POUR L'INTEGRATION EN TEMPS
C Rho    : DENSITE DE MASSE  Integrale Rho P2 P2 dX
C G      : ACCELERATION DE LA PESANTEUR TERRESTRE SELON l'AXE -Y
C CoBOUS : COEFFICIENT dans l'APPROXIMATION DE BOUSSINESQ
C          Rho0 - Rho = Rho CoBOUS (Temper-Temper0)
C          CoBOUS est le COEFFICIENT DE DILATATION THERMIQUE du FLUIDE
C TEMPER : TEMPERATURE au TEMPS tn+1 DES NOEUDS du MAILLAGE
C TEMPER0: TEMPERATURE au TEMPS tn   DES NOEUDS du MAILLAGE

C SORTIE :
C --------
C VE     : DL ELEMENTAIRES DU SECOND MEMBRE EN VITESSE(1:6,1:2)
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
      DOUBLE PRECISION  TP2P2(6,6), TEMPER(*), TEMPER0(*)
      INTEGER           NONOTR(6), I, J, M
      DOUBLE PRECISION  DELTAT, Rho, G, CoBOUS
      DOUBLE PRECISION  VE(6,2)
      DOUBLE PRECISION  DELTAe, X21, X31, Y21, Y31
      DOUBLE PRECISION  X1, Y1, D, S
      INTRINSIC         ABS

      IF( CoBOUS .EQ. 0D0 ) GOTO 9999

C     CALCUL DE LA MATRICE JACOBIENNE ET DE SON DETERMINANT
      M = NONOTR(1)
      X1  = XYZNOE(1,M)
      X21 = XYZNOE(1,NONOTR(2)) - X1
      X31 = XYZNOE(1,NONOTR(3)) - X1

      Y1  = XYZNOE(2,M)
      Y21 = XYZNOE(2,NONOTR(2)) - Y1
      Y31 = XYZNOE(2,NONOTR(3)) - Y1

C     CALCUL DU DETERMINANT DE LA MATRICE JACOBIENNE
      DELTAe = ABS( X21*Y31 - X31*Y21 )
      IF( DELTAe .LE. 0D0 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            PRINT *,'f2boup2p1: ATTENTION EF',
     %              ' de SURFACE*2=',DELTAe,' NON PRIS EN COMPTE'
         ELSE
            PRINT *,'f2boup2p1: ATTENTION FE',
     %              ' of SURFACE*2=',DELTAe,' NOT USED'
         ENDIF
         GOTO 9999
      ENDIF

C    -Integrale tP2 P2 dx  dt Rho G CoBOUS (Temper(tn+1)-Temper(t0))
C     g SELON L'AXE -Y   2-EME COMPOSANTE de VE
      D = DELTAe * DELTAT * Rho * G
      DO I=1,6
         DO J=1,6
            M = NONOTR(J)
            S = D * CoBOUS * ( TEMPER(M) - TEMPER0(M) )
cccC           COMPOSANTE X DE LA VITESSE
ccc            VE(I,1) = VE(I,1) + TP2P2(I,J) * S
C           COMPOSANTE Y DE LA VITESSE
            VE(I,2) = VE(I,2) + TP2P2(I,J) * S
         ENDDO
      ENDDO

 9999 RETURN
      END
