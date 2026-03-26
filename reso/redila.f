      SUBROUTINE REDILA( NYOBJT, NUOBJT, XPI, YPI, ZPI, MNDILA,
     %                   DILATA, TEMPIN )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++012
C BUT : REMPLIR LA VALEUR Alpha DU COEFFICIENT DE DILATATION THERMIQUE
C -----            (L - L0)/L0 = Alpha (Temp - Temp0)  en 1/Kelvin
C               LA VALEUR DE LA TEMPERATURE INITIALE DE REFERENCE
C               AU POINT XPI,YPI,ZPI
C
C      (en Anglais: DILATATION=THERMAL EXPANSION COEFFICIENT)
C       Ce COEFFICIENT PEUT AUSSI ETRE PRIS COMME COEFFICIENT dans
C       l'APPROXIMATION DE BOUSSINESQ SUR LA DENSITE DE MASSE
C       D'UN FLUIDE CHAUD (Rho0-Rho)/Rho = CoBOUS (Temp-Temp0)

C ENTREES :
C ---------
C NYOBJT : NUMERO DU TYPE (3:SURFACE, 4:VOLUME)
C NUOBJT : NUMERO DE SURFACE OU VOLUME
C XPI,YPI,ZPI : LES 3 COORDONNEES DU POINT D'INTEGRATION
C MNDILA : ADRESSE MCN DU TABLEAU 'DILATATION'
C
C SORTIE :
C --------
C DILATA : LE COEFFICIENT DE DILATATION THERMIQUE HOMOGENE ISOTROPE
C TEMPIN : TEMPERATURE INITIALE
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++012
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     JUILLET 1989
C....6...............................................................012
      IMPLICIT INTEGER(W)
      include"./incl/a___dilatation.inc"
      include"./incl/ctemps.inc"
      include"./incl/cthet.inc"
      include"./incl/cnonlin.inc"
C
      include"./incl/pp.inc"
      COMMON           MCN(MOTMCN)
      REAL             RMCN(1)
      EQUIVALENCE     (MCN(1),RMCN(1))
C
      DOUBLE PRECISION XPI,YPI,ZPI,DILATA,TEMPIN,XYZ(7)
C
C     LE TYPE DES DONNEES DE LA DILATATION
      LTDILA = MCN( MNDILA + WTDILA )
C
C     REMPLISSAGE SELON LE TYPE
      IF( LTDILA .EQ. 1 ) THEN
C
C        MATERIAU HOMOGENE ISOTROPE
C        ==========================
         DILATA = RMCN( MNDILA + WILATA )
         TEMPIN = RMCN( MNDILA + WEMPIN )
C
      ELSE
C
C        FONCTION UTILISATEUR
C        ====================
         XYZ(1) = TEMPS
         XYZ(2) = XPI
         XYZ(3) = YPI
         XYZ(4) = ZPI
         XYZ(5) = NYOBJT
         XYZ(6) = NUOBJT
         XYZ(7) = TEMPEL
C        LE COEFFICIENT DE DILATATION THERMIQUE
         CALL FONVAL( MCN(MNDILA+WFDILA), 7, XYZ,  NCODEV, DILATA )
C        LA TEMPERATURE DE REFERENCE POUR LA DILATATION THERMIQUE
         CALL FONVAL( MCN(MNDILA+WFTEMI), 6, XYZ,  NCODEV, TEMPIN )
C
      ENDIF
      END
