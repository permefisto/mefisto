      SUBROUTINE REDOGP( NDIM,   NBPOLY, X,      NBJEUX, JEU,
     %                   NOOBSF, NUMISU, NUMASU, LTDESU,
     %                   Rho,    Omega,  Alfa,   Beta, Force )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : RETROUVER LES COEFFICIENTS SUPPOSES CONSTANTS DU PROBLEME DE
C ----- GROSS-PITAEVSKII POUR UNE TRIANGULATION P1 (2P1D)

C ENTREES:
C --------
C NDIM   : DIMENSION DE L'ESPACE (2 ou 3)
C NBPOLY : NOMBRE DE SOMMETS DES EF (3 ou 4)
C X      : X(NBPOLY,NDIM) COORDONNEES DES NBPOLY SOMMETS DE L'EF
C NBJEUX : NOMBRE DE JEUX DE DONNEES
C JEU    : NUMERO DU JEU  DE DONNEES POUR CE CALCUL DE LA MATRICE ELEMENTAIRE

C NOOBSF : NUMERO DE L'OBJET SURFACE DE CET ELEMENT
C NUMISU : NUMERO MINIMAL DES OBJETS SURFACES
C NUMASU : NUMERO MAXIMAL DES OBJETS SURFACES
C LTDESU : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES DE THERMIQUE
C          DES OBJETS SURFACES
C
C SORTIES:
C --------
C Omega  : VITESSE ANGULAIRE DE LA ROTATION

C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & Saint Pierre du Perray   Mai 2014
C23456---------------------------------------------------------------012
      include"./incl/donthe.inc"
      include"./incl/cthet.inc"
      include"./incl/ctemps.inc"
      include"./incl/cnonlin.inc"

      REAL              X(NBPOLY, NDIM)
      INTEGER           LTDESU( 1:MXDOTH, 1:NBJEUX, NUMISU:NUMASU )
      DOUBLE PRECISION  Rho, Omega(3), Alfa, Beta, Force(2)
      DOUBLE PRECISION  COND(6), XYZ(3)

      TEMPEL = 0D0
      ONDEPI = 0D0
      XYZ(1) = X(1,1)
      XYZ(2) = X(1,2)
      XYZ(3) = 0D0

C     INITIALISATION DE LA VITESSE ANGULAIRE Omega/Z SUPPOSEE CONSTANTE
C     -----------------------------------------------------------------
      MN = LTDESU(LPVIANT,JEU,NOOBSF)
      IF( MN .GT. 0 ) THEN
C        RECUPERATION DE LA VITESSE ANGULAIRE Omega SUPPOSEE CONSTANTE
C        VECTEUR(3) DE VITESSE ANGULAIRE AU 1-ER SOMMET DU TRIANGLE
         XYZ(1) = X(1,1)
         XYZ(2) = X(1,2)
         CALL REVIAN( 3, NOOBSF, XYZ(1), XYZ(2), 0D0,
     %                   LTDESU(LPVIANT,JEU,NOOBSF), Omega )
C        Si SEUL OMEGA(1) EST DONNE C'EST SUPPOSE ETRE Omega(3)/axeZ
      ELSE
C        PAS DE ROTATION
         Omega(1) = 0D0
         Omega(2) = 0D0
         Omega(3) = 0D0
      ENDIF

C     DENSITE DE MASSE SUPPOSEE CONSTANTE
C     -----------------------------------
      MN = LTDESU(LPMAST,JEU,NOOBSF)
      IF( MN .GT. 0 ) THEN
         CALL REMASS( 3, NOOBSF, 3, XYZ, MN, Rho )
      ELSE
         Rho = 1D0
      ENDIF

C     COEFFICIENT DU LAPLACIEN = CONDUCTIVITE(1) SUPPOSEE CONSTANTE
C     -------------------------------------------------------------
      MN = LTDESU(LPCOND,JEU,NOOBSF)
      IF( MN .GT. 0 ) THEN
         CALL RECOND( 3, NOOBSF, 3, XYZ,
     %                LTDESU(LPCOND,JEU,NOOBSF), COND )
         Alfa = COND(1)
      ELSE
         Alfa = 0.5D0
      ENDIF

C     Beta COEFFICIENT DU TERME V**2+W**2 SUPPOSEE CONSTANT
C     -----------------------------------------------------
      MN = LTDESU(LPCOET,JEU,NOOBSF)
      IF( MN .GT. 0 ) THEN
         CALL RENLSE( 3, NOOBSF, 3, XYZ, TEMPS, 0D0, 0D0, MN,
     %                Beta )
      ELSE
         Beta = 0D0
      ENDIF

C     FORCE R et I AU SECOND MEMBRE SUPPOSEE CONSTANTE
C     ------------------------------------------------
      MN = LTDESU(LPSOUR,JEU,NOOBSF)
      IF( MN .GT. 0 ) THEN
         CALL REFORC( 3, NOOBSF, 2, XYZ(1),XYZ(2),XYZ(3),
     %                              0D0,0D0,0D0,
     %                MN, FORCE )
      ELSE
         FORCE(1) = 0D0
         FORCE(2) = 0D0
      ENDIF

      RETURN
      END
