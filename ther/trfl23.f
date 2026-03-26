      SUBROUTINE TRFL23( NDIMES, NBELFI, NBPNFX, NBCAS,
     %                   COPNFX, FLUXNP, NCAS  , CMFLUX )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACER LES VECTEURS FLUX NORMAUX D'UN ENSEMBLE D'EF
C -----           OBJET 1D OU 2D OU 3D
C
C ENTREES :
C ---------
C NDIMES : DIMENSION DE L"ESPACE DE L'OBJET (1 2 OU 3 )
C NBELFI : NOMBRE D'ELEMENTS FINIS DU TYPE A TRAITER
C NBPNFX : NOMBRE DE POINTS DE CALCUL DES FLUX NORMAUX PAR EF
C NBCAS  : NOMBRE DE CAS TRAITES
C COPNFX : (*,*,*,1) 3 COORDONNEES DES POINTS DE CALCUL DES FLUX NORMAUX
C          (*,*,*,2) 3 COMPOSANTES DU VECTEUR NORMAL UNITAIRE AUX POINTS
C FLUXNP : LES FLUX NORMAUX AUX POINTS D'INTEGRATION DES FACES DES EF
C NCAS   : NUMERO DU CAS A TRAITER
C CMFLUX : NOMBRE DE CM POUR L'UNITE DE GRADIENT
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     OCTOBRE 1994
C23456---------------------------------------------------------------012
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      include"./incl/ponoel.inc"
      REAL              COPNFX( 1:3, 1:NBPNFX, 1:NBELFI, 1:2 )
      DOUBLE PRECISION  FLUXNP( 1:NBPNFX, 1:NBELFI, 1:NBCAS )
      REAL              CXYZCM(3)
      EQUIVALENCE      (CXYZCM(1),CONXCM), (CXYZCM(2),CONYCM),
     %                 (CXYZCM(3),CONZCM)
C
C     TRACE DES FLUX NORMAUX AUX POINTS D'INTEGRATION DES SURFACES DES EF
C     -------------------------------------------------------------------
      IF( NDIMES .EQ. 1 ) THEN
C
C        DIMENSION 1
         DO 15 K=1,NBELFI
            DO 10 L=1,NBPNFX
C              LA DIRECTION DU VECTEUR
C              L'INTENSITE DU FLUX NORMAL EN CE POINT RAMENE EN CM
               FL = REAL( FLUXNP( L, K, NCAS ) * CMFLUX )
C              LONGUEUR CM DE LA FLECHE DANS LA DIRECTION NORMALE
               CONXCM = COPNFX( 1, L, K, 2 ) * FL
               CONCM = ABS( CONXCM )
C              LE TRACE DE LA FLECHE DU FLUX NORMAL
               CALL T2FLEC( NCOUFL, COPNFX(1,L,K,1), 0.0,
     %                      CONCM, CONXCM, 0.0 )
 10         CONTINUE
 15      CONTINUE
C
      ELSE IF( NDIMES .EQ. 2 ) THEN
C
C        DIMENSION 2
         DO 25 K=1,NBELFI
            DO 20 L=1,NBPNFX
C              LA DIRECTION DU VECTEUR
C              L'INTENSITE DU FLUX NORMAL EN CE POINT RAMENE EN CM
               FL = REAL( FLUXNP( L, K, NCAS ) * CMFLUX )
C              LONGUEUR CM DE LA FLECHE DANS LA DIRECTION NORMALE
               CONXCM = COPNFX( 1, L, K, 2 ) * FL
               CONYCM = COPNFX( 2, L, K, 2 ) * FL
               CONCM  = SQRT( CONXCM * CONXCM + CONYCM * CONYCM )
C              LE TRACE DE LA FLECHE DU FLUX NORMAL
               CALL T2FLEC( NCOUFL, COPNFX(1,L,K,1), COPNFX(2,L,K,1),
     %                      CONCM, CONXCM, CONYCM )
 20         CONTINUE
 25      CONTINUE
C
      ELSE IF( NDIMES .EQ. 3 ) THEN
C
C        DIMENSION 3
         DO 35 K=1,NBELFI
            DO 30 L=1,NBPNFX
C              LA DIRECTION DU VECTEUR
C              L'INTENSITE DU FLUX NORMAL EN CE POINT RAMENE EN CM
               FL = REAL( FLUXNP( L, K, NCAS ) * CMFLUX )
C              LONGUEUR CM DE LA FLECHE DANS LA DIRECTION NORMALE
               CONXCM = COPNFX( 1, L, K, 2 ) * FL
               CONYCM = COPNFX( 2, L, K, 2 ) * FL
               CONZCM = COPNFX( 3, L, K, 2 ) * FL
               CONCM  = SQRT( CONXCM**2 + CONYCM**2 + CONZCM**2 )
C              LE TRACE DE LA FLECHE DU FLUX NORMAL
               CALL T3FLEC( NCOUFL, COPNFX(1,L,K,1), CONCM, CXYZCM )
 30         CONTINUE
 35      CONTINUE
      ENDIF
C
      RETURN
      END
