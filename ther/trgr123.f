      SUBROUTINE TRGR123( NBELFI, NBPIEF, NDIMES, NBCAS,
     %                    COPIEF, GRATEM,
     %                    NCAS  , CMPGRA )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACER LES VECTEURS GRADIENTS DES EF D'UN TYPE DONNE
C -----    OBJET 1D ou 2D ou 3D
C
C ENTREES :
C ---------
C NBELFI : NOMBRE D'ELEMENTS FINIS DU TYPE A TRAITER
C NBPIEF : NOMBRE DE POINTS DE CALCUL DES GRADIENTS PAR ELEMENT
C NDIMES : ESPACE DE TRAVAIL 1 OU 2 OU 3
C NBCAS  : NOMBRE DE CAS TRAITES
C COPIEF : LES NDIMES COORDONNEES DES POINTS DE CALCUL DES GRADIENTS
C GRATEM : LES GRADIENTS DE TEMPERATURE AUX POINTS D'INTEGRATION
C XYZPOI : LES 3 COORDONNEES DES POINTS DE L'OBJET
C NUPGEL : LE NUMERO DES POINTS DE CHAQUE ELEMENT FINI DE CE TYPE
C NCAS   : NUMERO DU CAS A TRAITER
C CMPGRA : NOMBRE DE CM POUR L'UNITE DE GRADIENT
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     OCTOBRE 1994
C23456---------------------------------------------------------------012
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      REAL              COPIEF(1:NBELFI,1:NBPIEF,1:NDIMES)
      DOUBLE PRECISION  GRATEM(1:NBELFI,1:NDIMES,1:NBPIEF,1:NBCAS)
      REAL              XYZ(3)
      EQUIVALENCE      (XYZ(1),X), (XYZ(2),Y), (XYZ(3),Z)
      REAL              CXYZCM(3)
      EQUIVALENCE      (CXYZCM(1),CONXCM), (CXYZCM(2),CONYCM),
     %                 (CXYZCM(3),CONZCM)
C
C     TRACE DES GRADIENTS DE TEMPERATURE
C     ----------------------------------
      IF( NDIMES .EQ. 1 ) THEN
C
C        DIMENSION 1
         DO 15 K=1,NBELFI
            DO 10 L=1,NBPIEF
C              LONGUEUR CM DU CORPS DE LA FLECHE
               CONXCM = REAL( GRATEM(K,1,L,NCAS) * CMPGRA )
               CONCM  = ABS( CONXCM )
C              LE TRACE DE LA FLECHE GRADIENT
               CALL T2FLEC( NCOUFL, COPIEF(K,L,1), 0.0,
     %                      CONCM,  CONXCM, 0.0 )
 10         CONTINUE
 15      CONTINUE
C
         ELSE IF( NDIMES .EQ. 2 ) THEN
C
C        DIMENSION 2
         DO 25 K=1,NBELFI
            DO 20 L=1,NBPIEF
C              LONGUEUR CM DU CORPS DE LA FLECHE
               CONXCM = REAL( GRATEM(K,1,L,NCAS) * CMPGRA )
               CONYCM = REAL( GRATEM(K,2,L,NCAS) * CMPGRA )
               CONCM  = SQRT( CONXCM * CONXCM + CONYCM * CONYCM )
C              LE TRACE DE LA FLECHE GRADIENT
               CALL T2FLEC( NCOUFL, COPIEF(K,L,1), COPIEF(K,L,2),
     %                      CONCM,  CONXCM, CONYCM )
 20         CONTINUE
 25      CONTINUE
C
      ELSE IF( NDIMES .EQ. 3 ) THEN
C
C        DIMENSION 3
         DO 35 K=1,NBELFI
            DO 30 L=1,NBPIEF
C              LONGUEUR CM DU CORPS DE LA FLECHE
               CONXCM = REAL( GRATEM(K,1,L,NCAS) * CMPGRA )
               CONYCM = REAL( GRATEM(K,2,L,NCAS) * CMPGRA )
               CONZCM = REAL( GRATEM(K,3,L,NCAS) * CMPGRA )
               CONCM  = SQRT( CONXCM**2 + CONYCM**2 + CONZCM**2  )
C              COORDONNEES OBJET DU POINT DEPART DE LA FLECHE
               XYZ(1) = COPIEF(K,L,1)
               XYZ(2) = COPIEF(K,L,2)
               XYZ(3) = COPIEF(K,L,3)
C              LE TRACE DE LA FLECHE GRADIENT
               CALL T3FLEC( NCOUFL, XYZ, CONCM, CXYZCM )
 30         CONTINUE
 35      CONTINUE
      ENDIF
C
      RETURN
      END
