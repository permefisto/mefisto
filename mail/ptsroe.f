      SUBROUTINE PTSROE( P, PTXYZD, NONOUI )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  LE POINT P EST IL SUR L'UNE DES FACES DE L'OCTAEDRE ENGLOBANT?
C -----  LES POINTS 1 ET 2 DE PTXYZD SONT EXTREMITES D'UNE ARETE DE
C        L'OCTAEDRE ENGLOBANT
C
C ENTREES:
C --------
C P      : LES 3 COORDONNEES DU POINT
C PTXYZD : X Y Z DISTANCE SOUHAITEE DES POINTS
C
C SORTIES:
C --------
C NONOUI : 1 SI LE POINT EST DANS OU SUR L'UNE DES 8 FACES DE L'OCTAEDRE
C            ENGLOBANT DE SOMMETS 1 A 6 DANS PTXYZD
C          0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC       MARS 1993
C....................................................................012
      DOUBLE PRECISION  PTXYZD(4,*), P(3), COBARY(3), DISPLA, ARETE
      INTEGER           NOSOTR(3)
C
C     LA LONGUEUR DE L'ARETE DE L'OCTAEDRE ENGLOBANT
      ARETE  = SQRT( (PTXYZD(1,2)-PTXYZD(1,1))**2 +
     %               (PTXYZD(2,2)-PTXYZD(2,1))**2 +
     %               (PTXYZD(3,2)-PTXYZD(3,1))**2 )
      ARETE  = ARETE * 1D-8
      NONOUI = 0
C
C     LE POINT P EST IL SUR L'UNE DES FACES SUPERIEURES
      NOSOTR(1) = 1
      DO I = 2, 5
         NOSOTR(2) = I
         IF( I .EQ. 5 ) THEN
            NOSOTR(3) = 2
         ELSE
            NOSOTR(3) = I+1
         ENDIF
C        DISTANCE DU POINT P AU PLAN
         CALL DIPTPL( P,                   PTXYZD(1,NOSOTR(1)),
     %                PTXYZD(1,NOSOTR(2)), PTXYZD(1,NOSOTR(3)), DISPLA )
         IF( DISPLA .LT. ARETE ) THEN
C           LE POINT EST IL DANS OU SUR LES ARETES DU TRIANGLE
            CALL PTDSTR( P, PTXYZD, NOSOTR, COBARY, NONOUI )
            IF( NONOUI .GT. 0 ) RETURN
         ENDIF
      ENDDO
C
C     LE POINT P EST IL SUR L'UNE DES FACES INFERIEURES
      NOSOTR(1) = 6
      DO I = 2, 5
         NOSOTR(2) = I
         IF( I .EQ. 5 ) THEN
            NOSOTR(3) = 2
         ELSE
            NOSOTR(3) = I+1
         ENDIF
C        DISTANCE DU POINT P AU PLAN
         CALL DIPTPL( P,                   PTXYZD(1,NOSOTR(1)),
     %                PTXYZD(1,NOSOTR(2)), PTXYZD(1,NOSOTR(3)), DISPLA )
         IF( DISPLA .LT. ARETE ) THEN
C           LE POINT EST IL DANS OU SUR LES ARETES DU TRIANGLE
            CALL PTDSTR( P, PTXYZD, NOSOTR, COBARY, NONOUI )
            IF( NONOUI .GT. 0 ) RETURN
         ENDIF
      ENDDO

      RETURN
      END
