      SUBROUTINE INTRTE( STR, STE, LINTER )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    L'INTERSECTION D'UN TRIANGLE ET D'UN TETRAEDRE EST ELLE VIDE?
C -----  I.E. AUCUNE ARETE DU TRIANGLE INTERSECTE LES FACES DU TETRAEDRE
C          +  AUCUNE ARETE DU TETRAEDRE INTERSECTE LE TRIANGLE
C          +  AUCUN DES 3 SOMMETS DU TRIANGLE EST INTERNE AU TETRAEDRE?
C ENTREES:
C --------
C STR    : 3 COORDONNEES DES 3 SOMMETS DU TRIANGLE
C STE    : 3 COORDONNEES DES 4 SOMMETS DU TETRAEDRE
C
C SORTIE :
C --------
C LINTER : 0 PAS D'INTERSECTION
C          1 SI UNE ARETE DU TRIANGLE INTERSECTE UNE FACE  DU TETRAEDRE
C          2 SI LE TRIANGLE EST INTERSECTE PAR   UNE ARETE DU TETRAEDRE
C          3 SI UN DES SOMMETS DU TRIANGLE EST STRICTEMENT INTERNE
C            AU TETRAEDRE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY  DECEMBRE 2014
C MODIFS: ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY  AOUT     2015
C23456...............................................................012
      DOUBLE PRECISION  STR(3,3), STE(3,4), XYZINT(3), CBT(4), VOLUMT

      INTEGER  NOSOFATE(3,4), NOSOARTE(2,6)
      DATA     NOSOFATE / 1,3,2,  2,3,4, 3,1,4, 4,1,2 /
      DATA     NOSOARTE / 1,2, 2,3, 3,1, 1,4, 2,4, 3,4 /

      NSA1 = 3
      DO NSA2=1,3
         DO L=1,4
C           INTERSECTION DE L'ARETE NSA2 DU TRIANGLE AVEC LA FACE L DU TETRAEDRE
            CALL INARTR( STR(1,NSA1), STR(1,NSA2),
     %                   STE(1,NOSOFATE(1,L)),
     %                   STE(1,NOSOFATE(3,L)),
     %                   STE(1,NOSOFATE(2,L)),
     %                   LINTER, XYZINT, CBT )
            IF( LINTER .EQ. 1   .AND.
     %          1D-6 .LT. CBT(1) .AND. CBT(1) .LT. 0.999999D0  .AND.
     %          1D-6 .LT. CBT(2) .AND. CBT(2) .LT. 0.999999D0  .AND.
     %          1D-6 .LT. CBT(3) .AND. CBT(3) .LT. 0.999999D0  ) THEN
C               OUI: UN POINT D'INTERSECTION
                GOTO 9000
            ENDIF
         ENDDO
         NSA1 = NSA2
      ENDDO

      DO L=1,6
C        INTERSECTION DE L'ARETE L DU TETRAEDRE AVEC LE TRIANGLE
         CALL INARTR( STE(1,NOSOARTE(1,L)),
     %                STE(1,NOSOARTE(2,L)),
     %                STR(1,1), STR(1,2), STR(1,3),
     %                LINTER, XYZINT, CBT )
         IF( LINTER .EQ. 1   .AND.
     %       1D-6 .LT. CBT(1) .AND. CBT(1) .LT. 0.999999D0 .AND.
     %       1D-6 .LT. CBT(2) .AND. CBT(2) .LT. 0.999999D0 .AND.
     %       1D-6 .LT. CBT(3) .AND. CBT(3) .LT. 0.999999D0 ) THEN
C            OUI: UN POINT D'INTERSECTION
             LINTER = 2
             GOTO 9000
         ENDIF
      ENDDO

      DO K=1,3
C        COORDONNEES BARYCENTRIQUES DU SOMMET K DU TRIANGLE DANS LE TETRAEDRE
         CALL COBATET( STR(1,K), STE(1,1), STE(1,2), STE(1,3), STE(1,4),
     %                 VOLUMT, CBT, IERR )
         IF( IERR .EQ. 0     .AND.
     %       1D-6 .LT. CBT(1) .AND. CBT(1) .LT. 0.999999D0 .AND.
     %       1D-6 .LT. CBT(2) .AND. CBT(2) .LT. 0.999999D0 .AND.
     %       1D-6 .LT. CBT(3) .AND. CBT(3) .LT. 0.999999D0 .AND.
     %       1D-6 .LT. CBT(4) .AND. CBT(4) .LT. 0.999999D0 ) THEN
C            OUI: LE SOMMET K DU TRIANGLE EST INTERNE AU TETRAEDRE
             LINTER = 3
             XYZINT(1) = STR(1,K)
             XYZINT(2) = STR(2,K)
             XYZINT(3) = STR(3,K)
             GOTO 9000
         ENDIF
      ENDDO

C     PAS D'INTERSECTION
      LINTER = 0

cccC     TRACE EVENTUEL DE L'INTERSECTION
ccc 9000 CALL TRINTRTE( STR, STE, LINTER, XYZINT )

 9000 RETURN
      END
