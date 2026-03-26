      SUBROUTINE INTPP( NBPOID, POIDS, NBPOL, POLYP, PP )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  CALCUL DE L'INTEGRALE P P dx dy dz sur l'EF de REFERENCE
C -----  DES POLYNOMES LAGRANGE x DES POLYNOMES LAGRANGE
C      ( NECESSAIRE POUR L'EF DE TAYLOR HOOD )
C
C ENTREES:
C --------
C NBPOID : NOMBRE DE  POIDS DE LA FORMULE D'INTEGRATION NUMERIQUE
C POIDS  : VALEUR DES POIDS DE LA FORMULE D'INTEGRATION NUMERIQUE
C NBPOL  : NOMBRE DE POLYNOMES
C POLYP  : POLYP(NBPOLY,NBPOID)  ou  POLYP(I,L) = PI (XL,YL,ZL)
C
C SORTIE :
C --------
C PP     : INTEGRALE Pi Pj dX SUR L'EF UNITE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY      Mai 2010
C23456---------------------------------------------------------------012
      DOUBLE PRECISION POIDS(NBPOID), S
      DOUBLE PRECISION POLYP(NBPOL,NBPOID)
      DOUBLE PRECISION PP(NBPOL,NBPOL)
C
      DO J=1,NBPOL
         DO I=1,NBPOL
            S = 0D0
            DO M=1,NBPOID
               S = S + POIDS(M) * POLYP(I,M) * POLYP(J,M)
            ENDDO
            IF( ABS(S) .LT. 1D-14 ) S=0D0
            PP(I,J) = S
         ENDDO
      ENDDO
C
C     AFFICHAGE DES INTEGRALES
      PRINT 10000
10000 FORMAT(//'INTEGRALE P2 P2 dX pour l''EF de TAYLOR-HOOD')
      DO I=1,NBPOL
         PRINT 10001, (I,J,PP(I,J),J=1,NBPOL)
      ENDDO
10001 FORMAT(5(i3,i3,'=',D25.17))
C
      RETURN
      END
