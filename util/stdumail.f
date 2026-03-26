      SUBROUTINE STDUMAIL( NS, NBSOEF, NBEFOB, NOSOEF, NONOUI )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    SAVOIR SI LE SOMMET NS EST UN DES SOMMETS DES EF ACTIFS
C -----    DU MAILLAGE NOSOEF

C ENTREES:
C --------
C NS     : >0 NUMERO DU SOMMET
C NBSOEF : >0 NOMBRE DE SOMMETS PAR EF
C NBEFOB : NOMBRE D'EF DU MAILLAGE NOSOEF
C NOSOEF : NUMERO DES NBSOEF DES NBEFOB EF 

C SORTIE :
C --------
C NONOUI : =1 NS   EST     UN SOMMET DES EF ACTIFS DU MAILLAGE
C          =0 NS N'EST PAS UN SOMMET DES EF ACTIFS DU MAILLAGE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: PERRONNET ALAIN Saint PIERRE du PERRAY              Avril 2020
C2345X...............................................................012
      INTEGER   NOSOEF(NBSOEF,NBEFOB)

      DO N = 1, NBEFOB

         IF( NOSOEF(1,N) .GT. 0 ) THEN
C           EF ACTIF
            DO K = 1, NBSOEF
               IF( NS .EQ. NOSOEF(K,N) ) GOTO 1
            ENDDO
         ENDIF

      ENDDO

C     NS N'EST PAS UN SOMMET DES EF ACTIFS DU MAILLAGE
      NONOUI = 0
      GOTO 9999

C     NS   EST     UN SOMMET DES EF ACTIFS DU MAILLAGE
 1    NONOUI = 1

 9999 RETURN
      END
