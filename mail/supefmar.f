      SUBROUTINE SUPEFMAR( NBEFMA, NOEFMA, NBSOEF, NBEFOB, NOSOEF )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    SUPPRIMER LES NBEFMA EF MARQUES
C -----

C ENTREES:
C --------
C NBEFMA : NOMBRE D'EF MARQUES
C NOEFMA : NUMEROS NOSOEF DES NBEFMA EF MARQUES
C NBSOEF : NOMBRE DE SOMMETS DECLARABLES DANS NOSOEF PAR EF
C NBEFOB : NOMBRE D'EF ACTIFS DANS LE TABLEAU NOSOEF

C MODIFIE:
C --------
C NOSOEF : NUMERO XYZSOM DES NBSOEF SOMMETS DES NBEFOB EF
C          LES EF SUPPRIMES ONT LEURS NO XYZSOM DES SOMMETS NULS
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: PERRONNET ALAIN   Saint PIERRE du PERRAY              Mai 2020
C2345X...............................................................012
      INTEGER     NOEFMA(NBEFMA), NOSOEF(NBSOEF,NBEFOB)

      DO 10 N = 1, NBEFMA

         NEF = NOEFMA( N )
         IF( NEF .LE. 0 ) GOTO 10

C        LES NBSOEF SOMMETS SONT MIS A ZERO
         DO I = 1, NBSOEF
            NOSOEF( I, NEF ) = 0
         ENDDO

 10   ENDDO

      RETURN
      END

   
