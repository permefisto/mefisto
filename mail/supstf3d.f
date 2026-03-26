      SUBROUTINE SUPSTF3D( NOST, NBSOEF, NBEFOB, NOSOEF )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    SUPPRIMER LES EF DE SOMMET NOST
C -----

C ENTREES:
C --------
C NOST   : NUMERO XYZSOM DU SOMMET DE BARYCENTRE DES VOISINS A CALCULER
C NBSOEF : NOMBRE DE SOMMETS DECLARABLES DANS NOSOEF PAR EF
C NBEFOB : NOMBRE D'EF ACTIFS DANS LE TABLEAU NOSOEF

C MODIFIE:
C --------
C NOSOEF : NUMERO XYZSOM DES NBSOEF SOMMETS DES NBEFOB EF
C          LES EF SUPPRIMES ONT LEURS NO XYZSOM DES SOMMETS NULS
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: PERRONNET ALAIN   Saint PIERRE du PERRAY              Mai 2020
C2345X...............................................................012
      INTEGER     NOSOEF(NBSOEF,NBEFOB)

      NBF = 0
      IF( NOST .LE. 0 ) GOTO 9999

      DO 100 NEF = 1, NBEFOB

         IF( NOSOEF(1,NEF) .LE. 0 ) GOTO 100

C        RECHERCHE DU SOMMET NOST PARMI LES NBSOEF SOMMETS DE NEF
         DO I = 1, NBSOEF
            IF( NOSOEF( I, NEF ) .EQ. NOST ) GOTO 10
         ENDDO
         GOTO 100

C        LES NBSOEF SOMMETS SONT MIS A ZERO
 10      NBF = NBF + 1
         DO I = 1, NBSOEF
            NOSOEF( I, NEF ) = 0
         ENDDO

 100  ENDDO

 9999 RETURN
      END

   
