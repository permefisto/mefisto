      SUBROUTINE BARYSTVO( NOST, NBSOEF, NBEFOB, NOSOEF, XYZSOM, XYZBAR)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LES XYZ DU BARYCENTRE DES BARYCENTRES DES EF
C -----    DE SOMMET NOST

C ENTREES:
C --------
C NOST   : NUMERO XYZSOM DU SOMMET DE BARYCENTRE DES VOISINS A CALCULER
C NBSOEF : NOMBRE DE SOMMETS DECLARABLES DANS NOSOEF PAR EF
C NBEFOB : NOMBRE D'EF ACTIFS DANS LE TABLEAU NOSOEF
C NOSOEF : NUMERO XYZSOM DES NBSOEF SOMMETS DES NBEFOB EF
C XYZSOM : XYZ DES SOMMETS DU MAILLAGE

C SORTIE :
C --------
C XYZBAR : XYZ DU BARYCENTRE DES BARYCENTRES DES EF DE SOMMET NOST
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: PERRONNET ALAIN   Saint PIERRE du PERRAY              Mai 2020
C2345X...............................................................012
      INTEGER     NOSOEF(NBSOEF,NBEFOB)
      REAL        XYZSOM(3,*), XYZBAR(3), XYZB(3)

      IF( NOST .LE. 0 ) GOTO 9999

      DO K = 1, 3
         XYZBAR( K ) = 0
      ENDDO

      NBF = 0
      DO 100 NEF = 1, NBEFOB

         IF( NOSOEF(1,NEF) .LE. 0 ) GOTO 100

C        RECHERCHE DU SOMMET NOST
         DO I = 1, NBSOEF
            IF( NOSOEF( I, NEF ) .EQ. NOST ) GOTO 10
         ENDDO
         GOTO 100

C        CALCUL DES XYZ DU BARYCENTRE DE l'EF NEF
 10      NBF = NBF + 1
         DO K = 1, 3
            XYZB( K ) = 0
         ENDDO

         NBS = 0
         DO I = 1, NBSOEF
            NS = NOSOEF( I, NEF )
            IF( NS .GT. 0 ) THEN
               NBS = NBS + 1
               DO K = 1, 3
                  XYZB( K ) = XYZB( K ) + XYZSOM( K, NS )
               ENDDO
            ENDIF
         ENDDO

         DO K = 1, 3
            XYZBAR( K ) = XYZBAR( K ) + XYZB( K ) / NBS
         ENDDO

 100  ENDDO

C     CALCUL DES XYZ DU BARYCENTRE DES BARYCENTRE DES EF DE SOMMET NOST
      DO K = 1, 3
         XYZBAR( K ) = XYZBAR( K ) / NBF
      ENDDO


 9999 RETURN
      END

   
