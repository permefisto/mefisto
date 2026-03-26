      SUBROUTINE LEGCOULSO( VMIN, VMAX )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACE LA LEGENDE DES COULEURS SELON LA SOLUTION
C -----    COULEURS => VALEUR DE LA SOLUTION

C ENTREES:
C --------
C VMIN   : VALEUR MINIMALE de la SOLUTION
C VMAX   : VALEUR MAXIMALE de la SOLUTION
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET TEXAS A & M UNIVERSITY at QATAR  Fevrier 2012
C23456---------------------------------------------------------------012
      include"./incl/xvfontes.inc"
      include"./incl/trvari.inc"
      REAL           VMIN, VMAX
      CHARACTER*10   KNOM

C     SAUVEGARDE DU NUMERO DE LA FONTE DE CARACTERES ACTUELS
      NOFONT0 = NOFONT

C     CHOIX DE LA FONTE
      LHPXCA = 20
      CALL CHOIXFONTE( LHPXCA )

C     LE TRACE DE LA LEGENDE: COULEURS => VALEURS
      NBCOUL = NDCOUL - NDCORE
      IF( NOPACL .EQ. 17 ) THEN
C        1/2 PALETTE ARC EN CIEL + 1/2 PALETTE GRISE
         NBCOUL = NBCOUL / 2
      ENDIF
      NDCOU = NDCORE + NBCOUL

      NCPAS  = NBCOUL / 10
      VALPAS = (VMAX-VMIN) / 10
      VAL    =  VMIN

C     TRACE DE 11 VALEURS A DROITE EN BAS DE LA FENETRE DE TRACE
      NCOUL = N1COUL

C     COORDONNEES DANS LA FENETRE
      NX = LAPXFE - 150
      NY = LHPXFE - 30

      DO I = 0, 10
         CALL XVCOULEUR( NCOUL )
         CALL XVRECTANGLE( NX, NY, 35, 10 )
         WRITE( KNOM(1:10), '(G10.3)' ) VAL
         CALL XVTEXTE( KNOM(1:10), 10, NX+40, NY+10 )
         NCOUL = NCOUL + NCPAS
         IF( NCOUL .GT. NDCOU ) NCOUL = NDCOU
         VAL = VAL + VALPAS
         NY  = NY - 15
      ENDDO

C     RESTAURATION DU NUMERO DE LA FONTE DE CARACTERES ACTUELS
      CALL CHARGEFONTE( NOFONT0 )

      RETURN
      END
