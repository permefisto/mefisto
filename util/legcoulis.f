      SUBROUTINE LEGCOULIS( NBISO, VALISO )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACE LA LEGENDE DES COULEURS SELON LES ISO-SOLUTIONS
C -----    COULEURS => VALEUR DE LA SOLUTION
C          AFFICHAGE DE LA LEGENDE D'AU MAXIMUM 20 ISOS
C
C ENTREES:
C --------
C NBISO  : NOMBRE D'ISOVALEURS DE LA SOLUTION
C VALISO : VALEURS DES NBISO de la SOLUTION
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET TAXAS A & M UNIVERSITY at QATAR  Fevrier 2012
C23456---------------------------------------------------------------012
      include"./incl/trvari.inc"
      CHARACTER*20   KNOM
      REAL           VALISO(NBISO)
C
      NBCOUL = NDCOUL - N1COUL + 1
      ISOPAS = 1 + (NBISO-1) / 20
C
C     COORDONNEES DANS LA FENETRE
      NX = LAPXFE - 220
      NY = LHPXFE - 30
C
      NBFOIS = ( NDCOUL - N1COUL + 1 ) / 10 - 2
C
      DO I=NBISO,1,-ISOPAS
C
C        LA COULEUR DE L'ISOTHERME ET DU RECTANGLE
         NCOUL = N1COUL - 1 + MOD(I,NBCOUL)
         IF( NDIMLI .EQ. 3 ) THEN
C           COULEUR VIVE POUR LA LEGENDE
            NCOUL = N1COUL - 1 + I + 10 * NBFOIS
         ENDIF
         CALL XVCOULEUR( NCOUL )
         CALL XVRECTANGLE( NX, NY, 30, 10 )
C
C        TRACE DE LA VALEUR
         WRITE( KNOM(1:4), '(I4)' ) NBISO+1-I
         KNOM(5:7) = ' : '
         WRITE( KNOM(8:17), '(G10.3)' ) VALISO(NBISO+1-I)
         CALL XVTEXTE( KNOM(1:17), 17, NX+40, NY+10 )
         NY = NY - 15
C
      ENDDO
C
      RETURN
      END
