      SUBROUTINE SANSBL( KNOM, NBCAR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : SUPPRIME LES BLANCS D'UNE CHAINE DE CARACTERES EN CADRANT A GAUCHE
C -----

C MODIFIE:
C --------
C KNOM   : CHAINE AVANT ET APRES

C SORTIE :
C --------
C NBCAR  : NOMBRE DE CARACTERES NON BLANCS CADRES A GAUCHE DE KNOM
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS        JUIN 1994
C2345X7..............................................................012
      CHARACTER*(*) KNOM

C     NOMBRE DE CARACTERES DECLARES DE KNOM
      NBCKNOM = LEN( KNOM )

      NBCAR = 0
      DO I = 1,NBCKNOM
         IF( KNOM(I:I) .NE. ' ' ) THEN
            NBCAR = NBCAR + 1
            KNOM(NBCAR:NBCAR) = KNOM(I:I)
         ENDIF
      ENDDO

C     COMPLETION PAR DES BLANCS AU DELA DE NBCAR
      DO I = NBCAR+1, NBCKNOM
         KNOM( I : I ) = ' '
      ENDDO

      RETURN
      END
