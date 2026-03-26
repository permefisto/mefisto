      SUBROUTINE ECLAIR( NORMAL, RCOUL )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   LA COULEUR EST PROPORTIONNELLE AU COSINUS DE L'ANGLE
C -----   ENTRE LA NORMALE AU POINT ET LES DIRECTIONS DES ECLAIRAGES
C         DEFINIS PAR LES LAMPES ET LE POINT VU DU COMMON TRVARI
C
C ENTREE :
C --------
C NORMAL : 3 COORDONNEES DU VECTEUR NORMAL
C
C SORTIE :
C --------
C RCOUL  : VALEUR REELLE DE LA COULEUR ENTRE N1COUL ET NDCOUL
C          SERT DIRECTEMENT A APPELER LE SP TRIACOUL
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC  PARIS       AOUT 1998
C2345X7..............................................................012
      include"./incl/trvari.inc"
      REAL     NORMAL(3), RCOUL
C
      IF( NBLAMP .LE. 0 ) THEN
C
C        PAS DE LAMPE => L'OEIL SERT DE LAMPE
         RCOUL = ABS( DIREVI(1) * NORMAL(1)
     %              + DIREVI(2) * NORMAL(2)
     %              + DIREVI(3) * NORMAL(3) )
C
      ELSE
C
C        AU MOINS UNE LAMPE => COULEUR MOYENNE DES ECLAIRAGES
         RCOUL = 0
         DO 10 I=1,NBLAMP
            S = SQRT( (AXOLAM(1,I)-AXOPTV(1)) ** 2
     %              + (AXOLAM(2,I)-AXOPTV(2)) ** 2
     %              + (AXOLAM(3,I)-AXOPTV(3)) ** 2 )
            IF( S .EQ. 0 ) S=1.
CCC
CCC            RCOUL = RCOUL + ABS( (AXOLAM(1,I)-AXOPTV(1)) * NORMAL(1)
CCC     %                         + (AXOLAM(2,I)-AXOPTV(2)) * NORMAL(2)
CCC     %                         + (AXOLAM(3,I)-AXOPTV(3)) * NORMAL(3) )/S
CCC
            S = ABS(  (AXOLAM(1,I)-AXOPTV(1)) * NORMAL(1)
     %              + (AXOLAM(2,I)-AXOPTV(2)) * NORMAL(2)
     %              + (AXOLAM(3,I)-AXOPTV(3)) * NORMAL(3) ) / S
            RCOUL = MAX( RCOUL, S )
 10      CONTINUE
CCC         RCOUL = RCOUL / NBLAMP
C
      ENDIF
C
      RCOUL = N1COUL + (NDCOUL-N1COUL) * (1.-RCOUL)
      RETURN
      END
