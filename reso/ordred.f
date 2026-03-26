      SUBROUTINE ORDRED ( NC, TABL, MUTA )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    ORDONNER LES NC VALEURS DU TABLEAU TABL SUIVANT LEURS MODULES
C -----    CROISSANTS. LE TABLEAU MUTA CONTIENT LA PERMUTATION EFFECTUEE
C
C ENTREE :
C --------
C NC     : NOMBRE DE VALEURS A ORDONNER
C
C MODIFIE:
C --------
C TABL   : TABLEAU DES VALEURS AVANT ET APRES ORDONNANCEMENT
C
C SORTIE :
C --------
C MUTA   : ANCIENNNE POSITION DE LA VALEUR ACTUELLE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : BENAZETH.GOURDIN-SERVENIERE      LAN189 PARIS   DECEMBRE 1976
C          ALAIN PERRONNET                  LAN189 PARIS   JANVIER  1979
C ......................................................................
      DOUBLE PRECISION TABL(NC),S,T
      INTEGER          MUTA(NC)
C
C     INITIALISATION DE LA PERMUTATION
      DO 4 I=1,NC
         MUTA(I) = I
    4 CONTINUE
C
      DO 10 I=1,NC-1
C
C        RECHERCHE DU PLUS PETIT MODULE
         S  = TABL(I)
         IJ = I
         DO 2 J=I+1,NC
            T = TABL(J)
            IF( T .GE. S ) GOTO 2
            S  = T
            IJ = J
    2    CONTINUE
C
C        PERMUTATION DES VALEURS I ET IJ. MISE A JOUR DU TABLEAU MUTA
         IF( IJ .NE. I ) THEN
            S = TABL(IJ)
            TABL(IJ) = TABL(I)
            TABL(I)  = S
            J        = MUTA(IJ)
            MUTA(IJ) = MUTA(I)
            MUTA(I)  = J
         ENDIF
C
   10 CONTINUE
      END
