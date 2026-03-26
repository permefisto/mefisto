      SUBROUTINE VSLPOB( NBOBPR, MNOBPR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : TRIER LES PLSV PREMIERS D'UN OBJET SELON
C ----- LES VOLUMES, PUIS LES SURFACES, PUIS LES LIGNES,
C       PUIS LES POINTS
C
C ENTREES :
C ---------
C NBOBPR : NOMBRE DE PLSV PREMIERS
C
C MODIFIES:
C ---------
C MNOBPR : ADRESSE MCN DU TABLEAU OBPR
C          OBPR(1,I) = NUMERO DU TYPE DU PLSV PREMIER
C                      1:POINT, 2:LIGNE, 3:SURFACE, 4:VOLUME
C          OBPR(2,I) = NUMERO DU PLSV PREMIER DANS SON LEXIQUE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE UPMC PARIS    OCTOBRE 1999
C2345X7..............................................................012
      include"./incl/pp.inc"
      COMMON  MCN(MOTMCN)
C
C     TRI DES OBJETS PREMIERS SELON LE TYPE DECROISSANT V, S, L, P
      MN0 = MNOBPR
      DO 50 I=1,NBOBPR-1
C        LE PLSV I A T IL UN PLSV DE NUMERO DE TYPE SUPERIEUR?
 10      NUTY0 = MCN(MN0  )
         NUPL0 = MCN(MN0+1)
         MN1   = MN0 + 2
         DO 30 J=I+1,NBOBPR
C           LE TYPE DU PLSV J
            NUTY1 = MCN(MN1  )
            NUPL1 = MCN(MN1+1)
            IF( NUTY0 .LT. NUTY1 ) THEN
C              PERMUTATION CIRCULAIRE ENTRE CES 2 TYPES
               DO 20 MN=MN1-2,MN0,-2
C                 ECHANGE DES 2 PLSV
                  MCN(MN+2) = MCN(MN)
                  MCN(MN+3) = MCN(MN+1)
 20            CONTINUE
               MCN(MN0  ) = NUTY1
               MCN(MN0+1) = NUPL1
               GOTO 10
            ENDIF
            MN1 = MN1 + 2
 30      CONTINUE
C        PASSAGE AU PLSV SUIVANT
         MN0 = MN0 + 2
 50   CONTINUE
C
      RETURN
      END
