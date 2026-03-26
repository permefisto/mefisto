      SUBROUTINE T21SOM( NOMOBJ , NUOBJT , MNSOMM )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACER LES SOMMETS DU TABLEAU 'XYZSOMMET' DEFINI PAR
C -----    LE TABLEAU MS D'ADRESSE MNSOMM
C
C ENTREES:
C --------
C NOMOBJ : NOM DE L'OBJET A TRACER
C NUOBJT : NUMERO DE L'OBJET
C MNSOMM : ADRESSE MCN DU TABLEAU 'XYZSOMMET' A TRACER
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC PARIS  OCTOBRE 1988
C ....................................................................
      IMPLICIT INTEGER (W)
      include"./incl/a___xyzsommet.inc"
      include"./incl/trvari.inc"
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE       (MCN(1),RMCN(1))
      CHARACTER*(*)     NOMOBJ
C
C     LE NOMBRE DE SOMMETS A TRACER
      NBSOM = MCN( MNSOMM + WNBSOM )
C
C     L'ADRESSE AVANT LE DEBUT DU TABLEAU XYZSOM
      MN = MNSOMM + WYZSOM
C
C     LA BOUCLE SUR LES SOMMETS
C     =========================
      DO 10 I=1,NBSOM
C
C        TRACE DU POINT
         CALL ITEMP2( RMCN(MN) , NOMOBJ , NUOBJT )
C
C        L'ADRESSE DES COORDONNEES DU SOMMET SUIVANT
         MN = MN + 3
 10   CONTINUE
      END
