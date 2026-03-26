      SUBROUTINE T31SOM( NMOBJT , NUOBJT , MNSOMM )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACER LES SOMMETS DU TABLEAU 'XYZSOMMET' DEFINI PAR
C -----    LE TABLEAU MS D'ADRESSE MNSOMM
C
C ENTREES:
C --------
C NMOBJT : NOM DE L'OBJET A TRACER
C NUOBJT : NUMERO DE L'OBJET DANS SON LEXIQUE
C MNSOMM : ADRESSE MCN DU TABLEAU 'XYZSOMMET' A TRACER
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC PARIS  OCTOBRE 1988
C ......................................................................
      include"./incl/a___xyzsommet.inc"
      include"./incl/trvari.inc"
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      include"./incl/pp.inc"
      COMMON         MCN(MOTMCN)
      REAL          RMCN(1)
      EQUIVALENCE  (RMCN(1),MCN(1))

      CHARACTER*(*)     NMOBJT

C     LE NOMBRE DE SOMMETS A TRACER
      NBSOM = MCN( MNSOMM + WNBSOM )

C     L'ADRESSE AVANT LE DEBUT DU TABLEAU XYZSOM
      MN = MNSOMM + WYZSOM

C     LA BOUCLE SUR LES SOMMETS
C     =========================
      DO I=1,NBSOM

C        TRACE EVENTUEL DU POINT
         CALL ITEMP3( RMCN(MN), NMOBJT, NUOBJT )

C        L'ADRESSE DES COORDONNEES DU SOMMET SUIVANT
         MN = MN + 3

      ENDDO

      RETURN
      END
