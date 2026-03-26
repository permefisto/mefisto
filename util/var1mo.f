      REAL FUNCTION VAR1MO( NUMTYV )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : RETOURNER LE NOMBRE DE VARIABLES DE NUMERO DE TYPE NUMTYV
C ----- CONTENUES DANS 1 MOT MEMOIRE
C **********************************************************************
C ATTENTION : CETTE FONCTION EST DEPENDANTE MACHINE ICI VERSION IBM
C **********************************************************************
C ENTREE :
C --------
C NUMTYV : NUMERO COMPRIS ENTRE 1 ET 9 AVEC LA CORRESPONDANCE
C          LOGIQUE  <= 1 CARACTERE<= 2 ENTIER/2 <= 3 ENTIER   <= 4
C          REEL     <= 5 REEL2    <= 6 REEL4    <= 7 COMPLEXE <= 8
C          COMPLEXE2<= 9
C
C SORTIE :
C --------
C VAR1MO : NOMBRE DE VARIABLES CONTENUES DANS 1 MOT MEMOIRE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : PERRONNET ALAIN ANALYSE NUMERIQUE PARIS  NOVEMBRE 1983
C.......................................................................
      include"./incl/msvaau.inc"
C
C     AFFECTATION DIRECTE
      VAR1MO = VAD1MO( NUMTYV )
      END
