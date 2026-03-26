      SUBROUTINE LXPERE( NTFILS , NTPERE )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : RETOURNER LE NUMERO DE TAMS DU LEXIQUE PERE DU LEXIQUE DE TAMS
C ----- DE NUMERO NTFILS
C       REMONTER DANS L'ARBRE DES LEXIQUES
C
C ENTREE :
C -------
C NTFILS : NUMERO DU TAMS DU LEXIQUE FILS
C
C SORTIE :
C --------
C NTPERE : NUMERO DU TAMS DU LEXIQUE PERE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS  OCTOBRE 1985
C.......................................................................
      include"./incl/pp.inc"
      COMMON MCN(MOTMCN)
C
C     OUVERTURE DU LEXIQUE FILS
      CALL TAMSOU( NTFILS , MNFILS )
C
C     LE NUMERO DU TAMS DU LEXIQUE PERE
      NTPERE = MCN( MNFILS + 4 )
      END
