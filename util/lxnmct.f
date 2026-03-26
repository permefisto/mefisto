      SUBROUTINE LXNMCT( NTLX , KNOM , NVTAMS )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : REMPLACER LE NUMERO DU TAMS DE NOM KNOM DANS LE LEXIQUE NTLX
C ----- PAR LE NUMERO NVTAMS
C
C ENTREES :
C ---------
C NTLX   : NUMERO DU TABLEAU TAMS DU LEXIQUE
C KNOM   : NOM DU TAMS DE NUMERO A MODIFIER
C NVTAMS : NOUVEAU NUMERO DU TAMS DE NOM KNOM
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS  OCTOBRE 1985
C.......................................................................
      CHARACTER*(*) KNOM
      include"./incl/pp.inc"
      COMMON        MCN(MOTMCN)
C
C     RECHERCE DU NUMERO DU NOM KNOM DANS LE LEXIQUE
      CALL LXNMNO( NTLX , KNOM , NOKNOM , MNLX )
C
C     MISE A JOUR DU NUMERO DU TAMS DE NOM KNOM
      MCN(MNLX+MCN(MNLX+2)+2+MCN(MNLX)*NOKNOM) = NVTAMS
      END
