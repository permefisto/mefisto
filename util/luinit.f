      SUBROUTINE LUINIT
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : INITIALISER LES VARIABLES DU COMMON
C -----
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS    OCTOBRE 1989
C23456---------------------------------------------------------------012
      include"./incl/lu.inc"
      include"./incl/td.inc"
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      COMMON / MSSFTA / MSSF(28), NTADAM
C
      LHKLG  = 0
      NLPTV1 = 0
      NCPTV1 = 0
      NLPTV2 = 0
      NCPTV2 = 0
      LUIMPR = 1
      NBCTE  = 0
      NBPACH = 0
      NBVARU = 0
      NBVATS = 0
      NBETIQ = 0
      DO I=1,MXKLG
         KLG(I) = ' '
      ENDDO
C
C     DECLARATION DU TABLEAU DE LECTURE DES MOTS CLES
      CALL TNMCDC( 'REEL2', 1, MNMTCL )
C
cccC     DECLARATION DU TAMS  'LU'
ccc      CALL LXTNDC( NTADAM, 'LU', 'MOTS', 2 )
C
      RETURN
      END
