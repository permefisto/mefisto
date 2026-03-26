       SUBROUTINE COUOBJ( NUTYOB , NMOBJT , MNTRAC )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    FIXER LES COULEURS DE L'OBJET DE TYPE NUTYOB DE NOM NMOBJT
C -----    DANS LA MESURE OU SON TABLEAU TRACE EXISTE
C ENTREES:
C --------
C NUTYOB : NUMERO DU TYPE DE L'OBJET A TRACER
C NMOBJT : NOM    DE L'OBJET A TRACER
C
C SORTIE :
C --------
C MNTRAC : ADRESSE MCN DU TABLEAU DE TRACE
C          0 S'IL N'EXISTE PAS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS      AVRIL 1990
C23456---------------------------------------------------------------012
      IMPLICIT INTEGER(W)
      include"./incl/ntmnlt.inc"
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      include"./incl/a___trace.inc"
      include"./incl/pp.inc"
      COMMON          MCN(MOTMCN)
      REAL           RMCN(1)
      EQUIVALENCE    (MCN(1),RMCN(1))
      CHARACTER*(*)   NMOBJT
C
C     LE LEXIQUE DE L'OBJET
      CALL LXLXOU( NTMN(NUTYOB) , NMOBJT , NTLXOB , MNLXOB )
C
C     LE TABLEAU TRACE EXISTE T IL ?
      CALL LXTSOU( NTLXOB , 'TRACE' , NTTRAC , MNTRAC )
C
      IF( NTTRAC .GT. 0 ) THEN
C
C        LE TABLEAU TRACE EXISTE  LES ATTRIBUTS SONT IMPOSES
C        =======================
         CALL TRTATA( MCN(MNTRAC+WBRCOU) , LATTRI , MLATTR )
C
      ENDIF
      IF( NDCOUL .EQ. N1COUL+1 ) THEN
C        ECRAN NOIR ET BLANC
         NBRCOU = 0
      ELSE
C        ECRAN COULEUR
         NBRCOU = NDCOUL - N1COUL + 1
      ENDIF
      END
