      SUBROUTINE SDRES20( NTLXOB , MNAD , KNOMAT )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : EXTRAIRE LA DIAGONALE DE LA MATRICE LOCALE
C -----
C
C ENTREES :
C ---------
C NTLXOB : NUMERO  DU LEXIQUE DE L'OBJET SOUS-DOMAINE
C MNAD   : ADRESSE MCN DU TABLEAU DIAGONALE
C KNOMAT : NOM DE LA MATRICE A RECHERCHER ( RAIDEUR OU CONDUCTIVITE )
C
C SORTIES :
C ---------
C IERR   : 0 SI PAS D'ERREUR , NON NUL SINON
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY     ANALYSE NUMERIQUE UPMC PARIS  MAI 1990
C23456---------------------------------------------------------------012
      IMPLICIT          INTEGER (W)
      include"./incl/gsmenu.inc"
      include"./incl/a___profil.inc"
      include"./incl/a___morse.inc"
C
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      COMMON / SOUSDO / NORESO,NTDL,NDSM,NDIM,NPIMAX
      REAL              RMCN(1)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1),DMCN(1))
      CHARACTER*24      KNOM
      CHARACTER*(*)     KNOMAT
C
C     LE NOM DE LA MATRICE
C     --------------------
      KNOM = 'MORSE"' // KNOMAT
C     RECHERCHE DU 1-ER BLANC DANS LE NOM
      L0 = INDEX( KNOM , ' ' )
      IF( L0 .GT. 0 ) THEN
         L0 = L0 - 1
      ELSE
         L0 = LEN( KNOM )
      ENDIF
      CALL LXTSOU( NTLXOB , KNOM(1:L0) , NTMORS , MNMORS )
      LPMORS = MCN(MNMORS + WPMORS)
      MNAG   = MNMORS + LPMORS
      MNLPLI = MNMORS + WPLIGN
      MNLPCO = MNLPLI + NTDL + 1
C
C     LA DIAGONALE DE PRECONDITIONNEMENT
C     ----------------------------------
      IAG = ( MNAG - 1 ) / 2
      IAD = ( MNAD - 1 ) / 2
      DO K = 1 , NTDL
         DMCN(IAD+K) = DMCN(IAG+MCN(MNLPLI+K))
      ENDDO
C
      RETURN
      END
