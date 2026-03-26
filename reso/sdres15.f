      SUBROUTINE SDRES15( NTLXOB , MNR , MNZ , KNOMAC )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : PRECONDITIONNEMENT DU RESIDU (METHODE DES JOINTS) : DESCENTE
C -----
C
C ENTREE :
C ---------
C NTLXOB : NUMERO  DU LEXIQUE DE L'OBJET SOUS-DOMAINE
C MNR    : ADRESSE DU VECTEUR RESIDU
C
C SORTIES :
C ---------
C MNZ    : ADRESSE DU VECTEUR RESIDU PRECONDITIONNE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY     ANALYSE NUMERIQUE UPMC PARIS  MARS 1994
C23456---------------------------------------------------------------012
      IMPLICIT          INTEGER (W)
      include"./incl/gsmenu.inc"
      include"./incl/a___profil.inc"
      include"./incl/a___morse.inc"
C
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      COMMON / SOUSDO / NORESO,NTDL,NDSM,NDIM,NPIMAX
      CHARACTER*24      KNOM
      CHARACTER*(*)     KNOMAC
C
C     LE NOM DE LA MATRICE
C     --------------------
      KNOM = 'MORSE"' // KNOMAC
C     RECHERCHE DU 1-ER BLANC DANS LE NOM
      L0 = INDEX( KNOM , ' ' )
      IF( L0 .GT. 0 ) THEN
         L0 = L0 - 1
      ELSE
         L0 = LEN( KNOM )
      ENDIF
      CALL LXTSOU( NTLXOB , KNOM(1:L0) , NTMORS , MNMORS )
      LPMORS = MCN(MNMORS + WPMORS)
      MNAGC  = MNMORS + LPMORS
      MNLPLI = MNMORS + WPLIGN
      MNLPCO = MNLPLI + NTDL + 1
C
C     PRECONDITIONNEMENT : DESCENTE
C     -----------------------------
      CALL DECHGC(NTDL,MCN(MNLPLI),MCN(MNLPCO),MCN(MNAGC),
     S                 MCN(MNR),MCN(MNZ))
C
      RETURN
      END
