      SUBROUTINE SDRES1( NTLXOB , MNU , MNB , MNR , KNOMAT , IERR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCULER LE RESIDU DU SYSTEME LINEAIRE : R = B - A * X
C -----
C
C ENTREES :
C ---------
C NTLXOB : NUMERO  DU LEXIQUE DE L'OBJET SOUS-DOMAINE
C MNU    : ADRESSE DU VECTEUR SOLUTION
C MNB    : ADRESSE DU VECTEUR SECOND MEMBRE
C MNR    : ADRESSE DU VECTEUR RESIDU
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
CPU      DOUBLE PRECISION  DINFO,DCPU
      EQUIVALENCE      (MCN(1),RMCN(1),DMCN(1))
      CHARACTER*24      KNOM
      CHARACTER*(*)     KNOMAT
C
C     INITIALISATIONS
C     ---------------
CPU   DCPU   = DINFO( 'DELTA CPU' )
      MOREE2 = MOTVAR(6)
      IERR   = 0
      CALL TRTATA( MCN(MNB) , MCN(MNR) , NTDL * NDSM * MOREE2 )
C     LA MATRICE EST SYMETRIQUE
      NCODSA = 1
C
C     LA METHODE DE RESOLUTION
C     ========================
C
C     ***  ATTENTION : POUR RAISON D'EFFICACITE ( MULTIPLICATION D'UNE
C                      MATRICE PAR UNE VECTEUR ) , LA MATRICE UTILISEE
C                      DANS CE SP EST TOUJOURS STOCKEE SOUS FORME MORSE.
C                      VOIR AUSSI LES SP SDELA1 ET SDTHE1
C
CM      IF (NORESO.EQ.1) THEN
C
C        METHODE DE CHOLESKY :
C        -------------------
C        LE NOM DE LA MATRICE
CM         KNOM = 'PROFIL"' // KNOMAT
C        RECHERCHE DU 1-ER BLANC DANS LE NOM
CM         L0 = INDEX( KNOM , ' ' )
CM         IF( L0 .GT. 0 ) THEN
CM            L0 = L0 - 1
CM         ELSE
CM            L0 = LEN( KNOM )
CM         ENDIF
C        RECHERCHE DE LA MATRICE PROFIL SYMETRIQUE
CM         CALL LXTSOU( NTLXOB , KNOM(1:L0) , NTPROF , MNPROF )
CM         IAMUDL = MNPROF + WPDIAG
CM         LPPROF = MCN( MNPROF + WPPROF )
CM         MNAG   = MNPROF + LPPROF
CM         CALL PRBMAX(NCODSA,NDSM,NTDL,MCN(IAMUDL),MCN(MNAG),
CM     &              MCN(MNU),MCN(MNR))
C
CM      ELSE IF(NORESO.EQ.2) THEN
C
C        METHODE DU GRADIENT CONJUGUE
C        ----------------------------
C        LE NOM DE LA MATRICE
         KNOM = 'MORSE"' // KNOMAT
C        RECHERCHE DU 1-ER BLANC DANS LE NOM
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
         CALL MOBMAX(NCODSA,NDSM,NTDL,MCN(MNLPLI),MCN(MNLPCO),
     &                  MCN(MNAG),MCN(MNU),MCN(MNR))
CM      ENDIF
C
CPU   I = NINT( DINFO( 'DELTA CPU' ) )
CPU   WRITE(IMPRIM,12000) I
CCC12000 FORMAT(/'FIN NORMALE SDRES1 : COUT=',I12,' SECONDES CPU'/)
      RETURN
      END
