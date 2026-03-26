      SUBROUTINE LUFE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :      SAUVEGARDER LE COMMON / KLANUT / ET / ILANUT /
C -----      LES CONSTANTES ET VARIABLES DU LANGAGE UTILISATEUR
C            DANS LE TMS ADAM>LU
C            CF $MEFISTO/g/grammaire_lu DE DEFINITION DU LANGAGE UTILISATEUR
C            CF $MEFISTO/incl/lu.inc
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS  JANVIER 1990
C23456---------------------------------------------------------------012
      include"./incl/motmcg.inc"
      include"./incl/nbcamo.inc"
      include"./incl/lu.inc"
      include"./incl/a_lu.inc"
      include"./incl/gsmenu.inc"
      include"./incl/darete.inc"
      include"./incl/pp.inc"
      COMMON             MCN(MOTMCN)
      COMMON / MSSFTA /  MSSF(28), NTADAM
      COMMON / UNITES /  LECTEU, IMPRIM, INTERA, NUNITE(29)

C     LA TAILLE D'ARETE PAR DEFAUT DARETE EST SAUVEGARDEE DANS DCTE(NBCTE+1)
      IF( NBCTE .GT. 0 ) THEN
         NBCTE = NBCTE + 1
         DCTE( NBCTE ) = DARETE
      ENDIF
C
C     LE TABLEAU LU EXISTE T IL? SI OUI IL EST DETRUIT
      CALL LXTSOU( NTADAM, 'LU', NTLU, MNLU )
      IF( NTLU .GT. 0 ) THEN
C        DESTRUCTION DU TABLEAU LU
         CALL LXTSDS( NTADAM, 'LU' )
      ENDIF
C
C     DECLARATION DU TABLEAU ADAM>LU
C     ==============================
C     LA DIVISION DOIT ETRE EXACTE SINON PB
      MXCAVA = NBCAVA / NBCAMO
      MXCATS = NBCATS / NBCAMO
      MXCACH = NBCACH / NBCAMO
      MOTS = WDCTE + (NBCTE + NBVARU)*MOTVAR(6) +
     %       MXCAVA*NBVARU + MXCATS*NBVATS + MXCACH*NBPACH
      CALL TNMCMX( 'MOTS', MAXVAR )
      IF( MOTS .GT. MAXVAR ) GOTO 9999
      CALL LXTNDC( NTADAM, 'LU', 'MOTS', MOTS )
C
C     OUVERTURE DU TABLEAU LU
      CALL LXTSOU( NTADAM, 'LU', NTLU, MNLU )
      IF( NTLU .LE. 0 ) GOTO 9999
C
C     INITIALISATION DU TABLEAU LU
      MCN( MNLU + WNBCTE ) = NBCTE
      MCN( MNLU + WBVARU ) = NBVARU
      MCN( MNLU + WXCAVA ) = MXCAVA
      MCN( MNLU + WBVATS ) = NBVATS
      MCN( MNLU + WXCATS ) = MXCATS
      MCN( MNLU + WBPACH ) = NBPACH
      MCN( MNLU + WXCACH ) = MXCACH
C
C     LES CONSTANTES
      MN   = MNLU  + WDCTE
      MOTS = NBCTE * MOTVAR(6)
      CALL TRTATA( DCTE, MCN(MN), MOTS )
C
C     LES VALEURS DES VARIABLES UTILISATEUR
      MN   = MN + MOTS
      MOTS = NBVARU * MOTVAR(6)
      CALL TRTATA( DVARU, MCN(MN), MOTS )
C
C     LES MOTS DE NBCAMO CARACTERES DES NOMS DE VARIABLE UTILISATEUR
      MN = MN + MOTS
      DO 20 I = 1, NBVARU
         I1 = 1
         I2 = NBCAMO
         DO 10 MOTS = 1, MXCAVA
            MCN(MN) = ICHARX( KVARU(I)(I1:I2) )
            MN = MN + 1
            I1 = I2 + 1
            I2 = I2 + NBCAMO
 10      CONTINUE
 20   CONTINUE
C
C     LES MOTS DE NBCAMO CARACTERES DES NOMS DE VARIABLE TMS
      DO 40 I = 1, NBVATS
         I1 = 1
         I2 = NBCAMO
         DO 30 MOTS = 1, MXCATS
            MCN(MN) = ICHARX( KVATS(I)(I1:I2) )
            MN = MN + 1
            I1 = I2 + 1
            I2 = I2 + NBCAMO
 30      CONTINUE
 40   CONTINUE
C
C     LES MOTS DE NBCAMO CARACTERES DES NOMS DE PARAMETRES CHAINES
      DO 60 I = 1, NBPACH
         I1 = 1
         I2 = NBCAMO
         DO 50 MOPC = 1, MXCACH
            MCN(MN) = ICHARX( KPACH(I)(I1:I2) )
            MN = MN + 1
            I1 = I2 + 1
            I2 = I2 + NBCAMO
 50      CONTINUE
 60   CONTINUE
C
C     AJOUT DE LA DATE
      CALL ECDATE( MCN(MNLU) )
C     AJOUT DU NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNLU + MOTVAR(6) ) = NONMTD( '~>LU' )
      RETURN
C
C     ERREUR
 9999 NBLGRC(NRERR) = 2
      KERR(1) = 'IMPOSSIBLE CREER TMS LU'
      KERR(2) = 'ARRET EN CATASTROPHE'
      CALL LEREUR
      STOP 'ARRET DANS LUFE'
      END
