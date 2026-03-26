      SUBROUTINE LXIM0( MNLX )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    AFFICHER LES NOMS CONTENUS DANS LE LEXIQUE D'ADRESSE MCN MNLX
C -----    A CONCURRENCE DE LA PLACE MAXIMALE POSSIBLE DANS KERR
C
C ENTREE :
C --------
C MNLX   : ADRESSE MCN DU TABLEAU MS CONTENANT LE LEXIQUE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS          OCTOBRE 1984
C2345X---------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/msvaau.inc"
      include"./incl/gsmenu.inc"
      COMMON / MSSFTA / NOFISF,NBPASF,MOPASF,MGBUSF,NSFLIB,
     %                  M1FIMS,M2FIMS,MGFIMS,NSFIMS,LPFIMS,
     %                  M1TAMS,M2TAMS,MGBUTA,NBBUTA,NPTAMS,NATAMS,
     %                  NBCTMS,LLTAMS,LFTAMS,MGNPSF,NSFNPS,NPSNPS,
     %                  MGZLMG,MGZLMK,MGZLMN,MOTSMG,MOTSMK,MOTSMN,NTADAM
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      COMMON / MSKNOM / KNOM
      CHARACTER*160     KNOM
C
C     MNLX : ADRESSE MCN DU LEXIQUE SUPPOSE OUVERT
      IF( MNLX .LE. 0 .OR. MNLX .GT. MOTSMN ) RETURN
C     M1LX : NOMBRE D ENTIERS PAR NOM ET ATTRIBUTS DU LEXIQUE
      M1LX   = MCN( MNLX )
C     NBENNM : NOMBRE D'ENTIERS POUR STOCKER LES CARACTERES D'UN NOM
      NBENNM = MCN( MNLX + 2 )
C     NBCANM : NOMBRE DE CARACTERES D'UN NOM DU LEXIQUE
C     NBCANM = MCN( MNLX + 3 )
C
C     N NUMERO DU 1-ER NOM OCCUPE DU LEXIQUE
      N = MCN( MNLX + 5 )
      IF( N .LE. 0 ) THEN
C        LEXIQUE VIDE
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'lxim0: LEXIQUE VIDE'
         ELSE
            KERR(1) = 'lxim0: EMPTY LEXICON'
         ENDIF
         CALL LERESU
         RETURN
      ENDIF
C
C     PASSAGE A LA LIGNE
      NL0    = NBLGRC(NRERR)
      NCAR   = 0
      NL1    = NL0
      NL2    = 0
      NCAR1  = 0
      MXLGE1 = MXLGER -1
C
C     LE NOM EST IL OCCUPE ?
 10   IF( N .GT. 0 ) THEN
C        OUI.TRANSFORMATION DES NBENNM ENTIERS EN CARACTERES
         MN = MNLX + M1LX * N
         IF( MN .GT. MOTSMN ) RETURN
         CALL ENTNOM( NBENNM, MCN(MN), KNOM )
C        IMPRESSION DU NOM
 20      IF( NL1 .LT. MXLGE1 ) THEN
C           MXLGER-1 POUR AJOUTER LA LIGNE CLIQUER POUR CONTINUER...
            NL1 = NL1 + 1
            N   = NUDCNB( KNOM )
            IF( NCAR + N .GT. NBCAER ) THEN
C              TROP DE NOMS
               KERR(MXLGER)(NBCAER-22:NBCAER)='... LISTE ECOURTEE ...'
               GOTO 9000
            ENDIF
C           LE CARACTERE MAX ATTEINT
            NCAR1 = MAX( NCAR + N , NCAR1 )
            KERR( NL1 )(1+NCAR:NBCAER) = KNOM(1:N)
C           PASSAGE AU SUIVANT
            N = MCN( MN + NBENNM )
            GOTO 10
         ELSE
C           DEBORDEMENT EN LIGNES
            NL2   = MXLGER
            NL1   = NL0
            NCAR  = NCAR1 + 1
            NCAR1 = 0
            GOTO 20
         ENDIF
      ENDIF
C
C     L'AFFICHAGE
      IF( NL2 .EQ. 0 ) NL2 = NL1
 9000 NBLGRC(NRERR) = NL2
      CALL LERESU
      NBLGRC(NRERR) = 0

      RETURN
      END
