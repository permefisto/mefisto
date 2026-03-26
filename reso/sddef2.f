      SUBROUTINE SDDEF2( NUOBJE, NUTYOB, NUMOBJ, LISTE, LLISTE,
     %                   KPB,    KPBMTC, NBJEUX )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    LECTURE DES CARACTERISTIQUES PHYSIQUES DE CHAQUE MATERIAU
C -----    cad du VOLUME en 3D, de la SURFACE en 2D, de la LIGNE en 1D
C          ( THERMIQUE, ELASTICITE, FLUIDE )
C          INPUT DATA of the CHARACTERISTICS of EVERY MATERIAL

C ENTREE :
C --------
C NUOBJE : NUMERO DE L'OBJET DANS LE LEXIQUE DES OBJETS
C NUTYOB : NUMERO DU TYPE VOLUME/SURFACE
C NUMOBJ : NUMERO DANS LE LEXIQUE DES OBJETS DE CE TYPE
C LISTE  : LISTE DES MOTS CLES
C LLISTE : LONGUEUR DE LA LISTE DES MOTS CLES
C KPB    : PROVENANCE DU PROBLEME : ELASTICITE THERMIQUE NLSE ou FLUIDE
C KPMTC  : NOM DU TABLEAU DES MOTS CLES
C NBJEUX : NOMBRE DE JEUX DE DONNEES
C          1 EN ETUDE STANDARD, >1 POUR LES POLYNOMES DE VALEURS PROPRES
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL    JOLY  ANALYSE NUMERIQUE UPMC PARIS     FEVRIER 1990
C MODIFS : PERRONNET ALAIN TIMS NTU TAIPEI TAIWAN          NOVEMBRE 2009
C2345X7..............................................................012
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NFDOCU,NFFRAP,NUNITE(27)
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/a___young.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
C     LU A DECLARER     LU(LLISTE)
      LOGICAL           LU(9)
      CHARACTER*(*)     LISTE(1:LLISTE), KPB, KPBMTC
      CHARACTER*10      NMTYOB, KNOMTY, KTYPE
      CHARACTER*24      KNOMOB, KNOMTD
      CHARACTER*80      NOMTS
      CHARACTER*1       KNOJEU
      CHARACTER*2       KAJOUT
      CHARACTER*24      KLISTE

      IERR = 0
      IF( NBJEUX .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'NOMBRE DE JEUX DE DONNEES <=0'
         ELSE
            KERR(1) = 'NUMBER of INPUT DATA GAMES <=0'
         ENDIF
         CALL LEREUR
         RETURN
      ENDIF
      print*,
     %'PROVENANCE DU PROBLEME ELASTICITE THERMIQUE NLSE ou FLUIDE: ',KPB

C     LECTURE DES DONNEES
C     -------------------
C     LE TYPE DE L'OBJET
      KNOMTY = NMTYOB( NUTYOB )
      L1 = INDEX( KNOMTY, ' ' )
      IF( L1 .GT. 0 ) THEN
         L1 = L1 - 1
      ELSE
         L1 = LEN( KNOMTY )
      ENDIF
C     SON NOM
      CALL NMOBNU( KNOMTY(1:L1), NUMOBJ, KNOMOB )
      L2 = INDEX( KNOMOB, ' ' )
      IF( L2 .GT. 0 ) THEN
         L2 = L2 - 1
      ELSE
         L2 = LEN( KNOMOB )
      ENDIF

C     IMPRESSION D'UN MESSAGE
C     -----------------------
C     AJOUT DU NOM DU PLSV
      NBLGRC(NRHIST) = 1
      KHIST(1) = KNOMTY(1:L1) // ': ' // KNOMOB(1:L2)
      CALL LHISTO
      L3 = NUDCNB( KHIST(1) )
      WRITE( NFFRAP, * ) '{ ', KHIST(1)(1:L3), ' }'
      WRITE( IMPRIM, * )
      WRITE( IMPRIM, * ) KHIST(1)(1:L3)

C     BOUCLE SUR LES JEUX DE DONNEES
C     ==============================
      DO 100 NOJEU = 1, NBJEUX
         WRITE( KNOJEU, '(I1)' ) NOJEU

ccc         WRITE(IMPRIM,*)
ccc         NBLGRC(NRERR) = 1
ccc         IF( LANGAG .EQ. 0 ) THEN
ccc            KERR(1)='ATTENTION: JEU de DONNEES Numero='//KNOJEU
ccc         ELSE
ccc            KERR(1)='ATTENTION: INPUT DATA GAME Number='//KNOJEU
ccc         ENDIF
ccc         CALL LERESU

C        LECTURE AU NIVEAU VOLUME/SURFACE/LIGNE
C        --------------------------------------
         DO 10 NLI=1,LLISTE
            LU(NLI) = .FALSE.
 10      ENDDO

 30      CALL LIMTCL( KPBMTC, NMTCL )
         IF( NMTCL .GT. 0 ) THEN
C           LE NOM DU TABLEAU TMS A REMPLIR
            NOMTS = '~>' // KNOMTY(1:L1) // '>' // KNOMOB(1:L2) //
     %               '>' // LISTE(NMTCL)
            L3 = INDEX( NOMTS, ' ' )
            IF( L3 .GT. 0 ) THEN
               L3 = L3 - 1
            ELSE
               L3 = LEN( NOMTS )
            ENDIF

            KNOMTD = '~>>>' // LISTE(NMTCL)
            CALL MOTSTD( KNOMTD, NOMTS(1:L3), IERR )

C           LE TABLEAU EST IL CREE ?
            CALL LXLXOU( NTMN(NUTYOB), KNOMOB, NTOB, MNTS )
            CALL LXTSOU( NTOB, LISTE(NMTCL), NTTS, MNTS )
            IF( NTTS .GT. 0 ) THEN

               IF( NBJEUX .GT. 1 ) THEN
C                 AJOUT DE "NOJEU au NOM du TMS si NBJEUX>1
                  KAJOUT = '"' // KNOJEU
                  NOMTS(L3+1:L3+2) = KAJOUT
                  L3 = L3+2
C                 LE TABLEAU EST RENOMME AVEC "NOJEU
                  L4 = NUDCNB( LISTE(NMTCL) )
                  KLISTE = LISTE(NMTCL)(1:L4) // KAJOUT
                  CALL LXNMNM( NTOB, LISTE(NMTCL),  KLISTE )
               ELSE
                  KAJOUT = '  '
               ENDIF

               L4 = NUDCNB( LISTE(NMTCL) )
               KLISTE = LISTE(NMTCL)(1:L4) // KAJOUT
               L4 = NUDCNB( KLISTE )
               CALL LXTSOU( NTOB, KLISTE(1:L4), NTTS, MNTS )

C              LE TYPE DE DONNEES DES CARACTERISTIQUES PHYSIQUES
C              ATTENTION : TOUS LES TYPES DOIVENT ETRE A LA MEME ADRESSE
C                          DANS LE TABLEAU DESCRIPTEUR

               IF( IERR .NE. 0 .OR. MCN(MNTS+WTYOUN) .EQ. 0 ) THEN
C                 LE TYPE EST NUL => DESTRUCTION DU TABLEAU
                  CALL LXTSDS( NTOB, KLISTE(1:L4) )
                  NBLGRC(NRERR) = 2
                  KERR(1)(1:L3) = NOMTS(1:L3)
                  IF( LANGAG .EQ. 0 ) THEN
                  KERR(2)='Cette CARACTERISTIQUE PHYSIQUE est SUPPRIMEE'
                  ELSE
                     KERR(2) = 'This PHYSICAL CHARACTERISTIC is DELETED'
                  ENDIF
                  CALL LERESU
                  GOTO 30
               ENDIF
            ENDIF

C           AFFICHAGE DES VALEURS DU TABLEAU LU
            IF( IERR .EQ. 0 ) CALL AFTSTD( NOMTS(1:L3) )
            LU(NMTCL) = .TRUE.
            GOTO 30
         ENDIF
 100  ENDDO

C     ESSAI DE LECTURE AU NIVEAU DE L'OBJET  (1 seul JEU)
C     -------------------------------------
      NUTYSD = 5
C     LE NUMERO DE TMS DU LEXIQUE DE L'OBJET
      CALL LXNLOU( NTMN(NUTYSD), NUOBJE, NTSD, MNLXSD)
C     LE NUMERO DE TMS DU LEXIQUE DU VOLUME/SURFACE
      NTLXOB   = NTMN( NUTYOB )
      CALL LXNLOU( NTLXOB, NUMOBJ, NTOB, MNLXOB )
      DO 200 NLI = 1 , LLISTE
         IF (LU(NLI)) GO TO 200
         NOMTS =  LISTE(NLI)
         L = INDEX( NOMTS, ' ' )
         IF( L .GT. 0 ) THEN
            L = L - 1
         ELSE
            L = LEN( NOMTS )
         ENDIF
C        RECHERCHE DU TABLEAU DE NOM NOMTS DANS LE LEXIQUE DE L'OBJET
         CALL LXTSOU( NTSD, NOMTS(1:L), NTTAB1, MNTAB1 )
C        SI LE TABLEAU EXISTE ON LE RECOPIE
         IF( NTTAB1 .GT. 0 ) THEN
C           SA LONGUEUR
            CALL TAMSTV( NTTAB1, KTYPE, LTAB )
C           DECLARATION DU NOUVEAU TABLEAU DANS LE LEXIQUE
            CALL LXTSOU( NTOB, NOMTS(1:L),  NTTAB2, MNTAB2 )
            IF( NTTAB2 .GT. 0 ) THEN
C              DESTRUCTION AVANT REGENERATION
               CALL LXTSDS( NTOB, NOMTS(1:L) )
            ENDIF
            CALL LXTNDC( NTOB, NOMTS(1:L), 'MOTS' , LTAB )
            CALL LXTSOU( NTOB, NOMTS(1:L),  NTTAB2, MNTAB2 )
C           COPIE DU TABLEAU1 DANS LE TABLEAU2
            CALL TRTATA( MCN(MNTAB1), MCN(MNTAB2), LTAB )
C           MISE A JOUR DE LA DATE
            CALL ECDATE( MCN(MNTAB2) )
         ENDIF
 200  ENDDO

      RETURN
      END
