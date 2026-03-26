      SUBROUTINE SDDEF1( KNOMSD , LISTE , KPB , KPBMTC )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    LECTURE AU NIVEAU DE L'OBJET DES DONNEES INTERNES
C -----    DANS LE CAS D'UNE RESOLUTION CLASSIQUE OU PAR SOUS-DOMAINES
C
C ENTREE :
C --------
C KNOMSD : NOM DE L'OBJET ASSOCIE AUX DONNEES
C LISTE  : LISTE DES MOTS CLES
C KPB    : PROVENANCE DU PROBLEME : ELASTICITE OU THERMIQUE
C KPBMTC : NOM DU TABLEAU DES MOTS CLES
C
C EXEMPLE D'APPELS:
C CALL SDDEF1( KNOMOB , LISTE1 , 'ELASTICITE', 'coefelas' )
C CALL SDDEF1( KNOMOB , LISTE1 , 'THERMIQUE',  'coefther' )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY     ANALYSE NUMERIQUE UPMC PARIS     FEVRIER 1990
C2345X7..............................................................012
      IMPLICIT          INTEGER (W)
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/a___young.inc"
      include"./incl/a___contact.inc"
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NFDOCU,NFFRAP,NUNITE(27)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      CHARACTER*16      KNOMSD
      CHARACTER*(*)     LISTE(1:*),KPB,KPBMTC
      CHARACTER*24      KNOMTD
      CHARACTER*80      NOMTS
C
      IERR = 0
C
C     LE NOM DU PROBLEME
C     ------------------
C     RECHERCHE DU 1-ER BLANC DANS LE NOM
      L0 = INDEX( KPB , ' ' )
      IF( L0 .GT. 0 ) THEN
         L0 = L0 - 1
      ELSE
         L0 = LEN( KPB )
      ENDIF
      IF ( KPB(1:L0) .EQ. 'ELASTICITE' ) THEN
         WTMTCL = WTYOUN
      ELSE IF ( KPB(1:L0) .EQ. 'THERMIQUE' ) THEN
         WTMTCL = WTCONT
      ELSE
          NBLGRC(NRERR) = 2
          IF( LANGAG .EQ. 0 ) THEN
             KERR(1) = 'ERREUR A L''APPEL DE SDDFE1'
          ELSE
             KERR(1) = 'ERROR in the CALLING of SDDFE1'
          ENDIF
          KERR(2) = ' KPB =' // KPB
          CALL LEREUR
          RETURN
      ENDIF
C
C     LE NOM DE L'OBJET
C     -----------------
C     RECHERCHE DU 1-ER BLANC DANS LE NOM
      L1 = INDEX( KNOMSD , ' ' )
      IF( L1 .GT. 0 ) THEN
         L1 = L1 - 1
      ELSE
         L1 = LEN( KNOMSD )
      ENDIF
C     IMPRESSION DU NOM DE L'OBJET
      NBLGRC(NRHIST) = 1
      IF( LANGAG .EQ. 0 ) THEN
         KHIST(1) = 'OBJET: ' // KNOMSD
      ELSE
         KHIST(1) = 'OBJECT: ' // KNOMSD
      ENDIF
      CALL LHISTO
      L2 = NUDCNB( KHIST(1) )
      WRITE( NFFRAP, * ) '{ ' , KHIST(1)(1:L2), ' }'
      WRITE( IMPRIM, * )
      WRITE( IMPRIM, * ) KHIST(1)(1:L2)
C
C     LE NOM DU TABLEAU
C     -----------------
 110  CALL LIMTCL( KPBMTC , NMTCL )
      IF( NMTCL .LE. 0 ) RETURN
C
C     LE NOM DU TABLEAU TMS A REMPLIR
      NOMTS = '~>OBJET>' // KNOMSD(1:L1) // '>' // LISTE(NMTCL)
      L2 = INDEX( NOMTS , ' ' )
      IF( L2 .GT. 0 ) THEN
         L2 = L2 - 1
      ELSE
         L2 = LEN( NOMTS )
      ENDIF
      KNOMTD = '~>>>' // LISTE(NMTCL)
      CALL MOTSTD( KNOMTD , NOMTS(1:L2) , IERR )
C
C     LE TABLEAU DOIT IL ETRE SUPPRIME ?
      CALL LXLXOU( NTMN(5) , KNOMSD , NTOB , MNTS )
      CALL LXTSOU( NTOB , LISTE(NMTCL) , NTTS , MNTS )
C     LE TYPE DE DONNEES DES CARACTERISTIQUES PHYSIQUES
      IF( NTTS .GT. 0 ) THEN
C
C        ATTENTION : TOUS LES TYPES DOIVENT ETRE A LA MEME ADRESSE
C                    DANS LE TABLEAU DESCRIPTEUR
C
         IF( IERR .NE. 0 .OR. MCN(MNTS+WTMTCL) .EQ. 0 ) THEN
C            LE TYPE EST NUL => DESTRUCTION DU TABLEAU
             CALL LXTSDS( NTOB , LISTE(NMTCL) )
             NBLGRC(NRERR) = 2
             KERR(1) = LISTE(NMTCL)
             IF( LANGAG .EQ. 0 ) THEN
                KERR(2) = 'CARACTERISTIQUE PHYSIQUE SUPPRIMEE'
             ELSE
                KERR(2) = 'PHYSICAL CHARACTERISTIC DELETED'
             ENDIF
             CALL LEREUR
             GOTO 110
         ENDIF
      ENDIF
C
C     AFFICHAGE DES VALEURS DU TABLEAU LU
      CALL AFTSTD( NOMTS(1:L2) )
      GOTO 110
C
      END
