      SUBROUTINE LXLXDS( NTLX , KNOM )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : DETRUIRE LE LEXIQUE FILS DE NOM KNOM DANS LE LEXIQUE PERE
C ----- ET TOUS SES LEXIQUES,TABLEAUX,TABLEAUX DE TABLEAUX FILS OU
C       RELATIONS
C
C ENTREES :
C ---------
C NTLX   : NUMERO DU TABLEAU MS LEXIQUE PERE
C KNOM   : NOM DU LEXIQUE FILS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS    OCTOBRE 1985
C.......................................................................
      PARAMETER        (MXPILE=1024)
      include"./incl/motmcg.inc"
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      CHARACTER*(*)     KNOM
      CHARACTER*4       CHARX,K4,KTYPE
      CHARACTER*160     KNMLX
      COMMON / MSKNOM / KNMLX
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
C
C     RECHERCHE DU NOM KNOM PARMI LES NOMS DU LEXIQUE NTLX
      CALL LXNMNU( NTLX , KNOM , NONOM0 , NO , MNLX )
C
C     LE NOM A T IL ETE RETROUVE ?
      IF( NO .LE. 0 ) THEN
C
C        NON
C        ===
         NBLGRC(NRERR) = 2
         KERR(1) = 'LEXIQUE INCONNU '//KNOM
         KERR(2) = 'PARMI'
         CALL LEREUR
         CALL LXIM( NTLX )
         RETURN
      ENDIF
C
C     OUI
C     ===
      MN     = MNLX + MCN( MNLX ) * NO + MCN( MNLX+2 ) + 1
      KTYPE  = CHARX( MCN(MN) )
      IF( KTYPE .NE. 'LEXI' ) THEN
C        TYPE INCORRECT =/ LEXI
         NBLGRC(NRERR) = 1
         KERR(1) = KNOM
         KERR(2) = 'N''EST PAS UN LEXIQUE'
         CALL LEREUR
         RETURN
      ENDIF
C
C     LE LEXIQUE FILS NTLXLX EST PARCOURU
C     LES TAMS OU TATA SONT IMMEDIATEMENT DETRUITS
C     LE NO DES TAMS DES LEXIQUES SONT EMPILES
C     LEUR NO EST RENDU <0 LORSQU'ILS ONT ETE TRAITES MAIS
C     RESTENT A DETRUIRE
C
C     DECLARATION DE LA PILE
      MNPILE = -1
      CALL TNMCDC( 'ENTIER' , MXPILE , MNPILE )
      IF( MNPILE .LE. 0 ) RETURN
C
C     LE LEXIQUE NTLXLX=MCN(MN+1) EST EMPILE
      LHPILE = 1
      MCN(MNPILE) = MCN( MN + 1 )
C
C     TANT QUE LA PILE EST NON VIDE FAIRE
C     ===================================
C     LA PILE EST ELLE VIDE ?
 20   IF( LHPILE .GT. 0 ) THEN
C
C        NON.LECTURE DU HAUT DE LA PILE
         MNP = MNPILE + LHPILE - 1
C        LE NUMERO DU TAMS DU LEXIQUE A TRAITER
         NTL = MCN( MNP )
C
C        LE LEXIQUE RETROUVE A T IL DEJA ETE TRAITE ?
         IF( NTL .LT. 0 ) THEN
C           OUI.IL EST DETRUIT
            CALL TAMSDS( - NTL )
C           IL EST DEPILE
            LHPILE = LHPILE - 1
            GOTO 20
         ENDIF
C
C        LE LEXIQUE DOIT ETRE EXPLORE.IL EST MARQUE POUR ETRE
C        DETRUIT ENSUITE
         MCN(MNP) = -NTL
C
C        RECUPERATION DE SON ADRESSE
         CALL TAMSOU( NTL , MNL )
C
C        PARCOURS DES NOMS OCCUPES DU LEXIQUE (M1LX,0:M2LX)
         M1LX   = MCN( MNL )
         NBENNM = MCN( MNL + 2 )
C        NUMERO DU 1-ER NOM OCCUPE
         I      = MCN( MNL + 5 )
C
C        LA BOUCLE SUR LES NOMS OCCUPES
C        ------------------------------
 30      IF( I .GT. 0 ) THEN
C           LE NOM I EST OCCUPE. ADRESSE DE SES INFORMATIONS
            MN    = MNL + M1LX * I + NBENNM
C           LE NUMERO DU TMS ASSOCIE AU I-EME NOM
            NOTAB = MCN( MN + 2 )
C           QUEL EST LE TYPE DU TMS ?
            K4 = CHARX( MCN( MN + 1 ) )
            IF( K4 .EQ. 'TAMS' .OR. K4 .EQ. 'RELA' ) THEN
C              DESTRUCTION DU TABLEAU MS
               CALL TAMSDS( NOTAB )
            ELSE IF( K4 .EQ. 'TATA' ) THEN
C              DESTRUCTION DU TABLEAU DE TABLEAUX
C              REGENERATION DE SON NOM A PARTIR DU LEXIQUE
               CALL ENTNOM( NBENNM,MCN(MN-NBENNM),KNMLX )
               CALL LXTTDS( NTL,KNMLX )
            ELSE IF( K4 .EQ. 'LEXI' ) THEN
C              LE NOUVEAU LEXIQUE EST EMPILE
               IF( LHPILE .GE. MXPILE ) THEN
C                 DEBORDEMENT DE LA PILE
                  NBLGRC(NRERR) = 1
                  WRITE(KERR(MXLGER)(1:10),'(I10)') MXPILE
                  KERR(1) = 'LXLXDS:PILE SATUREE MALGRE '//
     %                       KERR(MXLGER)(1:10) //' MOTS'
                  CALL LEREUR
                  RETURN
               ELSE
C                 LE LEXIQUE EST EMPILE
                  IF( NOTAB .GT. 0 ) THEN
                     MCN( MNPILE + LHPILE ) = NOTAB
                     LHPILE = LHPILE + 1
                  ENDIF
               ENDIF
            ELSE
C              AUCUN DES 4 TYPES
               NBLGRC(NRERR) = 1
               KERR(1) = 'AUCUN TYPE RECONNU POUR '//KNOM
               CALL LEREUR
               RETURN
            ENDIF
C
C           PASSAGE AU SUIVANT
            I = MCN( MN )
            GOTO 30
         ENDIF
C
C        FIN DE LA BOUCLE SUR LES NOMS OCCUPES
C        -------------------------------------
         GOTO 20
      ENDIF
C
C     LA PILE EST VIDE ET DETRUITE
C     ============================
      CALL TNMCDS( 'ENTIER' , MXPILE , MNPILE )
C
C     MISE A JOUR DU CHAINAGE DANS LE LEXIQUE NTLX
      M1LX   = MCN( MNLX )
      NBENNM = MCN( MNLX + 2 )
      MN     = MNLX + M1LX * NO + NBENNM
      IF( NONOM0 .GT. 0 ) THEN
C        IL EXISTE UN PRECEDENT
         MCN( MNLX + M1LX * NONOM0 + NBENNM ) = MCN( MN )
      ELSE
C        NO EST TETE DE LISTE. SON SUIVANT DEVIENT TETE DE LISTE
         MCN( MNLX + 5 ) = MCN( MN )
      ENDIF
C
C     LE NOM KNOM DANS LE LEXIQUE REDEVIENT LIBRE
      MCN( MN     ) = MCN( MNLX + 6 )
      MCN( MN + 1 ) = ICHARX( 'V DE' )
      MCN( MN + 2 ) = 0
C
C     NO DEVIENT TETE DU CHAINAGE DES VIDES
      MCN( MNLX + 6 ) = NO
      END
