      SUBROUTINE LXLXFE( NTLX , KNOM )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : FERMER LE LEXIQUE FILS DE NOM KNOM DU LEXIQUE PERE
C -----        TOUS SES FILS C'EST A DIRE TOUS LES
C       LEXIQUES,TABLEAUX,TABLEAUX DE TABLEAUX ,RELATIONS  FILS
C
C ENTREES :
C ---------
C NTLX   : NUMERO DU TABLEAU MS DU LEXIQUE PERE
C KNOM   : NOM DU LEXIQUE FILS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS     OCTOBRE 1985
C.......................................................................
      PARAMETER       (MXPILE=1024)
      include"./incl/motmcg.inc"
      include"./incl/gsmenu.inc"
      CHARACTER*(*)    KNOM
      CHARACTER*4      CHARX,K4,KTYPE
      CHARACTER*160    KNMLX
      COMMON / MSKNOM/ KNMLX
      include"./incl/pp.inc"
      COMMON           MCN(MOTMCN)
      COMMON / UNITES/ LECTEU,IMPRIM,INTERA,NUNITE(29)
C
C     RECHERCHE DU NOM KNOM PARMI LES NOMS DU LEXIQUE NTLX
      CALL LXNMNO( NTLX , KNOM , NO , MNLX )
C
C     LE NOM A T IL ETE RETROUVE ?
      IF( NO .LE. 0 ) THEN
C        NON
         NBLGRC(NRERR) = 1
         KERR(1) = 'NOM INCONNU: ' // KNOM
         CALL LEREUR
         CALL LXIM( NTLX )
         RETURN
      ENDIF
C
C     OUI
      MN     = MNLX + MCN( MNLX ) * NO + MCN( MNLX+2 ) + 1
      KTYPE  = CHARX( MCN(MN) )
      IF( KTYPE .NE. 'LEXI' ) THEN
C        TYPE INCORRECT =/ LEXI
         NBLGRC(NRERR) = 2
         KERR(1) = KNOM
         KERR(2) = 'N''EST PAS UN LEXIQUE'
         CALL LEREUR
         RETURN
      ENDIF
C
C     LE LEXIQUE FILS DE NOM KNOM A POUR NUMERO DE TAMS NTLXLX
      NTLXLX = MCN( MN + 1 )
C
C     LE LEXIQUE FILS NTLXLX EST PARCOURU
C     LES TAMS OU TATA SONT IMMEDIATEMENT FERMES
C     LE NO DES TAMS DES LEXIQUES SONT EMPILES
C     LEUR NO EST RENDU <0 LORSQU'ILS ONT ETE TRAITES MAIS
C     RESTENT A FERMER
C
C     DECLARATION DE LA PILE
      MNPILE = 0
      CALL TNMCDC( 'ENTIER',MXPILE,MNPILE )
C
C     LE LEXIQUE NTLXLX EST EMPILE
      LHPILE = 1
      MCN(MNPILE) = NTLXLX
C
C     TANT QUE LA PILE EST NON VIDE FAIRE
C     ===================================
C     LA PILE EST ELLE VIDE ?
 20   IF( LHPILE .LE. 0 ) GOTO 9000
C
C        NON.LECTURE DU HAUT DE LA PILE
         MNP = MNPILE + LHPILE - 1
C        LE NUMERO DU TAMS DU LEXIQUE A TRAITER
         NTL = MCN( MNP )
C
C        LE LEXIQUE RETROUVE A T IL DEJA ETE TRAITE ?
         IF( NTL .LT. 0 ) THEN
C           OUI.IL EST FERME
            CALL TAMSFE( - NTL )
C           IL EST DEPILE
            LHPILE = LHPILE - 1
            GOTO 20
         ENDIF
C
C        LE LEXIQUE DOIT ETRE EXPLORE.IL EST MARQUE POUR ETRE
C        FERME ENSUITE
         MCN(MNP) = -NTL
C
C        RECUPERATION DE SON ADRESSE
         CALL TAMSOU( NTL , MNL )
C
C        PARCOURS DES NOMS OCCUPES DU LEXIQUE (M1LX,0:M2LX)
         M1LX   = MCN( MNL )
C         M2LX   = MCN( MNL + 1 )
         NBENNM = MCN( MNL + 2 )
C
C        ADRESSE MCN DU PREMIER NOM OCCUPE
         I      = MCN( MNL + 5 )
C
C        LA BOUCLE SUR LES NOMS OCCUPES
C        ------------------------------
 30      IF( I .LE. 0 ) GOTO 20
         MN     = MNL + M1LX * I + NBENNM
C              LE NOM I EST OCCUPE.QUEL EST SON TYPE ?
               K4 = CHARX( MCN( MN + 1 ) )
               IF( K4 .EQ. 'TAMS' .OR. K4 .EQ. 'RELA' ) THEN
C                 FERMETURE DU TABLEAU MS
                  CALL TAMSFE( MCN( MN + 2 ) )
               ELSE IF( K4 .EQ. 'TATA' ) THEN
C                 FERMETURE DU TABLEAU DE TABLEAUX
C                 REGENERATION DE SON NOM A PARTIR DU LEXIQUE
                  CALL ENTNOM( NBENNM,MCN(MN-NBENNM),KNMLX )
                  CALL LXTTFE( NTL,KNMLX )
               ELSE IF( K4 .EQ. 'LEXI' ) THEN
C                 LE NOUVEAU LEXIQUE EST EMPILE
                  IF( LHPILE .GE. MXPILE ) THEN
C                    DEBORDEMENT DE LA PILE
                     NBLGRC(NRERR) = 1
                     KERR(1) = 'LXLXFE:PILE SATUREE'
                     CALL LEREUR
                     RETURN
                  ELSE
C                    LE LEXIQUE EST IL OUVERT ?
C                    SI NON NE PAS L'EMPILER
                     NOTAB = MCN( MN + 2 )
                     CALL TAMSRE( NOTAB , MGTAMS )
                     MCTAMS = MCG( MGTAMS + 4 )
                     IF( MCTAMS .GT. 0 ) THEN
C                       IL EST OUVERT
C                       LE NO DE TAMS DU LEXIQUE
                        MCN( MNPILE + LHPILE ) = NOTAB
                        LHPILE = LHPILE + 1
                     ENDIF
                  ENDIF
               ELSE
C                 AUCUN DES 4 TYPES
                  NBLGRC(NRERR) = 1
                  KERR(1) = KNOM
                  KERR(2) = 'DE TYPE INCONNU'
                  CALL LEREUR
                  RETURN
               ENDIF
C
C        LA FIN DE LA BOUCLE SUR LES NOMS OCCUPES
C        ----------------------------------------
         I = MCN( MN )
         GOTO 30
C
C     LA PILE EST VIDE ET DETRUITE
C     ============================
 9000 CALL TNMCDS( 'ENTIER' , MXPILE , MNPILE )
      END
