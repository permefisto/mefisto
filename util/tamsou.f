      SUBROUTINE TAMSOU( NOTAMS , MCTAMS )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    OUVRIR LE TABLEAU MS NOTAMS
C -----    C'EST A DIRE LE DECLARER EN MC
C                       LE LIRE SUR SON FICHIER ET LE TRANSFERER EN MC
C
C ENTREE :
C --------
C NOTAMS : NUMERO DU TABLEAU MS A OUVRIR
C
C SORTIE :
C --------
C MCTAMS : ADRESSE EN MC DU 1-ER MOT DU TABLEAU MS APRES LECTURE
C                  EN MCN POUR UN TABLEAU NON CARACTERE
C                  EN MCK POUR UN TABLEAU CARACTERE
C          0 SI NOTAMS INCORRECT (<=0 OU >M2TAMS)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : PERRONNET ALAIN ANALYSE NUMERIQUE PARIS  DECEMBRE 1983
C.......................................................................
      include"./incl/langue.inc"
      include"./incl/motmcg.inc"
      include"./incl/msvaau.inc"
      include"./incl/gsmenu.inc"
      include"./incl/pp.inc"
      include"./incl/ppmck.inc"
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      COMMON /        / MCN(MOTMCN)
      COMMON / MSSFTA / NOFISF,NBPASF,MOPASS,MGBUSF,NSFLIB,
     %                  M1FIMS,M2FIMS,MGFIMS,NSFIMS,LPFIMS,
     %                  M1TAMS,M2TAMS,MGBUTA,NBBUTA,NPTAMS,NATAMS,
     %                  NBCTMS,LLTAMS,LFTAMS,MGNPSF,NSFNPS,NPSNPS,
     %                  MGZLMG,MGZLMK,MGZLMN,MOTSMG,MOTSMK,MOTSMN,NTADAM
      COMMON / MSIMTA / NOIMPR
      CHARACTER*9       TYPNUM
C
C     LE NUMERO DU TAMS EST IL VALIDE ?
      IF( NOTAMS .LE. 0 .OR. NOTAMS .GT. M2TAMS ) THEN
C        NON
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10000) NOTAMS,M2TAMS
         ELSE
            WRITE(IMPRIM,20000) NOTAMS,M2TAMS
         ENDIF
10000 FORMAT(' ERREUR TAMSOU: NO DE TMS=',I7,' NON COMPRIS ENTRE 1 ET',
     % I7)
20000 FORMAT(' ERROR TAMSOU: TMS NUMBER=',I7,' NOT BETWEEN 1 and',
     % I7)
         MCTAMS = 0
         RETURN
      ENDIF
C
C     TEMOIN DE SATURATION DE LA MEMOIRE CENTRALE
      MCSATU = 0
C
C     RECHERCHE DU DESCRIPTEUR RETAMS(1..6,NOTAMS)
 10   CALL TAMSRE( NOTAMS , MGTAMS )
C     RETAMS(1)=NTTAMS NO DU TYPE DU TABLEAU MS
C                      0 SI LE TABLEAU N EST PAS DECLARE
C     RETAMS(2)=NVTAMS NOMBRE DE VARIABLES DU TABLEAU MS
C     RETAMS(3)=NFTAMS NUMERO DU FICHIER MS (1 A M2FIMS) DU TABLEAU MS
C                      0 SI NON DECLARE
C     RETAMS(4)=N1TAMS NUMERO DE LA PREMIERE PAGE DU TABLEAU SUR
C                      LE FICHIER NFTAMS QUI LE SAUVEGARDE
C     RETAMS(5)=MCTAMS ADRESSE MC DU PREMIER MOT DU TABLEAU MS
C                      >0 LE TABLEAU EST OUVERT( VERROUILLE ) EN MC
C                         IL EST PRET A ETRE UTILISE
C                         IL NE PEUT ETRE SAUVEGARDE SUR SON FICHIER
C                         TANT QU IL N EST PAS FERME
C                      =0 LE TABLEAU N EST PAS EN MC
C                         IL EST DECLARE MAIS PAS OUVERT
C                      <0 LE TABLEAU EST FERME,SUSCEPTIBLE
C                         D ETRE SAUVEGARDE SUR SON FICHIER MS
C                         IL EST UTILISABLE SEULEMENT S IL EST REOUVERT
C     RETAMS(6)=LSTAMS CHAINAGE AVAL DES TABLEAUX MS LIBRES OU FERMES
C
C     LE TABLEAU EST-IL EN MEMOIRE CENTRALE MC ?
      MCTAMS = MCG( MGTAMS + 4 )
      IF( MCTAMS .GT. 0 ) THEN
C
C         IL EXISTE EN MC.IL EST DEJA OUVERT
C         =============== ------------------
          IF( MCTAMS .LE. MOTMCN ) GOTO 9000
C
C         ADRESSE MCN ERRONEE
C         TENTATIVE DE TABLEAU SUPPOSE NON EN MC
          MCG( MGTAMS + 4 ) = 0
          MCTAMS = 0
          GOTO 10
C
      ELSE IF( MCTAMS .LT. 0 ) THEN
C
C        IL EXISTE EN MC.IL ETAIT FERME .IL EST REOUVERT.
C        =============== --------------------------------
         MCTAMS = - MCTAMS
         IF( MCTAMS .GT. MOTMCN ) THEN
C           ADRESSE MCN ERRONEE
C           TENTATIVE DE TABLEAU SUPPOSE NON EN MC
            MCG( MGTAMS + 4 ) = 0
            MCTAMS = 0
            GOTO 10
         ENDIF
C
         MCG( MGTAMS + 4 ) =   MCTAMS
C        CETTE ADRESSE EST SAUVEGARDEE AUTOMATIQUEMENT AVEC LE BUFFER
C        LE NUMERO DE CE TABLEAU EST OTE DU CHAINAGE DES TABLEAUX FERMES
         LSTAMS = MCG( MGTAMS + 5 )
         MCG( MGTAMS + 5 ) = 0
         NOTAF0 = 0
         NOTAF1 = LFTAMS
C
 15      IF( NOTAF1 .NE. NOTAMS ) THEN
C           CE N EST PAS LE TABLEAU QUI PRECEDE CELUI A SUPPRIMER
C           DES TABLEAUX FERMES
            CALL TAMSRE( NOTAF1 , MGTAMS )
            NOTAF0 = NOTAF1
            NOTAF1 = MCG( MGTAMS + 5 )
            GOTO 15
C
         ELSE
C           C EST LE TABLEAU FERME QUI PRECEDE CELUI A SUPPRIMER
C           EST-IL LE PREMIER ?
            IF( NOTAF0 .EQ. 0 ) THEN
C              OUI.LE SUIVANT FERME DU TABLEAU NOTAMS DEVIENT LE PREMIER
               LFTAMS = LSTAMS
C
            ELSE
C              NON.LE PRECEDENT EST RELIE AU SUIVANT
               MCG( MGTAMS + 5 ) = LSTAMS
            ENDIF
         ENDIF
C
         GOTO 9000
C
      ELSE
C
C        IL N'EXISTE PAS EN MC. QUEL EST SON TYPE ?
C        =====================
         NOTYTA = ABS( MCG( MGTAMS ) )
C        ES-CE UN TYPE CARACTERE ?
         IF( NOTYTA .EQ. NOTYPK ) THEN
C           OUI
            MGZL   = MGZLMK
            NOTYPT = 1
         ELSE
C           NON
            MGZL   = MGZLMN
            NOTYPT = 0
         ENDIF
C
C        LE RETOUR EST FORCE EN CAS DE PLACE INSUFFISANTE DANS LE
C        SUPER-TABLEAU MCK OU MCN
         MCTAMS = 1
         CALL TAMCDC( TYPNUM(NOTYTA),MCG(MGTAMS+1),MCG(MGZL),MCTAMS )
C
C        LA DECLARATION EST-ELLE CORRECTE ?
         IF( MCTAMS .LE. 0 .OR. MCTAMS .GT. MOTMCN ) THEN
C
C           NON.LE TABLEAU NOTAMS N'A PAS DE PLACE EN MC
C           --------------------------------------------
C           LE PREMIER TABLEAU FERME ET SUPERIEUR AU BESOIN EST RECHERCH
C
C           LE NOMBRE DE MOTS DU TABLEAU A OUVRIR
            MOTABL = MOTTAB( NOTYTA , MCG(MGTAMS+1) )
C
C           LE PREDECESSEUR DU 1-ER TABLEAU MC FERME
            NOTAF0 = 0
C
C           LE PREMIER TABLEAU MC FERME
            NOTAF1 = LFTAMS
C
C           CE TABLEAU EST-IL LE DERNIER DANS LE CHAINAGE DES TABLEAUX FERMES ?
 20         IF( NOTAF1 .EQ. 0 ) THEN
C              LA MC ETAIT-ELLE DEJA SATUREE ?
               IF( MCSATU .NE. 0 ) THEN
C
C                 OUI: LA MC EST INSUFFISANTE POUR CONTENIR LES TABLEAUX
                  NBLGRC( NRERR ) = 4
                  IF( LANGAG .EQ. 0 ) THEN
                     KERR(1)='TAMSOU: DEMANDE OUVERTURE DE TMS MAIS'
                     KERR(2)='MEMOIRE CENTRALE (MCN) TOTALEMENT SATUREE'
                     KERR(3)='ARRET IMPOSE APRES LE CLIC SOURIS'
                     KERR(4)='AUGMENTER MOTMCN DANS /incl/pp.inc'
                  ELSE
                     KERR(1)='TAMSOU: OPENING of a TMS but'
                     KERR(2)='MAIN MEMORY (MCN) FULL SATURATED'
                     KERR(3)='IMPOSED STOP AFTER the CLICK on MOUSE'
                     KERR(4)='AUGMENT MOTMCN in /incl/pp.inc'
                  ENDIF
                  CALL LEREUR
C
C                 SAISIE D'UN POINT PAR CLIC DE LA SOURIS OU ENTREE D'UN CARACTE
                  CALL SAIPTC( NOTYEV, NX, NY, NOCHAR )
                  STOP
CCCC
CCCC                 SAUVEGARDE IMPOSSIBLE ICI CAR ELLE DEMANDE TAMSOU
CCC                  CALL ARRET( 100 )
C
               ELSE
C
C                 NON: TOUS LES TABLEAUX FERMES SONT TRANSFERES EN MS
C                 ===================================================
                  MCSATU = 1
C
C                 LE PREMIER TABLEAU FERME
                  NOTAF1 = LFTAMS
C
C                 ES-CE LE DERNIER DES TABLEAUX FERMES ?
 30               IF( NOTAF1 .NE. 0 ) THEN
C
C                    RECHERCHE DU DESCRIPTEUR DE CE TABLEAU FERME
                     CALL TAMSRE( NOTAF1 , MGTAMS )
C
C                    TRANSFERT DE CE TABLEAU SUR SON FICHIER
                     CALL TAMSEC( MCG(MGTAMS) )
C
C                    LE TYPE DU TABLEAU A DETRUIRE EN MC
                     NOTYTA = ABS( MCG(MGTAMS) )
C
C                    RECHERCHE DE SON REPERTOIRE DE ZONES LIBRES
                     IF( NOTYTA .EQ. NOTYPK ) THEN
C                       CARACTERE
                        MGZL = MGZLMK
                     ELSE
C                       NUMERIQUE
                        MGZL = MGZLMN
                     ENDIF
C
C                    LA DESTRUCTION EN MCK OU MCN DE CE TABLEAU FERME
                     MNMNMN = ABS( MCG(MGTAMS+4) )
                     CALL TAMCDS( TYPNUM(NOTYTA),MCG(MGTAMS+1),
     %                            MCG(MGZL), MNMNMN )
                     MCG(MGTAMS+4) = MNMNMN
C
C                    PASSAGE AU TABLEAU FERME SUIVANT
                     NOTAF1 = MCG( MGTAMS + 5 )
C
C                    LE TABLEAU N'EST PLUS FERME.CHAINAGE SUIVANT NUL
                     MCG( MGTAMS + 5 ) = 0
                     GOTO 30
C
                  ELSE
C
C                    LE CHAINAGE DES TABLEAUX FERMES EST VIDE
                     LFTAMS = 0
                     GOTO 10
                  ENDIF
               ENDIF
C
            ELSE
C              NON:RECHERCHE DU DESCRIPTEUR RETAMS DE CE TABLEAU NOTAF1
               CALL TAMSRE( NOTAF1 , MGTAFE )
C
C              LE NOMBRE DE SES MOTS
               MOTAFE = MOTTAB( MCG(MGTAFE) , MCG(MGTAFE+1) )
C
C              SON TYPE :  1 SI CARACTERE 0 SINON
               NOTYP1 = MCG( MGTAFE )
               IF( NOTYP1 .EQ. NOTYPK ) THEN
                  NOTYP1 = 1
               ELSE
                  NOTYP1 = 0
               ENDIF
C
C              EST-IL SUPERIEUR ET DE MEME TYPE ?
               IF( MOTAFE .LT. MOTABL .OR. NOTYPT .NE. NOTYP1 ) THEN
C
C                 NON:RECHERCHE DU SUIVANT
                  NOTAF0 = NOTAF1
                  NOTAF1 = MCG( MGTAFE + 5 )
                  GOTO 20
C
               ELSE
C
C                 OUI:LE TABLEAU NOTAF1 EST SAUVEGARDE SUR SON FICHIER
                  CALL TAMSEC( MCG( MGTAFE ) )
C
C                 DESTRUCTION SELON LE TYPE CARACTERE OU NUMERIQUE
                  IF( NOTYTA .EQ. NOTYPK ) THEN
C                    CARACTERE
                     MGZL = MGZLMK
                  ELSE
C                    NUMERIQUE
                     MGZL = MGZLMN
                  ENDIF
C
C                 LE TABLEAU NOTAF1 EST DETRUIT EN MC
                  MCTAFE = ABS( MCG(MGTAFE+4) )
                  CALL TAMCDS( TYPNUM(MCG(MGTAFE)),MCG(MGTAFE+1),
     %                         MCG(MGZL),MCTAFE                   )
                  MCG(MGTAFE+4) = MCTAFE
C
C                 SON ADRESSE A ETE MISE A ZERO ET
C                 LE CHAINAGE SUR LE SUIVANT EST MIS A JOUR
                  NOTAF1 = MCG( MGTAFE + 5 )
C
C                 LE TABLEAU EST NI LIBRE NI FERME EN MC
                  MCG( MGTAFE + 5 ) = 0
C
C                 NOTAF1 ETAIT-IL LE PREMIER TABLEAU FERME EN MC ?
                  IF( NOTAF0 .EQ. 0 ) THEN
C                    OUI.IL DEVIENT LE PREMIER DES TABLEAUX LIBRES
                     LFTAMS = NOTAF1
                  ELSE
C                    NON:LE DESCRIPTEUR DU PREDECESSEUR EST LU
                     CALL TAMSRE( NOTAF0 , MGTAF0 )
                     MCG( MGTAF0 + 5 ) = NOTAF1
                  ENDIF
C
C                 REMONTEE AU DEBUT CAR LE BUFFER DU TABLEAU
C                 A PEUT ETRE ETE MODIFIE
                  GOTO 10
               ENDIF
            ENDIF
C
         ELSE
C
C           OUI.LE TABLEAU NOTAMS EST DECLARE A L'ADRESSE MCTAMS
C           ----------------------------------------------------
C
C           SON ADRESSE MC EST MISE A JOUR
            MCG( MGTAMS + 4 ) = MCTAMS
C           LE TABLEAU EST NI LIBRE NI FERME
            MCG( MGTAMS + 5 ) = 0
C           LE TABLEAU EST EVENTUELLEMENT TRANSFERE EN MC
            IF( MCG( MGTAMS ) .GT. 0 ) THEN
               CALL TAMSLE( MCG(MGTAMS) )
            ELSE
               MCG( MGTAMS ) = - MCG( MGTAMS )
            ENDIF
         ENDIF
      ENDIF
C
 9000 IF( NOIMPR .NE. 0 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,19000) NOTAMS,MCTAMS
         ELSE
            WRITE(IMPRIM,29000) NOTAMS,MCTAMS
         ENDIF
      ENDIF
C
19000 FORMAT(' OUVERTURE   TABLEAU MS NO',I9,' ADRESSE MC',I10)
29000 FORMAT(' OPEN of the ARRAY MS Number',I9,' ADRESS MC',I10)
C
      RETURN
      END
