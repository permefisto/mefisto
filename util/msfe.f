      SUBROUTINE MSFE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : FERMER LA MEMOIRE SECONDAIRE
C ----- SAUVEGARDER TOUS LES TABLEAUX MS OUVERTS ET FERMES
C       FERMER TOUS LES FICHIERS NUMERIQUES ET LE SUPER-FICHIER
C
C       POUR REUTILISER LA MS FAIRE UN CALL MSOU(...)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : PERRONNET ALAIN ANALYSE NUMERIQUE PARIS    DECEMBRE 1983
C....................................................................012
      include"./incl/ppmck.inc"
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      include"./incl/motmcg.inc"
      include"./incl/xyzext.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      include"./incl/msvaau.inc"
      COMMON / MSSFTA / NOFISF,NBPASF,MOPASF,MGBUSF,NSFLIB,
     %                  M1FIMS,M2FIMS,MGFIMS,NSFIMS,LPFIMS,
     %                  M1TAMS,M2TAMS,MGBUTA,NBBUTA,NPTAMS,NATAMS,
     %                  NBCTMS,LLTAMS,LFTAMS,MGNPSF,NSFNPS,NPSNPS,
     %                  MGZLMG,MGZLMK,MGZLMN,MOTSMG,MOTSMK,MOTSMN,NTADAM
      COMMON / MSIMTA / NOIMPR
      REAL              RMCN(1)
      INTEGER           PAGE1(29)
      EQUIVALENCE       (PAGE1(1),NOFISF)
      EQUIVALENCE       (MCN(1),RMCN(1))
      CHARACTER*16      ERREUR
      DATA              ERREUR / 'ERREUR MSFE   :' /
C
C     LA MS EST-ELLE OUVERTE ?
      IF( MGBUSF .LE. 0 ) RETURN
C
C     OUI.TOUS LES TABLEAUX MS SONT SAUVEGARDES SUR LEUR FICHIER
C     ==========================================================
C     SAUVEGARDE DU TABLEAU COOEXT
      CALL LXTSOU( NTADAM , 'MIN_MAX_SOMMETS' , NTCOEX , MNCOEX )
      IF( NTCOEX .GT. 0 ) THEN
C        SAUVEGARDE DU TABLEAU COOEXT ET DES 2+4 MOTS SUIVANTS
         CALL TRTATA( COOEXT , RMCN(MNCOEX) , 6*2+2+4 )
      ENDIF
C
C     SAUVEGARDE DES TABLEAUX MS ACTUELLEMENT EN MC
C     =============================================
      DO 20 I=1,M2TAMS
C
C        LE TABLEAU I EST-IL EN MC ?
         CALL TAMSRE( I , IG )
         IF( MCG(IG+4) .EQ. 0 ) GOTO 20
C
C        OUI: SON ADRESSE EST-ELLE CORRECTE ?
         MCT = ABS( MCG(IG+4) )
         IF( MCG(IG) .EQ. NOTYPK ) THEN
C           TABLEAU DE CARACTERES
            MOTS = MOTSMK
         ELSE
C           TABLEAU NUMERIQUE
            MOTS = MOTSMN
         ENDIF
         IF( MCT .LE. 0 .OR. MCT .GT. MOTS ) THEN
C
C           NON: IMPRESSION ERREUR
            WRITE(IMPRIM,10020) ERREUR,I,MOTS,(MCG(IG+N),N=0,5)
10020 FORMAT(A16,'ADRESSE DU TABLEAU MS',I12,' NULLE OU EXCEDE LE MAXIMU
     %M DE MOTS =',I12/' RETAMS=',8I12)
C
         ELSE
C
C           ADRESSE CORRECTE: LE TABLEAU MS EST SAUVEGARDE
            CALL TAMSEC( MCG(IG) )
C
         ENDIF
C
C        L ADRESSE MC DU TABLEAU EST MISE A ZERO
         MCG( IG + 4 ) = 0
C        LE TABLEAU N'EST PLUS FERME
         MCG( IG + 5 ) = 0
C
 20   CONTINUE
C
C     IL N Y A PLUS DE TABLEAU FERME
      LFTAMS = 0
C
C     SAUVEGARDE DES PAGES RETAMS ENCORE EN MCG
C     =========================================
      IF( MGBUTA .LE. 0 ) GOTO 35
      IG = MGBUTA
      DO 30 I=1,NBBUTA
         IF( MCG(IG) .GT. 0 ) THEN
C           LE NO DE LA PAGE RETAMS ACTUELLEMENT DANS LE BUFFER
C           DE 1 A NPTAMS
            N = MCG( IG + 1 )
C           LE NO DE LA PAGE SF DE CETTE MEME PAGE D APRES NPSFTA
            N = MCG( MGNPSF - 1 + N )
C           LA SAUVEGARDE SUR LE SF DE CETTE PAGE BUFFER
            CALL FIPAEC( NOFISF,N,MOPASF,MCG( MCG(IG) )  )
         ENDIF
         IG = IG + 3
 30   CONTINUE
C
C     DESTRUCTION MCG DU TABLEAU MCGBUTA
C     ==================================
      MGBUTA = 0
C
C     LE TABLEAU NPSFTA NUMERO SUPER-FICHIER DES PAGES DE RETAMS
C     EST DETRUIT EN MCG MAIS N'EST PAS SAUVEGARDE CAR IL N'A PAS CHANGE
C     =================================================================
 35   MGNPSF = 0
C
C     SAUVEGARDE DU BUFFER RED1FI (NO PAGES SUIVANTES) DE CHAQUE FICHIER
C     NUMERIQUE ET FERMETURE DE CHAQUE FICHIER
C     ==================================================================
      IF( MGFIMS .LE. 0 ) GOTO 45
      IG = MGFIMS
      DO 40 I=1,M2FIMS
C
C        LE FICHIER EST-IL DECLARE?
         N = MCG(IG)
         IF( N .GT. 0 ) THEN
C           OUI:LE BUFFER DE RED1FI EST SAUVEGARDE DANS SA PAGE SF
            IF( MCG(IG+4) .GT. 0 )
     %      CALL FIPAEC(NOFISF,MCG(IG+5)+MCG(IG+4)-1,MOPASF,
     %                         MCG(MCG(IG+3)) )
C           LA PAGE BUFFER EST DETRUITE EN MCG
            IF( MCG(IG+3) .GT. 0 ) MCG(IG+3) = 0
C           LE NO DE LA PAGE BUFFER EN MCG EST ANNULE
            MCG(IG+4) = 0
C           LE FICHIER EST FERME
            CLOSE( UNIT=MCG(IG) , ERR=38 , STATUS='KEEP' ,
     %             IOSTAT=IOERR )
         ENDIF
         GOTO 39
C
C        TRAITEMENT DE L ERREUR A LA FERMETURE DU FICHIER I
 38      WRITE(IMPRIM,10038) ERREUR,I,IOERR
10038 FORMAT(A16,' FICHIER MS',I6,' ERREUR',I12,' A LA FERMETURE'/)
         CALL IMFIMS
C
 39      IG = IG + M1FIMS
 40   CONTINUE
C
C     SAUVEGARDE DU TABLEAU REFIMS
C     ============================
C     PROTECTION DU NO DE LA PREMIERE PAGE SF DU TABLEAU REFIMS
      N = NSFIMS
      CALL TAMCEC( NOFISF,N,MOPASF,M1FIMS*M2FIMS,MCG(MGFIMS) )
      MGFIMS = 0
C
C     DESTRUCTION DU BUFFER DU SUPER-FICHIER
C     ======================================
 45   MGBUSF = 0
C
C     LE NOM DU PROJET N'A PU ETRE MODIFIE.IL N'EST PAS SAUVEGARDE
C
C     LE COMMON / MSSFTA / EST SAUVEGARDE DANS LA 1-RE PAGE DU SF
C     ===========================================================
      CALL FIPAEC( NOFISF, 1, 29, PAGE1 )
C
C     LE SUPER-FICHIER EST FERME
C     ==========================
      CLOSE( UNIT=NOFISF, ERR=9900, STATUS='KEEP', IOSTAT=N )
C
C     LES 3 REPERTOIRES MGZLM C-G-K SONT DETRUITS EN MCG
C     =================================================
      MOTSMG = MOTSMG + 3 * ( MCG(MGZLMN)+MCG(MGZLMG)+MCG(MGZLMK)+3 )
      MGZLMN = 0
      MGZLMG = 0
      MGZLMK = 0
      IF( NOIMPR .NE. 0 ) WRITE(IMPRIM,19000)
19000 FORMAT(' FIN DE LA FERMETURE DE LA MEMOIRE SECONDAIRE')
      RETURN
C
C     LE SUPER-FICHIER NE PEUT ETRE FERME => ERREUR
C     ---------------------------------------------
 9900 WRITE(IMPRIM,19900) ERREUR,NOFISF,N
19900 FORMAT(A16,'LE SUPER-FICHIER DE NO SYSTEME=',I5,' NE PEUT ETRE FER
     %ME.CODE ERREUR=',I12/)
      STOP
      END
