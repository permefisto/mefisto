      SUBROUTINE MSIM( OPTION , NOPTIO )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : IMPRIMER LA MEMOIRE SECONDAIRE        SELON LE PARAMETRE
C ----  OPTION  : 'TOUT    ','PAGE1   ','PROJET  ',
C                 'FICHIER ','ZONE_LIB','TABLEAUX'
C       OPTION  = 'TOUT    ' => IMPRESSION TOTALE DE LA MS
C       OPTION  = 'FICHIER ' + NOPTIO => IMPRESSION DU REPERTOIRE
C                                        DES FICHIERS ET DU CHAINAGE
C                                        DES PAGES DU FICHIER MS NOPTIO
C       OPTION  = 'ZONE_LIB' + NOPTIO = 1 => IMPRESSION REZLMG
C                                     = 2 => IMPRESSION REZLMK
C                                     = 3 => IMPRESSION REZLMN
C       OPTION  = 'TABLEAUX' + NOPTIO = NEGATIF OU SUPERIEUR M2TAMS
C                                       ALORS IMPRESSION DE
C                                       TOUS LES TABLEAUX MS . SINON
C                                     = IMPRESSION DU TABLEAU NOPTIO
C
C ATTENTION:LA MS DOIT AUPARAVANT ETRE OUVERTE
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : PERRONNET ALAIN ANALYSE NUMERIQUE PARIS  DECEMBRE 1983
C.......................................................................
      include"./incl/motmcg.inc"
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      COMMON / MSSFTA / NOFISF,NBPASF,MOPASF,MGBUSF,NSFLIB,
     %                  M1FIMS,M2FIMS,MGFIMS,NSFIMS,LPFIMS,
     %                  M1TAMS,M2TAMS,MGBUTA,NBBUTA,NPTAMS,NATAMS,
     %                  NBCTMS,LLTAMS,LFTAMS,MGNPSF,NSFNPS,NPSNPS,
     %                  MGZLM(3)            ,MOTSMG,MOTSMK,MOTSMN,NTADAM
      include"./incl/msvaau.inc"
      include"./incl/nmproj.inc"
      CHARACTER*48      NMBLAN
      CHARACTER*3       MM(3)
      CHARACTER*8       OPTIOS(6),OPTION
      CHARACTER*16      ERREUR
      CHARACTER*1       KMC(2)
      DATA              KMC    / 'K','N' /
      DATA              MM     / 'MCG','MCK','MCN' /
      DATA              NMBLAN / '                             ' /
      DATA              ERREUR / 'ERREUR MSIM   :' /
      DATA              OPTIOS/ 'TOUT    ','PAGE1   ','PROJET  ',
     %                           'FICHIER ','ZONE_LIB','TABLEAUX'/
C
C     LA MS EST-ELLE OUVERTE ?
C     ========================
      IF( MGBUSF .LE. 0 ) THEN
         WRITE(IMPRIM,10000) ERREUR
         GOTO 9900
      ENDIF
10000 FORMAT(A16,'LA MEMOIRE SECONDAIRE N EST PAS OUVERTE.PAS D IMPRESSI
     %ON'/)
      WRITE(IMPRIM,10005) OPTION,NOPTIO,NMORCI
10005 FORMAT(1H0,79(1H=)/' IMPRESSION DE TABLEAUX DE LA MS SELON L OPTIO
     %N ',A8,' NO',I4,' VERSION ',A8/1X,79(1H=))
C
C     RECHERCHE DE L'OPTION CHOISIE
C     =============================
      DO 10 NOP=1,6
         IF( OPTION(1:4) .EQ. OPTIOS(NOP)(1:4) ) GOTO 20
 10   CONTINUE
C
C     OPTION NON RETROUVEE
C     --------------------
      WRITE(IMPRIM,10010) ERREUR,OPTION,OPTIOS
10010 FORMAT(A16,'OPTION ',A8,' NON RETROUVEE PARMI'/8(1X,A8))
      GOTO 9900
C
C     TRAITEMENT DE L'OPTION
C     ======================
 20   GOTO( 100 , 100 , 200 , 300 , 400 , 500 ) , NOP
C
C     IMPRESSION DE LA PREMIERE PAGE= COMMON / MSSFTA / PAGE1(29)
 100  CALL IMMSP1
      IF( NOP .NE. 1 ) GOTO 9900
C
C     IMPRESSION DU NOM DU PROJET
C     ===========================
 200  WRITE(IMPRIM,10200) NMPROJ
10200 FORMAT('0NOM DU PROJET : ',A48/)
      IF( NOP .NE. 1 ) GOTO 9900
C
C     IMPRESSION DU REPERTOIRE DES FICHIERS
C                DU REPERTOIRE DES PAGES SUIVANTES DES FICHIERS
C     =========================================================
 300  CALL IMFIMS
      DO 350 I=1,M2FIMS
         IF( NOP .EQ. 1 ) GOTO 310
         IF( NOPTIO .NE. I ) GOTO 350
C
C        IMPRESSION DU REPERTOIRE DES PAGES SUIVANTES DU FICHIER I
 310     IG     = MGFIMS + M1FIMS * ( I - 1 )
         IF( MCG(IG) .LE. 0 ) GOTO 350
         MGBUFI = MCG( IG + 3 )
         IF( MGBUFI .LE. 0 ) GOTO 350
C
C        LE NOMBRE DE PAGES DU FICHIER I
         NBPAFI = MCG( IG + 1 )
C        LE NO DE LA PAGE DANS LE BUFFER
         NOPAG  = MCG( IG + 4 )
C        LE NOMBRE DE PAGES DU REPERTOIRE RED1FI DES PAGES SUIVANTES
         N  = ( MCG(IG+1) - 1 ) / MOPASF + 1
C        LE NUMERO DE LA 1-RE PAGE SF DE RED1FI
         NSFRED = MCG(IG+5) - 1
         IF( NOPAG .GT. 0 .AND. NOPAG .LE. N ) THEN
C           LA PAGE DANS LE BUFFER N EST PAS LA 1-RE.ELLE EST SAUVEGARDE
            CALL FIPAEC( NOFISF , NSFRED + NOPAG , MOPASF , MCG(MGBUFI))
         ENDIF
C
         WRITE(IMPRIM,10315) I
10315 FORMAT('0LE REPERTOIRE DES PAGES SUIVANTES DU FICHIER MS',I6)
C
C        LA BOUCLE SUR LES PAGES DE RED1FI
         NP     = 0
         DO 320 J=1,N
C
C           LECTURE DE LA PAGE J DU REPERTOIRE
            CALL FIPALE( NOFISF , NSFRED+J , MOPASF , MCG(MGBUFI) )
C
C           IMPRESSION
            IF( NP+MOPASF .LE. NBPAFI ) THEN
               LL = MOPASF
            ELSE
               LL = NBPAFI - NP
            ENDIF
            WRITE(IMPRIM,10320) (NP+L,MCG(MGBUFI-1+L),L=1,LL)
10320 FORMAT(3(' PAGE',I7,'=>',I7))
            NP = NP + MOPASF
 320     CONTINUE
C
C        MISE A JOUR DU NO DE PAGE DANS LE BUFFER
         MCG( IG + 4 ) = N
C
 350  CONTINUE
C
C     LES ZONES LIBRES
      IF( NOP .NE. 1 ) GOTO 9900
      J = 1
      L = 3
      GOTO 410
C
C     IMPRESSION DES ZONES LIBRES SELON NOPTION
C     NOPTION : 1 => MCG , 2 => MCK , 3 => MCN
C     =========================================
 400  IF( NOPTIO .LE. 0 .OR. NOPTIO .GT. 3 ) GOTO 9900
      IF( MGZLM(NOPTIO) .LE. 0 ) GOTO 9900
      J = NOPTIO
      L = NOPTIO
C
 410  DO 420 I=J,L
         CALL IMZLMC( MM(I) , MCG( MGZLM(I) )   )
 420  CONTINUE
      IF( NOP .NE. 1 ) GOTO 9900
C
C     IMPRESSION DU REPERTOIRE DES TABLEAUX DE LA MS
C     ==============================================
 500  IF( MGBUTA .LE. 0 .OR. MGNPSF .LE. 0 ) GOTO 9900
      WRITE(IMPRIM,10500)
10500 FORMAT('0LE REPERTOIRE DES TABLEAUX DE LA MS :'/1X,79(1H=))
C
C     IMPRESSION DU DESCRIPTEUR DES BUFFERS DE RETAMS
C     ===============================================
      CALL IMBUTA
C
C     IMPRESSION DU TABLEAU NPSFTA NO DES PAGES SF DE RETABD
C     ======================================================
      WRITE(IMPRIM,10510)
10510 FORMAT('0LE NUMERO DES PAGES SF DU REPERTOIRE DES TABLEAUX MS :'/
     %1X,79(1H=))
      WRITE(IMPRIM,10512) (I,MCG(MGNPSF-1+I),I=1,NPTAMS)
10512 FORMAT(4(' PAGE SF',I4,'=',I6))
      IF( NOP .EQ. 1 ) GOTO 510
      IF( NOPTIO .LE. 0 .OR. NOPTIO .GT. M2TAMS ) GOTO 510
C
C     IMPRESSION DU TABLEAU MS NOPTIO
C     ===============================
      J = NOPTIO
      L = NOPTIO
      GOTO 520
C
C     IMPRESSION DE TOUS LES TABLEAUX DE LA MS
C     ========================================
 510  J = 1
      L = M2TAMS
C
 520  DO 550 I=J,L
C        RECHERCHE DU DESCRIPTEUR DU TABLEAU MS I
         CALL TAMSRE( I , MGTAMS )
C        SI LE TABLEAU MS N EST PAS DECLARE SAUT AU SUIVANT
         IF( MCG(MGTAMS) .LE. 0 ) GOTO 550
C        IMPRESSION DE CE DESCRIPTEUR APRES DETERMINATION DE MCK OU MCN
         IF( MCG(MGTAMS) .EQ. NOTYPK ) THEN
            LL = 1
         ELSE
            LL = 2
         ENDIF
         WRITE(IMPRIM,10550) I,(MCG(MGTAMS+NP),NP=0,3)
     %                        ,KMC(LL),(MCG(MGTAMS+NP),NP=4,M1TAMS-1)
C        IMPRESSION DES VALEURS DU TABLEAU MS I
         CALL TAMSIM( I , NMBLAN )
10550 FORMAT(' TABLEAU MS',I8,': TYPE=',I4,' VARIABLES=',I8,
     %' FICHIER MS=',I3,' NO PAGE1=',I8,' ADRESSE MC',A1,'=',I8,
     %' SUIVANT=',I8)
 550  CONTINUE
C
 9900 RETURN
      END
