      SUBROUTINE IMPXYZNSEF
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : IMPORTER UN FICHIER de CARACTERES ASCII xyznsef.plsv.nomplsv
C------ contenant les XYZ des SOMMETS et TANGENTES (tms 'XYZSOMMET') et
C       les NUMEROS DES SOMMETS DU MAILLAGE (tms 'NSEF')
C       d'un Point ou Ligne ou Surface ou Volume
C       et CREER le Point ou Ligne ou Surface ou Volume et son MAILLAGE
C       dans le PROJET MEFISTO
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : Alain PERRONNET Saint PIERRE du PERRAY           Janvier 2021
C23456---------------------------------------------------------------012
      include"./incl/gsmenu.inc"
      include"./incl/langue.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/a_volume__definition.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/xyzext.inc"
      include"./incl/trvari.inc"
      include"./incl/pp.inc"
      COMMON         MCN(MOTMCN)
      REAL          RMCN(1)
      EQUIVALENCE   (MCN(1),RMCN(1))
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)

      LOGICAL        LEXIST, LOPEN
      CHARACTER*10   NMTYOB
      CHARACTER*2    K2C
      CHARACTER*24   KNMPLSV
      CHARACTER*40   KNMFIC
      CHARACTER*80   KLIGNE
      CHARACTER*1    KSUFIX(4), KCHAR
      DATA  KSUFIX / 'p', 'l', 's', 'v' /

C     LECTURE DU NOM DU FICHIER xyznsef.PLSV.NomPLSV
      GOTO 5

 1    IF( LANGAG .EQ. 0 ) THEN
         PRINT *,'NOM INCORRECT du FICHIER: ',KNMFIC
         PRINT *,'IL DOIT ETRE du TYPE xyznsef.plsv.nomplsv '
      ELSE
        PRINT *,'INCORRECT FILE NAME: ',KNMFIC
        PRINT *,'It MUST be of TYPE xyznsef.plsv.nomplsv '
      ENDIF
      GOTO 5

 3    IF( LANGAG .EQ. 0 ) THEN
         PRINT *,'INCORRECT CONTENU du FICHIER ',KNMFIC
      ELSE
        PRINT *,'INCORRECT CONTAIN of FILE ',KNMFIC
      ENDIF

 5    CALL INVITE( 170 )
      NCVALS = 0
      CALL LIRCAR( NCVALS, KNMFIC )
      IF( NCVALS .EQ. -1 ) GOTO 9999

C     TOUTES LES LETTRES SONT MISES EN MINUSCULES
      CALL MINUSC( KNMFIC )

      PRINT*
      IF( LANGAG .EQ. 0 ) THEN
         PRINT *,'Lecture du fichier: ',KNMFIC
      ELSE
         PRINT *,'Reading of file: ',KNMFIC
      ENDIF

C     EXTRACTION DU TYPE P L S V du NOM du FICHIER
      K = INDEX( KNMFIC, '.' )
      IF( K .LE. 0 ) GOTO 1

      K = K+1
      KCHAR = KNMFIC(K:K)
      DO NUTYOB=1,4
         IF( KCHAR .EQ. KSUFIX(NUTYOB) ) GOTO 10
      ENDDO
      GOTO 1

C     LE NOM du P L S V
 10   CALL SANSBL( KNMFIC, L )
      KNMPLSV = KNMFIC(K+2:L)

C     LE FICHIER KNMFIC EXISTE T IL?
      INQUIRE( FILE=KNMFIC, EXIST=LEXIST, OPENED=LOPEN )
      IF( .NOT. LEXIST ) THEN
C        LE FICHIER KNMFIC N'EXISTE PAS
         IF( LANGAG .EQ. 0 ) THEN
            PRINT*,'le FICHIER ',KNMFIC,' N''EXISTE PAS'
         ELSE
            PRINT*,'the FILE ',KNMFIC,' DO''NT EXIST'
         ENDIF
         GOTO 5
      ENDIF

C     OUVERTURE DU FICHIER KNMFIC
      CALL TRUNIT( NF )
      OPEN( UNIT=NF, ERR=1, STATUS='OLD',
     %      FILE=KNMFIC, ACCESS='SEQUENTIAL', FORM='FORMATTED' )

C     LECTURE DU NO et NOM du TYPE D'OBJET en PREMIERE LIGNE du FICHIER
C     -----------------------------------------------------------------
      READ( NF, '(A)' ) KLIGNE
      PRINT *, KLIGNE

ccc      IF( LANGAG .EQ. 0 ) THEN
ccc         READ(NF,10050) NUTYOB, NMTYOB(NUTYOB)
ccc10050    FORMAT('{ ',I2,'; ',A10,'; { Numero et Type du PLSV } }')
ccc      ELSE
ccc         READ(NF,10051) NUTYOB, NMTYOB(NUTYOB)
ccc10051    FORMAT('{ ',I2,'; ',A10,'; { Number and PLSV Type } }')
ccc      ENDIF

C     LECTURE DU NOM DU PLSV en SECONDE LIGNE du FICHIER
C     --------------------------------------------------
      READ( NF, '(A)' ) KLIGNE
      PRINT *, KLIGNE

ccc      IF( LANGAG .EQ. 0 ) THEN
ccc         WRITE(NF,10059) KNMPLS
ccc10059    FORMAT('{ ',A24,'; { NOM du PLSV } }')
ccc      ELSE
ccc         WRITE(NF,10060) KNMPLS
ccc10060    FORMAT('{ ',A24,'; { PLSV''s NAME } }')
ccc      ENDIF

C     LECTURE DU NOMBRE DE SOMMETS DU PLSV en TROISIEME LIGNE du FICHIER
C     ------------------------------------------------------------------
      READ( NF, '(I8)' ) NBSOM
      PRINT 10061,       NBSOM
10061 FORMAT( I8, ';   {NBSOM  NOMBRE de SOMMETS DU PLSV}' )
      IF( NBSOM .LE. 0 ) GOTO 3

      READ( NF, '(I8)' ) NBTGS
      PRINT 10062,       NBTGS
10062 FORMAT( I8, ';   {NBTGS  NOMBRE de TANGENTES DU PLSV}' )

C     LECTURE DE NBCOOR
      READ( NF, '(I8)' ) NBCOOR
      PRINT 10063,       NBCOOR
10063 FORMAT( I8, ';   {NBCOOR NOMBRE de COORDONNEES D UN SOMMET}')

C     RECHERCHE DU NOM DU P L S V DANS LE LEXIQUE DES P L S V
C     -------------------------------------------------------
      CALL MAJUSC( KNMPLSV )
      CALL LXLXOU( NTMN(NUTYOB), KNMPLSV, NTLXOB, MNLXOB )
C     S'IL EXISTE ALORS IL EST DETRUIT
      IF( NTLXOB .GT. 0 ) THEN
         CALL LXLXDS( NTMN(NUTYOB), KNMPLSV )
      ENDIF

C     GENERATION DU P L S V dans LE REPERTOIRE de son TYPE
C     ----------------------------------------------------
C     LE P L S V EST CREE ET OUVERT
      CALL LXLXDC( NTMN(NUTYOB), KNMPLSV, 24, 8 )
      CALL LXLXOU( NTMN(NUTYOB), KNMPLSV, NTLXOB, MNLXOB )
      IF( NTLXOB .LE. 0 ) GOTO 3

C     GENERATION DU TABLEAU 'XYZSOMMET' DU P L S V
C     --------------------------------------------
      CALL LXTSOU( NTLXOB, 'XYZSOMMET', NTXYZS, MNXYZS )
      IF( NTXYZS .GT. 0 ) THEN
         CALL LXTSDS( NTLXOB, 'XYZSOMMET' )
      ENDIF
      NBMOTS = WYZSOM + NBCOOR * (NBSOM+NBTGS)
      CALL LXTNDC( NTLXOB, 'XYZSOMMET','MOTS',  NBMOTS )
      CALL LXTSOU( NTLXOB, 'XYZSOMMET', NTXYZS, MNXYZS )

C     LE NOMBRE DE SOMMETS DU P L S V
      MCN( MNXYZS + WNBSOM ) = NBSOM
C     LE NOMBRE DE TANGENTES
      MCN( MNXYZS + WNBTGS ) = NBTGS
C     LE NOMBRE DE COORDONNEES D'UN SOMMET
      MCN( MNXYZS + WBCOOR ) = NBCOOR

C     LECTURE DES XYZ DES NBSOM SOMMETS
10070 FORMAT( 6(E15.7,A2) )
C     ADRESSE DE LA 1-ERE COORDONNEE DU 1-ER SOMMET DU TMS 'XYZSOMMET'
      MNS = MNXYZS + WYZSOM -1
      DO I=1,NBSOM
         READ(NF,10070) (RMCN(MNS+K),K2C,K=1,NBCOOR)
         MNS = MNS + NBCOOR
      ENDDO

C     LECTURE DE XYZ DES NBTGS TANGENTES
C     MNS ADRESSE DE LA 1-ERE COORDONNEE DE LA 1-ERE TG DU TMS 'XYZSOMMET'
      IF( NBTGS .GT. 0 ) THEN
         DO I=1,NBTGS
            READ(NF,10070) (RMCN(MNS+K),K2C,K=1,3)
            MNS = MNS + 3
         ENDDO
      ENDIF

C     AJOUT DU NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNXYZS + MOTVAR(6) ) = NONMTD( '~>>>XYZSOMMET' )

      IF( NUTYOB .EQ. 1 ) THEN
C        UN POINT N'A PAS DE TMS 'NSEF'
         NTNSEF = 0
         MNNSEF = 0
         GOTO 9900
      ENDIF


C     GENERATION DU TABLEAU 'NSEF' DU P L S V
C     ---------------------------------------
C     LECTURE DES VALEURS DE NSEF
      READ(NF, '(I8)') NUTYO
      PRINT 10091,     NUTYO
10091 FORMAT( I8, '; {NUTYOB  type du PLSV(1 2 3 4)}')
      IF( NUTYOB .NE. NUTYO ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            PRINT *,'INCOHERENCE du TYPE du PLSV NUTYOB=',NUTYOB,
     %              ' et NUTYOB LU=',NUTYO
         ELSE
            PRINT *,'INCORECT TYPE of PLSV NUTYOB=',NUTYOB,
     %              ' and reading NUTYOB=',NUTYO
         ENDIF
         GOTO 3
      ENDIF

      READ(NF, '(I8)') NUTFMA
      PRINT 10092,     NUTFMA
10092 FORMAT( I8, '; {NUTFMA  type FERME(1) ou NON(0) ou INCONNU(-1)}')

      READ(NF,'(I8)') NBSOEF
      PRINT 10093,    NBSOEF
10093 FORMAT( I8, '; {NBSOEF  nombre de sommets par EF}')

      READ(NF,'(I8)') NBTGEF
      PRINT 10094,    NBTGEF
10094 FORMAT( I8, '; {NBTGEF  nombre de tangentes par EF}')

      READ(NF,'(I8)') NBEFOB
      PRINT 10095,    NBEFOB
10095 FORMAT( I8, '; {NBEFOB  nombre des EF du PLSV}')

      READ(NF,'(I8)') NBEFTG
      PRINT 10096,    NBEFTG
10096 FORMAT( I8, '; {NBEFTG  nombre des EF avec TG}')

      READ(NF,'(I8)') NBEFAP
      PRINT 10097,    NBEFAP
10097 FORMAT( I8, '; {NBEFAP  nombre des EF avec POINTEUR sur EF a TG}')

      READ(NF,'(I8)') NUTYMA
      PRINT 10098,    NUTYMA
10098 FORMAT( I8, '; {NUTYMA  type structure(1,...,7) ou non(0)}')

      CALL LXTSOU( NTLXOB, 'NSEF', NTNSEF, MNNSEF )
      IF( NTNSEF .GT. 0 ) THEN
         CALL LXTSDS( NTLXOB, 'NSEF' )
      ENDIF
      NBMOTS = WUSOEF + NBSOEF*NBEFOB + NBEFAP + NBEFTG + NBTGEF*NBEFTG
      CALL LXTNDC( NTLXOB, 'NSEF','MOTS',  NBMOTS )
      CALL LXTSOU( NTLXOB, 'NSEF', NTNSEF, MNNSEF )

      MCN(MNNSEF+WUTYOB) = NUTYOB
      MCN(MNNSEF+WUTFMA) = NUTFMA
      MCN(MNNSEF+WBSOEF) = NBSOEF
      MCN(MNNSEF+WBTGEF) = NBTGEF
      MCN(MNNSEF+WBEFOB) = NBEFOB
      MCN(MNNSEF+WBEFTG) = NBEFTG
      MCN(MNNSEF+WBEFAP) = NBEFAP
      MCN(MNNSEF+WUTYMA) = NUTYMA

C     LA BOUCLE SUR LES EF DU MAILLAGE
      MNSEF = MNNSEF + WUSOEF -1
      DO 100 NEF = 1, NBEFOB
C        LECTURE DES NUMEROS DES SOMMETS DE L'EF NEF
         READ(NF,11100) (MCN(MNSEF+K),K2C,K=1,NBSOEF)
         MNSEF = MNSEF + NBSOEF
 100  ENDDO
11100 FORMAT(  8(I8,A2) )
11110 FORMAT( 10(I8,A2) )

C     ADRESSE DERRIERE LES NUMEROS DES SOMMETS DANS NSEF
      IF( NBTGS .GT. 0 ) THEN
         LDAPEF = WUSOEF + NBSOEF * NBEFOB
         MNSEF  = MNNSEF + LDAPEF - 1

         IF( NBEFAP .GT. 0 ) THEN
C           LECTURE DES numero>0 de l'EF a TG sinon 0
            READ(NF,11110) (MCN(MNSEF+K),K2C,K=1,NBEFAP)
         ENDIF

         LDNGEF = LDAPEF + NBEFAP
         MNSEF  = MNNSEF + LDNGEF - 1
         IF( NBEFTG .GT. 0 ) THEN

C           LECTURE du numero geometrique des EF a TG
C          ( 0 : 'C1 degre 3' ,
C            1 : 'CERCLE',  2 : 'ELLIPSE' ,    3 : 'COURBE B-SPLINE' ,
C           11 : 'SPHERE', 12 : 'ELLIPSOIDE', 13 : 'CYLINDRE', 14 : 'CONE',
C           15 : 'TORE'  , 16 : 'TRANSFINI' , 17 : 'SURFACE B-SPLINE') ;
            READ(NF,11110) (MCN(MNSEF+K),K2C,K=1,NBEFTG)

            IF( NBTGEF .GT. 0 ) THEN
               LDTGEF = LDNGEF + NBEFTG
               MNSEF  = MNNSEF + LDTGEF - 1
               DO I=1,NBEFTG
C                 LECTURE du +-no des TANGENTES de l'EF a TG
                  READ(NF,11110) (MCN(MNSEF+K),K2C,K=1,NBTGEF)
                  MNSEF = MNSEF + NBTGEF
               ENDDO
            ENDIF

         ENDIF

      ENDIF

C     AJOUT DU NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNNSEF + MOTVAR(6) ) = NONMTD( '~>>>NSEF' )

C     LE TMS DEFINITION DU P L S V A UN TMS DEFINITION
C     A PARTIR DES TMS XYZSOMMET et NSEF
C     10: 'DONNEE DIRECTE DES TMS XYZSOMMET ET NSEF'
C     ------------------------------------------------
 9900 CALL LXTNDC( NTLXOB, 'DEFINITION', 'MOTS', WUTSSV+1 )
      CALL LXTSOU( NTLXOB, 'DEFINITION', NTDEFI, MNDEFI )
C     TRANSFORMATION (I pour IDENTITE)
      MCN( MNDEFI + WTYTRV ) = 1

C     TYPE DU P L S V
C     10: 'DONNEE DIRECTE DES TMS XYZSOMMET ET NSEF'
      MCN( MNDEFI + WUTYVO ) = 10

C     NUTSOV 'tableau MS xyzsommet' tms ~>>>XYZSOMMET
      MCN( MNDEFI + WUTSOV ) = NTXYZS

C     NUTSSV 'tableau MS no des sommets des EF' tms ~>>>NSEF
      MCN( MNDEFI + WUTSSV ) = NTNSEF

C     AJOUT DE LA DATE du TMS 'DEFINITION'
      CALL ECDATE( MCN(MNDEFI) )

C     AJOUT DE LA DATE du TMS 'XYZSOMMET' POSTERIEURE a DEFINITION
      CALL ECDATE( MCN(MNXYZS) )

C     MISE A JOUR DE LA DIMENSION DES COORDONNEES
      CALL DIMCOO( NBSOM, MCN(MNXYZS+WYZSOM), NDIMLI )

C     MISE A JOUR DU CADRE MINMAX DES NBCOOR COORDONNEES DES SOMMETS
      CALL MAJEXT( MNXYZS )

C     FERMETURE DU FICHIER
      CLOSE( UNIT=NF )

      IF( NUTYOB .GT. 1 ) THEN
C        AJOUT DE LA DATE du TMS 'NSEF' POSTERIEURE a DEFINITION
         CALL ECDATE( MCN(MNNSEF) )
      ENDIF

      IF( LANGAG .EQ. 0 ) THEN
         PRINT *,'Creation ',NMTYOB(NUTYOB),' ',KNMPLSV,' CORRECTE'
      ELSE
         PRINT *,'CORRECT CREATION of ',NMTYOB(NUTYOB),' ',KNMPLSV
      ENDIF
      PRINT*

C     TRACE du MAILLAGE du P L S V
      IF( INTERA .GE. 1 .AND. NUTYOB .GT. 0 .AND. MNXYZS .GT. 0 ) THEN
         INIEXT = 0
         NOTYVI = 0
         CALL MAJEXT( MNXYZS )
         CALL NUOBNM( NMTYOB(NUTYOB), KNMPLSV, NUPLSV )
         CALL T1SOBJ( NUTYOB, KNMPLSV, NUPLSV,
     %                NTNSEF, MNNSEF,  NTXYZS, MNXYZS )
      ENDIF

C     AFFICHAGE DE LA QUALITE DU MAILLAGE SI SURFACE ou VOLUME
      CALL IMPQUA( NUTYOB, KNMPLSV, MNNSEF, MNXYZS,
     %             NBEFMQ, QUAMIN, SURVOLEF ) 

 9999 RETURN
      END
