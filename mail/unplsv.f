      SUBROUTINE UNPLSV( NTYPEU, NUTYOB, NUOBUN, NBOBUN, NUMOBJ,
     %                   NTTSMA, MNTSMA, NTSOMM, MNSOMM,
     %                   NTUNIO, MNUNIO, IERR   )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  GENERER LES TABLEAUX 'XYZSOMMET' 'NSEF' ET 'UNION'
C -----  DU PLSV UNION DE NBOBUN PLSV APRES IDENTIFICATION DES
C        SOMMETS COMMUNS
C        PLSV DESIGNE UN POINT OU UNE LIGNE OU UNE SURFACE OU UN VOLUME
C        L'UNION SE FAIT ICI ENTRE PLSV DE MEME TYPE
C        UNE LIGNE UNION EST UNE UNION DE PLUSIEURS LIGNES
C        MAIS NE PEUT ETRE UNE UNION DE LIGNES ET DE POINTS ...
C
C        A PRIORI SI CONFORMITE LE RESULTAT EST C0-CONTINU
C        PUISQUE LES SOMMETS COMMUNS AUX 2 MAILLAGES SONT IDENTIFIES
C
C        SI SURFACE ET TYPE 52 ALORS LES TANGENTES AUX POINTS IDENTIFIES
C        COMMUNS SONT PROJETEES SUR UN MEME PLAN A DISTANCE MINIMALE
C        CE QUI ASSURE UNE GEOMETRIE SURFACIQUE G1-CONTINUE

C ENTREES :
C ---------
C NTYPEU : NUMERO DE TYPE DE L'UNION
C          51 => STANDARD C0-CONTINU
C          52 => SURFACE RENDUES G1-CONTINUES PAR TGS PROJETEES SUR UN MEME PLAN
C NUTYOB : NUMERO DU TYPE DES PLSV ( IL EST LE MEME POUR TOUS )
C          1:POINT 2:LIGNE 3:SURFACE 4:VOLUME 5:PLSV
C NUOBUN : NUMERO DU PLSV UNION DANS SON LEXIQUE

C ENTREES ET SORTIES:
C -------------------
C NBOBUN : NOMBRE DE  PLSV DE L'UNION
C NUMOBJ : NUMERO DES PLSV DE L'UNION
C          CES 2 PARAMETRES SONT MODIFIES SI UN NUMERO EST DOUBLE

C SORTIES :
C ---------
C NTTSMA : NUMERO      DU TMS 'NSEF' DU PLSV UNION
C MNTSMA : ADRESSE MCN DU TMS 'NSEF'
C NTSOMM : NUMERO      DU TMS 'XYZSOMMET'
C MNSOMM : ADRESSE MCN DU TMS 'XYZSOMMET'
C NTUNIO : NUMERO      DU TMS 'UNION'
C MNUNIO : ADRESSE MCN DU TMS 'UNION' DU PLSV UNION
C IERR   : 0 SI PAS D'ERREUR, >0 NOMBRE DES ERREURS RENCONTREES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     OCTOBRE 1988
C MODIFS : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 1995
C MODIFS : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS        JUIN 1998
C MODIFS : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS         MAI 1999
C MODIFS : ALAIN PERRONNET LJLL UPMC et St Pierre du Perray    Juin 2008
C23456---------------------------------------------------------------012
      PARAMETER        (NBIPAX=23,NBIPX1=24,
     %                  NBIPAY=23,NBIPY1=24,
     %                  NBIPAZ=23,NBIPZ1=24)
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/a___nsef.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___materiaux.inc"
      include"./incl/a___union.inc"
      include"./incl/trvari.inc"
      COMMON / TRTETR / STOPTE,TRACTE
      LOGICAL           STOPTE,TRACTE,TRACTE0
      COMMON / EPSSSS / EPZERO,EPSXYZ
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL             RMCN(1)
      EQUIVALENCE     (RMCN(1),MCN(1))
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)

      INTEGER           NUMOBJ(1:NBOBUN), NOSOEL(1:64)
      REAL              COIN(6,2)
      CHARACTER*10      NMTYOB,NAMTYOB
      CHARACTER*24      KNOM
      REAL              XX(8),YY(8),ZZ(8),XYZ(3,8)

C     INITIALISATIONS
      TRACTE0 = TRACTE
      IERR   = 0
      MNNSEF = 0
      MNXYOB = 0
      MNDESO = 0
      MNPAVE = 0
      MNCHAI = 0
      MNNBSU = 0

C     MAXIMUM DES EPSILON
      EPS    = MAX( EPZERO, EPSXYZ )
      UNMEPS = 1.0 - EPS

C     LE NUMERO DE TMS DU LEXIQUE DU TYPE DE PLSV
      IF( NUTYOB .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:4),'(I4)') NUTYOB
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'unplsv: TYPE PLSV INCORRECT '//KERR(MXLGER)(1:4)
         ELSE
            KERR(1) = 'unplsv: INCORRECT PLSV TYPE '//KERR(MXLGER)(1:4)
         ENDIF
         CALL LEREUR
         RETURN
      ENDIF
      CALL TAMSOU( NTMN(NUTYOB), MNLX )

C     POUR ACCELERER L'IDENTIFICATION UN PAVAGE DE L'ESPACE EST FAIT
C     RECHERCHE DU MIN_MAX DES SOMMETS DES NBOBUN PLSV
      R = RINFO( 'GRAND' )
      XMINS = R
      YMINS = R
      ZMINS = R
      XMAXS =-R
      YMAXS =-R
      ZMAXS =-R

C     LES TABLEAUX DE SAUVEGARDE DES ADRESSES MCN DES PLSV
      NBOBU0 = NBOBUN
      CALL TNMCDC( 'ENTIER', NBOBU0, MNNSEF )
      CALL TNMCDC( 'ENTIER', NBOBU0, MNXYOB )

C     LE TABLEAU POINTEUR SUR LE DERNIER SOMMET DE CHAQUE PLSV
      CALL TNMCDC( 'ENTIER', 1 + NBOBU0, MNDESO )
      MCN( MNDESO ) = 0

C     OUVERTURE DES TABLEAUX 'NSEF' ET 'XYZSOMMET' DE CHAQUE PLSV
      NBEFOB = 0
      NBSOEF = 0
      NBSOMM = 0
      NBQUST = 0
      NBTGUN = 0
      NBTGEU = 0
      NBEFTG = 0
      NBCOOR = 0

      DO I=1,NBOBU0

C        NOMBRE D'ERREURS POUR LE PLSV I
         IER = 0

C        LE NOM DU PLSV I DE L'UNION
         CALL NMOBNU( NMTYOB(NUTYOB), NUMOBJ(I), KNOM )

C        L'ADRESSE MCN DU NUMERO DU TMS LEXIQUE DU PLSV
         MN = MNLX + MCN( MNLX ) * NUMOBJ(I) + MCN( MNLX+2 ) + 2

C        LE NUMERO DE TMS DU LEXIQUE DU PLSV I DE L'UNION
         NTLXOB = MCN( MN )
         IF( NTLXOB .LE. 0 ) THEN
            NBLGRC(NRERR) = 2
            KERR(1) = 'unplsv: PLSV  ' // KNOM
            WRITE(KERR(MXLGER)(1:4),'(I4)') I
            IF( LANGAG .EQ. 0 ) THEN
               KERR(2) = 'unplsv: PLSV '//KERR(MXLGER)(1:4)//' INCONNU'
            ELSE
               KERR(2) = 'unplsv: UNKNOWN PLSV '//KERR(MXLGER)(1:4)
            ENDIF
            CALL LEREUR
            IER = IER + 1
            GOTO 3
         ENDIF

C        VERIFICATION DE NOMS NON REPETES DANS LA LISTE
         DO J=1,I-1
            IF( NUMOBJ(J) .EQ. NUMOBJ(I) ) THEN
C              NOM RETROUVE 2 FOIS
               NBLGRC(NRERR) = 3
               IF( LANGAG .EQ. 0 ) THEN
                  KERR(1) = 'UNION AVEC 2 MEMES NOMS PARMI LES PLSV'
                  KERR(2) = 'NOM DOUBLE: ' // KNOM
                  KERR(3) = 'LE SECOND NOM EST SUPPRIME'
               ELSE
                  KERR(1) = 'UNION WITH 2 SAME NAMES AMONG the PLSV'
                  KERR(2) = 'DOUBLE NAME: ' // KNOM
                  KERR(3) = 'The SECOND NAME is DELETED'
               ENDIF
               CALL LEREUR
C              LE NUMERO DU PLSV EST ANNULE
               NUMOBJ(I) = 0
C              LE NOUVEAU NOMBRE DE PLSV APRES IDENTIFICATION
               NBOBUN = NBOBUN - 1
               GOTO 3
            ENDIF
         ENDDO

C        LE TMS 'XYZSOMMET' DU PLSV I
         CALL LXTSOU( NTLXOB, 'XYZSOMMET', NT2, MNS )
         IF( MNS .LE. 0 ) THEN
            NBLGRC(NRERR) = 1
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) ='unplsv: PAS de TMS ''XYZSOMMET'' pour ' // KNOM
            ELSE
             KERR(1) ='unplsv: ' // KNOM // ' ''XYZSOMMET'' TMS UNKNOWN'
            ENDIF
            CALL LEREUR
            IER = IER + 1
            GOTO 3
         ENDIF

C        LE TMS 'XYZSOMMET' EXISTE
         MCN(MNXYOB-1+I) = MNS
         NBC = MCN(MNS+WBCOOR)
         IF( NBCOOR .EQ. 0 ) THEN
            NBCOOR = NBC
         ELSE IF( NBC .NE. NBCOOR ) THEN
            NBLGRC(NRERR) = 2
            WRITE(KERR(MXLGER)(1:4),'(I4)') NBC
            WRITE(KERR(MXLGER)(5:8),'(I4)') NBCOOR
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) ='unplsv: NB COORDONNEES DIFFERENTS'
     %       //KERR(MXLGER)(1:4)//' et '//KERR(MXLGER)(5:8)
               KERR(2) = 'A PARTIR DE ' // KNOM
            ELSE
               KERR(1) ='unplsv: DIFFERENT COORDINATE NUMBER'
     %       //KERR(MXLGER)(1:4)//' and '//KERR(MXLGER)(5:8)
               KERR(2) = 'FROM ' // KNOM
            ENDIF
            CALL LEREUR
            IER = IER + 1
            GOTO 3
         ENDIF

C        LE TMS 'NSEF' DU PLSV I
         IF( NUTYOB .GT. 1 ) THEN
            CALL LXTSOU( NTLXOB, 'NSEF', NT1, MNMA )
            IF( MNMA .LE. 0 ) THEN
               NBLGRC(NRERR) = 1
               IF( LANGAG .EQ. 0 ) THEN
                  KERR(1) ='unplsv: PAS de TMS ''NSEF'' pour ' // KNOM
               ELSE
                  KERR(1) ='unplsv: ' // KNOM //' ''NSEF'' TMS UNKNOWN'
               ENDIF
               CALL LEREUR
               IER = IER + 1
               GOTO 3
            ENDIF

C           LE TMS 'NSEF' EXISTE
            MCN(MNNSEF-1+I) = MNMA
            IF(MCN( MNMA + WUTYMA ).EQ.4) NBQUST = NBQUST+1
         ELSE
            NT1 = NT2
         ENDIF

C        SI TRAITEMENT DE POINTS,IL N'EXISTE PAS DE TABLEAU 'NSEF'
         IF( NT1 .LE. 0 .OR. NT2 .LE. 0 ) THEN
            NBLGRC(NRERR) = 2
            WRITE(KERR(MXLGER)(1:4),'(I4)') I
            IF( LANGAG .EQ. 0 ) THEN
              KERR(1) ='unplsv: PLSV '//KERR(MXLGER)(1:4)//' NON MAILLE'
              KERR(2) = KNOM // ' INCONNU'
            ELSE
              KERR(1) ='unplsv: PLSV '//KERR(MXLGER)(1:4)//' NOT MESHED'
               KERR(2) = KNOM // ' UNKNOWN'
            ENDIF
            CALL LEREUR
            IER = IER + 1
            GOTO 3
         ENDIF

C        LE NOMBRE DE SOMMETS DU PLSV I ET DE L'UNION
         NX = MCN( MNS + WNBSOM )
         IF( NX .LE. 0 ) THEN
            NBLGRC(NRERR) = 2
            WRITE(KERR(MXLGER)(1:12),'(I12)') NX
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) ='unplsv: PLSV ' // KNOM
               KERR(2) = 'A' // KERR(MXLGER)(1:12) // ' SOMMETS!'
            ELSE
               KERR(1) ='unplsv: PLSV ' // KNOM
               KERR(2) = 'HAS ' // KERR(MXLGER)(1:12) // ' VERTICES!'
            ENDIF
            CALL LEREUR
            IER = IER + 1
            GOTO 3
         ENDIF
         NBSOMM = NBSOMM + NX

C        LE DERNIER SOMMET DU PLSV I APRES SOMMATION
         MCN( MNDESO + I ) = NBSOMM

C        LE NOMBRE DE TANGENTES DU PLSV I ET DE L'UNION
         NBTGUN = NBTGUN + MCN( MNS + WNBTGS )

C        LE NOMBRE D'EF A TG DE L'UNION
         NBEFTG = NBEFTG + MCN( MNMA + WBEFTG )

C        LE CADRE EXTREME DES COORDONNEES
         CALL CADEXT( MNS, COIN )
         XMINS = MIN( XMINS, COIN(1,1) )
         YMINS = MIN( YMINS, COIN(2,1) )
         ZMINS = MIN( ZMINS, COIN(3,1) )
         XMAXS = MAX( XMAXS, COIN(1,2) )
         YMAXS = MAX( YMAXS, COIN(2,2) )
         ZMAXS = MAX( ZMAXS, COIN(3,2) )

C        LES PARAMETRES DES NO SOMMET DU MAILLAGE DU PLSV I
         IF( NUTYOB .GT. 1 ) THEN

C           UNION DE LSV
            CALL NSEFPA( MCN(MNMA),
     %                   NUTYMA, NBSOEL, NBSOS, NBTGEF,
     %                   LDAPEF, LDNGEF, LDTGEF, NBSSO,
     %                   NX  , NY  , NZ  ,
     %                   IERRR  )
            IER = IER + IERRR

         ELSE

C           UNION DE POINTS
            NBSOEL = 1
            NBSOS  = 1
            NBSSO  = MCN( MNS + WNBSOM )
            NBTGEF = 0
            IF( MCN(MNS+WNBTGS) .GT. 0 ) NBTGEF=1

         ENDIF

         IF( IER .NE. 0 ) THEN
            NBLGRC(NRERR) = 1
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = KNOM // ' a UN MAILLAGE INCORRECT'
            ELSE
               KERR(1) = KNOM // ' has an INCORRECT MESH'
            ENDIF
            CALL LEREUR
            GOTO 3
         ENDIF

         IF( I .GT. 1 .AND. NBSOS .NE. NBSOEF ) THEN
            NBLGRC(NRERR) = 2
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = KNOM // ' MAILLAGE INCORRECT'
               KERR(2) = 'NOMBRES DE SOMMETS PAR EF DIFFERENTS'
            ELSE
               KERR(1) =  KNOM // ' HAS an INCORRECT MESH'
               KERR(2) = 'DIFFERENT VERTEX NUMBERS by FE'
            ENDIF
            CALL LEREUR
            IER = IER + 1
            GOTO 3
         ENDIF
         NBSOEF = MAX( NBSOEF, NBSOS )

         IF( NBTGEU.GT.0 .AND. NBTGEF.GT.0 .AND. NBTGEU.NE.NBTGEF ) THEN
            NBLGRC(NRERR) = 2
            IF( LANGAG . EQ. 0 ) THEN
               KERR(1) = KNOM // ' MAILLAGE INCORRECT'
               KERR(2) = 'NOMBRES DE TANGENTES PAR EF DIFFERENTS'
            ELSE
               KERR(1) =  KNOM // ' Has an INCORRECT MESH'
               KERR(2) = 'DIFFERENT TANGENT NUMBERS by FE'
            ENDIF
            CALL LEREUR
            IER = IER + 1
            GOTO 3
         ENDIF
         NBTGEU = MAX( NBTGEU, NBTGEF )

C        LE NOMBRE TOTAL D'EF DU MAILLAGE DU PLSV UNION
         NBEFOB = NBEFOB + NBSSO

C        AFFICHAGE DES PLSV DE L'UNION
         NAMTYOB = NMTYOB(NUTYOB)
         IF( INDEX(NAMTYOB,'LIGNE') .GT. 0 ) NAMTYOB='LINE'
         IF( LANGAG . EQ. 0 ) THEN
            WRITE(IMPRIM,10003) NAMTYOB,KNOM,
     %                          MCN(MNS+WNBSOM),MCN(MNS+WNBTGS),
     %                          NBSSO,MCN(MNMA+WBEFTG)
         ELSE
            WRITE(IMPRIM,20003) NAMTYOB,KNOM,
     %                          MCN(MNS+WNBSOM),MCN(MNS+WNBTGS),
     %                          NBSSO,MCN(MNMA+WBEFTG)
         ENDIF

C        CUMUL DES ERREURS DE L'UNION
 3       IERR = IERR + IER

      ENDDO

      IF( IERR .NE. 0 ) GOTO 9000

10003 FORMAT(' UNION ',A,' ',A,' avec',I9,' ST',I9,' TG',
     %I9,' EF' , I9, ' EF a TG')
20003 FORMAT(' UNION ',A,' ',A,' with',I9,' VT',I9,' TG',
     %I9,' FE' , I9, ' FE with TG')

C     SUPPRESSION DES NUMEROS DE PLSV NULS CAR DOUBLES
      J = 0
      DO I=1,NBOBU0
         IF( NUMOBJ(I) .EQ. 0 ) THEN
C           NOM DOUBLE A SUPPRIMER
            J = J + 1
         ELSE
C           NOM SIMPLE A DECALER
            NUMOBJ(I-J) = NUMOBJ(I)
C           LES  ADRESSES DES TABLEAUX
            MCN(MNXYOB-1+I-J) = MCN(MNXYOB-1+I)
            IF(NUTYOB.GT.1) MCN(MNNSEF-1+I-J) = MCN(MNNSEF-1+I)
            MCN(MNDESO  +I-J) = MCN(MNDESO  +I)
         ENDIF
      ENDDO

      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10004) NBSOMM,NBTGUN
      ELSE
         WRITE(IMPRIM,20004) NBSOMM,NBTGUN
      ENDIF
10004 FORMAT(' unplsv: NOMBRE SOMMETS   AVANT IDENTIFICATION=',I9/
     %       ' unplsv: NOMBRE TANGENTES AVANT IDENTIFICATION=',I9)
20004 FORMAT(' unplsv: BEFORE IDENTIFICATION VERTEX  NUMBER=',I9/
     %       ' unplsv: BEFORE IDENTIFICATION TANGENT NUMBER=',I9)

      IF( NBCOOR .NE. 3 ) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:4),'(I4)') NBCOOR
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) ='unplsv: NB COORDONNEES='
     %      //KERR(MXLGER)(1:4)//' NON TRAITE'
         ELSE
            KERR(1) ='unplsv: COORDINATE NUMBER='
     %      //KERR(MXLGER)(1:4)//' NOT PROGRAMMED'
         ENDIF
         CALL LEREUR
         IERR = IERR + 1
         GOTO 9000
      ENDIF

C     INITIALISATION DU PAVAGE
C     ========================
C     DECLARATION ET OUVERTURE EN MC DU TABLEAU 'PAVES'
      MOPAVE = (NBIPX1+1) * (NBIPY1+1) * (NBIPZ1+1)
      CALL TNMCDC( 'ENTIER', MOPAVE, MNPAVE )
C     MISE A ZERO DU CHAINAGE DES PAVES
      CALL AZEROI( MOPAVE, MCN(MNPAVE) )

C     LE TABLEAU DES CHAINAGES DES SOMMETS DANS LES PAVES
      CALL TNMCDC( 'ENTIER', NBSOMM, MNCHAI )
      CALL AZEROI( NBSOMM, MCN(MNCHAI) )
      MNCHA1 = MNCHAI - 1

C     CALCUL DES 'ECHELLES' DANS LES 3 DIRECTIONS
      ECHPAX = XMAXS - XMINS
      IF( ECHPAX .LE. EPS ) THEN
         XMAXS  = XMAXS + 1.
         ECHPAX = 1.
      ENDIF
      XMINS  = XMINS - EPS * ECHPAX
      XMAXS  = XMAXS + EPS * ECHPAX
      ECHPAX = XMAXS - XMINS
      ECHPAX = NBIPX1 / ECHPAX

      ECHPAY = YMAXS - YMINS
      IF( ECHPAY .LE. EPS ) THEN
         YMAXS  = YMAXS + 1.
         ECHPAY = 1.
      ENDIF
      YMINS  = YMINS - EPS * ECHPAY
      YMAXS  = YMAXS + EPS * ECHPAY
      ECHPAY = YMAXS - YMINS
      ECHPAY = NBIPY1 / ECHPAY

      ECHPAZ = ZMAXS - ZMINS
      IF( ECHPAZ .LE. EPS ) THEN
         ZMAXS  = ZMAXS + 1.
         ECHPAZ = 1.
      ENDIF
      ZMINS  = ZMINS - EPS * ECHPAZ
      ZMAXS  = ZMAXS + EPS * ECHPAZ
      ECHPAZ = ZMAXS - ZMINS
      ECHPAZ = NBIPZ1 / ECHPAZ

C     DECLARATION OUVERTURE DU TABLEAU 'XYZSOMMET' DU PLSV UNION
      MN = MNLX + MCN( MNLX ) * NUOBUN + MCN( MNLX+2 ) + 2
      NTLXOB = MCN( MN )
      CALL LXTSOU( NTLXOB, 'XYZSOMMET', NTSOMM, MNSOMM )
      IF( NTSOMM .GT. 0 ) THEN
C        DESTRUCTION DU TABLEAU
         CALL LXTSDS( NTLXOB, 'XYZSOMMET' )
      ENDIF
      MOT = WYZSOM + ( NBSOMM + NBTGUN ) * 3
      CALL LXTNDC( NTLXOB, 'XYZSOMMET', 'MOTS', MOT )
      CALL LXTSOU( NTLXOB, 'XYZSOMMET', NTSOMM, MNSOMM )

C     DECLARATION OUVERTURE DU TABLEAU 'UNION' DU PLSV UNION
      CALL LXTSOU( NTLXOB, 'UNION', NTUNIO, MNUNIO )
      IF( NTUNIO .GT. 0 ) THEN
C        DESTRUCTION DU TABLEAU
         CALL LXTSDS( NTLXOB, 'UNION' )
      ENDIF
      MOT = WDSCOU + 1 + NBOBUN + NBSOMM + 1 + NBOBUN + NBTGUN
      CALL LXTNDC( NTLXOB, 'UNION', 'MOTS', MOT )
      CALL LXTSOU( NTLXOB, 'UNION', NTUNIO, MNUNIO )
      MNDETG = MNUNIO + WDSCOU + 1 + NBOBUN + NBSOMM
      MNNUTG = MNDETG + 1 + NBOBUN

C     INITIALISATIONS
C     LE NOMBRE DE SOMMETS DE L'UNION APRES IDENTIFICATION
      NBSOM = 0
C     LE NOMBRE DE TANGENTES DE L'UNION APRES IDENTIFICATION
      NUTG = 0
      MCN(MNDETG) = 0
C     ADRESSE -1 DE LA 1-ERE COORDONNEE DU 1-ER SOMMET DE L'UNION DE PLSV
      MNSU = MNSOMM + WYZSOM - 3

      PRINT*,'unplsv: EPZERO=',EPZERO,' EPSXYZ=',EPSXYZ

      DO 920 J=1,NBOBUN

C        ADRESSE DU TABLEAU 'XYZSOMMET' DU PLSV J
         MN = MCN( MNXYOB - 1 + J )

C        LE NOMBRE DE SOMMETS
         NBS = MCN( MN + WNBSOM )
C
C        ADRESSE DE LA 1-ERE COORDONNEE DU 1-ER SOMMET
         MNS = MN + WYZSOM - 3

C        LE NUMERO DE LA DERNIERE TG DU PLSV J
         NUTG = NUTG + MCN( MN + WNBTGS )
         MCN(MNDETG+J) = NUTG

C        ADRESSE -1 DU 1-ER NUMERO DE SOMMET DU PLSV UNION
         MNUN = MNUNIO+WDSCOU+NBOBUN+MCN(MNDESO+J-1)

         DO 900 I=1,NBS
C           LE NUMERO DES PAVES SUSCEPTIBLES DE CONTENIR CE SOMMET I
            MNS = MNS + 3

C           RECHERCHE DES INDICES DU PAVAGE EN X POUVANT CONTENIR LE SOMMET
            X   = RMCN(MNS)
            D   = ( X - XMINS ) * ECHPAX
            LX  = INT( D )
            D   = D - LX
            IF( D .LE. EPS ) THEN
               LX2 = LX
               LX1 = MAX( 0, LX - 1)
            ELSE IF( D .GE. UNMEPS ) THEN
               LX1 = LX
               LX2 = MIN( NBIPAX, LX+1 )
            ELSE
               LX1 = LX
               LX2 = LX
            ENDIF

C           RECHERCHE DES INDICES DU PAVAGE EN Y POUVANT CONTENIR LE SOMMET
            Y   = RMCN(MNS+1)
            D   = ( Y - YMINS ) * ECHPAY
            LY  = INT( D )
            D   = D - LY
            IF( D .LE. EPS ) THEN
               LY2 = LY
               LY1 = MAX( 0, LY - 1)
            ELSE IF( D .GE. UNMEPS ) THEN
               LY1 = LY
               LY2 = MIN( NBIPAY, LY+1 )
            ELSE
               LY1 = LY
               LY2 = LY
            ENDIF

C           RECHERCHE DES INDICES DU PAVAGE EN Z POUVANT CONTENIR LE SOMMET
            Z   = RMCN(MNS+2)
            D   = ( Z - ZMINS ) * ECHPAZ
            LZ  = INT( D )
            D   = D - LZ
            IF( D .LE. EPS ) THEN
               LZ2 = LZ
               LZ1 = MAX( 0, LZ - 1)
            ELSE IF( D .GE. UNMEPS ) THEN
               LZ1 = LZ
               LZ2 = MIN( NBIPAZ, LZ+1 )
            ELSE
               LZ1 = LZ
               LZ2 = LZ
            ENDIF

C           PARCOURS DES PAVES
            DO NZ = LZ1, LZ2
               MZ = NZ * NBIPY1
               DO NY = LY1, LY2
                  MY = ( NY + MZ ) * NBIPX1
                  DO NX = LX1, LX2
                     NUPAVE = NX + MY

C                    LE PREMIER SOMMET DU PAVE NUPAVE
                     NOSOMM = MCN( MNPAVE + NUPAVE )
C
 8                   IF( NOSOMM .GT. 0 ) THEN

C                       LE SOMMET NOSOMM EXISTE : COMPARAISON DES 3 COORDONNEES
                        MN1 = MNSU + 3 * NOSOMM
                        CALL XYZIDE( RMCN(MN1), RMCN(MNS), IDENTQ )
                        IF( IDENTQ .EQ. 0 ) GOTO 30

C                       LES 2 SOMMETS XYZ(I) DU PLSV J et NOSOMM SONT IDENTIFIES
                        R1 = RMCN(MNS  )
                        R2 = RMCN(MNS+1)
                        R3 = RMCN(MNS+2)
C                       AFFICHAGE SI LE 4-EME CHIFFRE SIGNIFICATIF EST DIFFERENT
                        IF( ABS(RMCN(MN1  )-R1) .GT. 5E-4*ABS(R1) .OR.
     %                      ABS(RMCN(MN1+1)-R2) .GT. 5E-4*ABS(R2) .OR.
     %                      ABS(RMCN(MN1+2)-R3) .GT. 5E-4*ABS(R3) ) THEN
                           IF( LANGAG .EQ. 0 ) THEN
                             WRITE(IMPRIM,10020)(RMCN(MNS+K),K=0,2),I,J,
     %                                        (RMCN(MN1+K),K=0,2),NOSOMM
                            ELSE
                             WRITE(IMPRIM,20020)(RMCN(MNS+K),K=0,2),I,J,
     %                                        (RMCN(MN1+K),K=0,2),NOSOMM
                            ENDIF
                        ENDIF
                        GOTO 800

C                       LE SOMMET N'EST PAS XYZ . PASSAGE AU SOMMET SUIVANT
  30                    NOSOMM = MCN( MNCHA1 + NOSOMM )
                        GOTO 8

                     ENDIF
                  ENDDO
               ENDDO
            ENDDO

10020 FORMAT(' unplsv: REMARQUE  :  SOMMET=',3G15.7,' NUMERO:',I9,
     %               ' du PLSV',I6/
     %       ' unplsv: IDENTIFIE AU SOMMET=',3G15.7,' NUMERO:',I9,
     %               ' FINAL'/)
20020 FORMAT(' unplsv: REMARK:       VERTEX=',3G15.7,' NUMBER:',I9,
     %               ' of PLSV',I6/
     %       ' unplsv: IDENTIFIED to VERTEX=',3G15.7,' NUMBER:',I9,
     %               ' FINAL'/)

C           SOMMET NON IDENTIFIE . AJOUT DE CE SOMMET
            NBSOM  = NBSOM + 1
            NOSOMM = NBSOM
C           CHAINAGE DU NOUVEAU SOMMET EN DEBUT DE PAVE
            NUPAVE = LX + NBIPX1 * ( LY + NBIPY1 * LZ )
            MCN( MNCHA1+NOSOMM ) = MCN( MNPAVE+NUPAVE )
            MCN( MNPAVE+NUPAVE ) = NOSOMM
C           LES 3 COORDONNEES
            MN1 = MNSU + 3 * NOSOMM
            RMCN( MN1     ) = RMCN( MNS     )
            RMCN( MN1 + 1 ) = RMCN( MNS + 1 )
            RMCN( MN1 + 2 ) = RMCN( MNS + 2 )

C           LE NUMERO DU SOMMET I DU PLSV J DANS L'UNION EST NOSOMM
 800        MCN( MNUN + I ) = NOSOMM

 900     ENDDO
 920  ENDDO

C     LE NOMBRE DE COORDONNEES PAR SOMMET
      MCN( MNSOMM + WBCOOR ) = 3

C     LE NOMBRE DE SOMMETS DE L'UNION
      MCN( MNSOMM + WNBSOM ) = NBSOM
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10920) NBSOM
      ELSE
         WRITE(IMPRIM,20920) NBSOM
      ENDIF
10920 FORMAT(' unplsv: NOMBRE SOMMETS   APRES IDENTIFICATION=',I9)
20920 FORMAT(' unplsv: AFTER  IDENTIFICATION VERTEX  NUMBER=',I9)

C     IDENTIFICATION DES TANGENTES
C     ============================
C     ADRESSE -1 DU 1-ER NUMERO DE SOMMET DE L'OBJET UNION
      MNDETG = MNUNIO + WDSCOU + 1+NBOBUN + NBSOMM
      MNNUTG = MNDETG + 1 + NBOBUN
      CALL IDETGS( NBOBUN, MNXYOB, MNSOMM, NBTGUN,
     %             MNDETG, MNNUTG, NBTGS )
C     LE NOMBRE DE TANGENTES DE L'OBJET UNION
      MCN( MNSOMM + WNBTGS ) = NBTGS

C     MISE A JOUR DU TABLEAU 'XYZSOMMET' DU PLSV UNION
C     ------------------------------------------------
C     LA DATE ET LE NUMERO DE TD
      CALL ECDATE( MCN(MNSOMM) )
      MCN( MNSOMM + MOTVAR(6) ) = NONMTD( '~>>>XYZSOMMET' )
C     LE TMS XYZSOMMET EST RACCOURCI A SA VALEUR EXACTE
      CALL TAMSRA( NTSOMM, WYZSOM + ( NBSOM + NBTGS ) * 3 )

C     MISE A JOUR DU TABLEAU 'UNION' DU PLSV UNION
C     --------------------------------------------
C     LE NOMBRE DE PLSV DE L'UNION
      MCN( MNUNIO + WBOBUN ) = NBOBUN
C     COPIE DU TABLEAU POINTEUR SUR LE DERNIER SOMMET DE CHAQUE PLSV
      CALL TRTATA( MCN(MNDESO), MCN(MNUNIO+WDSCOU), 1+NBOBUN )
C     LE NUMERO APRES IDENTIFICATION DE CHAQUE TANGENTE
C     A DEJA ETE CALCULE DANS CALL IDETGS
C     LA DATE ET LE NUMERO DE TD
      CALL ECDATE( MCN(MNUNIO) )
      MCN( MNUNIO + MOTVAR(6) ) = NONMTD( '~>>>UNION' )

C     MISE A JOUR DU TABLEAU 'NSEF' DU PLSV UNION
C     -------------------------------------------
C     DECLARATION OUVERTURE DU TABLEAU 'NSEF' DU PLSV UNION
      IF( NUTYOB .GT. 1 ) THEN
         CALL LXTSOU( NTLXOB, 'NSEF', NTTSMA, MNTSMA )
         IF( NTTSMA .GT. 0 ) THEN
C           DESTRUCTION DU TABLEAU
            CALL LXTSDS( NTLXOB, 'NSEF' )
         ENDIF
         MOT = WUSOEF + NBSOEF*NBEFOB
C        LE NOMBRE DE POINTEURS SUR LES EF POUR LES EF A TG
         IF( NBTGS .GT. 0 ) THEN
            NBEFAP = NBEFOB
            MOT    = MOT + NBEFOB + (1+NBTGEU) * NBEFTG
         ELSE
            NBEFAP = 0
         ENDIF
         CALL LXTNDC( NTLXOB, 'NSEF', 'MOTS', MOT )
         CALL LXTSOU( NTLXOB, 'NSEF', NTTSMA, MNTSMA )
C
C        LE TYPE DU PLSV UNION
         MCN( MNTSMA + WUTYOB ) = NUTYOB
C        LE TYPE DU MAILLAGE : ICI NON STRUCTURE
         MCN( MNTSMA + WUTYMA ) = 0
C        LE NOMBRE DE SOMMETS PAR EF
         MCN( MNTSMA + WBSOEF ) = NBSOEF
C        LE NOMBRE DE TANGENTES PAR EF DANS L'UNION
         MCN( MNTSMA + WBTGEF ) = NBTGEU
C        LE NOMBRE D'EF DU PLSV UNION
         MCN( MNTSMA + WBEFOB ) = NBEFOB
C        LE NOMBRE D'EF A TG DU PLSV UNION
         MCN( MNTSMA + WBEFTG ) = NBEFTG
C        LE NOMBRE DE POINTEURS SUR LES EF POUR LES EF A TG
         MCN( MNTSMA + WBEFAP ) = NBEFAP
C        LE TABLEAU DU NUMERO DES SOMMETS DE L'UNION
         MNMA = MNTSMA + WUSOEF - 1
C        LES TABLEAUX EVENTUELS SUR LES EF A TG DE L'UNION
         MNLDAP = MNTSMA + WUSOEF + NBSOEF * NBEFOB - 1
         MNLDNG = MNLDAP + NBEFAP
         MNLDTG = MNLDNG + NBEFTG - NBTGEU

C        SI UNION G1-CONTINUE ALORS
C           POUR CHAQUE SOMMET STOCKAGE DE -1 SI SOMMET APPARTENANT A PLUSIEURS
C                                          >=0 SI APPARTIENT A UN SEUL PLSV
         IF( NTYPEU .EQ. 52 ) THEN
C           DECLARATION DU TABLEAU NOMBRE PLSV DES SOMMETS
            CALL TNMCDC( 'ENTIER', 1+NBSOM, MNNBSU )
C           MISE A ZERO
            CALL AZEROI( 1+NBSOM, MCN(MNNBSU) )
         ENDIF

C        GENERATION DU TABLEAU 'MATERIAUX' DE L'UNION CAR PLUSIEURS MATERIAUX
C        --------------------------------------------------------------------
C        AJOUT 12 JUIN 2008
         CALL LXTSOU( NTLXOB, 'MATERIAUX', NTMATE, MNMATE )
         IF( NTMATE .GT. 0 ) CALL LXTSDS( NTLXOB, 'MATERIAUX' )
         MOTS2M = WUDMEF + NBEFOB
         CALL LXTNDC( NTLXOB, 'MATERIAUX', 'ENTIER', MOTS2M )
         CALL LXTSOU( NTLXOB, 'MATERIAUX',  NTMATE , MNMATE )
C        LES PREMIERES VARIABLES DU TABLEAU 'MATERIAUX'
C        NOMBRE DE MATERIAUX
         MCN( MNMATE + WNBDM ) = NBOBUN
C        NOMBRE D'ELEMENTS FINIS DU MAILLAGE
         MCN( MNMATE + WBDMEF ) = NBEFOB
C        ADRESSE MCN DE DEBUT DU TABLEAU NUDMEF(1..NBDMEF)
C        'NUMERO DU MATERIAU DE CHAQUE EF'
         MNUDMEF = MNMATE + WUDMEF

C        MISE A ZERO DU TABLEAU NOSOEL SI LECTURE INCOMPLETE
         CALL AZEROI( 64, NOSOEL )

C        NOMBRE ACTUEL DE TANGENTES DANS L'UNION
         NBTU   = 0
         NE     = 0
         NETG   = 0
         NBEFPB = 0
         DO 1200 J=1,NBOBUN

C           LE NO DE MATERIAU = NUMERO DU PLSV DANS SON LEXIQUE
            NOMATE = NUMOBJ( J )

C           LE NOMBRE DE TANGENTES DU PLSV J AVANT IDENTIFICATION
            NBT = MCN( MCN(MNXYOB-1+J) + WNBTGS )

C           L'ADRESSE DU TABLEAU 'NSEF' DU PLSV J
            MN = MCN( MNNSEF - 1 + J )

C           LES PARAMETRES DES NO SOMMET DU MAILLAGE DU PLSV J
            CALL NSEFPA( MCN(MN),
     %                   NUTYML, NBSOEL, NBSOEF, NBTGEF,
     %                   LDAPEF, LDNGEF, LDTGEF, NBEF,
     %                   NX   , NY   , NZ   ,
     %                   IER   )
            IERR = IERR + IER

C           LE NOMBRE D'EF A TG DU PLSV J DE L'UNION
            NBEFT = MCN( MN + WBEFTG )

C           ADRESSE -1 DU 1-ER NUMERO DE SOMMET DU PLSV UNION
            MNUN = MNUNIO+WDSCOU+NBOBUN+MCN(MNDESO+J-1)

            DO 1100 I=1,NBEF

C              UN ELEMENT FINI DE L'UNION EN PLUS
               NE = NE + 1

C              LE NO DE MATERIAU DE CET EF
               MCN( MNUDMEF-1+NE ) = NOMATE

C              LE NUMERO DES NBSOEF SOMMETS ET DES NBTGEF TGS DE L'EF I
               CALL NSEFNS( I     , NUTYML, NBSOEF, NBTGEF,
     %                      LDAPEF, LDNGEF, LDTGEF,
     %                      MN    , NX, NY, NZ  ,
     %                      NCOGEL, NUGEEF, NUEFTG, NOSOEL, IER )
               IERR = IERR + IER

C              EN SORTIE NBSOEF NUMEROS DE SOMMETS SONT INITIALISES
C                        LES DERNIERS PEUVENT ETRE NULS
C                        PAR EXEMPLE LE QUATRIEME D'UN TRIANGLE
C                                    LE CINQUIEME D'UN TETRAEDRE ...
C                        DE MEME POUR LES NUMEROS DES TANGENTES

C              TRANSFORMATION DU NUMERO DES SOMMETS DE CET EF I
C              DU PLSV J EN LES NUMEROS DE SOMMETS DU PLSV UNION
               LEFAPB = 0
               DO 1020 K=1,NBSOEF

C                 SON NUMERO DE SOMMET DANS SON PLSV AVANT UNION
                  NS = NOSOEL(K)
                  IF( NS .EQ. 0 ) GOTO 1015

C                 SON NUMERO DE SOMMET DANS L'UNION APRES IDENTIFICATION
                  NS = MCN( MNUN + NS )

                  IF( NS .LE. 0 .OR. NS .GT. NBSOM ) THEN
C                    ERREUR DETECTEE : NO DE SOMMET INCORRECT
                     NBLGRC(NRERR) = 5
C                    LE NOM DU PLSV
                     CALL NMOBNU( NMTYOB(NUTYOB), NUMOBJ(J), KNOM )
                     KERR(1) = 'PLSV ' // KNOM
                     WRITE(KERR(MXLGER)(1:6),'(I6)') I
                     WRITE(KERR(MXLGER-1)(1:9),'(I9)') NOSOEL(K)
                     WRITE(KERR(MXLGER-1)(11:12),'(I2)') K
                     WRITE(KERR(MXLGER-2)(1:9),'(I9)') NS
                     WRITE(KERR(MXLGER-3)(1:9),'(I9)') NBSOM
                     IF( LANGAG .EQ. 0 ) THEN
                        KERR(2) = 'EF '// KERR(MXLGER)(1:6) //
     %                            'AVEC NUMERO DE SOMMET INCORRECT'
                        KERR(3) = 'AVANT UNION NO SOMMET '
     %                         // KERR(MXLGER-1)(11:12)
     %                         // ' = ' // KERR(MXLGER-1)(1:9)
                        KERR(4) = 'APRES UNION NO SOMMET '
     %                         // KERR(MXLGER-1)(11:12)
     %                         // ' = ' // KERR(MXLGER-2)(1:9)
                        KERR(5) = 'ET FINALEMENT IMPOSE A '
     %                         // KERR(MXLGER-3)(1:9)
                     ELSE
                        KERR(2) = 'FE '// KERR(MXLGER)(1:6) //
     %                            'with INCORRECT VERTEX NUMBER'
                        KERR(3) = 'Before UNION VERTEX NUMBER '
     %                         // KERR(MXLGER-1)(11:12)
     %                         // ' = ' // KERR(MXLGER-1)(1:9)
                        KERR(4) = 'After  UNION VERTEX NUMBER '
     %                         // KERR(MXLGER-1)(11:12)
     %                         // ' = ' // KERR(MXLGER-2)(1:9)
                        KERR(5) = 'and finally IMPOSED to '
     %                         // KERR(MXLGER-3)(1:9)
                     ENDIF
                     CALL LEREUR
                     IERR = IERR + 1
C                    NUMERO DE SOMMET IMPOSE A NBSOM POUR MONTRER L'ERREUR
                     NS = NBSOM
                  ENDIF

C                 EXISTE T IL UN SOMMET DE MEME NUMERO NS
C                 => MAUVAISE IDENTIFICATION DE SOMMET
                  DO 1012 L=1,K-1
                     IF( NS .EQ. MCN(MNMA+L) ) THEN
C                       LE SOMMET L EST AUSSI NS => ERREUR
                        LEFAPB = 1
                     ENDIF
 1012             ENDDO

C                 STOCKAGE DU NOUVEAU NO DE SOMMET POUR L'UNION
 1015             MCN( MNMA + K ) = NS
C
                  IF( NTYPEU .EQ. 52 ) THEN
C                    LE NUMERO DU DERNIER PLSV DE L'UNION DE CE SOMMET NS
                     NUPLSV = MCN( MNNBSU + NS )
                     IF( NUPLSV .EQ. 0 ) THEN
C                       STOCKAGE DU NUMERO DU PLSV DANS L'UNION
                        MCN( MNNBSU + NS ) = J
                        GOTO 1020
                     ENDIF
                     IF( NUPLSV .NE. J ) THEN
C                       CE SOMMET APPARTIENT A PLUSIEURS PLSV DE L'UNION
                        MCN( MNNBSU + NS ) = -1
                     ENDIF
                  ENDIF

 1020          ENDDO

               IF( LEFAPB .NE. 0 ) THEN

C                 ERREUR DETECTEE: AU MOINS UN SOMMET DOUBLE POUR L'EF
                  IF( NBEFPB .EQ. 0 ) THEN
C                    INITIALISATION DU TRACE DE L'EF A PB
                     TRACTE = .TRUE.
C                    AUCUN ITEM SUR L'ECRAN
                     CALL EFFACE
                     CALL ITEMS0
C                    LES PARAMETRES DU CADRE MAXIMAL
                     CALL VISEE0
                  ENDIF

C                 UN EF A PB DE PLUS
                  NBEFPB = NBEFPB + 1

C                 LES COORDONNEES DES SOMMETS DE L'EF
                  MM  = MNSOMM + WYZSOM - 4
                  NBS = NBSOME( NCOGEL )
                  DO 1025 L=1,NBS
C                    LE NUMERO DU SOMMET ET SES XYZ SONT STOCKES DANS X Y Z
                     NOSM     = MCN( MNMA + L ) * 3
                     XX(L)    = RMCN(MM+NOSM+1)
                     XYZ(1,L) = XX(L)
                     YY(L)    = RMCN(MM+NOSM+2)
                     XYZ(2,L) = YY(L)
                     ZZ(L)    = RMCN(MM+NOSM+3)
                     XYZ(3,L) = ZZ(L)
 1025             ENDDO

C                 TRACE DE L'EF PERDU
                  IF( NDIMLI .EQ. 3 ) THEN
C
C                    TRACE EN DIMENSION 3
                     GOTO(8015, 8020, 8030,8030, 8045,8046,8047),NCOGEL
C
C                    UN POINT=SOMMET
 8015                IF( LANGAG .EQ. 0 ) THEN
                        CALL SYMBOLE3D(NCROUG, XYZ(1,1),'X PoinT Perdu')
                     ELSE
                        CALL SYMBOLE3D( NCROUG, XYZ(1,1),'X Lost Point')
                     ENDIF
                     GOTO 8090

C                    UNE ARETE
 8020                CALL  TRAIT3D( NCROUG, XYZ(1,1), XYZ(1,2) )
                     GOTO 8090

C                    UN TRIANGLE OU QUADRANGLE
 8030                CALL FAP13D( NCROUG, NCNOIR, 0, NCOGEL, XX,YY,ZZ )
                     GOTO 8090

C                    UN TETRAEDRE
 8045                CALL  TRAIT3D( NCROUG, XYZ(1,1), XYZ(1,2) )
                     CALL  TRAIT3D( NCROUG, XYZ(1,2), XYZ(1,3) )
                     CALL  TRAIT3D( NCROUG, XYZ(1,3), XYZ(1,1) )
                     CALL  TRAIT3D( NCROUG, XYZ(1,1), XYZ(1,4) )
                     CALL  TRAIT3D( NCROUG, XYZ(1,2), XYZ(1,4) )
                     CALL  TRAIT3D( NCROUG, XYZ(1,3), XYZ(1,4) )
                     GOTO 8090

C                    UN PENTAEDRE
 8046                CALL  TRAIT3D( NCROUG, XYZ(1,1), XYZ(1,2) )
                     CALL  TRAIT3D( NCROUG, XYZ(1,2), XYZ(1,3) )
                     CALL  TRAIT3D( NCROUG, XYZ(1,3), XYZ(1,1) )
                     CALL  TRAIT3D( NCROUG, XYZ(1,1), XYZ(1,4) )
                     CALL  TRAIT3D( NCROUG, XYZ(1,2), XYZ(1,5) )
                     CALL  TRAIT3D( NCROUG, XYZ(1,3), XYZ(1,6) )
                     CALL  TRAIT3D( NCROUG, XYZ(1,4), XYZ(1,5) )
                     CALL  TRAIT3D( NCROUG, XYZ(1,5), XYZ(1,6) )
                     CALL  TRAIT3D( NCROUG, XYZ(1,6), XYZ(1,4) )
                     GOTO 8090

C                    UN HEXAEDRE
 8047                CALL  TRAIT3D( NCROUG, XYZ(1,1), XYZ(1,2) )
                     CALL  TRAIT3D( NCROUG, XYZ(1,2), XYZ(1,3) )
                     CALL  TRAIT3D( NCROUG, XYZ(1,3), XYZ(1,4) )
                     CALL  TRAIT3D( NCROUG, XYZ(1,4), XYZ(1,1) )
                     CALL  TRAIT3D( NCROUG, XYZ(1,1), XYZ(1,5) )
                     CALL  TRAIT3D( NCROUG, XYZ(1,2), XYZ(1,6) )
                     CALL  TRAIT3D( NCROUG, XYZ(1,3), XYZ(1,7) )
                     CALL  TRAIT3D( NCROUG, XYZ(1,4), XYZ(1,8) )
                     CALL  TRAIT3D( NCROUG, XYZ(1,5), XYZ(1,6) )
                     CALL  TRAIT3D( NCROUG, XYZ(1,6), XYZ(1,7) )
                     CALL  TRAIT3D( NCROUG, XYZ(1,7), XYZ(1,8) )
                     CALL  TRAIT3D( NCROUG, XYZ(1,8), XYZ(1,5) )

                  ELSE

C                    TRACE EN DIMENSION 2
                     GOTO(8050, 8060, 8070, 8070, 8090,8090,8090),NCOGEL

C                    UN POINT=SOMMET
 8050                IF( LANGAG .EQ. 0 ) THEN
                      CALL SYMBOLE2D(NCROUG,XX(1),YY(1),'X PoinT PERDU')
                     ELSE
                       CALL SYMBOLE2D(NCROUG,XX(1),YY(1),'X Lost Point')
                     ENDIF
                     GOTO 8090

C                    UNE ARETE
 8060                CALL  TRAIT2D( NCROUG, XX(1), YY(1), XX(2), YY(2) )
                     GOTO 8090

C                    UN TRIANGLE OU QUADRANGLE
 8070                CALL FACE2D( NCROUG, NCNOIR, NCOGEL, XX, YY )

                  ENDIF

C                 IMPRESSION DES NO ET XYZ DES SOMMETS DE L'EF
 8090             WRITE(IMPRIM,18089)
                  IF( LANGAG .EQ. 0 ) THEN
                     WRITE(IMPRIM,18090) NE
                     WRITE(IMPRIM,18091)
     %                 (MCN(MNMA+L),XX(L),YY(L),ZZ(L),L=1,NBS)
                  ELSE
                     WRITE(IMPRIM,28090) NE
                     WRITE(IMPRIM,28091)
     %                 (MCN(MNMA+L),XX(L),YY(L),ZZ(L),L=1,NBS)
                  ENDIF
                  WRITE(IMPRIM,18089)
18089    FORMAT(80('*'))
18090    FORMAT(' EF',I9,' A PROBLEME AVEC AU MOINS UN SOMMET DOUBLE')
18091    FORMAT(' SOMMET',I9,' : X=',G15.7,' Y=',G15.7,' Z=',G15.7)
28090    FORMAT(' FE',I9,' with at LEAST a DOUBLE VERTEX')
28091    FORMAT(' VERTEX',I9,' : X=',G15.7,' Y=',G15.7,' Z=',G15.7)
C
C                 LE NOM DU PLSV
                  CALL NMOBNU( NMTYOB(NUTYOB), NUMOBJ(J), KNOM )
                  WRITE(KERR(MXLGER)(1:6),'(I6)') I
                  KERR(1) = 'PLSV ' // KNOM
                  IF( LANGAG .EQ. 0 ) THEN
                     KERR(2) = 'EF '// KERR(MXLGER)(1:6) //
     %                         ' AVEC AU MOINS UN SOMMET DOUBLE'
                     KERR(3) = 'REVOIR PRECISION pour IDENTIFIER POINTS'
                     KERR(4) = 'DE L''OPTION 20 du MENU DEBUT'
                  ELSE
                     KERR(2) = 'FE '// KERR(MXLGER)(1:6) //
     %                         ' with at LEAST a DOUBLE VERTEX'
                     KERR(3) = 'cf PRECISION to IDENTIFY POINTS'
                     KERR(4) = 'cf OPTION 20 of the DEBUT MENU'
                  ENDIF
                  NBLGRC(NRERR) = 4
                  CALL LEREUR
                  IERR = IERR + 1
               ENDIF

C              MISE A JOUR DE L'ADRESSE
               MNMA = MNMA + NBSOEF

               IF( NBTGS .GT. 0 ) THEN

C                 IL EXISTE DES EF A TG DANS LE PLSV UNION
C                 ========================================
                  IF( NBEFT .GT. 0 ) THEN
C                    IL EXISTE DES EF A TG POUR LE PLSV J DE L'UNION
C                    L'EF DE CE PLSV EST IL UN EF A TG?
                     IF( NUEFTG .GT. 0 ) THEN
C                       OUI EF A TG : IL EST AJOUTE
                        NETG = NETG + 1
C                       LE POINTEUR OU NUMERO DE L'EF A TG DANS L'UNION
                        MCN( MNLDAP + NE ) = NETG
C                       LE NUMERO GEOMETRIQUE DE CET EF A TG DANS L'UNION
                        MCN( MNLDNG + NETG ) = NUGEEF
C                       LE TABLEAU NUTGEF EST PRESENT DANS NSEF DE CE PLSV
C                       ADRESSE MCN DU NUMERO QUI PRECEDE LA PREMIERE TG DE CET
                        MNTG = MNLDTG + NETG * NBTGEU
                        DO 1030 K=1,NBTGEU
C                          LE NUMERO DE LA TANGENTE DANS SON PLSV
                           NT = NOSOEL( NBSOEF + K )
                           IF( NT .EQ. 0 ) THEN
                              LSIGNE = 0
                              GOTO 1029
                           ELSE IF( NT .LT. 0 ) THEN
                              LSIGNE = -1
                           ELSE
                              LSIGNE = 1
                           ENDIF
C                          LE NUMERO DE LA TANGENTE DANS L'UNION DES PLSV
                           NT = MCN( MNNUTG-1+NBTU+ABS(NT) )
                           IF( NT .LT. 0 ) THEN
                              LSIGNE = -LSIGNE
                           ENDIF
C                          STOCKAGE DU NUMERO DE LA TANGENTE DANS L'UNION
 1029                      MCN( MNTG + K ) =  LSIGNE * ABS( NT )
 1030                   ENDDO
                     ELSE
C                       EF SANS TG
C                       LE POINTEUR OU NUMERO DE L'EF A TG DANS L'UNION
                        MCN( MNLDAP + NE  ) = 0
                     ENDIF
                  ELSE
C                    PLSV SANS TG MAIS UNION A TG
C                    LE POINTEUR
                     MCN( MNLDAP + NE ) = 0
                  ENDIF
               ENDIF
C
 1100       ENDDO

C           MISE A JOUR DU NOMBRE ACTUEL DE TANGENTES DE L'UNION
            NBTU = MCN(MNDETG+J)

            TRACTE = TRACTE0
 1200    ENDDO

         IF( NBEFPB .NE. 0 ) IERR = IERR + 1
         IF( NBEFTG .NE. NETG ) THEN
            IF( LANGAG .EQ. 0 ) THEN
               WRITE(IMPRIM,*) 'unplsv: MAUVAIS CALCUL DES EF A TG'
            ELSE
               WRITE(IMPRIM,*) 'unplsv: BAD COMPUTATION of FE with TG'
            ENDIF
            IERR = IERR + 2
         ENDIF
         IF( IERR .NE. 0 ) GOTO 9000

C        LE TYPE INCONNU DE FERMETURE DU MAILLAGE
         MCN( MNTSMA + WUTFMA ) = -1

C        LA DATE ET LE NUMERO DE TD
         CALL ECDATE( MCN(MNTSMA) )
         MCN( MNTSMA + MOTVAR(6) ) = NONMTD( '~>>>NSEF' )

C        LA DATE DE CREATION DU TABLEAU 'MATERIAUX' DE L'UNION
         CALL ECDATE( MCN(MNMATE) )
C        AJOUT DU NUMERO DU TABLEAU DESCRIPTEUR
         MCN( MNMATE + MOTVAR(6) ) = NONMTD( '~>>>MATERIAUX' )
      ELSE
         NTTSMA = 0
         MNTSMA = 0
      ENDIF

      IF( NUTYOB .EQ. 2 ) THEN
C        LA LIGNE ACTUELLEMENT NON STRUCTUREE EST ELLE STRUCTURABLE ?
C        ATTENTION: SI LA LIGNE EST FERMEE (DONC NON STRUCTURABLE)
C        LES ARETES ET SOMMETS SONT REMIS DANS UN ORDRE DE PARCOURS FAVORABLE
C        --------------------------------------------------------------------
         CALL LIGSTR( NTLXOB, NTTSMA, MNTSMA, NTSOMM, MNSOMM, IER )
         IERR = IERR + IER
      ENDIF

      IF( NUTYOB .EQ. 3 ) THEN

C        TOUTE ARETE C1 TRANSFORME SON EF C0 EN C1
C        -----------------------------------------
         CALL NMOBNU( 'SURFACE', NUOBUN, KNOM )
         CALL SUC0C1( KNOM,   NTLXOB,
     %                NTTSMA, MNTSMA, NTSOMM, MNSOMM, IER )
         IERR = IERR + IER

         IF( NTYPEU .EQ. 52 ) THEN
C
C           DEMANDE DE SURFACE G1-CONTINUE PAR PROJECTION DES TGS D'UN SOMMET
C           SUR LE PLAN DE DISTANCE MINIMALE DES EXTREMITES DE SES TANGENTES
C           DES SEULS SOMMETS COMMUNS A PLUSIEURS PLSV DE L'UNION
C           -----------------------------------------------------------------
            CALL SUC0G1( NTLXOB, NBSOM,  MCN(MNNBSU),
     %                   NTTSMA, MNTSMA, NTSOMM, MNSOMM, IER )
            IERR = IERR + IER

         ENDIF

C        ASSURER L'ORIENTATION DU MAILLAGE DE LA SURFACE
C        PAR PARCOURS DES EF A TRAVERS LES ARETES COMMUNES
C        ET PERMUTATION DES SOMMETS 2-NBSTEF DU SECOND EF
C        D'UNE ORIENTATION DIFFERENTE DU PREMIER EF
         CALL SUORIENT( KNOM,
     %                  NTTSMA, MNTSMA, NTSOMM, MNSOMM, IER )
         IERR = IERR + IER

      ENDIF

CCCC     A METTRE A JOUR POUR LE CAS C1!
CCC      IF( (NUTYOB .EQ. 3) .AND. (NBQUST.EQ.NBOBUN) ) THEN
CCCC
CCCC        LA SURFACE UNION DE SURFACES QUADRANGLES STRUCTURES
CCCC        ACTUELLEMENT NON STRUCTUREE EST ELLE STRUCTURABLE ?
CCCC        ---------------------------------------------------
CCCC        DECLARATIONS DE TABLEAUX AUXILIAIRES
CCC         CALL TNMCDC( 'ENTIER', NBEFOB, MNINDX )
CCC         CALL TNMCDC( 'REEL' , 3*NBSOM, MNSOME )
CCC         CALL SUEY40( MCN(MNTSMA+WBSOEF), NBEFOB, MCN(MNTSMA+WUSOEF),
CCC     %                RMCN(MNSOMM+WYZSOM),
CCC     %                NBSOM, MCN(MNINDX), RMCN(MNSOME),
CCC     %                NBA1 , NBA2,IER )
CCC         IERR = IERR + IER
CCC         IF( IERR .EQ. 0 ) THEN
CCCC           LA SURFACE EST STRUCTUREE POUR SES SOMMETS
CCCC           **********************************************
CCCC           ATTENTION : ACTUELLEMENT PERTE DES TANGENTES !
CCCC           **********************************************
CCCC           MISE A JOUR DU TABLEAU  'NSEF'
CCC            CALL LXTSDS( NTLXOB, 'NSEF' )
CCC            CALL LXTNDC( NTLXOB, 'NSEF', 'MOTS', WBARYQ+1 )
CCC            CALL LXTSOU( NTLXOB, 'NSEF', NTTSMA, MNTSMA )
CCCC           LE TYPE DU PLSV UNION
CCC            MCN( MNTSMA + WUTYOB ) = 3
CCCC           LE TYPE NON-FERME DE FERMETURE DU MAILLAGE
CCC            MCN( MNTSMA + WUTFMA ) = 0
CCCC           LE TYPE DU MAILLAGE : ICI STRUCTURE
CCC            MCN( MNTSMA + WUTYMA ) = 4
CCCC           LE NOMBRE DE SOMMETS PAR FACE
CCC            MCN( MNTSMA + WBSOEF ) = 4
CCCC           LE NOMBRE DE TANGENTES STOCKEES PAR EF : SURFACE C0
CCC            MCN( MNTSMA + WBTGEF ) = 0
CCCC           LE NOMBRE D'EF A TG
CCC            MCN( MNTSMA + WBEFTG ) = 0
CCCC           PAS DE POINTEUR SUR LES EF A TG
CCC            MCN( MNTSMA + WBEFAP ) = 0
CCCC           LE NOMBRE D'EF DE LA SURFACE
CCC            MCN( MNTSMA + WBEFOB ) = NBA1 * NBA2
CCCC           LE NOMBRE DE SEGMENTS DANS LES DEUX DIRECTIONS
CCC            MCN( MNTSMA + WBARXQ ) = NBA1
CCC            MCN( MNTSMA + WBARYQ ) = NBA2
CCCC           LA DATE ET LE NUMERO DE TD
CCC            CALL ECDATE( MCN(MNTSMA) )
CCC            MCN( MNTSMA + MOTVAR(6) ) = NONMTD( '~>>>NSEF' )
CCCC
CCCC           MISE A JOUR DU TABLEAU 'XYZSOMMET' DU PLSV UNION
CCC            CALL LXTSDS( NTLXOB, 'XYZSOMMET' )
CCC            CALL LXTNDC( NTLXOB, 'XYZSOMMET', 'MOTS', WYZSOM+3*NBSOM)
CCC            CALL LXTSOU( NTLXOB, 'XYZSOMMET', NTSOMM, MNSOMM )
CCCC           LE NOMBRE DE SOMMETS DE L'UNION
CCC            MCN( MNSOMM + WNBSOM ) = NBSOM
CCCC           LE NOMBRE DE TANGENTES
CCC            MCN( MNSOMM + WNBTGS ) = 0
CCCC           COPIE DES 3 COORDONNEES DES NBSOM SOMMETS
CCC            CALL TRTATA( MCN(MNSOME), MCN(MNSOMM+WYZSOM), 3*NBSOM )
CCCC           LA DATE ET LE NUMERO DE TD
CCC            CALL ECDATE( MCN(MNSOMM) )
CCC            MCN( MNSOMM + MOTVAR(6) ) = NONMTD( '~>>>XYZSOMMET' )
CCC         ENDIF
CCCC        DESTRUCTION DES TABLEAUX AUXILIAIRES
CCC         CALL TNMCDS( 'ENTIER', NBEFOB, MNINDX )
CCC         CALL TNMCDS( 'REEL' , 3*NBSOM, MNSOME )
CCC      ENDIF

      IF( NUTYOB .EQ. 4 ) THEN

C        TOUTE FACE C1 TRANSFORME SON EF C0 EN C1
C        ----------------------------------------
         CALL NMOBNU( 'VOLUME', NUOBUN, KNOM )
         CALL VOC0C1( KNOM,   NTLXOB,
     %                NTTSMA, MNTSMA, NTSOMM, MNSOMM, IER )
         IERR = IERR + IER

      ENDIF

C     PAS D'ERREUR
      IERR = 0
      GOTO 9100

C     ERREUR => DESTRUCTION DU LEXIQUE DU PLSV UNION
C     ----------------------------------------------
 9000 NTTSMA = 0
      MNTSMA = 0
      NTSOMM = 0
      MNSOMM = 0
      NTUNIO = 0
      MNUNIO = 0

C     DESTRUCTION DES TABLEAUX TEMPORAIRES
C     LES TABLEAUX DE SAUVEGARDE DES ADRESSES MCN DES PLSV
 9100 IF( MNNSEF .GT. 0 ) CALL TNMCDS( 'ENTIER', NBOBU0, MNNSEF )
      IF( MNXYOB .GT. 0 ) CALL TNMCDS( 'ENTIER', NBOBU0, MNXYOB )
C     LE TABLEAU POINTEUR SUR LE DERNIER SOMMET DE CHAQUE PLSV
      IF( MNDESO .GT. 0 ) CALL TNMCDS( 'ENTIER', 1 + NBOBU0, MNDESO )
C     DESTRUCTION DU TABLEAU 'PAVES'
      IF( MNPAVE .GT. 0 ) CALL TNMCDS( 'ENTIER', MOPAVE, MNPAVE )
C     DESTRUCTION DU TABLEAU DES CHAINAGES DES SOMMETS DANS LES PAVES
      IF( MNCHAI .GT. 0 ) CALL TNMCDS( 'ENTIER', NBSOMM,  MNCHAI )
      IF( MNNBSU .GT. 0 ) CALL TNMCDS( 'ENTIER', 1+NBSOM, MNNBSU )

C     IDENTIFICATION DES TANGENTES MULTIPLES OU EGALES AUX ARETES
C     ===========================================================
      IF( IERR .EQ. 0 .AND. NUTYOB .GT. 1 ) THEN
         CALL MOINTG( NTTSMA, MNTSMA, NTSOMM, MNSOMM, IERR )
      ENDIF

      IF( IERR .EQ. 0 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,19100) MCN(MNSOMM+WNBTGS)
         ELSE
            WRITE(IMPRIM,29100) MCN(MNSOMM+WNBTGS)
         ENDIF
      ENDIF
19100 FORMAT(' NOMBRE de TANGENTES apres IDENTIFICATION=',I9)
29100 FORMAT(' AFTER IDENTIFICATION TANGENT NUMBER=',I9)

      TRACTE = TRACTE0
      RETURN
      END
