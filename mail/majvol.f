      SUBROUTINE MAJVOL( NUVOLU, IERR )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    METTRE A JOUR LES ELEMENTS FINIS DU VOLUME COMPOSE DE
C -----    PLUSIEURS MATERIAUX (DE NUMERO DANS LE TMS 'MATERIAUX')
C          POUR CHAQUE MATERIAU=VOLUME
C
C          PAR EXEMPLE UNE TETRAEDRISATION DE PLUSIEURS VOLUMES
C          GENERE UN VOLUME AVEC DES TETRAEDRES DANS PLUSIEURS MATERIAUX
C          CETTE SUBROUTINE REGROUPE DANS CHAQUE VOLUME=MATERIAU LES
C          ELEMENTS FINIS DE CE MATERIAU EN ECRASANT LES EF PRECEDENTS
C
C ENTREE :
C --------
C NUVOLU : NUMERO DU VOLUME MULTI-MATERIAUX DANS LE LEXIQUE DES VOLUMES
C
C SORTIE :
C --------
C IERR   : 0  SI PAS D'ERREUR,  >0 SI ERREUR
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC et St Pierre du Perray    Juin 2008
C2345X7..............................................................012
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/a_volume__definition.inc"
      include"./incl/a___nsef.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___materiaux.inc"
C
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
      CHARACTER*24      NMVOLU
C
      MNNUMA = 0
      MNLIPO = 0
      NBSTVL = 0
C
C     ADRESSE DU LEXIQUE DES VOLUMES
      CALL TAMSOU( NTVOLU, MNVOLU )
C
C     LE NUMERO DU VOLUME MULTI_MATERIAUX
      CALL NMOBNU( 'VOLUME', NUVOLU, NMVOLU )
C
C     NO DE LEXIQUE DE CE VOLUME MULTI-MATERIAUX
      CALL LXNLOU( NTVOLU, NUVOLU, NTLXVL, MNLXVL )
C
C     LE TABLEAU 'MATERIAUX' DU VOLUME MULTI_MATERIAUX
      CALL LXTSOU( NTLXVL, 'MATERIAUX', NTMATE, MNMATE )
 5    IF( NTMATE .LE. 0 ) THEN
C        PAS DE MATERIAUX ASSOCIES AUX ELEMENTS FINIS
         NBLGRC(NRERR) = 2
         KERR(1) = NMVOLU
         IF( LANGAG .EQ. 0 ) THEN
            KERR(2) = 'VOLUME MONO MATERIAU'
         ELSE
            KERR(2) = 'VOLUME WITH ONLY ONE MATERIAL'
         ENDIF
         CALL LEREUR
         IERR = 1
         GOTO 9990
      ENDIF
C
C     VRAI VOLUME MULTI-MATERIAUX
C     CHAQUE EF A SON NUMERO DE MATERIAU=VOLUME RANGE A L'ADRESSE MNUDMEF
C     NOMBRE DE MATERIAUX
      NBDM   = MCN( MNMATE + WNBDM )
C     NOMBRE D''ELEMENTS FINIS DU MAILLAGE AVEC NO DE MATERIAU
      NBDMEF = MCN( MNMATE + WBDMEF )
C     ADRESSE MCN DU NO DE MATERIAU DE CHAQUE EF
      MNUDMEF= MNMATE + WUDMEF
C
C     RETROUVER LE NO DES NBDM MATERIAUX=VOLUMES et NOMBRE D'EF DES MATERIAUX
      CALL TNMCDC( 'ENTIER', 2*NBDM, MNNUMA )
      IF( MNNUMA .LE. 0 ) GOTO 9990
      CALL AZEROI( NBDM, MCN(MNNUMA) )
      MNNBEF = MNNUMA + NBDM - 1
C
      WRITE(IMPRIM,*)
      NBDMATE = 0
      DO 20 NEF = 1, NBDMEF
C        NO DU MATERIAU DE L'EF NEF
         NOMATE = MCN( MNUDMEF - 1 + NEF )
C
C        CE NO A T IL ETE RECENSE?
         DO 10 NMA=1,NBDMATE
            IF( NOMATE .EQ.  MCN(MNNUMA-1+NMA) ) THEN
C              UN EF DE CE MATERIAU EN PLUS
               MCN(MNNBEF+NMA) = MCN(MNNBEF+NMA) + 1
               GOTO 20
            ENDIF
 10      CONTINUE
C
C        NOM DU NOUVEAU MATERIAU RETROUVE
         CALL NMOBNU( 'VOLUME', NOMATE, KERR(1) )
C        C'EST LE PREMIER EF DE CE MATERIAU
         MCN(MNNUMA+NBDMATE) = NOMATE
C        UN MATERIAU RECENSE DE PLUS
         NBDMATE = NBDMATE + 1
C        LE PREMIER EF DU MATERIAU
         MCN(MNNBEF+NBDMATE) = 1
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,*) 'MODIFICATION du MATERIAU',NBDMATE,
     %                      ': ',KERR(1)(1:24)
         ELSE
            WRITE(IMPRIM,*) 'MODIFICATION of MATERIAL',NBDMATE,
     %                      ': ',KERR(1)(1:24)
         ENDIF
C        NO DE LEXIQUE DE CE MATERIAU
         CALL LXNLOU( NTVOLU, NOMATE, NTLXMA, MNLXMA )
         IF( NTLXMA .LE. 0 ) THEN
            NBLGRC(NRERR) = 2
            IF( LANGAG .EQ. 0 ) THEN
               KERR(2) = 'N''EST PAS UN VOLUME'
            ELSE
               KERR(2) = 'NOT A VOLUME'
            ENDIF
            CALL LEREUR
            IERR = 2
            GOTO 9990
         ENDIF
         CALL LXTSOU( NTLXMA, 'XYZSOMMET', NTXYZM, MNXYZM )
         IF( NTXYZM .LE. 0 ) THEN
            NBLGRC(NRERR) = 2
            IF( LANGAG .EQ. 0 ) THEN
               KERR(2) = 'VOLUME SANS LE TMS XYZSOMMET'
            ELSE
               KERR(2) = 'VOLUME WITHOUT XYZSOMMET TMS'
            ENDIF
            CALL LEREUR
            IERR = 2
            GOTO 9990
         ENDIF
         CALL LXTSOU( NTLXMA, 'NSEF', NTNSEM, MNNSEM )
         IF( NTNSEM .LE. 0 ) THEN
            NBLGRC(NRERR) = 2
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'VOLUME SANS LE TMS NSEF'
            ELSE
               KERR(1) = 'VOLUME WITHOUT NSEF TMS'
            ENDIF
            CALL LEREUR
            IERR = 2
            GOTO 9990
         ENDIF
 20   CONTINUE
      IF( NBDMATE .NE. NBDM ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'DIFFERENTS NOMBRES DE MATERIAUX'
         ELSE
            KERR(1) = 'DIFFERENT NUMBERS OF MATERIALS'
         ENDIF
         CALL LEREUR
         IERR = 2
         GOTO 9990
      ENDIF
C
C     ICI TOUS LES MATERIAUX XYZSOMMET ET NSEF ONT ETE RETROUVES
C
C     LE TABLEAU 'XYZSOMMET' DU VOLUME MULTI-MATERIAUX
      CALL LXTSOU( NTLXVL, 'XYZSOMMET', NTXYZS, MNXYZS )
      IF( NTXYZS .LE. 0 ) THEN
         NBLGRC(NRERR) = 2
         KERR(1) = NMVOLU
         IF( LANGAG .EQ. 0 ) THEN
            KERR(2) = 'VOLUME SANS LE TMS XYZSOMMET'
         ELSE
            KERR(2) = 'VOLUME WITHOUT XYZSOMMET TMS'
         ENDIF
         CALL LEREUR
         IERR = 3
         GOTO 9990
      ENDIF
C     LE NOMBRE DE SOMMETS DU VOLUME MULTI-MATERIAUX
      NBSTVL = MCN( MNXYZS + WNBSOM )
C
C     MCN(MNLIPO) ( NBSTVL ) NO DU POINT AVANT ET APRES TRI
      MNLIPO = -1
      CALL TNMCDC( 'ENTIER', NBSTVL+1, MNLIPO )
      IF( MNLIPO .LE. 0 ) GOTO 9990
C
C     LE TABLEAU 'NSEF' DU VOLUME MULTI-MATERIAUX
      CALL LXTSOU( NTLXVL, 'NSEF', NTNSEF, MNNSEF )
      IF( NTNSEF .LE. 0 ) THEN
         NBLGRC(NRERR) = 2
         KERR(1) = NMVOLU
         IF( LANGAG .EQ. 0 ) THEN
            KERR(2) = 'VOLUME SANS LE TMS NSEF'
         ELSE
            KERR(2) = 'VOLUME WITHOUT NSEF TMS'
         ENDIF
         CALL LEREUR
         IERR = 4
         GOTO 9990
      ENDIF
C
C     LES CARACTERISTIQUES DES EF DU VOLUME MULTI-MATERIAUX
      CALL NSEFPA( MCN(MNNSEF) ,
     %             NUTYMA, NBSOEL, NBSOEF, NBTGEF,
     %             LDAPEF, LDNGEF, LDTGEF, NBEFVL,
     %             NX,     NY,     NZ,     IERR   )
      IF( IERR .NE. 0 ) GOTO 9990
C     NUTYMA : 'NUMERO DE TYPE DU MAILLAGE'    ENTIER
C              0 : 'NON STRUCTURE'       , 2 : 'SEGMENT    STRUCTURE',
C              3 : 'TETRAEDRE  STRUCTURE', 4 : 'TETRAEDRE  STRUCTURE',
C              5 : 'TETRAEDRE STRUCTURE' , 6 : 'PENTAEDRE  STRUCTURE',
C              7 : 'HEXAEDRE  STRUCTURE'
      IF( NUTYMA .NE. 0 ) THEN
C        TOUT VOLUME STRUCTURE EST MONO MATERIAU
         NTMATE = 0
         GOTO 5
      ENDIF
C     NBSOEL : NOMBRE DE SOMMETS DES EF (0 SI MAILLAGE NON STRUCTURE)
C     NBSOEF : NOMBRE DE SOMMETS DE STOCKAGE DE CHAQUE EF DE NSEF
C              ( TETRAEDRE NBSOEF=8 )
C     NBEFVL : NOMBRE DE EF DU MAILLAGE DU VOLUME MULTI-MATERIAUX
C     NX, NY, NZ : LE NOMBRE D'ARETES DANS LES DIRECTION X Y Z
C                CF LE TMS ~td/d/a___nsef
C
      IF( NBDMEF .NE. NBEFVL ) THEN
         NBLGRC(NRERR) = 2
         KERR(1) = NMVOLU
         IF( LANGAG .EQ. 0 ) THEN
            KERR(2) = 'TMS NSEF et MATERIAUX INCOMPATIBLES'
         ELSE
            KERR(2) = 'TMS NSEF and MATERIAUX are INCOMPATIBLE'
         ENDIF
         CALL LEREUR
         IERR = 5
         GOTO 9990
      ENDIF
C
C     BOUCLE SUR LES NBDM MATERIAUX DU VOLUME
      DO 100 NMA=1, NBDM
C
C        NO DU MATERIAU
         NOMATE = MCN(MNNUMA-1+NMA)
C
C        NO DE LEXIQUE DE CE MATERIAU
         CALL LXNLOU( NTVOLU, NOMATE, NTLXMA, MNLXMA )
C        PAS DE VERIFICATION CAR FAITES AVANT
C
C        DESTRUCTION DU TMS XYZSOMMET et NSEF DU MATERIAU
         CALL LXTSDS( NTLXMA, 'XYZSOMMET' )
         CALL LXTSDS( NTLXMA, 'NSEF' )
C
C        CONSTRUCTION DU TMS NSEF DU MATERIAU
         NBEFMA = MCN( MNNBEF + NMA )
         CALL LXTNDC( NTLXMA, 'NSEF', 'MOTS', WUSOEF+NBEFMA*NBSOEF )
         CALL LXTSOU( NTLXMA, 'NSEF', NTNSEM, MNNSEM )
C
C        variable NUTYOB 'numero de type de l''objet' entier
C        ( 1 : 'point' , 2 : 'ligne' , 3 : 'surface' , 4 : 'volume' )  ;
         MCN( MNNSEM + WUTYOB ) = 4
C
C        variable NUTFMA 'Ligne ou Surface fermee ou non-ferme V=> inconnu
C        ( -1 : 'inconnu' , 0 : 'non-ferme' , 1 : 'ligne ou surface fermee' )
         MCN( MNNSEM + WUTFMA ) = -1
C
C        variable NBSOEF 'nombre de sommets par EF' entier
C        ( 1 : 'sommet' , 2 : 'arete' , 4 : 'face' , 8 : 'cube' ) ;
         MCN( MNNSEM + WBSOEF ) = 8
C
C        variable NBTGEF 'nombre de tangentes par EF' entier
C        ( 0 : '0 tg par EF' , 1 : '1 tg par sommet' , 2 : '2 tg par arete' ,
C          8 : '8 tg par face' ,  24 : '24 tg par cube' )  ;
         MCN( MNNSEM + WBTGEF ) = 0
         MCN( MNNSEM + WBEFAP ) = 0
         MCN( MNNSEM + WBEFTG ) = 0
C
C        variable NBEFOB 'nombre des EF du PLSV'  entier ;
         MCN( MNNSEM + WBEFOB ) = NBEFMA
C
C        variable NUTYMA 'numero de type du maillage' entier
C        ( 0 : 'non structure' , ... )
         MCN( MNNSEM + WUTYMA ) = 0
C
C        RECENSEMENT DES EF DE CE MATERIAU NOMATE
         CALL NMOBNU( 'VOLUME', NOMATE, KERR(1) )
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,*) NBEFMA,' ELEMENTS FINIS pour le VOLUME ',
     %                      KERR(1)(1:24)
         ELSE
            WRITE(IMPRIM,*) NBEFMA,' FINITE ELEMENTS of the VOLUME ',
     %                      KERR(1)(1:24)
         ENDIF
C
         MNEFVL = MNNSEF + WUSOEF - 1
         MNEFMA = MNNSEM + WUSOEF - 1
         DO 50 NEF = 1, NBEFVL
C
            IF(  MCN( MNUDMEF - 1 + NEF ) .EQ. NOMATE ) THEN
C              REMPLISSAGE DU TABLEAU NSEF DU MATERIAU
               DO 40 K=1,NBSOEF
                  MCN( MNEFMA + K ) = MCN( MNEFVL + K )
 40            CONTINUE
               MNEFMA = MNEFMA + NBSOEF
            ENDIF
C
C           PASSAGE A L'EF SUIVANT DU VOLUME
            MNEFVL = MNEFVL + NBSOEF
 50      CONTINUE
C
C        RECHERCHE DU NOMBRE ET DES NUMEROS DES SOMMETS DE CE MATERIAU
         CALL AZEROI( NBSTVL+1, MCN(MNLIPO) )
C        PARCOURS DES SOMMETS DES NSEF DE CE MATERIAU
         MNS = MNNSEM + WUSOEF - 1
         DO 58 NEF=1,NBEFMA
            DO 55 K=1,NBSOEF
C              LE NUMERO DU SOMMET DANS LE VOLUME MULTI-MATERIAUX
               NS = MCN( MNS + K )
               IF( NS .LE. 0 ) GOTO 56
C              SOMMET NON NUL EST MARQUE
               MCN( MNLIPO + NS ) = 1
 55         CONTINUE
 56         MNS = MNS + NBSOEF
 58      CONTINUE
C
C        LES SOMMETS SONT RENUMEROTES DE 1 A NBS
         NBS = 0
         DO 60 K=1,NBSTVL
            IF( MCN(MNLIPO+K) .GT. 0 ) THEN
               NBS = NBS + 1
               MCN(MNLIPO+K) = NBS
            ENDIF
 60      CONTINUE
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,*) NBS,' SOMMETS pour le VOLUME ',KERR(1)(1:24)
         ELSE
            WRITE(IMPRIM,*) NBS,' VERTICES of the VOLUME ',KERR(1)(1:24)
         ENDIF
C
C        PARCOURS DES SOMMETS DES NSEF DE CE MATERIAU POUR
C        OBTENIR UNE NUMEROTATION LOCALE DES SOMMETS
         MNS = MNNSEM + WUSOEF - 1
         DO 64 NEF=1,NBEFMA
            DO 62 K=1,NBSOEF
C              LE NUMERO DU SOMMET DANS LE VOLUME PARTITION
               NS = MCN( MNS + K )
               IF( NS .LE. 0 ) GOTO 63
C              LE NOUVEAU NUMERO LOCAL
               MCN( MNS + K ) = MCN( MNLIPO + NS )
 62         CONTINUE
 63         MNS = MNS + NBSOEF
 64      CONTINUE
C        LA DATE DE CREATION
         CALL ECDATE( MCN(MNNSEM) )
C        AJOUT DU NUMERO DU TABLEAU DESCRIPTEUR
         MCN( MNNSEM + MOTVAR(6) ) = NONMTD( '~>>>NSEF' )
C
C        DECLARATION OUVERTURE DU TABLEAU 'XYZSOMMET' DU MATERIAU NOMATE
         MOTS = WYZSOM + 3 * NBS
         CALL LXTNDC( NTLXMA, 'XYZSOMMET', 'ENTIER', MOTS )
         CALL LXTSOU( NTLXMA, 'XYZSOMMET',  NTXYZM , MNXYZM )
C        LE NOMBRE DE SOMMETS
         MCN( MNXYZM + WNBSOM ) = NBS
C        LE NOMBRE DE TANGENTES
         MCN( MNXYZM + WNBTGS ) = 0
C        LE NOMBRE DE COORDONNEES PAR SOMMET
         MCN( MNXYZM + WBCOOR ) = 3
C        L'ADRESSE DES COORDONNEES DANS XYZSOMMET DU MATERIAU
         MNSC = MNXYZM + WYZSOM
C        L'ADRESSE DES COORDONNEES DANS XYZSOMMET DU VOLUME
         MNSP = MNXYZS + WYZSOM
C        LES COORDONNEES
         DO 70 K = 1, NBSTVL
            NS = MCN( MNLIPO + K )
            IF( NS .GT. 0 ) THEN
C              LES 3 COORDONNEES DU SOMMET DU MATERIAU
               RMCN( MNSC   ) = RMCN( MNSP   )
               RMCN( MNSC+1 ) = RMCN( MNSP+1 )
               RMCN( MNSC+2 ) = RMCN( MNSP+2 )
               MNSC = MNSC + 3
            ENDIF
            MNSP = MNSP + 3
 70      CONTINUE
C        LA DATE DE CREATION
         CALL ECDATE( MCN(MNXYZM) )
C        AJOUT DU NUMERO DU TABLEAU DESCRIPTEUR
         MCN( MNXYZM + MOTVAR(6) ) =  NONMTD( '~>>>XYZSOMMET' )
C
 100  CONTINUE
C
 9990 IF( MNLIPO .GT. 0 ) CALL TNMCDS( 'ENTIER', NBSTVL+1, MNLIPO )
      IF( MNNUMA .GT. 0 ) CALL TNMCDS( 'ENTIER', 2*NBDM, MNNUMA )
      RETURN
      END
