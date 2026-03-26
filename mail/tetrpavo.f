      SUBROUTINE TETRPAVO( NUVOLU, NTLXVL,  LADEFI, NBVOPA,
     %                     NBSOMM, PTXYZD,  MNSOFR, MNLIPO,
     %                     NUDTETR,NBTETR0, NOTETR, NF1VO, NVOLTE,
     %                     MXFACO, MNFACO,
     %                     NBSP,   NBTETT,  VOLUMT,
     %                     NTCUPA, MNCUPA,  NTSOMM, MNSOMM)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : REGENERATION DES VRAIES COORDONNEES DES POINTS DEPLACES
C ----- GENERATION DU TABLEAU 'MATERIAUX' DE LA PARTITION DE VOLUMES
C       DISTRUBUER LA TETRAEDRISATION     DE LA PARTITION DES VOLUMES
C       GENERER LES SURFACES FRONTIERES ET INTERFACES

C ENTREES:
C --------
C NUVOLU : NUMERO DU VOLUME PARTITION DANS LE LEXIQUE DES VOLUMES
C NTLXVL : NUMERO DU TABLEAU TS DU LEXIQUE DU VOLUME PARTITION
C LADEFI : TABLEAU DE DEFINITION DU VOLUME PARTITIONNE
C          CF '~/td/d/a_volume__definition'
C NBVOPA : NOMBRE DE VOLUMES DE LA PARTITION
C NBSOMM : NOMBRE ACTUEL DE POINTS DECLARES DANS PTXYZD
C PTXYZD : PAR POINT : X  Y  Z  DISTANCE_SOUHAITEE
C MNSOFR : ADRESSE MCN DU TABLEAU NPSOFR
C          1 000 000  + NUMERO DU POINT INTERNE DES POINTS AJOUTES
C         (1 000 000) * (NO SURFACE + 1) + NUMERO DU SOMMET DANS
C          LA NUMEROTATION DES SOMMETS DE LA SURFACE
C          -4 SI SOMMET DE LA GRILLE DES TETRAEDRES
C          -1 SI SOMMET DE LA GRILLE DES TETRAEDRES ELIMINE
C MNLIPO : TABLEAU(0:MXSOMM) D'ENTIERS UTILE POUR DES TRIS
C NUDTETR: NUMERO DU DERNIER TETRAEDRE ACTIF  DANS NOTETR
C NBTETR0: NOMBRE TOTAL DE TETRAEDRES STOCKES DANS NOTETR
C NOTETR : LISTE DES TETRAEDRES OCCUPES POINTEE PAR N1TETS
C                               VIDES   POINTEE PAR N1TEVI
C          1:4 SOMMET1   < SOMMET2    SOMMET3     SOMMET4
C          ORDRE CLASSIQUE 1 2 3 DIRECT EN BAS ET 4 EN HAUT
C          => VOLUME POSITIF SINON TETRAEDRE DEGENERE
C          SUITE INUTILE DANS CE SOUS PROGRAMME
C          5:8 TETRAEDRE1  TETREDRE2  TETRAEDRE3  TETREDRE4  VOISINS PAR
C              FACE1       FACE2D      FACE3       FACE4
C NF1VO  : NF1VO(1,I) = NUMERO DU 1-ER TETRAEDRE DU VOLUME I
C          NF1VO(3,I) = NUMERO DU VOLUME I de 1 a NBVOPA
C          NF1VO TABLEAU ENTIER( 1:3 , 1:NBVOPA )
C NVOLTE : NVOLTE(I) = NUMERO DU VOLUME DU TETRAEDRE DE 1 A NBVOPA
C MXFACO : NOMBRE MAXIMAL DE FACES DECLARABLES DANS LEFACO
C MNFACO : ADRESSE MCN DU TABLEAU LEFACO
C          FACES DU CONTOUR OU INTERFACES ENTRE VOLUMES
C          IL CONTIENT DANS CET ORDRE
C          1:   =0 POUR UNE FACE VIDE
C          123: NO (DANS PTXYZD) DU SOMMET 1, SOMMET 2, SOMMET 3
C          45:  NO (DANS NUVOPA 0 SINON) DU VOLUME1 , VOLUME2 DE LA FACE
C          678: NO (DANS LEFACO) DE LA FACE ADJACENTE PAR L'ARETE 1 2 3
C          9: ATTENTION: UNE ARETE PEUT APPARTENIR A PLUS DE 2 FACES
C             => CHAINAGE CIRCULAIRE DE CES FACES DANS LEFACO
C             LEFACO(9,*) -> FACE SUIVANTE (*=0:VIDE, *<>0:NON VIDE)
C          10: HACHAGE AVEC LA SOMME DES 3 SOMMETS MODULO MXFACO
C              LF = MOD( NOSOFA(1)+NOSOFA(2)+NOSOFA(3) , MXFACO ) + 1
C              NF = LEFACO( 10, LF ) LE NUMERO DE LA 1-ERE FACE DANS LEFACO
C              SI LA FACE NE CONVIENT PAS. PASSAGE A LA SUIVANTE
C              NF = LEFACO( 9, NF )  ...
C          11: >0  NO NOTETR D'UN TETRAEDRE AYANT CETTE FACE,
C              =0  SINON
C SORTIES:
C --------
C NBTETT : NOMBRE TOTAL DE TETRAEDRES DE LA PARTITION DE VOLUMES
C NBSP   : NOMBRE DE SOMMETS DES VOLUMES DE LA PARTITION
C VOLUMT : VOLUME TOTAL DES TETRAEDRES DE LA PARTITION DE VOLUMES
C NTCUPA : NUMERO      DU TS DES NO DES SOMMETS DES CUBES DU VOLUME PARTITION
C MNCUPA : ADRESSE MCN DU TS DES NO DES SOMMETS DES CUBES DU VOLUME PARTITION
C          CF '~/td/d/a___nsef'
C NTSOMM : NUMERO DU TABLEAU TS DES XYZ DES SOMMETS DE LA PARTITION
C MNSOMM : ADRESSE MCN DU    TS DES XYZ DES SOMMETS DE LA PARTITION
C          CF '~/td/d/a___xyzsommet'
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE DU PERRAY     Mars 2016
C MODIFS : ALAIN PERRONNET Saint PIERRE du PERRAY         Septembre 2018
C2345X7..............................................................012
      PARAMETER        (NBCNSU=12, MXNMSU=5)
      include"./incl/ntmnlt.inc"
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/a___materiaux.inc"
      include"./incl/a_surface__definition.inc"
      include"./incl/a_volume__definition.inc"
      include"./incl/a___nsef.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/xyzext.inc"
      include"./incl/trvari.inc"

      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)
      COMMON / TRTETR / STOPTE,TRACTE
      LOGICAL           STOPTE,TRACTE, TRACTE0

      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1),DMCN(1))

      DOUBLE PRECISION  PTXYZD(4,*)
      INTEGER           LADEFI(0:*), NOTETR(8,*),
     %                  NF1VO(3,NBVOPA), NVOLTE(1:*)

      INTEGER           NOSOCU(1:8), NOTETRA(5)
      CHARACTER*24      NMVOLU, NMSU, NMVO, NMNO
      DOUBLE PRECISION  V, VV, VOLTET, VOLUMT
      CHARACTER*80      KTITRE

      WRITE(IMPRIM,*)
      WRITE(IMPRIM,*)'tetrpavo: GENERATION des TMS xyzsommet et nsef des
     % VOLUMES de la PARTITION'


C     REGENERATION DES VRAIES COORDONNEES DES POINTS DEPLACES
C     -------------------------------------------------------
      MNSOF1 = MNSOFR - 1
      DO N=1,NBSOMM

C        SI MCN(MNSOF1+N) EST NEGATIF LES COORDONNEES DU POINT
C        ONT ETE MODIFIEES DANS LA TETRAEDRISATION ETOILEE
         NP = MCN( MNSOF1 + N )
         IF( NP .LT. 0 ) THEN

C           SES COORDONNEES ONT ETE MODIFIEES
C           RETOUR AUX ANCIENNES COORDONNEES
            NP = ABS( NP )
            IF( NP .GT. 1 000 000 .AND. NP .LT. 2 000 000 ) THEN

C              POINT INITIAL DE NUMERO NP1:( 1 000 000 + NP1 )
               NP = NP - 1 000 000
C              RECHERCHE DE SES COORDONNEES
               CALL LXNLOU( NTPOIN, NP       , NTPO, MNPO )
               CALL LXTSOU( NTPO, 'XYZSOMMET', NTSO, MNSO )
               MNSO = MNSO + WYZSOM + 3 * NP - 3
               PTXYZD(1,N) = RMCN( MNSO   )
               PTXYZD(2,N) = RMCN( MNSO+1 )
               PTXYZD(3,N) = RMCN( MNSO+2 )

            ELSE IF( NP .GT. 2 000 000 ) THEN

C              SOMMET D'UNE SURFACE INITIALE:(1 000 000*(NOSU+1) + NP1 )
               NOSU = NP / 1 000 000 - 1
               NP   = MOD( NP, 1 000 000 )
C              LE LEXIQUE DE LA SURFACE
               CALL LXNLOU( NTSURF, NOSU, NTLXSU, MNPO )
C              LE TABLEAU 'XYZSOMMET' DE CETTE SURFACE
               CALL LXTSOU( NTLXSU, 'XYZSOMMET', NTSOSU, MNSOSU )
C              LES 3 COORDONNEES DE CE SOMMET
               MNSO = MNSOSU + WYZSOM + 3 * NP - 3
               PTXYZD(1,N) = RMCN( MNSO   )
               PTXYZD(2,N) = RMCN( MNSO+1 )
               PTXYZD(3,N) = RMCN( MNSO+2 )

            ENDIF

         ENDIF
         MNS = MNS + 4
      ENDDO

C     SI AU MOINS 2 VOLUMES DANS LA PARTITION ALORS
C     GENERATION DU TABLEAU 'MATERIAUX' DE LA PARTITION DE VOLUMES
C     ============================================================
C     TABLEAU 'MATERIAUX' DE LA PARTITION TOTALE
C     LA PARTITION EST FORMEE DE NBVOPA VOLUMES
      NBVOPA = LADEFI( WBVOPA )

      IF( NBVOPA .LE. 1 ) THEN
         NTMATE = 0
         MNMATE = 0
         MNUDMEF= 0
         GOTO 10
      ENDIF

C     CONSTRUCTION DE 'MATERIAUX' CAR PLUSIEURS MATERIAUX
      CALL LXTSOU( NTLXVL, 'MATERIAUX', NTMATE, MNMATE )
      IF( NTMATE .GT. 0 ) CALL LXTSDS( NTLXVL, 'MATERIAUX' )
      MOTS2M = WUDMEF + NBTETR0
      CALL LXTNDC( NTLXVL, 'MATERIAUX', 'ENTIER', MOTS2M )
      CALL LXTSOU( NTLXVL, 'MATERIAUX',  NTMATE , MNMATE )

C     LES PREMIERES VARIABLES DU TABLEAU 'MATERIAUX'
C     LE NUMERO DU TABLEAU DESCRIPTEUR DE CE TMS
C     AJOUT DU NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNMATE + MOTVAR(6) ) = NONMTD( '~>>>MATERIAUX' )

C     NOMBRE DE MATERIAUX
      MCN( MNMATE + WNBDM ) = NBVOPA

C     NOMBRE D'ELEMENTS FINIS DU MAILLAGE
      MCN( MNMATE + WBDMEF ) = NBTETR0

C     ADRESSE MCN DE DEBUT DU TABLEAU NUDMEF(1..NBDMEF)
C     'NUMERO DU MATERIAU DE CHAQUE EF'
      MNUDMEF = MNMATE + WUDMEF


C     GENERATION DU TABLEAU 'NSEF' DES NBVOPA VOLUMES ET DE LA PARTITION
C     ==================================================================
 10   CALL LXTSOU( NTLXVL, 'NSEF', NTCUPA, MNCUPA )
      IF( NTCUPA .GT. 0 ) THEN
         CALL NMOBNU( 'VOLUME', NUVOLU, NMVOLU )
         CALL MAILDS( 'VOLUME', NMVOLU )
      ENDIF
      NBSOEF = 8
      MOTS   = WUSOEF + NBSOEF * NBTETR0
      CALL LXTNDC( NTLXVL, 'NSEF', 'ENTIER', MOTS )
      CALL LXTSOU( NTLXVL, 'NSEF',  NTCUPA , MNCUPA )
C
C     LES PREMIERES VARIABLES DU TABLEAU 'NSEF'
C
C     LE NUMERO DU TABLEAU DESCRIPTEUR DE CE TMS
      NOTADS = NONMTD( '~>>>NSEF' )
      NOTADC = NONMTD( '~>>>XYZSOMMET' )
C
C     variable NUTYOB        'numero de type de l''objet'    entier
C     ( 1 : 'point' , 2 : 'ligne' , 3 : 'surface' , 4 : 'volume' )  ;
      MCN( MNCUPA + WUTYOB ) = 4
C
C     variable NUTFMA 'Ligne ou Surface fermee ou non-ferme (P et V sont inconnu
C     ( -1 : 'inconnu' , 0 : 'non-ferme' , 1 : 'ligne ou surface fermee' ) ;
      MCN( MNCUPA + WUTFMA ) = -1
C
C     variable NBSOEF 'nombre de sommets par EF' entier
C     ( 1 : 'sommet' , 2 : 'arete' , 4 : 'face' , 8 : 'cube' ) ;
      MCN( MNCUPA + WBSOEF ) = 8
C
C     variable NBTGEF 'nombre de tangentes par EF' entier
C     ( 0 : '0 tg par EF' , 1 : '1 tg par sommet' , 2 : '2 tg par arete' ,
C       8 : '8 tg par face' ,  24 : '24 tg par cube' )  ;
      MCN( MNCUPA + WBTGEF ) = 0
      MCN( MNCUPA + WBEFAP ) = 0
      MCN( MNCUPA + WBEFTG ) = 0
C
C     variable NBEFOB 'nombre des EF du PLSV'  entier ;
      MCN( MNCUPA + WBEFOB ) = NBTETR0
C
C     variable NUTYMA 'numero de type du maillage' entier
C     ( 0 : 'non structure' , ... )
      MCN( MNCUPA + WUTYMA ) = 0
C
C     TABLEAU  NUSOEF(1..NBSOEF,1..NBEFOB)
C    'NUMERO DES NBSOEF SOMMETS DES NBEFOB NSEF DE L''OBJET'
C
C     BOUCLE SUR LES TETRAEDRES NSEF DES VOLUMES DE LA PARTITION
C     ----------------------------------------------------------
      NBTETT = 0
      VOLUMT = 0D0
C     LES SOMMETS 5 A 8 DU CUBE SONT NULS
      NOSOCU(5) = 0
      NOSOCU(6) = 0
      NOSOCU(7) = 0
      NOSOCU(8) = 0
      MNTEVP    = MNCUPA + WUSOEF - 1

      DO 50 NVP = 1, NBVOPA

C        LE VOLUME DU VOLUME NVP
         VV = 0D0

cccC        LE NUMERO DU 1-ER TETRAEDRE DU VOLUME NVP DE LA PARTITION
ccc         NT = NF1VO( 1, NVP )
ccc         IF( NT .LE. 0 ) GOTO 50

C        LE NO DE VOLUME DE 1 A NBVOPA
         NOVOLU = NF1VO( 3, NVP )

C        LE NUMERO DU VOLUME DANS LE LEXIQUE DES VOLUMES
         NOVOLX = LADEFI( WUVOPA - 1 + NOVOLU )

C        BOUCLE SUR LES TETRAEDRES DU VOLUME NOVOLU
         MNT0 = MNTEVP
         NBTT = 0

         DO NT = 1, NUDTETR

            IF( NOTETR(1,NT) .GT. 0 .AND. NVOLTE(NT) .EQ. NOVOLU ) THEN

C              LE VOLUME PARTITION
C              LE NUMERO EVENTUEL DE VOLUME DU TETRAEDRE DANS 'MATERIAUX'
               IF( MNUDMEF .GT. 0 ) THEN
                  MCN( MNUDMEF + NBTT ) = NOVOLX
               ENDIF
               NBTT = NBTT + 1

C              LE NUMERO DES 4 SOMMETS DU TETRAEDRE NT DE NOTETR
               DO I=1,4
                  NOSOCU(I) = NOTETR(I,NT)
               ENDDO
               DO I=1,4
                  IF( NOSOCU(I) .LE. 0 .OR. NOSOCU(I) .GT. NBSOMM ) THEN
C                    DETECTION D'UNE ERREUR
                     IF( LANGAG .EQ. 0 ) THEN
                        WRITE(IMPRIM,19321) NT,(NOSOCU(J),J=1,4),NBSOMM
                     ELSE
                        WRITE(IMPRIM,29321) NT,(NOSOCU(J),J=1,4),NBSOMM
                     ENDIF
                  ENDIF
               ENDDO
19321          FORMAT(' TETRAEDRE ',I9,' NO SOMMETS =',4I9,' >',I9,
     %         '=NBSOMM => INCORRECT.',
     %         ' UTILISER ENSUITE AMELIORATION d''UNE TETRAEDRISATION')
29321          FORMAT(' TETRAHEDRON ',I9,' VERTICES=',4I9,
     %         ' >',I9,'=NBSOMM => INCORRECT.',
     %         ' TETRAHEDRIZATION IMPROVEMENT is RECOMMANDED AFTER')

C              LE VOLUME EST-IL POSITIF ?
               V = VOLTET( PTXYZD( 1, NOSOCU(1) ),
     %                     PTXYZD( 1, NOSOCU(2) ),
     %                     PTXYZD( 1, NOSOCU(3) ),
     %                     PTXYZD( 1, NOSOCU(4) ) )
               IF( V .LE. 0.D0 ) THEN
                  IF( LANGAG .EQ. 0 ) THEN
                     WRITE(IMPRIM,19320) NT,(NOSOCU(J),J=1,4),V
                  ELSE
                     WRITE(IMPRIM,29320) NT,(NOSOCU(J),J=1,4),V
                  ENDIF
19320 FORMAT(' TETRAEDRE ',I9,' NO SOMMETS =',4I9,' de VOLUME ',
     %          G25.16,' =<0 INCORRECT.',
     %        ' UTILISER ENSUITE AMELIORATION d''UNE TETRAEDRISATION')
29320 FORMAT(' TETRAHEDRON ',I9,' VERTICES=',4I9,' HAS a VOLUME ',
     %          G25.16,' =<0 INCORRECT.',
     %        ' EXECUTE AFTER the TETRAHEDRIZATION IMPROVEMENT')

cccC                 SINON PERMUTATION DE 2 SOMMETS
CCC                  J         = NOSOCU(3)
CCC                  NOSOCU(3) = NOSOCU(4)
CCC                  NOSOCU(4) = J

C                 TRACE DES FACES TRIANGULAIRES DES NBTETRA TETRAEDRES
C                 + FACES TRIANGULAIRES DES TETRAEDRES OPPOSES
                  NBTETRA = 1
                  NOTETRA(1) = NT
                  DO M = 1, 4
                     NTOP = NOTETR(4+M,NT)
                     IF( NTOP .GT. 0 .AND. NOTETR(1,NTOP) .GT. 0 ) THEN
                        NBTETRA = NBTETRA + 1
                        NOTETRA(NBTETRA) = NTOP
                     ENDIF
                  ENDDO
                  KTITRE='tetrpavo: TETRAEDRE de VOLUME<=0'
                  TRACTE0 = TRACTE
                  CALL TRFETO13(KTITRE, PTXYZD, NBTETRA,NOTETRA, NOTETR)
                  TRACTE = TRACTE0
               ENDIF

C              VOLUME TOTAL
               VV = VV + ABS( V )

C              LE VOLUME DU TETRAEDRE EST POSITIF
C              L'INTEGRATION DE CE CUBE DANS LE TABLEAU 'NSEF'
               DO J=1,8
                  MCN( MNTEVP + J ) = NOSOCU( J )
               ENDDO
               MNTEVP = MNTEVP + 8

            ENDIF

         ENDDO

C        LE LEXIQUE DU VOLUME NOVOLX
C        ---------------------------
         CALL TAMSOU( NTVOLU, MNVOLU )
         NTLXVO = MCN( MNVOLU+MCN(MNVOLU)*NOVOLX+MCN(MNVOLU+2)+2 )
         CALL LXTSOU( NTLXVO, 'NSEF', NTCUVO, MNCUVO )
         IF( NTCUVO .GT. 0 ) THEN
C           LES TABLEAUX 'NSEF' 'XYZSOMMET' EXISTANT SONT DETRUITS
            CALL NMOBNU( 'VOLUME', NOVOLX, NMVO )
            CALL MAILDS( 'VOLUME', NMVO )
         ENDIF

C        DECLARATION OUVERTURE DU TABLEAU 'NSEF' DE CE VOLUME NOVOLX
C        -----------------------------------------------------------
C        LE NOMBRE DE TETRAEDRES DE CE VOLUME
         NBT  = ( MNTEVP - MNT0 ) / 8
         MOTS = WUSOEF + NBSOEF * NBT
         CALL LXTNDC( NTLXVO, 'NSEF', 'ENTIER', MOTS )
         CALL LXTSOU( NTLXVO, 'NSEF',  NTCUVO , MNCUVO )

C        COPIE DE LA LISTE DES NUMEROS DES CUBES
         CALL TRTATA( RMCN(MNT0+1), RMCN(MNCUVO+WUSOEF), 8*NBT )

C        LES PREMIERES VARIABLES DU TABLEAU 'NSEF'
C        variable NUTYOB        'numero de type de l''objet'    entier
C        ( 1 : 'point' , 2 : 'ligne' , 3 : 'surface' , 4 : 'volume' )  ;
         MCN( MNCUVO + WUTYOB ) = 4

C        variable NUTFMA 'Ligne ou Surface fermee ou non-ferme (P V inconnus)
C        ( -1 : 'inconnu' , 0 : 'non-ferme' , 1 : 'ligne ou surface fermee' ) ;
         MCN( MNCUVO + WUTFMA) = -1

C        variable NBSOEF 'nombre de sommets par EF' entier
C        ( 1 : 'sommet' , 2 : 'arete' , 4 : 'face' , 8 : 'cube' ) ;
         MCN( MNCUVO + WBSOEF) = 8

C        variable NBTGEF 'nombre de tangentes par EF' entier
C        ( 0 : '0 tg par EF' , 1 : '1 tg par sommet' , 2 : '2 tg par arete' ,
C          8 : '8 tg par face' ,  24 : '24 tg par cube' )  ;
         MCN( MNCUVO + WBTGEF ) = 0
         MCN( MNCUVO + WBEFAP ) = 0
         MCN( MNCUVO + WBEFTG ) = 0

C        variable NBEFOB 'nombre des EF du PLSV'  entier ;
         MCN( MNCUVO + WBEFOB ) = NBT

C        variable NUTYMA 'numero de type du maillage' entier
C        ( 0 : 'non structure' , ... )
         MCN( MNCUVO + WUTYMA ) = 0

C        RECHERCHE DU NOMBRE ET DES NUMEROS DES SOMMETS DE CE VOLUME
         CALL AZEROI( NBSOMM+1, MCN(MNLIPO) )
C        PARCOURS DES SOMMETS DES NSEF DE CE VOLUME
         MNS = MNCUVO + WUSOEF - 1
         DO J=1,NBT
            DO I=1,4
C              LE NUMERO DU SOMMET DANS LE VOLUME PARTITION
               NP = MCN( MNS + I )
C              IL EST MARQUE
               MCN( MNLIPO + NP ) = 1
            ENDDO
            MNS = MNS + 8
         ENDDO

C        LES SOMMETS SONT RENUMEROTES DE 1 A NBS
         NBS = 0
         DO I=1,NBSOMM
            IF( MCN(MNLIPO+I) .GT. 0 ) THEN
               NBS = NBS + 1
               MCN(MNLIPO+I) = NBS
            ENDIF
         ENDDO

C        PARCOURS DES SOMMETS DES NSEF DE CE VOLUME  POUR
C        OBTENIR UNE NUMEROTATION LOCALE DES SOMMETS
         MNS = MNCUVO + WUSOEF - 1
         DO J=1,NBT
            DO I=1,4
C              LE NUMERO DU SOMMET DANS LE VOLUME PARTITION
               NP = MCN( MNS + I )
C              LE NOUVEAU NUMERO LOCAL
               MCN( MNS + I ) = MCN( MNLIPO + NP )
            ENDDO
            MNS = MNS + 8
         ENDDO
C        LA DATE DE CREATION
         CALL ECDATE( MCN(MNCUVO) )
C        AJOUT DU NUMERO DU TABLEAU DESCRIPTEUR
         MCN( MNCUVO + MOTVAR(6) ) = NOTADS

C        DECLARATION OUVERTURE DU TABLEAU 'XYZSOMMET' DE CE VOLUME NOVOLX
C        ----------------------------------------------------------------
         MOTS = WYZSOM + 3 * NBS
         CALL LXTNDC( NTLXVO, 'XYZSOMMET', 'ENTIER', MOTS )
         CALL LXTSOU( NTLXVO, 'XYZSOMMET',  NTCUSO , MNCUSO )
C        LE NOMBRE DE SOMMETS
         MCN( MNCUSO + WNBSOM ) = NBS
C        LE NOMBRE DE TANGENTES
         MCN( MNCUSO + WNBTGS ) = 0
C        LE NOMBRE DE COORDONNEES PAR SOMMET
         MCN( MNCUSO + WBCOOR ) = 3
C        L'ADRESSE DES COORDONNEES DE CE VOLUME
         MNSC = MNCUSO + WYZSOM
C        LES COORDONNEES
         DO J=1,NBSOMM
            NP = MCN( MNLIPO + J )
            IF( NP .GT. 0 ) THEN
C              LES 3 COORDONNEES
               RMCN( MNSC   ) = REEL( PTXYZD(1,J) )
               RMCN( MNSC+1 ) = REEL( PTXYZD(2,J) )
               RMCN( MNSC+2 ) = REEL( PTXYZD(3,J) )
               MNSC = MNSC + 3
            ENDIF
         ENDDO
C        LA DATE DE CREATION
         CALL ECDATE( MCN(MNCUSO) )
C        AJOUT DU NUMERO DU TABLEAU DESCRIPTEUR
         MCN( MNCUSO + MOTVAR(6) ) = NOTADC

C        LE VOLUME GENERAL
         VOLUMT = VOLUMT + VV
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10049) NOVOLX,VV,NBTT,NBS
         ELSE
            WRITE(IMPRIM,20049) NOVOLX,VV,NBTT,NBS
         ENDIF
10049    FORMAT(' tetrpavo: VOLUME',I7,' de VOLUME=',G25.16,
     %' pour',I10,' TETRAEDRES et',I10,' SOMMETS')
20049    FORMAT(' tetrpavo: VOLUME',I7,' of VOLUME=',G25.16,
     %' for',I10,' TETRAHEDRA and',I10,' VERTICES')

         NBTETT = NBTETT + NBTT
 50   ENDDO


C     AFFICHAGE DU VOLUME TOTAL ET DU NOMBRE DE TETRAEDRES DE LA PARTITION
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10050) VOLUMT,NBTETT
      ELSE
         WRITE(IMPRIM,20050) VOLUMT,NBTETT
      ENDIF
10050 FORMAT(' tetrpavo: VOLUME du MAILLAGE TOTAL=',G25.16,
     %       ' pour un NOMBRE TOTAL de TETRAEDRES=',I10/)
20050 FORMAT(' tetrpavo: VOLUME  of   TOTAL MESH=',G25.16,
     %       ' for',I10,' TETRAHEDRA'/)

C     LA DATE DE CREATION DU TABLEAU 'NSEF' DU VOLUME PARTITION
      CALL ECDATE( MCN(MNCUPA) )
C     AJOUT DU NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNCUPA + MOTVAR(6) ) = NOTADS

      IF( MNUDMEF .GT. 0 ) THEN
C        LA DATE DE CREATION DU TABLEAU 'MATERIAUX' DU VOLUME PARTITION
         CALL ECDATE( MCN(MNMATE) )
C        AJOUT DU NUMERO DU TABLEAU DESCRIPTEUR
         MCN( MNMATE + MOTVAR(6) ) = NONMTD( '~>>>MATERIAUX' )
      ENDIF

C     DECLARATION OUVERTURE DU TABLEAU 'XYZSOMMET' DU VOLUME PARTITION
C     ================================================================
C     RECHERCHE DU NOMBRE ET DES NUMEROS DES SOMMETS DE CE VOLUME
      CALL AZEROI( NBSOMM+1, MCN(MNLIPO) )
C     PARCOURS DES SOMMETS DES NSEF DE CE VOLUME PARTITION
      MNS = MNCUPA + WUSOEF - 1
      DO J=1,NBTETT
         DO I=1,4
C           LE NUMERO DU SOMMET DANS LE VOLUME PARTITION
            NP = MCN( MNS + I )
C           IL EST MARQUE
            MCN( MNLIPO + NP ) = 1
         ENDDO
         MNS = MNS + 8
      ENDDO

C     LES SOMMETS DE LA PARTITION SONT RENUMEROTES DE 1 A NBS
      NBS = 0
      DO I=1,NBSOMM
         IF( MCN(MNLIPO+I) .GT. 0 ) THEN
            NBS = NBS + 1
            MCN(MNLIPO+I) = NBS
         ENDIF
      ENDDO
      NBSP = NBS

C     PARCOURS DES SOMMETS DES NSEF DE CE VOLUME POUR
C     OBTENIR UNE NUMEROTATION LOCALE DES SOMMETS
      MNS0 = MNCUPA + WUSOEF - 1
      MNS  = MNS0
      DO 55 J=1,NBTETT
         DO I=1,4
C           LE NUMERO DU SOMMET DANS LE VOLUME PARTITION
            NP = MCN( MNS + I )
C           LE NOUVEAU NUMERO LOCAL
            NEWP = MCN( MNLIPO + NP )
            IF( NEWP .LE. 0 ) THEN
               PRINT*,'tetrpavo: ANOMALIE EF',J,' SOMMET',NP,
     %                ' DEVIENT',NEWP
               GOTO 55
            ENDIF
            MCN( MNS + I ) = NEWP
C           LE NUMERO NUL DES SOMMETS AU DELA DE 4
            MCN( MNS + I + 4 ) = 0
         ENDDO
         MNS = MNS + 8
 55   ENDDO

C     LE NOMBRE TOTAL DE TETRAEDRES
      NBTETT = ( MNS - MNS0 ) / 8
      MCN( MNCUPA + WBEFOB ) = NBTETT

      PRINT*
      IF( LANGAG .EQ. 0 ) THEN
         PRINT*,'tetrpavo: NOMBRE TOTAL de SOMMETS      =',NBSP
         PRINT*,'tetrpavo: NOMBRE TOTAL d''ELEMENTS FINIS=',NBTETT
      ELSE
       PRINT*,'tetrpavo: VERTICES        TOTAL NUMBER=',NBSP
       PRINT*,'tetrpavo: FINITE ELEMENTS TOTAL NUMBER=',NBTETT
      ENDIF

      MOTS = WYZSOM + 3 * NBSP
      CALL LXTNDC( NTLXVL, 'XYZSOMMET', 'ENTIER', MOTS )
      CALL LXTSOU( NTLXVL, 'XYZSOMMET',  NTSOMM, MNSOMM )
C     LE NOMBRE DE SOMMETS DE LA PARTITION
      MCN( MNSOMM + WNBSOM ) = NBSP
C     LE NOMBRE DE TANGENTES
      MCN( MNSOMM + WNBTGS ) = 0
C     LE NOMBRE DE COORDONNEES PAR SOMMET
      MCN( MNSOMM + WBCOOR ) = 3
C     L'ADRESSE DES COORDONNEES DE CE VOLUME
      MNSC = MNSOMM + WYZSOM
C     LES COORDONNEES
      DO J=1,NBSOMM
         NP = MCN( MNLIPO + J )
         IF( NP .GT. 0 ) THEN
C           LES 3 COORDONNEES
            RMCN( MNSC   ) = REEL( PTXYZD(1,J) )
            RMCN( MNSC+1 ) = REEL( PTXYZD(2,J) )
            RMCN( MNSC+2 ) = REEL( PTXYZD(3,J) )
            MNSC = MNSC + 3
         ENDIF
      ENDDO
C     LA DATE DE CREATION
      CALL ECDATE( MCN(MNSOMM) )
C     AJOUT DU NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNSOMM + MOTVAR(6) ) = NOTADC

C     CREATION DES SURFACES INTERFACES ET FRONTIERES DE LA PARTITION
C     ==============================================================
C     DECLARATION DU TABLEAU FRIN POINTEUR SUR LA 1-ERE FACE
C     DE CHAQUE COUPLE ( NO VOLUME1, NO VOLUME2 )
C     AVEC NO VOLUME2 = 0 POSSIBLE POUR UNE FRONTIERE
      NBVOP1 = NBVOPA + 1
      MXFRIN = NBVOP1 * NBVOP1
      MNFRIN = -1
      CALL TNMCDC( 'ENTIER', MXFRIN, MNFRIN )
      IF( MNFRIN .LE. 0 ) THEN
         NBLGRC(NRERR) = 2
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1)='tetrpavo: SATURATION MEMOIRE du TABLEAU MNFRIN'
            KERR(2)='AUGMENTEZ la TAILLE du SUPER-TABLEAU MCN'
         ELSE
            KERR(1)='tetrpavo: ALL MEMORY is USED by the TMS'
            KERR(2)='AUGMENT the LENGTH of SUPER-ARRAY MCN'
         ENDIF
         CALL LEREUR
         IERR = 1
         GOTO 9999
      ENDIF

C     INITIALISATION A ZERO
      CALL AZEROI( MXFRIN, MCN(MNFRIN) )

C     BOUCLE SUR LES FACES DE FACO ( INTERFACES OU FRONTALIERES )
      MNF = MNFACO + 11
      DO I=1,MXFACO
C        LA FACE EXISTE-T-ELLE DANS LEFACO?
         IF( MCN(MNF) .GT. 0 ) THEN
C           LA FACE EST ACTIVE DANS LEFACO
C           ELLE EST CHAINEE DANS FRIN
            NV1   = MCN( MNF + 3 )
            NV2   = MCN( MNF + 4 )
C           LA FACE DEVIENT LA PREMIERE DU CHAINAGE (NV1,NV2)
            MNFI  = MNFRIN + NV1 + NBVOP1 * NV2
C           LA FACE SUIVANT I EST L'ANCIENNE PREMIERE DU CHAINAGE
            MCN(MNF+8) = MCN(MNFI)
C           LA FACE I DANS LEFACO EST LA NOUVELLE PREMIERE DU CHAINAGE
            MCN(MNFI) = I
         ENDIF
C        PASSAGE A LA FACE SUIVANTE
         MNF = MNF + 11
      ENDDO

C     LE 4-EME SOMMET DE CHAQUE FACE EST NUL
      NOSOCU(4) = 0

C     BOUCLE SUR LES COUPLES ( NO VOLUME1, 0 OU NO VOLUME 2 )
C     -------------------------------------------------------
      DO NV1 = 1, NBVOPA

C        LE NOM DU VOLUME NV1 DANS LE LEXIQUE DES VOLUMES
         CALL NMOBNU( 'VOLUME', LADEFI(WUVOPA-1+NV1), NMVO )
C        RECHERCHE DU 1-ER CARACTERE BLANC DE CE NOM
         N = INDEX( NMVO, ' ' )

C        LE 2-EME VOLUME EVENTUELLEMENT VIDE => FRONTIERE
         DO NV2 = 0, NBVOPA

C           LE COUPLE (NV1,NV2) EXISTE-T-IL DANS FRIN ?
            MNFI = MNFRIN + NV1 + NBVOP1 * NV2
            IF( MCN(MNFI) .GT. 0 ) THEN
C              OUI: IL EXISTE UNE FACE AU MOINS
C
C              DECLARATION DE LA SURFACE DANS LE LEXIQUE DES SURFACES
C              ------------------------------------------------------
C              FORMATION DU NOM DE CETTE SURFACE
C              NOM SUFFIXE .F SI ELLE EST FRONTIERE
C                          .I NUMERO DU VOLUME2 SI INTERFACE 2 VOLUMES
               I = N
               IF( NV2 .EQ. 0 ) THEN
C                 FRONTIERE DU VOLUME NV1
                  IF( I.EQ.0 .OR. I.GT.23 ) I = 23
                  NMSU = NMVO(1:I-1) // '.F'
               ELSE
C                 INTERFACE ENTRE LES VOLUMES NV1 ET NV2
                  CALL NUMCHA( LADEFI(WUVOPA-1+NV2), NBC, NMNO )
C                 LES 2 CARACTERES A LA PLACE DES BLANCS OU FINAUX
                  IF( I .EQ. 0 .OR. I .GT. 23-NBC ) I = 23-NBC
                  NMSU = NMVO(1:I-1) // '.I' // NMNO(1:NBC)
               ENDIF
C
C              LA SURFACE NMSU EXISTE-T-ELLE ?
               CALL LXLXOU( NTSURF, NMSU, NTLXSU, MNLXSU )
               IF( NTLXSU .GT. 0 ) THEN
C                 LA SURFACE EXISTE.ELLE EST DETRUITE
                  CALL MAILDS( 'SURFACE', NMSU )
                  CALL LXLXDS(  NTSURF  , NMSU )
               ENDIF
C              LA SURFACE EST DECLAREE PUIS OUVERTE
               CALL LXLXDC( NTSURF, NMSU, NBCNSU, MXNMSU )
C              ELLE EST OUVERTE
               CALL LXLXOU( NTSURF, NMSU, NTLXSU, MNLXSU )
C              IMPRESSION SIGNALANT LA CREATION
               IF( LANGAG .EQ. 0 ) THEN
                 WRITE(IMPRIM,*)'tetrpavo: SURFACE TRIANGULEE CREEE : ',
     %                            NMSU
               ELSE
                 WRITE(IMPRIM,*)'tetrpavo: CREATION of the TRIANGULATED
     %SURFACE: ',NMSU
               ENDIF

C              DECLARATION OUVERTURE INITIALISATION DU TABLEAU 'DEFINITION'
C              ICI SURFACE INTERSECTION DE 2 VOLUMES
C              ------------------------------------------------------------
               MOTS = 1 + WUVOL2
               CALL LXTNDC( NTLXSU, 'DEFINITION', 'ENTIER', MOTS )
               CALL LXTSOU( NTLXSU, 'DEFINITION',  NTDFSU , MNDFSU )
C              VARIABLE NTYTRS  ID POUR IDENTITE,...'
               MCN( MNDFSU + WTYTRS ) = 1
C              VARIABLE NUTYSU 'NUMERO DU TYPE DE LA SURFACE'  ENTIER
C              55 : 'SURFACE COMMUNE A 2 VOLUMES' ) ;
               MCN( MNDFSU + WUTYSU ) = 55
C              VARIABLE NUVOL1 'NOM DU VOLUME 1'  ^~>VOLUME ;
               MCN( MNDFSU + WUVOL1 ) = LADEFI( WUVOPA-1+NV1 )
C              VARIABLE NUVOL2 'NOM DU VOLUME 2'  ^~>VOLUME ;
               IF( NV2 .GT. 0 ) THEN
                  MCN( MNDFSU + WUVOL2 ) = LADEFI( WUVOPA-1+NV2 )
               ELSE
C                 SURFACE CONTOUR D'UN VOLUME
                  MCN( MNDFSU + WUVOL2 ) = LADEFI( WUVOPA-1+NV1 )
               ENDIF
C              LA DATE DE CREATION
               CALL ECDATE( MCN(MNDFSU) )
C              LE NUMERO DU TABLEAU DESCRIPTEUR DE CE TMS
               MCN( MNDFSU + MOTVAR(6) ) =
     %         NONMTD( '~>SURFACE>>DEFINITION' )

C              DECLARATION OUVERTURE INITIALISATION DU TABLEAU 'NSEF'
C              ------------------------------------------------------
C              CALCUL DU NOMBRE DE NSEF
               NBT   = 0
C              LE PREMIER SOUSOBJET=FACE DANS LEFACO
               NFACO = MCN( MNFI )
 60            IF( NFACO .GT. 0 ) THEN
                  NBT = NBT + 1
C                 LA FACE SUIVANTE
                  NFACO = MCN( MNFACO + 11 * NFACO + 8 )
                  GOTO 60
               ENDIF
               MOTS = WUSOEF + 4 * NBT
               CALL LXTNDC( NTLXSU, 'NSEF', 'ENTIER', MOTS )
               CALL LXTSOU( NTLXSU, 'NSEF',  NTFASU , MNFASU )
C
C              LE TYPE D'OBJET
               MCN( MNFASU + WUTYOB ) = 3
C              SURFACE FERMEE
               MCN( MNFASU + WUTFMA ) = 1
C              PAS DE TANGENTES
               MCN( MNFASU + WBTGEF ) = 0
               MCN( MNFASU + WBEFTG ) = 0
               MCN( MNFASU + WBEFAP ) = 0
C              LE TYPE DE MAILLAGE ICI NON STRUCTURE
               MCN( MNFASU + WUTYMA ) = 0
C              LE NOMBRE DE SOMMETS PAR FACE
               MCN( MNFASU + WBSOEF ) = 4
C              LE NOMBRE DE SES NSEF
               MCN( MNFASU + WBEFOB ) = NBT
C              ADRESSE DU 1-ER SOUS OBJET DE CETTE SURFACE
               MNSSS = MNFASU + WUSOEF - 1

C              LE NUMERO DE SES NBT FACES
               NOSOCU(4) = 0
C              LA PREMIERE FACE
               NFACO = MCN( MNFI )
 70            IF( NFACO .GT. 0 ) THEN
C                 L'ADRESSE DE LA FACE NF DANS NOFACE
                  MNFA  = MNFACO + 11 * NFACO - 1
                  DO K=1,3
C                    LE NUMERO DE SOMMET K DE LA FACE NFACO
                     NOSOCU(K) = MCN( MNFA + K )
                  ENDDO

C                 L'INTEGRATION DE CETTE FACE DANS LES NSEF
                  DO K=1,4
                     MCN( MNSSS + K ) = NOSOCU( K )
                  ENDDO
                  MNSSS = MNSSS + 4
C                 LA FACE SUIVANTE
                  NFACO = MCN( MNFA + 9 )
                  GOTO 70
               ENDIF

C              RECHERCHE DU NOMBRE ET DES NUMEROS DES SOMMETS
C              DE CETTE SURFACE
               CALL AZEROI( NBSOMM+1, MCN(MNLIPO) )
C              PARCOURS DES SOMMETS DES NSEF DE CETTE SURFACE
               MNS = MNFASU + WUSOEF - 1
               DO J=1,NBT
                  DO I=1,3
C                    LE NUMERO DU SOMMET DANS LE VOLUME PARTITION
                     NP = MCN( MNS + I )
C                    IL EST MARQUE
                     MCN( MNLIPO + NP ) = 1
                  ENDDO
                  MNS = MNS + 4
               ENDDO

C              LES SOMMETS SONT RENUMEROTES DE 1 A NBS
               NBS = 0
               DO I=1,NBSOMM
                  IF( MCN(MNLIPO+I) .GT. 0 ) THEN
                     NBS = NBS + 1
                     MCN(MNLIPO+I) = NBS
                  ENDIF
               ENDDO

C              PARCOURS DES SOMMETS DES NSEF DE CETTE SURFACE POUR
C              OBTENIR UNE NUMEROTATION LOCALE DES SOMMETS
               MNS = MNFASU + WUSOEF - 1
               DO J=1,NBT
                  DO I=1,3
C                    LE NUMERO DU SOMMET DANS LE VOLUME PARTITION
                     NP = MCN( MNS + I )
C                    LE NOUVEAU NUMERO LOCAL
                     MCN( MNS + I ) = MCN( MNLIPO + NP )
                  ENDDO
                  MNS = MNS + 4
               ENDDO

C              LE TYPE INCONNU DE FERMETURE DU MAILLAGE
               MCN( MNFASU + WUTFMA ) = -1
C              PAS DE TANGENTES STOCKEES
               MCN( MNFASU + WBTGEF ) = 0
               MCN( MNFASU + WBEFAP ) = 0
               MCN( MNFASU + WBEFTG ) = 0
C              LA DATE DE CREATION DU TABLEAU 'NSEF' DE CETTE SURFACE
               CALL ECDATE( MCN(MNFASU) )
C              LE NUMERO DU TABLEAU DESCRIPTEUR DE CE TMS
               MCN( MNFASU + MOTVAR(6) ) = NONMTD( '~>>>NSEF' )

C              DECLARATION OUVERTURE INITIALISATION
C              DU TABLEAU 'XYZSOMMET' DE CETTE SURFACE
C              ---------------------------------------
               MOTS = WYZSOM + 3 * NBS
               CALL LXTNDC( NTLXSU, 'XYZSOMMET', 'ENTIER', MOTS )
               CALL LXTSOU( NTLXSU, 'XYZSOMMET',  NTSUSO, MNSUSO )
C              LE NOMBRE DE SOMMETS
               MCN( MNSUSO + WNBSOM ) = NBS
C              LE NOMBRE DE TANGENTES
               MCN( MNSUSO + WNBTGS ) = 0
C              LE NOMBRE DE COORDONNEES PAR SOMMET
               MCN( MNSUSO + WBCOOR ) = 3
C              L'ADRESSE DES COORDONNEES DE CE VOLUME
               MNSC = MNSUSO + WYZSOM
C              LES COORDONNEES DES NBS SOMMETS
               DO J=1,NBSOMM
                  NP = MCN( MNLIPO + J )
                  IF( NP .GT. 0 ) THEN
C                    LES 3 COORDONNEES DU SOMMET NP
                     RMCN( MNSC   ) = REEL( PTXYZD(1,J) )
                     RMCN( MNSC+1 ) = REEL( PTXYZD(2,J) )
                     RMCN( MNSC+2 ) = REEL( PTXYZD(3,J) )
                     MNSC = MNSC + 3
                  ENDIF
               ENDDO
C              LA DATE DE CREATION DE 'XYZSOMMET' DE CETTE SURFACE
               CALL ECDATE( MCN(MNSUSO) )
C              AJOUT DU NUMERO DU TABLEAU DESCRIPTEUR
               MCN( MNSUSO + MOTVAR(6) ) = NOTADC
            ENDIF
         ENDDO
      ENDDO

 9999 RETURN
      END
