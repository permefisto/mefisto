      SUBROUTINE HAFAVO( NMVOLU, NBMOTS,
     %                   NTFAVO, MNFAVO, NBSOVL, IERR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CREER SI CE N'EST DEJA FAIT LE TMS DES FACES D'UN VOLUME
C -----    ET CHAINER LES FACES FRONTALIERES
C             CHAINER LES FACES INTERFACES ENTRE 2 MATERIAUX
C          cf ~/td/d/a___face

C ENTREES:
C --------
C NMVOLU : NOM DU VOLUME
C NBMOTS : NOMBRE DE MOTS D'INFORMATIONS A STOCKER EN PLUS PAR FACE(>=3)

C SORTIES:
C --------
C NTFAVO : NUMERO DU TMS DES FACES DU VOLUME  cf ~/td/d/a___face
C          =0 SI NON CREE PAR CAUSE D'ERREUR
C MNFAVO : ADRESSE MCN DU TABLEAU DES FACES DU VOLUME
C          =0 SI NON CREE PAR CAUSE D'ERREUR
C NBSOVL : NOMBRE DE SOMMETS DU MAILLAGE DU VOLUME
C IERR   : =0 SI PAS D'ERREUR
C          >0 SI ERREUR ET LE TMS DES FACES N'EXISTE PAS NTFAVO=MNFAVO=0

C EN SORTIE  LE TABLEAU LFACES CONTIENT :
C          LFACES(1,I)= NO DU 1-ER  SOMMET DE LA FACE
C          LFACES(2,I)= NO DU 2-EME SOMMET > 1-ER  SOMMET
C          LFACES(3,I)= NO DU 3-EME SOMMET DE LA FACE
C          LFACES(4,I)= NO DU 4-EME SOMMET DE LA FACE
C                       0 SI TRIANGLE
C          LFACES(5,I)= CHAINAGE HACHAGE SUR FACE SUIVANTE
C          SI SOMMET 2 < DERNIER SOMMET  => FACE   DIRECTE DANS L'ELEMENT FINI
C                      >                 => FACE INDIRECTE
C          UNE FACE DIRECTE EST VUE DE L EXTERIEUR DE L'ELEMENT FINI
C          SOUS LA FORME DIRECTE
C          LFACES(6,I)= NUMERO DU 1-ER ELEMENT FINI CONTENANT CETTE FACE
C                       >0 SI FACE   DIRECTE DANS CET ELEMENT FINI
C                       <0 SI FACE INDIRECTE DANS CET ELEMENT FINI

C          SI LA FACE APPARTIENT A 2 ELEMENTS FINIS ALORS
C          LFACES(7,I)= NUMERO DU 2-EME ELEMENT FINI CONTENANT CETTE FACE
C                       >0 SI FACE   DIRECTE DANS CET ELEMENT FINI
C                       <0 SI FACE INDIRECTE DANS CET ELEMENT FINI
C          SINON
C          LFACES(7,I)= NUMERO DE LA FACE FRONTALIERE SUIVANTE
C                       0 SI C'EST LA DERNIERE
C          L1FAFR = NUMERO DE LA PREMIERE FACE FRONTALIERE DANS LFACES
C                   CHAINAGE SUIVANT DANS LFACES(7,I)
C          L1FA2M = NUMERO DE LA PREMIERE FACE INTERFACE 2 MATERIAUX DANS LFACES
C                   CHAINAGE SUIVANT DANS LFACES(8,I)

C          SI VOLUME MULTI-MATERIAUX
C          LFACES(8,I)= NO DE LA FACE INTERFACE ENTRE 2 MATERIAUX SUIVANTE
C                       0 SI FACE SANS INFORMATION SUPPLEMENTAIRE
C          SINON VOLUME MONO-MATERIAU
C          LFACES(8,I)= NUMERO DE LA FACE A TANGENTE
C                       0 SI FACE SANS TANGENTE
C          Attention:   L'ENTIER 8 en cas de PLUSIEURS MATERIAUX
C                       est UTILISE car au 30/5/2008 Mefisto ne contient
C                       PAS de VOLUMES a FACES avec TANGENTES!
C                       Passer a 9 mots sinon

C          NUTGFA(1:8,1:NBFATG)= NUMERO DES 8 TANGENTES DES FACES A TG
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC PARIS    DECEMBRE 1996
C MODIF  : PERRONNET ALAIN LJLL UPMC ET ST PIERRE DU PERRAY     Mai 2008
C2345X7..............................................................012
      include"./incl/langue.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/a___face.inc"
      include"./incl/a___materiaux.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/gsmenu.inc"
      include"./incl/sotgfc.inc"
      include"./incl/trvari.inc"

      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
      CHARACTER*(*)     NMVOLU
      INTEGER           NOSOEL(64), NUSTFA(4), NUSTF1(4)
      INTEGER           NUTGFA(8)
      INTEGER           NBTYEF(0:9)
      LOGICAL           AVANT
      REAL              X(4),Y(4),Z(4)

      IERR   = 0
      NBF3EF = 0
      CALL AZEROI( 10, NBTYEF )

C     LE VOLUME INITIAL
C     =================
C     LE TABLEAU LEXIQUE DE CE VOLUME
      CALL LXLXOU( NTVOLU, NMVOLU, NTLXVL, MN )
      IF( NTLXVL .LE. 0 ) THEN
         NBLGRC(NRERR) = 2
         KERR(1) = 'hafavo:' // NMVOLU
         IF( LANGAG .EQ. 0 ) THEN
            KERR(2) = 'VOLUME INCONNU'
         ELSE
            KERR(2) = 'UNKNOWN VOLUME'
         ENDIF
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF

C     LE TABLEAU 'NSEF' DE CE VOLUME
      CALL LXTSOU( NTLXVL, 'NSEF', NTCUVL, MNCUVL )
      IF( NTCUVL .LE. 0 ) THEN
         NBLGRC(NRERR) = 2
         KERR(1) = NMVOLU
         IF( LANGAG .EQ. 0 ) THEN
            KERR(2) = 'VOLUME VOLUME SANS NSEF'
         ELSE
            KERR(2) = 'VOLUME WITHOUT NSEF'
         ENDIF
         CALL LEREUR
         IERR = 2
         RETURN
      ENDIF
      NBEFOB = MCN( MNCUVL + WBEFOB )

C     LE TABLEAU 'XYZSOMMET' DE CE VOLUME
      CALL LXTSOU( NTLXVL, 'XYZSOMMET', NTSOVL, MNSOVL )
      IF( NTSOVL .LE. 0 ) THEN
         NBLGRC(NRERR) = 2
         KERR(1) = NMVOLU
         IF( LANGAG .EQ. 0 ) THEN
            KERR(2) = 'VOLUME SANS XYZSOMMET'
         ELSE
            KERR(2) = 'VOLUME WITHOUT XYZSOMMET'
         ENDIF
         CALL LEREUR
         IERR = 3
         RETURN
      ENDIF
C     LE NOMBRE DE SOMMETS DU MAILLAGE DU VOLUME
      NBSOVL= MCN( MNSOVL + WNBSOM )
C     LE NOMBRE DE TANGENTES DU MAILLAGE DU VOLUME
      NBTGS = MCN( MNSOVL + WNBTGS )
C     LE TABLEAU 'XYZSOMMET' DE CE VOLUME

      CALL LXTSOU( NTLXVL, 'MATERIAUX', NTDMEF, MNDMEF )
      IF( NTDMEF .LE. 0 ) THEN
C        NOMBRE DE MATERIAUX
         NBDM = 0
C        NOMBRE D'ELEMENTS FINIS DU MAILLAGE
         NBDMEF = 0
C        ADRESSE MCN DU NUMERO DU MATERIAU DE CHAQUE ELEMENT FINI
         MNUDMEF = 0
      ELSE
C        NOMBRE DE MATERIAUX
         NBDM = MCN( MNDMEF + WNBDM )
C        NOMBRE D'ELEMENTS FINIS DU MAILLAGE
         NBDMEF = MCN( MNDMEF + WBDMEF )
         IF( NBEFOB .NE. NBDMEF ) THEN
            IF( LANGAG .EQ. 0 ) THEN
               NBLGRC(NRERR) = 2
               KERR(1) = NMVOLU
               IF( LANGAG .EQ. 0 ) THEN
            KERR(2) ='NOMBRE MATERIAUX DIFFERENTS MAILLAGE et MATERIAUX'
               ELSE
            KERR(2) ='DIFFERENT NUMBER of FINITE ELEMENTS in MESH and MA
     %TERIALS'
               ENDIF
               CALL LEREUR
               IERR = 4
               RETURN
            ENDIF
         ENDIF
C        ADRESSE-1 MCN DU TABLEAU NUDMEF : MATERIAU DE CHAQUE ELEMENT FINI
         MNUDMEF = MNDMEF + WUDMEF - 1
      ENDIF

C     LE TABLEAU 'FACE' DE CE VOLUME
C     ------------------------------
      CALL LXTSOU( NTLXVL, 'FACE', NTFAVO, MNFAVO )
      IF( NTFAVO .GT. 0 ) THEN
C        LES FACES SONT ELLES ANTERIEURES AUX SOMMETS ?
         IF( AVANT( MCN(MNSOVL), MCN(MNFAVO) ) ) THEN
C           LES SOMMETS SONT ANTERIEURS AUX FACES
C           LE TABLEAU FACE EST ACCEPTE
            print*,'hafavo: SOMMETS ANTERIEURS AUX FACES => TABLEAU FACE
     % ACCEPTE'
            RETURN
         ENDIF
C        LE TABLEAU FACE EST DETRUIT POUR ETRE ENSUITE REGENERE
         CALL LXTSDS( NTLXVL, 'FACE' )
         NTFAVO = 0
         MNFAVO = 0
      ENDIF

C     LES PARAMETRES DES NO SOMMET DU MAILLAGE DU VOLUME
      CALL NSEFPA( MCN(MNCUVL),
     %             NUTYMA, NBSOEL, NBSOEF, NBTGEF,
     %             LDAPEF, LDNGEF, LDTGEF, NBCUVO,
     %             NX    , NY    , NZ    ,
     %             IERR   )
      IF( IERR .NE. 0 ) RETURN

C     NOMBRE DE TETRAEDRES, PENTAEDRES, HEXAEDRES
      DO N=1,NBCUVO
C        LE NUMERO DES NBSOEF SOMMETS DU SOUSOBJET N
         CALL NSEFNS( N     , NUTYMA, NBSOEF, NBTGEF,
     %                LDAPEF, LDNGEF, LDTGEF,
     %                MNCUVL, NX, NY, NZ,
     %                NCOGEL, NUGEEF, NUEFTG, NOSOEL, IERR )
C        NCOGEL: CODE GEOMETRIQUE DE L ELEMENT
C                1:POINT 2:SEGMENT 3:TRIANGLE 4:QUADRANGLE
C                5:TETRAEDRE 6:PENTAEDRE 7:HEXAEDRE  8:6-CUBE 9:PYRAMIDE
C                >9:ERREUR
C        NOMBRE D'EF SELON LE CODE GEOMETRIQUE
         NBTYEF(NCOGEL) = NBTYEF(NCOGEL) + 1
      ENDDO

      NBTETR = NBTYEF( 5 )
      NBPENT = NBTYEF( 6 )
      NBHEXA = NBTYEF( 7 )

C     NOMBRE D'ENTIERS POUR UNE FACE
      MOFACE = 5 + NBMOTS

C     MAJORATION DU NOMBRE DE FACES : FORMULE A MODIFIER EVENTUELLEMENT
      IF( NBCUVO .LT. 10000 ) THEN
ccc         MXFACE = NBCUVO * 6
         MXFACE = 4 * NBTETR + 5 * NBPENT + 6 * NBHEXA
      ELSE
C        0.61 CORRECT SI PAS D'ERREUR POUR UN ELEMENT FINI SUBDIVISE EN N**3 SOUS-EF
ccc         MXFACE = INT( ( 6 * NBCUVO ) * 0.61 )
         MXFACE = INT( ( 4 * NBTETR + 5 * NBPENT + 6 * NBHEXA ) * 0.61 )
      ENDIF

C     LE MAXIMUM DE FACES A TG = 6 FOIS LE NOMBRE D'ELEMENTS FINIS A TG
C     LE NOMBRE D'ELEMENTS FINIS (TETRAEDRES,PENTAEDRES,HEXAEDRES) AVEC TG DU VOLUME
 5    MXFATG = 6 * MCN( MNCUVL + WBEFTG )

C     LE NOMBRE DE FACES AVEC TG
      NBFATG = 0

C     ADRESSAGE ET OUVERTURE DU TABLEAU FACES
C     ---------------------------------------
      CALL LXTNDC( NTLXVL, 'FACE', 'ENTIER',
     %             WFACES + MOFACE*MXFACE + 8*MXFATG )
      CALL LXTSOU( NTLXVL, 'FACE', NTFAVO, MNFAVO )

C     LE TABLEAU DES FACES EST INITIALISE A ZERO
      MNFACE = MNFAVO + WFACES
      CALL AZEROI( MOFACE*MXFACE, MCN(MNFACE) )

C     LE TABLEAU DES NUMEROS DES 8 TG DES FACES A TG
      MNTGFA = MNFACE + MOFACE*MXFACE

      DO 10 NCOGEL=5,9
C        LE NUMERO DES SOMMETS DES FACES DU TETRAEDRE, PENTAEDRE,
C        HEXAEDRE  ET PYRAMIDE
         IF( NCOGEL .EQ. 8 ) GOTO 10
         CALL SOFACU( NCOGEL, NBSTFA(1,NCOGEL), NOSTFA(1,1,NCOGEL) )
         IF( NBTGS .GT. 0 ) THEN
C           LE NUMERO DES TANGENTES DES FACES DU TETRAEDRE, PENTAEDRE, HEXAEDRE
            CALL TGFACU( NCOGEL, NBTGFA(1,NCOGEL), NOTGFA(1,1,NCOGEL) )
         ENDIF
 10   CONTINUE

C     LA 1-ERE FACE LIBRE
      LIBREF = MXFACE

C     LA BOUCLE SUR LES ELEMENTS FINIS DU MAILLAGE DU VOLUME
C     ------------------------------------------------------
      NBFAPB = 0
      DO 1000 N=1,NBCUVO
C
C        LE NUMERO DES SOMMETS DE L'ELEMENT FINI N
         CALL NSEFNS( N     , NUTYMA, NBSOEF, NBTGEF,
     %                LDAPEF, LDNGEF, LDTGEF,
     %                MNCUVL, NX, NY, NZ,
     %                NCOGEL, NUGEEF, NUEFTG, NOSOEL, IERR )
         IF( IERR .NE. 0 ) GOTO 1001

C        LE NOMBRE DE FACES DE L'ELEMENT FINI N
         NFACE = NBFACE( NCOGEL )

C        BOUCLE SUR LES FACES DE L'ELEMENT FINI
C        --------------------------------------
         DO 200 NF=1,NFACE

C           LE NOMBRE DE SOMMETS DE LA FACE
            NS = NBSTFA( NF, NCOGEL )

C           LE NUMERO DES NS SOMMETS DE LA FACE NF
            DO 130 K=1,NS
               NUSTFA(K) = NOSOEL( NOSTFA(K,NF,NCOGEL) )
 130        CONTINUE

C           PRESENTATION DES NUMEROS DES SOMMETS POUR LE HACHAGE
C           ----------------------------------------------------
C           LE PLUS PETIT NUMERO DE SOMMET EST MIS EN PREMIERE POSITION DE NUSTF
C           PUIS LE PLUS PETIT SOMMET PARMI LE SECOND ET LE NS-EME
            CALL NSFAHA( NS, NUSTFA, LESENS, NUSTF1, MINNUS )

C           RECHERCHE OU ADJONCTION DE LA FACE
C           ----------------------------------
            CALL HACHAG( NS, NUSTF1, MOFACE, MXFACE, MCN(MNFACE), 5,
     &                   LIBREF, NOFAC )
            IF( NOFAC .EQ. 0 ) THEN
C              LISTE SATUREE CAR MXFACE SOUS-ESTIMEE
               CALL LXTSDS( NTLXVL, 'FACE' )
C              MAJORATION SAUF SI DES ELEMENT FINI DISJOINTS
               MXFACE = NBCUVO * 6
               GOTO 5
            ENDIF

            MN = MNFACE + MOFACE * ( ABS(NOFAC) - 1 )

C           STOCKAGE ADRESSE DE L'ELEMENT FINI DANS LFACES(6,...;NOFAC)
C           -----------------------------------------------------------
            K = 4
 140        K = K + 1
            IF( K .LT. 7 ) THEN
C
C              RECHERCHE D'UN NUMERO NUL D'ELEMENT FINI PARMI LES 2 POSSIBLES
               IF( MCN(MN + K) .NE. 0 ) GOTO 140
C
C              NO ELEMENT FINI > 0 SI LA FACE EST   DIRECTE DANS L'ELEMENT FINI
C                              < 0 SI LA FACE EST INDIRECTE DANS L'ELEMENT FINI
               MCN( MN + K ) = N * LESENS
C
               IF( NBDMEF .NE. 0 .AND. K .EQ. 6 ) THEN
C                 FACE APPARTENANT A PLUSIEURS VOLUMES?
C                 => FACE A L'INTERFACE ENTRE 2 MATERIAUX
C                 NO DE L'ELEMENT FINI 1 DE CETTE FACE
                  NEF1 = ABS( MCN( MN + 5 ) )
C                 NO DE MATERIAU DE L'ELEMENT FINI 1 DE CETTE FACE
                  NOVOEF1 = MCN( MNUDMEF + NEF1)
C                 NO DE MATERIAU DE L'ELEMENT FINI 2 DE CETTE FACE
                  NOVOEF2 = MCN( MNUDMEF + N )
                  IF( NOVOEF1 .NE. NOVOEF2 ) THEN
C                    => FACE A L'INTERFACE ENTRE 2 MATERIAUX
C                    Z CETTE INSTRUCTION SUPPOSE DES FACES SANS TANGENTES!...
                     MCN( MN + 7 ) = NOVOEF1 + 10000 * NOVOEF2
                  ELSE
                     MCN( MN + 7 ) = 0
                  ENDIF
               ENDIF
C
            ELSE
C
C              ERREUR: IL Y A PLUS DE 2 ELEMENTS FINIS CONTENANT CETTE FACE
               IF( NBFAPB .EQ. 0 ) THEN
C
C                 INITIALISATION DU TRACE DES FACES A PROBLEME
C                 AUCUN ITEM SUR L'ECRAN
ccc                  CALL EFFACE
                  CALL EFFACEMEMPX
                  CALL ITEMS0
C                 LES PARAMETRES DU CADRE MAXIMAL
                  CALL VISEE0
C
               ENDIF
C
C              UNE FACE A PB DE PLUS
               NBFAPB = NBFAPB + 1

C              NOMBRE DE SOMMETS DE L'ELEMENT FINI
               NBS = NBSOME(NCOGEL)
               MM  = MNSOVL + WYZSOM - 4
               DO 145 L=1,NS
C                 LE NUMERO DU SOMMET ET SES XYZ SONT STOCKES DANS X Y Z
                  NOSM = NUSTF1(L) * 3
                  X(L) = RMCN(MM+NOSM+1)
                  Y(L) = RMCN(MM+NOSM+2)
                  Z(L) = RMCN(MM+NOSM+3)
 145           CONTINUE

C              TRACE DE LA FACE : UN TRIANGLE OU QUADRANGLE
               CALL FAP13D( NCROUG, NCNOIR, 0, NS, X, Y, Z )
               CALL TRFINS( 'hafavo: face dans au moins 3 EF!... ' )
               NBF3EF = NBF3EF + 1

               WRITE(IMPRIM,10200) N, (NOSOEL(L),L=1,NBS)
               WRITE(IMPRIM,10201) (NUSTF1(L),L=1,NS)
               WRITE(IMPRIM,10202) (NUSTF1(L),X(L),Y(L),Z(L),L=1,NS)

C              IMPRESSION DES NO ET XYZ DES SOMMETS DE LA FACE A PB
C              NO DE L'ELEMENT FINI 1 DE CETTE FACE
               NEF1 = ABS( MCN( MN + 5 ) )
C              NO DE L'ELEMENT FINI 2 DE CETTE FACE
               NEF2 = ABS( MCN( MN + 6 ) )

               WRITE(IMPRIM,10203) NEF1, NEF2

C              LE NUMERO DES SOMMETS DE L'ELEMENT FINI NEF1
               CALL NSEFNS( NEF1,   NUTYMA, NBSOEF, NBTGEF,
     %                      LDAPEF, LDNGEF, LDTGEF,
     %                      MNCUVL, NX, NY, NZ,
     %                      NCOGEL, NUGEEF, NUEFTG, NOSOEL, IERR )
               NBS = NBSOME(NCOGEL)
               WRITE(IMPRIM,10200) NEF1, (NOSOEL(L),L=1,NBS)

C              LE NUMERO DES SOMMETS DE L'ELEMENT FINI NEF2
               CALL NSEFNS( NEF2,   NUTYMA, NBSOEF, NBTGEF,
     %                      LDAPEF, LDNGEF, LDTGEF,
     %                      MNCUVL, NX, NY, NZ,
     %                      NCOGEL, NUGEEF, NUEFTG, NOSOEL, IERR )
               NBS = NBSOME(NCOGEL)
               WRITE(IMPRIM,10200) NEF2, (NOSOEL(L),L=1,NBS)

C              RETOUR A L'EF N
C              LE NUMERO DES SOMMETS DE L'ELEMENT FINI N
               CALL NSEFNS( N,      NUTYMA, NBSOEF, NBTGEF,
     %                      LDAPEF, LDNGEF, LDTGEF,
     %                      MNCUVL, NX, NY, NZ,
     %                      NCOGEL, NUGEEF, NUEFTG, NOSOEL, IERR )
               WRITE(IMPRIM,10204)

               NBLGRC(NRERR) = 1
               KERR(1)='PLUS DE 2 ELEMENTS FINIS CONTIENNENT CETTE FACE'
               CALL LEREUR
               IERR = 4
               GOTO 200

            ENDIF

C           ELEMENT FINI VOLUMIQUE AVEC TANGENTES?
            IF( NUEFTG .GT. 0 ) THEN
C
C              L'ELEMENT FINI VOLUMIQUE A DES TANGENTES
C              ----------------------------------------
C              PERMUTATION DES TANGENTES DE LA FACE SUIVANT CELLE
C              DES SOMMETS POUR AMENER LE PLUS PETIT NO DE SOMMET EN PREMIERE
C              POSITION DE NUTGFA ET LESENS POUR LA PERMUTATION EVENTUELLE
C              SOMMET 2 ET DERNIER SOMMET DE LA FACE
               CALL NU8TGF( NF, NCOGEL, MINNUS, LESENS, NBT, NUTGFA )

C              LES NUMEROS DES TANGENTES DANS LE MAILLAGE VOLUMIQUE
               DO 150 I=1,NBT
                  NUTGFA(I) = NOSOEL( 8 + NUTGFA(I) )
 150           CONTINUE

C              LA FACE A T ELLE DEJA DES TANGENTES ?
               NUFATG = MCN( MN + 7 )
               IF( NUFATG .EQ. 0 ) THEN

C                 NON. CELLES DE LA FACE NF DU VOLUME DEFINISSENT
C                      CELLES DE LA FACE A TANGENTES
                  K = 0
                  DO 155 L=1,NBT
                     K = K + ABS( NUTGFA(L) )
 155              CONTINUE
                  IF( K .EQ. 0 ) GOTO 200

C                 LA FACE A AU MOINS UNE TANGENTE NON NULLE
                  MNTF   = MNTGFA + 8 * NBFATG
C                 UNE FACE A TANGENTE DE PLUS
                  NBFATG = NBFATG + 1
                  MCN( MN + 7 ) = NBFATG
                  DO 160 L=1,NBT
                     MCN( MNTF + L ) = NUTGFA( L )
 160              CONTINUE
                  DO 165 L=NBT+1,8
                     MCN( MNTF + L ) = 0
 165              CONTINUE

               ELSE

C                 OUI. CELLES DE LA FACE NF DU VOLUME
C                      SONT REFONDUES AVEC CELLES PRECEDENTES
C                      AVEC PRIORITE A LA PRECEDENTE NON NULLE,
C                      PUIS, A LA NOUVELLE TG
                  MNTF  = MNTGFA + 8 * NUFATG - 8
                  DO 170 L=1,NBT
                     IF( NUTGFA(L) .NE. 0 ) THEN
                        IF( MCN( MNTF + L ) .EQ. 0 ) THEN
C                          AJOUT DE LA TG SI PAS DE PRECEDENTE
                           MCN( MNTF + L ) = NUTGFA( L )
                        ENDIF
                     ENDIF
 170              CONTINUE

               ENDIF

            ENDIF

 200     CONTINUE
10200 FORMAT(' hafavo: EF ',I9,' de SOMMETS:',T33,8I9)
10201 FORMAT('         La FACE de SOMMETS:',T33,4I9)
10202 FORMAT(' SOMMET',I9,' : X=',G15.7,' Y=',G15.7,' Z=',G15.7)
10203 FORMAT(' APPARTIENT DEJA aux 2 ELEMENTS FINIS',I9,I9,
     %       ' => ERREUR')
10204 FORMAT(' VERIFIER L''UNICITE DE CHAQUE VOLUME DANS LES UNIONS'/)

 1000 CONTINUE

C     BILAN
 1001 IF( IERR .NE. 0 ) THEN
C        DESTRUCTION DU TMS FACE
         CALL LXTSDS( NTLXVL, 'FACE' )
         NTFAVO = 0
         MNFAVO = 0
         GOTO 9999
      ENDIF

C     CHAINAGE DES FACES FRONTALIERES C-A-D APPARTENANT A UN SEUL ELEMENT FINI
      NBFACES= 0
      MN     = MNFACE
      L1FAFR = 0
      NBFAFR = 0
      L1FA2M = 0
      NBFA2M = 0
      DO 2000 N=1,MXFACE

C        ELIMINATION DES FACES NON UTILISEES DANS MNFACE
         IF( MCN( MN ) .EQ. 0 ) GOTO 1999

C        FACE ACTIVE
         NBFACES = NBFACES + 1

         IF( MCN( MN + 6 ) .EQ. 0 ) THEN

C           LA FACE EST FRONTALIERE => ELLE EST CHAINEE AVEC LA PRECEDENTE
            MCN( MN + 6 ) = L1FAFR
            L1FAFR        = N
            NBFAFR        = NBFAFR + 1

         ELSE IF( MCN( MN + 7 ) .GT. 0 ) THEN

C           LA FACE EST A L'INTERFACE ENTRE 2 MATERIAUX
            MCN( MN + 7 ) = L1FA2M
            L1FA2M        = N
            NBFA2M        = NBFA2M + 1
C
         ENDIF
C
C        PASSAGE A LA FACE SUIVANTE
 1999    MN = MN + MOFACE
 2000 CONTINUE

C     LE TABLEAU 'FACE' EST COMPLETE
C     ==============================
C     LE NOMBRE D'ENTIERS PAR FACE
      MCN( MNFAVO + WOFACE ) = MOFACE
C     LA MAJORATION DU NOMBRE DE FACES
      MCN( MNFAVO + WXFACE ) = MXFACE
C     LE NUMERO DE LA PREMIERE FACE FRONTALIERE
      MCN( MNFAVO + W1FAFR ) = L1FAFR
C     LE NOMBRE DE FACES FRONTALIERES
      MCN( MNFAVO + WBFAFR ) = NBFAFR
C     LE NUMERO DE LA PREMIERE FACE INTERFACE
      MCN( MNFAVO + W1FA2M ) = L1FA2M
C     LE NOMBRE DE FACES INTERFACES
      MCN( MNFAVO + WBFA2M ) = NBFA2M
C     LE NOMBRE DE FACES AVEC TANGENTES
      MCN( MNFAVO + WBFATG ) = NBFATG
C     AJOUT DE LA DATE
      CALL ECDATE( MCN(MNFAVO) )
C     AJOUT DU NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNFAVO + MOTVAR(6) ) = NONMTD( '~>>>FACE' )

C     REDUCTION EVENTUELLE
      IF( NBFATG .LT. MXFATG ) THEN
         CALL TAMSRA( NTFAVO, WFACES + MOFACE*MXFACE + 8*NBFATG )
      ENDIF

      print*,'hafavo: Fin de CONSTRUCTION des',NBFACES,
     %' FACES des',NBEFOB,' EF du volume ',NMVOLU

 9999 IF( NBF3EF .GT. 0 ) THEN
         PRINT*,'hafavo: Probleme',NBF3EF,
     %          ' FACES COMMUNES a au moins 3 TETRAEDRES'
      ENDIF

      RETURN
      END
