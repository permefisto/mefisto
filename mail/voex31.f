        SUBROUTINE VOEX31( NUMVOL, NTLXVF, LADEFI,
     %                     NTTETA, MNTETA, NTSOTT, MNSOTT, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : GENERER LA TETRAEDRISATION D'UNE TETRA-PYRA-PENTA-HEXAEDRISATION
C -----
C
C ENTREES:
C --------
C NUMVOL : NUMERO DU VOLUME DANS LE LEXIQUE DES VOLUMES
C NTLXVF : NUMERO DU TABLEAU TMS DU VOLUME FINAL
C LADEFI : TABLEAU DE DEFINITION DU VOLUME FINAL
C          CF '~td/d/a_volume__definition'
C
C SORTIES:
C --------
C NTTETA : NUMERO      DU TMS 'NSEF' DES NUMEROS DES TETRAEDRES
C MNTETA : ADRESSE MCN DU TMS 'NSEF' DES NUMEROS DES TETRAEDRES
C          CF '~td/d/a___nsef'
C NTSOTT : NUMERO      DU TMS 'XYZSOMMET' DU VOLUME
C MNSOTT : ADRESSE MCN DU TMS 'XYZSOMMET' DU VOLUME
C          CF '~td/d/a___xyzsommet'
C IERR   : 0 SI PAS D'ERREUR
C        > 0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS     Janvier 1990
C MODIFS : PERRONNET ALAIN Laboratoire J-L. LIONS UPMC PARIS   Mars 2007
C MODIFS : PERRONNET ALAIN LJLL UPMC et St Pierre du Perray    Juin 2008
C2345X7..............................................................012
      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/a___materiaux.inc"
      include"./incl/a_volume__definition.inc"
      include"./incl/a___face.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
      INTEGER           LADEFI(0:*)
      DOUBLE PRECISION  D,VOLTER
      CHARACTER*24      NMVOLU
      INTEGER           NOSOQU0(3), NOSOQU1(3), NOSOQU2(3)
C
      IERR    = 0
      NBEF    = 0
      MNIDMEF = 0
      MNFDMEF = 0
C
C     LE NUMERO DU VOLUME
      NUVOAS = LADEFI( WUVOAS )
C
C     OUVERTURE DU VOLUME
      CALL LXNLOU( NTVOLU, NUVOAS, NTLXVI, MNLXVI )
      IF( NTLXVI .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'VOLUME INEXISTANT'
         ELSE
            KERR(1) = 'UNKNOWN VOLUME'
         ENDIF
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF
C
C     EXISTENCE OU NON D'UN TMS 'MATERIAUX' POUR CE VOLUME
      CALL LXTSOU( NTLXVI, 'MATERIAUX', NTMATI, MNMATI )
      NBDM   = 0
      NBDMEF = 0
      IF( NTMATI .GT. 0 ) THEN
C        NOMBRE DE MATERIAUX DU VOLUME INITIAL
         NBDM   = MCN( MNMATI + WNBDM )
C        NOMBRE D''ELEMENTS FINIS DU MAILLAGE AVEC NO DE MATERIAU
         NBDMEF = MCN( MNMATI + WBDMEF )
C        ADRESSE MCN DU NO DE MATERIAU DE CHAQUE EF
         MNIDMEF= MNMATI + WUDMEF
      ENDIF
C
C     TRANSFORMATION DU MAILLAGE STRUCTURE OU NON EN NON STRUCTURE
C     ============================================================
C     LES SOMMETS COPIE DES SOMMETS DU MAILLAGE INITIAL
      CALL LXTSOU( NTLXVI, 'XYZSOMMET', NTSOMI, MNSOMI )
      IF( NTSOMI .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'VOLUME SANS SOMMETS'
         ELSE
            KERR(1) = 'VOLUME WITHOUT VERTICES'
         ENDIF
         CALL LEREUR
         IERR = 4
         RETURN
      ENDIF
C
C     LE TABLEAU 'NSEF'
      CALL LXTSOU( NTLXVI, 'NSEF', NTPEHX, MNPEHX )
      IF( NTPEHX .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'VOLUME SANS NO DES SOMMETS DES EF'
         ELSE
            KERR(1) = 'VOLUME WITHOUT FINITE ELEMENT VERTEX NUMBERS'
         ENDIF
         CALL LEREUR
         IERR = 3
         RETURN
      ENDIF
C
C     LES CARACTERISTIQUES DES EF DE CE MAILLAGE INITIAL
      CALL NSEFPA( MCN(MNPEHX),
     %             NUTYMA, NBSOEL, NBSOEF, NBTGEF,
     %             LDAPEF, LDNGEF, LDTGEF, NBPEHX,
     %             NX    , NY    , NZ    ,
     %             IERR   )
C     NUTYMA : 'NUMERO DE TYPE DU MAILLAGE'    ENTIER
C              0 : 'NON STRUCTURE'      , 2 : 'SEGMENT    STRUCTURE',
C              3 : 'TRIANGLE  STRUCTURE', 4 : 'QUADRANGLE STRUCTURE',
C              5 : 'TETRAEDRE STRUCTURE', 6 : 'PENTAEDRE  STRUCTURE',
C              7 : 'HEXAEDRE  STRUCTURE'
C     NBSOEL : NOMBRE DE SOMMETS DES EF
C              0 SI MAILLAGE NON STRUCTURE
C     NBSOEF : NOMBRE DE SOMMETS DE STOCKAGE DES EF
C              ( TETRAEDRE NBSOEL=4  NBSOEF=8 )
C     NBPEHX : NOMBRE DES EF DU MAILLAGE
C     NX, NY, NZ : LE NOMBRE D'ARETES DANS LES DIRECTION X Y Z
C                CF LE TMS ~td/d/a___nsef
C
      IF( NUTYMA .GE. 1 .AND. NUTYMA .LE. 5 ) THEN

         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'VOLUME IMPOSSIBLE A TETRAEDRISER'
         ELSE
            KERR(1) = 'IMPOSSIBLE VOLUME to TETRAHEDRIZE'
         ENDIF
         CALL LEREUR
         IERR = 2
         RETURN

      ELSE IF( NUTYMA .EQ. 0 ) THEN

C        MAILLAGE NON STRUCTURE
         MNSTTE = MNPEHX + WUSOEF
C        VERIFICATION DU VOLUME POSITIF DES PYRAMIDES PENTAEDRES ET HEXAEDRES
         CALL VOLPLUS( MNSOMI, MNPEHX )

      ELSE IF( NUTYMA .GE. 6 ) THEN

C        MAILLAGE STRUCTURE
         CALL TNMCDC( 'ENTIER', 8*NBPEHX, MNSTTE )
         MN = MNSTTE
         DO 100 NUELEM = 1, NBPEHX
C           LE NUMERO DES SOMMETS DE L'ELEMENT NUELEM
            CALL NSEFNS( NUELEM, NUTYMA, NBSOEF, NBTGEF,
     %                   LDAPEF, LDNGEF, LDTGEF,
     %                   MNTETA, NX, NY, NZ,
     %                   NCOGEL, NUGEEF, NUEFTG, MCN(MN), IERR )
            MN = MN + 8
            IF( IERR .NE. 0 ) RETURN
 100     CONTINUE

      ENDIF
C
C     LE TABLEAU DES SOMMETS DES TETRAEDRES
C     MAJORANT DU NOMBRE DE TETRAEDRES
      NBTEMX = 6 * NBPEHX
      CALL TNMCDC( 'ENTIER', NBTEMX*5, MNTETR )
      CALL AZEROI( NBTEMX*5, MCN(MNTETR) )
C
      IF( NBDM .GT. 0 ) THEN
C        CONSTRUCTION DE 'MATERIAUX' DU VOLUME FINAL CAR PLUSIEURS MATERIAUX
         CALL LXTSOU( NTLXVF, 'MATERIAUX', NTMATF, MNMATF )
         IF( NTMATF .GT. 0 ) CALL LXTSDS( NTLXVL, 'MATERIAUX' )
         MOTS2M = WUDMEF + NBTEMX
         CALL LXTNDC( NTLXVF, 'MATERIAUX', 'ENTIER', MOTS2M )
         CALL LXTSOU( NTLXVF, 'MATERIAUX',  NTMATF , MNMATF )
C        ADRESSE MCN DU NO DE MATERIAU DE CHAQUE EF
         MNFDMEF = MNMATF + WUDMEF
         IF( NBDMEF .NE. NBPEHX ) THEN
            NBLGRC(NRERR) = 1
            IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'VOLUME et MATERIAUX avec des  EF DIFFERENTS'
            ELSE
        KERR(1)='DIFFERENT FINITE ELEMENTS between VOLUME and MATERIALS'
            ENDIF
            CALL LEREUR
            IERR = 3
            RETURN
         ENDIF
      ENDIF
C
C     TETRAEDRISATION DES PENTAEDRES ET HEXAEDRES (PAS DES PYRAMIDES)
      CALL TETRPH( NBPEHX, MCN(MNSTTE),
     %             NBDM,   MCN(MNIDMEF), MCN(MNFDMEF),
     %             NBTETR, NBPYRA, MCN(MNTETR) )
C
C     LES SOMMETS COPIE DES SOMMETS DU MAILLAGE INITIAL
      MOTS = WYZSOM + MCN(MNSOMI+WNBSOM) * 3
      CALL LXTNDC( NTLXVF, 'XYZSOMMET', 'ENTIER', MOTS )
      CALL LXTSOU( NTLXVF, 'XYZSOMMET',  NTSOTT, MNSOTT )
      CALL TRTATA( MCN(MNSOMI), MCN(MNSOTT), MOTS )
      CALL ECDATE( MCN(MNSOTT) )
C
C     LE NOMBRE DE TETRAEDRES DE CE VOLUME
      MOTS = WUSOEF + 8 * (NBTETR+NBPYRA)
      CALL LXTNDC( NTLXVF, 'NSEF', 'ENTIER', MOTS )
      CALL LXTSOU( NTLXVF, 'NSEF',  NTTETA , MNTETA )
C
C     LES PREMIERES VARIABLES DU TABLEAU 'NSEF'
C     VARIABLE NUTYOB 'NUMERO DE TYPE DE L''OBJET' ENTIER
      MCN( MNTETA + WUTYOB ) = 4
C     VARIABLE NUTYMA 'NUMERO DE TYPE DU MAILLAGE' ENTIER
      MCN( MNTETA + WUTYMA ) = 0
C     VARIABLE NBSOEF 'NOMBRE DE SOMMETS PAR EF' ENTIER
      MCN( MNTETA + WBSOEF ) = 8
C     VARIABLE NBEFOB 'NOMBRE DE EF DE L''OBJET' ENTIER
      MCN( MNTETA + WBEFOB ) = NBTETR
C
C     LE NUMERO DES SOMMETS DES TETRAEDRES FINAUX
      NOMATE = 0
      MNT    = MNTETR
      MNS    = MNTETA + WUSOEF
      MNXYZ  = MNSOTT + WYZSOM - 3
      DO 200 I=1,NBTETR
C
         IF( MCN( MNT + 4 ) .EQ. 0 ) THEN
C
C           TETRAEDRE
            DO 190 J=0,3
               MCN( MNS + J ) = MCN( MNT + J )
               MCN( MNS + J + 4 ) = 0
 190        CONTINUE
C           SI VOLUME NEGATIF PERMUTATION DES SOMMETS 2 ET 3
            D = VOLTER( RMCN(MNXYZ+3*MCN(MNS  )),
     %                  RMCN(MNXYZ+3*MCN(MNS+1)),
     %                  RMCN(MNXYZ+3*MCN(MNS+2)),
     %                  RMCN(MNXYZ+3*MCN(MNS+3)) )
            IF( D .LT. 0D0 ) THEN
               MOTS           = MCN( MNS + 2 )
               MCN( MNS + 2 ) = MCN( MNS + 3 )
               MCN( MNS + 3 ) = MOTS
            ENDIF

         ELSE
C
C           PYRAMIDE
            DO 195 J=0,4
               MCN( MNS + J ) = MCN( MNT + J )
 195        CONTINUE
            MCN( MNS + 5 ) = 0
C           STOCKAGE DU NO DE MATERIAU
            IF( NBDM .GT. 0 ) NOMATE = MCN(MNFDMEF-1+I)
            MCN( MNS + 6 ) = NOMATE
            MCN( MNS + 7 ) = 0

         ENDIF
C
         MNS = MNS + 8
         MNT = MNT + 5
C
 200  CONTINUE
C
C     LE TYPE INCONNU DE FERMETURE DU MAILLAGE
      MCN( MNTETA + WUTFMA ) = -1
C     PAS DE TANGENTES STOCKEES
      MCN( MNTETA + WBTGEF ) = 0
      MCN( MNTETA + WBEFAP ) = 0
      MCN( MNTETA + WBEFTG ) = 0
C
      CALL ECDATE( MCN(MNTETA) )
      MCN(MNTETA+MOTVAR(6)) = NONMTD( '~>>>NSEF' )
C
      IF( NBPYRA .GT. 0 ) THEN
C
C        TRAITEMENT DES PYRAMIDES PAR RECHERCHE DE LA FACE QUADRANGULAIRE
C        RECHERCHE DES 2 TRIANGLES OPPOSES ET DECOUPAGE DE LA PYRAMIDE EN
C        2 TETRAEDRES A LA FIN DU TABLEAU MCN(MNTETR)
C
C        RECHERCHE DU NOM NMVOLU DE CE VOLUME
C        ADRESSE DU LEXIQUE DES VOLUMES
         CALL TAMSOU( NTMN(4) ,MNLXVF )
C        LE NOMBRE D'ENTIERS PAR NOM DE VOLUME
         NBENN = MCN( MNLXVF + 2 )
C        ADRESSE MCN DU LEXIQUE DU VOLUME DE NUMERO NUMVOL
         MN = MNLXVF + MCN(MNLXVF) * NUMVOL
C        LE NOM NMVOLU DE L'OBJET CONVERTI D'ENTIERS EN CARACTERES
         CALL ENTNOM( NBENN, MCN(MN), NMVOLU )
C
         IF( NBDM .GT. 0 ) THEN
C           NOMBRE DE MATERIAUX DU VOLUME ACTUEL SANS PYRAMIDES
            MCN( MNMATF + WNBDM ) = NBDM
C           NOMBRE D''ELEMENTS FINIS DU MAILLAGE AVEC NO DE MATERIAU
            MCN( MNMATF + WBDMEF ) = NBTETR
         ENDIF
C
C        GENERATION EVENTUELLE PAR HACHAGE DES FACES DES CUBES
C        CHAINAGE DES FACES FRONTALIERES EN POSITION 7
C        AVEC UN LIEN NEGATIF POUR LES FACES FRONTALIERES
C        =====================================================
         CALL HAFAVO( NMVOLU, 3, NTFAVO, MNFAVO, NBSOVL, IERR )
         IF( IERR .NE. 0 ) THEN
            NBLGRC(NRERR) = 2
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'VOLUME: ' // NMVOLU
               KERR(2) = 'IMPOSSIBLE DE CREER SES FACES'
            ELSE
               KERR(1) = 'VOLUME: ' // NMVOLU
               KERR(2) = 'IMPOSSIBLE to CREATE its FACES'
            ENDIF
            CALL LEREUR
            RETURN
         ENDIF
C
C        LE NOMBRE D'ENTIERS PAR FACE
         MOFACE = MCN( MNFAVO + WOFACE )
C        LA MAJORATION DU NOMBRE DE FACES
         MXFACE = MCN( MNFAVO + WXFACE )
C        ADRESSE DU TABLEAU LFACES
         MNFACE = MNFAVO + WFACES
C
C        PARCOURS DES PYRAMIDES DU MAILLAGE
         NBEF = NBTETR
         MNS  = MNTETA + WUSOEF
         MNT  = MNS + 8 * NBEF
         DO 300 I=1,NBTETR
C           LE NO DU 5-EME SOMMET
            NS5 = MCN( MNS + 4 )
            IF( NS5 .GT. 0 ) THEN
C
C              PYRAMIDE: RECUPERATION DU NO DE MATERIAU
               NOMATE = MCN( MNS + 6 )
C
C              LE NO DES 4 SOMMETS DE LA FACE QUADRANGULAIRE
C              SONT LES 4 PREMIERS DE LA PYRAMIDE
               NOSOQU0(1) = MCN( MNS     )
               NOSOQU0(2) = MCN( MNS + 1 )
               NOSOQU0(3) = MCN( MNS + 2 )
C
C              PRESENTATION DES NUMEROS DES SOMMETS POUR LE HACHAGE
C              PLUS PETIT NUMERO DE SOMMET EN PREMIERE POSITION DE NOSOQU2
C              PUIS LE PLUS PETIT SOMMET PARMI LE SECOND ET LE NS-EME
               CALL NSFAHA( 3, NOSOQU0, LESENS, NOSOQU2, MINNUS )
C
C              RECHERCHE DU TRIANGLE 1 2 3 DU QUADRANGLE DANS LES FACES
               CALL HACHAR( 3, NOSOQU2, MOFACE, MXFACE, MCN(MNFACE), 5,
     %                      NOFAC )
C              NOFAC: NO DE LA COLONNE DU TABLEAU VALEUR RETROUVEE
C                     > 0 SI LE TABLEAU VALEUR A ETE RETROUVE
C                     =<0 SI LE TABLEAU VALEUR N'A PAS ETE RETROUVE
               IF( NOFAC .GT. 0 ) THEN
C                 VERIFICATION 2-EME TRIANGLE 1 3 4 EST AUSSI UNE FACE
                  NOSOQU1(1) = MCN( MNS     )
                  NOSOQU1(2) = MCN( MNS + 2 )
                  NOSOQU1(3) = MCN( MNS + 3 )
C
C                 PRESENTATION DES NUMEROS DES SOMMETS POUR LE HACHAGE
C                 PLUS PETIT NUMERO DE SOMMET EN PREMIERE POSITION DE NOSOQU2
C                 PUIS LE PLUS PETIT SOMMET PARMI LE SECOND ET LE NS-EME
                  CALL NSFAHA( 3, NOSOQU1, LESENS, NOSOQU2, MINNUS )
C                 RECHERCHE DU TRIANGLE 1 3 4 DU QUADRANGLE DANS LES FACES
                  CALL HACHAR(3, NOSOQU2, MOFACE, MXFACE, MCN(MNFACE),5,
     %                        NOFAC )
                  IF( NOFAC .GT. 0 ) THEN
C                    IL EXISTE 2 FACES TRIANGULAIRES 123 et 134
C                    LA SUBDIVISION DE LA PYRAMIDE EN 2 TETRAEDRES
C                    PEUT ETRE FAITE
C                    TETRAEDRE 1235
                     MCN( MNS     ) = NOSOQU0(1)
                     MCN( MNS + 1 ) = NOSOQU0(2)
                     MCN( MNS + 2 ) = NOSOQU0(3)
                     MCN( MNS + 3 ) = NS5
                     MCN( MNS + 4 ) = 0
                     MCN( MNS + 5 ) = 0
                     MCN( MNS + 6 ) = 0
                     MCN( MNS + 7 ) = 0
C                    NO DU MATERIAU DEJA INITIALISE DANS TETRPH
C
C                    TETRAEDRE 1345
                     MCN( MNT     ) = NOSOQU1(1)
                     MCN( MNT + 1 ) = NOSOQU1(2)
                     MCN( MNT + 2 ) = NOSOQU1(3)
                     MCN( MNT + 3 ) = NS5
                     MCN( MNT + 4 ) = 0
                     MCN( MNT + 5 ) = 0
                     MCN( MNT + 6 ) = 0
                     MCN( MNT + 7 ) = 0
                     MNT = MNT + 8
C                    NO DU MATERIAU
                     IF( NBDM .GT. 0 ) MCN(MNFDMEF+NBEF) = NOMATE
                     NBEF = NBEF + 1
                     GOTO 290
                  ELSE
                     WRITE(IMPRIM,*)
     %              'ANOMALIE UN TRIANGLE 123 ET PAS L''AUTRE 134!'
                     GOTO 290
                  ENDIF
               ELSE
C                 AUTRE POSSIBILITE TRIANGLE 124 DU QUADRANGLE 1234
                  NOSOQU0(1) = MCN( MNS     )
                  NOSOQU0(2) = MCN( MNS + 1 )
                  NOSOQU0(3) = MCN( MNS + 3 )
C
C                 PRESENTATION DES NUMEROS DES SOMMETS POUR LE HACHAGE
C                 PLUS PETIT NUMERO DE SOMMET EN PREMIERE POSITION DE NOSOQU2
C                 PUIS LE PLUS PETIT SOMMET PARMI LE SECOND ET LE NS-EME
                  CALL NSFAHA( 3, NOSOQU0, LESENS, NOSOQU2, MINNUS )
C
C                 RECHERCHE DU TRIANGLE 1 2 4 DU QUADRANGLE DANS LES FACES
                  CALL HACHAR(3, NOSOQU2, MOFACE, MXFACE, MCN(MNFACE),5,
     %                        NOFAC )
C                 NOFAC: NO DE LA COLONNE DU TABLEAU VALEUR RETROUVEE
C                        > 0 SI LE TABLEAU VALEUR A ETE RETROUVE
C                        =<0 SI LE TABLEAU VALEUR N'A PAS ETE RETROUVE
                  IF( NOFAC .GT. 0 ) THEN
C                    VERIFICATION 2-EME TRIANGLE 234 EST AUSSI UNE FACE
                     NOSOQU1(1) = MCN( MNS + 1 )
                     NOSOQU1(2) = MCN( MNS + 2 )
                     NOSOQU1(3) = MCN( MNS + 3 )
C
C                    PRESENTATION DES NUMEROS DES SOMMETS POUR LE HACHAGE
C                    PLUS PETIT NUMERO DE SOMMET EN PREMIERE POSITION DE NOSOQU2
C                    PUIS LE PLUS PETIT SOMMET PARMI LE SECOND ET LE NS-EME
                     CALL NSFAHA( 3, NOSOQU1, LESENS, NOSOQU2, MINNUS )
C                    RECHERCHE DU TRIANGLE 1 3 4 DU QUADRANGLE DANS LES FACES
                     CALL HACHAR(3, NOSOQU2,MOFACE,MXFACE,MCN(MNFACE),5,
     %                           NOFAC )
                     IF( NOFAC .GT. 0 ) THEN
C                       IL EXISTE 2 FACES TRIANGULAIRES 123 et 134
C                       LA SUBDIVISION DE LA PYRAMIDE EN 2 TETRAEDRES
C                       PEUT ETRE FAITE
C                       TETRAEDRE 1245
                        MCN( MNS     ) = NOSOQU0(1)
                        MCN( MNS + 3 ) = NOSOQU0(2)
                        MCN( MNS + 2 ) = NOSOQU0(3)
                        MCN( MNS + 3 ) = NS5
                        MCN( MNS + 4 ) = 0
                        MCN( MNS + 5 ) = 0
                        MCN( MNS + 6 ) = 0
                        MCN( MNS + 7 ) = 0
C
C                       TETRAEDRE 2345
                        MCN( MNT     ) = NOSOQU1(1)
                        MCN( MNT + 1 ) = NOSOQU1(2)
                        MCN( MNT + 2 ) = NOSOQU1(3)
                        MCN( MNT + 3 ) = NS5
                        MCN( MNT + 4 ) = 0
                        MCN( MNT + 5 ) = 0
                        MCN( MNT + 6 ) = 0
                        MCN( MNT + 7 ) = 0
                        MNT  = MNT + 8
C                       NO DU MATERIAU
                        IF( NBDM .GT. 0 ) MCN(MNFDMEF+NBEF) = NOMATE
                        NBEF = NBEF + 1
                        GOTO 290
                     ELSE
                        WRITE(IMPRIM,*)
     %                 'ANOMALIE UN TRIANGLE 124 ET PAS L''AUTRE 234!'
                        GOTO 290
                     ENDIF
                  ELSE
                     WRITE(IMPRIM,*)'ANOMALIE NI TRIANGLE 123 NI 124!'
                     GOTO 290
                  ENDIF
               ENDIF
            ENDIF
C
 290        MNS = MNS + 8
 300     CONTINUE
C
C        DESTRUCTION DU TMS FACE INTERMEDIAIRE DU VOLUME
         CALL LXTSDS( NTLXVF, 'FACE' )
C
C        MISE A JOUR DU TMS NSEF DU VOLUME TETRAEDRISE
C        VARIABLE NBEFOB 'NOMBRE DE EF DE L''OBJET'  ENTIER
         MCN( MNTETA + WBEFOB ) = NBEF
         CALL ECDATE( MCN(MNTETA) )
      ENDIF
C
      IF( NBDM .GT. 0 ) THEN
C        MISE A JOUR DU TMS 'MATERIAUX' DU VOLUME TETRAEDRISE
C        NOMBRE DE MATERIAUX
         MCN( MNMATF + WNBDM ) = NBDM
C        NOMBRE D'ELEMENTS FINIS DU MAILLAGE
         MCN( MNMATF + WBDMEF ) = NBEF
C        LA DATE DE CREATION
         CALL ECDATE( MCN(MNMATF) )
C        LE NUMERO DU TABLEAU DESCRIPTEUR DE CE TMS 'MATERIAUX'
         MCN( MNMATF + MOTVAR(6) ) = NONMTD( '~>>>MATERIAUX' )
      ENDIF
C
C     FIN DU TRAITEMENT
      IF ( NUTYMA .GE. 6 ) THEN
C        MAILLAGE STRUCTURE
         CALL TNMCDS( 'ENTIER', 8*NBPEHX, MNSTTE )
      ENDIF
      CALL TNMCDS( 'ENTIER', NBTEMX*5, MNTETR )
C
C     MISE A JOUR DES ELEMENTS FINIS DES VOLUMES MATERIAUX
C     SI CE VOLUME EST MULTI-MATERIAUX
C     ----------------------------------------------------
      IF( NBDM .GT. 1 ) THEN
         CALL MAJVOL( NUMVOL, IERR )
C        SI ERREUR => PAS DE MODIFICATION DES VOLUMES
C                     PAS DE PROBLEME CREE => ERREUR ANNULEE
         IERR = 0
      ENDIF
C
      RETURN
      END
