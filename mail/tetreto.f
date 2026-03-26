      SUBROUTINE TETRETO( QUAMINEX, QTEAME, KTITRE,
     %                    MXSOMM, NBSOMM0,NBSOMM, PTXYZD,NPSOFR, VMOYEN,
     %                    MXTETR, N1TEVI, NOTETR, NUDTETR, N1TETS,
     %                    NBTRCF, NOTRCF, INFACO, MXFACO, LEFACO,NO0FAR,
     %                    N1FEVI, MXFETO, NFETOI,
     %                    MXSSET, NBSSET0,N1SSET, NBTECFPR, NOTECFPR,
     %                    MXTECR, NBTECR, NOTECR, VOLTECR, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : TETRAEDRISATION DES NBSSET0 SOUS-ETOILES DE L'ETOILE DES FACES
C ----- DEFINIES DANS LES TABLEAUX N1SSET et NFETOI
C       A L'APPEL LES NBTRCF PREMIERES FACES SONT LES NBTRCF FACES DU CF
C       SI CONFIGURATION EN TORE, ELLES SONT REPETEES 2 FOIS
C       LE NOMBRE DE SOUS-ETOILES EVOLUE AU COURS DE LA TETRAEDRISATION
C       DANS CHAQUE SOUS-ETOILE LES FACES SONT ORIENTEES DE TELLE SORTE
C       QUE LEUR VECTEUR NORMAL SOIT DIRIGE VERS L'INTERIEUR

C REMARQUE:
C ARETGR : LONGUEUR DE L'ARETE SOUHAITEE POUR L'ENSEMBLE DES TETRAEDRES
C          EST INITIALISE DANS DARETE et EXPLOITE dans CALL TAILIDEA
C NOFOTI : =NO DE LA FONCTION UTILISATEUR 'TAILLE_IDEALE'
C          =0 PAS D'EXISTENCE DE LA FONCTION UTILISATEUR 'TAILLE_IDEALE'
C          EST STOCKE dans le COMMON /DEARET/ et OBTENU par
C          NOFOTI = NOFOTIEL()

C ENTREES:
C --------
C QUAMINEX: QUALITE MINIMALE EXIGEE POUR CREER UN TETRAEDRE
C           C-A-D AU DESSOUS DE LAQUELLE UN TETRAEDRE N'EST PAS CREE
C           C'EST AUSSI LE SEUIL DE MEDIOCRITE DE LA QUALITE D'UN TETRAEDRE
C           = QUAMED
C QTEAME : QUALITE DES TETRAEDRES AU DESSOUS DE LAQUELLE UNE
C          AMELIORATION DE LA QUALITE DES TETRAEDRES EST DEMANDEE
C KTITRE : TITRE DU TRACE
C MXSOMM : NOMBRE MAXIMAL DE SOMMETS DECLARABLES DANS PTXYZD NPSOFR
C MXTETR : NOMBRE DE TETRAEDRES DECLARABLES DANS NOTETR
C MXFETO : NOMBRE DE FACES DECLARABLES DANS NFETOI
C MXTECR : NOMBRE MAXIMAL DE TETRAEDRES CREABLES DANS NOTECR
C VMOYEN : VOLUME MOYEN DES TETRAEDRES INITIAUX DE L'ETOILE
C          UTILE POUR ELIMINER LES TETRAEDRES DE VOLUME PROCHE DE ZERO
C NBSOMM0: NOMBRE DE SOMMETS DANS PTXYZD AVANT TRAITEMENT DE CETTE FACE PERDUE

C MODIFIES:
C ---------
C NBSOMM : NUMERO DU DERNIER SOMMET CREE DANS PTXYZD
C PTXYZD : TABLEAU DES COORDONNEES DES POINTS
C          PAR POINT : X  Y  Z DISTANCE_SOUHAITEE
C NPSOFR : =  0 SI POINT AJOUTE NON SUR LA SURFACE DE L'OBJET
C          =  1 SI LE POINT EST AJOUTE SUR UNE ARETE DE LA SURFACE
C          =  2 SI LE POINT EST AJOUTE COMME BARYCENTRE SUR LA SURFACE
C          =  1 000 000  + NUMERO DU  POINT INTERNE UTILISATEUR
C          = (1 000 000) * (NO SURFACE + 1) + NUMERO DU SOMMET DANS
C              LA NUMEROTATION DES SOMMETS DE LA SURFACE
C          = -4 SI LE POINT EST SOMMET D'OT NON TROP PROCHE PT OU FACE
C          = -1 SI LE POINT EST SOMMET D'OT TROP PROCHE PT OU FACE
C          = -3 SI LE POINT EST SOMMET D'OT REFUSE DANS LA TETRAEDRISATION
C          = -NPSOFR(I) SI POINT DEPLACE SUR LA SURFACE
C            DANS LEUR SURFACE FERMEE INITIALE ou NO DE POINT INTERNE
C N1TEVI : NO DU PREMIER TETRAEDRE VIDE DANS NOTETR
C NOTETR : LISTE DES TETRAEDRES
C          SOMMET1,    SOMMET2,    SOMMET3,    SOMMET4,
C          TETRAEDRE1, TETRAEDRE2, TETRAEDRE3, TETRAEDRE4
C          DE L'AUTRE COTE DE LA FACE
C          1: 123      2: 234      3: 341      4: 412
C NUDTETR: PLUS GRAND NUMERO DANS NOTETR DE TETRAEDRE OCCUPE
C N1TETS : N1TETS(I) NUMERO NOTETR D'UN TETRAEDRE AYANT POUR SOMMET I

C NBTRCF : NOMBRE DE TRIANGLES PERDUS FORMANT LE CF
C NOTRCF : SI NOTRCF(*)>0 NUMERO LEFACO DU TRIANGLE
C                      <0 NUMERO NO0FAR DU TRIANGLE AJOUTE AU CF
C INFACO : = 0 PAS DE TABLEAU LEFACO NI DE CONSERVATION DES
C              FACES FRONTIERE ( NO DE VOLUME CONNU PAR NVOLTE )
C          = 1 EXISTENCE DU TABLEAU LEFACO ET CONSERVATION DES
C              FACES DE LA FRONTIERE
C MXFACO : MAX DE FACES DECLARABLES DANS LEFACO
C LEFACO : FACE=TRIANGLE DE LA PEAU OU DES INTERFACES ENTRE VOLUMES
C          IL CONTIENT DANS CET ORDRE
C          NUMERO (DANS PTXYZD) DU SOMMET 1, SOMMET 2, SOMMET 3
C          NUMERO (DANS NUVOPA 0 SINON) DU VOLUME1, VOLUME2 DE LA FACE
C          NUMERO (DANS LEFACO) DE LA FACE ADJACENTE PAR L'ARETE 1 2 3
C
C          ATTENTION: UNE ARETE PEUT APPARTENIR A PLUS DE 2 FACES
C          => CHAINAGE CIRCULAIRE DE CES FACES DANS LEFACO
C          LEFACO(9,*)  -> FACE SUIVANTE (*=0:VIDE, *<>0:NON VIDE)
C          HACHAGE AVEC LA SOMME DES 3 SOMMETS MODULO MXFACO
C          LF = MOD( NOSOFA(1)+NOSOFA(2)+NOSOFA(3) , MXFACO ) + 1
C          NF = LEFACO( 10, LF ) LE NUMERO DE LA 1-ERE FACE DANS LEFACO
C          SI LA FACE NE CONVIENT PAS. PASSAGE A LA SUIVANTE
C          NF  = LEFACO( 9, NF )  ...
C          LEFACO(10,*) PREMIERE FACE DANS LE HACHAGE
C          LEFACO(11,.) = NO NOTETR D'UN TETRAEDRE AYANT CETTE FACE, 0 SINON
cccC          LEFACO(12,.) = NO FACEOC DE 1 A NBFACES D'OC
C NO0FAR : NUMERO DES 3 SOMMETS DE LA FACE AJOUTEE AU CF
C          NORMALE VERS L'INTERIEUR DU TETRAEDRE LA CONTENANT

C N1FEVI : POINTEUR SUR LA PREMIERE FACE VIDE DE L'ETOILE
C          CHAINAGE SUIVANT DANS NFETOI(5,*)
C NFETOI : VERSION 2 LES FACES TRIANGULAIRES DE L'ETOILE
C          1: NUMERO DU TETRAEDRE DANS NOTETR OPPOSE A CETTE FACE
C          2: NUMERO PTXYZD DU SOMMET 1 DE LA FACE
C          3: NUMERO PTXYZD DU SOMMET 2 DE LA FACE
C          4: NUMERO PTXYZD DU SOMMET 3 DE LA FACE
C             S1S2xS1S3 EST LA NORMALE DIRIGEE VERS L'INTERIEUR
C             DE L'ETOILE
C          5: CHAINAGE SUIVANT DES FACES OCCUPEES ET VIDES

C MXSSET : NOMBRE MAXIMAL DE SOUS ETOILES DECLARABLES DANS N1SSET
C NBSSET0: NOMBRE INITIAL DE SOUS-ETOILES STOCKEES DANS N1SSET-NFETOI
C N1SSET : NUMERO NFETOI DE LA 1-ERE FACE DES NBSSET0 SOUS ETOILES

C MODIFIES:
C ---------
C NBTECFPR: NOMBRE DE TETRAEDRES INITIAUX CONSTITUANTS L'ETOILE
C           DE NUMERO NOTETR SAUVEGARDES DANS NOTECFPR
C NOTECFPR: NUMERO DANS NOTETR DES NBTECFPR TETRAEDRES SAUVEGARDES

C SORTIES:
C --------
C NBTECR : NOMBRE DE TETRAEDRES CREES
C NOTECR : NUMERO NOTETR DES NBTECR TETRAEDRES CREES
C VOLTECR: VOLUME DES NBTECR TETRAEDRES CREES DANS L'ETOILE
C IERR   : 0 SI PAS D'ERREUR
C          1 ETOILE AVEC MOINS DE 4 FACES
C          2 SATURATION D'UN TABLEAU
C          3 3 SOMMETS NON DANS UN TETRAEDRE
C          4 NOMBRE INCORRECT DE FACES DANS L'ETOILE
C          5 SOMMETS DES 4 FACES DE L'ETOILE INCORRECTS
C          6 UNE ARETE DES FACES DE L'ETOILE APPARTIENT A MOINS DE 2
C            OU PLUS DE 2 FACES DE L'ETOILE
C          7 ou 8 ou 9 ALGORITHME DEFAILLANT A AMELIORER
C         10 PLUSIEURS ARETES A PROBLEME
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & VEULETTES SUR MER       Aout 2015
C MODIFS : ALAIN PERRONNET LJLL UPMC & VEULETTES SUR MER  Septembre 2017
C MODIFS : ALAIN PERRONNET Saint Pierre du Perray          Decembre 2019
C2345X7..............................................................012
      PARAMETER         ( QUAMED=0.005 )
      PARAMETER         ( MXVPSI=2048, MXCIAS=64, MXASFVP=512 )
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/trvari.inc"
      include"./incl/darete.inc"

      COMMON / TRTETR / STOPTE, TRACTE
      LOGICAL           STOPTE, TRACTE, TRACTE0

      DOUBLE PRECISION  PTXYZD(1:4,1:NBSOMM)
      INTEGER           NPSOFR(MXSOMM), NOTETR(8,MXTETR),
     %                  N1TETS(1:MXSOMM), NOTRCF(NBTRCF),
     %                  LEFACO(11,0:*),
     %                  NO0FAR(3,*), NFETOI(5,MXFETO),
     %                  NOTECFPR(*), NOTECR(MXTECR),
     %                  N1SSET(MXSSET)

      CHARACTER*(*)     KTITRE
      DOUBLE PRECISION  CBTR(3), VOLTECR, V, VMAX, VOLMAX, VMOYEN,
     %                  VOLUMT, VOLMOY, VOLTEC, XYZ(4), VETXYZ, RAVEVM,
     %                  ARMIN, ARMAX, SURFTR(4)
      REAL              QUATET, QUAMIN, Q

      INTEGER           NFVPSI0(MXVPSI), NFVPSI(MXVPSI),
     %                  N1CIAS(MXCIAS), NSASFVP(4,MXASFVP)

      INTEGER           NOSOTE(4), NOSOTR2(3), NOSOTR3(3)
      INTEGER           NSARP2F(1)
      EQUIVALENCE      (NOSOTE(1),NS1), (NOSOTE(2),NS2),
     %                 (NOSOTE(3),NS3), (NOSOTE(4),NS4)
      INTEGER           NTEOPF(4)
      EQUIVALENCE      (NTEOPF(1),NTEOPF1), (NTEOPF(2),NTEOPF2),
     %                 (NTEOPF(3),NTEOPF3), (NTEOPF(4),NTEOPF4)
      INTEGER           NOSOFATE(3,4)
      DATA              NOSOFATE / 1,3,2,  2,3,4,  3,1,4,  4,1,2 /
C     => LA NORMALE A LA FACE EST ORIENTEE VERS L'EXTERIEUR DU TETRAEDRE
      INTEGER           NOSOARTE(2,6)
      DATA              NOSOARTE/ 1,2, 2,3, 3,1, 1,4, 2,4, 3,4 /

      TRACTE0 = TRACTE
      LORBITE = 1


ccc      if( nbsomm0 .eq. 172078 ) then
ccc         print *,'tetreto: debut avec nbsomm0=',nbsomm0
ccc         nte=865332
ccc         print*,'tetreto: notetr(',nte,'): st',(notetr(k,nte),k=1,4),
ccc     %          ' teop',(notetr(k,nte),k=5,8)
ccc         nte=865333
ccc         print*,'tetreto: notetr(',nte,'): st',(notetr(k,nte),k=1,4),
ccc     %          ' teop',(notetr(k,nte),k=5,8)
ccc         tracte = .true.
ccc      endif


C     VERIFICATIONS SUR LES NBTECFPR TETRAEDRES DU TABLEAU NOTECFPR
C     -------------------------------------------------------------
ccc      PRINT 10000
ccc10000 FORMAT( 200('_') )

      PRINT *,'tetreto: Debut avec',NBSSET0,
     %' SOUS-ETOILES A TETRAEDRISER VMOYEN=',VMOYEN,' QTEAME=',QTEAME
      PRINT *,'tetreto:',NBTECFPR,' NO NOTETR des TETRAEDRES INITIAUX'
ccc      PRINT 10001,(NOTECFPR(K),K=1,NBTECFPR)
ccc10001 FORMAT( 20I9 )

ccc      DO K = 1, NBTECFPR
ccc         NTE = NOTECFPR(K)
ccc         PRINT*,'tetreto: TETRAEDRE(',NTE,') St:',(NOTETR(L,NTE),L=1,4)
ccc     %         ,' NteOp:',(NOTETR(L,NTE),L=5,8)
ccc      ENDDO

C     VERIFIER L'ABSENCE DE TETRAEDRE  OPPOSE  NEGATIF (INCONNU)
C     VERIFIER L'ABSENCE DE TETRAEDRES OPPOSES DOUBLES
C     VERIFIER L'ABSENCE DE FACE COMMUNE A AU MOINS 3 TETRAEDRES
C     VERIFIER L'OPPOSITION DES TETRAEDRES 
C     DANS LA LISTE DES NBTECFPR TETRAEDRES NOTECFPR
      CALL VEOPTE( NBTECFPR, NOTECFPR, NOTETR, PTXYZD, NBFANR )

C     TRACE DES NBSSET0 SOUS-ETOILES INITIALES
C     ----------------------------------------
      CALL TRFETO7( PTXYZD, 0, NFETOI, NBSSET0, N1SSET, 0, NSARP2F )

C     NOMBRE DE TRAITEMENTS DE L'ETOILE
      NOPASS6 = 0
      IERR    = 0
      MAXTETR = MIN( MAX(6,NBTRCF) * NBTECFPR, MXTECR )
ccc      MAXSOMM = 16 * NBTRCF
      MAXSOMM = 16

      NS4MAX00= -1
      NS4MAX0 = 0

      NBSOMM = NBSOMM0

C     NOMBRE DE SOUS-ETOILES DE L'ETOILE INITIALE
      NBSSET = NBSSET0

C     NOMBRE INITIAL DE TETRAEDRES CREES DANS L'ETOILE INITIALE
      NBTECR = 0

C     NOMBRE DE FACES AVEC UN TETRAEDRE OPPOSE INCONNU
      NBFANR = 0

C     VOLUME DES NBTECR TETRAEDRES CREES
      VOLTECR = 0D0

C     ======================================================================
C     TETRAEDRISATION DE LA PILE DES NBSSET0 SOUS-ETOILES DE L'ETOILE NFETOI
C     ======================================================================

C     TETRAEDRISATION DE LA SOUS-ETOILE NBSSET EN HAUT DE PILE N1SSET(NBSSET)
C     -----------------------------------------------------------------------
 5    IF( NBSSET .LE. 0 ) GOTO 9000

cccC     VERIFIER SI TOUT TETRAEDRE VIDE A SON 1-ER NUMERO NOTETR NUL
ccc      CALL VETEVIDE( N1TEVI, MXTETR, NOTETR, IERR )

C     RAPPORT DU VOLUME ACTUEL SUR LE VOLUME MOYEN D'UNE ETOILE DE TETRAEDRES
      RAVEVM = 1D0
      QUVEVM = QUAMINEX

C     NOMBRE DE TETRAEDRES CREES ACTUELLEMENT AVANT TRAITEMENT
C     DE LA SOUS-ETOILE NBSSET
      NBTECR0 = NBTECR
      MODIFS  = 0
      NOCHOIXYZ = 0
      NBVPSI = 0

C     1-ERE FACE NFETOI VERSION 2 DE LA SOUS-ETOILE NBSSET DE NFETOI A TETRAEDRISER
      N1FEOC = N1SSET( NBSSET )

C     VERIFIER QUE TOUTE ARETE DE LA SOUS ETOILE NBSSET APPARTIENT
C     SEULEMENT A 2n FACES DE L'ETOILE avec n=1,2
      CALL VETAFET( NOTETR, N1FEOC, NFETOI, NBFETO, NBARPB )

      IF( NBARPB .GT. 0 ) THEN
C        TRACE DE LA SOUS-ETOILE N1FEOC A TETRAEDRISER
         TRACTE = .TRUE.
         KTITRE='tetreto: ATTENTION LA SOUS-ETOILE A TETRAEDRISER avec A
     %RETES A PROBLEME'
         CALL SANSDBL( KTITRE, L )
         PRINT*,KTITRE(1:L),' NoSSET=',NBSSET,' NBARPB=',NBARPB
         CALL TRFETO4( KTITRE(1:L), PTXYZD, N1FEOC, NFETOI,
     %                 NBTRCF, NOTRCF, LEFACO, NO0FAR,
     %                 NBTECR, NOTECR, NOTETR )
         IERR = 10
         GOTO 9999
      ENDIF

C     NBFETO NOMBRE DE FACES DE LA SOUS-ETOILE NBSSET A TETRAEDRISER
      NBFETO = 0
      NF00   = 0
      NF0    = N1FEOC
 10   IF( NF0 .GT. 0 ) THEN
C        UNE FACE DE PLUS
         NBFETO = NBFETO + 1
C        NO DU 1-ER SOMMET DE LA FACE EST RENDU POSITIF
         NFETOI(2,NF0) = ABS( NFETOI(2,NF0) )
C        LA FACE SUIVANTE DE NF0 DANS NFETOI
         NF00 = NF0
         NF0  = NFETOI(5,NF0)
         GOTO 10
      ENDIF

ccc      PRINT *
ccc      PRINT *,'tetreto: Tetraedrisation de la SOUS-ETOILE',NBSSET
ccc     %       ,' de',NBFETO,' FACES. Face initiale N1FEOC=',N1FEOC

      N = NBSOMM - NBSOMM0
      IF( NBTECR .GT. MAXTETR .OR. N .GT. MAXSOMM )THEN
         print*,'?????????????????????????????????????????????????????'
         print*,'PB tetreto:    TROP DE FACES dans l''ETOILE',NBFETO
         print*,'PB tetreto: ou TROP DE TETRAEDRES CREES',NBTECR
         print*,'PB tetreto: ou TROP DE SOMMETS    CREES',N
         print*,'PB tetreto: NOMBRE DE TETRAEDRES INITIAUX',NBTECFPR
         print*,'?????????????????????????????????????????????????????'
         IERR = 8
         GOTO 9999
      ENDIF

C     LA SOUS-ETOILE EST DEPILEE CAR EN COURS DE TETRAEDRISATION
C     ET LA TETRAEDRISATION PEUT A SON TOUR CREER DES SOUS-ETOILES...
C     LA PREMIERE FACE SIMPLE DE L'ETOILE EST N1FEOC DANS NFETOI
      NBSSET = NBSSET - 1

C     ----------------------------------------------------------------
C     MODE DE CREATION DES TETRAEDRES SELON LE NOMBRE DE FACES SIMPLES
C     ----------------------------------------------------------------
      IF( NBFETO .LE. 0 ) THEN
C        PLUS AUCUNE FACE DANS LA SOUS-ETOILE NBSSET+1
C        PASSAGE A LA SOUS-ETOILE EN HAUT DE PILE
         GOTO 5
      ENDIF

      IF( NBFETO .LE. 3 ) THEN

C        ----------------------------------------------
C        MOINS DE 4 FACES DANS LA SOUS-ETOILE -> ERREUR
C        ----------------------------------------------
         IF( NBFETO .EQ. 2 ) THEN
C           CES 2 FACES SONT ELLES IDENTIQUES? (CAS 1 SEULE DEMI-ETOILE)
            NF0 = N1FEOC
            NF1 = NFETOI(5,NF0)
            DO K=1,3
               NOSOTR2(K) = NFETOI(1+K,NF0)
               NOSOTR3(K) = NFETOI(1+K,NF1)
            ENDDO
            CALL TRI3NO( NOSOTR2, NOSOTR2 )
            CALL TRI3NO( NOSOTR3, NOSOTR3 )
            IF( NOSOTR2(1) .EQ. NOSOTR3(1) .AND.
     %          NOSOTR2(2) .EQ. NOSOTR3(2) .AND.
     %          NOSOTR2(3) .EQ. NOSOTR3(3) ) THEN
               PRINT *,'tetreto: Attention SOUS-ETOILE',NBSSET+1,
     %                 ' AVEC 2 FACES IDENTIQUES RETIREES', NOSOTR2
               GOTO 5
            ENDIF
         ENDIF

         PRINT*
         PRINT*,'tetreto: Probleme NOMBRE DE FACES=',NBFETO,
     %          ' INCORRECT POUR FORMER UN TETRAEDRE!!!'
         KTITRE='PB tetreto: NBFETO=       '
         WRITE(KTITRE(21:24),'(I4)') NBFETO
         DO K=1,NBSSET+1
            PRINT*
            PRINT*,'SOUS ETOILE',K
            NF = N1SSET( K )
            CALL AFNFETOI( NF, NFETOI, NOTETR )
         ENDDO

C        TRACE DE LA SOUS-ETOILE NBSSET+1 A PROBLEME ISSUE DE N1FEOC
         TRACTE = .TRUE.
         KTITRE='tetreto: ATTENTION LA SOUS-ETOILE       N''A PAS ASSEZ 
     %de FACES       '
         WRITE(KTITRE(35:37),'(I3)') NBSSET+1
         WRITE(KTITRE(65:69),'(I5)') NBFETO
         CALL SANSDBL( KTITRE, L )
         PRINT*,KTITRE(1:L)
         CALL TRFETO4( KTITRE(1:L), PTXYZD, N1FEOC, NFETOI,
     %                 NBTRCF, NOTRCF, LEFACO, NO0FAR,
     %                 NBTECR, NOTECR, NOTETR )

C        TENTATIVE DE JOINDRE CETTE SOUS-ETOILE INCORRECTE A
C        LA SOUS ETOILE QUI LA PRECEDE
         IF( NBSSET .GT. 0 ) THEN
            PRINT*,'tetreto: JONCTION des SOUS-ETOILES',NBSSET,NBSSET+1
C           LA DERNIERE FACE DE NBSSET A POUR SUIVANTE
C           LA PREMIERE DE NBSSET+1
            NFETOI(5,NF00) = N1SSET( NBSSET )
            N1SSET( NBSSET ) = N1FEOC
            GOTO 5
         ENDIF 

C        FIN DE SOUS ETOILE NBSSET+1 AVEC ERREUR
         IERR = 4
         GOTO 9010

      ENDIF

      IF( NBFETO .EQ. 4 ) THEN
C        ---------------------------------------------------------
C        LES 4 SEULES FACES DE LA SOUS-ETOILE DOIVENT FORMER UN
C        TETRAEDRE SINON ERREUR
C        FIN de CREATION des TETRAEDRES de la SOUS-ETOILE NBSSET+1
C        ---------------------------------------------------------
C        TETRAEDRE OPPOSE A LA FACE N1FEOC
         NTEOPF(1) = NFETOI(1,N1FEOC)

C        LE NUMERO PTXYZD DES 4 SOMMETS DU TETRAEDRE NOSOTE
         NS1 = NFETOI(2,N1FEOC)
         NS2 = NFETOI(3,N1FEOC)
         NS3 = NFETOI(4,N1FEOC)
C        LE 4-EME SOMMET EST CELUI DES FACES DIFFERENT DES 3 DE LA FACE 1
C        LA FACE 2 DE LA SOUS-ETOILE DE 4 FACES
         NF2 = NFETOI(5,N1FEOC)
         DO K=1,3
            NS4 = NFETOI(1+K,NF2)
            IF( NS4 .NE. NS1 .AND. NS4 .NE. NS2 .AND. NS4 .NE. NS3 )
     %          GOTO 12
         ENDDO

C        RECHERCHE DES TETRAEDRES OPPOSES AUX 3 DERNIERES FACES
C        VERIFICATION: LES 3 DERNIERES FACES DOIVENT AVOIR LEURS 3
C        SOMMETS PARMI LES 4 DU TETRAEDRE
 12      DO 20 K=2,4

C           LES 3 SOMMETS DE LA FACE K DU TETRAEDRE NS1 NS2 NS3 NS4
            DO L=1,3
               NOSOTR2(L) = NOSOTE( NOSOFATE(L,K)  )
            ENDDO
C           TRI CROISSANT DES 3 NUMEROS DE SOMMETS DE LA FACE K
            CALL TRI3NO( NOSOTR2, NOSOTR2 )

            NF = NF2
C           LES 3 SOMMETS DE LA FACE NF DE LA SOUS-ETOILE
 16         DO L=1,3
               NOSOTR3(L) = NFETOI(1+L,NF)
            ENDDO
C           TRI CROISSANT DES 3 NUMEROS DE SOMMETS DE LA FACE NF
            CALL TRI3NO( NOSOTR3, NOSOTR3 )

            IF( NOSOTR2(1) .EQ. NOSOTR3(1) .AND.
     %          NOSOTR2(2) .EQ. NOSOTR3(2) .AND.
     %          NOSOTR2(3) .EQ. NOSOTR3(3) ) THEN
C               LA FACE K DU TETRAEDRE EST LA FACE NF DE LA SOUS-ETOILE
C               LE NO DU TETRAEDRE OPPOSE A LA FACE K DU TETRAEDRE
                NTEOPF(K) = NFETOI(1,NF)
                GOTO 20
            ENDIF

C           FACE SUIVANTE
            NF = NFETOI(5,NF)
            IF( NF .NE. 0 ) GOTO 16

            PRINT*,'PB tetreto: FACE',K,'DU TETRAEDRE de St:',NOSOTE,
     %' N''EST PAS UNE DES 4 SEULES FACES de l''ETOILE => ETOILE INCORRE
     %CTE'
            IERR  = 5
            NOCAS = 0
            KTITRE='PB tetreto: FACE NON TROUVEE DANS LE TETRAEDRE    '
            WRITE( KTITRE(49:92),'(4I11)') (NOSOTE(L),L=1,4)
            GOTO 9010

 20      ENDDO

C        VERIFICATION DE TETRAEDRES OPPOSES DOUBLES OU TRIPLES ...
         DO K=1,3
            NTEOP = NTEOPF(K)
            IF( NTEOP .GT. 0 ) THEN
               DO L=K+1,4
                  IF( NTEOPF(L) .EQ. NTEOP ) THEN
                  PRINT*,'PB tetreto: Pour les 4 DERNIERES FACES de l''E
     %TOILE de SOMMETS',NOSOTE
                  PRINT*,'PB tetreto: TETRAEDRES OPPOSES',NTEOPF,
     %                   ' avec des DOUBLES...'
                  IERR  = 5
                  NOCAS = 0
                  KTITRE='PB tetreto: TETRAEDRES OPPOSES MULTIPLES     '
                  WRITE( KTITRE(42:85),'(4I11)') NTEOPF
                  GOTO 9010
                  ENDIF
               ENDDO
            ENDIF
         ENDDO

C        AJOUT DU TETRAEDRE FORME DES 4 SEULES FACES DE LA SOUS-ETOILE
C        LES 3 SOMMETS DE LA FACE N1FEOC SONT DANS LE SENS DIRECT
C        DU TETRAEDRE CAR SA NORMALE EST DIRIGEE VERS L'INTERIEUR

ccc         NS1F0 = NFETOI(2,N1FEOC)
ccc         NS2F0 = NFETOI(3,N1FEOC)
ccc         NS3F0 = NFETOI(4,N1FEOC)
ccc         CALL AJTEET( NS1F0,   NS2F0,   NS3F0,   NS4,

         CALL AJTEET( NS1,     NS2,     NS3,     NS4,
     %                NTEOPF1, NTEOPF2, NTEOPF3, NTEOPF4,
     %                NBSOMM,  PTXYZD,  QUVEVM,
     %                MXTETR,  N1TEVI,  NOTETR,  NUDTETR,
     %                N1TETS,  MXTECR,  NBTECR,  NOTECR,  VOLTECR,
     %                VOLTEC,  QUATET,  IERR )
         IF( IERR .NE. 0 ) GOTO 9999

C        LE TETRAEDRE NBTECR EST CREE
         NTE = NOTECR( NBTECR )

C        LE VOLUME et LA QUALITE DU TETRAEDRE NTE
         CALL QUATETD( PTXYZD( 1, NOTETR(1,NTE) ),
     %                 PTXYZD( 1, NOTETR(2,NTE) ),
     %                 PTXYZD( 1, NOTETR(3,NTE) ),
     %                 PTXYZD( 1, NOTETR(4,NTE) ),
     %                 ARMIN, ARMAX, SURFTR, VOLTEC, QUATET )

         IF( VOLTEC .LT. 0D0 ) THEN

            PRINT*,'tetreto: NOTETR(',NTE,'): St',(NOTETR(k,NTE),k=1,4),
     %             ' Tetra Op',(NOTETR(k,NTE),k=5,8),
     %             ' V=',VOLTEC,' Q=',QUATET
            PRINT*,'tetreto:',(' NPSOFR(',NOTETR(k,NTE),')=',
     %              NPSOFR(NOTETR(k,NTE)),k=1,4)

C           PROBLEME: TRACE DES ARETES DU TETRAEDRE NTE, DE SES 4 TETRAEDRES
C           VOISINS, ET DES TETRAEDRES VOISINS DES TETRAEDRES VOISINS
ccc            tracte = .true.
      KTITRE='tetreto: 4FACES DONNENT le TETRAEDRE          de QUALITE   
     %                 et VA ETRE TRAITE'
            WRITE(KTITRE(38:45),'(I8)') NTE
            WRITE(KTITRE(58:71),'(G14.6)') QUATET
            CALL TRFETO9( KTITRE, PTXYZD, NTE, NOTETR )
            tracte = tracte0

C           TENTATIVE DE DEPLACER LE SOMMET CENTRAL DU TETRAEDRE PLAT NTEQM
C           DE L'AUTRE COTE DE LA FACE POUR RENDRE LE VOLUME POSITIF
C           ---------------------------------------------------------------
            CALL DEPSCTVN( NTE, NOTETR, N1TETS, PTXYZD, NPSOFR, INFACO,
     %                     VOLTEC, QUATET, MODIF )

            IF( MODIF .NE. 0 ) THEN
C              LES XYZ D'UN SOMMET DU TETRAEDRE NTE ONT ETE MODIFIES
C              POUR RENDRE LE VOLUME ET LA QUALITE POSITIFS
               MODIFS = MODIFS + MODIF
            ENDIF

         ENDIF

         NOCAS = 1
         print*,'tetreto: NOCAS=',NOCAS,': 4Faces=> +TETRAEDRE',NTE,'='
     %    ,(NOTETR(KK,NTE),KK=1,8),' N1TEVI=',N1TEVI,
     %     ' V=',VOLTEC,' Q=',QUATET

         IF( VOLTEC .LE. 0D0 .OR. QUATET .LE. 0 ) THEN

C           PROBLEME: LE TRAITEMENT CI-DESSUS N'A PAS MODIFIE LE VOLUME
C           TRACE DES ARETES DU TETRAEDRE NTE, DE SES 4 TETRAEDRES
C           VOISINS, ET DES TETRAEDRES VOISINS DES TETRAEDRES VOISINS
ccc            tracte = .true.
      KTITRE='tetreto: 4FACES DONNENT le TETRAEDRE          de QUALITE   
     %                VOLUME=               MALGRE le TRAITEMENT!'
            WRITE(KTITRE(38:45),'(I8)') NTE
            WRITE(KTITRE(58:71),'(G14.6)') QUATET
            WRITE(KTITRE(80:93),'(G14.6)') VOLTEC
            CALL TRFETO9( KTITRE, PTXYZD, NTE, NOTETR )
            tracte = tracte0

         ENDIF

C        AJOUTER EVENTUELLEMENT LES 4 FACES DU TETRAEDRE NTE
C        DANS LE TABLEAU LEFACO
         CALL AJTELEFA( NTE, NOTETR, INFACO, MXFACO, LEFACO )

C        UN TETRAEDRE OPPOSE AU TETRAEDRE NTE EST IL DOUBLE?
         CALL TETOPDOU( NTE, NOTETR, NF1, NF2 )
         GOTO 800

      ENDIF


C     ----------------------------------------------------------------------------
C     1) LA SOUS-ETOILE ISSUE DE N1FEOC A PLUS DE 4 FACES=TRIANGLES SIMPLES NFETOI
C        ESSAI DE TETRAEDRISER LES TRIANGLES DE LA SOUS-ETOILE SANS AJOUT DE POINT
C        RECHERCHE DANS LA SOUS-ETOILE DE 3 FACES ADJACENTES FORMANT UN TETRAEDRE
C        DE VOLUME>0 ET D'ARETES SANS INTERSECTION AVEC LES AUTRES FACES
C     ----------------------------------------------------------------------------
      NF0 = N1FEOC
 51   IF( NF0 .GT. 0 ) THEN

C        BOUCLE SUR LES 3 ARETES DE LA FACE NF0
         DO K=1,3

C           L'ARETE K DE LA FACE NF0
C           NUMERO DES 3 SOMMETS AVEC LES 2 DE L'ARETE EN PREMIER
            IF( K .EQ. 3 ) THEN
               KK = 1
            ELSE
               KK = K + 1
            ENDIF
C           LES 2 SOMMETS DE L'ARETE K DE LA FACE NF0
            NS1F0 = NFETOI(1+K ,NF0)
            NS2F0 = NFETOI(1+KK,NF0)
C           LE SOMMET 3 DE LA FACE NF0
            IF( KK .EQ. 3 ) THEN
               KKK = 1
            ELSE
               KKK = KK + 1
            ENDIF
            NS3F0 = NFETOI(1+KKK,NF0)

C           RECHERCHE DE LA 2-EME FACE DE LA SOUS-ETOILE NON NF0 et D'ARETE K
            NF1 = N1FEOC
 52         IF( NF1 .GT. 0 ) THEN
               IF( NF1 .NE. NF0 ) THEN
                  DO L=1,3
C                    L'ARETE L DE LA FACE NF1
                     IF( L .EQ. 3 ) THEN
                        LL = 1
                     ELSE
                        LL = L + 1
                     ENDIF
C                    LES 2 SOMMETS DE L'ARETE L DE LA FACE NF1
                     NS1F1 = NFETOI(1+L ,NF1)
                     NS2F1 = NFETOI(1+LL,NF1)
                     IF( NS1F0 .EQ. NS2F1 .AND. NS2F0 .EQ. NS1F1 ) THEN
ccc     %          .OR. NS1F0 .EQ. NS1F1 .AND. NS2F0 .EQ. NS2F1 ) THEN
C                       L'ARETE K DE NF0 EST L'ARETE L DE NF1
C                       LE SOMMET 3 DE LA FACE NF1
                        IF( LL .EQ. 3 ) THEN
                           LLL = 1
                        ELSE
                           LLL = LL + 1
                        ENDIF
                        NS3F1 = NFETOI(1+LLL,NF1)

                        IF( NS3F1 .EQ. NS3F0 ) THEN
C                          CAS DE 2 FACES DU CF IDENTIQUES DANS L'ENCOCHE
                       print*,'tetreto: 2FACES de l''etoile avec 3 MEME
     %S SOMMETS (ENCOCHE POSSIBLE) NF0=',NF0,' ST:',NS1F0,NS2F0,NS3F0,
     %                         ' NF1=',NF1,' ST:',NS1F1,NS2F1,NS3F1
                           GOTO 69
                        ENDIF

C                       LE VOLUME FACE NF0-SOMMET NS3F1 EST IL POSITIF?
C                       SON VOLUME ET SA QUALITE
                        CALL QUATETD( PTXYZD(1,NFETOI(2,NF0)),
     %                                PTXYZD(1,NFETOI(3,NF0)),
     %                                PTXYZD(1,NFETOI(4,NF0)),
     %                                PTXYZD(1,NS3F1),
     %                                ARMIN, ARMAX, SURFTR, V, Q )

                        IF( Q .LE. QUVEVM ) THEN
C                          TETRAEDRE DE VOLUME<=0 ou QUALITE TROP FAIBLE
C                           => RECHERCHE D'UNE AUTRE FACE NF1
                           GOTO 69
                        ENDIF

C                       RECHERCHE DE LA 3-eme FACE NON NF0 et NON NF1
                        NF2 = N1FEOC
 53                     IF( NF2 .GT. 0 ) THEN
                           IF( NF2 .NE. NF0 .AND. NF2 .NE. NF1 ) THEN

                              DO 65 M=1,3

C                                LES 2 SOMMETS DE L'ARETE M DE LA FACE NF2
                                 NS1F2 = NFETOI(1+M ,NF2)
                                 IF( NS1F2 .NE. NS3F1 ) GOTO 65

C                                SOMMETS NS1F2=NS3F1
                                 IF( M .EQ. 3 ) THEN
                                    MM = 1
                                 ELSE
                                    MM = M + 1
                                 ENDIF
                                 NS2F2 = NFETOI(1+MM,NF2)

                                 IF( MM .EQ. 3 ) THEN
                                    MMM = 1
                                 ELSE
                                    MMM = MM + 1
                                 ENDIF
C                                LE SOMMET 3 DE LA FACE NF2
                                 NS3F2 = NFETOI(1+MMM,NF2)

                                 IF( NS2F2 .EQ. NS1F0 ) THEN
                                    IF( NS3F2 .EQ. NS3F0 ) THEN

C                                      NOCAS 2:  TETRAEDRE NF0-NS3F1 avec
C                                      NS1F2=NS3F1 et NS2F2=NS1F0 et NS3F2=NS3F0
C                                      -----------------------------------------
                                       NOCAS = 2

C                                      UNE DES 3 ARETES DE NF2 INTERSECTE T ELLE
C                                      LA FACE NS3F0 NS2F0 NS3F1?
                                       N0 = NFETOI(4,NF2)
                                       DO N=1,3
                                          N1 = NFETOI(1+N,NF2)
                                          CALL INARTR( PTXYZD(1,N0),
     %                                                 PTXYZD(1,N1),
     %                                                 PTXYZD(1,NS3F0),
     %                                                 PTXYZD(1,NS2F0),
     %                                                 PTXYZD(1,NS3F1),
     %                                                 LINTER,XYZ,CBTR)
                                          IF( LINTER  .EQ. 1 .AND.
     %                                     CBTR(1) .LT. 0.999999D0 .AND.
     %                                     CBTR(2) .LT. 0.999999D0 .AND.
     %                                     CBTR(3) .LT. 0.999999D0 )THEN
C                                          OUI: UN POINT D'INTERSECTION
C                                          ABANDON DE CE CHOIX DE 2 FACES
                           print *,'tetreto: bizarre nocas=2 cbtr=',cbtr
                                           GOTO 69
                                          ENDIF
                                          N0 = N1
                                       ENDDO

C                                      RECUPERER LE TETRAEDRE OPPOSE A LA
C                                      FACE NS2F0-NS3F0-NS3F1
C                                      DU TETRAEDRE NS1F0+NS2F0+NS3F0+NS3F1
                                       NTEOPF( 2 )= -1
                                       NOSOTR2(1) = NS2F0
                                       NOSOTR2(2) = NS3F0
                                       NOSOTR2(3) = NS3F1
                                       CALL TRI3NO( NOSOTR2, NOSOTR2 )

                                       NF3 = N1FEOC
 54                                    IF( NF3 .GT. 0 ) THEN
                                      CALL TRI3NO(NFETOI(2,NF3),NOSOTR3)

                                       IF(NOSOTR2(1).EQ.NOSOTR3(1) .AND.
     %                                    NOSOTR2(2).EQ.NOSOTR3(2) .AND.
     %                                    NOSOTR2(3).EQ.NOSOTR3(3) )THEN
C                                         TETRAEDRE OPPOSE A LA FACE 2
C                                         DU TETRAEDRE NS1F0-NS2F0-NS3F0-NS3F1
                                          NTEOPF( 2 ) = NFETOI( 1, NF3 )
                                          GOTO 55
                                       ENDIF
                                       NF3 = NFETOI(5,NF3)
                                       GOTO 54
                                       ENDIF

C                                      TETRAEDRE OPPOSE A LA FACE 3
C                                      DU TETRAEDRE NS1F0-NS2F0-NS3F0-NS3F1
 55                                    NTEOPF( 3 ) = NFETOI( 1, NF2 )
                                       NTEOPF( 4 ) = NFETOI( 1, NF1 )
                                       GOTO 60
                                    ENDIF
                                 ENDIF

C                                AUTRE POSSIBILITE
                                 IF( NS2F2 .EQ. NS3F0 ) THEN
                                    IF( NS3F2 .EQ. NS2F0 ) THEN

C                                      NOCAS 3:  TETRAEDRE NF0-NS3F1 avec
C                                      NS1F2=NS3F1 et NS2F2=NS3F0 et NS3F2=NS2F0
C                                      -----------------------------------------
                                       NOCAS = 3

C                                      UNE DES 3 ARETES DE NF2 INTERSECTE T ELLE
C                                      LA FACE NS3F0 NS3F1 NS1F0?
                                       N0 = NFETOI(4,NF2)
                                       DO N=1,3
                                          N1 = NFETOI(1+N,NF2)
                                          CALL INARTR( PTXYZD(1,N0),
     %                                                 PTXYZD(1,N1),
     %                                                 PTXYZD(1,NS3F0),
     %                                                 PTXYZD(1,NS3F1),
     %                                                 PTXYZD(1,NS1F0),
     %                                                 LINTER,XYZ,CBTR)
                                          IF( LINTER  .EQ. 1 .AND.
     %                                      CBTR(1).LT.0.999999D0 .AND.
     %                                      CBTR(2).LT.0.999999D0 .AND.
     %                                      CBTR(3).LT.0.999999D0)THEN
C                                           OUI: UN POINT D'INTERSECTION
C                                           ABANDON DE CE CHOIX DE 2 FACES
                           print *,'tetreto: bizarre nocas=3 cbtr=',cbtr
                                            GOTO 69
                                          ENDIF
                                          N0 = N1
                                       ENDDO

C                                      RECUPERER LE TETRAEDRE OPPOSE A LA
C                                      FACE NS1F0-NS3F0-NS3F1
C                                      DU TETRAEDRE NS1F0+NS2F0+NS3F0+NS3F1
                                       NTEOPF( 3) = -1
                                       NOSOTR2(1) = NS1F0
                                       NOSOTR2(2) = NS3F0
                                       NOSOTR2(3) = NS3F1
                                       CALL TRI3NO( NOSOTR2, NOSOTR2 )

                                       NF3 = N1FEOC
 57                                    IF( NF3 .GT. 0 ) THEN
                                      CALL TRI3NO(NFETOI(2,NF3),NOSOTR3)
                                      IF(NOSOTR2(1) .EQ. NOSOTR3(1).AND.
     %                                   NOSOTR2(2) .EQ. NOSOTR3(2).AND.
     %                                   NOSOTR2(3) .EQ. NOSOTR3(3))THEN
                                         NTEOPF(3) = NFETOI( 1, NF3 )
                                         GOTO 59
                                      ENDIF
                                          NF3 = NFETOI(5,NF3)
                                          GOTO 57
                                       ENDIF

C                                      TETRAEDRE OPPOSE A LA FACE 3
 59                                    NTEOPF( 2 ) = NFETOI( 1, NF2 )
                                       NTEOPF( 4 ) = NFETOI( 1, NF1 )
                                       GOTO 60
                                    ENDIF
                                 ENDIF
                                 GOTO 65

C                                TETRAEDRE OPPOSE PAR LA 4-EME FACE DE NTE
 60                              NTEOPF(4) = NFETOI(1,NF1)

C                                LES 3 FACES NF0 NF1 NF2 FORMENT UN TETRAEDRE
C                                NS1F0-NS2F0-NS3F0-NS3F1 A AJOUTER
C                                CREATION DU TETRAEDRE NS1F0-NS2F0-NS3F0-NS3F1
                                 NFETOI(2,NF0) = NS1F0
                                 NFETOI(3,NF0) = NS2F0
                                 NFETOI(4,NF0) = NS3F0
                                 NTEOPF1 = NFETOI(1,NF0)

                                 CALL AJTEET(NS1F0, NS2F0, NS3F0, NS3F1,
     %                               NTEOPF1, NTEOPF2, NTEOPF3, NTEOPF4,
     %                               NBSOMM,  PTXYZD, QUVEVM,
     %                               MXTETR,  N1TEVI, NOTETR,
     %                               NUDTETR, N1TETS,
     %                               MXTECR,  NBTECR, NOTECR, VOLTECR,
     %                               VOLTEC,  QUATET, IERR )
                                 IF( IERR .NE. 0 ) GOTO 9999

                                 NTE = NOTECR( NBTECR )

C                                AJOUTER EVENTUELLEMENT LES 4 FACES DU
C                                TETRAEDRE NTE DANS LE TABLEAU LEFACO
                                 CALL AJTELEFA( NTE,    NOTETR,
     %                                          INFACO, MXFACO, LEFACO )

                                 print *,'tetreto: NOCAS=',NOCAS,
     %                           ': 3Faces=> +TETRAEDRE',NTE,'=',
     %                           (NOTETR(KK,NTE),KK=1,8),' N1TEVI=',
     %                            N1TEVI,' V=',VOLTEC,' Q=',QUATET

C                                EXAMEN DE LA SOUS-ETOILE
C                                APRES AJOUT DU TETRAEDRE NTE
                                 GOTO 800

 65                           CONTINUE
                           ENDIF

C                          LA FACE SUIVANTE DE NF2
                           NF2 = NFETOI(5,NF2)
                           GOTO 53

                        ENDIF
                     ENDIF
                  ENDDO
               ENDIF

C              LA FACE SUIVANTE DE NF1
 69            NF1 = NFETOI(5,NF1)
               GOTO 52
            ENDIF

         ENDDO

C        LA FACE SUIVANTE DE NF0 DANS NFETOI
         NF0 = NFETOI(5,NF0)
         GOTO 51

      ENDIF


C     --------------------------------------------------------------------
C     2) RECHERCHE D'UN POINT XYZ SOMMET D'UNE GRILLE UN CUBE DE CENTRE
C        LE BARYCENTRE DU BARYCENTRE DES FACES DE LA SOUS-ETOILE
C        ET DES FACES DU CONTOUR FERME QUI JOINT AUX FACES SIMPLES DONNENT
C        TOUTES DES TETRAEDRES DE VOLUME>0 ET SANS INTERSECTION AVEC ELLES
C     --------------------------------------------------------------------
      CALL SIPTGRBA( KTITRE, PTXYZD, N1FEOC, NFETOI,
     %               NBFETO, XYZ,    MXVPSI, NBVPSI, NFVPSI,
     %               QUAMIN, RAVEVM )

C     RAVEVM  : RAPPORT DU VOLUME DE LA SOUS-ETOILE ACTUELLE AU VOLUME DE
C               NBFETO FOIS LE VOLUME MOYEN D'UN TETRAEDRE
C               REDUIT A 1 SI PLUS GRAND QUE 1
      RAVEVM = MIN( 1D0, RAVEVM )

cccC     SI LA SOUS-ETOILE EST ECRASEE ALORS LA QUALITE MINIMALE DES TETRAEDRES
cccC     A CREER EST REDUITE D'AUTANT
ccc      QUVEVM = REAL( QUVEVM * RAVEVM )

      IF( NBVPSI .EQ. NBFETO ) THEN
         NOCHOIXYZ = 1
C        TETRAEDRISATION DE TOUTES LES FACES DE LA SOUS-ETOILE ETOILEES PAR XYZ
         PRINT*,'tetreto: NBFETO=NBVPSI=',NBFETO,
     %   ' TOUTES les FACES de la SOUS-ETOILE sont ETOILEES PAR XYZ',XYZ
         GOTO 500
      ENDIF


C     ----------------------------------------------------------------------------
C     3) ESSAI DE TETRAEDRISER LES TRIANGLES DE LA SOUS-ETOILE SANS AJOUT DE POINT
C        RECHERCHE DANS LA SOUS-ETOILE DE 2 FACES ADJACENTES FORMANT UN TETRAEDRE
C        DE VOLUME>0 ET DONT L'ARETE OPPOSEE INTERSECTE AUCUNE DES FACES
C     ----------------------------------------------------------------------------
      NF0 = N1FEOC
 150  IF( NF0 .GT. 0 ) THEN

C        BOUCLE SUR LES 3 ARETES DE LA FACE NF0
         DO K=1,3

C           L'ARETE K DE LA FACE NF0
            IF( K .EQ. 3 ) THEN
               KK = 1
            ELSE
               KK = K + 1
            ENDIF
C           LES 2 SOMMETS DE L'ARETE K DE LA FACE NF0
            NS1F0 = NFETOI(1+K ,NF0)
            NS2F0 = NFETOI(1+KK,NF0)

C           LE SOMMET 3 DE LA FACE NF0
            IF( KK .EQ. 3 ) THEN
               KKK = 1
            ELSE
               KKK = KK + 1
            ENDIF
            NS3F0 = NFETOI(1+KKK,NF0)

C           RECHERCHE DE LA 2-EME FACE NF1 DE LA SOUS-ETOILE NON NF0 et D'ARETE K
            NF1 = N1FEOC
 152        IF( NF1 .GT. 0 ) THEN
               IF( NF1 .NE. NF0 ) THEN
                  DO L=1,3
C                    L'ARETE L DE LA FACE NF1
                     IF( L .EQ. 3 ) THEN
                        LL = 1
                     ELSE
                        LL = L + 1
                     ENDIF
C                    LES 2 SOMMETS DE L'ARETE L DE LA FACE NF1
                     NS1F1 = NFETOI(1+L ,NF1)
                     NS2F1 = NFETOI(1+LL,NF1)
                     IF( NS1F0 .EQ. NS2F1 .AND. NS2F0 .EQ. NS1F1 ) THEN
c     %          .OR. NS1F0 .EQ. NS1F1 .AND. NS2F0 .EQ. NS2F1 ) THEN
C                       L'ARETE K DE NF0 EST L'ARETE L DE NF1
C                       LE SOMMET 3 DE LA FACE NF1
                        IF( LL .EQ. 3 ) THEN
                           LLL = 1
                        ELSE
                           LLL = LL + 1
                        ENDIF
                        NS3F1 = NFETOI(1+LLL,NF1)

                        IF( NS3F1 .EQ. NS3F0 ) THEN
C                          CAS DE 2 FACES DU CF IDENTIQUES DANS L'ENCOCHE
                       print*,'Pb tetreto: 2FACES avec 3 MEMES SOMMETS'
     %                        ,' NF0=',NF0,' NF1=',NF1,
     %                         ' ST:',NS1F0,NS2F0,NS3F0
                          GOTO 159
                        ENDIF

C                       NS1F0-NS2F0-NS3F0-NS3F1 TETRAEDRE DE VOLUME>0?
C                       SON VOLUME ET SA QUALITE
                        CALL QUATETD( PTXYZD(1,NS1F0),
     %                                PTXYZD(1,NS2F0),
     %                                PTXYZD(1,NS3F0),
     %                                PTXYZD(1,NS3F1),
     %                                ARMIN, ARMAX, SURFTR, V, Q )

                        IF( Q .LE. QUVEVM ) THEN
C                          TETRAEDRE DE VOLUME<=0 OU PROCHE 0
C                          ou QUALITE TROP FAIBLE
C                          => RECHERCHE D'UNE AUTRE FACE NF1
                           GOTO 159
                        ENDIF

C                       L'ARETE NS3F0-NS3F1 INTERSECTE T ELLE UNE FACE
C                       DE LA SOUS-ETOILE?
                        NF2 = N1FEOC
 153                    IF( NF2 .GT. 0 ) THEN
                           IF( NF2 .NE. NF0 .AND. NF2 .NE. NF1 ) THEN

C                             CETTE FACE NF2 EST ELLE INTERSECTEE PAR 
C                             L'ARETE NS3F0-NS3F1?
                              CALL INARTRNS( PTXYZD, NS3F0, NS3F1,
     %                        NFETOI(2,NF2),NFETOI(3,NF2),NFETOI(4,NF2),
     %                                      LINTER, XYZ, CBTR )
C                             LINTER: -2 SI S1=S2
C                             -1 SI S1-S2 PARALLELE AU PLAN DU TRIANGLE
C                              0 SI S1-S2 N'INTERSECTE PAS LE TRIANGLE
C                              1 SI S1-S2   INTERSECTE     LE TRIANGLE ET ENTRE S1-S2
C                            XYZ:  3 COORDONNEES DU POINT D'INTERSECTION SI LINTER=1
C                            CBTR: 3 COORDONNEES BARYCENTRIQUES DE XYZ DANS LE TRIANGLE
C                              2 LE SOMMET NS3F0 EST UN SOMMET DU TRIANGLE ET
C                                LE SOMMET NS3F1 VOIT LE TRIANGLE PAR DERRIERE (VOLUME<0)
C                                ou
C                                LE SOMMET NS3F1 EST UN SOMMET DU TRIANGLE ET
C                                LE SOMMET NS3F0 VOIT LE TRIANGLE PAR DERRIERE (VOLUME<0)

                              IF( LINTER .EQ. 2 ) THEN
C                                 LA FACE EST VUE PAR DERRIERE PAR
C                                 UN DES 2 POINTS DE L'ARETE
C                                 ABANDON DE CE CHOIX DE 2 FACES
                                  GOTO 159
                              ENDIF

                              IF( LINTER  .EQ. 1 .AND.
     %                            CBTR(1) .LT. 0.999999D0 .AND.
     %                            CBTR(2) .LT. 0.999999D0 .AND.
     %                            CBTR(3) .LT. 0.999999D0 ) THEN
C                                 OUI: UN POINT D'INTERSECTION
C                                 ABANDON DE CE CHOIX DE 2 FACES
                                  GOTO 159
                              ENDIF

C                             UNE DES 3 ARETES DE NF2 INTERSECTE T ELLE
C                             LA FACE NS3F0 NS2F0 NS3F1?
                              N0 = NFETOI(4,NF2)
                              DO M=1,3
                                 N1 = NFETOI(1+M,NF2)
                                 CALL INARTRNS( PTXYZD, N0,    N1,
     %                                          NS3F0,  NS2F0, NS3F1,
     %                                          LINTER, XYZ,   CBTR )

                                 IF( LINTER .EQ. 2 ) THEN
C                                    LA FACE EST VUE PAR DERRIERE PAR
C                                    UN DES 2 POINTS DE L'ARETE
C                                    ABANDON DE CE CHOIX DE 2 FACES
                                    GOTO 159
                                 ENDIF

                                 IF( LINTER  .EQ. 1 .AND.
     %                               CBTR(1) .LT. 0.999999D0 .AND.
     %                               CBTR(2) .LT. 0.999999D0 .AND.
     %                               CBTR(3) .LT. 0.999999D0 ) THEN
C                                    OUI: UN POINT D'INTERSECTION
C                                    ABANDON DE CE CHOIX DE 2 FACES
                                     GOTO 159
                                 ENDIF
                                 N0 = N1
                              ENDDO

C                             UNE DES 3 ARETES DE NF2 INTERSECTE T ELLE
C                             LA FACE NS1F0 NS3F0 NS3F1?
                              N0 = NFETOI(4,NF2)
                              DO M=1,3
                                 N1 = NFETOI(1+M,NF2)
                                 CALL INARTRNS( PTXYZD, N0,    N1,
     %                                          NS1F0,  NS3F0, NS3F1,
     %                                          LINTER, XYZ,   CBTR )

                                 IF( LINTER .EQ. 2 ) THEN
C                                    LA FACE EST VUE PAR DERRIERE
C                                    PAR UN DES 2 POINTS DE L'ARETE
C                                    ABANDON DE CE CHOIX DE 2 FACES
                                    GOTO 159
                                 ENDIF

                                 IF(  LINTER  .EQ. 1 .AND.
     %                                CBTR(1) .LT. 0.999999D0 .AND.
     %                                CBTR(2) .LT. 0.999999D0 .AND.
     %                                CBTR(3) .LT. 0.999999D0 ) THEN
C                                     OUI: UN POINT D'INTERSECTION
C                                     ABANDON DE CE CHOIX DE 2 FACES
                                    GOTO 159
                                 ENDIF
                                 N0 = N1
                              ENDDO

                           ENDIF

C                          LA FACE SUIVANTE DE NF2
                           NF2 = NFETOI(5,NF2)
                           GOTO 153

                        ENDIF

C                       CREATION DU TETRAEDRE NS1F0-NS2F0-NS3F0-NS3F1
C                       AVEC LES 2 FACES NF0-NF1 DE LA SOUS-ETOILE
C                       ---------------------------------------------
                        NOCAS = 4

C                       TETRAEDRE INCONNU OPPOSE AUX 2 AUTRES FACES DE NF0-NF1
                        NTEOPF2 = -1
                        NTEOPF3 = -1
                        NTEOPF4 = NFETOI(1,NF1)

C                       CREATION DU TETRAEDRE NS1F0-NS2F0-NS3F0-NS3F1
                        NFETOI(2,NF0) = NS1F0
                        NFETOI(3,NF0) = NS2F0
                        NFETOI(4,NF0) = NS3F0
                        NTEOPF1 = NFETOI(1,NF0)

                        CALL AJTEET( NS1F0,   NS2F0,   NS3F0,   NS3F1,
     %                               NTEOPF1, NTEOPF2, NTEOPF3, NTEOPF4,
     %                               NBSOMM,  PTXYZD,  QUVEVM,
     %                          MXTETR, N1TEVI, NOTETR, NUDTETR,
     %                          N1TETS, MXTECR, NBTECR, NOTECR, VOLTECR,
     %                          VOLTEC, QUATET, IERR )
                        IF( IERR .NE. 0 ) GOTO 9999

                        NTE = NOTECR( NBTECR )
                        print *,'tetreto: NOCAS=',NOCAS,
     %                  ': 2Faces=> +TETRAEDRE',NTE,'=',
     %                  (NOTETR(KK,NTE),KK=1,8),' N1TEVI=',N1TEVI,
     %                  ' V=',VOLTEC,' Q=',QUATET

C                       AJOUTER EVENTUELLEMENT LES 4 FACES DU
C                       TETRAEDRE NTE DANS LE TABLEAU LEFACO
                        CALL AJTELEFA( NTE,    NOTETR,
     %                                 INFACO, MXFACO, LEFACO )

C                       EXAMEN DE LA SOUS-ETOILE APRES AJOUT DU TETRAEDRE NTE
                        GOTO 800

                     ENDIF
                  ENDDO
               ENDIF

C              LA FACE SUIVANTE DE NF1
 159           NF1 = NFETOI(5,NF1)
               GOTO 152
            ENDIF

         ENDDO

C        LA FACE SUIVANTE DE NFETOI
         NF0 = NFETOI(5,NF0)
         GOTO 150

      ENDIF

      IF( N1FEOC .LE. 0 ) THEN
C        TOUTE LA SOUS-ETOILE A ETE TETRAEDRISEE
         GOTO 800
      ENDIF


C     -----------------------------------------------------------------------
C     4) RECHERCHE D'UNE FACE FORMANT AVEC UN SOMMET DE LA SOUS-ETOILE PROCHE
C        UN TETRAEDRE DE VOLUME POSITIF ET SANS INTERSECTION AVEC LES FACES
C     -----------------------------------------------------------------------
      NS4MAX = 0
      NF0MAX = 0
      VOLMAX = -1D123
      QUAMAX = 0.0
      NF0 = N1FEOC
 200  IF( NF0 .GT. 0 ) THEN

C        RECHERCHE D'UN SOMMET NS4 D'UNE FACE NF1 NON ADJACENTE A NF0
C        ET OPPOSEE DANS LA SOUS-ETOILE
C        ET SANS INTERSECTION AVEC LA SOUS-ETOILE
         NF1 = N1FEOC
 210     IF( NF1 .GT. 0 ) THEN
            IF( NF1 .NE. NF0 ) THEN

               DO 250 L=1,3

C                 LE SOMMET L DE LA FACE NF1
                  NS4 = NFETOI(1+L,NF1)
                  IF( NS4 .EQ. NFETOI(2,NF0)  .OR.
     %                NS4 .EQ. NFETOI(3,NF0)  .OR.
     %                NS4 .EQ. NFETOI(4,NF0) ) GOTO 250

C                 NS4 EST SOMMET DE NF1 MAIS PAS SOMMET DE NF0
C                 LE VOLUME DU TETRAEDRE FACE NF0 + NS4 EST-IL POSITIF?
                  CALL QUATETD( PTXYZD(1,NFETOI(2,NF0) ),
     %                          PTXYZD(1,NFETOI(3,NF0) ),
     %                          PTXYZD(1,NFETOI(4,NF0) ),
     %                          PTXYZD(1,NS4),
     %                          ARMIN, ARMAX, SURFTR, V, Q )

                  IF( Q .LE. QUVEVM ) GOTO 250

C                 LE TETRAEDRE NF0-NS4 EST IL SANS INTERSECTION AVEC LA SOUS-ETOILE?
                  NF2 = N1FEOC
 230              IF( NF2 .GT. 0 ) THEN
                     IF( NF2 .NE. NF0 ) THEN

C                       INTERSECTION TRIANGLE NF2-TETRAEDRE NF0+NS4?
                        CALL INTRITET( PTXYZD(1,NFETOI(2,NF2)),
     %                                 PTXYZD(1,NFETOI(3,NF2)),
     %                                 PTXYZD(1,NFETOI(4,NF2)),
     %                                 PTXYZD(1,NFETOI(2,NF0)),
     %                                 PTXYZD(1,NFETOI(3,NF0)),
     %                                 PTXYZD(1,NFETOI(4,NF0)),
     %                                 PTXYZD(1,NS4),
     %                                 LINTER )
                        IF( LINTER  .NE. 0 ) THEN
C                          OUI: UN POINT D'INTERSECTION=>ABANDON DE CE SOMMET NS4
                           GOTO 250
                        ENDIF

                     ENDIF

C                    LA FACE SUIVANTE DE NF2
                     NF2 = NFETOI(5,NF2)
                     GOTO 230

                  ENDIF

C                 TETRAEDRE FACE NF0 + NS4  de V>0 ET SANS INTERSECTION
ccc                  IF( V .GT. VOLMAX ) THEN  9/9/17
                  IF( Q .GT. QUAMAX ) THEN
C                    TETRAEDRE DE MEILLEURE QUALITE
                     VOLMAX = V
                     QUAMAX = Q
                     NF0MAX = NF0
                     NS4MAX = NS4
                     CALL QUATETD( PTXYZD(1,NFETOI(2,NF0)),
     %                             PTXYZD(1,NFETOI(3,NF0)),
     %                             PTXYZD(1,NFETOI(4,NF0)),
     %                             PTXYZD(1,NS4MAX),
     %                             ARMIN, ARMAX, SURFTR, V, Q )
                     NOCAS = 5
ccc                    print*,'tetreto: NOCAS=',NOCAS,': 1F+1StPr NS4MAX='
ccc     %                    ,NS4MAX,' +NF0=',NF0,' V=',V,' Q=',Q,
ccc     %                     ' EST POSSIBLE'

                  ENDIF
 250           ENDDO
            ENDIF

C           LA FACE SUIVANTE DE NF1
            NF1 = NFETOI(5,NF1)
            GOTO 210

         ENDIF

C        LA FACE SUIVANTE DE NFETOI
         NF0 = NFETOI(5,NF0)
         GOTO 200

      ENDIF

      IF( VOLMAX .EQ. -1D123 ) THEN
C        AUCUN SOMMET NS4 CONVIENT
         GOTO 300
      ENDIF

C     LE TETRAEDRE NF0MAX-NS4MAX EST CREE A TRAVERS NFVPSI
C     S'IL N'A PAS ETE LE POINT CHOISI PRECEDANT
      IF( NS4MAX .EQ. NS4MAX0 .OR. NS4MAX .EQ. NS4MAX00 ) THEN
C        ABANDON POUR EVITER DES CALCULS INUTILES
         PRINT*,'tetreto: RUPTURE de la BOUCLE sur 2 SOMMETS NS4MAX00='
     %         ,NS4MAX00,' NS4MAX0=',NS4MAX0,' et NS4MAX=',NS4MAX
         GOTO 300
      ENDIF

C     LE TETRAEDRE NF0MAX-NS4MAX EST CREE A TRAVERS NFVPSI
C     ----------------------------------------------------
      NBVPSI    = 1
      NFVPSI(1) = NF0MAX
      NS4MAX00  = NS4MAX0
      NS4MAX0   = NS4MAX
      NOCAS     = 5
      print*,'tetreto: NOCAS=',NOCAS,': 1F+1St NS4MAX=',NS4MAX,
     %' +NF0=',NF0MAX,' de St',(NFETOI(M,NF0MAX),M=2,4),
     %  ' V=',VOLMAX,' Q=',QUAMAX,' EST CREE avec QUVEVM=',QUVEVM
      GOTO 510


C     -----------------------------------------------------------------
C     5) RECHERCHE DE 2 ARETES DOUBLES DES FACES SIMPLES DE L'ETOILE
C        FORMANT LE TETRAEDRE DE LA MEILLEURE QUALITE POSSIBLE ET
C        SANS INTERSECTION AVEC LES FACES SIMPLES DE L'ETOILE ET
C        AUCUNE DE SES 4 FACES NE DOIT ETRE UNE FACE SIMPLE DE L'ETOILE
C     -----------------------------------------------------------------
 300  IF( NBFETO .GE. 0 ) GOTO 400
CCC       DESSOUS INACTIF ... test sans le cas 2 aretes


      QMAX = -2.
      VMAX = 0D0
      NF0  = N1FEOC

 310  IF( NF0 .GT. 0 ) THEN

C        PARCOURS DES ARETES DOUBLES DES FACES SIMPLES DE L'ETOILE
         DO 350 M1 = 1, 3

            NS1 = NFETOI( 1+M1, NF0 )
            IF( M1 .EQ. 3 ) THEN
               M2 = 1
            ELSE
               M2 = M1 + 1
            ENDIF
            NS2 = NFETOI( 1+M2, NF0 )

C           PARCOURS DES ARETES DOUBLES DES FACES SIMPLES DE L'ETOILE
C           AU DELA DE LA FACE NF0
            NF1 = NFETOI( 5, NF0 )
 320        IF( NF1 .GT. 0 ) THEN
               DO 340 N1 = 1, 3

C                 ARETE DANS L'ORDRE NS4-NS3 POUR OBTENIR UN VOLUME
C                 AVEC L'INTERIEUR DE L'ETOILE
                  NS4 = NFETOI( 1+N1, NF1 )
                  IF( NS4 .EQ. NS1 .OR. NS4 .EQ. NS2 ) GOTO 340

                  IF( N1 .EQ. 3 ) THEN
                     N2 = 1
                  ELSE
                     N2 = N1 + 1
                  ENDIF
                  NS3 = NFETOI( 1+N2, NF1 )

                  IF( NS3 .EQ. NS1 .OR. NS3 .EQ. NS2 ) GOTO 340

C                 LE TETRAEDRE DE SOMMETS PTXYZD NS1 NS2 NS3 NS4
C                 SOIT DE VOLUME>0 ET SA QUALITE>=0
                  CALL QUATETD( PTXYZD(1,NS1), PTXYZD(1,NS2),
     %                          PTXYZD(1,NS3), PTXYZD(1,NS4),
     %                          ARMIN, ARMAX, SURFTR, V, Q )
                  IF( V .LT. 0D0 ) THEN
C                    TETRAEDRE DE VOLUME<0 : NE CONVIENT PAS
C                    NON INTERIEUR A L'ETOILE
                     GOTO 340
                  ENDIF

                  IF( Q .LE. QMAX ) GOTO 340

C                 LA QUALITE DU TETRAEDRE EST MEILLEURE QUE QMAX
C                 LE TETRAEDRE INTERSECTE T IL UNE DES FACES DE L'ETOILE NFETOI 
C                 DE PREMIERE FACE N1FEOC?
                  CALL INTTESET( PTXYZD, NS1, NS2, NS3, NS4,
     %                           N1FEOC, NFETOI, LINTER )
C                 LINTER : 0 PAS D'INTERSECTION
C                 1 SI UNE ARETE D'UNE FACE NFETOI INTERSECTE UNE FACE DU TETRAEDRE
C                 2 SI UNE FACE NFETOI EST INTERSECTEE PAR UNE ARETE DU TETRAEDRE
C                 3 SI UN DES SOMMETS D'UNE FACE NFETOI EST STRICTEMENT INTERNE
C                   AU TETRAEDRE
                  IF( LINTER .EQ. 0 ) THEN

C                    NON: UNE FACE SIMPLE DE L'ETOILE EST ELLE
C                         UNE DES FACES DU TETRAEDRE NS1 NS2 NS3 NS4=NOSOTE?
                     NF3 = N1FEOC
 330                 IF( NF3 .GT. 0 ) THEN

                        CALL TRI3NO( NFETOI(2,NF3), NOSOTR2 )
C                       NS1 NS2 NS3 NS4 EN EQUIVALENCE AVEC NOSOTE
                        CALL NO1F1T( NOSOTR2, NOSOTE, NFL )
                        IF( NFL .GT. 0 ) THEN
C                          LE TETRAEDRE A UNE FACE SIMPLE DE L'ETOILE => REJET
                           GOTO 340
                        ENDIF

C                       UNE DES ARETES 2 3 4 5 DU TETRAEDRE NOSOTE
C                       EST ELLE UNE ARETE DOUBLE DES FACES SIMPLES DE L'ETOILE?
                        DO 336 NA=2,5
                           NSA1 = NOSOTE( NOSOARTE(1,NA) )
C                          NSA1 EST IL UN SOMMET DE NF3?
                           DO I=2,4
                              NSF = NFETOI( I, NF3 )
                              IF( NSF .EQ. NSA1 ) GOTO 333
                           ENDDO
                           GOTO 336

C                          NSA1 EST UN SOMMET DE NF3
 333                       NSA2 = NOSOTE( NOSOARTE(2,NA) )
                           DO I=2,4
                              NSF = NFETOI( I, NF3 )
                              IF( NSF .EQ. NSA2 ) THEN
C                                OUI: NSA1 NSA2 EST UNE ARETE DE NF3 => REJET
                                 GOTO 340
                              ENDIF
                           ENDDO
 336                    ENDDO

C                       NON: LES 4 ARETES DE NOSOTE NE SONT PAS DES ARETES
C                            DE L'ETOILE

C                       PASSAGE A LA FACE SIMPLE SUIVANTE DE L'ETOILE
                        NF3 = NFETOI( 5, NF3 )
                        GOTO 330

                     ENDIF

C                    ICI: LE TETRAEDRE NS1 NS2 NS3 NS4 DE MEILLEURE QUALITE
C                    INTERSECTE AUCUNE FACE SIMPLE DE L'ETOILE
C                    AUCUNE DE SES FACES EST UNE FACE SIMPLE DE L'ETOILE
C                    AUCUNE DE SES 4 AUTRES ARETES N'EST UNE ARETE DE L'ETOILE
C                    => UN MEILLEUR CHOIX POUR L'AJOUTER AUX TETRAEDRES
C                       DE L'ETOILE
C                    ----------------------------------------------------------
                     QMAX   = Q
                     VMAX   = V
                     NS1MAX = NS1
                     NS2MAX = NS2
                     NS3MAX = NS3
                     NS4MAX = NS4

                  ENDIF

 340           ENDDO

C              PASSAGE A LA FACE SIMPLE SUIVANTE DE L'ETOILE
               NF1 = NFETOI( 5, NF1 )
               GOTO 320

            ENDIF

 350     ENDDO

C        PASSAGE A LA FACE SIMPLE SUIVANTE DE L'ETOILE
         NF0 = NFETOI( 5, NF0 )
         GOTO 310

      ENDIF

      IF( QMAX .GT. 0.0 .AND. VMAX .GT. 0D0 ) THEN

C        IL EXISTE UN TETRAEDRE DE SOMMETS NS1MAX NS2MAX NS3MAX NS4MAX
C        EXTREMITES DE 2 ARETES DOUBLES DES FACES SIMPLES DE LA SOUS-ETOILE
C        DE VOLUME>0 et de QUALITE QMAX SANS INTERSECTION AVEC LES FACES
C        SIMPLES DE L'ETOILE -> IL EST CREE et AJOUTE A LA SOUS-ETOILE
C        ------------------------------------------------------------------
         NOCAS = 6

C        LES 4 FACES ONT DES TETRAEDRES OPPOSES INCONNUS
         NTEOPF1 = -1
         NTEOPF2 = -1
         NTEOPF3 = -1
         NTEOPF4 = -1

C        CREATION DU TETRAEDRE NS1MAX NS2MAX NS3MAX NS4MAX
         CALL AJTEET( NS1MAX,  NS2MAX,  NS3MAX,  NS4MAX,
     %                NTEOPF1, NTEOPF2, NTEOPF3, NTEOPF4,
     %                NBSOMM,  PTXYZD,  QUVEVM,
     %                MXTETR,  N1TEVI,  NOTETR, NUDTETR,
     %                N1TETS,  MXTECR,  NBTECR, NOTECR, VOLTECR,
     %                VOLTEC,  QUATET,  IERR )
         IF( IERR .NE. 0 ) GOTO 9999

         NTE = NOTECR( NBTECR )
         print *,'tetreto: NOCAS=',NOCAS,
     %           ': 2Aretes=> 1TETRAEDRE',NTE,':',
     %           (NOTETR(KK,NTE),KK=1,8),' N1TEVI=',N1TEVI,
     %           ' V=',VOLTEC,' Q=',QUATET

C        AJOUTER EVENTUELLEMENT LES 4 FACES DU
C        TETRAEDRE NTE DANS LE TABLEAU LEFACO
         CALL AJTELEFA( NTE, NOTETR, INFACO, MXFACO, LEFACO )

C        TRACE DE LA SOUS-ETOILE N1FEOC A TETRAEDRISER
         KTITRE='tetreto: Ajout 1 TETRAEDRE avec 2 ARETES DES FACES SIMP
     %LES'
         CALL SANSDBL( KTITRE, L )
         PRINT*,KTITRE(1:L),' NoSSET=',NBSSET,' NBARPB=',NBARPB
         CALL TRFETO4( KTITRE(1:L), PTXYZD, N1FEOC, NFETOI,
     %                 NBTRCF, NOTRCF, LEFACO, NO0FAR,
     %                 NBTECR, NOTECR, NOTETR )
         CALL TRFETO7( PTXYZD, N1FEOC, NFETOI, NBSSET, N1SSET,
     %                 0, NSARP2F )

C        EXAMEN DE LA SOUS-ETOILE APRES AJOUT DU TETRAEDRE NTE
         GOTO 800

      ENDIF


C     --------------------------------------------------------------
C     6) LA CONFIGURATION DES FACES DE LA SOUS-ETOILE NE PERMET PAS 
C        AVEC 3 OU 2 OU 1 FACES OU 2 ARETES DE FORMER UN TETRAEDRE
C     => UN  NOUVEAU POINT EST A LOCALISER PUIS AJOUTER ET
C        TETRAEDRISER AVEC LES FACES DE LA SOUS-ETOILE

C     CHOISIR UN POINT XYZ A AJOUTER ET SIMULER LA TETRAEDRISATION
C     DE TOUT OU UNE PARTIE DE LA SOUS-ETOILE => NBVPSI FACES NFVPSI
C     DE LA SOUS-ETOILE A TETRAEDRISER AVEC LE POINT XYZ
C     --------------------------------------------------------------
 400  NOPASS6 = NOPASS6 + 1

      IF( NOPASS6 .GE. 8 ) THEN
ccc      IF( NOPASS6 .GE. 0 ) THEN
C        ABANDON POUR EVITER DES CALCULS INUTILES
         IERR = 7
ccc         TRACTE = .TRUE.
         GOTO 9999
      ENDIF

C     SIMULER LA TETRAEDRISATION DES FACES SIMPLES DE L'ETOILE
C     NFETOI ISSUES DE N1FEOC A PARTIR D'UN POINT XYZ CALCULE
C     SUR UNE GRILLE 3D CENTRE AU MILIEU DE LA PLUS GRANDE ARETE
C     DES FACES DE L'ETOILE POUR MAXIMISER LE NOMBRE DE TETRAEDRES
C     FACE+XYZ DE VOLUMES POSITIFs ET SANS INTERSECTION AVEC LES FACES
C     ----------------------------------------------------------------
      CALL SIPTGRAR( KTITRE, PTXYZD, N1FEOC, NFETOI,
     %               NBFETO, XYZ,    MXVPSI, NBVPSI, NFVPSI, QUAMIMX )

      IF( NBVPSI .EQ. NBFETO ) THEN
         NOCHOIXYZ = 1
C        TETRAEDRISATION DE TOUTES LES FACES DE LA SOUS-ETOILE
C        ETOILEES PAR XYZ
         PRINT*,'tetreto: NBFETO=NBVPSI=',NBVPSI,
     %   ' TOUTES les FACES de l''ETOILE sont ETOILEES PAR XYZ',XYZ
         GOTO 500
      ENDIF

C     MAXFAC: NOMBRE MINIMUM DE FACES A VOLUME POSITIF AU DELA DUQUEL
C     ON SE CONTENTE DU MAX TROUVE POUR JOINDRE LE POINT XYZ AUX FACES
ccc      MAXFAC = MAX( NBFETO/2+1, 3 )
ccc      MAXFAC = NBFETO
ccc      MAXFAC = MAX( NBFETO/3+1, 3 )  2/3/2018
ccc      MAXFAC = MAX( 3, MIN( NBFETO/3+1, 7 ) )
ccc      MAXFAC = MAX( NBFETO/4, 3 )
ccc      MAXFAC = 2
      MAXFAC = MAX( 3, NBFETO/4+1 )

      IF( NBVPSI .GE. MAXFAC ) GOTO 500


cccC     SIMULER LA TETRAEDRISATION DES FACES SIMPLES DE L'ETOILE
cccC     NFETOI ISSUES DE N1FEOC A PARTIR D'UN POINT XYZ CALCULES POUR
cccC     MAXIMISER LE NOMBRE DE TETRAEDRES FACE+XYZ DE VOLUMES POSITIFS
cccC     ET SANS INTERSECTION AVEC LES FACES
cccC     --------------------------------------------------------------
ccc      CALL SIPTETO2( NOPASS6, KTITRE, PTXYZD, NPSOFR,
ccc     %               N1FEOC,  NFETOI, NBFETO, MAXFAC, XYZ,
ccc     %               MXVPSI,  NBVPSI,  NFVPSI, NFVPSI0,
ccc     %               MXCIAS,  NBCIAS,  N1CIAS,
ccc     %               MXASFVP, NBASFVP, NSASFVP,
ccc     %               QTEAME,  QUAMIN,  MIARSICI,
ccc     %               NOCHOIXYZ, VETXYZ, RAVEVM )
cccC     MIARSICI : NOMBRE MIN D'ARETES SIMPLES DES CYCLES DES FACES NFVPSI
cccC     NOCHOIXYZ: NUMERO DU CAS DONNANT LE NOMBRE MAXIMAL DE FACES A VOLUME>0
cccC                SANS INTERSECTION ET OBTENUES A TRAVERS LES ARETES ADJACENTES
ccc      IF( NBVPSI .GT. 0 ) THEN
cccC        TRACE DES NBVPSI FACES SIMPLES+XYZ A VOLUME>0 ET 0 INTERSECTION
ccc      KTITRE='Fin sipteto2:        FACES+XYZ TETRAEDRES V>0 pour        
ccc     %FACES NOCHOIXYZ=   QMIN=          VETXYZ=            '
ccc         WRITE(KTITRE(14:18), '(I5)'   ) NBVPSI
ccc         WRITE(KTITRE(51:55), '(I5)'   ) NBFETO
ccc         WRITE(KTITRE(74:75), '(I2)'   ) NOCHOIXYZ
ccc         WRITE(KTITRE(82:89), '(F8.5)' ) QUAMIN
ccc         WRITE(KTITRE(99:111),'(G13.5)') VETXYZ
ccc         CALL SANSDBL( KTITRE, L )
ccc         PRINT*, KTITRE(1:L)
ccc         CALL TRFETO3( KTITRE(1:L), XYZ, XYZ, PTXYZD, N1FEOC, NFETOI,
ccc     %                 NBVPSI, NFVPSI)
ccc      ENDIF


      IF( NBVPSI .LE. 0 ) THEN
C        ==============================================================
C        ALGORITHME DE SIPTETO2 INSUFFISANT ...
C        LA SOUS-ETOILE EST ECRASEE, PLATE ...
C        ABANDON POUR EVITER DES CALCULS INUTILES
C        ==============================================================
         IERR = 7
ccc         TRACTE = .TRUE.
C        TRACE DES NBVPSI FACES SIMPLES+XYZ A VOLUME>0 ET 0 INTERSECTION
         KTITRE='PB  tetreto:        FACES+XYZ TETRAEDRES V>0 pour     
     %   FACES NOCHOIXYZ=   QMIN=          VETXYZ=            '
         WRITE(KTITRE(14:18), '(I5)'   ) NBVPSI
         WRITE(KTITRE(51:55), '(I5)'   ) NBFETO
         WRITE(KTITRE(74:75), '(I2)'   ) NOCHOIXYZ
         WRITE(KTITRE(82:89), '(F8.5)' ) QUAMIN
         WRITE(KTITRE(99:111),'(G13.5)') VETXYZ
         CALL SANSDBL( KTITRE, L )
         PRINT*, KTITRE(1:L)
         CALL TRFETO3( KTITRE(1:L), XYZ, XYZ, PTXYZD, N1FEOC, NFETOI,
     %                 NBVPSI, NFVPSI)
         GOTO 9999
      ENDIF


C     EXISTE T IL DES CYCLES d'ARETES SIMPLES, BORD DES NBVPSI>0 FACES
C     TRIANGULABLES POUR REDUIRE LE NOMBRE DE FACES RENTRANTES
C     GENEREES PAR LES TETRAEDRES FACE VPSI + XYZ POINT
C     ET TRIANGULATION DU CYCLE POUR FERMER LA TETRAEDRISATION DES
C     FACES NFVPSI
C     ----------------------------------------------------------------
      NBVPSI0 = NBVPSI
      CALL TETRCYCL( QUVEVM, KTITRE, NOFOTI,
     %               MXSOMM, NBSOMM, PTXYZD, NPSOFR,
     %               MXFETO, N1FEVI, N1FEOC, NFETOI,
     %               XYZ,    MXVPSI, NBVPSI, NFVPSI0, NFVPSI,
     %               MXCIAS, N1CIAS, MXASFVP, NSASFVP,
     %               IERR )
      IF( IERR .NE. 0 ) GOTO 9999

      IF( NBVPSI0 .LT. NBVPSI ) GOTO 510


C     =========================================================
C     TETRAEDRISATION DES NBVPSI FACES NFVPSI DE LA SOUS-ETOILE
C     ISSUES DE N1FEOC A PARTIR DU POINT XYZ COMME 4-EME SOMMET
C     =========================================================
 500  NOCAS = 7 + NOCHOIXYZ
      IF( NBSOMM .GE. MXSOMM ) THEN
         PRINT*,'tetreto: MAXIMUM de SOMMETS MXSOMM=',MXSOMM,
     %          ' A AUGMENTER'
         IERR = 2
         GOTO 9999
      ENDIF

C     CREATION DU NOUVEAU SOMMET XYZ DE LA TETRAEDRISATION
C     ----------------------------------------------------
C     XYZ COORDONNEES DU 4-EME SOMMET DES TETRAEDRES A FORMER
      CALL TAILIDEA( NOFOTI, XYZ, NCODEV, XYZ(4) )

      NBSOMM = NBSOMM + 1
      DO L=1,4
         PTXYZD(L,NBSOMM) = XYZ(L)
      ENDDO

      PRINT*,'tetreto: NOCAS=',NOCAS,': Pt AJOUTE NBSOMM=',NBSOMM,
     %       ' XYZD=',XYZ,' avec les', NBVPSI,' FACES+Pt=Volume>0'

C     POINT AJOUTE NON SUR LA SURFACE DE L'OBJET
      NPSOFR( NBSOMM ) = 0

C     NBSOMM EST LE SOMMET CENTRAL DES NBVPSI TETRAEDRES
      NS4MAX = NBSOMM


C     CONSTRUCTION DES TETRAEDRES NBSOMM+FACES NFVPSI DE LA SOUS-ETOILE
C     -----------------------------------------------------------------
 510  DO 550 K = 1, NBVPSI

C        LA FACE A TETRAEDRISER
         NF0 = ABS( NFVPSI( K ) )
         IF( NF0 .LE. 0 ) GOTO 550

C        LA FACE NF0 + POINT NS4MAX FORMENT UN NOUVEAU TETRAEDRE
         NS1F0 = NFETOI(2,NF0)
         NS2F0 = NFETOI(3,NF0)
         NS3F0 = NFETOI(4,NF0)

cccC        PAS DE VERIFICATION SUR LE VOLUME ou la QUALITE DU TETRAEDRE
cccC        CELA A DU ETRE FAIT AUPARAVANT

C        LE VOLUME FACE NF0+XYZ DEVRAIT ETRE POSITIF?
         CALL QUATETD( PTXYZD(1,NS1F0), PTXYZD(1,NS2F0),
     %                 PTXYZD(1,NS3F0), PTXYZD(1,NS4MAX),
     %                 ARMIN, ARMAX, SURFTR, VOLTEC, Q )

         IF( Q .LE. QUVEVM ) THEN
C           TETRAEDRE DE VOLUME<=0 ou QUALITE TROP FAIBLE
            print *,'tetreto: NOCAS=',NOCAS,' Bizarre Face',NF0,' St',
     %              (NFETOI(L,NF0),L=2,4),' +XYZ',NS4MAX,
     %              ' => V=',VOLTEC,' Q=',Q
ccc            GOTO 550
         ENDIF

C        NUMERO NOTETR DES TETRAEDRES OPPOSES AUX 4 FACES
         NTEOPF1 = NFETOI(1,NF0)
         NTEOPF2 = -1
         NTEOPF3 = -1
         NTEOPF4 = -1

         CALL AJTEET( NS1F0,   NS2F0,   NS3F0,   NS4MAX,
     %                NTEOPF1, NTEOPF2, NTEOPF3, NTEOPF4,
     %                NS4MAX,  PTXYZD,  QUVEVM,
     %                MXTETR,  N1TEVI,  NOTETR, NUDTETR,
     %                N1TETS,  MXTECR,  NBTECR, NOTECR, VOLTECR,
     %                VOLTEC,  QUATET,  IERR )

         NTE = NOTECR(NBTECR)
         print *,'tetreto: NOCAS=',NOCAS,
     %           ': 1Face + XYZ4MAX => TETRAEDRE',NTE,'=',
     %            (NOTETR(KK,NTE),KK=1,8 ),' V=',VOLTEC,
     %           ' Q=',QUATET

C        AJOUTER EVENTUELLEMENT LES 4 FACES DU
C        TETRAEDRE NTE DANS LE TABLEAU LEFACO
         CALL AJTELEFA( NTE, NOTETR, INFACO, MXFACO, LEFACO )

 550  ENDDO


C     SUPPRESSION DE NFETOI DES FACES AJOUTEES NON STANDARD DE LA
C     TRIANGULATION DES CYCLES
C     -----------------------------------------------------------
      N      = NBVPSI
      NBVPSI = 0
      DO K = 1, N
         NF0 = NFVPSI(K)
         IF( NF0 .LT. 0 ) THEN
            NF0 = -NF0
            NFETOI( 5, NF0 ) = N1FEVI
            N1FEVI = NF0
         ELSE IF( NF0 .GT. 0 ) THEN
            NBVPSI = NBVPSI + 1
            NFVPSI( NBVPSI ) = NF0
         ENDIF
      ENDDO


C     MISE A JOUR DES TETRAEDRES OPPOSES DES FACES DES TETRAEDRES
C     CREES AVEC LES NBVPSI FACES NFVPSI
C     -----------------------------------------------------------
ccc      CALL MJOPTE( NBTECR, NOTECR, N1TETS, NOTETR, MXTETR,
ccc     %             N1TEVI, PTXYZD, NBFANR )
ccc      IF( NBFANR .EQ. 0 ) GOTO 610
cccC     AU DESSOUS VERSION PLUS FIABLE et RAPIDE ....?

      DO NTC = NBTECR0+1, NBTECR-1

         NTE = NOTECR( NTC )

C        LES 4 FACES DU TETRAEDRE NTE SONT RECHERCHEES
C        DANS LES AUTRES TETRAEDRES DE LA SOUS-ETOILE
         DO 600 L=1,4

C           LA FACE L DU TETRAEDRE NTE
            IF( NOTETR(4+L,NTE) .LE. 0 ) THEN

C              FACE SANS TETRAEDRE OPPOSE
C              LES 3 SOMMETS DE LA FACE L DE NTE
               DO K=1,3
                  NOSOTR3(K) = NOTETR( NOSOFATE(K,L), NTE )
               ENDDO
               CALL TRI3NO( NOSOTR3, NOSOTR3 )

C              PARCOURS DES TETRAEDRES AU DELA DE NTC
               DO NTC1 = NTC+1, NBTECR

                  NTE1 = NOTECR( NTC1 )
C                 LA FACE L DE NTE NOSOTR3 EST ELLE UNE FACE LL DE NTE1
                  CALL NUFATRTE( NOSOTR3, NOTETR(1,NTE1), LL )

                  IF( LL .GT. 0 ) THEN
C                    LA FACE L  DU TETRAEDRE NTE  EST 
C                    LA FACE LL DU TETRAEDRE NTE1
                     NOTETR(4+L ,NTE ) = NTE1
                     NOTETR(4+LL,NTE1) = NTE
                     GOTO 600
                  ENDIF

               ENDDO

            ENDIF

 600     ENDDO

      ENDDO

cccC     VERIFICATION DES TETRAEDRES OPPOSES  
ccc 610  NTC = NBTECR-NBTECR0
ccc      CALL VEOPTE( NTC, NOTECR(NBTECR0+1), NOTETR, PTXYZD, NBFANR )
ccc      IF( NBFANR .GT. 0 ) THEN
ccc         PRINT*,'tetreto: NOMBRE TETRAEDRES OPPOSES INCONNUS=',NBFANR
ccc      ENDIF


C     ================================================================
C     A PARTIR DE L'AJOUT DES NBTECR-NBTECR0 TETRAEDRES CREES
C     MISE A JOUR DES FACES NFETOI DE LA SOUS-ETOILE NBSSET+1 RESTANTE
C     MISE A JOUR DU NOMBRE ET CONTENU DES SOUS-ETOILES
C     ================================================================
C     NOMBRE DE FACES AJOUTEES A LA SOUS-ETOILE PAR L'AJOUT DES TETRAEDRES
 800  NBFAAJ = 0
C     NOMBRE DE FACES RETIREES DE LA SOUS-ETOILE
      NBFARE = 0

      DO NTC = NBTECR0+1, NBTECR

C        NUMERO NOTETR DU TETRAEDRE NTC
         NTE = NOTECR( NTC )

C        LES 4 FACES DU TETRAEDRE NTE SONT AJOUTEES A LA SOUS-ETOILE NFETOI
C         => 2 ou 3 ou 4 FACES DE NTE Y SONT EN FAIT RETIREES
C            1 ou 2 ou 3 ou 4 FACES PEUVENT LUI ETRE AJOUTEES
C        REMARQUE: UNE FACE DE LA TRIANGULATION D'UN CYCLE A ETE CONSTRUITE
C                  UTILISEE AVEC XYZ POUR CREER UN TETRAEDRE
C                  SUPPRIMEE DE NFETOI ET AU DESSOUS ELLE REAPPARAIT COMME
C                  FACE NFETOI A PARTIR DE LA FACE DU TETRAEDRE...
         DO K=1,4
C           LA FACE K DU TETRAEDRE NTE
            CALL AJFAET2( NTE, K, NOTETR, N1FEOC, N1FEVI, NFETOI,
     %                    NFA, NBFARE, NBFAAJ )
C           NFA : >0 NUMERO DE CETTE FACE AJOUTEE DANS NFETOI
C                 =0 SUPPRESSION DE LA FACE DE LA SOUS-ETOILE
         ENDDO

      ENDDO

C     TRACE DES FACES DE LA SOUS-ETOILE ISSUE DE N1FEOC
ccc      CALL TRFETO7( PTXYZD, N1FEOC, NFETOI, NBSSET, N1SSET,
      CALL TRFETO7( PTXYZD, N1FEOC, NFETOI, 0, N1SSET,
     %              0, NSARP2F )

cccC     NBFETO NOMBRE DE FACES DE LA SOUS-ETOILE ACTUELLE RESTANTE
ccc      NBFETO = 0
ccc      NF0    = N1FEOC
ccc 1010 IF( NF0 .GT. 0 ) THEN
cccC        UNE FACE DE PLUS
ccc         NBFETO = NBFETO + 1
cccC        LA FACE SUIVANTE DE NF0 DANS NFETOI
ccc         NF0  = NFETOI(5,NF0)
ccc         GOTO 1010
ccc      ENDIF
ccc      KTITRE='tetreto: LES       FACES VPSI de SOMMET NBSOMM=          
ccc     %NBTECR0=       NBTECR=       SOUS-ETOILE=       NBFETO=     '
ccc      WRITE( KTITRE(15:19), '(I5)' ) NBVPSI
ccc      WRITE( KTITRE(49:57), '(I9)' ) NBSOMM
ccc      WRITE( KTITRE(67:71), '(I5)' ) NBTECR0
ccc      WRITE( KTITRE(81:85), '(I5)' ) NBTECR
ccc      WRITE( KTITRE(100:103),'(I4)') NBSSET+1
ccc      WRITE( KTITRE(114:117),'(I4)') NBFETO
ccc      CALL SANSDBL( KTITRE, L )
ccc      PRINT*, KTITRE(1:L)
ccc      CALL TRFETO3( KTITRE(1:L), PTXYZD(1,NBSOMM), PTXYZD(1,NBSOMM),
ccc     %              PTXYZD, N1FEOC, NFETOI, NBVPSI, NFVPSI )

cccC     VERIFICATION QUE TOUTE ARETE DE L'ETOILE RESTANTE APPARTIENT
cccC     SEULEMENT A 2n FACES avec n=1,2
cccC     ------------------------------------------------------------
ccc      CALL VETAFET( NOTETR, N1FEOC, NFETOI, NBFETO, NBARPB )
ccc      IF( NBARPB .GT. 0 ) THEN
ccc         CALL AFETOI( N1FEOC, NFETOI )
ccc      ENDIF

ccc      PRINT*,'tetreto:',NBFARE,' FACES RETIREES',
ccc     %        NBFAAJ,' FACES AJOUTEES et l''ETOILE a maintenant',
ccc     %        NBFETO,' FACES'

cccC Verification visuelle si les tetraedres opposes aux faces sont corrects
cccc     a supprimer ensuite
ccc      do ntc=1,nbtecr
ccc      nte = notecr(ntc)
ccc      print *,'tetreto:',ntc,' notetr(',nte,')=',(notetr(kk,nte),kk=1,8)
ccc      enddo
ccc      print *,'tetreto:',NBFETO,' FACES DANS L''ETOILE'

cccC     LE TITRE DU TRACE DES FACES DE LA SOUS-ETOILE ACTUELLE
cccC     TRACE EN NOIR DES 6 ARETES DU DERNIER TETRAEDRE CREE
ccc      KTITRE='tetreto:       TETRAEDRES CREES et        FACES DE L''ETO
ccc     %ILE NOCAS=        '
cccC     AJOUT DANS LE TITRE DU NOMBRE DE FACES DE LA SOUS-ETOILE ET
cccC     DE TETRAEDRES CREES
ccc      WRITE(KTITRE(10:14),'(I5)') NBTECR
ccc      WRITE(KTITRE(37:41),'(I5)') NBFETO
ccc      WRITE(KTITRE(69:70),'(I2)') NOCAS
ccc      CALL SANSDBL( KTITRE, L )
ccc      CALL TRFETO4( KTITRE(1:L), PTXYZD, N1FEOC, NFETOI,
ccc     %              NBTRCF, NOTRCF, LEFACO, NO0FAR,
ccc     %              NBTECR, NOTECR, NOTETR )


      IF( N1FEOC .GT. 0 ) THEN

C        TOUTE LA SOUS-ETOILE N'A PAS ENCORE ETE TETRAEDRISEE
C        => UN NOUVEAU TRAITEMENT DE LA SOUS-ETOILE NBSSET+1
C        ----------------------------------------------------
C        LA TETRAEDRISATION A T ELLE CREEE PLUSIEURS SOUS-SOUS-ETOILES?
         NBSS0 = NBSSET
         CALL SOUSETO( PTXYZD, NOTETR, N1FEOC, NFETOI,
     %                 MXSSET, NBSSET, N1SSET, IERR )

         IF( NBSS0 .LT. NBSSET ) THEN
C           TRACE DES FACES ISSUES DES NBSSET SOUS-ETOILES
            CALL TRFETO7( PTXYZD, 0, NFETOI, NBSSET, N1SSET,
     %                    0, NSARP2F )
         ENDIF

         IF( IERR .NE. 0 ) GOTO 9999

      ENDIF

C     RETOUR A LA TETRAEDRISATION DE LA SOUS ETOILE EN HAUT DE PILE
      GOTO 5


C     MISE A JOUR DES TETRAEDRES OPPOSES DES FACES DES NBTECR TETRAEDRES
C     CREES POUR TOUTES LES SOUS-ETOILES de 1 A NBSSET0
C     ------------------------------------------------------------------
 9000 CALL MJOPTE( NBTECR, NOTECR, N1TETS, NOTETR, MXTETR,
     %             N1TEVI, PTXYZD, NBFANR )

C     VERIFIER L'ABSENCE DE TETRAEDRE  OPPOSE  NEGATIF (INCONNU)
C     VERIFIER L'ABSENCE DE TETRAEDRES OPPOSES DOUBLES
C     VERIFIER L'ABSENCE DE FACE COMMUNE A AU MOINS 3 TETRAEDRES
C     VERIFIER L'OPPOSITION DES TETRAEDRES 
C     DANS LA LISTE DES NBTECR TETRAEDRES NOTECR
C     ----------------------------------------------------------
      CALL VEOPTE( NBTECR, NOTECR, NOTETR, PTXYZD, NBFANR )


C     CALCUL DES QUALITES MIN MOYENNE et VOLUMES des NBTECR
C     TETRAEDRES NOTECR CREES
C     -----------------------------------------------------
      CALL QUALGRTE( PTXYZD, MXTETR, NOTETR, NBTECR, NOTECR,
     %               QUAMED, NBTMED,
     %               QUAMIN, QUAMOY, VOLUMT, VOLMOY )

      PRINT*,'tetreto: Fin',NBTECR,
     %       ' TETRAEDRES CREES de VOLUME',VOLUMT,
     %       ' VOLUME MOYEN=',VOLMOY,
     %       ' QUALITE MOYENNE=',QUAMOY,
     %       ' QUALITE MIN=',QUAMIN,' et',
     %        NBTMED,' TETRAEDRES de QUALITE<',QUAMED

ccc      DO K = 1, NBTECR
ccc         NTE = NOTECR(K)
ccc         PRINT*,'tetreto: TETRAEDRE(',NTE,') St:',(NOTETR(L,NTE),L=1,4)
ccc     %         ,' NteOp:',(NOTETR(L,NTE),L=5,8)
ccc      ENDDO

C     TRACE FINAL DES NBTECR TETRAEDRES CREES FORMANT L'ETOILE TOTALE
      KTITRE = 'Fin de tetreto:            TETRAEDRES CREES dans l''ETOI
     %LE TOTALE'
      WRITE(KTITRE(18:25),'(I8)') NBTECR

 9010 CALL SANSDBL( KTITRE, L )
      CALL TRFETO4( KTITRE(1:L), PTXYZD, N1FEOC, NFETOI,
     %              NBTRCF, NOTRCF, LEFACO, NO0FAR,
     %              NBTECR, NOTECR, NOTETR )

 9999 PRINT*,'tetreto: Fin avec IERR=',IERR,' NoPass6=',NOPASS6,NBTECR,
     %       ' TETRAEDRES CREES de Volume',VOLTECR,' NBSSET=',NBSSET,
     %       ' NBSOM0=',NBSOMM0,' NBSOMM=',NBSOMM,' NBFANR=',NBFANR


      IF( IERR .NE. 0 .AND. NBTECR .GT. 0 ) THEN

C        ABANDON de la TETRAEDRISATION et RETOUR AUX TETRAEDRES INITIAUX
C        ---------------------------------------------------------------
ccc         tracte = .true.
ccc         PRINT 10000
         PRINT*,'Fin tetreto: PROBLEME',NBTECR,' TETRAEDRES NE RECOUVRAN
     %T PAS l''ETOILE => ABANDONNEE. NBTRCF=',NBTRCF

      KTITRE='Probleme tetreto:            TETRAEDRES CREES dans l ETOIL
     %E ABANDONNEE. NBTRCF=      '
         WRITE(KTITRE(19:26),'(I8)') NBTECR
         WRITE(KTITRE(80:84),'(I5)') NBTRCF
         CALL SANSDBL( KTITRE, L )
         CALL TRFETO4( KTITRE(1:L), PTXYZD, N1FEOC, NFETOI,
     %                 NBTRCF, NOTRCF, LEFACO, NO0FAR,
     %                 NBTECR, NOTECR, NOTETR )
        
      KTITRE='Probleme tetreto:            TETRAEDRES INITIAUX RESTAUREN
     %T l''ETOILE INITIALE NBTRCF=      '
         WRITE(KTITRE(19:26),'(I8)') NBTECFPR
         WRITE(KTITRE(86:90),'(I5)') NBTRCF
         CALL SANSDBL( KTITRE, L )
         CALL TRFETO4( KTITRE(1:L), PTXYZD, 0, NFETOI,
     %                 NBTRCF, NOTRCF, LEFACO, NO0FAR,
     %                 NBTECFPR, NOTECFPR, NOTETR )

C        RESTAURATION DES NBTECFPR TETRAEDRES SAUVEGARDES
C        DESTRUCTION  DES NBTECR   TETRAEDRES NOTECF CREES DANS TETRETO
C        MISE A JOUR  DE N1TEVI, NOTETR, N1TETS, TETRAEDRES OPPOSES
         CALL RSTETRETO( PTXYZD, 
     %                   MXTETR,   N1TEVI,   NOTETR, NUDTETR, N1TETS,
     %                   INFACO,   MXFACO,   LEFACO,
     %                   NBTECFPR, NOTECFPR, NBTECR, NOTECR )

C        LES SOMMETS CREES PRECEDAMMENT SONT SUPPRIMES
         NBSOMM = NBSOMM0

      ENDIF

ccc      PRINT 10000
      TRACTE = TRACTE0
      RETURN
      END
