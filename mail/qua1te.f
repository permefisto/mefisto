      SUBROUTINE QUA1TE( NT,     QUTEMN, GRAND,  VOLMOYTE,
     %                   MOFACE, MXFACE, LFACES, NBFAPB,
     %                   NBSOMM, MXSOMM, XYZSOM, NPSOFR, NBDM,   NUDMEF,
     %                   NBSOTE, MXTETR, MXTETA, N1TEVI, NSTETR,
     %                   VOLUMT, QUALIT, NO1TSO, MXTESO, N1TESO, NOTESO,
     %                   MXFAET, N1FEOC, N1FEVI, NFETOI, VOETOI, QUETOI,
     %                   NOTEDS, VOLET0, QUAET0, VOLET1, QUAET1, IERR  )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    AMELIORER LA QUALITE DU TETRAEDRE NT DU TABLEAU NSTETR
C -----
C
C ENTREES:
C --------
C NT     : NUMERO NSTETR DU TETRAEDRE DE QUALITE A AMELIORER
C QUTEMN : QUALITE MINIMALE AU DESSOUS DE LAQUELLE UN TETRAEDRE DOIT ETRE
C          AMELIORE
C GRAND  : PLUS GRAND REEL STOCKABLE
C VOLMOYTE: VOLUME MOYEN D'UN TETRAEDRE DE LA TETREDRISATION
C MXSOMM : NOMBRE MAXIMAL DE SOMMETS DECLARABLES DANS XYZSOM ET NPSOFR
C XYZSOM : COORDONNEES X Y Z DES NBSOMM SOMMETS DES TETRAEDRES
C NBDM   : 0 SI 1 MATERIAU=VOLUME, SINON NOMBRE DE MATERIAUX DU VOLUME
C NBSOTE : VALEUR MAXIMALE DE DECLARATION DU PREMIER INDICE DE NSTETR(>3)
C MXTETR : NOMBRE MAXIMAL DE TETRAEDRES DECLARABLES DANS NSTETR
C MXTETA : NUMERO MAXIMUM DES TETRAEDRES ACTUELS
C MXTESO : NOMBRE MAXIMAL DE NUMERO DE TETRAEDRES DES SOMMETS
C MXFAET : NOMBRE MAXIMAL DE FACES DECLARABLES DANS NFETOI
C
C TABLEAUX AUXILIAIRES :
C ----------------------
C NFETOI(5,MXFAET)  DES ENTIERS
C VOETOI(MXFAET)    DES REELS
C QUETOI(MXFAET)    DES REELS
C NOTEDS(MXTETR)    DES ENTIERS
C
C MODIFIES :
C ----------
C NBSOMM : NOMBRE ACTUEL DE SOMMETS DE LA TETRAEDRISATION
C NPSOFR : NUMERO 0 SI SOMMET INTERNE
C                 1 SI SOMMET SUR LA FRONTIERE
C                 2 SI SOMMET SUR L'INTERFACE ENTRE 2 MATERIAUX
C                -1 SI SOMMET SUPPRIME LORS DE L'AMELIORATION
C NUDMEF : NUMERO DE MATERIAU DE CHAQUE TETRAEDRE DU MAILLAGE
C          ATTENTION: CE TABLEAU EXISTE SEULEMENT SI NBDM>0
C N1TEVI : NUMERO DANS NSTETR DE LA PREMIERE PLACE VIDE
C          CHAINAGE SUIVANT DANS NSTETR(2,*) ET DERNIER A ZERO
C N1TESO : NUMERO DANS NOTESO DE LA PREMIERE PLACE VIDE
C          CHAINAGE SUIVANT DANS NOTESO(2,*) ET DERNIER A ZERO
C NSTETR : NUMERO DES 4 SOMMETS DE CHAQUE TETRAEDRE
C N1FEOC : NUMERO NFETOI DE LA PREMIERE FACE OCCUPEE
C N1FEVI : NUMERO NFETOI DE LA PREMIERE FACE VIDE
C
C SORTIES:
C --------
C NO1TSO : NUMERO DU 1-ER TETRAEDRE DANS NOTESO DE CHAQUE SOMMET
C NOTESO : NOTESO(1,*) NUMERO DANS NSTETR DU TETRAEDRE
C          NOTESO(2,*) NUMERO DANS NOTESO DU TETRAEDRE SUIVANT
C                      0 SI C'EST LE DERNIER
C QUALIT : QUALITE DES TETRAEDRES DE LA TETRAEDRISATION
C VOLUMT : VOLUME  DES TETRAEDRES DE LA TETRAEDRISATION

C VOLET0 : VOLUME  DE L'ETOILE INITIALE
C QUAET0 : QUALITE DE L'ETOILE INITIALE
C VOLET1 : VOLUME  DE LA DERNIERE ETOILE CALCULEE
C QUAET1 : QUALITE DE LA DERNIERE ETOILE CALCULEE
C IERR   : =0 SI PAS D'ERREUR
C          >0 EN CAS DE SATURATION D'UN TABLEAU
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC   Novembre 1993
C MODIFS : ALAIN PERRONNET LJLL UPMC et St Pierre du Perray    Juin 2008
C MODIFS : ALAIN PERRONNET LJLL UPMC et St Pierre du Perray    Aout 2012
C MODIFS : ALAIN PERRONNET LJLL UPMC et St Pierre du Perray Fevrier 2016
C2345X7..............................................................012
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      COMMON / TRTETR / STOPTE,TRACTE
      LOGICAL           STOPTE,TRACTE,TRACTE0

      INTEGER           NPSOFR(MXSOMM), LFACES(NBSOTE,MXFACE),
     %                  NUDMEF(MXTETR),
     %                  NSTETR(NBSOTE,MXTETR),
     %                  NO1TSO(MXSOMM),
     %                  NOTESO(2,MXTESO) ,
     %                  NFETOI(5,MXFAET),
     %                  NOTEDS(MXTETR)
      REAL              XYZSOM(3,MXSOMM), VOLUMT(MXTETR),
     %                  QUALIT(MXTETR), VOETOI(MXFAET), QUETOI(MXFAET)

      REAL              P(3), B(3), VOLET0, VOLET1, QUAET0, QUAET1
      REAL              ARMIN, ARMAX, SURFTR(4), LONARE(6)
      INTEGER           NSOP(4), NTOP(4)
      INTEGER           NOTETRA(64)
      CHARACTER*120     KTITRE

      IF( NSTETR(1,NT) .LE. 0 ) RETURN

      TRACTE0= TRACTE
      NBCOOR = 3
      VOLET0 = 0.0
      VOLET1 = 0.0
      NBTETRA = 0

C     CALCUL DU VOLUME ET DE LA QUALITE DU TETRAEDRE NT ACTUEL
C     --------------------------------------------------------
      CALL QUATET( XYZSOM(1,NSTETR(1,NT)),
     %             XYZSOM(1,NSTETR(2,NT)),
     %             XYZSOM(1,NSTETR(3,NT)),
     %             XYZSOM(1,NSTETR(4,NT)),
     %             ARMIN, ARMAX, SURFTR, VOLUMT(NT), QUALIT(NT) )
      VNT = VOLUMT(NT)
 
C     QUAMIN0: QUALITE INITIALE DU TETRAEDRE NT
      QNT     = QUALIT( NT )
      QUAMIN0 = QUALIT( NT )

C     LA QUALITE DU TETRAEDRE NT
      QUAET0 = QUALIT( NT )
      IF( QUAET0 .GT. QUTEMN ) GOTO 9900

C     LONGUEUR DES 6 ARETES DU TETRAEDRE NT
      CALL LON6AR( XYZSOM(1,NSTETR(1,NT)),
     %             XYZSOM(1,NSTETR(2,NT)),
     %             XYZSOM(1,NSTETR(3,NT)),
     %             XYZSOM(1,NSTETR(4,NT)), LONARE )

C     LONGUEUR DE L'ARETE MOYENNE DU TETRAEDRE NT
      AREMOY = ( LONARE(1)+LONARE(2)+LONARE(3)
     %         + LONARE(4)+LONARE(5)+LONARE(6) ) / 6

C     NOMBRE NFR ET NUMERO DES TETRAEDRES OPPOSES AUX 4 FACES DE NT
      CALL NOTSOPTE( NT,  NO1TSO, NOTESO, NBSOTE, NSTETR,
     %               NFR, NTOP, NSOP )
ccc      PRINT *,'qua1te: QUALIT(',NT,')=',QUALIT(NT),' NFR=',NFR,
ccc     %        ' NTOP=',NTOP,' NSOP=',NSOP

C     CALCUL DU VOLUME MOYEN         DES TETRAEDRES OPPOSES AUTOUR DE NT
cccC            DU MINIMUM DES QUALITES DES TETRAEDRES OPPOSES AUTOUR DE NT
ccc      QUALIV = GRAND
      VOLALE = 0
      DO I=1,4
         IF( NTOP(I) .GT. 0 ) THEN
C           LE VOLUME MOYEN DES TETRAEDRES OPPOSES
            VOLALE = VOLALE + ABS( VOLUMT( NTOP(I) ) )
cccC           LA QUALITE MINIMALE DES TETRAEDRES OPPOSES
ccc            QUALIV = MIN( QUALIV, QUALIT( NTOP(I) ) )
         ENDIF
      ENDDO
      IF( NFR .GT. 0 ) VOLALE = VOLALE / NFR

C     TRACE DU TETRAEDRE ET DE SES NFR TETRAEDRES OPPOSES
C     ---------------------------------------------------
      IF( QUAMIN0 .LT. 0.001 ) THEN
         TRACTE = .TRUE.
         KTITRE = 'qua1te: TETRAEDRE           +   TETRAEDRES OPPOSES et
     % OPPOSES avant TRAITEMENT. V=                Q=         '
         WRITE(KTITRE(19:27),'(I9)' ) NT
         WRITE(KTITRE(31:31),'(I1)' ) NFR
         WRITE(KTITRE(83:97)  ,'(G14.6)' ) VNT
         WRITE(KTITRE(101:115),'(G14.6)' ) QNT
         CALL SANSDBL( KTITRE, L )
         PRINT*,KTITRE(1:L)
         NOPASS = 0
         CALL TR1TEVV( NOPASS, KTITRE(1:L), XYZSOM, NT, NBSOTE, NSTETR,
     %                 NO1TSO, NOTESO, NBTETRA, NOTETRA, VOLET0 )
         TRACTE = TRACTE0
      ENDIF

      IF( VOLUMT(NT) .LE. AREMOY**3 * 0.06 .OR.
     %    QUALIT(NT) .LT. 0.1 ) THEN

C        TRAITEMENT DU TETRAEDRE NT DE VOLUME NUL
C        ========================================
         CALL QTEVNUL( NT,     NFR,    NTOP,   NSOP,   SURFTR, VOLMOYTE,
     %                 MOFACE, MXFACE, LFACES, NBFAPB,
     %                 NBSOMM, MXSOMM, XYZSOM, GRAND,  NBDM,   NUDMEF,
     %                 NBSOTE, MXTETR, N1TEVI, NSTETR, NPSOFR,
     %                 VOLUMT, QUALIT, NO1TSO, MXTESO, N1TESO, NOTESO,
     %                 MXFAET, N1FEOC, N1FEVI, NFETOI, VOETOI, QUETOI,
     %                 NBTEDS, MXTETR, NOTEDS, MXTETA,
     %                 VOLET1, QUAET1, IERR  )
         IF( IERR .EQ. -1 ) THEN
            IERR = 0
            GOTO 9900
         ENDIF

         IF( IERR .GT. 1 ) GOTO 9900
         IERR = 0

cccC        DETECTION DES FACES APPARTENANT A 3 TETRAEDRES
ccc         CALL CRHAFAVO( NBSOTE, MXTETR, NSTETR, NUDTETR,
ccc     %                  MOFACE, MXFACE, LFACES, NBFAPB, IERR)

ccc         IF( NBFAPB .GT. 0 ) THEN
ccc            PRINT *
ccc            PRINT *,'qua1te 0: 1 NT=',NT,' avec',
ccc     %               NBFAPB,' FACES DE 3 TETRAEDRES'
ccc            TRACTE = .TRUE.
ccc         ENDIF

      ENDIF

C     BOUCLE SUR LES 6 ARETES DU TETRAEDRE NT
C     =======================================
ccc      TRACTE = .TRUE.
ccc      KTITRE = 'qua1te: TETRAEDRE           +   TETRAEDRES OPPOSES et OP
ccc     %POSES avant TRAITEMENT'
ccc      WRITE(KTITRE(19:27),'(I9)' ) NT
ccc      WRITE(KTITRE(31:31),'(I1)' ) NFR
ccc      CALL SANSDBL( KTITRE, L )
ccc      NOPASS = 0
ccc      CALL TR1TEVV( NOPASS, KTITRE(1:L), XYZSOM, NT, NBSOTE, NSTETR,
ccc     %              NO1TSO, NOTESO, NBTETRA, NOTETRA, VOLET0 )

      DO IAR=1,6

C        LES 2 SOMMETS NS4 NS5 DE L'ARETE IAR DU TETRAEDRE NT
         IF( IAR .LE. 3 ) THEN
C           NUMEROTATION SELON LE SP LON6AR DE UTIL
            NS4 = IAR
            IF( IAR .EQ. 3 ) THEN
               NS5 = 1
            ELSE
               NS5 = IAR + 1
            ENDIF
         ELSE
            NS4 = 4
            NS5 = IAR - 3
         ENDIF
         NS4 = NSTETR(NS4,NT)
         NS5 = NSTETR(NS5,NT)

C        NO DES TETRAEDRES D'ARETE NS4-NS5
         CALL NOTEAR( NS4,    NS5,
     %                NBSOTE, NSTETR, NO1TSO, NOTESO,
     %                NBTEDS, NOTEDS )

         IF( NBTEDS .EQ. 3 ) THEN
C
C           ESSAI 3T => 2T
C           --------------
            CALL QU3T2T( NT,     IAR,
     %                   XYZSOM, NBSOTE, N1TEVI, NSTETR, NBDM, NUDMEF,
     %                   NO1TSO, N1TESO, NOTESO,
     %                   NBTEDS, NOTEDS,
     %                   MXTETA, VOLUMT, QUALIT, IERR )
            IF( IERR .EQ. -1 ) THEN
C              REUSSITE
               WRITE(IMPRIM,*)'qua1te 1: NT=',NT,' 3T=>2T   ARETE ',IAR
               IERR = 0

cccC        DETECTION DES FACES APPARTENANT A 3 TETRAEDRES
ccc         CALL CRHAFAVO( NBSOTE, MXTETR, NSTETR, NUDTETR,
ccc     %                  MOFACE, MXFACE, LFACES, NBFAPB, IERR)

ccc         IF( NBFAPB .GT. 0 ) THEN
ccc            PRINT *
ccc            PRINT *,'qua1te 1: 1 NT=',NT,' avec',
ccc     %               NBFAPB,' FACES DE 3 TETRAEDRES'
ccc            TRACTE = .TRUE.
ccc         ENDIF


               GOTO 9900
            ENDIF
            IF( IERR .NE. 0 ) GOTO 9900
C
         ELSE IF( NBTEDS .GT. 3 ) THEN
C
C           ESSAI MT => 2(M-2) T
C           --------------------
            CALL QUMTNT( NT,     IAR,
     %                   XYZSOM, NBSOTE, N1TEVI, NSTETR, NBDM, NUDMEF,
     %                   NO1TSO, N1TESO, NOTESO,
     %                   NBTEDS, NOTEDS, QUETOI,
     %                   MXTETA, VOLUMT, QUALIT, IERR )
            IF( IERR .EQ. -1 ) THEN
C              REUSSITE
               WRITE(IMPRIM,*)'qua1te 2: NT=',NT,NBTEDS,'T=>',
     %                        2*(NBTEDS-2),'T   ARETE',IAR
               IERR = 0

cccC        DETECTION DES FACES APPARTENANT A 3 TETRAEDRES
ccc         CALL CRHAFAVO( NBSOTE, MXTETR, NSTETR, NUDTETR,
ccc     %                  MOFACE, MXFACE, LFACES, NBFAPB, IERR)

ccc         IF( NBFAPB .GT. 0 ) THEN
ccc            PRINT *
ccc            PRINT *,'qua1te 2: 1 NT=',NT,' avec',
ccc     %               NBFAPB,' FACES DE 3 TETRAEDRES'
ccc            TRACTE = .TRUE.
ccc         ENDIF


               GOTO 9900
            ENDIF
            IF( IERR .NE. 0 ) GOTO 9900
         ENDIF
      ENDDO


C     BOUCLE SUR LES 6 GRANDES ARETES DECROISSANTES
C     =============================================
      IF( NFR .LT. 4 ) GOTO 95
      DO IAR=1,6

C        ESSAI AVEC LA PLUS GRANDE ARETE
         D = 0
         DO J=1,6
            IF( LONARE(J) .GT. D ) THEN
               I = J
               D = LONARE(J)
            ENDIF
         ENDDO

C        L'ARETE I EST ACTUELLEMENT MAXIMALE
C        SIMULATION DE TETRAEDRISATION ETOILEE
C        PAR LE MILIEU DE CETTE ARETE
         CALL SITEMI( I,      NT,     XYZSOM,
     %                NBSOTE, NSTETR, NO1TSO, NOTESO,
     %                N1FEOC, N1FEVI, NFETOI, VOETOI, QUETOI,
     %                NBTEDS, MXTETR, NOTEDS,
     %                VOLUMT, QUALIT, VOLET0, QUAET0, VOLET1, QUAET1, P)
         IF( QUAET1 .LE. 0 ) GOTO 43

         IF( QUAET1 .GT. QUAET0 ) THEN
C
C           GENERATION EFFECTIVE DE CETTE TETRAEDRISATION ETOILEE
C           AVEC AJOUT DU MILIEU DE LA PLUS GRANDE ARETE
C           -----------------------------------------------------
            CALL GETEPT( P,      NBSOMM, MXSOMM, XYZSOM, NBDM, NUDMEF,
     %                   NBSOTE, N1TEVI, NSTETR,
     %                   NO1TSO, N1TESO, NOTESO,
     %                   N1FEOC, N1FEVI, NFETOI, NBTEDS, NOTEDS,
     %                   MXTETA, VOLUMT, QUALIT,
     %                   VOLET0, QUAET0, VOLET1, QUAET1, IERR )
            IF( IERR .NE. 0 ) GOTO 9900
            WRITE(IMPRIM,*) 'qua1te 3: NT=',NT,
     %                      ' MILIEU ARETE AJOUTE ',IAR

cccC        DETECTION DES FACES APPARTENANT A 3 TETRAEDRES
ccc         CALL CRHAFAVO( NBSOTE, MXTETR, NSTETR, NUDTETR,
ccc     %                  MOFACE, MXFACE, LFACES, NBFAPB, IERR)

ccc         IF( NBFAPB .GT. 0 ) THEN
ccc            PRINT *
ccc            PRINT *,'qua1te 3: 1 NT=',NT,' avec',
ccc     %               NBFAPB,' FACES DE 3 TETRAEDRES'
ccc            TRACTE = .TRUE.
ccc         ENDIF

            GOTO 9900
C
         ELSE
C
C           LES FACES OCCUPEES DE NFETOI DEVIENNENT VIDES
            CALL ETOIOCVI( N1FEOC, N1FEVI, NFETOI )
C
         ENDIF
C
C        POUR EVITER DE RETROUVER LA MEME ARETE
 43      LONARE(I) = -LONARE(I)

      ENDDO
C
      IF( QUALIT(NT) .GT. QUTEMN ) GOTO 9900

      IF( NFR .LT. 4 ) GOTO 95
C
C     RECHERCHE DE LA PLUS PETITE HAUTEUR DU TETRAEDRE DE NT
C     ET SA POSITION PAR RAPPORT A LA FRONTIERE
C     ======================================================
      D1   = 0
      DMIN = GRAND
      I    = 0
      DO 60 J=1,4

         D = VOLUMT(NT) / SURFTR(J)
         IF( D .LT. DMIN ) THEN
C           LE SOMMET NS DE LA HAUTEUR EST IL DEPLACABLE ?
C           LA HAUTEUR J DE SOMMET I0 EST ACTUELLEMENT MINIMALE
            IF( J .EQ. 1 ) THEN
               I0 = 4
            ELSE
               I0 = J - 1
            ENDIF
C           LE SOMMET DE LA HAUTEUR MINIMALE DU TETRAEDRE NT
            NS = NSTETR(I0,NT)
            IF( NPSOFR( NS ) .NE. 0 ) GOTO 60
C           LE SOMMET NS EST DEPLACABLE
            I    = J
            DMIN = D
         ENDIF
         IF( D .GT. D1 ) D1 = D

 60   CONTINUE
C
      IF( I .GT. 0 ) THEN
C
C        ESSAI DE RELEVER LE SOMMET DE LA PLUS PETITE HAUTEUR
C        LA HAUTEUR I DE SOMMET I0 EST MINIMALE  H = 3 * V / S
         DMIN = DMIN * 3.0
C
C        IL EXISTE UN SOMMET I0 DEPLACABLE
         IF( I .EQ. 1 ) THEN
            I0 = 4
         ELSE
            I0 = I - 1
         ENDIF
C        LE SOMMET DE LA HAUTEUR MINIMALE EST DEPLACABLE
         NS = NSTETR(I0,NT)
         IF( I .EQ. 4 ) THEN
            I1 = 1
         ELSE
            I1 = I + 1
         ENDIF
         IF( I1 .EQ. 4 ) THEN
            I2 = 1
         ELSE
            I2 = I1 + 1
         ENDIF
C        DIRECTION DE LA NORMALE B A CETTE FACE I
         CALL NORFA3( XYZSOM(1,NSTETR(I ,NT)),
     %                XYZSOM(1,NSTETR(I1,NT)),
     %                XYZSOM(1,NSTETR(I2,NT)), B, NDT )
         IF( NDT .NE. 0 ) THEN
            print *,'qua1te: norfa3 avec 2 SOMMETS IDENTIQUES'
            print *,'qua1te: nstetr(',nt,')=',(nstetr(kk,nt),kk=1,4)
            print *,'qua1te: i=',i,' i1=',i1,' i2=',i2
            GOTO 95
         ENDIF

C        ESSAI DE DEPLACER LE SOMMET NS SUR CETTE HAUTEUR
C        LA DIRECTION
         P(1) = DMIN * B(1)
         P(2) = DMIN * B(2)
         P(3) = DMIN * B(3)
         CALL QTEMXD( NS,     P,
     %                XYZSOM, NBSOTE, NSTETR, NO1TSO, NOTESO,
     %                VOLUMT, QUALIT,
     %                VOLET1, QUAET1, NTQMI1, NBTENS )

cccC        DETECTION DES FACES APPARTENANT A 3 TETRAEDRES
ccc         CALL CRHAFAVO( NBSOTE, MXTETR, NSTETR, NUDTETR,
ccc     %                  MOFACE, MXFACE, LFACES, NBFAPB, IERR)

ccc         IF( NBFAPB .GT. 0 ) THEN
ccc            PRINT *
ccc            PRINT *,'qua1te 4: 1 NT=',NT,' avec',
ccc     %               NBFAPB,' FACES DE 3 TETRAEDRES'
ccc            TRACTE = .TRUE.
ccc         ENDIF


         IF( NBTENS .LE. 0 ) GOTO 9900

      ENDIF

C     ESSAI DU BARYCENTRE COMME CENTRE DE L'ETOILE DU TETRAEDRE NT
C     ============================================================
 95   IF( QUALIT(NT) .GT. QUTEMN ) GOTO 9900
      DO J=1,3
         P(J) = ( XYZSOM(J,NSTETR(1,NT)) +
     %            XYZSOM(J,NSTETR(2,NT)) +
     %            XYZSOM(J,NSTETR(3,NT)) +
     %            XYZSOM(J,NSTETR(4,NT)) ) * 0.25
      ENDDO
      NBS = 0
      DO J=1,4
         IF( NPSOFR( NSTETR(J,NT) ) .GT. 0 ) NBS = NBS + 1
      ENDDO

C     SIMULATION DE TETRAEDRISATION ETOILEE PAR CE BARYCENTRE DE NT
      CALL SITEPT( P,      NT,     XYZSOM,
     %             NBSOTE, NSTETR, NO1TSO, NOTESO,
     %             N1FEOC, N1FEVI, NFETOI, VOETOI, QUETOI,
     %             NBTEDS, MXTETR, NOTEDS, VOLUMT, QUALIT,
     %             VOLET0, QUAET0, VOLET1, QUAET1 )
      IF( QUAET1 .LE. 0 ) GOTO 9900
      IF( QUAET1 .GT. QUAET0 ) THEN

C        GENERATION EFFECTIVE DE CETTE TETRAEDRISATION ETOILEE
         CALL GETEPT( P,      NBSOMM, MXSOMM, XYZSOM, NBDM, NUDMEF,
     %                NBSOTE, N1TEVI, NSTETR,
     %                NO1TSO, N1TESO, NOTESO,
     %                N1FEOC, N1FEVI, NFETOI, NBTEDS, NOTEDS,
     %                MXTETA, VOLUMT, QUALIT,
     %                VOLET0, QUAET0, VOLET1, QUAET1, IERR )
         WRITE(IMPRIM,*) 'qua1te 5: NT=',NT,' getept execute'


cccC        DETECTION DES FACES APPARTENANT A 3 TETRAEDRES
ccc         CALL CRHAFAVO( NBSOTE, MXTETR, NSTETR, NUDTETR,
ccc     %                  MOFACE, MXFACE, LFACES, NBFAPB, IERR)

ccc         IF( NBFAPB .GT. 0 ) THEN
ccc            PRINT *
ccc            PRINT *,'qua1te 5: 1 NT=',NT,' avec',
ccc     %               NBFAPB,' FACES DE 3 TETRAEDRES'
ccc            TRACTE = .TRUE.
ccc         ENDIF


         GOTO 9900

      ENDIF

ccc      PRINT *,'qua1te: TETRAEDRE',NT,' NON MODIFIE... St',
ccc     %         (NSTETR(K,NT),K=1,4),' QUALITE=',QUALIT(NT),
ccc     %        ' NTOP=',NTOP,' NSOP=',NSOP

ccc      PRINT *,'qua1te: QUALITE MINIMALE DES TETRAEDRES OPPOSES=',QUALIV

C     FIN DU TRAITEMENT DU TETRAEDRE NT
C     =================================

C     LES FACES OCCUPEES DE NFETOI DEVIENNENT VIDES
 9900 CALL ETOIOCVI( N1FEOC, N1FEVI, NFETOI )


ccc      KTITRE = 'qua1te: TETRAEDRE           +   TETRAEDRES OPPOSES et OP
ccc     %POSES APRES TRAITEMENT'
ccc      WRITE(KTITRE(19:27),'(I9)' ) NT
ccc      WRITE(KTITRE(31:31),'(I1)' ) NFR
ccc      CALL SANSDBL( KTITRE, L )
ccc      NOPASS = 1
ccc      TRACTE = .TRUE.
ccc      CALL TR1TEVV( NOPASS, KTITRE(1:L), XYZSOM, NT, NBSOTE, NSTETR,
ccc     %              NO1TSO, NOTESO, NBTETRA, NOTETRA, VOLET1 )
ccc      IF( VOLET0 .NE. 0 .AND.
ccc     %    ABS(VOLET1-VOLET0) .GT. VOLET0*0.001 ) THEN
ccc         PRINT *,'qua1te: VOLUMES DIFFERENTS V0=',VOLET0,' V1=',VOLET1,
ccc     %           'GAIN=',ABS((VOLET1-VOLET0)/VOLET0)*100,'%'
ccc      ENDIF
ccc      TRACTE = .FALSE.

cccC     CONTROLE DE REGRESSION DE LA QUALITE MINIMALE DU MAILLAGE
ccc      CALL QUAMESH( 'FIN Tetra ', NT, GRAND,  XYZSOM,
ccc     %               NBSOTE, MXTETR, NSTETR,
ccc     %               VOLUMT, QUALIT, QUAMIN0,
ccc     %               NBTU,   VOLT,   QUAMOY, QUAMIN, NBPB )
ccc      nbpb = 0
ccc      if( nbpb .ne. 0 ) then
ccc         print *,'qua1te: Regression qualite avec NT=',NT,
ccc     %           ' QUAMIN=',QUAMIN
ccc         QUAMIN0 = QUAMIN
ccc      endif

ccc      if( QUAET1 .LE. QUAMIN0 ) then
ccc         print *,'qua1te: REGRESSION de la QUALITE avec NT=',NT,
ccc     %           ' QUAMIN0=',QUAMIN0,' QUAET1=',QUAET1,' NTOP=',NTOP
ccc      KTITRE = 'qua1te: Regression de la qualite au tetraedre         '
ccc      WRITE(KTITRE(47:55),'(I9)' ) NT
ccc      CALL SANSDBL( KTITRE, L )
ccc      NOPASS = 1
ccc      TRACTE = .TRUE.
ccc      CALL TR1TEVV( NOPASS, KTITRE(1:L), XYZSOM, NT, NBSOTE, NSTETR,
ccc     %              NO1TSO, NOTESO, NBTETRA, NOTETRA, VOLET1 )
ccc      TRACTE = .FALSE.
ccc      endif

      TRACTE = TRACTE0

      RETURN
      END
