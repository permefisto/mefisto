      SUBROUTINE QTEVNUL( NT,     NFR,    NTOP,   NSOP, SURFTR,VOLMOYTE,
     %                    MOFACE, MXFACE, LFACES, NBFAPB,
     %                    NBSOMM, MXSOMM, XYZSOM, GRAND,  NBDM,  NUDMEF,
     %                    NBSOTE, MXTETR, N1TEVI, NSTETR, NPSOFR,
     %                    VOLUMT, QUALIT,
     %                    NO1TSO, MXTESO, N1TESO, NOTESO,
     %                    MXFAET, N1FEOC, N1FEVI, NFETOI, VOETOI,QUETOI,
     %                    NBTEDS, MXTEDS, NOTEDS, MXTETA,
     %                    VOLET1, QUAET1, IERR  )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   AMELIORER LA QUALITE DU TETRAEDRE NT DE VOLUME QUASI NUL ou <0
C -----   ou de QUALITE<0.12

C ENTREES:
C --------
C NT     : NUMERO DANS NSTETR DU TETRAEDRE
C NFR    : NOMBRE DE TETRAEDRES OPPOSES PAR LES FACES
C NTOP   : NUMERO DES 4 TETRAEDRES OPPOSES PAR LES 4 FACES
C          0 SI FACE FRONTALIERE
C NSOP   : NUMERO DU SOMMET OPPOSE DES TETRAEDRES OPPOSES AUX 4 FACES DE NT
C SURFTR : SURFACE DES 4 FACES DU TETRAEDRE
C VOLMOYTE: VOLUME MOYEN D'UN TETRAEDRE DE LA TETREDRISATION
C MXSOMM : NOMBRE MAXIMAL DE SOMMETS DECLARABLES DANS XYZSOM ET NPSOFR
C XYZSOM : COORDONNEES X Y Z DES NBSOMM SOMMETS DES TETRAEDRES
C GRAND  : LE PLUS GRAND REEL STOCKABLE
C NBDM   : 0 SI 1 MATERIAU=VOLUME, SINON NOMBRE DE MATERIAUX DU VOLUME
C NBSOTE : VALEUR MAXIMALE DE DECLARATION DU PREMIER INDICE DE NSTETR(>3)
C MXTETR : NOMBRE MAXIMAL DE TETRAEDRES DECLARABLES DANS NSTETR
C N1TEVI : NUMERO 1-ER TETRAEDRE VIDE DE NSTETR
C MXTESO : NOMBRE MAXIMAL DE NUMERO DE TETRAEDRES DES SOMMETS
C NPSOFR : NUMERO 0 SI SOMMET INTERNE
C                 1 SI SOMMET SUR LA FRONTIERE
C                 2 SI SOMMET SUR INTERFACE ENTRE 2 MATERIAUX
C                -1 SI SOMMET SUPPRIME LORS DE L'AMELIORATION

C MODIFIES :
C ----------
C NBSOMM : NOMBRE DE SOMMETS AVANT ET APRES
C MXTETA : NUMERO MAXIMAL DES TETRAEDRES ACTUELS
C NUDMEF : NUMERO DE MATERIAU DE CHAQUE TETRAEDRE DU MAILLAGE
C          ATTENTION: CE TABLEAU EXISTE SEULEMENT SI NBDM>0
C NSTETR : NUMERO DES 4 SOMMETS DE CHAQUE TETRAEDRE

C TABLEAUX AUXILIAIRES :
C ----------------------
C MXFAET : NOMBRE MAXIMAL DE FACES DECLARABLES DANS NFETOI VOETOI QUETOI
C NFETOI(5,MXFAET)  DES ENTIERS
C VOETOI(MXFAET)    DES REELS
C QUETOI(MXFAET)    DES REELS
C MXTEDS : NOMBRE MAXIMAL DE TETRAEDRES RANGES DANS NOTEDS
C NOTEDS : NUMERO DANS NSTETR DES TETRAEDRES A DETRUIRE

C SORTIES:
C --------
C NBTEDS : NOMBRE DE TETRAEDRES RANGES DANS NOTEDS
C NO1TSO : NUMERO DU 1-ER TETRAEDRE DANS NOTESO DE CHAQUE SOMMET
C NOTESO : NOTESO(1,*) NUMERO DANS NSTETR DU TETRAEDRE
C          NOTESO(2,*) NUMERO DANS NOTESO DU TETRAEDRE SUIVANT
C                      0 SI C'EST LE DERNIER
C QUALIT : QUALITE DES TETRAEDRES DE LA TETRAEDRISATION
C VOLUMT : VOLUME  DES TETRAEDRES DE LA TETRAEDRISATION
C VOLET1 : VOLUME  DE LA DERNIERE ETOILE CALCULEE
C QUAET1 : QUALITE DE LA DERNIERE ETOILE CALCULEE
C IERR   : = 0 SI PAS D'ERREUR ET PAS DE MODIFICATION
C          =-1 SI AMELIORATION DE NT EFFECTUEE CORRECTEMENT
C          > 0 EN CAS DE SATURATION D'UN TABLEAU
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC    JANVIER 1992
C MODIFS : ALAIN PERRONNET LJLL UPMC & ST PIERRE DU PERRAY  AOUT    2010
C MODIFS : ALAIN PERRONNET LJLL UPMC et St Pierre du Perray Fevrier 2016
C23456...............................................................012
      PARAMETER         (MXSSET=256)
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      COMMON / TRTETR / STOPTE,TRACTE
      LOGICAL           STOPTE,TRACTE,TRACTE0

      INTEGER           NSTETR(NBSOTE,MXTETR), LFACES(NBSOTE,MXFACE),
     %                  NPSOFR(MXSOMM),
     %                  NUDMEF(MXTETR),
     %                  NO1TSO(MXSOMM),
     %                  NOTESO(2,MXTESO),
     %                  NOTEDS(MXTEDS),
     %                  NFETOI(5,MXFAET)
      REAL              XYZSOM(3,MXSOMM), VOLUMT(MXTETR),
     %                  QUALIT(MXTETR), VOETOI(MXFAET), QUETOI(MXFAET)

      REAL              SURFTR(4), SURFT1(4),
     %                  P(3), PMX(3), BAR(3), PB(3), VN1(3), VN2(3),
     %                  VOLET0, VOLET1, QUAET0, QUAET1
      DOUBLE PRECISION  PTXYZD(4,3), PD(3), COBARY(3),
     %                  VMOYEN, VOLTER, V
      INTEGER           N1SSET(MXSSET), NSF(3), NSOP(4), NTOP(4)
      INTEGER           NOTETRA(32)

      INTEGER           NOSOTR(3)
      EQUIVALENCE      (NOSOTR(1),NS2),(NOSOTR(2),NS3),(NOSOTR(3),NS4)

C     LES 3 SOMMETS DE LA FACE DU TETRAEDRE SONT ORIENTES POUR QUE
C     LE VOLUME AVEC LE SOMMET MANQUANT DE LA FACE SOIT POSITIF
      INTEGER           NOSOFATE(3,4)
      DATA              NOSOFATE / 1,2,3, 2,4,3, 3,4,1, 4,2,1 /

      INTEGER           NOSOARTE(2,6)
      DATA              NOSOARTE / 1,2, 2,3, 3,1, 4,1, 4,2, 4,3 /
      CHARACTER*128     KTITRE

      PRINT*
      PRINT*,'qtevnul: Entree TETRAEDRE',NT,': St',(NSTETR(I,NT),I=1,4),
     %       ' V=',VOLUMT(NT),' Q=',QUALIT(NT),' <---------------------'

C     TRACE DU TETRAEDRE NT et SES TETRAEDRES VOISINS et VOISINS VOISINS
C     ------------------------------------------------------------------
      TRACTE0 = TRACTE
      KTITRE='qtevnul: TETRAEDRE            V=                Q=        
     %       '
      WRITE(KTITRE(20:28),'(I9)' ) NT
      WRITE(KTITRE(33:47),'(G15.6)' ) VOLUMT(NT)
      WRITE(KTITRE(51:65),'(G15.6)' ) QUALIT(NT)
      CALL SANSDBL( KTITRE, L )
      PRINT*,KTITRE(1:L)
      NOPASS = 0

      IF( QUALIT(NT) .LT. 0.01 ) THEN
         TRACTE = .TRUE.
         CALL TR1TEVV( NOPASS, KTITRE(1:L), XYZSOM, NT, NBSOTE, NSTETR,
     %                 NO1TSO, NOTESO, NBTETRA, NOTETRA, VOLET0 )
         TRACTE = TRACTE0
      ENDIF
            
C     RECHERCHE DU CAS DE DEGENERESCENCE DU TETRAEDRE
C     CAS 1 (TYPE QUADRANGLE) OU 2 (TRIANGLE AVEC BARYCENTRE)?
C     ========================================================
      IERR = 0
      NUMATE  = 0
      NBTETRA = 0
      NBTEDS  = 0

C     AIRE MIN ET MAX DES 4 FACES DU TETREDRE NT
      NFSFMAX = 1
      SFMAX   = SURFTR(1)
      DO I=2,4
         IF( SURFTR(I) .GT. SFMAX ) THEN
            SFMAX   = SURFTR(I)
            NFSFMAX = I
         ENDIF
      ENDDO

C     LE SOMMET OPPOSE A LA PLUS GRANDE FACE NFSFMAX DE NT
      IF( NFSFMAX .GT. 1 ) THEN
         IS = NFSFMAX - 1
      ELSE
         IS = 4
      ENDIF
      S3 = SURFTR(1)+SURFTR(2)+SURFTR(3)+SURFTR(4) - SFMAX

C     NS1 LE NUMERO DU SOMMET OPPOSE A LA PLUS GRANDE FACE NFSFMAX DE NT
      NS1 = NSTETR(IS,NT)

C     LES 3 SOMMETS DE LA PLUS GRANDE FACE NFSFMAX DE NT
C     DE NORMALE DIRIGEE VERS L'INTERIEUR DU TETRAEDRE NT
C     SONT NS2 NS3 NS4 => VOLUME (NS1234) >=0
      NS2 = NOSOFATE( 1, NFSFMAX )
      NS2 = NSTETR( NS2, NT )

      NS3 = NOSOFATE( 2, NFSFMAX )
      NS3 = NSTETR( NS3, NT )

      NS4 = NOSOFATE( 3, NFSFMAX )
      NS4 = NSTETR( NS4, NT )

C     PROJECTION DU SOMMET OPPOSE SUR LA FACE NFSFMAX
      CALL PRPTPL( XYZSOM(1,NS1),
     %             XYZSOM(1,NS2), XYZSOM(1,NS3), XYZSOM(1,NS4),
     %             P, IERR )
      IF( IERR .NE. 0 ) THEN
         IERR = 0
         GOTO 9999
      ENDIF
C
C     LE POINT PROJECTION DU SOMMET NS1 EST IL DANS LE TRIANGLE?
      DO I=1,3
         PD(I)       = P(I)
         PTXYZD(I,1) = XYZSOM(I,NS2)
         PTXYZD(I,2) = XYZSOM(I,NS3)
         PTXYZD(I,3) = XYZSOM(I,NS4)
         NSF(I) = I
      ENDDO
      CALL PTDSTR( PD, PTXYZD, NSF, COBARY, NONOUI )
C     NONOUI: 1 SI LE POINT EST DANS OU SUR LE TRIANGLE
C             0 SI LE POINT P EST EXTERIEUR AU TRIANGLE
C            -1 SI LE TRIANGLE EST DEGENERE


      IF( NONOUI .EQ. 0 .OR. ABS(S3-SFMAX) .GT. SFMAX*0.1 ) THEN

C        ====================================
C        TYPE 1: TETRAEDRE DE TYPE QUADRANGLE
C        ====================================
         IF( NFR .EQ. 0 ) THEN

C           0 TETRAEDRE OPPOSE => SUPPRESSION DU TETRAEDRE NT
C           =================================================
            CALL DS1TET( NT,     NBSOTE, N1TEVI, NSTETR,
     %                   NO1TSO, N1TESO, NOTESO )
            WRITE(IMPRIM,*) 'qtevnul 1: ds1tet REUSSI NT=',NT,
     %          ' TETRAEDRE ISOLE DE VOLUME NUL DETRUIT'
            IERR = -1
            QUALIT( NT ) = GRAND
            VOLUMT( NT ) = 0
            IF( NBDM .GT. 0 ) NUDMEF( NT ) = 0
C           SUPPRESSION DU NO DE MATERIAU
            GOTO 9999


         ELSE IF( NFR .EQ. 1 ) THEN

C           NT DE VOLUME QUASI NUL ET 1 TETRAEDRE OPPOSE
C           ============================================
            WRITE(IMPRIM,*) 'qtevnul 2: TETRAEDRE QUADRANGLE et 1 SEUL T
     %ETRAEDRE OPPOSE NT=',NT,' CAS IMPROBABLE?...'
            WRITE(IMPRIM,*) 'TETRAEDRE',NT,': ST',(NSTETR(I,NT),I=1,4)
            WRITE(IMPRIM,*) 'NTOP=',NTOP,' NSOP=',NSOP
            WRITE(IMPRIM,*) 'RETOUR SANS MODIFICATION a ANALYSER...'

            KTITRE='qtevnul: TETRAEDRE            avec 1 SEUL TETRAEDRE 
     %OPPOSE. O MODIF'
            WRITE(KTITRE(19:27),'(I9)' ) NT
            CALL SANSDBL( KTITRE, L )
            NOPASS = 0
            TRACTE = .TRUE.
          CALL TR1TEVV( NOPASS, KTITRE(1:L), XYZSOM, NT, NBSOTE, NSTETR,
     %                  NO1TSO, NOTESO, NBTETRA, NOTETRA, VOLET0 )
            TRACTE = .FALSE.
            GOTO 9999


         ELSE IF( NFR .EQ. 2 ) THEN
C
C           NT DE VOLUME QUASI NUL ET 2 TETRAEDRES OPPOSES
C           ==============================================
C           SI UN SOMMET DOUBLE POUR LES 2 TETRAEDRES OPPOSES ALORS
C           GENERER 2 TETRAEDRES A PARTIR DES 3 TETRAEDRES
            CALL Q3T2S2( NT,     NTOP,   NSOP,   GRAND,
     %                   XYZSOM, NBSOTE, N1TEVI, NSTETR, NBDM, NUDMEF,
     %                   NO1TSO, N1TESO, NOTESO,
     %                   MXTETA, VOLUMT, QUALIT, IERR )
            IF( IERR .EQ. -1 ) THEN
               WRITE(IMPRIM,*)'qtevnul 3: q3t2s2 REUSSI NT=',NT,
     %                        ' 2 TETRA OPPOSES NTOP=',NTOP
               GOTO 9999
            ENDIF

C           ESSAYER LE BARYCENTRE DES TETRAEDRES D'ARETE CELLE
C           COMMUNE AUX 2 TETRAEDRES OPPOSES AU TETRAEDRE NT
C           SI MEILLEURE QUALITE ALORS
C           TETRAEDRISER L'ETOILE DES TETRAEDRES DE CETTE ARETE
            CALL QTE2TOP( NBSOMM, MXSOMM, XYZSOM, NBDM, NUDMEF,
     %                    NT,     NTOP,   NBSOTE, N1TEVI, NSTETR,
     %                    NO1TSO, N1TESO, NOTESO,
     %                    N1FEOC, N1FEVI, NFETOI, VOETOI, QUETOI,
     %                    MXTETA, VOLUMT, QUALIT, NOTEDS,
     %                    VOLET0, QUAET0, VOLET1, QUAET1, IERR  )
            IF( IERR .EQ. -1 ) THEN
               WRITE(IMPRIM,*) 'qtevnul 4: qte2top REUSSI AVEC NT=',NT,
     %                         ' 2 TETRA OPPOSES NTOP=',NTOP
               GOTO 9999
            ENDIF

C           RECHERCHE DE L'ARETE COMMUNE AUX TETRAEDRES OPPOSES
C           ESSAI DE SUPPRESSION DE CETTE ARETE NAR DU TETRAEDRE NT
            DO I1=1,4
               IF( NTOP(I1) .NE. 0 ) GOTO 12
            ENDDO
 12         DO I2=I1+1,4
               IF( NTOP(I2) .NE. 0 ) GOTO 20
            ENDDO
C           NUMERO DE L'ARETE DE NT COMMUNE AUX 2 TETRAEDRES OPPOSES
 20         CALL NAC2FATE( I1, I2, NAR )
C           ESSAI DE SUPPRESSION DE CETTE ARETE NAR DU TETRAEDRE NT
            CALL QUMTNT( NT,     NAR,
     %                   XYZSOM, NBSOTE, N1TEVI, NSTETR, NBDM,NUDMEF,
     %                   NO1TSO, N1TESO, NOTESO,
     %                   NBTEDS, NOTEDS, QUETOI,
     %                   MXTETA, VOLUMT, QUALIT, IERR )
            IF( IERR .EQ. -1 ) THEN
C              REUSSITE
              WRITE(IMPRIM,*)'qtevnul 5: qumtnt REUSSI QUADRANGLE 2 TETR
     %A OPPOSES NT=',NT,' NAR=',NAR,' NBTEDS=',NBTEDS
               GOTO 9999
            ENDIF

C           CONSTRUIRE L'ETOILE DES TETRAEDRES D'ARETE LA PLUS LONGUE
C           ARETE DU TETRAEDRE NT
C           CHERCHER LE MEILLEUR BARYCENTRE DE CES TETRAEDRES
C           CREER CETTE TETRAEDRISATION SI VOLUME RESPECTE ET
C           QUALITE AMELIOREE DE L'ETOILE
            CALL QTEN4T( NBSOMM, MXSOMM, XYZSOM, NBDM,   NUDMEF,
     %                   NT,     NBSOTE, N1TEVI, NSTETR,
     %                   NO1TSO, N1TESO, NOTESO,
     %                   N1FEOC, N1FEVI, NFETOI, VOETOI, QUETOI,
     %                   MXTETA, VOLUMT, QUALIT, NOTEDS,
     %                   VOLET0, QUAET0, VOLET1, QUAET1, IERR  )
            IF( IERR .EQ. -1 ) THEN
               WRITE(IMPRIM,*) 'qtevnul 6: qten4t REUSSI AVEC NT=',NT,
     %                         ' 2 TETRA OPPOSES NTOP=',NTOP
               GOTO 9999
            ENDIF

C           CONSTRUIRE L'ETOILE DES TETRAEDRES DE CHACUNE DES 6 ARETES
C           DU TETRAEDRE NT SI ELLE EST INTERNE AU VOLUME
C           CHERCHER LE MEILLEUR BARYCENTRE DE CES TETRAEDRES
C           CREER CETTE TETRAEDRISATION SI VOLUME RESPECTE ET
C           QUALITE AMELIOREE DE L'ETOILE
            CALL QTE6AR( NBSOMM, MXSOMM, XYZSOM, NPSOFR, NBDM, NUDMEF,
     %                   NT,     NBSOTE, N1TEVI, NSTETR,
     %                   NO1TSO, N1TESO, NOTESO,
     %                   N1FEOC, N1FEVI, NFETOI, VOETOI, QUETOI,
     %                   MXTETA, VOLUMT, QUALIT, NOTEDS,
     %                   VOLET0, QUAET0, VOLET1, QUAET1, IERR  )
            IF( IERR .EQ. -1 ) THEN
               WRITE(IMPRIM,*) 'qtevnul 7: qte6ar REUSSI AVEC NT=',NT,
     %                         ' 2 TETRA OPPOSES NTOP=',NTOP
               GOTO 9999
            ENDIF

ccc            IF( QUALIT( NT ) .LT. 0.05 ) THEN
cccC              QUADRANGLE PLAT => SUPPRESSION DU TETRAEDRE NT
ccc               CALL DS1TET( NT,     NBSOTE, N1TEVI, NSTETR,
ccc     %                      NO1TSO, N1TESO, NOTESO )
ccc              WRITE(IMPRIM,*)'qtevnul 8: ds1tet REUSSI QUADRANGLE 2 TETR
ccc     %A OPPOSES NT=',NT,' TETRAEDRE DETRUIT car Qualite=',QUALIT(NT),
ccc     %                 ' Volume=',VOLUMT(NT)
ccc               IERR = -1
ccc               QUALIT( NT ) = GRAND
ccc               VOLUMT( NT ) = 0
ccc               IF( NBDM .GT. 0 ) NUDMEF( NT ) = 0
cccC              SUPPRESSION DU NO DE MATERIAU
ccc               GOTO 9999
ccc            ENDIF

C           PAS D'AMELIORATION....
           WRITE(IMPRIM,*)'qtevnul 9: PAS D''AMELIORATION. QUADRANGLE 2 
     %TETRA OPPOSES NT=',NT,' Qualite=',QUALIT(NT),' Volume=',VOLUMT(NT)
 
            KTITRE ='qtevnul: TETRAEDRE           + 2 TETRAEDRE OPPOSES 
     %et PAS d''AMELIORATION Qualite=           '
            WRITE(KTITRE(19:27),'(I9)' ) NT
            WRITE(KTITRE(85:95),'(F11.7)' ) QUALIT(NT)
            CALL SANSDBL( KTITRE, L )
            NOPASS = 0
            TRACTE = .TRUE.
            CALL TR1TEVV(NOPASS, KTITRE(1:L), XYZSOM, NT, NBSOTE,NSTETR,
     %                   NO1TSO, NOTESO, NBTETRA, NOTETRA, VOLET0 )
            TRACTE = TRACTE0
            IERR = 0
            GOTO 1010


         ELSE IF( NFR .EQ. 3 ) THEN
C
C           NT DE VOLUME QUASI NUL ET 3 TETRAEDRES OPPOSES
C           ==============================================
C           GENERER 2 TETRAEDRES A PARTIR DE 3 TETRAEDRES ET UN SOMMET
C           DOUBLE OPPOSE A UN TETRAEDRE DE MAUVAISE QUALITE?
C           ..........................................................
            CALL Q3T2S2( NT,     NTOP,   NSOP,   GRAND,
     %                   XYZSOM, NBSOTE, N1TEVI, NSTETR, NBDM, NUDMEF,
     %                   NO1TSO, N1TESO, NOTESO,
     %                   MXTETA, VOLUMT, QUALIT, IERR )
            IF( IERR .EQ. -1 ) THEN
               WRITE(IMPRIM,*)'qtevnul 10: q3t2s2 REUSSI NT=',NT,
     %                        ' 3 TETRA OPPOSES NTOP=',NTOP
               GOTO 9999
            ENDIF

C           RECHERCHE DE LA FACE FRONTALIERE I1 DE NT
            DO 22 I1=1,4
               IF( NTOP(I1) .EQ. 0 ) GOTO 24
 22         CONTINUE

C           RECHERCHE DES 3 ARETES N'APPARTENANT PAS A LA FACE FRONTALIERE
 24         IF( I1 .EQ. 1 ) THEN
C              LES ARETES N'APPARTENANT PAS A LA FACE 1
               NSF(1) = 4
               NSF(2) = 5
               NSF(3) = 6
            ELSE IF( I1 .EQ. 2 ) THEN
C              LES ARETES N'APPARTENANT PAS A LA FACE 2
               NSF(1) = 1
               NSF(2) = 3
               NSF(3) = 4
            ELSE IF( I1 .EQ. 3 ) THEN
C              LES ARETES N'APPARTENANT PAS A LA FACE 3
               NSF(1) = 1
               NSF(2) = 2
               NSF(3) = 5
            ELSE
C              LES ARETES N'APPARTENANT PAS A LA FACE 4
               NSF(1) = 2
               NSF(2) = 3
               NSF(3) = 6
            ENDIF

            DO I=1,3

C              ESSAI DE SUPPRESSION DE CETTE ARETE NSF(I) DU TETRAEDRE NT
C              ..........................................................
               CALL QUMTNT( NT,     NSF(I),
     %                      XYZSOM, NBSOTE, N1TEVI, NSTETR, NBDM,NUDMEF,
     %                      NO1TSO, N1TESO, NOTESO,
     %                      NBTEDS, NOTEDS, QUETOI,
     %                      MXTETA, VOLUMT, QUALIT, IERR )
               IF( IERR .EQ. -1 ) THEN
C                 REUSSITE
                 WRITE(IMPRIM,*)'qtevnul 11: qumtnt REUSSI 3 TETRA OPPOS
     %ES I=',I,' NT=',NT,' 3 TETRA OPPOSES NSF(',I,')=',NSF(I)
                  GOTO 9999
               ENDIF

            ENDDO

C           IL FAUT FORCER LE BARYCENTRE DU TETRAEDRE NT?
C            ............................................
            CALL QTEN4T( NBSOMM, MXSOMM, XYZSOM, NBDM, NUDMEF,
     %                   NT,     NBSOTE, N1TEVI, NSTETR,
     %                   NO1TSO, N1TESO, NOTESO,
     %                   N1FEOC, N1FEVI, NFETOI, VOETOI, QUETOI,
     %                   MXTETA, VOLUMT, QUALIT, NOTEDS,
     %                   VOLET0, QUAET0, VOLET1, QUAET1, IERR  )
            IF( IERR .EQ. -1 ) THEN
               WRITE(IMPRIM,*) 'qtevnul 12: qten4t REUSSI AVEC NT=',NT,
     %                         ' 3 TETRA OPPOSES NTOP=',NTOP
               GOTO 9999
            ENDIF

C           CONSTRUIRE L'ETOILE DES TETRAEDRES DE CHACUNE DES 6 ARETES
C           DU TETRAEDRE NT SI ELLE EST INTERNE AU VOLUME
C           CHERCHER LE MEILLEUR BARYCENTRE DE CES TETRAEDRES
C           CREER CETTE TETRAEDRISATION SI VOLUME RESPECTE ET
C           QUALITE AMELIOREE DE L'ETOILE
C           ..........................................................
            CALL QTE6AR( NBSOMM, MXSOMM, XYZSOM, NPSOFR, NBDM, NUDMEF,
     %                   NT,     NBSOTE, N1TEVI, NSTETR,
     %                   NO1TSO, N1TESO, NOTESO,
     %                   N1FEOC, N1FEVI, NFETOI, VOETOI, QUETOI,
     %                   MXTETA, VOLUMT, QUALIT, NOTEDS,
     %                   VOLET0, QUAET0, VOLET1, QUAET1, IERR  )
            IF( IERR .EQ. -1 ) THEN
               WRITE(IMPRIM,*) 'qtevnul 13: qte6ar REUSSI AVEC NT=',NT,
     %                         ' 3 TETRA OPPOSES NTOP=',NTOP
               GOTO 9999
            ENDIF


         ELSE
C
C           NT DE VOLUME QUASI NUL ET 4 TETRAEDRES OPPOSES
C           ==============================================
C           GENERER 2 TETRAEDRES A PARTIR DE 3 TETRAEDRES ET UN SOMMET
C           DOUBLE OPPOSE A UN TETRAEDRE DE MAUVAISE QUALITE
C           ..........................................................
            CALL Q3T2S2( NT,     NTOP,   NSOP,   GRAND,
     %                   XYZSOM, NBSOTE, N1TEVI, NSTETR, NBDM, NUDMEF,
     %                   NO1TSO, N1TESO, NOTESO,
     %                   MXTETA, VOLUMT, QUALIT, IERR )
            IF( IERR .EQ. -1 ) THEN
               WRITE(IMPRIM,*)'qtevnul 14: q3t2s2 REUSSI NT=',NT,
     %         ' 4 TETRA OPPOSES NTOP=',NTOP,' NSOP=',NSOP
               GOTO 9999
            ENDIF

            IF( IERR .EQ. 1 ) THEN
C              PAS DE SOMMET DOUBLE OPPOSE A 2 FACES DE NT
               IERR = 0

ccc      KTITRE = 'qtevnul: TETRAEDRE           + 4 TETRAEDRES OPPOSES et O
ccc     %PPOSES apres q3t2s2'
ccc      WRITE(KTITRE(19:27),'(I9)' ) NT
ccc      CALL SANSDBL( KTITRE, L )
ccc      NOPASS = 0
ccc      TRACTE = .TRUE.
ccc      CALL TR1TEVV( NOPASS, KTITRE(1:L), XYZSOM, NT, NBSOTE, NSTETR,
ccc     %              NO1TSO, NOTESO, NBTETRA, NOTETRA, VOLET0 )

C              RECHERCHE DES 2 DIAGONALES DU QUADRANGLE FACE MAX
C              ET DU POINT PROJETE P OU PLUS PRECISEMENT
C              RECHERCHE DE L'ARETE DU TRIANGLE FACE MAX
C              DU MEME COTE QUE LE POINT P DE PROJECTION
C              -------------------------------------------------
C              NORMALE A LA FACE DE NT DE PLUS GRANDE SURFACE
               CALL NORFA3( XYZSOM(1,NOSOTR(1)),
     %                      XYZSOM(1,NOSOTR(2)),
     %                      XYZSOM(1,NOSOTR(3)),
     %                      VN1, IERR )

C              NORMALE A P + ARETE K DE LA FACE SE SURFACE MAX
               DO K=1,3
                  IF( K .EQ. 3 ) THEN
                     K1 = 1
                  ELSE
                     K1 = K+1
                  ENDIF

                  NS1D = NOSOTR(K)
                  NS2D = NOSOTR(K1)

                  CALL NORFA3( P, XYZSOM(1,NS1D), XYZSOM(1,NS2D),
     %                         VN2, IERR )

C                 PRODUIT SCALAIRE DES 2 NORMALES
                  IF( VN1(1)*VN2(1)+VN1(2)*VN2(2)+VN1(3)*VN2(3) .LT. 0 )
     %                GOTO 31
               ENDDO
               GOTO 36

C              P EST DE L'AUTRE COTE DE L'ARETE K DE LA FACE MAX
C              => L'ARETE K NS1D-NS2D EST UNE DIAGONALE DU QUADRANGLE
C              NUMERO DE L'ARETE DANS LE TETRAEDRE?
 31            DO NAR=1,6
C                 NUMERO DES 2 SOMMETS DE L'ARETE I D'UN TETRAEDRE
                  NS1A = NOSOARTE( 1,  NAR )
                  NS1A = NSTETR( NS1A, NT )
                  NS2A = NOSOARTE( 2,  NAR )
                  NS2A = NSTETR( NS2A, NT )
                  IF( ( NS1A .EQ. NS1D .AND. NS2A .EQ. NS2D ) .OR.
     %                ( NS1A .EQ. NS2D .AND. NS2A .EQ. NS1D ) ) GOTO 34
               ENDDO
               GOTO 36

cccC              LA 1-ERE DIAGONALE NS1D-NS2D EST L'ARETE NAR DU TETRAEDRE NT
cccC              ESSAI D'ECHANGER DANS L'ETOILE UNE DIAGONALE PAR L'AUTRE
cccC              SANS AJOUTER DE NOUVEAU SOMMET A LA TETRAEDRISATION C-A-D
cccC              TETRAEDRISER EN SUPPRIMANT LA DIAGONALE
cccC              ------------------------------------------------------------
ccc 33            CALL SITEDIAG( NAR,    NT,     XYZSOM,
ccc     %                        NBSOTE, NSTETR, NO1TSO, NOTESO,
ccc     %                        NBTEDS, MXTETR, NOTEDS,
ccc     %                        QUALIT, VOLUET0, QUAET0, VOLUET1, QUAET1,
ccc     %                        P,      NESSAI )
ccc               IF( QUAET1 .GT. 0 .AND. QUAET1 .GT. QUAET0 ) THEN
cccC                 GENERATION EFFECTIVE DE CETTE TETRAEDRISATION ETOILEE
cccC                 AVEC SUPPRESSION DE L'UNE DES 2 DIAGONALES
cccC                 -----------------------------------------------------
ccc                  CALL GETEDIAG( NESSAI, NBSOMM, XYZSOM, NBDM,   NUDMEF,
ccc     %                           NAR,    NT,     NBSOTE, N1TEVI, NSTETR,
ccc     %                           NO1TSO, N1TESO, NOTESO,
ccc     %                           MXTEDS, NBTEDS, NOTEDS,
ccc     %                           MXTETA, VOLUMT, QUALIT,
ccc     %                           VOLET0, QUAET0, VOLET1, QUAET1, IERR )
ccc                  IF( IERR .EQ. 0 ) THEN
ccc                     IERR = -1
ccc                     GOTO 9999
ccc                  ENDIF
ccc                  IERR = 0
ccc               ENDIF

C              SIMULATION DE LA TETRAEDRISATION ETOILEE PAR LE BARYCENTRE
C              DU TETRAEDRE NT A PARTIR DE L'ARETE DIAGONALE I ET
C              LES TETRAEDRES DE L'AUTRE DIAGONALE
C              ..........................................................
 34            CALL SITEBA( NAR,    NT,     XYZSOM,
     %                      NBSOTE, NSTETR, NO1TSO, NOTESO,
     %                      N1FEOC, N1FEVI, NFETOI, VOETOI, QUETOI,
     %                      NBTEDS, MXTETR, NOTEDS,
     %                      QUALIT, VOLET0, QUAET0, VOLET1, QUAET1,
     %                      P,      NESSAI )

ccc               if( nbteds .GT. 0 ) then
ccc        KTITRE ='qtevnul: TETRAEDRE           + 4 TETRAEDRES OPPOSES ESS
ccc     %AI siteba'
ccc               WRITE(KTITRE(19:27),'(I9)' ) NT
ccc               CALL SANSDBL( KTITRE, L )
ccc               NOPASS = 1
ccc               TRACTE = .TRUE.
ccc               CALL TRNTETRA( NOPASS, KTITRE(1:L), XYZSOM, NBSOTE,
ccc     %                        NBTEDS, NOTEDS, NSTETR )
ccc               endif

               IF( QUAET1 .GT. 0 .AND. QUAET1 .GT. QUAET0 ) THEN

C                 GENERATION EFFECTIVE DE CETTE TETRAEDRISATION ETOILEE
C                 AVEC AJOUT DU BARYCENTRE DU TETRAEDRE NT
C                 .....................................................
                  CALL GETEPT( P,  NBSOMM, MXSOMM,XYZSOM,NBDM,NUDMEF,
     %                         NBSOTE, N1TEVI, NSTETR,
     %                         NO1TSO, N1TESO, NOTESO,
     %                         N1FEOC, N1FEVI, NFETOI, NBTEDS,NOTEDS,
     %                         MXTETA, VOLUMT, QUALIT,
     %                         VOLET0, QUAET0, VOLET1, QUAET1, IERR )
                  IF( IERR .EQ. 0 ) THEN
                     IERR = -1
                     GOTO 9999
                  ENDIF

               ENDIF

C              ESSAI DE FORCER LE BARYCENTRE DU TETRAEDRE NT
C              .............................................
 36            CALL QTEN4T( NBSOMM, MXSOMM, XYZSOM, NBDM, NUDMEF,
     %                      NT,     NBSOTE, N1TEVI, NSTETR,
     %                      NO1TSO, N1TESO, NOTESO,
     %                      N1FEOC, N1FEVI, NFETOI, VOETOI, QUETOI,
     %                      MXTETA, VOLUMT, QUALIT, NOTEDS,
     %                      VOLET0, QUAET0, VOLET1, QUAET1, IERR  )
               IF( IERR .EQ. -1 ) THEN
                  WRITE(IMPRIM,*) 'qtevnul 15: qten4t REUSSI NT=',NT,
     %                            ' 4 TETRA OPPOSES NTOP=',NTOP
                  GOTO 9999
               ENDIF

C              CONSTRUIRE L'ETOILE DES TETRAEDRES DE CHACUNE DES 6 ARETES
C              DU TETRAEDRE NT SI ELLE EST INTERNE AU VOLUME
C              CHERCHER LE MEILLEUR BARYCENTRE DE CES TETRAEDRES
C              CREER CETTE TETRAEDRISATION SI VOLUME RESPECTE ET
C              QUALITE AMELIOREE DE L'ETOILE
C              ..........................................................
               CALL QTE6AR( NBSOMM, MXSOMM, XYZSOM, NPSOFR, NBDM,NUDMEF,
     %                      NT,     NBSOTE, N1TEVI, NSTETR,
     %                      NO1TSO, N1TESO, NOTESO,
     %                      N1FEOC, N1FEVI, NFETOI, VOETOI, QUETOI,
     %                      MXTETA, VOLUMT, QUALIT, NOTEDS,
     %                      VOLET0, QUAET0, VOLET1, QUAET1, IERR  )
               IF( IERR .EQ. -1 ) THEN
                  WRITE(IMPRIM,*)'qtevnul 16: qte6ar REUSSI AVEC NT=',NT
     %                          ,' 4 TETRA OPPOSES NTOP=',NTOP
                  GOTO 9999
               ENDIF

            ENDIF

         ENDIF

C        FIN du TRAITEMENT TETRAEDRE NT DE TYPE QUADRANGLE
         GOTO 9999

      ENDIF


C     ==================================================
C     TYPE 2: TETRAEDRE DE TYPE TRIANGLE AVEC BARYCENTRE
C     ==================================================
ccc 40   
      IF( VOLUMT(NT) .LT. 0 ) THEN

C        TETRAEDRE NT DE VOLUME NEGATIF
C        DEPLACEMENT DU SOMMET NS1 DE L'AUTRE COTE DE SA PLUS GRANDE FACE
C        P EST LE POINT PROJETE DU SOMMET OPPOSE SUR LA FACE NFSFMAX
C        ----------------------------------------------------------------
C        SAUVEGARDE DANS PB
         DO K = 1, 3
            XYZ = XYZSOM( K, NS1 )
            PB( K ) = XYZ
            XYZSOM( K, NS1 ) = XYZ + 1.25 * ( P( K ) - XYZ )
         ENDDO

         CALL QUATET( XYZSOM(1,NSTETR(1,NT)),
     %                XYZSOM(1,NSTETR(2,NT)),
     %                XYZSOM(1,NSTETR(3,NT)),
     %                XYZSOM(1,NSTETR(4,NT)),
     %                ARMIN, ARMAX, SURFTR, VOLUMT(NT), QUALIT(NT) )
        PRINT*,'qtevnul 40: SOMMET NS1=',NS1,
     %         ' DEPLACE POUR RENDRE LE VOLUME',
     %          VOLUMT(NT),' >0  QUALITE=', QUALIT(NT)
        PRINT*,'qtevnul 40: NS1=',NS1,' AVANT',PB
        PRINT*,'qtevnul 40: NS1=',NS1,' APRES',(XYZSOM(K,NS1),K=1,3)

        IF( QUALIT(NT) .GE. 0.01 ) GOTO 9999

      ENDIF


      IF( NFR .EQ. 0 ) THEN

C           NT DE VOLUME QUASI NUL ET 0 TETRAEDRE OPPOSE
C           ============================================
C           => SUPPRESSION DU TETRAEDRE NT
C           ..............................
            WRITE(IMPRIM,*)'qtevnul 41: ds1tet NT=',NT,
     %      ' TETRAEDRE de VOLUME',VOLUMT(NT),' EST DETRUIT'
            CALL DS1TET( NT,     NBSOTE, N1TEVI, NSTETR,
     %                   NO1TSO, N1TESO, NOTESO )
            QUALIT( NT ) = GRAND
            VOLUMT( NT ) = 0
            IF( NBDM .GT. 0 ) NUDMEF( NT ) = 0
C           SUPPRESSION DU NO DE MATERIAU DE NT
            IERR = -1
            GOTO 9999


         ELSE IF( NFR .EQ. 1 ) THEN

C           NT DE VOLUME QUASI NUL ET 1 TETRAEDRE OPPOSE
C           ============================================
C           ESSAI 2341+2435  => 5124+5143+5132?
            IF( NTOP(NFSFMAX) .NE. 0 ) THEN

C              TETRAEDRE OPPOSE A LA PLUS GRANDE FACE
               NS5 = NSOP( NFSFMAX )
               IF( NS5 .EQ. NS1 .OR. NS5 .EQ. NS2 .OR. NS5 .EQ. NS3 .OR.
     %             NS5 .EQ. NS4 ) THEN
C                  LE SOMMET OPPOSE EST UN SOMMET DE NT
                  IERR = 0
                  GOTO 9999
               ENDIF

C              LES 2 TETRAEDRES SONT ILS DANS LE MEME MATERIAU? 
               NUMATE = 0
               IF( NBDM .GT. 0 ) THEN
                  NUMATE = NUDMEF(NT)
                  IF( NUMATE .NE. NUDMEF( NTOP(NFSFMAX) ) ) THEN
                     IERR = 0
                     GOTO 1000
                  ENDIF
               ENDIF

C              LA QUALITE DE 2431+2345 EST ELLE > QUALITE 5142+5134+5123 ?
C              CALCUL DU VOLUME ET QUALITE DU TETRAEDRE DU A LA FACE
               CALL QUATET( XYZSOM(1,NS5),
     %                      XYZSOM(1,NS1),
     %                      XYZSOM(1,NS4),
     %                      XYZSOM(1,NS2),

     %                      ARMIN, ARMAX, SURFT1, V1, Q1 )
ccc               IF( V1 .LE. 0 ) GOTO 44

               CALL QUATET( XYZSOM(1,NS5),
     %                      XYZSOM(1,NS1),
     %                      XYZSOM(1,NS3),
     %                      XYZSOM(1,NS4),
     %                      ARMIN, ARMAX, SURFT1, V2, Q2 )
ccc               IF( V2 .LE. 0 ) GOTO 44

               CALL QUATET( XYZSOM(1,NS5),
     %                      XYZSOM(1,NS1),
     %                      XYZSOM(1,NS2),
     %                      XYZSOM(1,NS3),
     %                      ARMIN, ARMAX, SURFT1, V3, Q3 )
ccc               IF( V3 .LE. 0 ) GOTO 44

               IF( MIN( QUALIT( NT ), QUALIT( NTOP(NFSFMAX) ) ) .GT.
     %             MIN( Q1, Q2, Q3 ) ) GOTO 44

C              SUPPRESSION DES 2 TETRAEDRES NT:2341 ET NT1:2345
               CALL DS1TET( NT,     NBSOTE, N1TEVI, NSTETR,
     %                      NO1TSO, N1TESO, NOTESO )
               QUALIT( NT ) = GRAND
               VOLUMT( NT ) = 0
               IF( NBDM .GT. 0 ) NUDMEF( NT ) = 0

               CALL DS1TET( NTOP(NFSFMAX),  NBSOTE, N1TEVI, NSTETR,
     %                      NO1TSO, N1TESO, NOTESO )
               QUALIT( NTOP(NFSFMAX) ) = GRAND
               VOLUMT( NTOP(NFSFMAX) ) = 0
               IF( NBDM .GT. 0 ) NUDMEF( NTOP(NFSFMAX) ) = 0

C              LES TETRAEDRES 2431+2345 DEVIENNENT 5142+5134+5123
               CALL Q05S3T( NS2, NS3, NS4, NS5, NS1, NUMATE, NUDMEF,
     %                      XYZSOM, NBSOTE, N1TEVI, NSTETR,
     %                      NO1TSO, N1TESO, NOTESO,
     %                      MXTETA, VOLUMT, QUALIT, IERR )

               WRITE(IMPRIM,*)'qtevnul 18: q05s3t REUSSI TYPE TRIANGLE 1
     % TETRA OP NT=',NT, ' St:',NS2, NS3, NS4, NS1, NS5


C              DETECTION DES FACES APPARTENANT A 3 TETRAEDRES (A SUPPRIMER...)
               CALL CRHAFAVO( NBSOTE, MXTETR, NSTETR, NUDTETR,
     %                        MOFACE, MXFACE, LFACES, NBFAPB, IERR )
               IF( NBFAPB .GT. 0 ) THEN
                  PRINT *
                  PRINT *,'qtevnul 18: q05s3t ATTENTION 1 NT=',NT,
     %                    ' avec',NBFAPB,' FACES DE 3 TETRAEDRES'
                  TRACTE = .TRUE.
               ENDIF

               IF( IERR .EQ. 0 ) THEN
                  IERR = -1
                  GOTO 9999
               ENDIF

            ENDIF

 44         IF( VOLUMT(NT) .LT. VOLMOYTE * 1E-3 ) THEN
C              => SUPPRESSION DU TETRAEDRE NT
C              ..............................
               PRINT*,'qtevnul 42: ds1tet NT=',NT,
     %            ' TETRAEDRE de VOLUME',VOLUMT(NT),' est DETRUIT'
               CALL DS1TET( NT,     NBSOTE, N1TEVI, NSTETR,
     %                      NO1TSO, N1TESO, NOTESO )
               QUALIT( NT ) = GRAND
               VOLUMT( NT ) = 0
C              SUPPRESSION DU NO DE MATERIAU de NT
               IF( NBDM .GT. 0 ) NUDMEF( NT ) = 0
               IERR = -1
               GOTO 9999
            ENDIF

         ELSE IF( NFR .EQ. 2 ) THEN

C          NT DE VOLUME QUASI NUL ET 2 TETRAEDRES OPPOSES
C          ==============================================

C           SI UN SOMMET OPPOSE DOUBLE POUR LES 2 TETRAEDRES OPPOSES ALORS
C           GENERER 2 TETRAEDRES A PARTIR DES 3 TETRAEDRES
C           ..............................................................
            CALL Q3T2S2( NT,     NTOP,   NSOP,   GRAND,
     %                   XYZSOM, NBSOTE, N1TEVI, NSTETR, NBDM, NUDMEF,
     %                   NO1TSO, N1TESO, NOTESO,
     %                   MXTETA, VOLUMT, QUALIT, IERR )
            IF( IERR .EQ. -1 ) THEN
               WRITE(IMPRIM,*)'qtevnul 19: q3t2s2 REUSSI NT=',NT,
     %                        ' 3 TETRA OPPOSES NTOP=',NTOP
               GOTO 9999
            ENDIF

C           ESSAYER LE BARYCENTRE DES TETRAEDRES D'ARETE CELLE
C           COMMUNE AUX 2 TETRAEDRES OPPOSES A NT
C           SI MEILLEURE QUALITE ALORS
C           TETRAEDRISER L'ETOILE DES TETRAEDRES DE CETTE ARETE
C           ...................................................
            CALL QTE2TOP( NBSOMM, MXSOMM, XYZSOM, NBDM, NUDMEF,
     %                    NT,     NTOP,   NBSOTE, N1TEVI, NSTETR,
     %                    NO1TSO, N1TESO, NOTESO,
     %                    N1FEOC, N1FEVI, NFETOI, VOETOI, QUETOI,
     %                    MXTETA, VOLUMT, QUALIT, NOTEDS,
     %                    VOLET0, QUAET0, VOLET1, QUAET1, IERR  )
            IF( IERR .EQ. -1 ) THEN
               WRITE(IMPRIM,*) 'qtevnul 20: qte2top REUSSI AVEC NT=',NT,
     %                         ' 2 TETRA OPPOSES NTOP=',NTOP
               GOTO 9999
            ENDIF

C           ESSAI DE SUPPRESSION DE L'ARETE COMMUNE
C           .......................................
            DO I1=1,4
               IF( NTOP(I1) .NE. 0 ) GOTO 45
            ENDDO
 45         DO I2=I1+1,4
               IF( NTOP(I2) .NE. 0 ) GOTO 50
            ENDDO
C           NUMERO DE L'ARETE DE NT COMMUNE AUX 2 TETRAEDRES OPPOSES
 50         CALL NAC2FATE( I1, I2, NAR )
C           ESSAI DE SUPPRESSION DE CETTE ARETE
            CALL QUMTNT( NT,     NAR,
     %                   XYZSOM, NBSOTE, N1TEVI, NSTETR, NBDM, NUDMEF,
     %                   NO1TSO, N1TESO, NOTESO,
     %                   NBTEDS, NOTEDS, QUETOI,
     %                   MXTETA, VOLUMT, QUALIT, IERR )
            IF( IERR .EQ. -1 ) THEN
C              REUSSITE
               WRITE(IMPRIM,*)'qtevnul 21: qumtnt TRIANGLE 2 TETRA OPPOS
     %ES NT=',NT,' NAR=',NAR
               GOTO 9999
            ENDIF

C           CONSTRUIRE L'ETOILE DES TETRAEDRES DE CHACUNE DES 6 ARETES
C           DU TETRAEDRE NT SI ELLE EST INTERNE AU VOLUME
C           CHERCHER LE MEILLEUR BARYCENTRE DE CES TETRAEDRES
C           CREER CETTE TETRAEDRISATION SI VOLUME RESPECTE ET
C           QUALITE AMELIOREE DE L'ETOILE
C           ..........................................................
            CALL QTE6AR( NBSOMM, MXSOMM, XYZSOM, NPSOFR, NBDM, NUDMEF,
     %                   NT,     NBSOTE, N1TEVI, NSTETR,
     %                   NO1TSO, N1TESO, NOTESO,
     %                   N1FEOC, N1FEVI, NFETOI, VOETOI, QUETOI,
     %                   MXTETA, VOLUMT, QUALIT, NOTEDS,
     %                   VOLET0, QUAET0, VOLET1, QUAET1, IERR  )
            IF( IERR .EQ. -1 ) THEN
               WRITE(IMPRIM,*) 'qtevnul 22: qte6ar REUSSI AVEC NT=',NT,
     %                         ' 2 TETRA OPPOSES NTOP=',NTOP
               GOTO 9999
            ENDIF


         ELSE IF( NFR .EQ. 3 ) THEN
C
C           NT DE VOLUME QUASI NUL ET 3 TETRAEDRES OPPOSES
C           ==============================================

C           UN GRAND TETRAEDRE CONTIENT IL NT ET SES 3 TETRAEDRES OPPOSES?
C           LES 3 SOMMETS OPPOSES SONT ILS IDENTIQUES?
            K = 0
 52         K = K + 1
            IF( NSOP(K) .EQ. 0 ) GOTO 52
C           1-ER SOMMET OPPOSE NON NUL
            NS4 = NSOP(K)
            NB  = 0
            DO I=1,4
               IF( NSOP(I) .EQ. NS4 ) NB = NB + 1
            ENDDO

            IF( NB .EQ. 3 ) THEN

C              UN GRAND TETRAEDRE CONTIENT NT ET SES 3 TETRAEDRES OPPOSES
C              => LES 4 TETRAEDRES INTERNES SONT DETRUITS
C                 LE GRAND TETRAEDRE EST CREE
C              ...........................................................
               DO K=1,4
                  IF( NTOP(K) .EQ. 0 ) GOTO 53
               ENDDO

C              LES SOMMETS DE LA FACE FRONTIERE DE NT
 53            NS1T = NSTETR( NOSOFATE(1,K), NT )
               NS2T = NSTETR( NOSOFATE(2,K), NT )
               NS3T = NSTETR( NOSOFATE(3,K), NT )

C              LES 1+3 TETRAEDRES SONT ILS DANS LE MEME MATERIAU?
               NUMATE = 0
               IF( NBDM .GT. 0 ) THEN
                  NUMATE = NUDMEF(NT)
                  DO I=1,4
                     NTO = NTOP(I)
                     IF( NTO .GT. 0 ) THEN
                        IF( NUMATE .NE. NUDMEF( NTO ) ) THEN
                           IERR = 0
                           GOTO 9999
                        ENDIF
                     ENDIF
                  ENDDO
               ENDIF

               DO I=1,4
                  NTO = NTOP(I)
                  IF( NTO .GT. 0 ) THEN
C                    DESTRUCTION DU TETRAEDRE NTO
                     CALL DS1TET( NTO,    NBSOTE, N1TEVI, NSTETR,
     %                            NO1TSO, N1TESO, NOTESO )
                     QUALIT( NTO ) = GRAND
                     VOLUMT( NTO ) = 0
                     IF( NBDM .GT. 0 ) NUDMEF( NTO ) = 0
                  ENDIF
               ENDDO

C              DESTRUCTION DU TETRAEDRE NT
               CALL DS1TET( NT,     NBSOTE, N1TEVI, NSTETR,
     %                      NO1TSO, N1TESO, NOTESO )
               QUALIT( NT ) = GRAND
               VOLUMT( NT ) = 0
               IF( NBDM .GT. 0 ) NUDMEF( NT ) = 0

C              LE NOUVEAU TETRAEDRE NS 1234
               CALL Q04S1T( NS1T, NS2T, NS3T, NS4, NUMATE, NUDMEF,
     %                      XYZSOM, NBSOTE, N1TEVI, NSTETR,
     %                      NO1TSO, N1TESO, NOTESO,
     %                      MXTETA, VOLUMT, QUALIT, IERR )

C              CONSTRUCTION CORRECTE
               IF( IERR .EQ. 0 ) THEN
                  IERR = -1
        WRITE(IMPRIM,*)'qtevnul 23: q04s1t TETRAEDRE NT=',NT,' SUPPRIME'
     %            ,' TETRAEDRE ',NS1T, NS2T, NS3T, NS4,' AJOUTE'
               ENDIF
               GOTO 9999

            ENDIF

C           SI UN SOMMET OPPOSE DOUBLE POUR LES 3 TETRAEDRES OPPOSES ALORS
C           GENERER 2 TETRAEDRES A PARTIR DES 3 TETRAEDRES
C           ..............................................................
            CALL Q3T2S2( NT,     NTOP,   NSOP,   GRAND,
     %                   XYZSOM, NBSOTE, N1TEVI, NSTETR, NBDM, NUDMEF,
     %                   NO1TSO, N1TESO, NOTESO,
     %                   MXTETA, VOLUMT, QUALIT, IERR )
            IF( IERR .EQ. -1 ) THEN
               WRITE(IMPRIM,*)'qtevnul 24: q3t2s2 REUSSI NT=',NT,
     %                        ' 3 TETRA OPPOSES NTOP=',NTOP
               GOTO 9999
            ENDIF

C           CONSTRUIRE L'ETOILE DES TETRAEDRES DE CHACUNE DES 6 ARETES
C           DU TETRAEDRE NT SI ELLE EST INTERNE AU VOLUME
C           CHERCHER LE MEILLEUR BARYCENTRE DE CES TETRAEDRES
C           CREER CETTE TETRAEDRISATION SI VOLUME RESPECTE ET
C           QUALITE AMELIOREE DE L'ETOILE
C           ...........................................................
            CALL QTE6AR( NBSOMM, MXSOMM, XYZSOM, NPSOFR, NBDM, NUDMEF,
     %                   NT,     NBSOTE, N1TEVI, NSTETR,
     %                   NO1TSO, N1TESO, NOTESO,
     %                   N1FEOC, N1FEVI, NFETOI, VOETOI, QUETOI,
     %                   MXTETA, VOLUMT, QUALIT, NOTEDS,
     %                   VOLET0, QUAET0, VOLET1, QUAET1, IERR  )
            IF( IERR .EQ. -1 ) THEN
               WRITE(IMPRIM,*) 'qtevnul 25: qte6ar REUSSI AVEC NT=',NT,
     %                         ' 2 TETRA OPPOSES NTOP=',NTOP
               GOTO 9999
            ENDIF

C           EXISTE T IL UN SOMMET DE NT APPARTENANT AUX 3 TETRAEDRES OPPOSES?
            DO K=1,4

C              SOMMET K DU TETRAEDRE NT
               NS = NSTETR(K,NT)
               NB = 0
               DO 54 L=1,4
                  NT1 = NTOP(L)
                  IF( NT1 .GT. 0 ) THEN
                     DO I=1,4
                        IF( NSTETR(I,NT1) .EQ. NS ) THEN
                           NB = NB + 1
                           GOTO 54
                        ENDIF
                     ENDDO
                  ENDIF
 54            ENDDO

               IF( NB .EQ. 3 ) THEN

C                 LE SOMMET NS VA ETRE DEPLACE AU BARYCENTRE ET 6 VOISINS
C                 DES TETRAEDRES NT ET 3 OPPOSES POUR TESTER
C                 SI LA QUALITE DE CE SOMMET EST AMELIOREE

C                 QUALITE DU SOMMET NS AVANT DEPLACEMENT
                  CALL QUALST( NS, XYZSOM, NBSOTE,NSTETR, NO1TSO,NOTESO,
     %                         VOLUMT, QUALIT,
     %                         VOLET0, QUAET0, NTQMIN, NBTENS )
                  IF( NBTENS .LE. 0 ) GOTO 65

                  IF( NPSOFR(NS) .NE. 0 ) GOTO 65

                  DO M=1,3
                     BAR(M) = 0
                  ENDDO
                  NB = 0

C                 LES 3 TETRAEDRES OPPOSES
                  DO L=1,4
                     NT1 = NTOP(L)
                     IF( NT1 .GT. 0 ) THEN
                        DO I=1,4
                           NST = NSTETR(I,NT1)
                           IF( NST .NE. NS ) THEN
                              NB = NB + 1
                              DO M=1,3
                                 BAR(M) = BAR(M) + XYZSOM(M,NST)
                              ENDDO
                           ENDIF
                        ENDDO
                     ENDIF
                  ENDDO

C                 LE TETRAEDRE NT
                  DO I=1,4
                     NST = NSTETR(I,NT)
                     IF( NST .NE. NS ) THEN
                        NB = NB + 1
                        DO M=1,3
                           BAR(M) = BAR(M) + XYZSOM(M,NST)
                        ENDDO
                     ENDIF
                  ENDDO

                  DO M=1,3

C                    LE BARYCENTRE
                     BAR(M) = BAR(M) / NB

C                    SAUVEGARDE DES XYZ DU SOMMET NS
                     P(M) = XYZSOM(M,NS)
                     PMX(M) = P(M)

C                    DEPLACEMENT DU SOMMET NS AU BARYCENTRE
                     XYZSOM(M,NS) = BAR(M)

                  ENDDO

C                 DIFFERENTES POSITIONS DU SOMMET NS
                  QUAET1MX = QUAET0

                  NBALFA = 0
ccc                  ALFA   = 0.4
                  ALFA   = 0.2


 55               NBALFA = NBALFA + 1
ccc                  ALFA   = ALFA * 4
                  ALFA   = ALFA * 2

                  PB(1) = ( BAR(1) - P(1) ) * ALFA
                  PB(2) = ( BAR(2) - P(2) ) * ALFA
                  PB(3) = ( BAR(3) - P(3) ) * ALFA

                  BARP  = SQRT( PB(1)**2 + PB(2)**2 + PB(3)**2 )

ccc                  print *,'qtevnul: St',NS,' XYZ=',P
ccc                  print *,'qtevnul: St',NS,' BAR=',BAR
ccc                  print *,'qtevnul: St',NS,' PB =',PB
ccc                  print *

                  NBDENS = 0

C                 NOMBRE DE DEPLACEMENT DU SOMMET NS
 56               NBDENS = NBDENS + 1
                  GOTO( 63, 57, 58, 59, 60, 61, 62, 64 ),NBDENS

 57               XYZSOM(1,NS) = BAR(1) + BARP
                  XYZSOM(2,NS) = BAR(2)
                  XYZSOM(3,NS) = BAR(3)
                  GOTO 63

 58               XYZSOM(1,NS) = BAR(1)
                  XYZSOM(2,NS) = BAR(2) + BARP
                  XYZSOM(3,NS) = BAR(3)
                  GOTO 63

 59               XYZSOM(1,NS) = BAR(1)
                  XYZSOM(2,NS) = BAR(2)
                  XYZSOM(3,NS) = BAR(3) + BARP
                  GOTO 63

 60               XYZSOM(1,NS) = BAR(1) - BARP
                  XYZSOM(2,NS) = BAR(2)
                  XYZSOM(3,NS) = BAR(3)
                  GOTO 63

 61               XYZSOM(1,NS) = BAR(1)
                  XYZSOM(2,NS) = BAR(2) - BARP
                  XYZSOM(3,NS) = BAR(3)
                  GOTO 63

 62               XYZSOM(1,NS) = BAR(1)
                  XYZSOM(2,NS) = BAR(2)
                  XYZSOM(3,NS) = BAR(3) - BARP
                  GOTO 63

C                 QUALITE DU SOMMET NS DEPLACE?
 63               CALL QUALST( NS, XYZSOM, NBSOTE,NSTETR, NO1TSO,NOTESO,
     %                         VOLUMT, QUALIT,
     %                         VOLET1, QUAET1, NTQMIN, NBTENS )

                  IF( QUAET1 .GT. QUAET1MX ) THEN
                     QUAET1MX = QUAET1
                     PMX(1) = XYZSOM(1,NS)
                     PMX(2) = XYZSOM(2,NS)
                     PMX(3) = XYZSOM(3,NS)
                  ENDIF

                  GOTO 56

 64               IF( NBALFA .LT. 3 ) GOTO 55

C                 BILAN: UN POINT DE MEILLEURE QUALITE?
                  IF( QUAET1MX .GT. QUAET0 ) THEN

C                    AMELIORATION DE LA QUALITE
C                    NS DEPLACE AU POINT DE QUALITE MAXIMALE
                     XYZSOM(1,NS) = PMX(1)
                     XYZSOM(2,NS) = PMX(2)
                     XYZSOM(3,NS) = PMX(3)
                     CALL QUALST( NS,     XYZSOM, NBSOTE, NSTETR,
     %                            NO1TSO, NOTESO,
     %                            VOLUMT, QUALIT,
     %                            VOLET1, QUAET1, NTQMIN, NBTENS )
                     IERR = -1
                     WRITE(IMPRIM,*)'qtevnul 26: REUSSI DEPLACEMENT du S
     %OMMET',NS,' avec la Qualite=',QUAET1MX
                     GOTO 9999
                  ENDIF

C                 ABANDON DE LA RECHERCHE D'UN MEILLEUR POINT
C                 SOMMET NS REMIS A SA PLACE INITIALE P
                  DO M=1,3
C                    RESTAURATION DES XYZ DU SOMMET NS
                     XYZSOM(M,NS) = P(M)
                  ENDDO

C                 QUALITE DU SOMMET NS AVANT DEPLACEMENT
                  CALL QUALST( NS, XYZSOM, NBSOTE,NSTETR, NO1TSO,NOTESO,
     %                         VOLUMT, QUALIT,
     %                         VOLET0, QUAET0, NTQMIN, NBTENS )

               ENDIF

            ENDDO

C           TENTATIVE DE SUPPRESSION DU SOMMET NS1 OPPOSE A LA PLUS GRANDE FACE
C           DU TETRAEDRE NT POUR AMELIORER LA QUALITE DES TETRAEDRES QUI
C           RECOUVRENT L'ETOILE DES TETRAEDRES AYANT NS1 COMME SOMMET
C           -------------------------------------------------------------------
 65         NOPASS = 0
            NBTEDS = 0
            NUMATE = 0

C           FORMATION DES NBTEDS TETRAEDRES AYANT NS1 COMME SOMMET
            IF( NPSOFR( NS1 ) .NE. 0 ) GOTO 72

C           POSITION DANS NOTESO DU 1-ER TETRAEDRE DE SOMMET NS1
            N = NO1TSO( NS1 )
            VMOYEN = 0D0

C           ETOILE DES TETRAEDRES DE SOMMET NS1 A L'AIDE DE NOTESO
C           TANT QU'IL EXISTE UN TETRAEDRE DE SOMMET NS1 FAIRE
 66         IF( N .GT. 0 ) THEN

C              LE NUMERO DU TETRAEDRE DANS NSTETR
               NT1 = NOTESO(1,N)

C              LE STOCKAGE DU TETRAEDRE DE SOMMET NS1
               NBTEDS = NBTEDS + 1
               NOTEDS( NBTEDS ) = NT1

C              LE TETRAEDRE SUIVANT
               N = NOTESO(2,N)
               GOTO 66

            ENDIF

C           VOLUME ET QUALITES DES NBTEDS TETRAEDRES DE L'ETOILE
 68         VOLET0 = 0
            QUAET0 = GRAND
            VMOYEN = 0D0
            DO N = 1, NBTEDS

C              NUMERO NSTETR DU TETRAEDRE N DE SOMMET NS1
               NT1 = NOTEDS( N )

C              CALCUL DU VOLUME ET QUALITE DU TETRAEDRE NT1
               CALL QUATET( XYZSOM(1,NSTETR(1,NT1)),
     %                      XYZSOM(1,NSTETR(2,NT1)),
     %                      XYZSOM(1,NSTETR(3,NT1)),
     %                      XYZSOM(1,NSTETR(4,NT1)),
     %                      ARMIN, ARMAX, SURFT1, V1, Q1 )
               VOLET0 = VOLET0 + ABS( V1 )
               QUAET0 = MIN( QUAET0, Q1 )

            ENDDO

C           VOLUME MOYEN DES TETRAEDRES DE L'ETOILE
            VMOYEN = VOLET0 / NBTEDS

C           LES NBTEDS TETRAEDRES SONT ILS DANS LE MEME MATERIAU?
            NUMATE = 0
            IF( NBDM .GT. 0 ) THEN
               NUMATE = NUDMEF( NOTEDS( 1 ) )
               DO N = 2, NBTEDS
                  NT1 = NOTEDS( N )
                  IF( NT1 .GT. 0 ) THEN
                     IF( NUMATE .NE. NUDMEF( NT1 ) ) THEN
                        IERR = 0
                        GOTO 72
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF

C           TRACE EVENTUEL DE L'ETOILE DES NBTEDS TETRAEDRES A TETRAEDRISER
            KTITRE='qtevnul: TETRAEDRE           + SUPPRESSION du SOMMET
     %                  ?'
            WRITE(KTITRE(19:27),'(I9)' ) NT
            WRITE(KTITRE(53:61),'(I9)' ) NS1
            CALL SANSDBL( KTITRE, L )
            N = 0
ccc            TRACTE = .TRUE.
            CALL TRNTETRA( N,      KTITRE, XYZSOM, NBSOTE,
     %                     NBTEDS, NOTEDS, NSTETR )

C           CONSTRUCTION DE L'ETOILE DES NBTEDS TETRAEDRES
C           REINITIALISATION A VIDE DES FACES DE NFETOI
            CALL ETOIOCVI( N1FEOC, N1FEVI, NFETOI )

            DO N = 1, NBTEDS

C              NUMERO NSTETR DU TETRAEDRE DE SOMMET NS1
               NT1 = NOTEDS( N )

C              AJOUT DES 4 FACES DU TETRAEDRE NT1 AJOUTE A L'ETOILE
               DO I=1,4
C                 SI    ( LA FACE I DU TETRAEDRE NT1 N'APPARTIENT PAS
C                         AUX FACES DE L'ETOILE NFETOI )
C                 ALORS ELLE EST AJOUTEE    A L'ETOILE NFETOI VERSION 2
C                 SINON ELLE EST SUPPRIMEE DE L'ETOILE NFETOI
                  CALL AJFACE( 1, NT1, I, NBSOTE, NSTETR,
     %                         N1FEOC, N1FEVI, NFETOI,  NF )
                  IF( NF .LT. 0 ) THEN
C                    SATURATION DES FACES DE L'ETOILE
                     GOTO 72
                  ENDIF
               ENDDO

            ENDDO

C           SIMULATION DE LA TETRAEDRISATION DE L'ETOILE DES TETRAEDRES
C           (DE NS1 SANS NS1) ET SANS AJOUT D'UN NOUVEAU POINT INTERNE
C           ...........................................................
            NBSSET0 = 1
            N1SSET( 1 ) = N1FEOC
            CALL SIETSTSU( NS1,    XYZSOM, NBDM,   NUDMEF, NUMATE,
     %                     NBSOTE, N1TEVI, NSTETR, VMOYEN,VOLUMT,QUALIT,
     %                     N1FEVI, NFETOI, MXSSET, NBSSET0, N1SSET,
     %                     NBTETC, MXTEDS-NBTEDS, NOTEDS(NBTEDS+1),
     %                     VOLET1, QUAET1, IERR )

            IF( NBTETC .LE. 0 ) THEN

C              PAS D'AMELIORATION DE LA QUALITE
C              LES NBTETC TETRAEDRES CREES REDEVIENNENT VIDES
               DO K=1,-NBTETC
C                 LE TETRAEDRE NT1 DE NSTETR EST CHAINE VIDE
                  NT1 = NOTEDS( NBTEDS+K )
                  NSTETR( 1, NT1 ) = 0
                  NSTETR( 2, NT1 ) = N1TEVI
                  N1TEVI = NT1
               ENDDO

               IF( NOPASS .EQ. 0 ) THEN

C                TENTATIVE D'AMELIORATION DE L'ETOILE DES TETRAEDRES DE NS1
C                EN L'ENRICHISSANT DES TETRAEDRES DES ARETES SANS SOMMET DE NT
C                ET COMMUNES A 2 FACES DE L'ETOILE FORMANT UN ANGLE > Pi
C                .............................................................

C                RECHERCHE DES ARETES DES FACES DE L'ETOILE TELLES
C                2 SOMMETS NE SOIENT PAS DES SOMMETS DU TETRAEDRE INITIAL NT
C                RAPPEL NFETOI : 1,2,3: NUMERO DES SOMMETS DE LA FACE
C                                DANS LE SENS DIRECT VU DE L'INTERIEUR
                 NF1 = N1FEOC
 680             IF( NF1 .GT. 0 ) THEN

C                  PARCOURS DES 3 ARETES DE LA FACE NF1 DE L'ETOILE
                   DO M1=1,3

                     NSA1 = NFETOI(M1,NF1)
                     IF( NSA1 .NE. NSTETR(1,NT) .AND.
     %                   NSA1 .NE. NSTETR(2,NT) .AND.
     %                   NSA1 .NE. NSTETR(3,NT) .AND.
     %                   NSA1 .NE. NSTETR(4,NT) ) THEN

                     DO 710 M2=M1+1,3

                       NSA2 = NFETOI(M2,NF1)
                       IF( NSA2 .NE. NSTETR(1,NT) .AND.
     %                     NSA2 .NE. NSTETR(2,NT) .AND.
     %                     NSA2 .NE. NSTETR(3,NT) .AND.
     %                     NSA2 .NE. NSTETR(4,NT) ) THEN

C                        NSA1-NSA2 ARETE DE NF1 SANS SOMMET DE NT
C                        RECHERCHE DE LA FACE NF2 COMMUNE DE SOMMETS NSA1 NSA2
                         NF2 = N1FEOC

 690                     IF( NF2 .GT. 0 ) THEN
                           IF( NF2 .NE. NF1 ) THEN
                             DO M3=1,3

                               NSA3 = NFETOI(M3,NF2)
                               IF( NSA3 .EQ. NSA1 ) THEN
                                 DO M4=1,3

                                   NSA4 = NFETOI(M4,NF2)
                                   IF( NSA4 .EQ. NSA2 ) THEN

C                                    L'ARETE M1-M2 DE NF1 EST
C                                    L'ARETE M3-M4 DE NF2
C                                    SIGNE DU VOLUME NF1-NF2 
C                                    POUR DES FACES AVEC NORMALE
C                                    VERS L'INTERIEUR DE L'ETOILE
C                                    RECHERCHE DU 3-EME SOMMET DE NF2
                                     DO M=1,3
                                       IF( M.NE.M3  .AND.
     %                                     M.NE.M4 ) GOTO 700
                                     ENDDO
 700                                 V = VOLTER(
     %                                   XYZSOM(1,NFETOI(1,NF1)),
     %                                   XYZSOM(1,NFETOI(2,NF1)),
     %                                   XYZSOM(1,NFETOI(3,NF1)),
     %                                   XYZSOM(1,NFETOI(M,NF2)) )

                                     IF( V .LT. 0D0 ) THEN
C                                      FACES CONCAVES
C                                      AJOUT A L'ETOILE DES
C                                      TETRAEDRES D'ARETE NSA1 NSA2
                                       CALL AJTE1AR( NSA1,   NSA2,
     %                                               NBSOTE, NSTETR,
     %                                               NO1TSO, NOTESO,
     %                                               MXTEDS, NBTEDS,
     %                                               NOTEDS, IERR )
                                       GOTO 710
                                     ENDIF
                                   ENDIF

                                 ENDDO
                               ENDIF
                             ENDDO

                           ENDIF

C                          PASSAGE A LA FACE NF2 SUIVANTE
                           NF2 = NFETOI(5,NF2)
                           GOTO 690

                         ENDIF
                       ENDIF

C                      FIN SOMMET M2 NSA2 DE NF1
 710                 ENDDO
                     ENDIF

C                    FIN SOMMET M1 NSA1 DE NF1
                   ENDDO

C                  PASSAGE A LA FACE NF1 SUIVANTE
                   NF1 = NFETOI(5,NF1)
                   GOTO 680

                 ENDIF

                 NOPASS = 1
C                FIN TRAITEMENT NOPASS=0 => 1
                 GOTO 68

               ENDIF

C              ABANDON DE CETTE TENTATIVE
               NBTETC = 0
               TRACTE = .FALSE.
               GOTO 72

            ENDIF

C           BILAN: QUALITE MEILLEURE DES NBTETC NOUVEAUX TETRAEDRES?
C           --------------------------------------------------------
            KTITRE='qtevnul: TETRAEDRE           + SUPPRESSION du SOMMET
     %            FAITE'
            WRITE(KTITRE(19:27),'(I9)' ) NT
            WRITE(KTITRE(53:61),'(I9)' ) NS1
            CALL SANSDBL( KTITRE, L )
            N = 1
            CALL TRNTETRA( N, KTITRE, XYZSOM, NBSOTE,
     %                     ABS(NBTETC), NOTEDS(NBTEDS+1), NSTETR )

            IF( IERR .NE. 0 .OR.
     %          ABS(VOLET1-VOLET0) .GT. VOLET0 * 0.001 .OR.
     %          QUAET0 .GE. QUAET1 ) THEN

C              NON: PAS D'AMELIORATION DE LA QUALITE
C              LES NBTETC TETRAEDRES CREES REDEVIENNENT VIDES
               DO K=1,NBTETC
C                 LE TETRAEDRE NT1 DE NSTETR EST CHAINE VIDE
                  NT1 = NOTEDS( NBTEDS+K )
                  NSTETR( 1, NT1 ) = 0
                  NSTETR( 2, NT1 ) = N1TEVI
                  N1TEVI = NT1
               ENDDO

C              PASSAGE AU TRAITEMENT SUIVANT DU TETRAEDRE NT
               GOTO 72

            ENDIF

C           LA NOUVELLE TETRAEDRISATION AMELIORE L'ANCIENNE
C           -----------------------------------------------
C           DESTRUCTION DES NBTEDS NOTEDS() TETRAEDRES DE L'ETOILE INITIALE
            DO N = 1, NBTEDS
C              NUMERO NSTETR DU TETRAEDRE N A DETRUIRE
               NT1 = NOTEDS( N )
C              DETRUIRE LE TETRAEDRE NT DU TABLEAU NSTETR ET NOTESO
               CALL DS1TET( NT1,    NBSOTE, N1TEVI, NSTETR,
     %                      NO1TSO, N1TESO, NOTESO )
            ENDDO

C           MISE A JOUR DES TABLEAUX AUTRES QUE NSTETR
C           DES NBTETC TETRAEDRES QUE sietstsu A CREER
            DO N = 1, NBTETC
C              NUMERO NSTETR DU TETRAEDRE N A DETRUIRE
               NT1 = NOTEDS( NBTEDS+N )
               CALL CHTESO( NT1, NBSOTE, NSTETR,NO1TSO,N1TESO,NOTESO,
     %                      IERR )
C              NUMERO DE MATERIAU
               IF( NBDM .GT. 0 ) NUDMEF( NT1 ) = NUMATE
            ENDDO

            WRITE(IMPRIM,*) 'qtevnul 27: sietstsu REUSSI AVEC NT=',NT,
     %             ' 3 TETRA OPPOSES. SUPPRESSION du SOMMET NS1=',NS1
            IERR = -1
            GOTO 9999


C           LES 4 TETRAEDRES NE FORMENT PAS UN GRAND TETRAEDRE
C           => SUPPRESSION DU TETRAEDRE NT SI LA QUALITE MINIMALE
C              DES 3 TETRAEDRES OPPOSES AUX FACES EST MEILLEURE
C           DEPLACEMENT DU 4-EME SOMMET AU BARYCENTRE DE LA FACE FRONTIERE
C           ..............................................................
 72         IF( NPSOFR( NS1 ) .EQ. 0 ) THEN

C              QUALITE DU SOMMET NS1 SANS DEPLACEMENT
               CALL QUALST( NS1, XYZSOM, NBSOTE, NSTETR, NO1TSO, NOTESO,
     %                      VOLUMT, QUALIT,
     %                      VOLET0, QUAET0, NTQMIN, NBTENS )

C              LE SOMMET NS1 DEVIENT LE BARYCENTRE DE LA FACE OPPOSEE
               L = 0
               DO K=1,3
                  S = 0.0
                  DO I=1,4
                     NS = NSTETR(I,NT)
                     IF( NS .NE. NS1 ) THEN
                        S  = S + XYZSOM(K,NS)
                        IF( NPSOFR(NS) .NE. 0 ) L=L+1
                     ENDIF
                  ENDDO
                  P(K)          = XYZSOM(K,NS1)
                  XYZSOM(K,NS1) = S / 3.0
               ENDDO
               IF( L .EQ. 9 ) NPSOFR(NS1) = 1
C
C              QUALITE DU SOMMET NS1 APRES DEPLACEMENT
               CALL QUALST( NS1, XYZSOM, NBSOTE, NSTETR, NO1TSO, NOTESO,
     %                      VOLUMT, QUALIT,
     %                      VOLET1, QUAET1, NTQMIN, NBTENS )
               IF( NBTENS .LE. 0 ) GOTO 73
C
               IF( QUAET1 .GT. QUAET0 ) THEN

C                 DESTRUCTION DU TETRAEDRE NT
                  CALL DS1TET( NT,     NBSOTE, N1TEVI, NSTETR,
     %                         NO1TSO, N1TESO, NOTESO )
                  QUALIT( NT ) = GRAND
                  VOLUMT( NT ) = 0
                  IF( NBDM .GT. 0 ) NUDMEF( NT ) = 0

                  WRITE(IMPRIM,*)'qtevnul 28: ds1tet REUSSI NT=',NT,
     %          ' TETRAEDRE=TRIANGLE 3 TETRA OPPOSES TUE SOMMET',
     %            NS1,' DEPLACE'
                  IERR = -1
                  GOTO 9999
C
               ENDIF
C
C              REGENERATION DU SOMMET NS1
 73            DO K=1,3
                  XYZSOM(K,NS1) = P(K)
               ENDDO

C              QUALITE DU SOMMET NS1 SANS DEPLACEMENT
               CALL QUALST( NS1, XYZSOM, NBSOTE, NSTETR, NO1TSO, NOTESO,
     %                      VOLUMT, QUALIT,
     %                      VOLET0, QUAET0, NTQMIN, NBTENS )
            ENDIF
C
         ELSE
C
C           NT DE VOLUME QUASI NUL ET 4 TETRAEDRES OPPOSES
C           ==============================================
C           EXISTE T IL 3 SOMMETS OPPOSES IDENTIQUES?
            DO K=1,4
               NS4 = NSOP(K)
               NB = 0
               DO I=1,4
                  IF( NSOP(I) .EQ. NS4 ) NB = NB + 1
               ENDDO
               IF( NB .EQ. 3 ) GOTO 74
            ENDDO

C           PAS 3 SOMMETS OPPOSES IDENTIQUES
            GOTO 80

C           UN GRAND TETRAEDRE CONTIENT NT ET 3 TETRAEDRES OPPOSES
C           => LES 4 TETRAEDRES INTERNES SONT DETRUITS
C              LE GRAND TETRAEDRE EST CREE
C           ......................................................
 74         DO K=1,4
               IF( NSOP(K) .NE. NS4 ) GOTO 75
            ENDDO

C           LES SOMMETS DE LA FACE FRONTIERE DE NT
 75         NS1 = NSTETR( NOSOFATE(1,K), NT )
            NS2 = NSTETR( NOSOFATE(2,K), NT )
            NS3 = NSTETR( NOSOFATE(3,K), NT )

C           LES 4 TETRAEDRES SONT ILS DANS LE MEME MATERIAU?
            NUMATE = 0
            IF( NBDM .GT. 0 ) THEN
               NUMATE = NUDMEF(NT)
               DO I=1,4
                  NTO = NTOP(I)
                  IF( I .NE. K ) THEN
                     IF( NUMATE .NE. NUDMEF( NTO ) ) THEN
                        GOTO 80
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF

            DO I=1,4
               IF( I .NE. K ) THEN
                  NTO = NTOP(I)
C                 DESTRUCTION DU TETRAEDRE NTO
                  CALL DS1TET( NTO,    NBSOTE, N1TEVI, NSTETR,
     %                         NO1TSO, N1TESO, NOTESO )
                  QUALIT( NTO ) = GRAND
                  VOLUMT( NTO ) = 0
                  IF( NBDM .GT. 0 ) NUDMEF( NTO ) = 0
               ENDIF
            ENDDO

C           DESTRUCTION DU TETRAEDRE NT
            CALL DS1TET( NT,  NBSOTE, N1TEVI, NSTETR,
     %                      NO1TSO, N1TESO, NOTESO )
            QUALIT( NT ) = GRAND
            VOLUMT( NT ) = 0
            IF( NBDM .GT. 0 ) NUDMEF( NT ) = 0

C           NT EST LE NOUVEAU TETRAEDRE
            CALL Q04S1T( NS1, NS2, NS3, NS4, NUMATE, NUDMEF,
     %                   XYZSOM, NBSOTE, N1TEVI, NSTETR,
     %                   NO1TSO, N1TESO, NOTESO,
     %                   MXTETA, VOLUMT, QUALIT, IERR )

C           CONSTRUCTION CORRECTE
            IF( IERR .EQ. 0 ) IERR = -1
            WRITE(IMPRIM,*) 'qtevnul 29: q04s1t REUSSI NT=',NT,
     %               ' SOMMETS',NS1, NS2, NS3, NS4

            GOTO 9999


C           LE TETRAEDRE OPPOSE A DECOMPOSER EN 3 TETRAEDRES
C           ................................................
 80         NT1 = NTOP( NFSFMAX )
            IF( NT1 .LE. 0 ) THEN
              WRITE(IMPRIM,*) 'qtevnul 30: ANOMALIE 80 A TRAITER NT=',NT
               IERR = 1
               GOTO 9999
            ENDIF
C
C           LE SOMMET OPPOSE DANS NT1
            NS5 = NSOP( NFSFMAX )
            IF( NS5 .EQ. NS1 .OR. NS5 .EQ. NS2 .OR. NS5 .EQ. NS3 .OR.
     %          NS5 .EQ. NS4 ) THEN
C               LE SOMMET OPPOSE EST UN SOMMET DE NT
               IERR = 0
               GOTO 9999
            ENDIF
C
C           LES 2 TETRAEDRES SONT ILS DANS LE MEME MATERIAU?
            NUMATE = 0
            IF( NBDM .GT. 0 ) THEN
               NUMATE = NUDMEF(NT)
               IF( NUMATE .NE. NUDMEF( NT1 ) ) THEN
                  IERR = 0
                  GOTO 1000
               ENDIF
            ENDIF

C           LA QUALITE DE 2431+2345 EST ELLE > QUALITE 5124+5143+5132?
C           ..........................................................
C           CALCUL DU VOLUME ET QUALITE DU TETRAEDRE DU A LA FACE
            CALL QUATET( XYZSOM(1,NS5),
     %                   XYZSOM(1,NS1),
     %                   XYZSOM(1,NS2),
     %                   XYZSOM(1,NS4),
     %                   ARMIN, ARMAX, SURFT1, V1, Q1 )
            IF( V1 .LE. 0 ) GOTO 9999

            CALL QUATET( XYZSOM(1,NS5),
     %                   XYZSOM(1,NS1),
     %                   XYZSOM(1,NS4),
     %                   XYZSOM(1,NS3),
     %                   ARMIN, ARMAX, SURFT1, V2, Q2 )
            IF( V2 .LE. 0 ) GOTO 9999

            CALL QUATET( XYZSOM(1,NS5),
     %                   XYZSOM(1,NS1),
     %                   XYZSOM(1,NS3),
     %                   XYZSOM(1,NS2),
     %                   ARMIN, ARMAX, SURFT1, V3, Q3 )
            IF( V3 .LE. 0 ) GOTO 9999

            IF( MIN( QUALIT( NT ), QUALIT( NTOP(NFSFMAX) ) ) .GE.
     %          MIN( Q1, Q2, Q3 ) ) GOTO 9999
C
C           SUPPRESSION DES 2 TETRAEDRES NT:1234 ET NT1:2345
            CALL DS1TET( NT,     NBSOTE, N1TEVI, NSTETR,
     %                   NO1TSO, N1TESO, NOTESO )
            QUALIT( NT ) = GRAND
            VOLUMT( NT ) = 0
            IF( NBDM .GT. 0 ) NUDMEF( NT ) = 0
C
            CALL DS1TET( NT1,    NBSOTE, N1TEVI, NSTETR,
     %                   NO1TSO, N1TESO, NOTESO )
            QUALIT( NT1 ) = GRAND
            VOLUMT( NT1 ) = 0
            IF( NBDM .GT. 0 ) NUDMEF( NT1 ) = 0

C           LES TETRAEDRES 2431+2345 DEVIENNENT 5124+5143+5132
            CALL Q05S3T( NS2, NS3, NS4, NS1, NS5, NUMATE, NUDMEF,
     %                   XYZSOM, NBSOTE, N1TEVI, NSTETR,
     %                   NO1TSO, N1TESO, NOTESO,
     %                   MXTETA, VOLUMT, QUALIT, IERR )

           WRITE(IMPRIM,*) 'qtevnul 31: q05s3t REUSSI NT=',NT,':',
     %                      NS2,NS4,NS3,NS1,' +',NS5
           WRITE(IMPRIM,*) 'qtevnul 31: DEVIENT',NS5,NS1,NS2,NS4,' +',
     %                      NS5,NS1,NS4,NS3,' +',NS5,NS1,NS3,NS2


C           DETECTION DES FACES APPARTENANT A 3 TETRAEDRES (A SUPPRIMER...)
            CALL CRHAFAVO( NBSOTE, MXTETR, NSTETR, NUDTETR,
     %                     MOFACE, MXFACE, LFACES, NBFAPB, IERR)
            IF( NBFAPB .GT. 0 ) THEN
               PRINT *
               PRINT *,'qtevnul 31: q05s3t ATTENTION  1 NT=',NT,' avec',
     %                  NBFAPB,' FACES DE 3 TETRAEDRES'
               TRACTE = .TRUE.
            ENDIF

            IF( IERR .EQ. 0 ) THEN
               IERR = -1
               GOTO 9999
            ENDIF

      ENDIF


C     EN DESESPOIR LE TETRAEDRE NT DE MAUVAISE QUALITE EST SUPPRIME
C     =============================================================
 1000 IF( QUALIT(NT) .GE. 0.1 ) GOTO 9999

 1010 KTITRE='qtevnul: TETRAEDRE            V=                Q=        
     %        est DETRUIT'
      WRITE(KTITRE(20:28),'(I9)' ) NT
      WRITE(KTITRE(33:47),'(G15.6)' ) VOLUMT(NT)
      WRITE(KTITRE(51:65),'(G15.6)' ) QUALIT(NT)
      CALL SANSDBL( KTITRE, L )
      PRINT*,KTITRE(1:L)

      TRACTE = .TRUE.
      CALL TR1TEVV( NOPASS, KTITRE(1:L), XYZSOM, NT, NBSOTE, NSTETR,
     %                 NO1TSO, NOTESO, NBTETRA, NOTETRA, VOLET0 )
      TRACTE = TRACTE0
 
      CALL DS1TET( NT,     NBSOTE, N1TEVI, NSTETR,
     %                NO1TSO, N1TESO, NOTESO )
      WRITE(IMPRIM,*) 'qtevnul 1000: Destruction de NT=',NT,
     %   ' NTOP=',NTOP,' V=',VOLUMT(NT),' Q=',QUALIT(NT)
      IERR = -1
      QUALIT( NT ) = GRAND
      VOLUMT( NT ) = 0
      IF( NBDM .GT. 0 ) NUDMEF( NT ) = 0
C     SUPPRESSION DU NO DE MATERIAU DE NT


C     FIN DU TRAITEMENT DU TETRAEDRE NT
C     =================================
C     LES FACES OCCUPEES DE NFETOI REDEVIENNENT VIDES
C     INITIALISATION A VIDE DES MXFAET FACES DU TABLEAU NFETOI
 9999 CALL ETOIVIDE( MXFAET, N1FEOC, N1FEVI, NFETOI )


cccC     TRACE DE FIN DE TRAITEMENT
ccc      IF( NBTEDS .GT. 0 ) THEN
ccc         KTITRE ='qtevnul: TETRAEDRE           + 4 TETRAEDRES OPPOSES av
ccc     %ec ou sans POINT ETOILANT OK'
ccc         WRITE(KTITRE(19:27),'(I9)' ) NT
ccc         CALL SANSDBL( KTITRE, L )
ccc         NOPASS = 0
ccc         TRACTE = .TRUE.
ccc         CALL TRNTETRA( NOPASS, KTITRE(1:L), XYZSOM, NBSOTE,
ccc     %                  NBTEDS, NOTEDS, NSTETR )
ccc         GOTO 9090
ccc      ENDIF

ccc      IF( NBTETRA .GT. 0 ) THEN
ccc         KTITRE = 'qtevnul: TETRAEDRE           + 4 TETRAEDRES OPPOSES e
ccc     %t OPPOSES apres TRAITEMENT'
ccc         WRITE(KTITRE(19:27),'(I9)' ) NT
ccc         CALL SANSDBL( KTITRE, L )
ccc         NOPASS = 1
ccc         TRACTE = .TRUE.
ccc         CALL TR1TEVV( NOPASS, KTITRE(1:L), XYZSOM, NT, NBSOTE, NSTETR,
ccc     %                 NO1TSO, NOTESO, NBTETRA, NOTETRA, VOLET1 )

ccc         IF( VOLET0 .NE. 0D0 .AND. VOLET1 .NE. 0D0 .AND.
ccc     %      ABS(VOLET1-VOLET0) .GT. VOLET0*0.001 ) THEN
ccc            PRINT *,'qtevnul 32: VOLUMES DIFFERENTS V0=',VOLET0,
ccc     %              ' V1=',VOLET1,
ccc     %              ' GAIN=',ABS((VOLET1-VOLET0)/VOLET0)*100,'%'
ccc      ENDIF

ccc      TRACTE = .FALSE.
ccc      ENDIF


cccC        DETECTION DES FACES APPARTENANT A 3 TETRAEDRES
ccc         CALL CRHAFAVO( NBSOTE, MXTETR, NSTETR, NUDTETR,
ccc     %                  MOFACE, MXFACE, LFACES, NBFAPB, IERR)

ccc         IF( NBFAPB .GT. 0 ) THEN
ccc            PRINT *
ccc            PRINT *,'qtevnul fin: 1 NT=',NT,' avec',
ccc     %               NBFAPB,' FACES DE 3 TETRAEDRES'
ccc            TRACTE = .TRUE.
ccc         ENDIF

      PRINT*,'qtevnul: Sortie TETRAEDRE',NT,': St',(NSTETR(I,NT),I=1,4),
     %       ' V=',VOLUMT(NT),' Q=',QUALIT(NT),' --------------------->'
      TRACTE = TRACTE0

      RETURN
      END
