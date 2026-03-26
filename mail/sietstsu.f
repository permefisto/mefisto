      SUBROUTINE SIETSTSU( NST1,   XYZSOM, NBDM,   NUDMEF, NUMATE,
     %                     NBSOTE, N1TEVI, NSTETR, VMOYEN,VOLUMT,QUALIT,
     %                     N1FEVI, NFETOI, MXSSET, NBSSET0, N1SSET,
     %                     NBTETC, MXTETC, NOTETC, VOLET1,QUAET1, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : SIMULATION D'UNE TETRAEDRISATION DES NBSSET0 SOUS-ETOILES
C ----- DE L'ETOILE DES FACES DEFINIES DANS N1SSET et NFETOI
C       SANS POSSIBILITE DE CREER UN POINT INTERNE ETOILANT LES FACES

C       LE NOMBRE DE SOUS-ETOILES EVOLUE AU COURS DE LA TETREDRISATION
C       DANS CHAQUE SOUS-ETOILE LES FACES SONT A NORMALE VERS L'INTERIEUR

C ENTREES:
C --------
C NST1   : NUMERO PTXYZD DU SOMMET EN RECHERCHE D'ELIMINATION DANS L'ETOILE
C XYZSOM : TABLEAU DES COORDONNEES X Y Z DES SOMMETS
C NBDM   : 0 SI 1 MATERIAU=VOLUME, SINON NOMBRE DE MATERIAUX DU VOLUME
C NBSOTE : VALEUR MAXIMALE DE DECLARATION DU PREMIER INDICE DE NSTETR(>3)
C VMOYEN : VOLUME MOYEN DES TETRAEDRES DE L'ETOILE
C          UTILE POUR ELIMINER LES TETRAEDRES DE VOLUME PROCHE DE ZERO
C MXTETC : NOMBRE MAXIMAL DE TETRAEDRES CREABLES DANS NOTETC

C MODIFIES:
C ---------
C NUDMEF : NUMERO DE MATERIAU DE CHAQUE TETRAEDRE DU MAILLAGE
C          ATTENTION: CE TABLEAU EXISTE SEULEMENT SI NBDM>0
C NUMATE : NUMERO DE MATERIAU DE L'ETOILE A REMPLIR

C N1TEVI : NO DU PREMIER TETRAEDRE VIDE DANS NSTETR
C NSTETR : LISTE DU NUMERO DES 4 SOMMETS DES TETRAEDRES
C          NSTETR(1,NT)=0 => TETRAEDRE VIDE ET 
C                            CHAINAGE SUIVANT PAR NSTETR(2,NT)
C VOLUMT : VOLUME  DES TETRAEDRES
C QUALIT : QUALITE DES TETRAEDRES

C N1FEVI : NUMERO NFETOI DE LA PREMIERE FACE VIDE
C NFETOI : LES FACES TRIANGULAIRES DE L'ETOILE VERSION 3
C          1: NUMERO XYZSOM DU SOMMET 1 DE LA FACE
C          2: NUMERO XYZSOM DU SOMMET 2 DE LA FACE
C          3: NUMERO XYZSOM DU SOMMET 3 DE LA FACE
C             S1S2xS1S3 EST LA NORMALE DIRIGEE VERS L'INTERIEUR
C             DE L'ETOILE
C          4: NUMERO DU TETRAEDRE DE CETTE FACE
C          5: CHAINAGE SUIVANT DES FACES OCCUPEES ET VIDES

C MXSSET : NOMBRE MAXIMAL DE SOUS ETOILES DECLARABLES DANS N1SSET
C NBSSET0: NOMBRE INITIAL DE SOUS-ETOILES STOCKEES DANS N1SSET-NFETOI
C N1SSET : NUMERO NFETOI DE LA 1-ERE FACE DES NBSSET0 SOUS ETOILES

C SORTIES:
C --------
C NBTETC : >0 NOMBRE DE TETRAEDRES CREES
C          <0 IMPOSSIBLE DE CREER LES TETRAEDRES SANS AJOUTER UN POINT
C NOTETC : NUMERO NSTETR DES NBTETC TETRAEDRES CREES
C VOLET1 : VOLUME DES NBTETC TETRAEDRES CREES DANS L'ETOILE
C IERR   : 0 SI PAS D'ERREUR TETRAEDRISATION REALISEE
C          1 ETOILE AVEC MOINS DE 4 FACES
C          2 SATURATION D'UN TABLEAU
C          3 3 SOMMETS NON DANS UN TETRAEDRE
C          4 NOMBRE INCORRECT DE FACES DANS L'ETOILE
C          5 SOMMETS DES 4 FACES DE L'ETOILE INCORRECTS
C          6 UNE ARETE DES FACES DE L'ETOILE APPARTIENT A MOINS DE 2
C            OU PLUS DE 2 FACES DE L'ETOILE
C          7 ou 8 ALGORITHME DEFAILLANT A AMELIORER
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE DU PERRAY  Fevrier 2016
C2345X7..............................................................012
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)
      COMMON / TRTETR / STOPTE,TRACTE
      LOGICAL           STOPTE,TRACTE

      REAL              XYZSOM(3,*), VOLUMT(*), QUALIT(*), VOLET1
      REAL              ARMIN, ARMAX, SURFTR(4), V, Q, XYZ(3)

      INTEGER           NSTETR(NBSOTE,*), NUDMEF(NBDM),
     %                  NFETOI(5,*), NOTETC(MXTETC),
     %                  N1SSET(MXSSET)

      DOUBLE PRECISION  CBTR(3), VOLTER, VD, VOLMAX, VMOYEN

      CHARACTER*128     KTITRE

      INTEGER           NOSOTR(4), NOSOTR2(3), NOSOTR3(3)
      EQUIVALENCE      (NOSOTR(1),NS1), (NOSOTR(2),NS2),
     %                 (NOSOTR(3),NS3), (NOSOTR(4),NS4)

      INTEGER           NOSOFA(3,4)
      DATA              NOSOFA / 1,3,2,  2,3,4,  3,1,4,  4,1,2 /
C     => LA NORMALE A LA FACE EST ORIENTEE VERS L'EXTERIEUR DU TETRAEDRE

      PRINT*,'sietstsu: Sommet',NST1,' dans',NBSSET0,
     %       ' SOUS-ETOILES A TETRAEDRISER Volume MOYEN=', VMOYEN

C     NOMBRE DE TRAITEMENTS DE L'ETOILE
      NBPASS = 0
      IERR   = 0

C     NOMBRE DE SOUS-ETOILES DE L'ETOILE
      NBSSET = NBSSET0

C     NOMBRE DE TETRAEDRES CREES DANS L'ETOILE
      NBTETC = 0

C     VOLUME DES NBTETC TETRAEDRES A CREER
      VOLET1 = 0D0
C     QUALITE MINIMALE DES NBTETC TETRAEDRES A CREER
      QUAET1 = 2.0

C     ======================================================================
C     TETRAEDRISATION DE LA PILE DES NBSSET0 SOUS-ETOILES DE L'ETOILE NFETOI
C     ======================================================================
 5    IF( NBSSET .LE. 0 ) GOTO 9000

C     TRACE DE L'ETOILE ACTUELLE NFETOI ET DES NBTETC TETRAEDRES CREES
      KTITRE='          TETRAEDRES CREES et FACES DE L''ETOILE'
      WRITE(KTITRE(1:9),'(I9)') NBTETC
      CALL SANSDBL( KTITRE, L )
      CALL TRFETOV3( KTITRE(1:L), XYZSOM, NBSSET, N1SSET, NFETOI,
     %               NBTETC, NOTETC, NBSOTE, NSTETR )

C     TETRAEDRISATION DE LA SOUS-ETOILE NBSSET EN HAUT DE PILE N1SSET
C     ===============================================================
C     1-ERE FACE NFETOI VERSION 3 DE LA SOUS-ETOILE NBSSET
      N1FEOC = N1SSET( NBSSET )

C     LA SOUS-ETOILE EST DEPILEE CAR EN TETRAEDRISATION
      NBSSET = NBSSET - 1

      IF( N1FEOC .LE. 0 ) GOTO 5

C     NOMBRE DE TETRAEDRES CREES ACTUELLEMENT DANS L'ETOILE
      NBTETC0 = NBTETC
      NBPASS  = NBPASS + 1
      NOAPPEL = 0

      IF( NBTETC .GT. MXTETC ) THEN
         print*,'?????????????????????????????????????????????????????'
         print*,'PB sietstsu: TROP DE TETRAEDRES CREES  ',NBTETC
         print*,'PB sietstsu: BOUCLE INFINIE?'
         print*,'?????????????????????????????????????????????????????'
         IERR = 8
         GOTO 9020
      ENDIF

C     NBFETO NOMBRE DE FACES DE LA SOUS-ETOILE NBSSET
C     -----------------------------------------------
      NBFETO = 0
      NF0    = N1FEOC

 10   IF( NF0 .GT. 0 ) THEN
C        UNE FACE DE PLUS
         NBFETO = NBFETO + 1
C        NO DU 1-ER SOMMET DE LA FACE EST RENDU POSITIF
         NFETOI(1,NF0) = ABS( NFETOI(1,NF0) ) 
C        LA FACE SUIVANTE DE NFETOI
         NF0 = NFETOI(5,NF0)
         GOTO 10
      ENDIF

ccc      PRINT *,'sietstsu: SOUS-ETOILE',NBSSET,' N1FEOC=',N1FEOC,
ccc     %        ' NB FACES de la SOUS-ETOILE=',NBFETO


C     MODE DE CREATION DES TETRAEDRES SELON LE NOMBRE DE FACES
C     --------------------------------------------------------
      IF( NBFETO .EQ. 0 ) THEN
C        PLUS AUCUNE FACE DANS LA SOUS-ETOILE NBSSET+1
C        PASSAGE A LA SOUS-ETOILE EN HAUT DE PILE
         GOTO 5
      ENDIF

      IF( NBFETO .LE. 3 ) THEN
C        MOINS DE 4 FACES DANS LA SOUS-ETOILE => INCORRECT
         PRINT*,'PB sietstsu: NOMBRE DE FACES=',NBFETO,' INCORRECT'
         IERR = 4
         KTITRE='PB sietstsu: NBFETO=      '
         WRITE(KTITRE(21:22),'(I2)') NBFETO
         GOTO 9010
      ENDIF

      IF( NBFETO .EQ. 4 ) THEN

C        ---------------------------------------------------------
C        0) LES 4 SEULES FACES DE LA SOUS-ETOILE DOIVENT FORMER LE
C           DERNIER TETRAEDRE DE LA SOUS-ETOILE SINON ERREUR
C        ---------------------------------------------------------
C        LES 3-ERS SOMMETS DU TETRAEDRE = CEUX DE LA FACE N1FEOC
         NS1 = NFETOI( 1, N1FEOC )
         NS2 = NFETOI( 2, N1FEOC )
         NS3 = NFETOI( 3, N1FEOC )
C        LE 4-EME SOMMET EST CELUI DES FACES DIFFERENT DES 3 DE LA FACE 1
         NF0 = NFETOI( 5, N1FEOC )
         DO K=1,3
            NS4 = NFETOI( K, NF0 )
            IF( NS4 .NE. NS1 .AND. NS4 .NE. NS2 .AND. NS4 .NE. NS3 )
     %          GOTO 12
         ENDDO

C        VERIFICATION: LES 3 DERNIERES FACES DOIVENT AVOIR LEURS
C                      3 SOMMETS PARMI LES 4 DU TETRAEDRE
 12      DO 20 K=2,4

C           LES 3 SOMMETS DE LA FACE K DU TETRAEDRE
            DO L=1,3
               NOSOTR2(L) = NOSOTR( NOSOFA(L,K)  )
            ENDDO
C           TRI CROISSANT DES 3 NUMEROS DE SOMMETS DE LA FACE K
            CALL TRI3NO( NOSOTR2, NOSOTR2 )

            NF = NF0
C           LES 3 SOMMETS DE LA FACE NF DE LA SOUS-ETOILE
 16         DO L=1,3
               NOSOTR3(L) = NFETOI(L,NF)
            ENDDO
C           TRI CROISSANT DES 3 NUMEROS DE SOMMETS DE LA FACE NF
            CALL TRI3NO( NOSOTR3, NOSOTR3 )
            IF( NOSOTR2(1) .EQ. NOSOTR3(1) .AND.
     %          NOSOTR2(2) .EQ. NOSOTR3(2) .AND.
     %          NOSOTR2(3) .EQ. NOSOTR3(3) ) THEN
C               LA FACE K DU TETRAEDRE EST LA FACE NF DE LA SOUS-ETOILE
                GOTO 20
            ENDIF
C           FACE SUIVANTE
            NF = NFETOI(5,NF)
            IF( NF .NE. 0 ) GOTO 16

            PRINT*,'PB sietstsu: FACE',K,'DU TETRAEDRE',NOSOTR,
     %   ' N''EST PAS UNE DES 4 FACES DE L''ETOILE => ETOILE INCORRECTE'
            IERR  = 5
            NOCAS = 0
            KTITRE='PB sietstsu: FACE NON TROUVEE DANS LE TETRAEDRE    '
            WRITE(KTITRE(49:57),'(I9)') NOSOTR(1)
            WRITE(KTITRE(59:67),'(I9)') NOSOTR(2)
            WRITE(KTITRE(69:77),'(I9)') NOSOTR(3)
            WRITE(KTITRE(79:87),'(I9)') NOSOTR(4)
            GOTO 9010

 20      ENDDO

C        AJOUT DU TETRAEDRE FORME DES 4 SEULES FACES DE LA SOUS-ETOILE
         NF = N1FEOC
         CALL AJ1TEET( XYZSOM, NBSOTE, N1TEVI, NSTETR,
     %                 NBDM,   NUDMEF, NUMATE,
     %                 NS4,    NF,     N1FEOC, N1FEVI, NFETOI,
     %                 VOLUMT, QUALIT,
     %                 MXTETC, NBTETC, NOTETC, VOLET1, QUAET1, IERR )

C        LE TETRAEDRE FACE N1FEOC + SOMMET NS4 NBTETC EST CREE
         NOCAS = 1
         NTE   = NOTETC( NBTETC )
        print*,'sietstsu NOCAS=',NOCAS,': 4Faces =>1 TETRAEDRE',NTE,'=',
     % (NSTETR(KK,NTE),KK=1,4),' N1TEVI=',N1TEVI,' QualTet=',QUALIT(NTE)

C        FIN TETRAEDRISATION
         GOTO 8000

      ENDIF

C     ----------------------------------------------------------------------------
C     1) ESSAI DE TETRAEDRISER LES TRIANGLES DE LA SOUS-ETOILE SANS AJOUT DE POINT
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
            NS1F0 = NFETOI(K ,NF0)
            NS2F0 = NFETOI(KK,NF0)
C           LE SOMMET 3 DE LA FACE NF0
            IF( KK .EQ. 3 ) THEN
               KKK = 1
            ELSE
               KKK = KK + 1
            ENDIF
            NS3F0 = NFETOI(KKK,NF0)

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
                     NS1F1 = NFETOI(L ,NF1)
                     NS2F1 = NFETOI(LL,NF1)
                     IF( NS1F0 .EQ. NS2F1 .AND. NS2F0 .EQ. NS1F1 ) THEN
ccc     %          .OR. NS1F0 .EQ. NS1F1 .AND. NS2F0 .EQ. NS2F1 ) THEN
C                       L'ARETE K DE NF0 EST L'ARETE L DE NF1
C                       LE SOMMET 3 DE LA FACE NF1
                        IF( LL .EQ. 3 ) THEN
                           LLL = 1
                        ELSE
                           LLL = LL + 1
                        ENDIF
                        NS3F1 = NFETOI(LLL,NF1)

                        IF( NS3F1 .EQ. NS3F0 ) THEN
C                          CAS DE 2 FACES DU CF IDENTIQUES DANS L'ENCOCHE
                       print*,'PB sietstsu: 2FACES avec 3 MEMES SOMMETS'
     %                        ,' NF0=',NF0,' NF1=',NF1,
     %                         ' ST:',NS1F0,NS2F0,NS3F0
                           GOTO 69
                        ENDIF

C                       LE VOLUME FACE NF0-SOMMET NS3F1 EST IL POSITIF?
                        VD = VOLTER( XYZSOM(1,NFETOI(1,NF0)),
     %                               XYZSOM(1,NFETOI(2,NF0)),
     %                               XYZSOM(1,NFETOI(3,NF0)),
     %                               XYZSOM(1,NS3F1) )
                        IF( VD .LE. 0D0 ) THEN
C                          TETRAEDRE DE VOLUME<=0
C                           => RECHERCHE D'UNE AUTRE FACE NF1
                           GOTO 69
                        ENDIF

C                       RECHERCHE DE LA 3-eme FACE NON NF0 et NON NF1
                        NF2 = N1FEOC
 53                     IF( NF2 .GT. 0 ) THEN
                           IF( NF2 .NE. NF0 .AND. NF2 .NE. NF1 ) THEN

                              DO 65 M=1,3

C                                LES 2 SOMMETS DE L'ARETE M DE LA FACE NF2
                                 NS1F2 = NFETOI(M ,NF2)
                                 IF( NS1F2 .NE. NS3F1 ) GOTO 65

C                                SOMMETS NS1F2=NS3F1
                                 IF( M .EQ. 3 ) THEN
                                    MM = 1
                                 ELSE
                                    MM = M + 1
                                 ENDIF
                                 NS2F2 = NFETOI(MM,NF2)

                                 IF( MM .EQ. 3 ) THEN
                                    MMM = 1
                                 ELSE
                                    MMM = MM + 1
                                 ENDIF
C                                LE SOMMET 3 DE LA FACE NF2
                                 NS3F2 = NFETOI(MMM,NF2)

                                 IF( NS2F2 .EQ. NS1F0 ) THEN
                                    IF( NS3F2 .EQ. NS3F0 ) THEN

C                                      NOCAS 2:  TETRAEDRE NF0-NS3F1 avec
C                                      NS1F2=NS3F1 et NS2F2=NS1F0 et NS3F2=NS3F0
C                                      -----------------------------------------
                                       NOCAS = 2

C                                      UNE DES 3 ARETES DE NF2 INTERSECTE T ELLE
C                                      LA FACE NS3F0 NS2F0 NS3F1?
                                       N0 = NFETOI(3,NF2)
                                       DO N=1,3
                                          N1 = NFETOI(N,NF2)
                                          CALL INARTRR( XYZSOM(1,N0),
     %                                                  XYZSOM(1,N1),
     %                                                  XYZSOM(1,NS3F0),
     %                                                  XYZSOM(1,NS2F0),
     %                                                  XYZSOM(1,NS3F1),
     %                                                  LINTER,XYZ,CBTR)
                                          IF( LINTER  .EQ. 1 .AND.
     %                                     CBTR(1) .LT. 0.999999D0 .AND.
     %                                     CBTR(2) .LT. 0.999999D0 .AND.
     %                                     CBTR(3) .LT. 0.999999D0 )THEN
C                                          OUI: UN POINT D'INTERSECTION
C                                          ABANDON DE CE CHOIX DE 2 FACES
                                    print *,'bizarre nocas=2 cbtr=',cbtr
                                           GOTO 69
                                          ENDIF
                                          N0 = N1
                                       ENDDO

C                                      FACE NS2F0-NS3F0-NS3F1
C                                      DU TETRAEDRE NS1F0+NS2F0+NS3F0+NS3F1
                                       NOSOTR2(1) = NS2F0
                                       NOSOTR2(2) = NS3F0
                                       NOSOTR2(3) = NS3F1
                                       CALL TRI3NO( NOSOTR2, NOSOTR2 )

                                       NF3 = N1FEOC
 54                                   IF( NF3 .GT. 0 ) THEN
                                      CALL TRI3NO(NFETOI(1,NF3),NOSOTR3)

                                      IF(NOSOTR2(1).EQ.NOSOTR3(1) .AND.
     %                                   NOSOTR2(2).EQ.NOSOTR3(2) .AND.
     %                                   NOSOTR2(3).EQ.NOSOTR3(3) )THEN
C                                        FACE 2 DU TETRAEDRE NS1F0-NS2F0-NS3F0-NS3F1
                                         GOTO 60
                                       ENDIF

                                       NF3 = NFETOI(5,NF3)
                                       GOTO 54
                                       ENDIF

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
                                       N0 = NFETOI(3,NF2)
                                       DO N=1,3
                                          N1 = NFETOI(N,NF2)
                                          CALL INARTRR( XYZSOM(1,N0),
     %                                                  XYZSOM(1,N1),
     %                                                  XYZSOM(1,NS3F0),
     %                                                  XYZSOM(1,NS3F1),
     %                                                  XYZSOM(1,NS1F0),
     %                                                  LINTER,XYZ,CBTR)
                                          IF( LINTER  .EQ. 1 .AND.
     %                                      CBTR(1).LT.0.999999D0 .AND.
     %                                      CBTR(2).LT.0.999999D0 .AND.
     %                                      CBTR(3).LT.0.999999D0)THEN
C                                           OUI: UN POINT D'INTERSECTION
C                                           ABANDON DE CE CHOIX DE 2 FACES
                                    print *,'bizarre nocas=3 cbtr=',cbtr
                                            GOTO 69
                                          ENDIF
                                          N0 = N1
                                       ENDDO

C                                      FACE NS1F0-NS3F0-NS3F1
C                                      DU TETRAEDRE NS1F0+NS2F0+NS3F0+NS3F1
                                       NOSOTR2(1) = NS1F0
                                       NOSOTR2(2) = NS3F0
                                       NOSOTR2(3) = NS3F1
                                       CALL TRI3NO( NOSOTR2, NOSOTR2 )

                                       NF3 = N1FEOC
 57                                   IF( NF3 .GT. 0 ) THEN
                                      CALL TRI3NO(NFETOI(1,NF3),NOSOTR3)
                                      IF(NOSOTR2(1) .EQ. NOSOTR3(1).AND.
     %                                   NOSOTR2(2) .EQ. NOSOTR3(2).AND.
     %                                   NOSOTR2(3) .EQ. NOSOTR3(3))THEN
                                         GOTO 60
                                       ENDIF
                                       NF3 = NFETOI(5,NF3)
                                       GOTO 57
                                       ENDIF

                                       GOTO 60
                                    ENDIF
                                 ENDIF
                                 GOTO 65

C                                LES 3 FACES NF0 NF1 NF2 FORMENT UN TETRAEDRE
C                                NS1F0-NS2F0-NS3F0-NS3F1 A AJOUTER
C                                CREATION DU TETRAEDRE NS1F0-NS2F0-NS3F0-NS3F1
 60                              NFETOI(1,NF0) = NS1F0
                                 NFETOI(2,NF0) = NS2F0
                                 NFETOI(3,NF0) = NS3F0

C                                AJOUT DU TETRAEDRE FORME DE LA FACE NF0
C                                                         + SOMMET NS3F1
                                 CALL AJ1TEET( XYZSOM, NBSOTE, N1TEVI,
     %                                         NSTETR,
     %                                         NBDM,   NUDMEF, NUMATE,
     %                                         NS3F1,  NF0,
     %                                         N1FEOC, N1FEVI, NFETOI,
     %                                         VOLUMT, QUALIT,
     %                                         MXTETC, NBTETC, NOTETC,
     %                                         VOLET1, QUAET1, IERR )
                                 IF( IERR .NE. 0 ) GOTO 9020

C                                TETRAEDRE FACE NF0+SOMMET NS3F1 NBTETC EST CREE
                                 NTE = NOTETC( NBTETC )
                                 print *,'sietstsu NOCAS=',NOCAS,
     %                           ': 3Faces=> TETRAEDRE+',NTE,'=',
     %                           (NSTETR(KK,NTE),KK=1,4),' N1TEVI=',
     %                            N1TEVI,' QualTet=',QUAET1

C                                NOUVELLE ETOILE APRES AJOUT DU TETRAEDRE NTE
                                 GOTO 8000

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

C        LA FACE SUIVANTE DE NFETOI
         NF0 = NFETOI(5,NF0)
         GOTO 51

      ENDIF

C     -------------------------------------------------------------------------
C     2) ESSAI DE TETRAEDRISER LES TRIANGLES DE LA SOUS-ETOILE SANS AJOUT DE POINT
C        RECHERCHE DANS LA SOUS-ETOILE DE 2 FACES ADJACENTES FORMANT UN TETRAEDRE
C        DE VOLUME>0 ET DONT L'ARETE OPPOSEE INTERSECTE AUCUNE FACE DE LA SOUS-ETOILE
C     -------------------------------------------------------------------------
      NF0 = N1FEOC
 151  IF( NF0 .GT. 0 ) THEN

C        BOUCLE SUR LES 3 ARETES DE LA FACE NF0
         DO K=1,3

C           L'ARETE K DE LA FACE NF0
            IF( K .EQ. 3 ) THEN
               KK = 1
            ELSE
               KK = K + 1
            ENDIF
C           LES 2 SOMMETS DE L'ARETE K DE LA FACE NF0
            NS1F0 = NFETOI(K ,NF0)
            NS2F0 = NFETOI(KK,NF0)

C           LE SOMMET 3 DE LA FACE NF0
            IF( KK .EQ. 3 ) THEN
               KKK = 1
            ELSE
               KKK = KK + 1
            ENDIF
            NS3F0 = NFETOI(KKK,NF0)

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
                     NS1F1 = NFETOI(L ,NF1)
                     NS2F1 = NFETOI(LL,NF1)
                     IF( NS1F0 .EQ. NS2F1 .AND. NS2F0 .EQ. NS1F1 ) THEN
c     %          .OR. NS1F0 .EQ. NS1F1 .AND. NS2F0 .EQ. NS2F1 ) THEN
C                       L'ARETE K DE NF0 EST L'ARETE L DE NF1
C                       LE SOMMET 3 DE LA FACE NF1
                        IF( LL .EQ. 3 ) THEN
                           LLL = 1
                        ELSE
                           LLL = LL + 1
                        ENDIF
                        NS3F1 = NFETOI(LLL,NF1)

                        IF( NS3F1 .EQ. NS3F0 ) THEN
C                          CAS DE 2 FACES DU CF IDENTIQUES DANS L'ENCOCHE
                       print*,'Pb sietstsu: 2FACES avec 3 MEMES SOMMETS'
     %                        ,' NF0=',NF0,' NF1=',NF1,
     %                         ' ST:',NS1F0,NS2F0,NS3F0
                          GOTO 159
                        ENDIF

C                       NS1F0-NS2F0-NS3F0-NS3F1 TETRAEDRE DE VOLUME>0?
                        VD = VOLTER( XYZSOM(1,NS1F0),
     %                               XYZSOM(1,NS2F0),
     %                               XYZSOM(1,NS3F0),
     %                               XYZSOM(1,NS3F1) )
                        IF( VD .LE. VMOYEN*1D-3 ) THEN
C                          TETRAEDRE DE VOLUME<=0 OU PROCHE 0
C                          => RECHERCHE D'UNE AUTRE FACE NF1
                           GOTO 159
                        ENDIF

C                       L'ARETE NS3F0-NS3F1 INTERSECTE T ELLE UNE FACE DE LA SOUS-ETOILE?
                        NF2 = N1FEOC
 153                    IF( NF2 .GT. 0 ) THEN
                           IF( NF2 .NE. NF0 .AND. NF2 .NE. NF1 ) THEN

C                             CETTE FACE NF2 EST ELLE INTERSECTEE PAR 
C                             L'ARETE NS3F0-NS3F1?
                              CALL INARTRR( XYZSOM(1,NS3F0),
     %                                     XYZSOM(1,NS3F1),
     %                                     XYZSOM(1,NFETOI(1,NF2)),
     %                                     XYZSOM(1,NFETOI(2,NF2)),
     %                                     XYZSOM(1,NFETOI(3,NF2)),
     %                                     LINTER, XYZ, CBTR )
C                             LINTER : -2 SI S1=S2
C                             -1 SI S1-S2 PARALLELE AU PLAN DU TRIANGLE
C                              0 SI S1-S2 N'INTERSECTE PAS LE TRIANGLE
C                              1 SI S1-S2   INTERSECTE     LE TRIANGLE ET ENTRE S1-S2
C                            XYZ : 3 COORDONNEES DU POINT D'INTERSECTION SI LINTER=1
C                            CBTR: 3 COORDONNEES BARYCENTRIQUES DE XYZ DANS LE TRIANGLE
                              IF(  LINTER  .EQ. 1 .AND.
     %                             CBTR(1) .LT. 0.999999D0 .AND.
     %                             CBTR(2) .LT. 0.999999D0 .AND.
     %                             CBTR(3) .LT. 0.999999D0 ) THEN
C                                OUI: UN POINT D'INTERSECTION
C                                ABANDON DE CE CHOIX DE 2 FACES
                                 GOTO 159
                              ENDIF

C                             UNE DES 3 ARETES DE NF2 INTERSECTE T ELLE
C                             LA FACE NS3F0 NS2F0 NS3F1?
                              N0 = NFETOI(3,NF2)
                              DO M=1,3
                                 N1 = NFETOI(M,NF2)
                                 CALL INARTRR( XYZSOM(1,N0),
     %                                        XYZSOM(1,N1),
     %                                        XYZSOM(1,NS3F0),
     %                                        XYZSOM(1,NS2F0),
     %                                        XYZSOM(1,NS3F1),
     %                                        LINTER, XYZ, CBTR )
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

C                             UNE DES 3 ARETES DE NF2 INTERSECTE T ELLE
C                             LA FACE NS3F0 NS3F1 NS1F0?
                              N0 = NFETOI(3,NF2)
                              DO M=1,3
                                 N1 = NFETOI(M,NF2)
                                 CALL INARTRR( XYZSOM(1,N0),
     %                                        XYZSOM(1,N1),
     %                                        XYZSOM(1,NS3F0),
     %                                        XYZSOM(1,NS3F1),
     %                                        XYZSOM(1,NS1F0),
     %                                        LINTER, XYZ, CBTR )
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

C                       CREATION DU TETRAEDRE NS1F0-NS2F0-NS3F0-NS3F1
                        NFETOI(1,NF0) = NS1F0
                        NFETOI(2,NF0) = NS2F0
                        NFETOI(3,NF0) = NS3F0

C                       AJOUT DU TETRAEDRE FORME DE LA FACE NF0 + SOMMET NS3F1
                        CALL AJ1TEET( XYZSOM, NBSOTE, N1TEVI, NSTETR,
     %                                NBDM,   NUDMEF, NUMATE,
     %                                NS3F1,  NF0,
     %                                N1FEOC, N1FEVI, NFETOI,
     %                                VOLUMT, QUALIT,
     %                                MXTETC, NBTETC, NOTETC,
     %                                VOLET1, QUAET1, IERR )
                        IF( IERR .NE. 0 ) GOTO 9020

C                       LE TETRAEDRE FACE NF0 + SOMMET NS3F1 NBTETC EST CREE
                        NTE = NOTETC( NBTETC )
                        print *,'sietstsu NOCAS=',NOCAS,
     %                  ': 2Faces=> TETRAEDRE+',NTE,'=',
     %                  (NSTETR(KK,NTE),KK=1,4),' N1TEVI=',N1TEVI,
     %                  ' QualTet=',QUAET1

C                       NOUVELLE ETOILE APRES AJOUT DU TETRAEDRE NTE
                        GOTO 8000

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
         GOTO 151

      ENDIF

      IF( N1FEOC .LE. 0 ) THEN
C        TOUTE LA SOUS-ETOILE A ETE TETRAEDRISEE
         GOTO 8000
      ENDIF


C     -----------------------------------------------------------------------
C     3) RECHERCHE D'UNE FACE FORMANT AVEC UN SOMMET DE LA SOUS-ETOILE PROCHE
C        UN TETRAEDRE DE VOLUME POSITIF ET SANS INTERSECTION AVEC LES FACES
C     -----------------------------------------------------------------------
      VOLMAX = -1D123
      NF0 = N1FEOC
 200  IF( NF0 .GT. 0 ) THEN

C        RECHERCHE D'UN SOMMET NS4 D'UNE FACE NON ADJACENTE ET OPPOSEE DANS LA SOUS-ETOILE
C        ET SANS INTERSECTION AVEC LA SOUS-ETOILE
         NF1 = N1FEOC
 210     IF( NF1 .GT. 0 ) THEN
            IF( NF1 .NE. NF0 ) THEN

               DO 250 L=1,3

C                 LE SOMMET L DE LA FACE NF1
                  NS4 = NFETOI(L,NF1)
                  IF( NS4 .EQ. NFETOI(1,NF0)  .OR.
     %                NS4 .EQ. NFETOI(2,NF0)  .OR.
     %                NS4 .EQ. NFETOI(3,NF0) ) GOTO 250

C                 LE VOLUME DU TETRAEDRE FACE NF0 + NS4 EST-IL POSITIF?
                  VD = VOLTER( XYZSOM(1,NFETOI(1,NF0) ),
     %                         XYZSOM(1,NFETOI(2,NF0) ),
     %                         XYZSOM(1,NFETOI(3,NF0) ),
     %                         XYZSOM(1,NS4) )
                  IF( VD .LE. VMOYEN*1D-3 ) GOTO 250

C                 LE TETRAEDRE NF0-NS4 EST IL SANS INTERSECTION AVEC LA SOUS-ETOILE?
                  NF2 = N1FEOC
 230              IF( NF2 .GT. 0 ) THEN
                     IF( NF2 .NE. NF0 ) THEN

C                       INTERSECTION TRIANGLE NF2-TETRAEDRE NF0+NS4?
                        CALL INTRITETR( XYZSOM(1,NFETOI(1,NF2)),
     %                                  XYZSOM(1,NFETOI(2,NF2)),
     %                                  XYZSOM(1,NFETOI(3,NF2)),
     %                                  XYZSOM(1,NFETOI(1,NF0)),
     %                                  XYZSOM(1,NFETOI(2,NF0)),
     %                                  XYZSOM(1,NFETOI(3,NF0)),
     %                                  XYZSOM(1,NS4),
     %                                  LINTER )
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
                  IF( V .GT. VOLMAX ) THEN
                     VOLMAX = V
                     NF0MAX = NF0
                     NS4MAX = NS4
                     CALL QUATET( XYZSOM(1,NFETOI(1,NF0)),
     %                            XYZSOM(1,NFETOI(2,NF0)),
     %                            XYZSOM(1,NFETOI(3,NF0)),
     %                            XYZSOM(1,NS4),
     %                            ARMIN, ARMAX, SURFTR, V, Q )
                     print*,'sietstsu nocas=5: VOLMAX=',VOLMAX,
     %                      ' qualite du tetra=',Q
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

C     ------------------------------------------------------------
C     4) IL EXISTE UNE CONFIGURATION DE FACES DE LA SOUS-ETOILE
C        NE POUVANT PAS AVEC 3 OU 2 OU 1 FACE FORMER UN TETRAEDRE
C        => UN  NOUVEAU POINT SERAIT A AJOUTER ET A TETRAEDRISER
C           AVEC LES FACES DE LA SOUS-ETOILE
C        MAIS ICI, UN TEL POINT CHERCHE A ETRE DETRUIT => ABANDON
C     ------------------------------------------------------------
 300  GOTO 9990


C     AFFICHAGE ET TRACE DE LA SOUS-ETOILE ACTUELLE
C     ---------------------------------------------
C     NOMBRE DE SES FACES
 8000 NBFETO = 0
      NF0    = N1FEOC
 8010 IF( NF0 .GT. 0 ) THEN
         NBFETO = NBFETO + 1
C        LA FACE SUIVANTE DE NFETOI
         NF0 = NFETOI(5,NF0)
         GOTO 8010
      ENDIF

      IF( N1FEOC .GT. 0 ) THEN
C        TOUTE LA SOUS-ETOILE N'A PAS ENCORE ETE TETRAEDRISEE
C        => UN NOUVEAU TRAITEMENT DE LA SOUS-ETOILE
C        LA TETRAEDRISATION A T ELLE CREEE PLUSIEURS SOUS-SOUS-ETOILES?
         CALL SOUSETOR(XYZSOM, N1FEOC,NFETOI, MXSSET,NBSSET,N1SSET,IERR)
         IF( IERR .NE. 0 ) GOTO 9020
      ENDIF

C     NOUVEAU TRAITEMENT DE LA SOUS-ETOILE NBSSET
      GOTO 5


C     TRACE FINAL DES TETRAEDRES CREES FORMANT L'ETOILE
C     =================================================
 9000 KTITRE = '      TETRAEDRES CREES en FIN de sietstsu'
      WRITE(KTITRE(1:5),'(I5)') NBTETC

 9010 CALL SANSDBL( KTITRE, L )
      CALL TRFETOV3( KTITRE(1:L), XYZSOM, NBSSET, N1SSET, NFETOI,
     %               NBTETC, NOTETC, NBSOTE, NSTETR )

 9020 PRINT*,'sietstsu: NB TETRAEDRES CREES=',NBTETC,' NOCAS=',NOCAS,
     %       ' IERR=',IERR,' N1TEVI=',N1TEVI
      GOTO 9999


C     ABANDON DE LA TETRAEDRISATION DE L'ETOILE
C     LES NBTETC TETRAEDRES CREES REDEVIENNENT VIDES
C     ==============================================
 9990 NBTETC = - NBTETC

 9999 RETURN
      END
