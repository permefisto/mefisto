      SUBROUTINE TRTR3RE1( LSIGNE, NUIT2,
     %                     NBS1,   XYZS1,  NBTRS1,  NUSTS1, PSGISF1,
     %                     NBS2,   XYZS2,  NBTRS2,  NUSTS2,
     %                     MXSGI,  NSTSGI, MXCHSGI, LCHSGI,
     %                     MXARET, NBARET, LCHARET, MX1ARET, LC1ARET,
     %                     MXPTA,  NBPTA,  XYZPTA,  NBSSGI,
     %                     MXTRA,  NBTRA,  NUSTRA,  IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    SOUS-TRIANGULER LES TRIANGLES INTERSECTES DES 2 SURFACES
C -----    SELON DES SEGMENTS D'INTERSECTION (SGI) ET L'ENTIER LSIGNE
C          ISSU DE L'OPERATEUR LOGIQUE NUOPER
C
C ENTREES:
C --------
C LSIGNE : -1 OU +1 SELON L'OPERATEUR INTERSECTION, ADDITION OU SOUSTRACTION
C NUIT2  : NUMERO DE L'INDICE DANS NSTSGI 3 ou 4 DU NUMERO DE TRIANGLE DE LA
C          SURFACE 2 (1-ER APPEL NUIT2=4, 2-EME APPEL NUIT2=3)
C NBS1   : NOMBRE  DE     SOMMETS DE LA SURFACE 1
C XYZS1  : 3 XYZ DES NBS1 SOMMETS DE LA SURFACE 1
C NBTRS1 : NOMBRE DE TRIANGLES    DE LA SURFACE 1
C NUSTS1 : (4,NBTRS1) NUMERO DANS XYZS1 DES 3 SOMMETS + 0 EN POSITION 4
C PSGISF1: (NBTRS1) POINTEUR SUR 1-ER SEGMENT INTERSECTION DE CHAQUE
C          TRIANGLE DE LA SURFACE 1
C
C NBS2   : NOMBRE  DE     SOMMETS DE LA SURFACE 2
C XYZS2  : 3 XYZ DES NBS2 SOMMETS DE LA SURFACE 2
C NBTRS2 : NOMBRE DE TRIANGLES    DE LA SURFACE 2
C NUSTS2 : (4,NBTRS2) NUMERO DANS XYZS2 DES 3 SOMMETS + 0 EN POSITION 4
C
C MXSGI  : NOMBRE MAXIMAL DE SGI DECLARABLES DANS NSTSGI
C NSTSGI : NUMERO DANS XYZPTA DES 2 POINTS EXTREMITES DU SEGMENT INTERSECTION
C          NUMERO DU TRIANGLE DANS NUSTS1 ET DANS NUSTS2
C
C MXCHSGI: NOMBRE MAXIMAL DE SGI DECLARABLES DANS LCHSGI
C LCHSGI : NUMERO DU SGI et CHAINAGE SUR LE SUIVANT
C
C MXPTA  : NOMBRE MAXIMAL DE POINTS AJOUTES DECLARABLES DANS XYZPTA
C MXARET : NOMBRE MAXIMAL D'ARETES DECLARABLES DANS LCHARET
C MXSTAR : NOMBRE MAXIMAL DE SOMMETS DE SGI SUR LES 3 ARETES D'UN TRIANGLE
C MX1ARET: NOMBRE MAXIMAL DE CF DECLARABLES DANS LC1ARET
C
C AUXILIAIRES:
C ------------
C NBARET : NOMBRE D'ARETES D'UN CF
C LCHARET: (3,MXARET) NUMERO DANS XYZPTA DU SOMMET 1 ET 2 DU SEGMENT,
C          CHAINAGE SUR LE SUIVANT
C          TABLEAU CLONE DU CF POUR SA TRIANGULATION
C LC1ARET: (MX1ARET) INDICE DANS LCHARET DE LA PREMIERE ARETE DE CHAQUE CF
C NBSSGI : (MXPTA) POUR REPERER LES SOMMETS SIMPLES DES SGI
C
C MODIFIES:
C ---------
C NBPTA  : NOMBRE DE  POINTS AJOUTES, SOMMET DES ARETES D'INTERSECTION
C XYZPTA : 3 XYZ  DES POINTS AJOUTES, SOMMET DES ARETES D'INTERSECTION
C MXTRA  : NOMBRE MAXIMAL DE TRIANGLES AJOUTES DECLARABLES DANS NUSTRA
C NBTRA  : NOMBRE DE TRIANGLES AJOUTES
C NUSTRA : (3,NBTRA) NO XYZPTA DES 3 SOMMETS DES NBTRA TRIANGLES
C IERR   : =0 SI PAS D'ERREUR RENCONTREE
C          >0 SI LE CAS N'EST PAS TRAITE => ARRET SANS TRIANGULATION
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET Alain LJLL UPMC & St Pierre du Perray Novembre 2011
C2345X7..............................................................012
      include"./incl/gsmenu.inc"
      include"./incl/langue.inc"
      LOGICAL           TRATRI
      COMMON / DV2DCO / TRATRI
C
      DOUBLE PRECISION  XYZS1(3,NBS1), XYZS2(3,NBS2), XYZPTA(3,MXPTA)
      INTEGER           NSTSGI(4,MXSGI),  LCHSGI(2,MXCHSGI),
     %                  PSGISF1(NBTRS1),
     %                  NUSTS1(4,NBTRS1), NUSTS2(4,NBTRS2),
     %                  NUSTRA(3,MXTRA),  LCHARET(3,MXARET),
     %                  LC1ARET(MX1ARET)
C
C     NUMERO DES SOMMETS UNIQUES DES SGI D'UN TRIANGLE (CAS GENERAL=2)
      INTEGER           NBSSGI(MXPTA), L1COURB(0:32),
     %                  NUSU(64), NUSUCB(2,0:32), NBTRCB(0:32)
      INTEGER           NUARET(2+64), NUARCB(2,0:32)
      EQUIVALENCE      (NUARET(1),NUARCB(1,0))
C
      INTEGER           NPA1, NPA2, NPA3, NPA(3)
      EQUIVALENCE      (NPA(1),NPA1),(NPA(2),NPA2),(NPA(3),NPA3)
      DOUBLE PRECISION  VOLTET, D, VNT1(3), VOL,
     %                  CBINF, CBSUP,
     %                  CBPTTRA(3,64), CBPTCB(3,2,0:32)
      INTRINSIC         SQRT
C
C     TRACE OU NON DES SEGMENTS D'INTERSECTION
ccc      TRATRI = .TRUE.
      TRATRI = .FALSE.
      VOLSIG = 0
C
C     BOUCLE SUR LES TRIANGLES NT1 DE LA SURFACE 1 AYANT DES SGI
C     ==========================================================
      DO 100 NT1=1,NBTRS1
C
C        LA TETE DE CHAINAGE DES SGI DU TRIANGLE NT1
         LCH = PSGISF1( NT1 )
         IF( LCH .EQ. 0 ) GOTO 100
C
C        NOMBRE DE PERMUTATION POUR RENDRE LE TRIANGLE NT1 DECIDABLE
         NBDECI = 0
c
ccc         print *
ccc         print *,'trtr3re1: Les SGI initiaux du triangle NT1=',NT1,
ccc     %    '==================================================='
ccc 1       CALL SGIIMPR( PSGISF1(NT1), LCHSGI, NSTSGI )
C
C        NUMERO DES 3 SOMMETS DU TRIANGLE NT1 CONTENANT AU MOINS 1 SGI
 1       NS1T1 = NUSTS1(1,NT1)
         NS2T1 = NUSTS1(2,NT1)
         NS3T1 = NUSTS1(3,NT1)
C
C        SURFACE ET NORMALE DU TRIANGLE NT1
         CALL VECNOR3D( XYZS1(1,NS1T1), XYZS1(1,NS2T1), XYZS1(1,NS3T1),
     %                  VNT1 )
C
C        AJOUT DES 3 SOMMETS DE NT1 COMME POINTS AJOUTES A XYZPTA
C        --------------------------------------------------------
C        LE NUMERO DE POINT AJOUTE DU SOMMET K EST NPA(K)=NPAk
         IF( NBPTA+3 .GT. MXPTA ) GOTO 9200
         DO K=1,3
            NP = NUSTS1(K,NT1)
            CALL XYZIDEDS( XYZS1(1,NP), NBPTA, XYZPTA, NPA(K) )
            IF( NPA(K) .EQ. 0 ) THEN
C              SOMMET NON IDENTIFIE COMME POINT AJOUTE => IL EST AJOUTE
               NBPTA  = NBPTA + 1
               NPA(K) = NBPTA
               XYZPTA(1,NBPTA) = XYZS1(1,NP)
               XYZPTA(2,NBPTA) = XYZS1(2,NP)
               XYZPTA(3,NBPTA) = XYZS1(3,NP)
            ENDIF
         ENDDO
C
C        NBSGI: NOMBRE DE SGI DU TRIANGLE NT1
         NBSGI = 1
C        NUMERO DANS LCHSGI DU SGI SUIVANT
 2       LCH = LCHSGI( 2, LCH )
         IF( LCH .GT. 0 ) THEN
            NBSGI = NBSGI + 1
            GOTO 2
         ENDIF
C
C        RECHERCHE DES POINTS D'INTERSECTION NUSU APPARTENANT A UN SEUL SGI
C        EXTREMITES DES COURBES DES SGI
C        ------------------------------------------------------------------
         LCH = PSGISF1( NT1 )
C        NOMBRE DE SGI ARRIVANT EN CHAQUE SOMMET DE XYZPTA
         DO K=1,NBPTA
            NBSSGI(K)=0
         ENDDO
C
C        COMPTAGE DU NOMBRE DE FOIS LES SOMMETS SONT VUS DANS LES SGI
C        NUMERO DANS LCHSGI DU SGI
 4       NSGI = LCHSGI( 1, LCH )
C        NSG1 NUMERO DU SOMMET 1 DU SGI DANS XYZPTA
         NSG1 = NSTSGI( 1, NSGI )
         NBSSGI(NSG1) = NBSSGI(NSG1) + 1
C        NSG2 NUMERO DU SOMMET 2 DU SGI DANS XYZPTA
         NSG2 = NSTSGI( 2, NSGI )
         NBSSGI(NSG2) = NBSSGI(NSG2) + 1
C        NUMERO DANS LCHSGI DU SGI SUIVANT
         LCH = LCHSGI( 2, LCH )
         IF( LCH .GT. 0 ) GOTO 4
C
C        NBSU: NOMBRE DE SOMMETS UNIQUES DES SGI DU TRIANGLE NT1
C        VUS UNE SEULE FOIS ET NUMERO NUSU DES SOMMETS UNIQUES
         NBSU = 0
         DO K=1,NBPTA
            IF( NBSSGI(K) .EQ. 1 ) THEN
               NBSU = NBSU + 1
               NUSU(NBSU) = K
            ENDIF
         ENDDO
ccc         print *,'TRIANGLE NT1=',NT1,NBSU,' SOMMETS UNIQUES DES SGI'
ccc         print *,'SOMMETS UNIQUES=',(NUSU(M),M=1,NBSU)
C
         IF( NBSU .GT. 64 ) THEN
            print*,'trtr3re1: PLUS DE 64 SOMMETS UNIQUES NBSU=',NBSU
            print*,'CAS NON TRAITE'
            IERR = 1
            GOTO 9999
         ENDIF
C
         IF( MOD(NBSU,2) .NE. 0 ) THEN
            print*,'trtr3re1: NOMBRE IMPAIR DE SOMMETS UNIQUES NT1=',NT1
            print*,'CAS NON TRAITE'
            IERR = 1
            GOTO 9999
         ENDIF
C
         IF( NBSU .EQ. 0 ) THEN
         print*,'trtr3re1: AU MOINS UNE COURBE FERMEE de SGI dans NT1=',
     %           NT1
         print*,'UNE SURFACE INTERSECTE L''AUTRE A L''INTERIEUR D'' UN T
     %RIANGLE. DESOLE, CAS NON TRAITE'
            IERR = 1
            print*,'MAILLER PLUS FINEMENT LA PLUS GRANDE SURFACE'
            GOTO 9999
C
C            RESTE A AJOUTER ICI LA JONCTION DES COURBES FERMEES
C            INTERNES A UN POINT D'UNE ARETE DE NT1 ...  ?????????
C
         ENDIF
C
C        CALCUL NUARET DU NUMERO D'ARETE NT1 DES NBSU SOMMETS UNIQUES
C        ------------------------------------------------------------
         DO 7 M=1,NBSU
C           COORDONNEES BARYCENTRIQUES DU POINT NUSU(M)
            CALL CBPTTR( XYZS1(1,NS1T1), XYZS1(1,NS2T1),
     %                   XYZS1(1,NS3T1), XYZPTA(1,NUSU(M)),
     %                   CBPTTRA(1,M) )
C           NUSU(M) EST SUR L'ARETE K SI CBPTTRA(K-1)=0
            DO K=1,3
               IF( K .GT. 1 ) THEN
                  K2 = K-1
               ELSE
                  K2 = 3
               ENDIF
               IF( ABS(CBPTTRA(K2,M)) .LE. 1D-4 ) THEN
C                 NUSU(M) EST INTERNE A L'ARETE K DE NT1
                  NUARET(M) = K
ccc                  print *,'Triangle',NT1,' NUSU(',M,')=',NUSU(M),
ccc     %                    ' est sur arete',NUARET(M)
                  GOTO 7
               ENDIF
            ENDDO
 7       CONTINUE
C
C        LES NUSU SONT ORDONNES PAR ARETE 1 A 3 DU TRIANGLE NT1
C        ------------------------------------------------------
         NU = 0
         DO K = 1, 3
            DO M = 1, NBSU
               IF( NUARET(M) .EQ. K ) THEN
C                 PERMUTATIONS POUR FAIRE SE SUIVRE LES SOMMETS UNIQUES
C                 DE L'ARETE K
                  NU = NU + 1
                  L          = NUSU( NU )
                  NUSU( NU ) = NUSU( M )
                  NUSU( M  ) = L
                  L            = NUARET( NU )
                  NUARET( NU ) = NUARET( M )
                  NUARET( M  ) = L
               ENDIF
            ENDDO
         ENDDO
C
C        PERMUTATIONS POUR FAIRE SE SUIVRE DANS LE SENS DIRECT
C        DES 3 ARETES DU TRIANGLE NT1 LES NBSU SOMMETS UNIQUES NUSU DES SGI
C        ------------------------------------------------------------------
         DO K = 1, 3
            IF( K .LT. 3 ) THEN
               K1 = K+1
            ELSE
               K1 = 1
            ENDIF
            NB = 0
            DO 9 M = 1, NBSU
               IF( NUARET(M) .EQ. K ) THEN
                  NB = NB + 1
C                 COORDONNEES BARYCENTRIQUES DU POINT NUSU(M)
                  CALL CBPTTR( XYZS1(1,NS1T1), XYZS1(1,NS2T1),
     %                         XYZS1(1,NS3T1), XYZPTA(1,NUSU(M)),
     %                         CBPTTRA(1,M) )
C                 NUSU(M) EST SUR L'ARETE K
                  NU = M
                  DO LL = 2, NB
                     IF( CBPTTRA(K1,NU) .LT. CBPTTRA(K1,NU-1) ) THEN
C                       L'ABSCISSE CURVILIGNE ISSUE DU SOMMET K
C                       EST CBPTTRA(K1,.)
                        L          = NUSU( NU )
                        NUSU( NU ) = NUSU( NU-1 )
                        NUSU(NU-1) = L
                        L            = NUARET( NU )
                        NUARET( NU ) = NUARET( NU-1 )
                        NUARET(NU-1) = L
                        DO MM=1,3
                           D                = CBPTTRA(MM,NU)
                           CBPTTRA(MM,NU)   = CBPTTRA(MM,NU-1)
                           CBPTTRA(MM,NU-1) = D
                        ENDDO
                        NU = NU-1
                     ELSE
                        GOTO 9
                     ENDIF
                  ENDDO
               ENDIF
 9          CONTINUE
         ENDDO
C
C        CONSTRUCTION A PARTIR DES NUSU(1:NBSU) DE NBCB COURBES CHAINEES
C        ---------------------------------------------------------------
C        ATTENTION: LA FIN DE CHAINAGE DE CHAQUE SUITE
c                   DONNE LE SUIVANT AVEC UN CHAINAGE NEGATIF!...
C                   POUR SIGNALER LE SAUT AU DEBUT DE LA SUITE SUIVANTE
         CALL SGICOURB1( NBSGI, LCHSGI, MXSGI, NSTSGI, PSGISF1(NT1),
     %                   NBSU,  NUSU,   NBCB, L1COURB, NUSUCB )
C
C        NUARCB NUMERO D'ARETE DES 2 SU DE CHACUNE DES COURBES
         DO M=1,NBCB
            DO 10 L=1,2
C              COORDONNEES BARYCENTRIQUES DU POINT NUSUCB(L,M)
               CALL CBPTTR( XYZS1(1,NS1T1), XYZS1(1,NS2T1),
     %                      XYZS1(1,NS3T1), XYZPTA(1,NUSUCB(L,M)),
     %                      CBPTCB(1,L,M) )
C              NUSUCB(L,M) EST SUR L'ARETE K SI CBPTCB(K-1)=0
               DO K=1,3
                  IF( K .GT. 1 ) THEN
                     K2 = K-1
                  ELSE
                     K2 = 3
                  ENDIF
                  IF( ABS(CBPTCB(K2,L,M)) .LE. 1D-4 ) THEN
C                    NUSUCB(L,M) EST INTERNE A L'ARETE K DE NT1
                     NUARCB(L,M) = K
ccc                     print *,'Triangle',NT1,
ccc     %               ' NUSUCB(',L,',',M,')=',NUSUCB(L,M),
ccc     %                      ' est sur arete',NUARCB(L,M)
                     GOTO 10
                  ENDIF
               ENDDO
 10         CONTINUE
         ENDDO
C
C        AJOUT EVENTUEL DES COURBES ISSUES DE SOMMETS DOUBLES
C        SUR UNE ARETE DU TRIANGLE NT1
         NBICB = NBCB
         CALL SGICOURB2( NT1,    NUSTS1, XYZS1,   XYZPTA,
     %                   LCHSGI, NSTSGI,
     %                   NBCB,   L1COURB, NUSUCB, NUARCB )
C
C        REMISE DANS LE SENS DIRECT DES NBCB COURBES
         CALL SGICOURB3( NT1,    NUSTS1,  XYZS1,  XYZPTA,
     %                   LCHSGI, NSTSGI,
     %                   NBCB,   L1COURB, NUSUCB, NUARCB )
C
ccc         print *,'trtr3re1: TRIANGLE',NT1,' ST ajoutes=',NPA
C
C        NUARCB NUMERO D'ARETE DES 2 SU DE CHACUNE DES COURBES
C        ET CBPTCB APRES LE REORDONNANCEMENT DES COURBES
         DO M=1,NBCB
            DO 11 L=1,2
C              COORDONNEES BARYCENTRIQUES DU POINT NUSUCB(L,M)
               CALL CBPTTR( XYZS1(1,NS1T1), XYZS1(1,NS2T1),
     %                      XYZS1(1,NS3T1), XYZPTA(1,NUSUCB(L,M)),
     %                      CBPTCB(1,L,M) )
C              NUSUCB(L,M) EST SUR L'ARETE K SI CBPTCB(K-1)=0
               DO K=1,3
                  IF( K .GT. 1 ) THEN
                     K2 = K-1
                  ELSE
                     K2 = 3
                  ENDIF
                  IF( ABS(CBPTCB(K2,L,M)) .LE. 1D-4 ) THEN
C                    NUSUCB(L,M) EST INTERNE A L'ARETE K DE NT1
                     NUARCB(L,M) = K
ccc                     print *,'Triangle',NT1,
ccc     %               ' NUSUCB(',L,',',M,')=',NUSUCB(L,M),
ccc     %                      ' est sur arete',NUARCB(L,M)
                     GOTO 11
                  ENDIF
               ENDDO
 11         CONTINUE
         ENDDO
C
C        ================================================================
C        CONSTRUCTION DE LA COURBE 0 PARTANT DU SOMMET 1 DE NT1
C        ARRIVANT AU DERNIER SU DES DERNIERES COURBES et LE PLUS PROCHE
C        DU SOMMET 1 SUR L'ARETE 3, ET EN ABSENCE DE COURBES, AU SOMMET 3
C        ================================================================
C        LES ARETES CHAINEES DU CF POUR LA SOUS-TRIANGULATION
         NBARET = 0
C
C        RECHERCHE DU DERNIER SOMMET DES COURBES LE PLUS PROCHE DU SOMMET 1
         NUDERCB = 0
         CBSUP   = 0D0
         DO K=1,NBCB
            IF( NUARCB(2,K) .EQ. 3 ) THEN
               IF( CBPTCB(1,2,K) .GT. CBSUP ) THEN
C                 PLUS GRANDE COORDONNEE BARYCENTRIQUE 1 SUR L'ARETE 3
                  CBSUP   = CBPTCB(1,2,K)
C                 NUMERO DE LA DERNIERE COURBE ARRIVANT AUPRES DU SOMMET 1
                  NUDERCB = K
               ENDIF
            ENDIF
         ENDDO
         IF( NUDERCB .EQ. 0 ) THEN
C           LE SU FINAL DE LA COURBE 0 EST LE SOMMET 3 DE NT1
            NUSUDER = NPA3
            CBPTCB(1,2,0) = 0D0
            CBPTCB(2,2,0) = 0D0
            CBPTCB(3,2,0) = 1D0
         ELSE
C           LE SU FINAL DE LA COURBE 0 EST LE SOMMET FINAL DE LA COURBE NUDERCB
            NUSUDER = NUSUCB(2,NUDERCB)
            CBPTCB(1,2,0) = CBPTCB(1,2,NUDERCB)
            CBPTCB(2,2,0) = CBPTCB(2,2,NUDERCB)
            CBPTCB(3,2,0) = CBPTCB(3,2,NUDERCB)
         ENDIF
C
C        COORDONNEES BARYCENTRIQUES DU SOMMET 1 DE NT1 = 1ER SU COURBE 0
         CBPTCB(1,1,0) = 1D0
         CBPTCB(2,1,0) = 0D0
         CBPTCB(3,1,0) = 0D0
C
C        CARACTERISTIQUES DE LA PSEUDO-COURBE 0
         NUSUCB(1,0) = NPA1
         NUSUCB(2,0) = NUSUDER
C
C        NUASU1 NUMERO DE L'ARETE ACTIVE TRAITEE DE NT1
         NUASU1 = 1
         NUARCB(1,0) = NUASU1
C
C        NUASU2 NUMERO DE L'ARETE FINALE DE LA COURBE NUDERCB
         NUASU2 = 3
         NUARCB(2,0) = NUASU2
C
C        PAS DE POINTEUR SUR LES SGI DE LA COURBE 0
C        CAR IL N'Y A PAS DE SGI SUR LA COURBE 0
         L1COURB(0)=0
C
C        NBTRCB NOMBRE DE TRAITEMENTS DE CHAQUE COURBE DE SGI DE NT1
         DO M=1,NBCB
            NBTRCB(M) = 0
         ENDDO
C        POUR HOMOGENEITE, LA COURBE EST SUPPOSEE DEJA TRAITEE UNE FOIS
         NBTRCB(0)=1
C
C        NOMBRE-1 DE CF TRAITES
         NBCFTR = -1
C
C        ===============================================
C        EXISTE T IL ENCORE UNE COURBE NUCBTR A TRAITER?
C        ===============================================
C        REMARQUE UNE COURBE DOIT ETRE TRAITEE 2 FOIS
C        UNE FOIS DANS LE SENS DIRECT, UNE FOIS DANS LE SENS INDIRECT
 25      DO NUCBTR=0, NBCB
            IF( NBTRCB(NUCBTR) .LT. 2 ) THEN
               NBCFTR = NBCFTR + 1
               GOTO 30
            ENDIF
         ENDDO
C
C        ICI TOUTES LES COURBES DE NT1 ONT ETE TRAITEES
C        FIN DES SOUS-TRIANGULATIONS DU TRIANGLE NT1
         GOTO 100
C
C        ICI IL RESTE AU MOINS LA COURBE NUCBTR A TRAITER
C        ------------------------------------------------
C        NUMERO D'ARETE NT1 DU 1-ER SU DE LA COURBE NUCBTR
 30      NUASU1 = NUARCB(1,NUCBTR)
         NUARK  = NUASU1
C
C        NUMERO D'ARETE NT1 DU DERNIER SU DE LA COURBE NUCBTR
         NUASU2 = NUARCB(2,NUCBTR)
C
C        SOMMET UNIQUE SU DE DEPART DU CHAINAGE
         NUSUI = NUSUCB(1,NUCBTR)
C
C        SOMMET UNIQUE FIN DE LA COURBE NUCBTR
         NUSUF = NUSUCB(2,NUCBTR)
C
C        NUMERO DE LA DERNIERE COURBE INTERNE TRAITEE
         NUDECB = NUCBTR
C
C        COORDONNEE BARYCENTRIQUE DU DERNIER SU
         IF( NUASU1 .LT. 3 ) THEN
            NUCOBA = NUASU1 + 1
         ELSE
            NUCOBA = 1
         ENDIF
         CBINF = CBPTCB(NUCOBA,1,NUCBTR)
         CBSUP = 1D0
C
C        ================================================================
C        TRAITEMENT DE LA COURBE NUCBTR QUI PART DE L'ARETE NUASU1 DE NT1
C        ================================================================
         DO 50 NUAR = NUASU1, NUASU2
C
            IF( NUAR .LT. NUARK ) GOTO 50
C
C           EXISTE T IL SUR L'ARETE NUAR DES COURBES INTERNES
C           A LA COURBE NUCBTR?
C           -------------------------------------------------
C           NO DE LA COORDONNEE BARYCENTRIQUE SUR L'ARETE NUAR
            IF( NUAR .LT. 3 ) THEN
               NUCOBA = NUAR + 1
            ELSE
               NUCOBA = 1
            ENDIF
            IF( NUAR .EQ. NUASU2 ) THEN
               CBSUP = CBPTCB(NUCOBA,2,NUCBTR)
            ENDIF
C
C           TRAITEMENT DES COURBES K AU DELA DE NUDECB
            K = NUDECB
C
 45         K = K + 1
            IF( K .LE. NBCB ) THEN
C
C              LA COURBE K PART ELLE ET ARRIVE T ELLE SUR L'ARETE NUAR?
               IF(NUARCB(1,K) .EQ. NUAR .AND. NUARCB(2,K) .EQ. NUAR)THEN
C
C                 OUI: LES 2 SU DE LA COURBE K SONT SUR L'ARETE NUAR
                  IF( CBPTCB(NUCOBA,1,K) .GE. CBINF  .AND.
     %                CBPTCB(NUCOBA,2,K) .LE. CBSUP ) THEN
C
C                    LA COURBE K EST INTERNE => ELLE EST CHAINEE SENS DIRECT
                     IF( NUSUI .NE. NUSUCB(1,K) ) THEN
C                       CREATION DE LA SOUS-ARETE NUSUI-NUSUCB(1,K)
                        NBARET = NBARET + 1
C                       NUMERO XYZPTA SOMMET 1 DE NT1
                        LCHARET(1,NBARET) = NUSUI
C                       NUMERO XYZPTA SOMMET 2 DE L'ARETE
                        LCHARET(2,NBARET) = NUSUCB(1,K)
C                       CHAINAGE LCHARET DE L'ARETE SUIVANTE
                        LCHARET(3,NBARET) = NBARET + 1
                     ENDIF
C
C                    CHAINAGE DES SGI NUSUCB(1,K)-NUSUCB(2,K)
C                    COMME SOUS ARETES DU CF
                     CALL SGIARET( L1COURB(K), LCHSGI, MXSGI, NSTSGI,
     %                             MXARET, NBARET, LCHARET )
C
C                    UN TRAITEMENT DE PLUS DE LA COURBE K
                     NBTRCB(K) = NBTRCB(K) + 1
C
C                    NUMERO DE LA DERNIERE COURBE INTERNE TRAITEE
                     NUDECB = K
C
C                    SOMMET UNIQUE FIN DE LA COURBE K
                     NUSUI = NUSUCB(2,K)
C
C                    NUMERO DE L'ARETE A TRAITER ENSUITE
                     NUARK = NUAR
C
C                    COORDONNEE BARYCENTRIQUE DU DERNIER SU DE LA COURBE K
                     CBINF = CBPTCB(NUCOBA,2,K)
C
                  ENDIF
C
               ELSE
C
C                 AUTRE COURBE K PARTANT MAIS ARRIVANT SUR NUAR+1 +2 ...
                  IF( NUARCB(1,K) .EQ. NUAR ) THEN
C
C                    LA COURBE K PART ELLE DE L'ARETE K AVANT CBSUP?
                     IF( CBPTCB(NUCOBA,1,K) .GE. CBINF .AND.
     %                   CBPTCB(NUCOBA,1,K) .LE. CBSUP ) THEN
C
C                       OUI: AJOUT DES SGI DE LA CORBE K AU CF
                        IF( NUSUI .NE. NUSUCB(1,K) ) THEN
C                          CREATION DE LA SOUS-ARETE NUSUI-NUSUCB(1,K)
                           NBARET = NBARET + 1
C                          NUMERO XYZPTA SOMMET 1 DE NT1
                           LCHARET(1,NBARET) = NUSUI
C                          NUMERO XYZPTA SOMMET 2 DE L'ARETE
                           LCHARET(2,NBARET) = NUSUCB(1,K)
C                          CHAINAGE LCHARET DE L'ARETE SUIVANTE
                           LCHARET(3,NBARET) = NBARET + 1
                        ENDIF
C
C                       LA COURBE K ARRIVE SUR NUARCB(2,K) DIFFERENTE DE NUAR
C                       CHAINAGE DES SGI NUSUCB(1,K)-NUSUCB(2,K)
C                       COMME SOUS ARETES DU CF
                        CALL SGIARET( L1COURB(K), LCHSGI, MXSGI, NSTSGI,
     %                                MXARET, NBARET, LCHARET )
C
C                       UN TRAITEMENT DE PLUS DE LA COURBE K
                        NBTRCB(K) = NBTRCB(K) + 1
C
C                       NUMERO DE LA DERNIERE COURBE TRAITEE
                        NUDECB = K
C
C                       SOMMET UNIQUE FIN DE LA COURBE K
                        NUSUI = NUSUCB(2,K)
C
C                       NUMERO DE L'ARETE A TRAITER ENSUITE
                        NUARK = NUARCB(2,K)
C
C                       COORDONNEE BARYCENTRIQUE DU DERNIER SU DE LA COURBE K
                        IF( NUARK .LT. 3 ) THEN
                           NUCOBA = NUARK + 1
                        ELSE
                           NUCOBA = 1
                        ENDIF
                        CBINF = CBPTCB(NUCOBA,2,K)
                        GOTO 50
                     ENDIF
C
                  ENDIF
C
               ENDIF
C
C              PASSAGE A LA COURBE SUIVANTE K
               GOTO 45
C
            ENDIF
C
C           AUCUNE AUTRE COURBE K PARTANT DE L'ARETE NUAR DE NT1
            IF( NUAR .LT. 3 ) THEN
               NUCOBA = NUAR + 1
            ELSE
               NUCOBA = 1
            ENDIF
C
C           LE DERNIER SU AVANT LE CHAINAGE INVERSE DE LA COURBE NUCBTR
            IF( NUAR .EQ. NUASU2 ) THEN
C              PASSAGE DU CHAINAGE PAR LE 2-EME SOMMET DE LA COURBE NUCBTR
               NUSUF = NUSUCB(2,NUCBTR)
               CBINF = CBPTCB(NUCOBA,2,NUCBTR)
               CBSUP = 1D0
            ELSE
C              PASSAGE DU CHAINAGE PAR LE SOMMET NUAR+1
               NUSUF = NPA(NUCOBA)
               CBINF = 0D0
               CBSUP = 1D0
            ENDIF
C
            IF( NUSUI .NE. NUSUF ) THEN
C              CREATION DE LA SOUS-ARETE NUSUI-NUSUF
               NBARET = NBARET + 1
C              NUMERO XYZPTA SOMMET 1 DE NT1
               LCHARET(1,NBARET) = NUSUI
C              NUMERO XYZPTA SOMMET 2 DE L'ARETE
               LCHARET(2,NBARET) = NUSUF
C              CHAINAGE LCHARET DE L'ARETE SUIVANTE
               LCHARET(3,NBARET) = NBARET + 1
            ENDIF
C
C           SOMMET UNIQUE A CHAINER
            NUSUI = NUSUF
C
C           NUMERO D'ARETE A TRAITER ENSUITE
            NUARK = NUCOBA
C
 50      CONTINUE
C
C        CHAINAGE FINAL DU CF POUR FERMETURE ET TRIANGULATION SI NECESSAIRE
C        ==================================================================
         IF( NUCBTR .EQ. 0 ) THEN
C
C           TRAITEMENT (SPECIAL) DE LA COURBE 0
C           -----------------------------------
C           CREATION DU SOUS-SEGMENT NUSUCB(2,0)-NPA1 SOMMET 1 DE NT1
            IF( NUSUCB(2,0) .NE. NPA1 ) THEN
C              SOUS-SEGMENT NUSUCB(2,0)-NPA1
               NBARET = NBARET + 1
C              NUMERO XYZPTA SOMMET 1 DE NT1
               LCHARET(1,NBARET) = NUSUCB(2,0)
C              NUMERO XYZPTA SOMMET 2 DE L'ARETE
               LCHARET(2,NBARET) = NPA1
C              CHAINAGE LCHARET DE L'ARETE SUIVANTE
               LCHARET(3,NBARET) = NBARET + 1
            ENDIF
C
C           TRAITEMENT DE LA COURBE 0 SUPPOSE TRAITE 2 FOIS
            NBTRCB(0) = 2
C
C           POUR SAVOIR SI NT1 EST AU DESSUS OU AU DESSOUS DU TRIANGLE NT2
C           DE LA SURFACE 2 CHOIX DU SOMMET NST1 DE NT1 EN DEBUT D'ARETE
C           DU SOMMET INITIAL DE LA COURBE 1
C           --------------------------------------------------------------
            NST1 = NPA( NUARCB(1,1) )
C
C           RECHERCHE DU 1-ER SGI POUR DEFINIR NT2 TRIANGLE DE LA SURFACE 2
            NSGI = LCHSGI( 1, L1COURB(1) )
C           NSG1 NUMERO DU SOMMET 1 DU SGI DANS XYZPTA
            NSG1 = NSTSGI( 1, NSGI )
            IF( NSG1 .EQ. NST1 ) THEN
C
C              NST1 EST UN SOMMET DE NT1 ET LE PREMIER SOMMET DE LA COURBE 1
C              => LA PARTIE SOLIDE AU SOMMET NST1 EST INDECIDABLE!
C              => PERMUTATION DES SOMMETS DU TRIANGLE NT1 et RECALCULS
               NBDECI = NBDECI + 1
               IF( NBDECI .GE. 3 ) THEN
C                 LES 3 SOMMETS SONT INDECIDABLES
                  NBLGRC(NRERR) = 3
                  IF( LANGAG .EQ. 0 ) THEN
        KERR(1)='LES SEGMENTS D''INTERSECTION PASSENT PAR LES 3 SOMMETS'
                  KERR(2)='IMPOSSIBLE DE SAVOIR OU EST LA PARTIE SOLIDE'
                  KERR(3)='REMAILLER UNE SURFACE'
                  ELSE
                   KERR(1)='INTERSECTION SEGMENTS PASSES TO 3 VERTICES'
                   KERR(2)='IMPOSSIBLE TO KNOW WHERE IS THE SOLID'
                   KERR(3)='REMESH a SURFACE'
                  ENDIF
                  CALL LEREUR
                  IERR = 1
                  GOTO 9999
               ELSE
C                 PERMUTATION DES SOMMETS DU TRIANGLE NT1 et RECALCULS
                  NS1T1 = NUSTS1(1,NT1)
                  NUSTS1(1,NT1) = NUSTS1(2,NT1)
                  NUSTS1(2,NT1) = NUSTS1(3,NT1)
                  NUSTS1(3,NT1) = NS1T1
                  GOTO 1
               ENDIF
            ENDIF
C
C           NUMERO DANS NUSTS2 DU TRIANGLE NT2 DU SGI NSGI
            NT2  = NSTSGI( NUIT2, NSGI )
C
C           CALCUL DU DETERMINANT(NST1 S1T2,NST1 S2T2,NST1 S3T2)
C           ou 6 FOIS LE VOLUME SIGNE DU TETRAEDRE NST1 S1T2 S2T2 S3T2
            VOL = VOLTET( XYZPTA(1,NST1),
     %                    XYZS2(1,NUSTS2(1,NT2)),
     %                    XYZS2(1,NUSTS2(2,NT2)),
     %                    XYZS2(1,NUSTS2(3,NT2)) ) * 6D0
C
C           SI LE SOMMET NST1 DU TRIANGLE NT1 EST AU DESSUS DU TRIANGLE NT2
C           VOL = DETERMINANT(NST1 S1T2,NST1 S2T2,NST1 S3T2) >0
C           ALORS LES SOUS-TRIANGLES DU CF SONT CREES DANS  LE TRIANGLE NT1
C           LA COURBE 0 DECIDE DU CF A TRIANGULER OU NON
C           ENSUITE, ALTERNANCE DES SOLIDES (PLEIN-VIDE-PLEIN-VIDE-...)
C           ---------------------------------------------------------------
C           1ER CF A TRIANGULER OU NON
            VOLSIG = REAL( - VOL * LSIGNE )
C
ccc            print *,'trtr3re1: NUARCB(1,1)=',NUARCB(1,1),
ccc     %            ' NSTSGI(',NSGI,')=',(NSTSGI(KM,NSGI),KM=1,4)
ccc            print *,'trtr3re1: NT1=',NT1,' de sommet NST1=',NST1,
ccc     %            ' NT2=',NT2,' NUSTS2=',(NUSTS2(KM,NT2),KM=1,3)
ccc            print *,'VOL=',VOL,' LSIGNE=',LSIGNE
C
            GOTO 70
C
         ENDIF
C
C        SECOND TRAITEMENT DE LA COURBE NUCBTR PARCOURUE EN SENS INVERSE
C        CHAINAGE INVERSE DES SGI NUSUCB(2,NUCBTR)-NUSUCB(1,NUCBTR)
C        ---------------------------------------------------------------
         CALL SGIARINV( L1COURB(NUCBTR), LCHSGI, MXSGI, NSTSGI,
     %                  MXARET, NBARET, LCHARET )
C
C        UN TRAITEMENT DE PLUS DE LA COURBE NUCBTR
         NBTRCB(NUCBTR) = NBTRCB(NUCBTR) + 1
C
C        FERMETURE DU CF
C        ---------------
C        CHAINAGE LCHARET DE L'ARETE SUIVANTE
C        C'EST LA DERNIERE ARETE DU CF
 70      LCHARET(3,NBARET) = 0
C
C        LA PREMIERE ARETE DU CF
         LC1ARET(1) = 1
C
cccC        AFFICHAGE DES SOUS-ARETES CHAINEES DU CF
ccc         print *
ccc         LCH = LC1ARET(1)
ccc 90      IF( LCH .GT. 0 ) THEN
ccc            print *,'SOUS-ARETE du CF', LCHARET(1,LCH), LCHARET(2,LCH),
ccc     %              ' suivante',LCHARET(3,LCH)
ccc            LCH = LCHARET(3,LCH)
ccc            GOTO 90
ccc         ENDIF
C
C        LA COURBE 0 DECIDE DU CF A TRIANGULER OU NON
C        ENSUITE, ALTERNANCE DES SOLIDES (PLEIN-VIDE-PLEIN-VIDE-...)
         VOLSIG = - VOLSIG
ccc         print *,'NT1=',NT1,' NUCBTR=',NUCBTR,' VOLSIG=',VOLSIG
C
         IF( VOLSIG .LT. 0.D0 ) THEN
C
C           SOUS-TRIANGULATION DU CF DU TRIANGLE NT1
C           ----------------------------------------
            NBCF   = 1
            NBTRA0 = NBTRA
            CALL TRIANGCF( VNT1,  NBCF,   LC1ARET, MXARET, LCHARET,
     %                     MXPTA, XYZPTA, MXTRA,   NBTRA,  NUSTRA,
     %                     IERR )
            IF( IERR .EQ. 1 ) GOTO 9400
            IF( IERR .EQ. 2 ) GOTO 9999
C
C           TRACE DES NOUVEAUX SOUS-TRIANGLES DU TRIANGLE NT1
            IF( TRATRI .AND. NBTRA0 .LT. NBTRA ) THEN
               CALL TRATRCF( NT1,    NBTRS1,   NUSTS1, XYZS1,
     %                       XYZPTA, NBTRA0+1, NBTRA,  NUSTRA )
            ENDIF
C
         ENDIF
C
C        RETOUR A LA BOUCLE SUR LES COURBES NUCBTR NON 2 FOIS TRAITEES
         NBARET = 0
         GOTO 25
C
C        FIN DE BOUCLE SUR LES TRIANGLES NT1
 100  CONTINUE
C
      GOTO 9999
C
C     PAS ASSEZ DE POINTS SOMMETS DES SGI DECLARABLES
 9200 NBLGRC(NRERR) = 2
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = 'TRTR3RE1: PAS ASSEZ DE SOMMETS DES SGI DECLARES'
         KERR(2) = 'AUGMENTER mxpta'
      ELSE
         KERR(1) = 'TRTR3RE1: NOT ENOUGH VERTICES OF SGI DECLARED'
         KERR(2) = 'AUGMENT mxpta'
      ENDIF
      CALL LEREUR
      IERR = 2
      GOTO 9999
C
C     PAS ASSEZ DE TRIANGLES AJOUTABLES
 9400 NBLGRC(NRERR) = 1
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = 'TRTR3RE1: PAS ASSEZ DE TRIANGLES AJOUTES DECLARES'
         KERR(2) = 'AUGMENTER mxtra dans oplor3.f'
      ELSE
         KERR(1) = 'TRTR3RE1: NOT ENOUGH DECLARED ADDING TRIANGLES'
         KERR(2) = 'AUGMENT mxtra IN oplor3.f'
      ENDIF
      CALL LEREUR
      IERR = 4
C
 9999 RETURN
      END
