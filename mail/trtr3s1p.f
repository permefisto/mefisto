      SUBROUTINE TRTR3S1P( NBS1,    XYZS1,   NBTRS1,  NUSTS1, PSGISF1,
     %                     L1ARET1, L2ARET1, NARET1,
     %                     MXCHSTA, LCHSTA,  PCHSTA,
     %                     MXSGI,   NSTSGI,  MXCHSGI, LCHSGI,
     %                     MXARET,  NBARET,  LCHARET, MX1ARET, LC1ARET,
     %                     MXPILE,  LAPILE,
     %                     MXPTA,   NBPTA,   XYZPTA,
     %                     MXTRA,   NBTRA,   NUSTRA,
     %                     IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CONSTRUCTION DES TRIANGLES DE SURFACE1
C -----    ADJACENTS PAR UNE ARETE AUX TRIANGLES AVEC SGI DE SURFACE1
C
C ENTREES:
C --------
C NBS1   : NOMBRE  DE     SOMMETS DE LA SURFACE 1
C XYZS1  : 3 XYZ DES NBS1 SOMMETS DE LA SURFACE 1
C NBTRS1 : NOMBRE DE TRIANGLES    DE LA SURFACE 1
C NUSTS1 : (4,NBTRS1) NUMERO DANS XYZS1 DES 3 SOMMETS + 0 EN POSITION 4
C PSGISF1: (NBTRS1) POINTEUR SUR 1-ER SEGMENT INTERSECTION DE CHAQUE
C          TRIANGLE DE LA SURFACE 1
C L1ARET1: NOMBRE DE MOTS PAR ARETE DU TABLEAU NARET1
C L2ARET1: NOMBRE DE ARETES DU TABLEAU NARET1
C NARET1 : NARET1(1,I)= NO DU 1-ER  SOMMET DE L'ARETE
C          NARET1(2,I)= NO DU 2-EME SOMMET > 1-ER  SOMMET
C          NARET1(3,I)= CHAINAGE HACHAGE SUR ARETE SUIVANTE
C          NARET1(4,I)= +-NUMERO DU 1-ER EF CONTENANT CETTE ARETE
C                       + SI ARETE DIRECTE DANS L'ELEMENT, - SINON
C                       0 SI PAS DE 1-ER EF
C          NARET1(5,I)= +-NUMERO DU 2-EME EF CONTENANT CETTE ARETE
C                       + SI ARETE DIRECTE DANS L'ELEMENT, - SINON
C                       0 SI PAS DE 2-EME EF
C          NARET1(6,I)= NO DANS LCHSTA DU 1-ER POINT DE L'ARETE I
C                       0 SI PAS DE POINT
C
C MXCHSTA: NOMBRE MAXIMAL DE POINTEURS SUR LES SOMMETS DES ARETES
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
C MXILE  : NOMBRE MAXIMAL D'ARETES EMPILABLES DANS LAPILE
C MX1ARET: NOMBRE MAXIMAL DE CF DECLARABLES DANS LC1ARET
C
C AUXILIAIRES:
C ------------
C LCHSTA : (2,MXCHSTA) NO SOMMET SUR UNE ARETE et SUIVANT DANS LCHSTA
C PCHSTA : MAX(NBEF1,NBEF2) POINTEUR LCHSTA POUR LES TRIANGLES DE SF1 ou SF2
C NBARET : NOMBRE D'ARETES D'UN CF
C LCHARET: (3,MXARET) NUMERO DANS XYZST DU SOMMET 1 ET 2 DU SEGMENT,
C          CHAINAGE SUR LE SUIVANT
C          TABLEAU CLONE DU CF POUR SA TRIANGULATION
C LC1ARET: (MX1ARET) INDICE DANS LCHARET DE LA PREMIERE ARETE DE CHAQUE CF
C LAPILE : (MXPILE) NO D'ARETE DANS NARET1 DE TRIANGLES A TRIANGULER
C
C MODIFIES:
C ---------
C NBPTA  : NOMBRE DE  POINTS AJOUTES, SOMMET DES ARETES D'INTERSECTION
C XYZPTA : 3 XYZ  DES POINTS AJOUTES, SOMMET DES ARETES D'INTERSECTION
C MXTRA  : NOMBRE MAXIMAL DE TRIANGLES AJOUTES DECLARABLES DANS NUSTRA
C
C SORTIES:
C --------
C NBTRA  : NOMBRE DE TRIANGLES AJOUTES
C NUSTRA : (3,NBTRA) NO XYZPTA DES 3 SOMMETS DES NBTRA TRIANGLES
C IERR   : =0 SI PAS D'ERREUR RENCONTREE, >0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR :PERRONNET Alain LJLL UPMC & St Pierre du Perray Septembre 2011
C2345X7..............................................................012
      include"./incl/gsmenu.inc"
      include"./incl/langue.inc"
      LOGICAL           TRATRI
      COMMON / DV2DCO / TRATRI
C
      DOUBLE PRECISION  XYZS1(3,NBS1), XYZPTA(3,MXPTA)
      INTEGER           NSTSGI(4,MXSGI),  LCHSGI(2,MXCHSGI),
     %                  PSGISF1(NBTRS1),
     %                  NUSTS1(4,NBTRS1),
     %                  NARET1(L1ARET1,L2ARET1),
     %                  LCHSTA(3,MXCHSTA), PCHSTA(*),
     %                  NUSTRA(3,MXTRA),  LCHARET(3,MXARET),
     %                  LC1ARET(MX1ARET),
     %                  LAPILE(MXPILE)
C
      DOUBLE PRECISION  CBPTTRA(3),CBPTTR2(3), VNT(3)
      INTEGER           NST(2), NBARE(3), NOARTR(3)
      EQUIVALENCE      (NBARE(1),NBARE1),(NBARE(2),NBARE2),
     %                 (NBARE(3),NBARE3)
      INTEGER           NPA(3), NPA1, NPA2, NPA3
      EQUIVALENCE      (NPA(1),NPA1),(NPA(2),NPA2),(NPA(3),NPA3)
C
C     TRACE OU NON DES SEGMENTS D'INTERSECTION
ccc      TRATRI = .TRUE.
      TRATRI = .FALSE.
C
C     MISE A ZERO DU POINTEUR SUR LES POINTS INTERNES AUX ARETES
      LIBRE1 = L2ARET1
      DO K = 1, L2ARET1
         NARET1(L1ARET1,K) = 0
      ENDDO
C
C     LE CHAINAGE SUIVANT DES SOMMETS DES ARETES
      NBCHSTA = 0
      DO K = 1, MXCHSTA
         LCHSTA(2,K) = 0
      ENDDO
      DO K = 1, NBTRS1
         PCHSTA(K) = 0
      ENDDO
C
C     BOUCLE SUR LES TRIANGLES DE LA SURFACE 1
C     RECHERCHE DES ARETES DE TRIANGLES NT1 AVEC SGI AYANT DES SOMMETS
C     A L'INTERIEUR DES ARETES DE NT1 ET DE TRIANGLE ADJACENT NT0 SANS SGI
C     MARQUAGE DE L'ARETE D'ADJACENCE DANS LES ARETES DE LA SURFACE AVEC NT0
C     ======================================================================
      NBTRA0 = NBTRA
      NBA    = 0
      DO 100 NT1=1,NBTRS1
C
C        LA TETE DE CHAINAGE DES SGI DU TRIANGLE NT1
         IF( PSGISF1( NT1 ) .EQ. 0 ) GOTO 100
C
C        NT1 EST UN TRIANGLE AVEC SGI
         DO 20 K = 1, 3
C
C           K1 ARETE SUIVANTE DE L'ARETE K DE NT1
            IF( K .LT. 3 ) THEN
               K1 = K+1
            ELSE
               K1 = 1
            ENDIF
C           K2 ARETE PRECEDENTE DE L'ARETE K DE NT1
            IF( K .GT. 1 ) THEN
               K2 = K-1
            ELSE
               K2 = 3
            ENDIF
C
C           LES 2 SOMMETS DE NUMEROS CROISSANTS DE L'ARETE K DE NT1
            NST(1) = NUSTS1(K ,NT1)
            NST(2) = NUSTS1(K1,NT1)
            IF( NST(1) .GT. NST(2) ) THEN
               M      = NST(1)
               NST(1) = NST(2)
               NST(2) = M
            ENDIF
C
C           NOARE NO DE L'ARETE HACHEE DANS NARET1 ARETES DE LA SURFACE
C           PAR ARETE:st1,st2,lien,ef1,ef2,POINTEUR STA
            CALL HACHAG( 2, NST, L1ARET1, L2ARET1, NARET1, 3,
     %                   LIBRE1, NOARE )
            IF( NOARE .LT. 0 ) THEN
               NBLGRC(NRERR) = 1
               WRITE(KERR(5)( 1:10),'(I10)') NST(1)
               WRITE(KERR(5)(11:20),'(I10)') NST(2)
               IF( LANGAG .EQ. 0 ) THEN
                  KERR(1) ='ARETE'//KERR(5)(1:20)//' NON RETROUVEE PARMI
     % LES ARETES'
               ELSE
                  KERR(1) ='EDGE'//KERR(5)(1:20)//' NOT FOUND AMONG THE
     %EDGES'
               ENDIF
               CALL LEREUR
               IERR = 9
               GOTO 100
            ENDIF
C
C           LE TRIANGLE NT0 ADJACENT A L'ARETE K DE NT1
C           -------------------------------------------
            NT0 = ABS( NARET1(5,NOARE) )
            IF( NT0 .EQ. NT1 ) NT0 = ABS( NARET1(4,NOARE) )
            IF( NT0 .EQ. 0 ) GOTO 20
C
C           NT0 A-T-IL DES SGI?
            IF( PSGISF1( NT0 ) .GT. 0 ) GOTO 20
C
C           NON: NT0 N'A PAS DE SGI
C           L'ARETE K DE NT1 SUPPORTE T ELLE DES SOMMETS DE SGI?
C           PARCOURS DES SGI DU TRIANGLE NT1
C           LE SGI LCH A T IL UN SOMMET SUR L'ARETE K DE NT1?
            LCH = PSGISF1( NT1 )
            NB  = 0
C
 10         NSGI = LCHSGI( 1, LCH )
            DO 15 NS=1,2
C
C              NSG NUMERO DU SOMMET NS DU SGI DANS XYZPTA
               NSG = NSTSGI( NS, NSGI )
C
C              COORDONNEES BARYCENTRIQUES DU POINT NSG
               CALL CBPTTR( XYZS1(1,NUSTS1(1,NT1)),
     %                      XYZS1(1,NUSTS1(2,NT1)),
     %                      XYZS1(1,NUSTS1(3,NT1)),
     %                      XYZPTA(1,NSG), CBPTTRA )
C
C              NSG EST SUR L'ARETE K SI CBPTTRA(K2)=0
               IF( ABS(CBPTTRA(K2)) .LE. 1D-4 ) THEN
C
C                 NSG EST INTERNE A L'ARETE K DE NT1
C                 NSG EST IL UN SOMMET DE l'ARETE de NT1?
                  IF( CBPTTRA(K)  .GE. 1D-4 .AND.
     %                CBPTTRA(K1) .LE. 1D0-1D-4 ) THEN
C
C                    NSG EST UN SOUS-SOMMET INTERNE A L'ARETE K DE NT1
                     IF( NBCHSTA .EQ. 0 .OR.
     %                   LCHSTA(2,NBCHSTA) .NE. NSG ) THEN
C                       C'EST SA PREMIERE VISION (EVITER 2 FOIS LE MEME NSG)
                        IF( NBCHSTA .GE. MXCHSTA ) GOTO 9600
                        NBCHSTA = NBCHSTA + 1
C                       RECHERCHE DU NUMERO L DE L'ARETE DE NSG DANS NT0
                        DO L = 1, 3
                           IF( L .LT. 3 ) THEN
                              L1 = L+1
                           ELSE
                              L1 = 1
                           ENDIF
C                          LE NO DES 2 SOMMETS DE L'ARETE L DE NT0
                           M1 = NUSTS1(L ,NT0)
                           M2 = NUSTS1(L1,NT0)
                           IF((M1 .EQ. NST(2) .AND. M2 .EQ. NST(1)) .OR.
     %                        (M1 .EQ. NST(1) .AND. M2 .EQ. NST(2)))THEN
C
C                             CHAINAGE DE L'ARETE L DE NT0 SUPPORT
C                             D'AU MOINS UN SOMMET DE SGI DE NT1
C                             L EST NO DE L'ARETE DANS NT0 DU SOMMET NSG
                              LCHSTA(1,NBCHSTA) = L
C                             NUMERO DU SOMMET NSG
                              LCHSTA(2,NBCHSTA) = NSG
C                             LA SUIVANTE EST LA TETE DE LISTE
                              LCHSTA(3,NBCHSTA) = ABS( PCHSTA( NT0 ) )
C                             CHAINAGE POUR DEVENIR NOUVELLE TETE DE LISTE
                              PCHSTA( NT0 ) = NBCHSTA
                              NB = NB + 1
C
ccc                       print *,'LCH=',LCH,' SOMMET1',NSG,' SUR ARETE',K,
ccc     %              ' de NT1=',NT1,' ARETE',L,' de NT0=',NT0,' NST=',NST
C
                              GOTO 15
                           ENDIF
                        ENDDO
                     ENDIF
                  ENDIF
               ENDIF
 15         ENDDO
C
C           SGI SUIVANT DE NT1
            LCH  = LCHSGI( 2, LCH )
            IF( LCH .GT. 0 ) GOTO 10
C
            IF( NB .EQ. 0 ) THEN
C              ICI NT1 AVEC SGI MAIS PAS DE SOMMET SGI SUR L'ARETE K
               IF( PCHSTA(NT0) .GE. 0 ) THEN
C
C                 NT0 N'A PAS ETE SOUS-TRIANGULE
C                 L'ARETE K DE NT1 TRIANGLE AVEC SGI, SANS SOMMET DE SGI,
C                 EST ELLE UNE ARETE DES SOUS-TRIANGLES DE NT1?
C                 LES 2 SOMMETS DOIVENT ETRE DES SOMMETS AJOUTES
                  CALL XYZIDEDS( XYZS1(1,NST(1)), NBPTA, XYZPTA, NPA1 )
                  IF( NPA1 .GT. 0 ) THEN
C                    LE SOMMET 1 EST UN SOMMET AJOUTE
                     CALL XYZIDEDS(XYZS1(1,NST(2)), NBPTA, XYZPTA, NPA2)
                     IF( NPA2 .GT. 0 ) THEN
C                       LES 2 SOMMETS SONT DES SOMMETS AJOUTES
C                       CELA NE SUFFIT PAS POUR DIRE QUE L'ARETE NPA1-NPA2
C                       EST UNE ARETE DE LA SOUS-TRIANGULATION
                        CALL REARTRI( NPA1, NPA2, NBTRA, NUSTRA,
     %                                NUTRAR, NUAT )
C                       NUTRAR=0 SI L'ARETE NPA1-NPA2 N'EST PAS UNE ARETE
C                              NO DANS NUSTRA DU TRIANGLE AVEC CETTE ARETE
C                              => NUAT=NO ARETE DANS LE TRIANGLE NUTRAR
                        IF( NUTRAR .GT. 0 ) THEN
C
C                          L'ARETE K DE NT1 EST UNE ARETE DES
C                          SOUS-TRIANGLES DE NT1
C                          L'ARETE NOARE EST MARQUEE DANS NARET1 AVEC NT0
C                          AFIN ENSUITE D'AJOUTER LE TRIANGLE NT0 A NUSTRA
                           NARET1(L1ARET1,NOARE) = NT0
ccc                           print *,'SOUS-TRIANGULER NT0=',NT0,
ccc     %                             ' NO ARETE=',NOARE,' St:',NST
                           NBA = NBA + 1
C
                        ENDIF
                     ENDIF
                  ENDIF
C
               ENDIF
            ENDIF
C
 20      CONTINUE
C
C        NT1 EST MARQUE COMME ETANT SOUS-TRIANGULE
         PCHSTA(NT1) = -ABS( PCHSTA(NT1) )
C
 100  CONTINUE
ccc      print *,NBA,' ARETES DE NT1 AVEC SGI DE TRIANGLE ADJACENT A SOUS-T
ccc     %RIANGULER dans BOUCLE 100'
C
C     SOUS-TRIANGULATION DES TRIANGLES NT0 ADJACENTS AUX TRIANGLES AVEC SGI
C     ET DONT AU MOINS UNE ARETE A AU MOINS UN SOMMET DE SGI SUR SON OPPOSE
C     CONSTRUCTION DU CF DU TRIANGLE NT0 ET SOUS-TRIANGULATION
C     MARQUAGE DES ARETES DE NT0 SANS SOMMETS SGI POUR LES TRAVERSER ENSUITE
C     ======================================================================
      NBA = 0
      DO 200 NT0 = 1, NBTRS1
C
C        LA TETE DE CHAINAGE DES SGI DU TRIANGLE NT0
         IF( PCHSTA( NT0 ) .LE. 0 ) GOTO 200
C
C        IL EXISTE AU MOINS UN SOMMET DE SGI SUR UNE ARETE DE NT0
C        ET CE TRIANGLE NT0 EST A SOUS-TRIANGULER
C
C        CONSTRUCTION DU CF A PARTIR DES 3 SOMMETS DE NT0
C        AJOUT DES 3 SOMMETS DE NT0 COMME POINTS AJOUTES SI CE N'EST DEJA FAIT
         DO 111 K=1,3
            NP = NUSTS1(K,NT0)
            CALL XYZIDEDS( XYZS1(1,NP), NBPTA, XYZPTA, NPA(K) )
            IF( NPA(K) .GT. 0 ) GOTO 111
            IF( NBPTA  .GE. MXPTA ) GOTO 9200
C           AJOUT DU SOMMET K AUX SOMMETS AJOUTES
            NBPTA  = NBPTA + 1
            NPA(K) = NBPTA
            XYZPTA(1,NBPTA) = XYZS1(1,NP)
            XYZPTA(2,NBPTA) = XYZS1(2,NP)
            XYZPTA(3,NBPTA) = XYZS1(3,NP)
 111     CONTINUE
C
C        TEMOIN DE SOMMETS DE SGI INTERNE A CHAQUE ARETE DU TRIANGLE NT0
         NOARTR(1) = 0
         NOARTR(2) = 0
         NOARTR(3) = 0
C
C        INITIALEMENT, LE CF EST Sk-Sk1-Sk2
         NBCF = 1
         LC1ARET(1) = 1
         NBARET = 0
         IF( NBARET+3 .GT. MXARET ) GOTO 9100
C
C        Sk-Sk1
         NBARET = NBARET + 1
         NBARE1 = NBARET
C        NUMERO XYZPTA SOMMET 1 DE L'ARETE
         LCHARET(1,NBARET) = NPA1
C        NUMERO XYZPTA SOMMET 2 DE L'ARETE
         LCHARET(2,NBARET) = NPA2
C
C        Sk1-Sk2
         NBARET = NBARET + 1
         NBARE2 = NBARET
C        NUMERO XYZPTA SOMMET 1 DE L'ARETE
         LCHARET(1,NBARET) = NPA2
C        NUMERO XYZPTA SOMMET 2 DE L'ARETE
         LCHARET(2,NBARET) = NPA3
C        CHAINAGE LCHARET DE L'ARETE PRECEDENTE SUR CELLE CI
         LCHARET(3,NBARE1) = NBARE2
C
C        Sk2-Sk
         NBARET = NBARET + 1
         NBARE3 = NBARET
C        NUMERO XYZPTA SOMMET 1 DE L'ARETE
         LCHARET(1,NBARET) = NPA3
C        NUMERO XYZPTA SOMMET 2 DE L'ARETE
         LCHARET(2,NBARET) = NPA1
C        CHAINAGE LCHARET DE L'ARETE PRECEDENTE SUR CELLE CI
         LCHARET(3,NBARE2) = NBARE3
C
C        LA DERNIERE ARETE DU CF DES 3 SOMMETS
         LCHARET(3,NBARE3) = 0
C
C        PARCOURS DES SOMMETS DES SGI DES ARETES DE NT0
C        ----------------------------------------------
         LCH = ABS( PCHSTA( NT0 ) )
C
C        NO DE L'ARETE L DANS NT0 DU SOMMET NSG
 210     L = LCHSTA(1,LCH)
C        L'ARETE L DE NT0 EST TRAITEE
         NOARTR(L) = L
C        NUMERO DU SOMMET NSG SUR L'ARETE L DE NT0
         NSG = LCHSTA(2,LCH)
C
C        OU SE TROUVE NSG PAR RAPPORT AUX SOMMETS DES SGI SUR L'ARETE L?
C        CALCUL DES COORDONNEES BARYCENTRIQUES DE NSG SUR NT0
         CALL CBPTTR( XYZS1(1,NUSTS1(1,NT0)),
     %                XYZS1(1,NUSTS1(2,NT0)),
     %                XYZS1(1,NUSTS1(3,NT0)),
     %                XYZPTA(1,NSG), CBPTTRA )
C        CBPTTRA(L1) EST LA COORDONNEE BARYCENTRIQUE SUR Sl-Sl1
         IF( L .LT. 3 ) THEN
            L1 = L+1
         ELSE
            L1 = 1
         ENDIF
C
C        NO DU SOMMET 2 DE LA PREMIERE SOUS ARETE DE L'ARETE L
         LCHA = NBARE(L)
C
 222     NSG2 = LCHARET(2,LCHA)
C        CALCUL DES COORDONNEES BARYCENTRIQUES DE NSG2 SUR NT0
         CALL CBPTTR( XYZS1(1,NUSTS1(1,NT0)),
     %                XYZS1(1,NUSTS1(2,NT0)),
     %                XYZS1(1,NUSTS1(3,NT0)),
     %                XYZPTA(1,NSG2), CBPTTR2 )
         IF( CBPTTRA(L1) .GT. CBPTTR2(L1) ) THEN
C           ORDRE Sl-NSG-NSG2-Sl1
C           IL FAUT PLACER NSG AU DELA DE NSG2
            LCHA = LCHARET(3,LCHA)
            IF( LCHA .GT. 0 ) GOTO 222
         ENDIF
C
C        ORDRE Sl-...-NSG2-NSG-Sl1 SUR L'ARETE L DE NT0
         IF( NBARET .GE. MXARET ) GOTO 9100
C        AJOUT DE NSG SUR L'ARETE NBARET DE NT0
         NBARET = NBARET + 1
C        NUMERO XYZPTA DES SOMMETS DE L'ARETE SGI
         LCHARET(1,NBARET) = NSG
         LCHARET(2,NBARET) = LCHARET(2,LCHA)
         LCHARET(2,LCHA)   = NSG
C        CHAINAGE LCHARET DE L'ARETE PRECEDENTE SUR CELLE CI
         LCHARET(3,NBARET) = LCHARET(3,LCHA)
         LCHARET(3,LCHA)   = NBARET
C
C        CHAINAGE SUR LE SUIVANT
         LCH = LCHSTA(3,LCH)
         IF( LCH .GT. 0 ) GOTO 210
C
C        SOUS-TRIANGULATION DU CF DU TRIANGLE NT0
C        SURFACE ET NORMALE DU TRIANGLE NT0
ccc         print *
ccc         print *,'SOUS-TRIANGULATION DU CF DU TRIANGLE NT0=',NT0
         CALL VECNOR3D( XYZS1(1,NUSTS1(1,NT0)),
     %                  XYZS1(1,NUSTS1(2,NT0)),
     %                  XYZS1(1,NUSTS1(3,NT0)),
     %                  VNT )
         CALL TRIANGCF( VNT,   NBCF,   LC1ARET, MXARET, LCHARET,
     %                  MXPTA, XYZPTA, MXTRA,   NBTRA,  NUSTRA,
     %                  IERR )
         IF( NBTRA .GT. NBTRA0 ) THEN
            CALL TRATRCF( NT0,    NBTRS1, NUSTS1, XYZS1,
     %                    XYZPTA, 1,      NBTRA, NUSTRA )
         ENDIF
C
C        NT0 EST MARQUE COMME ETANT SOUS-TRIANGULE
         PCHSTA(NT0) = -ABS( PCHSTA(NT0) )
C
C        QUELLES SONT LES ARETES DE NT0 A TRAVERSER ENSUITE?
C        ---------------------------------------------------
         DO K = 1, 3
C
            IF( NOARTR(K) .EQ. 0 ) THEN
C              ARETE K A TRAVERSER CAR SANS SOMMET DE SGI
C
C              K1 ARETE SUIVANTE DE L'ARETE K DE NT0
               IF( K .LT. 3 ) THEN
                  K1 = K+1
               ELSE
                  K1 = 1
               ENDIF
C
C              LES 2 SOMMETS CROISSANTS DE L'ARETE K DE NT0
               NST(1) = NUSTS1(K ,NT0)
               NST(2) = NUSTS1(K1,NT0)
               IF( NST(1) .GT. NST(2) ) THEN
                  M      = NST(1)
                  NST(1) = NST(2)
                  NST(2) = M
               ENDIF
C
C              NOARE NO DE L'ARETE HACHEE DANS NARET1
C              PAR ARETE:st1,st2,lien,ef1,ef2,POINTEUR STA
               CALL HACHAG( 2, NST, L1ARET1, L2ARET1, NARET1, 3,
     %                      LIBRE1, NOARE )
C
C              LE TRIANGLE NT1 ADJACENT A L'ARETE K DE NT0
               NT1 = ABS( NARET1(5,NOARE) )
               IF( NT0 .EQ. NT1 ) NT1 = ABS( NARET1(4,NOARE) )
C
               IF( PCHSTA(NT1) .GE. 0 .AND. PSGISF1(NT1) .LE. 0 ) THEN
C
C                 LE TRIANGLE NT1 ADJACENT A NT0 PAR L'ARETE NOARE
C                 EST SANS SGI. L'ARETE NOARE EST A TRAVERSER
                  NARET1(L1ARET1,NOARE) = NT1
ccc                  print *,'SOUS-TRIANGULER NT1=',NT1,' NOARE=',NOARE
                  NBA = NBA + 1
C
               ENDIF
C
            ENDIF
         ENDDO
C
 200  CONTINUE
ccc      print *,NBA,' ARETES DE NT0 SANS SGI A TRAVERSER dans BOUCLE 200'
C
C     TRACE DES SGI, DES TRIANGLES DE LA SURFACE1 ET DES TRIANGLES AJOUTES
      IF( TRATRI .AND. NBTRA .GT. NBTRA0 ) THEN
         CALL TRASGIS1A( NBS1,  XYZS1, NBTRS1, NUSTS1, MXPTA,  XYZPTA,
     %                   MXTRA, NBTRA, NUSTRA  )
      ENDIF
C
C     EMPILAGE DES
C     ARETES DES TRIANGLES QUI VIENNENT D'ETRE SOUS-TRIANGULES
C     ARETES DE NT1 AVEC SGI SANS SOMMET INTERNE A L'ARETE
C     ARETES DE NT0 SANS SOMMET INTERNE A L'ARETE
C     ARETES DES TRIANGLES ADJACENTS NON NT1 SGI ET NON NT0
C     ========================================================
C     RECHERCHE D'UNE ARETE NOAR A TRAVERSER
      NBA = 0
      DO 400 NOAR = 1, L2ARET1
C
C        L'EVENTUEL TRIANGLE NT0 A SOUS-TRIANGULER DE L'ARETE NOA
         NT0 = NARET1(L1ARET1,NOAR)
         IF( NT0 .LE. 0 ) GOTO 400
C        L'ARETE NOAR EST A TRAVERSER
C
         IF( PCHSTA(NT0) .LT. 0 ) GOTO 400
C        LE TRIANGLE NT0 N'EST PAS ENCORE SOUS-TRIANGULE
C
         IF( PSGISF1(NT0) .GT. 0 ) GOTO 400
C        ICI TRIANGLE NT0 SANS SGI A SOUS-TRIANGULER DANS NUSTRA
C
C        FORMATION DE LA PILE DES ARETES AVEC CETTE PREMIERE ARETE NOAR
C        LA PREMIERE ARETE
         LAPILE(1) = NOAR
         LHPILE = 1
C
C        TRAITEMENT DE L'ARETE EN HAUT DE LA PILE
 250     IF( LHPILE .GT. 0 ) THEN
C
C           L'ARETE EN HAUT DE LA PILE
            NOA = LAPILE( LHPILE )
C           L'ARETE EST DEPILEE
            LHPILE = LHPILE - 1
C
            NT0 = NARET1(L1ARET1,NOA)
C           NT0 NO DU TRIANGLE A SOUS-TRIANGULER
            IF( NT0 .LE. 0 ) GOTO 250
C           TRIANGLE NT0 NON ENCORE SOUS-TRIANGULE
C
            IF( PCHSTA(NT0) .LT. 0 ) GOTO 250
C           TRIANGLE NT0 NON ENCORE SOUS-TRIANGULE
C
            IF( PSGISF1(NT0) .GT. 0 ) GOTO 250
C           TRIANGLE NT0 SANS SGI  NON ENCORE SOUS-TRIANGULE
ccc            print *,'SOUS-TRIANGULATION DU HAUT DE PILE NT0=',NT0,
ccc     %              ' NOA=',NOA
C
C           NT0 EST UN TRIANGLE A SOUS-TRIANGULER
C           NT0 EST MARQUE COMME ETANT SOUS-TRIANGULE
            PCHSTA(NT0) = -NT0
C
C           L'ARETE EST MARQUEE TRAVERSEE
            NARET1(L1ARET1,NOA) = -NT0
C
C           AJOUT DU TRIANGLE NT0 A NUSTRA
            IF( NBTRA .GE. MXTRA ) GOTO 9400
            NBTRA = NBTRA + 1
C
C           AJOUT EVENTUEL DES 3 SOMMETS DE NT0 COMME SOMMETS AJOUTES
            DO 301 K=1,3
               NP = NUSTS1(K,NT0)
               CALL XYZIDEDS( XYZS1(1,NP), NBPTA, XYZPTA, NPA(K) )
               IF( NPA(K) .GT. 0 ) GOTO 301
               IF( NBPTA .GE. MXPTA ) GOTO 9200
C              AJOUT DU SOMMET NP COMME SOMMET AJOUTE
               NBPTA  = NBPTA + 1
               NPA(K) = NBPTA
               XYZPTA(1,NBPTA) = XYZS1(1,NP)
               XYZPTA(2,NBPTA) = XYZS1(2,NP)
               XYZPTA(3,NBPTA) = XYZS1(3,NP)
 301        CONTINUE
C
C           NT0 EST LE SOUS TRIANGLE NBTRA
            NUSTRA(1,NBTRA) = NPA1
            NUSTRA(2,NBTRA) = NPA2
            NUSTRA(3,NBTRA) = NPA3
ccc            print *,'TRIANGLE NUSTRA Sous-Triangle',NBTRA,' NPA=',NPA
C
C           TRACE DU SOUS-TRIANGLE NBTRA DE NUSTRA
            IF( TRATRI ) THEN
               CALL TRATRCF( NT0,    NBTRS1, NUSTS1, XYZS1,
     %                       XYZPTA, NBTRA,  NBTRA,  NUSTRA )
            ENDIF
C
C           AJOUT DES ARETES NON TRAVERSEES DE NT0
C           QUELLES SONT LES ARETES DE NT0 A TRAVERSER?
C           -------------------------------------------
            DO L = 1, 3
C
C              L1 ARETE SUIVANTE DE L'ARETE L DE NT1
               IF( L .LT. 3 ) THEN
                  L1 = L+1
               ELSE
                  L1 = 1
               ENDIF
C
C              LES 2 SOMMETS CROISSANTS DE L'ARETE L DE NT0
               NST(1) = NUSTS1(L ,NT0)
               NST(2) = NUSTS1(L1,NT0)
               IF( NST(1) .GT. NST(2) ) THEN
                  M      = NST(1)
                  NST(1) = NST(2)
                  NST(2) = M
               ENDIF
C
C              NOARE NO DE L'ARETE HACHEE DANS NARET1
C              PAR ARETE:st1,st2,lien,ef1,ef2,POINTEUR STA
               CALL HACHAG( 2, NST, L1ARET1, L2ARET1, NARET1, 3,
     %                      LIBRE1, NOARE )
C
C              LE TRIANGLE NT1 ADJACENT A L'ARETE L DE NT0 (NOARE DANS NARET1)
               NT1 = ABS( NARET1(5,NOARE) )
               IF( NT0 .EQ. NT1 ) NT1 = ABS( NARET1(4,NOARE) )
C
               IF( PCHSTA(NT1) .EQ. 0 ) THEN
C                 LE TRIANGLE NT1 N'EST PAS ENCORE TRAVERSE
                  IF( PSGISF1( NT1 ) .LE. 0 ) THEN
C
C                    LE TRIANGLE NT1 N'A PAS DE SGI
C                    SON ARETE NOARE EST A TRAVERSER
                     NARET1(L1ARET1,NOARE) = NT1
ccc                     print *,'SOUS-TRIANGULER NT1=',NT1,' NOARE=',NOARE,
ccc     %               ' St',NST
C
C                    L'ARETE NOARE EST EMPILEE POUR ETRE TRAVERSEE ENSUITE
                     IF( LHPILE .GE. MXPILE ) GOTO 9700
                     LHPILE = LHPILE + 1
                     LAPILE(LHPILE) = NOARE
                     NBA = NBA + 1
C
                  ENDIF
C
               ENDIF
C
            ENDDO
C
            GOTO 250
C
         ENDIF
C        LA PILE DES ARETES A TRAVERSER EST VIDE
C
 400  CONTINUE
ccc      print *,NBA,' ARETES TRAVERSEES DANS LA BOUCLE 400'
C
C     TRACE DES SGI, DES TRIANGLES DE LA SURFACE 1 ET DES TRIANGLES AJOUTES
      IF( TRATRI ) THEN
         CALL TRASGIS1A( NBS1,  XYZS1, NBTRS1, NUSTS1, MXPTA, XYZPTA,
     %                   MXTRA, NBTRA, NUSTRA  )
      ENDIF
C
      GOTO 9999
C
C     PAS ASSEZ DE SEGMENTS D'ARETES DECLARABLES
 9100 NBLGRC(NRERR) = 2
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) ='TRTR3S1P: PAS ASSEZ D''ARETES LCHARET DECLAREES'
         KERR(2) ='AUGMENTER mxaret'
      ELSE
         KERR(1) ='TRTR3S1P: NOT ENOUGH DECLARED LCHARET EDGES'
         KERR(2) ='AUGMENT mxaret'
      ENDIF
      CALL LEREUR
      IERR = 1
      GOTO 9999
C
C     PAS ASSEZ DE POINTS SOMMETS DES SGI DECLARABLES
 9200 NBLGRC(NRERR) = 2
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = 'TRTR3S1P: PAS ASSEZ DE SOMMETS DECLARES'
         KERR(2) = 'AUGMENTER mxpta'
      ELSE
         KERR(1) = 'TRTR3S1P: NOT ENOUGH DECLARED VERTICES'
         KERR(2) = 'AUGMENT mxpta'
      ENDIF
      CALL LEREUR
      IERR = 2
      GOTO 9999
C
C     PAS ASSEZ DE TRIANGLES AJOUTABLES
 9400 NBLGRC(NRERR) = 1
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = 'TRTR3S1P: PAS ASSEZ DE TRIANGLES AJOUTES DECLARES'
         KERR(2) = 'AUGMENTER mxtra dans oplor3.f'
      ELSE
         KERR(1) = 'TRTR3S1P: NOT ENOUGH DECLARED ADDING TRIANGLES'
         KERR(2) = 'AUGMENT mxtra IN oplor3.f'
      ENDIF
      CALL LEREUR
      IERR = 4
      GOTO 9999
C
C     PAS ASSEZ DE SOMMETS DE SGI SUR LES 3 ARETES D'UN TRIANGLE
 9600 NBLGRC(NRERR) = 1
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = 'TRTR3S1P: PAS ASSEZ DE SOMMETS SGI SUR ARETES'
         KERR(2) = 'AUGMENTER mxchsta dans oplor3.f'
      ELSE
         KERR(1) = 'TRTR3S1P: NOT ENOUGH DECLARED ADDING SGI VERTICES'
         KERR(2) = 'AUGMENT mxchsta IN oplor3.f'
      ENDIF
      CALL LEREUR
      IERR = 6
      GOTO 9999
C
C     PAS ASSEZ DE HAUTEUR DE PILE
 9700 NBLGRC(NRERR) = 1
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = 'TRTR3S1P: PAS ASSEZ DE HAUTEUR DE PILE'
         KERR(2) = 'AUGMENTER mxstar dans oplor3.f'
      ELSE
         KERR(1) = 'TRTR3S1P: NOT ENOUGH HIGH STACK'
         KERR(2) = 'AUGMENT mxstar IN oplor3.f'
      ENDIF
      CALL LEREUR
      IERR = 7
C
 9999 RETURN
      END
