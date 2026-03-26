      SUBROUTINE TRTR3S2T( NBPTSG, MXPT,   XYZPT,   NBTRSG, MXTR, NUSTR,
     %                     MXSGI,  NSTSGI, MXCHSGI, LCHSGI,
     %                     QUALIT, NOANCT, NONOPT,  NODEST, NODETR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    AMELIORER LA QUALITE DES TRIANGLES AUTOUR DES SGI
C -----    ASSURER LA CONTINUITE DES NUMEROS DES SOMMETS 1 A MXPT
C          ET DES TRIANGLES 1 A MXTR DU TABLEAU NUSTR
C
C ENTREES:
C --------
C NBPTSG : NOMBRE DE POINTS SOMMETS DES SGI
C MXPT   : MAXIMUM DE  SOMMETS DES MXTR TRIANGLES
C XYZPT  : 3 XYZ   DES SOMMETS DES MXTR TRIANGLES
C NBTRSG : NOMBRE  DE TRIANGLES DE SOMMETS SGI
C MXTR   : MAXIMUM DE TRIANGLES A TRAITER
C
C MXSGI  : NOMBRE MAXIMAL DE SGI DECLARABLES DANS NSTSGI
C NSTSGI : NUMERO DANS XYZPTA DES 2 POINTS EXTREMITES DU SEGMENT INTERSECTION
C          NUMERO DU TRIANGLE DANS NUSTS1 ET DANS NUSTS2
C MXCHSGI: NOMBRE MAXIMAL DE SGI DECLARABLES DANS LCHSGI
C
C MODIFIES:
C ---------
C NUSTR  : (3,MXTR) NO XYZPT DES 3 SOMMETS DES MXTR TRIANGLES
C LCHSGI : NUMERO DU SGI et CHAINAGE SUR LE SUIVANT
C
C SORTIES:
C --------
C QUALIT : QUALITE DES NODETR TRIANGLES DANS L'ORDRE CROISSANT
C NOANCT : NUMERO ANCIEN OU POSITION DE LA QUALITE DU TRIANGLE
C          NOANCT(1)=NO POSITION POINTEUR SUR QUALIT(1) PLUS PETITE QUALITE
C NONOPT : NOUVEAU NUMERO DE CHAQUE POINT (-NS NOUVEAU SI SUPPRIME)
C NODEST : NUMERO DU DERNIER SOMMET DES TRIANGLES APRES SUPPRESSIONS
C NODETR : NUMERO DU DERNIER TRIANGLE APRES SUPPRESSIONS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET Alain LJLL UPMC & St Pierre du Perray Decembre 2011
C2345X7..............................................................012
      include"./incl/gsmenu.inc"
      include"./incl/langue.inc"
      include"./incl/darete.inc"
      LOGICAL           TRATRI
      COMMON / DV2DCO / TRATRI
C
      DOUBLE PRECISION  XYZPT(3,MXPT)
      INTEGER           NSTSGI(4,MXSGI),  LCHSGI(2,MXCHSGI)
      INTEGER           NUSTR(3,MXTR), NOANCT(MXTR), NONOPT(MXPT),
     %                  NOSOTR(3), NOSOTR1(3), NTRAR(2,2)
      INTEGER           NOTR(16)
      DOUBLE PRECISION  QUALIT(MXTR), D, DG, DM,
     %                  QUALT0, QUALT1, QUALT2
      DOUBLE PRECISION  SURTRD, AIRE1, AIRE2, AIRE3, AIRE4, AIRE(4)
      EQUIVALENCE       (AIRE(1),AIRE1), (AIRE(2),AIRE2),
     %                  (AIRE(3),AIRE3), (AIRE(4),AIRE4)
C
      print *
      print *,'===================================================='
      print *,'DEBUT trtr3s2t: NBPTSG=',NBPTSG,' MXPT=',MXPT,
     %        ' NBTRSG=',NBTRSG,' MXTR=',MXTR,' MXSGI=',MXSGI
C
C     TRACE OU NON DES SEGMENTS D'INTERSECTION
      TRATRI = .TRUE.
ccc      TRATRI = .FALSE.
C
C     NOFOTI NUMERO DANS LX DES FONCTIONS DE LA FONCTION TAILLE_IDEALE(X,Y,Z)
      NOFOTI = NOFOTIEL()
C
C     IDENTITE DU NOUVEAU NO DE SOMMET
      DO M=1,MXPT
         NONOPT( M ) = M
      ENDDO
C
C     RECHERCHE DES 2 TRIANGLES ADJACENTS DE CHAQUE SGI
C     CALCUL DE LONGUEUR SGI/DARETE ou TAILLE_IDEALE du SGI
C     TRI CROISSANT DES SGI SELON CE CRITERE
C     SUPPRESSION DES 2 TRIANGLES et UN SOMMET DU SGI EST DETRUIT
C     -----------------------------------------------------------
      L1CHSGI = 0
      NBSGI = 0
      DO 5 NSGI = 1, MXSGI
         IF( NSTSGI(1,NSGI) .LE. 0 ) GOTO 5
         NBSGI = NBSGI + 1
C        NUMERO DU SGI
         LCHSGI(1,NBSGI) = NSGI
C        CHAINAGE SUR LE SUIVANT
         LCHSGI(2,NBSGI) = L1CHSGI
         L1CHSGI = NBSGI
C
C        CALCUL DE LONGUEUR SGI/DARETE ou TAILLE_IDEALE du SGI
C        NSG NUMERO DU SOMMET NS DU SGI DANS XYZPT
         NSG1 = NSTSGI( 1, NSGI )
         NSG2 = NSTSGI( 2, NSGI )
C
C        TAILLE DE L'ARETE IDEALE COMME MOYENNE DE CELLE AUX 2 SOMMETS
         IF( NOFOTI .GT. 0 ) THEN
            CALL FONVAL( NOFOTI, 3, XYZPT(1,NSG1), NCODEV, DM )
            IF( NCODEV .LE. 0 ) THEN
               IERR = 1
               DM   = 0D0
            ENDIF
C
            CALL FONVAL( NOFOTI, 3, XYZPT(1,NSG2), NCODEV, D )
            IF( NCODEV .LE. 0 ) THEN
               IERR = 1
               D    = 0D0
            ENDIF
C
C           LONGUEUR MOYENNE
            D = ( D + DM ) / 2D0
C
            IF( D .LE. 0D0 ) THEN
               IF( DARETE .LE. 0.D0 ) THEN
                  D=1D0
               ELSE
                  D=DARETE
               ENDIF
            ENDIF
         ELSE
            D = DARETE
         ENDIF
C
C        CALCUL DE LONGUEUR SGI/DARETE ou TAILLE_IDEALE du SGI
         QUALIT(NBSGI) = SQRT( (XYZPT(1,NSG2)-XYZPT(1,NSG1)) **2
     %                       + (XYZPT(2,NSG2)-XYZPT(2,NSG1)) **2
     %                       + (XYZPT(3,NSG2)-XYZPT(3,NSG1)) **2 ) / D
C
         NOANCT(NBSGI) = NBSGI
C
 5    CONTINUE
C
C     TRI CROISSANT DE LONGUEUR SGI/DARETE ou TAILLE_IDEALE du SGI
      CALL TRITRD( NBSGI, QUALIT, NOANCT )
C     NOANCT NUMERO ANCIEN POSITION DE LONGUEUR SGI/DARETE
C     NOANCT(1)=NO POINTEUR SUR QUALIT(1) LA PLUS PETIT LONGUEUR SGI/DARETE
C
      DO NOSGI=1,NBSGI
C
         IF( QUALIT(NOSGI) .LT. 0.2D0 ) GOTO 30
C        NUMERO DANS QUALIT DU PLUS PETIT RAPPORT
         LSGI = NOANCT(NOSGI)
C        NUMERO DU SGI DANS NSTSGI
         NSGI = LCHSGI( 1, LSGI )
C        NSG NUMERO DU SOMMET NS DU SGI DANS XYZPT
         NSG1 = NSTSGI( 1, NSGI )
         NSG2 = NSTSGI( 2, NSGI )
C
C        RECHERCHE DES 2 TRIANGLES ADJACENTS NTR0 NTR1 PAR L'ARETE NSG1-NSG2
         NTR0 = 0
         NTR1 = 0
         DO 10 NT=1,MXTR
            IF( NUSTR(1,NT) .GT. 0 ) THEN
               DO N1=1,3
                  IF( N1 .EQ. 3 ) THEN
                     N2 = 1
                  ELSE
                     N2 = N1+1
                  ENDIF
                  NS1 = NUSTR(N1,NT)
                  NS2 = NUSTR(N2,NT)
                  IF(( NSG1.EQ.NS1 .AND. NSG2.EQ.NS2 ) .OR.
     %               ( NSG2.EQ.NS1 .AND. NSG1.EQ.NS2 ) ) THEN
                     IF( NTR0 .EQ. 0 ) THEN
                        NTR0 = NT
                        GOTO 10
                     ELSE
                        NTR1 = NT
                        GOTO 20
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
 10      CONTINUE
C
C        TRACE DES 2 TRIANGLES AVANT DESTRUCTION
 20      IF( TRATRI ) THEN
            NOTR(1)=NTR0
            NOTR(2)=NTR1
            CALL TRATRINT12( NBPTSG, NONOPT, XYZPT,
     %                       1, NBTRSG, NUSTR, 2, NOTR )
         ENDIF
C
C        LE SIGNE - POUR INDIQUER LA DESTRUCTION DU POINT
         IF( NSG2 .GT. NSG1 ) THEN
            NONOPT( NSG2 ) = -NSG1
            print *,'2T=>0T  2 SOMMETS',NSG1,NSG2,' => st',NSG1,
     %           ' Destruction des 2 triangles',NTR0,NTR1
         ELSE
            NONOPT( NSG1 ) = -NSG2
            print *,'2T=>0T  2 SOMMETS',NSG2,NSG1,' => st',NSG2,
     %           ' Destruction des 2 triangles',NTR0,NTR1
         ENDIF
C
C        DESTRUCTION DU TRIANGLE NTR1
         NUSTR(1,NTR1) = 0
C
C        DESTRUCTION DU TRIANGLE NTR0
         NUSTR(1,NTR0) = 0
C
C        TRACE DES TRIANGLES 1 A NBTRSG
         IF( TRATRI ) THEN
            CALL TRATRINT13( NBPTSG, NONOPT, XYZPT, 1, NBTRSG, NUSTR )
         ENDIF
C
      ENDDO
C
C     SUPPRESSION DU TRIANGLE DE MAUVAISE QUALITE ET DU TRIANGLE
C     ADJACENT PAR SA PLUS PETITE ARETE AVEC 2 SOMMETS DE SGI
C     2T=>0T  => CES 2 SOMMETS SONT CONFONDUS EN L'UN DES DEUX
C     ----------------------------------------------------------
 30   DO 300 NOITER=1,4
C
C     TRI CROISSANT DE LA QUALITE DES TRIANGLES
C     -----------------------------------------
      DO NT=1,MXTR
C        QUALITE DU TRIANGLE NT
         IF( NUSTR(1,NT) .GT. 0 ) THEN
            CALL QUATRID( NUSTR(1,NT), XYZPT, QUALIT(NT) )
         ELSE
            QUALIT(NT) = 2D0
         ENDIF
         NOANCT(NT) = NT
      ENDDO
      CALL TRITRD( MXTR, QUALIT, NOANCT )
C     NOANCT NUMERO ANCIEN POSITION DE LA QUALITE DU TRIANGLE
C     NOANCT(1)=NO POSITION POINTEUR SUR QUALIT(1) LA PLUS PETITE QUALITE
      print *
      print *,'QUALITE DES 30 PLUS MAUVAIS TRIANGLES INITIAUX ITER=',
     % NOITER,' -----------------------------------------------'
      DO NT=1,30
         NTR0 = NOANCT( NT )
         print *,'iter',NOITER,' QUALITE TRIANGLE',NTR0,'=',QUALIT(NT),
     %   ' NUSTR(',NTR0,')=',(NUSTR(L,NTR0),L=1,3)
      ENDDO
cccC     TRACE DES TRIANGLES 1 A NBTRSG
ccc      print *,'ITERATION',NOITER
ccc      print *,'AVANT SUPPRESSION 2TRIANGLES COTE MIN et 2 SOMMETS SGI'
ccc      IF( TRATRI ) THEN
ccc         CALL TRATRINT( NBPTSG, NONOPT, XYZPT, 1, NBTRSG, NUSTR )
ccc      ENDIF
C
C     BOUCLE SUR LES TRIANGLES DE PLUS MAUVAISES QUALITES
      L1 = 0
      DO 40 NT=1,NBTRSG
C
C        NUMERO DU TRIANGLE DE PLUS MAUVAISE QUALITE
         NTR0 = NOANCT( NT )
         IF( NUSTR(1,NTR0) .NE. 0 ) THEN
C
C           NO ACTUEL DES SOMMETS DU TRIANGLE NTR0
            DO K=1,3
               NOSOTR(K) = NOUVNOPT( NUSTR(K,NTR0), NONOPT )
               NUSTR(K,NTR0) = NOSOTR(K)
            ENDDO
C
            CALL QUATRID( NOSOTR, XYZPT, QUALIT(NT) )
C
C           RECHERCHE DU PLUS PETIT ET PLUS GRAND COTE DU TRIANGLE NTR0
            DG = 0D0
            KM = 0
            DM = 1D100
            DO K0=1,3
               IF( K0 .LT. 3 ) THEN
                  K1 = K0+1
               ELSE
                  K1 = 1
               ENDIF
               NS0 = NOSOTR(K0)
               NS1 = NOSOTR(K1)
               D = (XYZPT(1,NS1)-XYZPT(1,NS0)) ** 2
     %           + (XYZPT(2,NS1)-XYZPT(2,NS0)) ** 2
     %           + (XYZPT(3,NS1)-XYZPT(3,NS0)) ** 2
               IF( D .LT. DM ) THEN
                  DM = D
                  KM = K0
               ENDIF
               IF( D .GT. DG ) THEN
                  DG = D
                  KG = K0
               ENDIF
            ENDDO
C
C           RAPPORT LONGUEUR COTE MIN / LONGUEUR COTE MAX
            IF( 3D0*DM .GT. DG ) GOTO 40
            IF( QUALIT(NT) .GT. 0.1D0+0.025*NOITER ) GOTO 40
C
C           LE NUMERO DES 2 SOMMETS DE L'ARETE LA PLUS COURTE
            IF( KM .LT. 3 ) THEN
               K1 = KM+1
            ELSE
               K1 = 1
            ENDIF
            NS0 = NOSOTR(KM)
            NS1 = NOSOTR(K1)
C           LES 2 SOMMETS DE L'ARETE LA PLUS COURTE SONT ILS SGI?
            IF( NS0 .GT. NBPTSG .OR. NS1 .GT. NBPTSG ) GOTO 40
C
C           OUI: RECHERCHE DU TRIANGLE OPPOSE A CE COTE KM DE NTR0
            DO 35 NTR1=1,NBTRSG
               IF( NTR1 .EQ. NTR0 .OR. NUSTR(1,NTR1) .EQ. 0 ) GOTO 35
C
C              NO ACTUEL DES SOMMETS DU TRIANGLE NTR1
               DO L0=1,3
                  NOSOTR1(L0) = NOUVNOPT( NUSTR(L0,NTR1), NONOPT )
                  NUSTR(L0,NTR1) = NOSOTR1(L0)
               ENDDO
C
               DO L0=1,3
                  IF( L0 .LT. 3 ) THEN
                     L1 = L0+1
                  ELSE
                     L1 = 1
                  ENDIF
                  NS2 = NOSOTR1(L0)
                  NS3 = NOSOTR1(L1)
                  IF( (NS0 .EQ. NS3 .AND. NS1 .EQ. NS2) .OR.
     %                (NS0 .EQ. NS2 .AND. NS1 .EQ. NS3) ) THEN
C
C                    TRACE DES 2 TRIANGLES AVANT DESTRUCTION
                     print *
                     print *,'2T=>0T de NUSTR(',NTR0,')=',NOSOTR,
     %                       ' Qualite=', QUALIT(NT)
                     print *,'2T=>0T de NUSTR(',NTR1,')=',NOSOTR1
                     IF( TRATRI ) THEN
                        NOTR(1)=NTR0
                        NOTR(2)=NTR1
                        CALL TRATRINT12( NBPTSG, NONOPT, XYZPT,
     %                                   1, NBTRSG, NUSTR, 2, NOTR )
                     ENDIF
C
C                    L'ARETE K0 DE NTR0 EST L'ARETE L0 DE NTR1
C                    LE SIGNE - POUR INDIQUER LA DESTRUCTION DU POINT
                     IF( NS1 .GT. NS0 ) THEN
                        NONOPT( NS1 ) = -NS0
                        print *,'2T=>0T   2 SOMMETS',NS0,NS1,' en ',NS0,
     %                  ' Destruction des 2 triangles',NTR0,NTR1
                     ELSE
                        NONOPT( NS0 ) = -NS1
                        print *,'2T=>0T   2 SOMMETS',NS1,NS0,' en ',NS1,
     %                  ' Destruction des 2 triangles',NTR0,NTR1
                     ENDIF
C
C                    DESTRUCTION DU TRIANGLE NTR1
                     NUSTR(1,NTR1) = 0
C
C                    DESTRUCTION DU TRIANGLE NTR0
                     NUSTR(1,NTR0) = 0
C
C     TRACE DES TRIANGLES 1 A NBTRSG
      IF( TRATRI ) THEN
         CALL TRATRINT13( NBPTSG, NONOPT, XYZPT, 1, NBTRSG, NUSTR )
      ENDIF
C
                     GOTO 40
C
                  ENDIF
               ENDDO
 35         CONTINUE
         ENDIF
 40   CONTINUE
C
C     TRACE DES TRIANGLES 1 A NBTRSG
      IF( TRATRI ) THEN
         CALL TRATRINT( NBPTSG, NONOPT, XYZPT, 1, NBTRSG, NUSTR )
      ENDIF
C
C     SUPPRESSION DU TRIANGLE DE MAUVAISE QUALITE SI LES
C     2 TRIANGLES ADJACENTS PAR SES 2 PLUS COURTES ARETES
C     FORMENT AVEC LUI UN TRIANGLE GLOBAL
C     ---------------------------------------------------
      DO NT=1,NBTRSG
C        QUALITE DU TRIANGLE NT
         IF( NUSTR(1,NT) .NE. 0 ) THEN
            DO K=1,3
C              MISE A JOUR DU NUMERO DES SOMMETS
               NOSOTR(K)   = NOUVNOPT( NUSTR(K,NT), NONOPT )
               NUSTR(K,NT) = NOSOTR(K)
            ENDDO
            CALL QUATRID( NOSOTR, XYZPT, QUALIT(NT) )
         ELSE
            QUALIT(NT) = 2D0
         ENDIF
         NOANCT(NT) = NT
      ENDDO
      CALL TRITRD( NBTRSG, QUALIT, NOANCT )
C     NOANCT NUMERO ANCIEN POSITION DE LA QUALITE DU TRIANGLE
C     NOANCT(1)=NO POSITION POINTEUR SUR QUALIT(1) LA PLUS PETITE QUALITE
      print *
      print *,'QUALITE DES PLUS MAUVAIS TRIANGLES AVANT 3t=>1T'
      DO NT=1,30
         NTR0 = NOANCT( NT )
         print *,'QUALITE du TRIANGLE',NOANCT(NT),'=',QUALIT(NT),
     %   ' NUSTR(',NTR0,')=',(NOUVNOPT(NUSTR(L,NTR0),NONOPT),L=1,3)
      ENDDO
      print *
      print *,'QUALITE DES PLUS MAUVAIS TRIANGLES AVANT 3T => 1T'
C
      DO 275 NT=1,NBTRSG
C
C        NUMERO DU TRIANGLE DE PLUS MAUVAISE QUALITE
         NTR0 = NOANCT( NT )
         IF( NUSTR(1,NTR0) .LE. 0 ) GOTO 275
C
C        QUALITE DU TRIANGLE NTR0
         DO K0=1,3
            NOSOTR(K0) = NOUVNOPT( NUSTR(K0,NTR0), NONOPT )
            NUSTR(K0,NTR0) = NOSOTR(K0)
         ENDDO
         CALL QUATRID( NOSOTR, XYZPT, QUALIT(NT) )
C
         IF( QUALIT(NT) .GT. 0.025D0*NOITER ) GOTO 275
         print *,'iter=',NOITER,' 3t=>1T? de NUSTR(',NTR0,')=',
     %            NOSOTR,' QUAL=', QUALIT(NT)
C
C        RECHERCHE DE LA PLUS GRANDE ARETE KM DU TRIANGLE NTR0
         KM = 0
         DM = 0D0
         DO K0=1,3
            IF( K0 .LT. 3 ) THEN
               K1 = K0+1
            ELSE
               K1 = 1
            ENDIF
            NS0 = NOSOTR(K0)
            NS1 = NOSOTR(K1)
            D = (XYZPT(1,NS1)-XYZPT(1,NS0)) ** 2
     %        + (XYZPT(2,NS1)-XYZPT(2,NS0)) ** 2
     %        + (XYZPT(3,NS1)-XYZPT(3,NS0)) ** 2
            IF( D .GT. DM ) THEN
               DM = D
               KM = K0
            ENDIF
         ENDDO
C
C        SOMMETS DE L'ARETE SUIVANT LA PLUS GRANDE ARETE DE NTR0
         IF( KM .LT. 3 ) THEN
            K0 = KM+1
         ELSE
            K0 = 1
         ENDIF
C
         DO N=1,2
C
C           RECHERCHE DU TRIANGLE OPPOSE A L'ARETE COURTE K0 DE NTR0
C           POUR SUPPRIMER LE POINT INTERNE A L'UNION DES 3 TRIANGLES
C           NTR0+2 TRIANGLES OPPOSES SI ELLE FORME UN TRIANGLE
            IF( K0 .LT. 3 ) THEN
               K1 = K0+1
            ELSE
               K1 = 1
            ENDIF
            NS0 = NOSOTR(K0)
            NS1 = NOSOTR(K1)
C
            DO NTR1=1,NBTRSG
               IF( NTR1 .NE. NTR0 .AND. NUSTR(1,NTR1) .NE. 0 ) THEN
                  DO L0=1,3
                     NOSOTR1(L0) = NOUVNOPT( NUSTR(L0,NTR1), NONOPT )
                     NUSTR(L0,NTR1) = NOSOTR1(L0)
                  ENDDO
                  DO L0=1,3
                     IF( L0 .LT. 3 ) THEN
                        L1 = L0+1
                     ELSE
                        L1 = 1
                     ENDIF
                     NS2 = NOSOTR1(L0)
                     NS3 = NOSOTR1(L1)
                     IF( (NS0 .EQ. NS3 .AND. NS1 .EQ. NS2) .OR.
     %                   (NS0 .EQ. NS2 .AND. NS1 .EQ. NS3) ) THEN
C
C                       L'ARETE K0 DE NTR0 EST L'ARETE L0 DE NTR1
C                       NUMERO DU TRIANGLE
                        NTRAR(1,N) = NTR1
C                       NUMERO DU SOMMET NON SUR L'ARETE NS0-NS1
                        IF( L0 .NE. 1 ) THEN
                           L1 = L0-1
                        ELSE
                           L1 = 3
                        ENDIF
                        NTRAR(2,N) = NOSOTR1(L1)
                        GOTO 255
C
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
C
 255        K0 = K1
C
         ENDDO
C
C        CES 2 TRIANGLES NTRAR + NTR0 FORMENT ILS UN TRIANGLE?
C        C-A-D PARTAGENT ILS UN MEME SOMMET NON SUR NS0-NS1?
         IF( NTRAR(2,1) .EQ. NTRAR(2,2) ) THEN
C
C           OUI: SUPPRESSION DU POINT CENTRAL NS2
C                DESTRUCTION DES 3 SOUS-TRIANGLES NTRAR
C                CONSTRUCTION DU TRIANGLE NS0-NS1-NS3
C
C          trace des 3 triangles avant destruction
            IF( TRATRI ) THEN
               NOTR(1)=NTRAR(1,1)
               NOTR(2)=NTRAR(1,2)
               NOTR(3)=NTR0
               CALL TRATRINT12( NBPTSG, NONOPT, XYZPT,
     %                          1, NBTRSG, NUSTR, 3, NOTR )
            ENDIF
            IF( KM .LT. 3 ) THEN
               K1 = KM+1
            ELSE
               K1 = 1
            ENDIF
            IF( K1 .LT. 3 ) THEN
               K2 = K1+1
            ELSE
               K2 = 1
            ENDIF
C           SOMMETS DU TRIANGLE NTR0 DANS LA NOUVELLE NUMEROTATION DES POINTS
            NS0 = NOSOTR(KM)
            NS1 = NOSOTR(K1)
            NS2 = NOSOTR(K2)
C           LE SOMMET OPPOSE DANS NTR1
            if( l1 .Le. 0 ) print *,'PB trtr3s2t: L1=',l1
            NS3 = NOSOTR1(L1)
C
C           DESTRUCTION DU TRIANGLE ET CONSTRUCTION DU TRIANGLE GLOBAL
            NUSTR(1,NTR0) = NS0
            NUSTR(2,NTR0) = NS1
            NUSTR(3,NTR0) = NS3
C
C           LE SIGNE - POUR INDIQUER LA DESTRUCTION DU POINT NS2
            NONOPT( NS2 ) = -NS2
C
C           DESTRUCTION DU TRIANGLE NTRAR(1,1)
            NUSTR(1,NTRAR(1,1)) = 0
C
C           DESTRUCTION DU TRIANGLE NTRAR(1,2)
            NUSTR(1,NTRAR(1,2)) = 0
C
            print *,'3t=>1T SUPPRESSION DU SOMMET',NS2,' et 3Triangles',
     %      NTRAR(1,1),NTRAR(1,2),NTR0,' => 1 Triangle Englobant',NTR0
            print *,'Modification NUSTR(',NTR0,')=',NOSOTR
            IF( TRATRI ) THEN
               CALL TRATRINT13( NBPTSG, NONOPT, XYZPT,
     %                          1, NBTRSG, NUSTR )
            ENDIF
C
         ENDIF
 275  CONTINUE
C
C     TRACE DES TRIANGLES 1 A NBTRSG
      IF( TRATRI ) THEN
         CALL TRATRINT( NBPTSG, NONOPT, XYZPT, 1, NBTRSG, NUSTR )
      ENDIF
C
 300  CONTINUE
C
C     SUPPRESSION DU TRIANGLE DE MAUVAISE QUALITE ET DU TRIANGLE
C     ADJACENT PAR SA PLUS PETITE ARETE SI AU MOINS 1 SOMMET DE SGI
C     2t=>0t  => CES 2 SOMMETS SONT CONFONDUS EN L'UN DES DEUX
C     -------------------------------------------------------------
C     test a supprimer apres mise au point du dessous
      IF( NBTRSG .GT. 0 ) GOTO 477
C
      DO NT=1,MXTR
C        QUALITE DU TRIANGLE NT
         IF( NUSTR(1,NT) .NE. 0 ) THEN
            DO K=1,3
C              MISE A JOUR DU NUMERO DES SOMMETS
               NOSOTR(K)   = NOUVNOPT( NUSTR(K,NT), NONOPT )
               NUSTR(K,NT) = NOSOTR(K)
            ENDDO
            CALL QUATRID( NOSOTR, XYZPT, QUALIT(NT) )
         ELSE
            QUALIT(NT) = 2D0
         ENDIF
         NOANCT(NT) = NT
      ENDDO
      CALL TRITRD( MXTR, QUALIT, NOANCT )
C     NOANCT NUMERO ANCIEN POSITION DE LA QUALITE DU TRIANGLE
C     NOANCT(1)=NO POSITION POINTEUR SUR QUALIT(1) LA PLUS PETITE QUALITE
      print *
      print *,'QUALITE DES PLUS MAUVAIS TRIANGLES AVANT 2t=>0t'
      DO NT=1,30
         NTR0 = NOANCT( NT )
         print *,'QUALITE du TRIANGLE',NOANCT(NT),'=',QUALIT(NT),
     %   ' NUSTR(',NTR0,')=',(NOUVNOPT(NUSTR(L,NTR0),NONOPT),L=1,3)
      ENDDO
C
C     BOUCLE SUR LES TRIANGLES DE PLUS MAUVAISES QUALITES
      DO 420 NT=1,NBTRSG
C
C        NUMERO DU TRIANGLE DE PLUS MAUVAISE QUALITE
         NTR0 = NOANCT( NT )
         IF( NUSTR(1,NTR0) .NE. 0 ) THEN
C
C           NO ACTUEL DES SOMMETS DU TRIANGLE NTR0
            DO K=1,3
               NOSOTR(K) = NOUVNOPT( NUSTR(K,NTR0), NONOPT )
               NUSTR(K,NTR0) = NOSOTR(K)
            ENDDO
C
            CALL QUATRID( NOSOTR, XYZPT, QUALIT(NT) )
ccc            IF( QUALIT(NT) .GT. 0.3D0 ) GOTO 420
C
C           RECHERCHE DU PLUS PETIT COTE KM DU TRIANGLE NTR0
            DG = 0D0
            KM = 0
            DM = 1D100
            DO K0=1,3
               IF( K0 .LT. 3 ) THEN
                  K1 = K0+1
               ELSE
                  K1 = 1
               ENDIF
               NS0 = NOSOTR(K0)
               NS1 = NOSOTR(K1)
               D = (XYZPT(1,NS1)-XYZPT(1,NS0)) ** 2
     %           + (XYZPT(2,NS1)-XYZPT(2,NS0)) ** 2
     %           + (XYZPT(3,NS1)-XYZPT(3,NS0)) ** 2
               IF( D .LT. DM ) THEN
                  DM = D
                  KM = K0
               ENDIF
               IF( D .GT. DG ) THEN
                  DG = D
                  KG = K0
               ENDIF
            ENDDO
C
cccC           RAPPORT LONGUEUR COTE MIN / LONGUEUR COTE MAX
ccc            RAP = DM / DG
ccc            IF( 3D0*DM .GT. DG ) GOTO 420
            print *,'DEBUT 2t=>0t? de NUSTR(',NTR0,')=',
     %              NOSOTR,' QUAL=', QUALIT(NT)
C
C           LE NUMERO DES 2 SOMMETS DE L'ARETE LA PLUS COURTE
            IF( KM .LT. 3 ) THEN
               K1 = KM+1
            ELSE
               K1 = 1
            ENDIF
            NS0 = NOSOTR(KM)
            NS1 = NOSOTR(K1)
C           LES 2 SOMMETS DE L'ARETE LA PLUS COURTE NE SONT ILS PAS SGI?
            IF( NS0 .GT. NBPTSG .AND. NS1 .GT. NBPTSG ) GOTO 420
C
C           AU MOINS UN DES SOMMETS DE CETTE COURTE ARETE EST SGI
C           RECHERCHE DU TRIANGLE OPPOSE A CE COTE KM DE NTR0
            DO 410 NTR1=1,NBTRSG
               IF( NTR1 .EQ. NTR0 .OR. NUSTR(1,NTR1) .EQ. 0 ) GOTO 410
C
C              NO ACTUEL DES SOMMETS DU TRIANGLE NTR1
               DO L0=1,3
                  NOSOTR1(L0) = NOUVNOPT( NUSTR(L0,NTR1), NONOPT )
                  NUSTR(L0,NTR1) = NOSOTR1(L0)
               ENDDO
C
               DO L0=1,3
                  IF( L0 .LT. 3 ) THEN
                     L1 = L0+1
                  ELSE
                     L1 = 1
                  ENDIF
                  NS2 = NOSOTR1(L0)
                  NS3 = NOSOTR1(L1)
                  IF( (NS0 .EQ. NS3 .AND. NS1 .EQ. NS2) .OR.
     %                (NS0 .EQ. NS2 .AND. NS1 .EQ. NS3) ) THEN
C
C                    L'ARETE K0 DE NTR0 EST L'ARETE L0 DE NTR1
                     CALL QUATRID( NOSOTR1, XYZPT, QUALIT(NTR1) )
                     IF( MIN(QUALIT(NT),QUALIT(NTR1)).GT.0.2D0 )GOTO 420
C
C                 trace des 2 triangles avant destruction
                     IF( TRATRI ) THEN
                        NOTR(1)=NTR0
                        NOTR(2)=NTR1
                        CALL TRATRINT12( NBPTSG, NONOPT, XYZPT,
     %                                   1, NBTRSG, NUSTR, 2, NOTR )
                     ENDIF
C
C                    LE SIGNE - POUR INDIQUER LA DESTRUCTION DU POINT
                     IF( NS0 .LE. NBPTSG ) THEN
                        NONOPT( NS1 ) = -NS0
                        print *,'2t=>0t  2 SOMMETS',NS0,NS1,' en ',NS0,
     %                  ' Destruction des 2 triangles',NTR0,NTR1
                     ELSE
                        NONOPT( NS0 ) = -NS1
                        print *,'2t=>0t 2 SOMMETS',NS1,NS0,' en ',NS1,
     %                  ' Destruction des 2 triangles',NTR0,NTR1
                     ENDIF
C
C                    DESTRUCTION DU TRIANGLE NTR1
                     NUSTR(1,NTR1) = 0
C
C                    DESTRUCTION DU TRIANGLE NTR0
                     NUSTR(1,NTR0) = 0
            IF( TRATRI ) THEN
               CALL TRATRINT13( NBPTSG, NONOPT, XYZPT,
     %                          1, NBTRSG, NUSTR )
            ENDIF
C
                     GOTO 420
C
                  ENDIF
               ENDDO
 410        CONTINUE
         ENDIF
 420  CONTINUE
C
C     TRACE DES TRIANGLES 1 A NBTRSG
      IF( TRATRI ) THEN
         CALL TRATRINT( NBPTSG, NONOPT, XYZPT, 1, NBTRSG, NUSTR )
      ENDIF
C
C     RENUMEROTATION DES XYZ DES SOMMETS COMPTE TENU DES SUPPRESSIONS DE POINTS
C     =========================================================================
 477  NODEST = 0
      DO K=1,MXPT
         IF( K .EQ. NBPTSG ) THEN
            NBPTSG1 = NODEST
         ENDIF
         IF( NONOPT(K) .GT. 0 ) THEN
C
C           SOMMET NON SUPPRIME. 1 SOMMET DE PLUS
            NODEST = NODEST + 1
            XYZPT(1,NODEST) = XYZPT(1,K)
            XYZPT(2,NODEST) = XYZPT(2,K)
            XYZPT(3,NODEST) = XYZPT(3,K)
C
C           NOUVEAU NUMERO DU SOMMET K
            NOANCT(K) = NODEST
C
         ENDIF
      ENDDO
      print *,'APRES RENUMEROTATION le NOUVEAU NBPTSG=',nbptsg1
C
C     RENUMEROTATION DES SOMMETS DES MXTR TRIANGLES
C     =============================================
      NODETR = 0
      DO NT=1,MXTR
         IF( NT .EQ. NBTRSG ) THEN
C           NOMBRE DE TRIANGLES SGI APRES RENUMEROTATION
            NBTRSG1 = NODETR
         ENDIF
         IF( NUSTR(1,NT) .NE. 0 ) THEN
C
C           TRIANGLE NON SUPPRIME. 1 TRIANGLE DE PLUS
            NODETR = NODETR + 1
            DO L=1,3
C              LE NUMERO ANCIEN DU SOMMET
               NS = NUSTR(L,NT)
C              LE NOUVEAU NUMERO DU SOMMET
               NS = NOUVNOPT( NS, NONOPT )
               NP = NOANCT(NS)
               NUSTR(L,NODETR) = NP
            ENDDO
C
C           VERIFICATION: LE TRIANGLE A T IL 2 SOMMETS IDENTIQUES?
            IF( NUSTR(1,NODETR) .EQ. NUSTR(2,NODETR) .OR.
     %          NUSTR(2,NODETR) .EQ. NUSTR(3,NODETR) .OR.
     %          NUSTR(3,NODETR) .EQ. NUSTR(1,NODETR) ) THEN
               print *,'NUSTR(',NODETR,')=',(NUSTR(L,NODETR),L=1,3),
     %                 ' EST SUPPRIME'
C              LE TRIANGLE EST SUPPRIME
               NODETR = NODETR - 1
            ENDIF
C
         ENDIF
C
      ENDDO
      print *,'APRES RENUMERATION LE NOUVEAU NBTRSG=',NBTRSG1
C
C     CHANGEMENT DE DIAGONALES D'UN TRIANGLE DE MAUVAISE QUALITE
C     ET SON TRIANGLE ADJACENT PAR SA PLUS LONGUE ARETE ET
cccC     SEULEMENT SI UN DES SOMMETS EST SGI ET
C     SEULEMENT SI LE MINIMUM DES QUALITES DES 2 TRIANGLES CROIT
C     ----------------------------------------------------------
      NBECHA = 0
      DO NOITER=1,3
C
      DO NT=1,NODETR
C        QUALITE DU TRIANGLE NT
         CALL QUATRID( NUSTR(1,NT), XYZPT, QUALIT(NT) )
         NOANCT(NT) = NT
      ENDDO
C
      CALL TRITRD( NODETR, QUALIT, NOANCT )
C     NOANCT NUMERO ANCIEN POSITION DE LA QUALITE DU TRIANGLE
C     NOANCT(1)=NO POSITION POINTEUR SUR QUALIT(1) LA PLUS PETITE QUALITE
C
      print *
      print *,'APRES RENUMEROTATION QUALITE DES PLUS MAUVAIS TRIANGLES'
      DO NT=1,30
         NTR0 = NOANCT( NT )
         print *,'QUALITE du TRIANGLE',NOANCT(NT),'=',QUALIT(NT),
     %   ' NUSTR(',NTR0,')=',(NUSTR(L,NTR0),L=1,3)
      ENDDO
C
      DO 500 NT=1,NODETR
C
C        NUMERO DU TRIANGLE DE PLUS MAUVAISE QUALITE
         NTR0 = NOANCT( NT )
C
C        QUALITE DU TRIANGLE NTR0
         CALL QUATRID( NUSTR(1,NTR0), XYZPT, QUALIT(NT) )
C
         IF( QUALIT(NT) .GT. 0.1D0+0.025D0*NOITER ) GOTO 500
C
cccC        AU MOINS UN DES 3 SOMMETS DOIT ETRE SGI
ccc         DO K0=1,3
ccc            IF( NUSTR(K0,NTR0) .LE. NBPTSG1 ) GOTO 480
ccc         ENDDO
ccc         GOTO 500
C
         print *,'2T=>2T? QUALITE=',QUALIT(NT),
     %   ' NUSTR(',NTR0,')=',(NUSTR(L,NTR0),L=1,3)
C
C        RECHERCHE DU PLUS GRAND COTE DU TRIANGLE NTR0
         QUALT0 = QUALIT(NT)
         KM = 0
         DM = 0D0
         DO K0=1,3
            IF( K0 .LT. 3 ) THEN
               K1 = K0+1
            ELSE
               K1 = 1
            ENDIF
            NS0 = NUSTR(K0,NTR0)
            NS1 = NUSTR(K1,NTR0)
            D = (XYZPT(1,NS1)-XYZPT(1,NS0)) ** 2
     %        + (XYZPT(2,NS1)-XYZPT(2,NS0)) ** 2
     %        + (XYZPT(3,NS1)-XYZPT(3,NS0)) ** 2
            IF( D .GT. DM ) THEN
               DM = D
               KM = K0
            ENDIF
         ENDDO
C        SOMMETS DE L'ARETE LA PLUS GRANDE
         IF( KM .LT. 3 ) THEN
            K1 = KM+1
         ELSE
            K1 = 1
         ENDIF
         NS0 = NUSTR(KM,NTR0)
         NS1 = NUSTR(K1,NTR0)
C
C        RECHERCHE DU TRIANGLE OPPOSE A CE COTE KM DE NTR0
C        POUR FAIRE L'ECHANGE DES DIAGONALES
         DO 490 NTR1=1,NODETR
            IF( NTR1 .NE. NTR0 ) THEN
               DO L0=1,3
                  IF( L0 .LT. 3 ) THEN
                     L1 = L0+1
                  ELSE
                     L1 = 1
                  ENDIF
                  NS2 = NUSTR(L0,NTR1)
                  NS3 = NUSTR(L1,NTR1)
                  IF( NS0 .EQ. NS3 .AND. NS1 .EQ. NS2 ) THEN
C
C                    L'ARETE K0 DE NTR0 EST L'ARETE L0 DE NTR1
C                    LES 2 AUTRES SOMMETS DE CES 2 TRIANGLES
                     IF( K1 .LT. 3 ) THEN
                        K2 = K1+1
                     ELSE
                        K2 = 1
                     ENDIF
                     NS2 = NUSTR(K2,NTR0)
C
                     IF( L1 .LT. 3 ) THEN
                        L2 = L1+1
                     ELSE
                        L2 = 1
                     ENDIF
                     NS3 = NUSTR(L2,NTR1)
C
C                    CALCUL DE LA QUALITE DU TRIANGLE NTR1
                     CALL QUATRID( NUSTR(1,NTR1), XYZPT, QUALT1 )
C
C                    LA QUALITE MINIMALE DES 2 TRIANGLES NTR0 NTR1
                     QUALT0 = MIN( QUALT0, QUALT1 )
C
C                    LA QUALITE DES 2 TRIANGLES DE L'AUTRE SUBDIVISION
C                    LE TRIANGLE NS0-NS5-NS2
                     NOSOTR(1) = NS0
                     NOSOTR(2) = NS3
                     NOSOTR(3) = NS2
C                    CALCUL DE LA QUALITE DU TRIANGLE NS0-NS3-NS2
                     CALL QUATRID( NOSOTR, XYZPT, QUALT2 )
                     IF( QUALT2 .LE. QUALT0 ) GOTO 500
C
C                    LE TRIANGLE NS3-NS1-NS2
                     NOSOTR(1) = NS3
                     NOSOTR(2) = NS1
                     NOSOTR(3) = NS2
C                    CALCUL DE LA QUALITE DU TRIANGLE NS3-NS1-NS2
                     CALL QUATRID( NOSOTR, XYZPT, QUALT2 )
                     IF( QUALT2 .LE. QUALT0 ) GOTO 500
C
C                    CALCUL DE L'AIRE1 DU TRIANGLE NTR0
                     AIRE1 = SURTRD( XYZPT(1,NS0), XYZPT(1,NS1),
     %                               XYZPT(1,NS2) )
C
C                    CALCUL DE L'AIRE2 DU TRIANGLE NTR1
                     AIRE2 = SURTRD( XYZPT(1,NS0), XYZPT(1,NS3),
     %                               XYZPT(1,NS1) )
C
C                    CALCUL DE L'AIRE3 DU TRIANGLE NS0-NS3-NS2
                     AIRE3 = SURTRD( XYZPT(1,NS0), XYZPT(1,NS3),
     %                               XYZPT(1,NS2) )
C
C                    CALCUL DE L'AIRE4 DU TRIANGLE NS3-NS1-NS2
                     AIRE4 = SURTRD( XYZPT(1,NS3), XYZPT(1,NS1),
     %                               XYZPT(1,NS2) )
C
C                    L'UNION DES 2 TRIANGLES EST ELLE UN QUADRANGLE CONVEXE?
                     AIRE2 = AIRE1 + AIRE2
                     AIRE3 = AIRE3 + AIRE4
                     print *,'TRIANGLES',ntr0,ntr1,' AIRES',AIRE2,AIRE3,
     %                       ABS(AIRE2-AIRE3)/AIRE3
                     IF( ABS( AIRE2 - AIRE3 ) .LE. AIRE3*1D-3 ) THEN
C
c                       QUADRANGLE CONVEXE
C
C                 trace des 2 triangles avant destruction
                     IF( TRATRI ) THEN
                        NOTR(1)=NTR0
                        NOTR(2)=NTR1
                        CALL TRATRINT12( NBPTSG, NONOPT, XYZPT,
     %                                   1, NBTRSG, NUSTR, 2, NOTR )
                     ENDIF
C                       LES TRIANGLES NS0-NS3-NS2  NS3-NS1-NS2
C                       ONT UNE MEILLEURE MINIMALE QUALITE
C                       => ECHANGE NTR0 NTR1 EN  NS0-NS3-NS2 NS3-NS1-NS2
                        NUSTR(1,NTR0)=NS0
                        NUSTR(2,NTR0)=NS3
                        NUSTR(3,NTR0)=NS2
C
                        NUSTR(1,NTR1)=NS3
                        NUSTR(2,NTR1)=NS1
                        NUSTR(3,NTR1)=NS2
C
                       print *,'2T=>2T TRIANGLES',NTR0,NTR1,
     %                   ' ECHANGE DIAGONALE',NS0,NS1,' PAR',NS3,NS2
                        NBECHA = NBECHA + 1
C
                     ENDIF
            IF( TRATRI ) THEN
               CALL TRATRINT13( NBPTSG, NONOPT, XYZPT,
     %                          1, NBTRSG, NUSTR )
            ENDIF
C
                     GOTO 500
C
                  ENDIF
               ENDDO
            ENDIF
 490     CONTINUE
C
         IF( TRATRI ) THEN
            print *
            print *,'FIN trtr3s2t: TRACE FINAL DES',NBTRSG1,
     %   ' PREMIERS TRIANGLES sur',NODETR,' TRIANGLES'
            DO K=1,NODEST
               NONOPT(K) = K
            ENDDO
            CALL TRATRINT( NODEST, NONOPT, XYZPT, 1, NBTRSG1, NUSTR )
         ENDIF
C
 500  CONTINUE
C
      print *,NBECHA,' ECHANGES DES DIAGONALES 2T => 2T'
      ENDDO
C
      RETURN
      END
