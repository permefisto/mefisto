      SUBROUTINE TRTR3SGI( NBS1,    XYZS1,   NBTRS1,  NUSTS1, PSGISF1,
     %                     NBS2,    XYZS2,   NBTRS2,  NUSTS2, PSGISF2,
     %                     MXSGI,   NBSGI,   NSTSGI,
     %                     MXCHSGI, NBCHSGI, LCHSGI,
     %                     MXPTA,   NBPTA,   XYZPTA,  IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DES SEGMENTS D'INTERSECTION ENTRE LES TRIANGLES NT1
C -----    DE LA SURFACE 1 ET LES TRIANGLES NT2 DE LA SURFACE 2
C
C ENTREES:
C --------
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
C PSGISF2: (NBTRS2) POINTEUR SUR 1-ER SEGMENT INTERSECTION DE CHAQUE
C          TRIANGLE DE LA SURFACE 2
C
C MXSGI  : NOMBRE MAXIMAL DE SGI DECLARABLES DANS NSTSGI
C MXCHSGI: NOMBRE MAXIMAL DE SGI DECLARABLES DANS LCHSGI
C MXPTA  : NOMBRE MAXIMAL DE POINTS AJOUTES DECLARABLES DANS XYZPTA
C
C SORTIES:
C --------
C NBSGI  : NOMBRE DE SEGMENTS D'INTERSECTION
C NSTSGI : NUMERO DANS XYZPTA DES 2 POINTS EXTREMITES DU SEGMENT INTERSECTION
C          NUMERO DU TRIANGLE DANS NUSTS1 ET DANS NUSTS2
C NBCHSGI: NOMBRE DE CHAINAGE DE SGI DANS LCHSGI
C LCHSGI : NUMERO DU SGI et CHAINAGE SUR LE SUIVANT
C NBPTA  : NOMBRE DE  POINTS AJOUTES, SOMMET DES ARETES D'INTERSECTION
C XYZPTA : 3 XYZ  DES POINTS AJOUTES, SOMMET DES ARETES D'INTERSECTION
C IERR   : =0 SI PAS D'ERREUR RENCONTREE, >0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET Alain LJLL UPMC& St Pierre du Perray SEPTEMBRE 2011
C2345X7..............................................................012
      include"./incl/gsmenu.inc"
      include"./incl/langue.inc"
      LOGICAL           TRATRI
      COMMON / DV2DCO / TRATRI
C
      DOUBLE PRECISION  XYZS1(3,NBS1), XYZS2(3,NBS2), XYZPTA(3,MXPTA)
      INTEGER           NSTSGI(4,MXSGI),   LCHSGI( 2, MXCHSGI ),
     %                  PSGISF1( NBTRS1 ), PSGISF2( NBTRS2 ),
     %                  NUSTS1(4,NBTRS1),  NUSTS2(4,NBTRS2)
C
      INTEGER           NOCODE(3), NBPTI(3), KK(3)
      DOUBLE PRECISION  D, PTI(3,3), CBPTIA, CBPTIT1(3,3), CBPTIT2(3,3),
     %                  DISPLA1, DISPLA2, STSGT2(3,3), STSGI(3,2)
      INTRINSIC         SQRT
C
C     TRACE OU NON DES SEGMENTS D'INTERSECTION
ccc      TRATRI = .TRUE.
      TRATRI = .FALSE.
C
      NBPTA   = 0
      NBSGI   = 0
      NBCHSGI = 0
C
C     MISE A ZERO DES POINTEURS SUR LES SGI DE CHAQUE TRIANGLE
      DO K=1,NBTRS1
         PSGISF1(K) = 0
      ENDDO
      DO K=1,NBTRS2
         PSGISF2(K) = 0
      ENDDO
C
      DO 200 NT1=1,NBTRS1
C
C        NUMERO DES 3 SOMMETS DU TRIANGLE NT1
         NS1T1 = NUSTS1(1,NT1)
         NS2T1 = NUSTS1(2,NT1)
         NS3T1 = NUSTS1(3,NT1)
C
         DO 100 NT2=1,NBTRS2
C
C           NUMERO DES 3 SOMMETS DU TRIANGLE NT2
            NS1T2 = NUSTS2(1,NT2)
            NS2T2 = NUSTS2(2,NT2)
            NS3T2 = NUSTS2(3,NT2)
C
ccc            if( nt1 .eq. 563 .and. nt2 .eq. 561 ) then
ccc              print *,'nt1=',nt1,' nt2=',nt2,' st',NS1T2,NS2T2,NS3T2
ccc            endif
C
C           RECHERCHE DES EVENTUELS 3 POINTS D'INTERSECTION DES 3 ARETES
C           DU TRIANGLE NT1 AVEC LE PLAN DU TRIANGLE NT2
            DO K=1,3
C
C              NUMERO DE SOMMET 1 DE L'ARETE K DU TRIANGLE NT1
               NS1AK1 = NUSTS1(K,NT1)
               IF( K .EQ. 3 ) THEN
                  KP = 1
               ELSE
                  KP = K + 1
               ENDIF
C              NUMERO DE SOMMET 2 DE L'ARETE K DU TRIANGLE NT1
               NS2AK1 = NUSTS1(KP,NT1)
C
C              CALCUL DU POINT D'INTERSECTION DE L'ARETE K DU
C              TRIANGLE NT1 AVEC LE PLAN DU TRIANGLE NT2
               CALL INDRPL( XYZS1(1,NS1AK1), XYZS1(1,NS2AK1),
     %                   XYZS2(1,NS1T2), XYZS2(1,NS2T2), XYZS2(1,NS3T2),
     %                   PTI(1,K), NOCODE(K) )
C              NOCODE  : 1 SI LA DROITE EST PARALLELE AU PLAN
C                        2 SI S1=S2  SUPPOSE ICI IMPOSSIBLE
C                        0 SI LE POINT D'INTERSECTION A ETE CALCULE
               IF( NOCODE(K) .EQ. 1 ) THEN
C
C                 LA DROITE NS1AK1-NS2AK1 // EST ELLE DANS LE PLAN DE NT2?
                  D = SQRT( (XYZS1(1,NS2AK1)-XYZS1(1,NS1AK1))**2
     %                    + (XYZS1(2,NS2AK1)-XYZS1(2,NS1AK1))**2
     %                    + (XYZS1(3,NS2AK1)-XYZS1(3,NS1AK1))**2 )
                  CALL DIPTPL( XYZS1(1,NS1AK1),
     %                         XYZS2(1,NS1T2), XYZS2(1,NS2T2),
     %                         XYZS2(1,NS3T2), DISPLA1 )
                  DISPLA1 = DISPLA1 / D
                  CALL DIPTPL( XYZS1(1,NS2AK1),
     %                         XYZS2(1,NS1T2), XYZS2(1,NS2T2),
     %                         XYZS2(1,NS3T2), DISPLA2 )
                  DISPLA2 = DISPLA2 / D
                  IF( DISPLA1 .LT. 1D-5 .AND. DISPLA2 .LT. 1D-5 ) THEN
C
C                    ARETE K DU TRIANGLE NT1 EST DANS LE PLAN DU TRIANGLE NT2
                     NOCODE(K) = -1
C
C                    LES 2 SOMMETS DE L'ARETE SONT 2 POINTS D'INTERSECTION
C                    S1T1 = PTI(K)
                     PTI(1,K) = XYZS1(1,NS1AK1)
                     PTI(2,K) = XYZS1(2,NS1AK1)
                     PTI(3,K) = XYZS1(3,NS1AK1)
C
C                    S2T1 = PTI(K+1)
                     PTI(1,KP) = XYZS1(1,NS2AK1)
                     PTI(2,KP) = XYZS1(2,NS2AK1)
                     PTI(3,KP) = XYZS1(3,NS2AK1)
C
ccc               ELSE
C
C                    MARQUEUR D'ARETE SANS INTERSECTION // PLAN
C                    ET NON DANS LE PLAN
ccc                  NOCODE(K) = 1
C
                  ENDIF
C
               ENDIF
C
            ENDDO
C
C           CALCUL DES COORDONNEES BARYCENTRIQUES DES PTI DANS NT1 et NT2
            DO K = 1, 3
C
               IF( NOCODE(K) .NE. 1 ) THEN
C                 QUELLE EST LA POSITION DU POINT D'INTERSECTION PTI
C                 PAR RAPPORT AU TRIANGLE NT1?
C                 CALCUL DES COORDONNEES BARYCENTRIQUES DE PTI DANS NT1
                  CALL CBPTTR( XYZS1(1,NS1T1), XYZS1(1,NS2T1),
     %                         XYZS1(1,NS3T1), PTI(1,K), CBPTIT1(1,K) )
ccc               print *,'NT1=',NT1,'k=',k,' CBPTI=',(CBPTIT1(M,K),M=1,3)
C
C                 QUELLE EST LA POSITION DU POINT D'INTERSECTION PTI
C                 PAR RAPPORT AU TRIANGLE NT2?
C                 CALCUL DES COORDONNEES BARYCENTRIQUES DE PTI DANS NT2
                  CALL CBPTTR( XYZS2(1,NS1T2), XYZS2(1,NS2T2),
     %                         XYZS2(1,NS3T2), PTI(1,K), CBPTIT2(1,K) )
ccc               print *,'NT2=',NT2,'k=',k,' CBPTI=',(CBPTIT2(M,K),M=1,3)
               ENDIF
C
            ENDDO
C
C           BILAN DES POINTS D'INTERSECTION ENTRE LES 3 ARETES
C           DU TRIANGLE NT1 ET LE PLAN DU TRIANGLE NT2
C           NOCODE(K) = 1 SI LA DROITE EST PARALLELE AU PLAN, NON DANS LE PLAN
C                         => PAS DE POINT D'INTERSECTION
C                      -1 SI LA DROITE EST PARALLELE AU PLAN, ET  DANS LE PLAN
C                         PTI(K)=S1T1  PTI(K+1)=S2T1 => 2 POINTS D'INTERSECTION
C                       2 SI S1=S2  SUPPOSE ICI IMPOSSIBLE
C                       0 SI LE POINT D'INTERSECTION A ETE CALCULE
C                         => UN POINT D'INTERSECTION PTI(K)
C
C           NOMBRE D'ARETES // PLAN DE NT2 NON DANS LE PLAN
            KA    = 0
            NBA1P = 0
            DO K = 1, 3
               IF( NOCODE(K) .EQ. 1 ) THEN
                  NBA1P = NBA1P + 1
                  KA    = K
               ENDIF
            ENDDO
C
            IF( NBA1P .GE. 2 ) THEN
C
C              CAS 1: AU MOINS 2 ARETES PARALLELES AU PLAN NT2
C                     QUI NE SONT PAS DANS LE PLAN NT2
C                     PAS DE POINT D'INTERSECTION DE NT1 AVEC NT2
               GOTO 100
C
            ELSE IF( NBA1P .EQ. 1 ) THEN
C
C              1 ARETE PARALLELE AU PLAN NT2 NON DANS LE PLAN
C              LES 2 ARETES SUIVANTES K2 K3 INTERSECTENT LE PLAN NT2
               IF( KA .EQ. 1 ) THEN
                  K2 = 2
                  K3 = 3
               ELSE IF( KA .EQ. 2 ) THEN
                  K2 = 3
                  K3 = 1
               ELSE
                  K2 = 1
                  K3 = 2
               ENDIF
C
C              ARETE K2 DE NT1 NON // A NT2
               NS1 = NUSTS1(K2,NT1)
               NS2 = NUSTS1(K3,NT1)
               CALL CBPTAR( XYZS1(1,NS1), XYZS1(1,NS2),
     %                      PTI(1,K2), CBPTIA )
               IF( CBPTIA .LE. -1D-14 .OR. CBPTIA .GE. 1D0 ) THEN
C                 LE POINT D'INTERSECTION PTI(K2) EST EXTERNE A L'ARETE K2
C                 =>PAS DE POINT D'INTERSECTION DE NT1 AVEC NT2
                  GOTO 100
               ENDIF
C
            ELSE
C
C              3 POINTS D'INTERSECTION = AUCUNE ARETE DE NT1 // PLAN DE NT2
C              RECHERCHE DES CAS OU 2 ARETES K2 K3 INTERSECTENT LE PLAN NT2
C              EN UN POINT INTERIEUR AUX ARETES DE NT1
               DO K = 1, 3
                  NS1 = NUSTS1(K,NT1)
                  IF( K .EQ. 3 ) THEN
                     KP = 1
                  ELSE
                     KP = K + 1
                  ENDIF
                  NS2 = NUSTS1(KP,NT1)
                  CALL CBPTAR( XYZS1(1,NS1), XYZS1(1,NS2),
     %                         PTI(1,K), CBPTIA )
                  IF( CBPTIA .GE. -1D-14 .AND. CBPTIA .LE. 1D0 ) THEN
C                    LE POINT D'INTERSECTION EST INTERNE A L'ARETE K DE NT1
                     KK(K) = K
                  ELSE
C                    LE POINT D'INTERSECTION N'EST PAS INTERNE A L'ARETE K DE NT
                     KK(K) = 0
                  ENDIF
               ENDDO
               IF( KK(1)+KK(2)+KK(3) .EQ. 0 ) THEN
C                 LES 3 POINTS D'INTERSECTION SONT EXTERNES A NT1
C                 PAS DE POINT D'INTERSECTION DE NT1 AVEC NT2
                  GOTO 100
               ELSE IF( KK(1) .EQ. 0 ) THEN
                  K2 = 2
                  K3 = 3
               ELSE IF( KK(2) .EQ. 0 ) THEN
                  K2 = 3
                  K3 = 1
               ELSE
                  K2 = 1
                  K3 = 2
               ENDIF
C
            ENDIF
C
C           ICI LES ARETES K2 et K3 DE NT1 COUPENT LE PLAN DE NT2
C           EN 2 POINTS PTI(K2) PTI(K3) INTERNES (SOMMETS COMPRIS) A NT1
C           CES 2 POINTS PEUVENT ETRE CONFONDUS EN UN SEUL SOMMET DE NT1
            CALL XYZIDED( PTI(1,K2), PTI(1,K3), IDENTQ )
C           IDENTQ=1 SI LES 2 POINTS SONT JUGES IDENTIQUES, 0 SINON
            IF( IDENTQ .EQ. 1 ) THEN
C              LES 2 ARETES INTERSECTENT LE PLAN NT2 EN UN SEUL POINT
C              QUI EST DONC LE SOMMET K3 DE NT1
C              => PAS DE SEGMENT D'INTERSECTION DE NT1 AVEC NT2
               GOTO 100
            ENDIF
C
C           ICI LES ARETES K2 et K3 DE NT1 COUPENT LE PLAN DE NT2
C           EN 2 POINTS DISTINCTS PTI(K2) et PTI(K3) INTERNES AUX ARETES
C           (SOMMETS DE NT1 COMPRIS) MAIS NON REDUIT AU SOMMET COMMUN
C           DES 2 ARETES K2 K3 DE NT1
C           COMMENT EST PLACE PTI(K2) PTI(K3) PAR RAPPORT AU TRIANGLE NT2?
C           RECHERCHE DANS PTI(K2)-PTI(K3) DU SEGMENT INTERNE A NT2
            IF(   CBPTIT2(1,K2).GE.-1D-14 .AND. CBPTIT2(1,K2).LE.1D0
     %      .AND. CBPTIT2(2,K2).GE.-1D-14 .AND. CBPTIT2(2,K2).LE.1D0
     %      .AND. CBPTIT2(3,K2).GE.-1D-14 .AND. CBPTIT2(3,K2).LE.1D0)
     %      THEN
C
C              POINT K2 INTERNE OU FRONTALIER AU TRIANGLE NT2
               IF(  CBPTIT2(1,K3).GE.-1D-14 .AND. CBPTIT2(1,K3).LE.1D0
     %         .AND.CBPTIT2(2,K3).GE.-1D-14 .AND. CBPTIT2(2,K3).LE.1D0
     %         .AND.CBPTIT2(3,K3).GE.-1D-14 .AND. CBPTIT2(3,K3).LE.1D0)
     %         THEN
C
C                 POINTS K2 K3 SONT INTERNES OU FRONTALIERS DU TRIANGLE NT2
C                 AJOUT DU SEGMENT PTI(1,K2) PTI(1,K3) DANS NT2 et NT1
C                 LES 2 POINTS D'INTERSECTION EXTREMITES DU SEGMENT
C                 PTI(K2)=STSGI(1,1)->PTI(K3)=STSGI(1,2)
C                 EXCEPTE SI PTI(K2) EST UN SOMMET DE NT2
C                         ET PTI(K3) EST UN SOMMET DE NT2 (DISTINCT OU NON)
                  IF( CBPTIT2(1,K2) .EQ. 1D0 .OR.
     %                CBPTIT2(2,K2) .EQ. 1D0 .OR.
     %                CBPTIT2(3,K2) .EQ. 1D0 ) THEN
                     IF( CBPTIT2(1,K3) .EQ. 1D0 .OR.
     %                   CBPTIT2(2,K3) .EQ. 1D0 .OR.
     %                   CBPTIT2(3,K3) .EQ. 1D0 ) THEN
C                       PTI(K2) EST UN SOMMET DE NT2 ET
C                       PTI(K3) EST UN SOMMET DISTINCT DE NT2
C                       L'ARETE DE NT2 NE SERA PAS A TRIANGULER
                        GOTO 100
                     ENDIF
                  ENDIF
C                 LE SGI EST PTI(K2)-PTI(K3) 2 PT INTERNES A NT2
                  STSGI(1,1) = PTI(1,K2)
                  STSGI(2,1) = PTI(2,K2)
                  STSGI(3,1) = PTI(3,K2)
                  STSGI(1,2) = PTI(1,K3)
                  STSGI(2,2) = PTI(2,K3)
                  STSGI(3,2) = PTI(3,K3)
                  GOTO 50
C
               ELSE
C
C                 POINT K2 INTERNE ET K3 EXTERNE A NT2
C                 RECHERCHE DES AU PLUS 3 POINTS D'INTERSECTION
C                 DU SEGMENT PTI(K2)-PTI(K3) AVEC LES 3 ARETES
C                 DU TRIANGLE NT2
                  CALL INTARTRPL( PTI(1,K2), PTI(1,K3),
     %                   XYZS2(1,NS1T2), XYZS2(1,NS2T2), XYZS2(1,NS3T2),
     %                   NBPTI, STSGT2 )
C                 NBPTI(K)=NB DE POINT INTERSECTION DE PTI(K2)-PTI(K3)
C                          ET DE L'ARETE K (EXTREMITES INCLUSES) DE NT2
C                         =0 LE POINT D'INTERSECTION EST EXTERNE A
C                            PTI(K2)-PTI(K3) ET/OU L'ARETE K DE P1P2P3
C                         =1 LE POINT D'INTERSECTION EST INTERNE A
C                            PTI(K2)-PTI(K3) ET A L'ARETE K DE NT2
C                            (SOMMETS EXTREMITES INCLUS)
                  NBPTIT = NBPTI(1) + NBPTI(2) + NBPTI(3)
                  IF( NBPTIT .EQ. 3 ) THEN
                     print *,'trtr3sgi nbptit=3 est un CAS IMPOSSIBLE!'
                     print *,'nt1=',nt1,' nt2=',nt2
                  ELSE IF( NBPTIT .EQ. 2 ) THEN
C                    LES 2 POINTS INTERNES D'INTERSECTION SONT
C                    UN SOMMET DE NT2 ET L'ARETE EST EXTERNE
C                    PAS DE SEGMENT D'INTERSECTION INTERNE A NT2
                     DO K=1,3
                        IF( NBPTI(K) .EQ. 1 ) GOTO 11
                     ENDDO
 11                  STSGI(1,1) = STSGT2(1,K)
                     STSGI(2,1) = STSGT2(2,K)
                     STSGI(3,1) = STSGT2(3,K)
                     DO K3=K+1,3
                        IF( NBPTI(K3) .EQ. 1 ) GOTO 12
                     ENDDO
 12                  STSGI(1,2) = STSGT2(1,K3)
                     STSGI(2,2) = STSGT2(2,K3)
                     STSGI(3,2) = STSGT2(3,K3)
                     GOTO 50
                  ELSE IF( NBPTIT .EQ. 1 ) THEN
C                    UNE SEULE ARETE D'INTERSECTION INTERNE PTI(K2)-STSGT2(?)
                     DO K=1,3
                        IF( NBPTI(K) .EQ. 1 ) GOTO 14
                     ENDDO
C                    LES 2 POINTS D'INTERSECTION EXTREMITES DU SEGMENT
C                    STSGI(1,1)=STSGT2(1,K) et PTI(K2)=STSGT2(1,2)
 14                  STSGI(1,1) = STSGT2(1,K)
                     STSGI(2,1) = STSGT2(2,K)
                     STSGI(3,1) = STSGT2(3,K)
                     STSGI(1,2) = PTI(1,K2)
                     STSGI(2,2) = PTI(2,K2)
                     STSGI(3,2) = PTI(3,K2)
                     GOTO 50
                  ELSE
C                    INTERSECTION NT1 NT2 REDUITE A UN POINT SUR UNE ARETE
C                    DE NT1 ET NT2 => PAS DE SGI
C                    TRACE DES 2 TRIANGLES NT1 et NT2 et
C                    des SEGMENTS INTERSECTION
                     IF( TRATRI ) THEN
               CALL TRASGI( NBS1,  XYZS1,  NBTRS1, NUSTS1, PSGISF1, NT1,
     %                      NBS2,  XYZS2,  NBTRS2, NUSTS2, PSGISF2, NT2,
     %                      MXSGI, NSTSGI, MXCHSGI,LCHSGI, MXPTA,XYZPTA,
     %                      3, PTI)
                     ENDIF
                     GOTO 100
                  ENDIF
C
               ENDIF
C
            ELSE
C
C              POINT K2 EXTERNE AU TRIANGLE NT2
               IF(   CBPTIT2(1,K3).GE.-1D-14 .AND. CBPTIT2(1,K3).LE.1D0
     %         .AND. CBPTIT2(2,K3).GE.-1D-14 .AND. CBPTIT2(2,K3).LE.1D0
     %         .AND. CBPTIT2(3,K3).GE.-1D-14 .AND. CBPTIT2(3,K3).LE.1D0)
     %         THEN
C
C                 POINT K2 EXTERNE K3 INTERNE OU FRONTALIER AU TRIANGLE NT2
C                 RECHERCHE DES AU PLUS 3 POINTS D'INTERSECTION
C                 DU SEGMENT PTI(K2)-PTI(K3) AVEC LES 3 ARETES
C                 DU TRIANGLE NT2
                  CALL INTARTRPL( PTI(1,K2), PTI(1,K3),
     %                   XYZS2(1,NS1T2), XYZS2(1,NS2T2), XYZS2(1,NS3T2),
     %                   NBPTI, STSGT2 )
C                 NBPTI(K)=NOMBRE DE POINT INTERNE INTERSECTION DE S1S2
C                          ET DE L'ARETE K (EXTREMITES INCLUSES)
C                         =0 LE POINT D'INTERSECTION EST EXTERNE A
C                            S1S2 ET/OU L'ARETE K DE P1P2P3
C                         =1 LE POINT D'INTERSECTION EST INTERNE A
C                            S1S2 ET A L'ARETE K DE NT2
C                            (SOMMETS EXTREMITES INCLUS)
                  NBPTIT = NBPTI(1) + NBPTI(2) + NBPTI(3)
                  IF( NBPTIT .EQ. 3 ) THEN
                     print *,'trtr3sgi nbptit=3 est un CAS IMPOSSIBLE!'
                     print *,'nt1=',nt1,' nt2=',nt2
                  ELSE IF( NBPTIT .EQ. 2 ) THEN
C                    LES 2 POINTS INTERNES D'INTERSECTION SONT
C                    UN SOMMET DE NT2 ET L'ARETE EST EXTERNE
C                    PAS DE SEGMENT D'INTERSECTION INTERNE A NT2
                     DO K=1,3
                        IF( NBPTI(K) .EQ. 1 ) GOTO 21
                     ENDDO
 21                  STSGI(1,1) = STSGT2(1,K)
                     STSGI(2,1) = STSGT2(2,K)
                     STSGI(3,1) = STSGT2(3,K)
                     DO K3=K+1,3
                        IF( NBPTI(K3) .EQ. 1 ) GOTO 22
                     ENDDO
 22                  STSGI(1,2) = STSGT2(1,K3)
                     STSGI(2,2) = STSGT2(2,K3)
                     STSGI(3,2) = STSGT2(3,K3)
                     GOTO 50
                  ELSE IF( NBPTIT .EQ. 1 ) THEN
C                    UNE SEULE ARETE D'INTERSECTION INTERNE PTI(K3)-STSGT2(?)
                     DO K=1,3
                        IF( NBPTI(K) .EQ. 1 ) GOTO 24
                     ENDDO
C                    LES 2 POINTS D'INTERSECTION EXTREMITES DU SEGMENT
C                    STSGI(1,1)=STSGT2(1,K) et PTI(K3)=STSGT2(1,2)
 24                  STSGI(1,1) = STSGT2(1,K)
                     STSGI(2,1) = STSGT2(2,K)
                     STSGI(3,1) = STSGT2(3,K)
                     STSGI(1,2) = PTI(1,K3)
                     STSGI(2,2) = PTI(2,K3)
                     STSGI(3,2) = PTI(3,K3)
                     GOTO 50
                  ELSE
C                    INTERSECTION NT1 NT2 REDUITE A UN POINT SUR UNE ARETE
C                    DE NT1 ET NT2 => PAS DE SGI
C                    TRACE DES 2 TRIANGLES NT1 et NT2 et
C                    des SEGMENTS INTERSECTION
                     IF( TRATRI ) THEN
               CALL TRASGI( NBS1,  XYZS1,  NBTRS1, NUSTS1, PSGISF1, NT1,
     %                      NBS2,  XYZS2,  NBTRS2, NUSTS2, PSGISF2, NT2,
     %                      MXSGI, NSTSGI, MXCHSGI,LCHSGI, MXPTA,XYZPTA,
     %                      3, PTI)
                     ENDIF
                     GOTO 100
                  ENDIF
               ELSE
C
C                 POINT K2 EXTERNE ET K3 EXTERNE AU TRIANGLE NT2
C                 RECHERCHE DES AU PLUS 3 POINTS D'INTERSECTION
C                 DU SEGMENT PTI(K2)-PTI(K3) AVEC LES 3 ARETES
C                 DU TRIANGLE NT2
                  CALL INTARTRPL( PTI(1,K2), PTI(1,K3),
     %                   XYZS2(1,NS1T2), XYZS2(1,NS2T2), XYZS2(1,NS3T2),
     %                   NBPTI, STSGT2 )
C                 NBPTI(K)=NOMBRE DE POINT INTERNE INTERSECTION DE S1S2
C                          ET DE L'ARETE K (EXTREMITES INCLUSES)
C                         =0 LE POINT D'INTERSECTION EST EXTERNE A
C                            S1S2 ET/OU L'ARETE K DE P1P2P3
C                         =1 LE POINT D'INTERSECTION EST INTERNE A
C                            S1S2 ET A L'ARETE K DE NT2 (SOMMETS EXTREMITES INCL
                  NBPTIT = NBPTI(1) + NBPTI(2) + NBPTI(3)
                  IF( NBPTIT .EQ. 2 ) THEN
C                    LES 2 POINTS INTERNES D'INTERSECTION SONT LES 2 SOMMETS
C                    DU SEGMENT D'INTERSECTION INTERNE A NT2
C                    STSGI(1,1)=STSGT2(1,K) et PTI(K3)=STSGT2(1,K+1ou2)
                     DO K=1,3
                        IF( NBPTI(K) .EQ. 1 ) GOTO 26
                     ENDDO
 26                  STSGI(1,1) = STSGT2(1,K)
                     STSGI(2,1) = STSGT2(2,K)
                     STSGI(3,1) = STSGT2(3,K)
                     DO KP=K+1,3
                        IF( NBPTI(KP) .EQ. 1 ) GOTO 28
                     ENDDO
 28                  STSGI(1,2) = STSGT2(1,KP)
                     STSGI(2,2) = STSGT2(2,KP)
                     STSGI(3,2) = STSGT2(3,KP)
                     GOTO 50
                  ELSE
C                    UN SEUL POINT INTERNE SOMMET DE NT2 OU 0 POINT INTERNE
C                    PAS D'ARETE D'INTERSECTION
                     GOTO 100
                  ENDIF
               ENDIF
            ENDIF
            print *,'trtr3sgi: NT1=',NT1,' NT2=',NT2,
     %        ' nbpti=',nbpti,' CAS NON TRAITE. A programmer...'
C
C           TRACE DES 2 TRIANGLES NT1 et NT2 et des SEGMENTS INTERSECTION
               CALL TRASGI( NBS1,  XYZS1,  NBTRS1, NUSTS1, PSGISF1, NT1,
     %                      NBS2,  XYZS2,  NBTRS2, NUSTS2, PSGISF2, NT2,
     %                      MXSGI, NSTSGI, MXCHSGI,LCHSGI, MXPTA,XYZPTA,
     %                      3, PTI)
            GOTO 100
C
C           AJOUT DU SEGMENT D'INTERSECTION STSGI(1,1)->STSGI(1,2)
C           DANS LE CHAINAGE DES SEGMENTS DES TRIANGLES NT1 ET NT2
C
C           IDENTIFICATION AUX SOMMETS DES SEGMENTS DEJA AJOUTES
 50         CALL  XYZIDEDS( STSGI(1,1), NBPTA, XYZPTA, NP1 )
            IF( NP1 .GT. 0 ) GOTO 60
C           AJOUT D'UN NOUVEAU PT-ST D'UN SGI
            IF( NBPTA .GE. MXPTA ) GOTO 9200
            NBPTA = NBPTA + 1
            XYZPTA(1,NBPTA) = STSGI(1,1)
            XYZPTA(2,NBPTA) = STSGI(2,1)
            XYZPTA(3,NBPTA) = STSGI(3,1)
            NP1 = NBPTA
C
 60         CALL  XYZIDEDS( STSGI(1,2), NBPTA, XYZPTA, NP2 )
            IF( NP2 .GT. 0 ) GOTO 70
C           AJOUT D'UN NOUVEAU PT-ST D'UN SGI
            IF( NBPTA .GE. MXPTA ) GOTO 9200
            NBPTA = NBPTA + 1
            XYZPTA(1,NBPTA) = STSGI(1,2)
            XYZPTA(2,NBPTA) = STSGI(2,2)
            XYZPTA(3,NBPTA) = STSGI(3,2)
            NP2 = NBPTA
C
 70         IF( NP1 .EQ. NP2 ) THEN
C              SGI REDUIT A UN SEUL POINT
               GOTO 100
            ENDIF
C
C           AJOUT DU SEGMENT D'INTERSECTION STSGI(1,1)->STSGI(1,2)
            IF( NBSGI .GE. MXSGI ) GOTO 9100
            NBSGI = NBSGI + 1
            NSTSGI(1,NBSGI) = NP1
            NSTSGI(2,NBSGI) = NP2
            NSTSGI(3,NBSGI) = NT1
            NSTSGI(4,NBSGI) = NT2
C
C           CHAINAGE DU SGI EN TETE DU CHAINAGE DE NT1
            IF( NBCHSGI .GE. MXCHSGI ) GOTO 9300
            NBCHSGI = NBCHSGI + 1
C           LE PREMIER SGI DU TRIANGLE NT1
            L1CH = PSGISF1( NT1 )
C           NUMERO DANS NSTSGI DU SGI  => LE NO DES 2 SOMMETS DANS PTA
            LCHSGI( 1, NBCHSGI ) = NBSGI
C           NUMERO DANS LCHSGI DU SGI SUIVANT
            LCHSGI( 2, NBCHSGI ) = L1CH
C           LE NOUVEAU SGI EN TETE DU CHAINAGE
            PSGISF1( NT1 ) = NBCHSGI
C
C           CHAINAGE DU SGI EN TETE DU CHAINAGE DE NT2
            IF( NBCHSGI .GE. MXCHSGI ) GOTO 9300
            NBCHSGI = NBCHSGI + 1
C           LE PREMIER SGI DU TRIANGLE NT2
            L1CH = PSGISF2( NT2 )
C           NUMERO DANS NSTSGI DU SGI  => LE NO DES 2 SOMMETS DANS PTA
            LCHSGI( 1, NBCHSGI ) = NBSGI
C           NUMERO DANS LCHSGI DU SGI SUIVANT
            LCHSGI( 2, NBCHSGI ) = L1CH
C           LE NOUVEAU SGI EN TETE DU CHAINAGE
            PSGISF2( NT2 ) = NBCHSGI
C
C           TRACE DES 2 TRIANGLES NT1 et NT2 et des SEGMENTS INTERSECTION
            IF( TRATRI ) THEN
               CALL TRASGI( NBS1,  XYZS1,  NBTRS1, NUSTS1, PSGISF1, NT1,
     %                      NBS2,  XYZS2,  NBTRS2, NUSTS2, PSGISF2, NT2,
     %                      MXSGI, NSTSGI, MXCHSGI,LCHSGI, MXPTA,XYZPTA,
     %                      0, PTI)
            ENDIF
C
C           PASSAGE AU TRIANGLE NT2 SUIVANT
 100     CONTINUE
C
 200  CONTINUE
C
C     TRACE DES SGI ET DES TRIANGLES DES 2 SURFACES
      TRATRI = .TRUE.
      IF( TRATRI ) THEN
         CALL TRASGIS12(  NBS1,  XYZS1,  NBTRS1,  NUSTS1, PSGISF1,
     %                    NBS2,  XYZS2,  NBTRS2,  NUSTS2, PSGISF2,
     %                    MXSGI, NSTSGI, MXCHSGI, LCHSGI, MXPTA, XYZPTA)
      ENDIF
      GOTO 9999
C
C     PAS ASSEZ DE SEGMENTS D'INTERSECTION DECLARABLES
 9100 NBLGRC(NRERR) = 1
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) ='TRTR3SGI: PAS ASSEZ DE SGI DECLARES. AUGMENTER mxsgi'
      ELSE
         KERR(1) ='TRTR3SGI: NOT ENOUGH SGI DECLARED. AUGMENT mxsgi'
      ENDIF
      CALL LEREUR
      IERR = 1
      GOTO 9999
C
C     PAS ASSEZ DE POINTS SOMMETS DES SGI DECLARABLES
 9200 NBLGRC(NRERR) = 2
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = 'TRTR3SGI: PAS ASSEZ DE SOMMETS DES SGI DECLARES'
         KERR(2) = 'AUGMENTER mxpta'
      ELSE
         KERR(1) = 'TRTR3SGI: NOT ENOUGH VERTICES OF SGI DECLARED'
         KERR(2) = 'AUGMENT mxpta'
      ENDIF
      CALL LEREUR
      IERR = 2
      GOTO 9999
C
C     PAS ASSEZ DE CHAINAGES DES SGI DECLARABLES
 9300 NBLGRC(NRERR) = 2
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = 'TRTR3SGI: PAS ASSEZ DE CHAINAGES DES SGI DECLARES'
         KERR(2) = 'AUGMENTER mxchsgi'
      ELSE
         KERR(1) = 'TRTR3SGI: NOT ENOUGH CHAINS of SGI DECLARED'
         KERR(2) = 'AUGMENT mxchsgi'
      ENDIF
      CALL LEREUR
      IERR = 2
      GOTO 9999
C
 9999 RETURN
      END
