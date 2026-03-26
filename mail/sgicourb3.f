      SUBROUTINE SGICOURB3( NT1,    NUSTS1,  XYZS1,  XYZPTA,
     %                      LCHSGI, NSTSGI,
     %                      NBCB,   L1COURB, NUSUCB, NUARCB )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    REORDONNER LES COURBES SELON LE SENS DIRECT DES 3 ARETES DU
C -----    TRIANGLE NT1
C
C ENTREES:
C --------
C NT1    : NUMERO DU TRIANGLE
C NUSTS1 : (4,NBTRS1) NUMERO DANS XYZS1 DES 3 SOMMETS + 0 EN POSITION 4
C XYZS1  : 3 XYZ DES NBS1 SOMMETS DE LA SURFACE 1
C XYZPTA : 3 XYZ DES SOMMETS AJOUTES PAR LES SGI
C NBCB   : NOMBRE FINAL DE COURBES DE SGI
C
C MODIFIES:
C----------
C LCHSGI : LCHSGI(1,.) NUMERO DU SGI DANS NSTSGI
C          LCHSGI(2,.) CHAINAGE SUR LE SGI SUIVANT DANS LCHSGI
C NSTSGI : NSTSGI(1,.) NUMERO DANS XYZPTA DU PREMIER POINT DU SGI
C          NSTSGI(2,.) NUMERO DANS XYZPTA DU SECOND  POINT DU SGI
C          NSTSGI(3,.) NUMERO DU TRIANGLE DANS LA TRIANGULATION NUSTS1
C          NSTSGI(4,.) NUMERO DU TRIANGLE DANS LA TRIANGULATION NUSTS2
C
C SORTIES:
C --------
C L1COURB: L1COURB(k) NUMERO DANS LCHSGI DU PREMIER SGI DE LA COURBE k
C          DE SOMMET INITIAL NUSUCB(1,k) et FINAL NUSUCB(2,k)
C NUSUCB : NUSUCB(1:2,k) NO DU PREMIER ET DERNIER SOMMET DANS XYZPTA
C          DE LA COURBE k
C NUARCB : NUARCB(L,M) NUMERO DE L'ARETE NT1 DU SU L DE LA COURBE M
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET Alain LJLL UPMC & St Pierre du Perray Novembre 2011
C2345X7..............................................................012
      INTEGER  NUSTS1(4,*), LCHSGI(2,*), NSTSGI(4,*),
     %         L1COURB(0:NBCB), NUSUCB(2,0:NBCB), NUARCB(2,0:NBCB)
      DOUBLE PRECISION  XYZS1(3,*), XYZPTA(3,*),
     %                  CBPTTR1(3), CBPTTR0(3), CB0, CB
C
C     NUMERO DES 3 SOMMETS DU TRIANGLE NT1 CONTENANT AU MOINS 1 SGI
      NS1T1 = NUSTS1(1,NT1)
      NS2T1 = NUSTS1(2,NT1)
      NS3T1 = NUSTS1(3,NT1)
C
C     PARCOURS DES NBCB COURBES ACTUELLES
      DO 20 NUCB=1,NBCB
C
C        INVERSION DE LA COURBE NUCB SI LE NUMERO D'ARETE DES SU
C        EXTREMITES N'EST PAS CROISSANT
C        -------------------------------------------------------
         IF( NUARCB(1,NUCB) .GT. NUARCB(2,NUCB) ) THEN
C
C           PERMUTATION DES EXTREMITES ET SGI DE LA COURBE NUCB
            N              = NUSUCB(1,NUCB)
            NUSUCB(1,NUCB) = NUSUCB(2,NUCB)
            NUSUCB(2,NUCB) = N
C
            N              = NUARCB(1,NUCB)
            NUARCB(1,NUCB) = NUARCB(2,NUCB)
            NUARCB(2,NUCB) = N
C
C           INVERSION DES SGI DE LA COURBE NUCB
C           LA PREMIERE ARETE EST LE DERNIER DES SGI
            LCHINI = L1COURB(NUCB)
            LCHDER = LCHINI
C
C           RECHERCHE DU DERNIER SGI DE LA COURBE
 5          LCH = LCHSGI( 2, LCHDER )
            IF( LCH .GT. 0 ) THEN
               LCHDER = LCH
               GOTO 5
            ENDIF
C           LCHDER EST LE DERNIER SGI QUI DEVIENT LE PREMIER DE LA COURBE
            L1COURB(NUCB) = LCHDER
C           LE NUMERO DU SU DE DEPART DE LA COURBE
            NSG0 = NUSUCB(1,NUCB)
C
C           INVERSION DU DERNIER SGI
C           NUMERO DANS LCHSGI DU SGI
 8          NSGI = LCHSGI( 1, LCHDER )
C           NSG1 NUMERO DU SOMMET 1 DU SGI DANS XYZPTA
            NSG1 = NSTSGI( 1, NSGI )
C           NSG2 NUMERO DU SOMMET 2 DU SGI DANS XYZPTA
            NSG2 = NSTSGI( 2, NSGI )
C
            IF( NSG2 .EQ. NSG0 ) THEN
C              INVERSION DES 2 SOMMETS DU SGI
               NSTSGI( 1, NSGI ) = NSG2
               NSTSGI( 2, NSGI ) = NSG1
            ENDIF
C
C           LE DERNIER SOMMET DU SGI ACTUEL
            NSG0 = NSTSGI( 2, NSGI )
C
C           RECHERCHE DU PRECEDENT DU NOUVEAU DERNIER SGI DE LA COURBE
            IF( LCHDER .NE. LCHINI ) THEN
               LCH0 = LCHINI
 10            LCH = LCHSGI( 2, LCH0 )
               IF( LCH .NE. LCHDER ) THEN
                   LCH0 = LCH
                   GOTO 10
               ENDIF
C
C              LCH0 PRECEDE LCHDER
               LCHSGI( 2, LCHDER ) = LCH0
C
C              CE SGI DEVIENT LE PROCHAIN DERNIER A TRAITER
               LCHDER = LCH0
               GOTO 8
C
            ENDIF
C
C           CHAINAGE FINAL DE LA COURBE NUCB
            LCHSGI( 2, LCHDER ) = 0
C
         ENDIF
C
 20   CONTINUE
C
C     REORDONNER LES COURBES POUR LE DEPART CROISSANT SELON
C     LE SENS DIRECT DES ARETES 1 2 3 DE NT1
C     -----------------------------------------------------
C     PARCOURS DES ARETES DE NT1
      DO 80 NUAR= 1, 3
C
C        NUMERO DE LA COORDONNEE BARYCENTRIQUE SUR L'ARETE
         IF( NUAR .LT. 3 ) THEN
            NUCOBA = NUAR+1
         ELSE
            NUCOBA = 1
         ENDIF
C
C        PARCOURS DES NBCB COURBES ACTUELLES
         NUCB0 = 0
         CB0   = 0D0
C
         DO 50 NUCB=1,NBCB
C
C           COURBE NUCB A TRAITER
            IF( NUARCB(1,NUCB) .NE. NUAR ) GOTO 50
C
C           CALCUL DE LA COORDONNEE BARYCENTRIQUE SUR L'ARETE NUAR
            CALL CBPTTR( XYZS1(1,NS1T1), XYZS1(1,NS2T1),
     %                   XYZS1(1,NS3T1), XYZPTA(1,NUSUCB(1,NUCB)),
     %                   CBPTTR1 )
            CB = CBPTTR1(NUCOBA)
C
            IF( CB .GT. CB0 ) GOTO 40
C
            IF( CB .EQ. CB0 ) THEN
C
C              ORDRE SELON LE NUMERO D'ARETE FINALE DECROISSANT
               NUA  = NUARCB(2,NUCB)
               NUA0 = NUARCB(2,NUCB0)
               IF( NUA0 .GT. NUA ) GOTO 40
C
               IF( NUA0 .EQ. NUA ) THEN
C
C                 CALCUL DES COORDONNEES BARYCENTRIQUES DE NUSUCB(2,NUCB0)
                  CALL CBPTTR( XYZS1(1,NS1T1),XYZS1(1,NS2T1),
     %                         XYZS1(1,NS3T1),XYZPTA(1,NUSUCB(2,NUCB0)),
     %                         CBPTTR0 )
C
C                 CALCUL DES COORDONNEES BARYCENTRIQUES DE NUSUCB(2,NUCB)
                  CALL CBPTTR( XYZS1(1,NS1T1), XYZS1(1,NS2T1),
     %                         XYZS1(1,NS3T1), XYZPTA(1,NUSUCB(2,NUCB)),
     %                         CBPTTR1 )
C
C                 NO DE LA COORDONNEE BARYCENTRIQUE DE L'ARETE NUA
                  IF( NUA .LT. 3 ) THEN
                     N = NUA+1
                  ELSE
                     N = 1
                  ENDIF
C                 COMPARAISON DE LA COORDONNEE BARYCENTRIQUE SUR L'ARETE NUA
                  IF( CBPTTR0(N) .GT. CBPTTR1(N) ) GOTO 40
C                 ORDRE SELON LA COORDONNEE BARYCENTRIQUE DECROISSANTE
C
               ENDIF
C
            ENDIF
C
C           PERMUTATION DES COURBES NUCB ET NUCB0
            N                = L1COURB( NUCB0 )
            L1COURB( NUCB0 ) = L1COURB( NUCB  )
            L1COURB( NUCB  ) = N
C
            N               = NUSUCB(1,NUCB0)
            NUSUCB(1,NUCB0) = NUSUCB(1,NUCB )
            NUSUCB(1,NUCB ) = N
C
            N               = NUSUCB(2,NUCB0)
            NUSUCB(2,NUCB0) = NUSUCB(2,NUCB )
            NUSUCB(2,NUCB ) = N
C
            N               = NUARCB(1,NUCB0)
            NUARCB(1,NUCB0) = NUARCB(1,NUCB )
            NUARCB(1,NUCB ) = N
C
            N               = NUARCB(2,NUCB0)
            NUARCB(2,NUCB0) = NUARCB(2,NUCB )
            NUARCB(2,NUCB ) = N
C
C           DERNIERE COURBE TRAITEE
 40         CB0   = CB
            NUCB0 = NUCB
C
 50      CONTINUE
C
 80   CONTINUE
C
cccC     AFFICHAGE DES NBCB COURBES DE SGI
ccc      print *
ccc      print *, ('L1COURB(',M,')=',L1COURB(M),M=1,NBCB)
ccc      DO M=1,NBCB
ccc         print *,'SGICOURB3:  NUSUCB(1,',M,')=',NUSUCB(1,M),
ccc     %                      ' NUSUCB(2,',M,')=',NUSUCB(2,M)
ccc         print *,'SGICOURB3:  NUARCB(1,',M,')=',NUARCB(1,M),
ccc     %                      ' NUARCB(2,',M,')=',NUARCB(2,M)
ccc      ENDDO
ccc      CALL SGIIMPR( L1COURB(1), LCHSGI, NSTSGI )
C
      RETURN
      END
