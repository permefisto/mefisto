      SUBROUTINE XTN1TNBF( XYZPOI, NBEF,   NBNOEF, NUNOTE,
     %                     MOFACE, MXFACE, LFACES,
     %                     NDDLNO, NTDLHB, VXYZPNtn, VITMXtn,
     %                     tn,     tn1,  NTE1,  XYZ1,  V1,
     %                     DELTAe0,      NTE0,  XYZ0,  V0, IERR )
Cte ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCUL DU POINT X(tn;tn+1,XYZ1)=XYZ0 PAR INTEGRATION RETROGRADE
C ----- LE LONG DE LA CARACTERISTIQUE DANS DES TETRAEDRES BREZZI-FORTIN
C       ou LA VITESSE EST UN POLYNOME LAGRANGE DE DEGRE 1 +BULLE P4
C         LA PRESSION EST UN POLYNOME LAGRANGE DE DEGRE 1
C       LE TETRAEDRE ELEMENT FINI: e ref -> e EST P1 POUR CHAQUE COORDONNEE
C       cf Olivier PIRONNEAU

C ENTREES:
C --------
C XYZPOI : 3 COORDONNEES DES SOMMETS=POINTS DES TETRAEDRES
C NBEF   : NOMBRE DE TETRAEDRES DU MAILLAGE
C NBNOEF : NOMBRE DE NOEUDS D'UN TETRAEDRE
C NUNOTE : NUMERO DES NBNOEF NOEUDS DES TETRAEDRES

C MOFACE : NOMBRE DE MOTS DE CHAQUE FACE DU TABLEAU LFACES
C MXFACE : NOMBRE MAXIMAL de FACES DU TABLEAU LFACES
C LFACES : LFACES(1,I)= NO DU 1-ER  SOMMET DE LA FACE
C          LFACES(2,I)= NO DU 2-EME SOMMET > 1-ER  SOMMET
C          LFACES(3,I)= NO DU 3-EME SOMMET DE LA FACE
C          LFACES(4,I)= NO DU 4-EME SOMMET DE LA FACE
C                       0 SI TRIANGLE
C          LFACES(5,I)= CHAINAGE HACHAGE SUR FACE SUIVANTE
C          SI SOMMET 2 < DERNIER SOMMET  => FACE   DIRECTE DANS L'ELEMENT FINI
C                      >                 => FACE INDIRECTE
C          UNE FACE DIRECTE EST VUE DE L EXTERIEUR DE L'ELEMENT FINI
C          SOUS LA FORME DIRECTE
C          LFACES(6,I)= NUMERO DU 1-ER ELEMENT FINI CONTENANT CETTE FACE
C                       >0 SI FACE   DIRECTE DANS CET ELEMENT FINI
C                       <0 SI FACE INDIRECTE DANS CET ELEMENT FINI

C          SI LA FACE APPARTIENT A 2 ELEMENTS FINIS ALORS
C          LFACES(7,I)= NUMERO DU 2-EME ELEMENT FINI CONTENANT CETTE FACE
C                       >0 SI FACE   DIRECTE DANS CET ELEMENT FINI
C                       <0 SI FACE INDIRECTE DANS CET ELEMENT FINI
C          SINON
C          LFACES(7,I)= NUMERO DE LA FACE FRONTALIERE SUIVANTE
C                       0 SI C'EST LA DERNIERE
C          L1FAFR = NUMERO DE LA PREMIERE FACE FRONTALIERE DANS LFACES
C                   CHAINAGE SUIVANT DANS LFACES(7,I)
C          L1FA2M = NUMERO DE LA PREMIERE FACE INTERFACE 2 MATERIAUX DANS LFACES
C                   CHAINAGE SUIVANT DANS LFACES(8,I)
C          SI VOLUME MULTI-MATERIAUX
C          LFACES(8,I)= NO DE LA FACE INTERFACE ENTRE 2 MATERIAUX SUIVANTE
C                       0 SI FACE SANS INFORMATION SUPPLEMENTAIRE
C          SINON VOLUME MONO-MATERIAU
C          LFACES(8,I)= NUMERO DE LA FACE A TANGENTE
C                       0 SI FACE SANS TANGENTE

C NDDLNO : TABLEAU DES POINTEURS SUR LE DERNIER DL DE CHAQUE NOEUD FLUIDE
C          CE TABLEAU EST DIMENSIONNE 0:NBNOEU
C NTDLHB : NOMBRE TOTAL DE DEGRES DE LIBERTE DE L'OBJET HORS BARYCENTRES
C          (NDIM+1) NBNOHB (SOMMETS SANS les BARYCENTRES)
C VXYZPNtn: VECTEUR GLOBAL des DL VITESSES-PRESSIONS (1:NTDLHB) au TEMPS tn
C          SUIVI DES 3 VITESSES AU BARYCENTRE DES TETRAEDRES
C VITMXtn: NORME DE LA VITESSE MAXIMALE AU TEMPS tn

C tn     : TEMPS INITIAL DE L'INTERVALLE de TEMPS
C tn1    : TEMPS FINAL   DE L'INTERVALLE de TEMPS
C NTE1   : NUMERO du TETRAEDRE CONTENANT LE POINT XYZ1 au TEMPS tn+1
C XYZ1   : POINT XYZ1( tn+1 ) INITIAL DE LA CARACTERISTIQUE RETROGRADE

C SORTIES:
C --------
C V1     : VITESSE( tn+1, XYZ1 ) selon l'INTERPOLATION Brezzi-FORTIN

C DELTAe0: 6*VOLUME du TETRAEDRE NTE0 CONTENANT XYZ0 de VITESSE V0
C NTE0   : >0 NUMERO DU TETRAEDRE CONTENANT XYZ0=X(tn;tn+1,XYZ1)
C             CORRECTEMENT CALCULE ou ARRETE AU DERNIER TETRAEDRE PARCOURU
C XYZ0   : POINT XYZ0=X(tn;tn+1,XYZ1) FINAL DE LA CARACTERISTIQUE RETROGRADE
C V0     : VITESSE( tn, XYZ0 ) selon l'INTERPOLATION Brezzi-FORTIN

C IERR   : =0 V1 XYZ0 V0 CALCULES CORRECTEMENT ou
C             V1 VITESSE TRES FAIBLE  XYZ0=XYZ1 et V0=V1
C          =1 VITESSE SUR NTE1 TROP GRANDE > VITMXtn  XYZ0=XYZ1 et V0=V1
C          =2 XYZ0 EST EXTERNE aux TETRAEDRES
C          =3 DELTAe<=0 ou XYZ1 NON INTERNE A NTE1
C                       ou XYZ0 NON INTERNE A NTE0
C             -> XYZ0=XYZ1 et V0=V1=0
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Saint Pierre du Perray           Octobre 2020
C23456---------------------------------------------------------------012
      IMPLICIT  NONE
      include"./incl/langue.inc"
      REAL               tn, tn1, DTNTN1, XYZPOI(3,*), LONARE(6)

      INTEGER            NTE0, NTE1, NTE, NTES,
     %                   NBEF, NBNOEF, NUNOTE(NBEF,NBNOEF),
     %                   MOFACE, MXFACE, LFACES(MOFACE,MXFACE),
     %                   NTDLHB, NDDLNO(0:*), IERR,
     %                   I, K, NDL, NS, NONOUI

      DOUBLE PRECISION   VXYZPNtn(1:*), VITMXtn, VMXNTE1, ARETMXNTE1,
     %                   TEDDT, DDT, XYZ0(3), XYZ1(3), XYZ(3),
     %                   D, PRCBL1(5), PRCBL0(5), PRCBLT,
     %                   V1NORM, V1(3), V0(3), V(3),
     %                   DELTAe0, DELTAe, BULLE, CBTE(4), CB1234
ccc     %                  ,VNORM

      INTEGER            NOSOFATE(3,4)
      DATA               NOSOFATE / 1,3,2,  2,3,4,  3,1,4,  4,1,2 /

      INTEGER            NOSOARTE(2,6)
      DATA               NOSOARTE / 1,2,  2,3,  3,1,  1,4,  2,4,  3,4 /

      INTEGER            NOAR2FAC(4,4)
      DATA               NOAR2FAC / 0, 2, 3, 1,
     %                              2, 0, 6, 5,
     %                              3, 6, 0, 4,
     %                              1, 5, 4, 0 /

      INTRINSIC          SQRT, DBLE

C     VERIFICATION QUE XYZ1 EST INTERNE AU TETRAEDRE NTE1
      CALL XYZD1TE( NBEF, NBNOEF, NUNOTE, XYZPOI,
     %              XYZ1, NTE1,   DELTAe, CBTE,   NONOUI )
C     NONOUI= 1 OUI XYZ   EST     INTERNE AU TETRAEDRE NTE  CBTE CALCULE
C           = 0 NON XYZ N'EST PAS INTERNE AU TETRAEDRE NTE  CBTE CALCULE
C           =-1 TETRAEDRE NTE DE VOLUME<=0  CBTE NON CALCULE

      IF( NONOUI .EQ. -1 ) THEN
C        SORTIE avec XYZ0=XYZ1 V1=V0=0
         DELTAe0 = DELTAe
         GOTO 9900
      ENDIF

      CB1234 = ABS(CBTE(1)) + ABS(CBTE(2)) + ABS(CBTE(3)) + ABS(CBTE(4))
      IF( CB1234 .GT. 1D0 ) THEN
         CBTE(1) = CBTE(1) / CB1234
         CBTE(2) = CBTE(2) / CB1234
         CBTE(3) = CBTE(3) / CB1234
         CBTE(4) = CBTE(4) / CB1234
      ENDIF

C     INTERPOLATION de la VITESSE V1( tn, XYZ1 ) dans le TETRAEDRE Brezzi-Fortin
C     CALCUL DES POLYNOMES LAGRANGE DE DEGRE 1 + BULLE P4 sur TETRAEDRE UNITE
      BULLE = 64D0 * CBTE(1) * CBTE(2) * CBTE(3) * CBTE(4)
      PRCBL1(1) = CBTE(1) - BULLE
      PRCBL1(2) = CBTE(2) - BULLE
      PRCBL1(3) = CBTE(3) - BULLE
      PRCBL1(4) = CBTE(4) - BULLE
      PRCBL1(5) = 4D0 * BULLE
      PRCBLT = PRCBL1(1) + PRCBL1(2) + PRCBL1(3) + PRCBL1(4) + PRCBL1(5)

      IF( PRCBLT .GT. 1D0 ) THEN
ccc         PRINT*,'xtn1tnbf: NTE1=',NTE1,' PRCBL1=',PRCBL1,
ccc     %          ' PRCBLT1=',PRCBLT
         PRCBL1(1) = PRCBL1(1) / PRCBLT
         PRCBL1(2) = PRCBL1(2) / PRCBLT
         PRCBL1(3) = PRCBL1(3) / PRCBLT
         PRCBL1(4) = PRCBL1(4) / PRCBLT
         PRCBL1(5) = PRCBL1(5) / PRCBLT
      ENDIF

C     CALCUL DE LA VITESSE V1( tn, (CBTE(2),CBTE(3),CBTE(4)) ) dans l'EF NTE1
      V1(1) = 0D0
      V1(2) = 0D0
      V1(3) = 0D0
      DO I=1,4
C        NO DU DL I DANS L'EF NTE1 DE LA VITESSE AU TEMPS tn
C        AU POINT X(tn;tn+1,XYZ1)
         NDL  = NDDLNO( NUNOTE(NTE1,I) - 1 )
         V1(1) = V1(1) + PRCBL1(I) * VXYZPNtn( NDL+1 )
         V1(2) = V1(2) + PRCBL1(I) * VXYZPNtn( NDL+2 )
         V1(3) = V1(3) + PRCBL1(I) * VXYZPNtn( NDL+3 )
      ENDDO
C     LE NO DU DERNIER DL VITESSE DU BARYCENTRE DU TETRAEDRE NTE1
      NDL  = NTDLHB + 3 * NTE1
      V1(1) = V1(1) + PRCBL1(5) * VXYZPNtn( NDL-2 )
      V1(2) = V1(2) + PRCBL1(5) * VXYZPNtn( NDL-1 )
      V1(3) = V1(3) + PRCBL1(5) * VXYZPNtn( NDL   )
C     NORME DE LA VITESSE V1
      V1NORM = SQRT( V1(1)**2 + V1(2)**2 + V1(3)**2 )

      IF( V1NORM .GT. 1.001D0 * VITMXtn  ) THEN
ccc         print*,'xtn1tnbf: EF1=',NTE1,' VITMXtn=',VITMXtn,
ccc     %             ' < V1NORM=',V1NORM,' TROP GRAND!... V1=',V1
ccc         print*,'xtn1tnbf: CBTE=',CBTE
ccc         print*,'xtn1tnbf: PRCBL1(1:5)=',PRCBL1,' PRCBLT=',PRCBLT
ccc         print*,'xtn1tnbf: XYZ1=', XYZ1
ccc         print*
         D = V1NORM / VITMXtn
         V1(1) = V1(1) / D
         V1(2) = V1(2) / D
         V1(3) = V1(3) / D
      ENDIF

C     VITESSE MAXIMALE AUX 4 SOMMETS DU TETRAEDRE NTE1
      VMXNTE1 = 0D0
      DO I=1,4
C        NUMERO DU NOEUD I DU TETRAEDRE NTE1
         NS = NUNOTE(NTE1,I)
C        NO DU DL I DANS L'EF NTE1 DE LA VITESSE
         NDL = NDDLNO( NS - 1 )
C        CARRE DE LA VITESSE AU NOEUD I
         D= VXYZPNtn(NDL+1)**2 + VXYZPNtn(NDL+2)**2 + VXYZPNtn(NDL+3)**2
         IF( D .GT. VMXNTE1 ) VMXNTE1 = D
      ENDDO

C     VITESSE AU BARYCENTRE DE L'EF NTE1
      NDL = NTDLHB + 3 * NTE1
C     CARRE DE LA VITESSE AU BARYCENTRE DE L'EF NTE1
      D = VXYZPNtn(NDL-2)**2 + VXYZPNtn(NDL-1)**2 + VXYZPNtn(NDL)**2
      IF( D .GT. VMXNTE1 ) VMXNTE1 = D

C     VITESSE MAXIMALE AUX 5 NOEUDS DU TETRAEDRE NTE1 Brezzi-Fortin
      VMXNTE1 = SQRT( VMXNTE1 )

ccc      IF( VMXNTE1 .GT. 1.1D0 * VITMXtn ) THEN
ccc         PRINT*,'xtn1tnbf: NTE1=',NTE1,' VMXNTE1=',VMXNTE1,
ccc     %          ' > VITMXtn=',VITMXtn,' ?...'
ccc      ENDIF

      IF( VMXNTE1 .LE. 1D-3 * VITMXtn ) THEN
C        VITESSE TRES FAIBLE
         GOTO 9910
      ENDIF

C     CALCUL DE LA LONGUEUR DES 6 ARETES ET DE L'ARETE MAXIMALE DE NTE1
      CALL LON6AR( XYZPOI( 1, NUNOTE(NTE1,1) ),
     %             XYZPOI( 1, NUNOTE(NTE1,2) ),
     %             XYZPOI( 1, NUNOTE(NTE1,3) ),
     %             XYZPOI( 1, NUNOTE(NTE1,4) ),  LONARE )
      ARETMXNTE1 = DBLE( MAX( LONARE(1), LONARE(2), LONARE(3),
     %                        LONARE(4), LONARE(5), LONARE(6) ) )

C     Remontee RETROGRADE de la CARACTERISTIQUE de tn+1 a tn: c-a-d
C     CALCUL du POINT XYZ0=X(tn;tn+1,XYZ1) au temps tn
C     QUI SERA au temps tn+1 AU POINT XYZ1.
C     En fin d'integration retrograde:
C     XYZ0 = Fe( CBTE(2),CBTE(3),CBTE(4) ) = X(tn;tn+1,XYZ1)
C     =============================================================

C     INITIALISATIONS du PAS de TEMPS DDT pour PASSER de tn+1 a tn:
C     -------------------------------------------------------------
C     AJUSTEMENT DU SOUS-PAS DE TEMPS de tn+1 a tn
C     PAS DE TEMPS ACTUEL
      DTNTN1 = tn1 - tn

C     16 SOUS PAS DE TEMPS DDT ENTRE tn1 et tn POUR PARCOURIR ARETMXNTE1
C     A LA VITESSE VITMXtn
      DDT = MIN( ARETMXNTE1 / ( 16 * VITMXtn ),
     %           ARETMXNTE1 / (  6 * VMXNTE1 ), DTNTN1/4 )
      IF( DDT .LE. 0D0 ) GOTO 9920


C     DEPART dans NTE1: COORDONNEES du POINT INITIAL XYZ0 de PARCOURS
      XYZ0(1) = XYZ1(1)
      XYZ0(2) = XYZ1(2)
      XYZ0(3) = XYZ1(3)

C     DEPART dans NTE1: La VITESSE V0=V1 en XYZ0=XYZ1 au TEMPS tn+1
      V0(1) = V1(1)
      V0(2) = V1(2)
      V0(3) = V1(3)

C     No DU 1-er TETRAEDRE de RECHERCHE DEVANT CONTENIR XYZ0
      NTE0 = NTE1
      DELTAe0 = DELTAe

C     TEMPS D'INTEGRATION RETROGRADE de tn1 a tn EFFECTUE
      TEDDT = 0D0


C     LA BOUCLE SUR LES PAS DE TEMPS -DDT de tn+1 a tn:
C     =================================================
 10   IF( TEDDT .GE. DTNTN1 ) GOTO 100

      TEDDT = TEDDT + DDT

C     LE NOUVEAU POINT XYZ INTEGRE RETROGRADE EN -DDT DE XYZ0
C     -------------------------------------------------------
      XYZ(1) = XYZ0(1) - DDT * V0(1)
      XYZ(2) = XYZ0(2) - DDT * V0(2)
      XYZ(3) = XYZ0(3) - DDT * V0(3)

C     RECHERCHE A PARTIR DE XYZ0 DU TETRAEDRE NTE0 DU TETRAEDRE NTE
C     CONTENANT LE NOUVEAU POINT RETROGRADE XYZ
C     -------------------------------------------------------------
      CALL XYZ1VXYZ( XYZPOI, NBEF,   NBNOEF, NUNOTE,
     %               MOFACE, MXFACE, LFACES,
     %               XYZ0, NTE0, XYZ,  NTES )
C     NTES: >0  NUMERO DU TETRAEDRE CONTENANT XYZ INTERNE
C           <0 -NUMERO DU TETRAEDRE CONTENANT XYZ SUR UNE FACE FRONTALIERE
C           =0  PAS DE TETRAEDRE CONTENANT XYZ

      NTE = ABS( NTES )
      IF( NTE .EQ. 0 ) THEN
C        PAS DE TETRAEDRE CONTENANT XYZ
         IERR = 2
         GOTO 9990
      ENDIF

C     CALCUL DE DELTAe, CBTE de XYZ
      CALL XYZD1TE( NBEF, NBNOEF, NUNOTE, XYZPOI,
     %              XYZ,  NTE,    DELTAe, CBTE,   NONOUI )
C     NONOUI= 1 OUI XYZ   EST     INTERNE AU TETRAEDRE NTE  CBTE CALCULE
C           = 0 NON XYZ N'EST PAS INTERNE AU TETRAEDRE NTE  CBTE CALCULE
C           =-1 TETRAEDRE NTE DE VOLUME<=0  CBTE NON CALCULE

      IF( NONOUI .EQ. -1 ) THEN
C        NTE TETRAEDRE DE VOLUME<=0
C        SORTIE avec XYZ0=XYZ1 V1=V0=0
         GOTO 9900
      ENDIF

C     XYZ EST INTERNE OU FRONTIERE DU TETRAEDRE NTE
      DO I=1,4
         IF( CBTE(I) .LT. 0D0 ) THEN
C           POUR EVITER UNE EXTRAPOLATION AU DELA DU TETRAEDRE
            CBTE(I) = 0D0
         ENDIF
      ENDDO

      CB1234 = ABS(CBTE(1)) + ABS(CBTE(2))
     %       + ABS(CBTE(3)) + ABS(CBTE(4))
      IF( CB1234 .GT. 1D0 ) THEN
         CBTE(1) = CBTE(1) / CB1234
         CBTE(2) = CBTE(2) / CB1234
         CBTE(3) = CBTE(3) / CB1234
         CBTE(4) = CBTE(4) / CB1234
      ENDIF

C     CALCUL DE LA VITESSE V( XYZ ) DANS NTE AVEC LES CBTE
      BULLE = 64D0 * CBTE(1) * CBTE(2) * CBTE(3) * CBTE(4)
      PRCBL0(1) = CBTE(1) - BULLE
      PRCBL0(2) = CBTE(2) - BULLE
      PRCBL0(3) = CBTE(3) - BULLE
      PRCBL0(4) = CBTE(4) - BULLE
      PRCBL0(5) = 4D0 * BULLE
      PRCBLT = PRCBL0(1) + PRCBL0(2) + PRCBL0(3) + PRCBL0(4) + PRCBL0(5)

ccc      IF( PRCBLT.GT.1.015D0) THEN
ccc         CB1234 = ABS(CBTE(1)) + ABS(CBTE(2))
ccc     %          + ABS(CBTE(3)) + ABS(CBTE(4))
ccc       PRINT*,'xtn1tnbf: NTE=',NTE,' CBTE  =',CBTE  ,' CBTE1234=',CB1234
ccc       PRINT*,'xtn1tnbf: NTE=',NTE,' PRCBL0=',PRCBL0,' PRCBLT0 =',PRCBLT
ccc      ENDIF

      IF( PRCBLT .GT. 1D0 ) THEN
         PRCBL0(1) = PRCBL0(1) / PRCBLT
         PRCBL0(2) = PRCBL0(2) / PRCBLT
         PRCBL0(3) = PRCBL0(3) / PRCBLT
         PRCBL0(4) = PRCBL0(4) / PRCBLT
         PRCBL0(5) = PRCBL0(5) / PRCBLT
      ENDIF

C     CALCUL DE LA VITESSE V(tn,(CBTE(2),CBTE(3),CBTE(4))) dans l'EF NTE
      V(1) = 0D0
      V(2) = 0D0
      V(3) = 0D0
      DO I=1,4
C        NO DU DL I DANS L'EF NTE DE LA VITESSE V
         NDL  = NDDLNO( NUNOTE(NTE,I) - 1 )
         V(1) = V(1) + PRCBL0(I) * VXYZPNtn( NDL+1 )
         V(2) = V(2) + PRCBL0(I) * VXYZPNtn( NDL+2 )
         V(3) = V(3) + PRCBL0(I) * VXYZPNtn( NDL+3 )
      ENDDO
C     CONTRIBUTION DE LA VITESSE AU BARYCENTRE DU TETRAEDRE NTE
      NDL  = NTDLHB + 3 * NTE
      V(1) = V(1) + PRCBL0(5) * VXYZPNtn( NDL-2 )
      V(2) = V(2) + PRCBL0(5) * VXYZPNtn( NDL-1 )
      V(3) = V(3) + PRCBL0(5) * VXYZPNtn( NDL   )

cccC     NORME DE LA VITESSE V en XYZ
ccc      VNORM = SQRT( V(1)**2 + V(2)**2 + V(3)**2 )
ccc      IF( VNORM .GT. 1.1D0 * VITMXtn ) THEN
ccc         print*,'xtn1tnbf: EF1=',NTE1,' A VOIR NTE0=',NTE0,' V=',V,
ccc     %          ' VNORM=',VNORM,' > VITMXtn=',VITMXtn,' trop grand!...'
ccc         print*,'xtn1tnbf: CBTE=',CBTE
ccc         print*,'xtn1tnbf: PRCBL0(1:5)=',PRCBL0,' PRCBLT=',PRCBLT
ccc         print*
ccc      ENDIF

C     XYZ V NTE SONT ACCEPTES  et  DEVIENNENT XYZ0 V0 ET NTE0
C     -------------------------------------------------------
      NTE0 = NTE
      DELTAe0 = DELTAe

      XYZ0( 1 ) = XYZ( 1 )
      XYZ0( 2 ) = XYZ( 2 )
      XYZ0( 3 ) = XYZ( 3 )

      V0( 1 ) = V( 1 )
      V0( 2 ) = V( 2 )
      V0( 3 ) = V( 3 )

      IF( NTES .LT. 0 ) THEN
C        XYZ0 EST SUR UNE FACE FRONTIERE DE NTE0.
C        L'INTEGRATION RETROGRADE EN TEMPS PAR PAS DDT S'ARRETE LA
         GOTO 100
      ENDIF

C     PASSAGE A L'INTERVALLE de TEMPS -DDT SUIVANT entre tn et tn+1
      GOTO 10


C     en FIN D'INTEGRATION entre tn+1 et tn:
C     ici XYZ0=X(tn;tn+1,XYZ1) et V0 = Vitesse( X(tn;tn+1,XYZ1) )
C     SONT CALCULES dans NTE0 avec DELTAe0
C     ===========================================================
 100  IERR = 0
      GOTO 9990


C     DELTAe<=0 ou XYZ1 NON INTERNE A NTE1 ou XYZ0 NON INTERNE A NTE0
 9900 DO K=1,3
         XYZ0(K) = XYZ1(K)
         V1(K)   = 0D0
         V0(K)   = 0D0
      ENDDO
      IERR = 3
      print*,'xtn1tnbf: NTE1=',NTE1,' XYZ1=',XYZ1,' DELTAe=',DELTAe,
     %       ' En SORTIE: XYZ0=XYZ1 V1=V0=0 IERR=3'
      GOTO 9999


C     V1 VITESSE TRES FAIBLE: XYZ0=XYZ1
 9910 DO K=1,3
         XYZ0(K) = XYZ1(K)
         V0(K)   = V1(K)
      ENDDO
      IERR = 0
ccc      print*,'xtn1tnbf: NTE1=',NTE1,' XYZ1=',XYZ1,' FAIBLE V1=',V1,
ccc     %       ' En SORTIE XYZ0=XYZ1 V0=V1 IERR=0'
      GOTO 9990


C     DDT<=0
 9920 DO K=1,3
         XYZ0(K) = XYZ1(K)
         V0(K)   = V1(K)
      ENDDO
      IERR = 1
      print*,'xtn1tnbf: NTE1=',NTE1,' XYZ1=',XYZ1,' DDT=',DDT,' V1=',V1,
     %       ' En SORTIE XYZ0=XYZ1 V0=V1 IERR=1'


C     PROJECTION EVENTUELLE SUR LA VITESSE MAX ANTERIEURE
 9990 D = SQRT( V0(1)**2 + V0(2)**2 + V0(3)**2 ) / VITMXtn
      IF( D .GT. 1D0 ) THEN
         V0(1) = V0(1) / D
         V0(2) = V0(2) / D
         V0(3) = V0(3) / D
      ENDIF

 9999 RETURN
      END
