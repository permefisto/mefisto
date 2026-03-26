      SUBROUTINE SUEX28( NTLXSU , LADEFI , RADEFI ,
     %                   NTFASU , MNFASU , NTSOFA , MNSOFA , IERR )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    GENERER LE MAILLAGE D'UNE SURFACE JOIGNANT DEUX LIGNES
C -----    DE TYPE "CONGE" (RACCORD DE DEUX SEGMENTS PAR TANGENCE)
C
C ENTREES:
C --------
C NTLXSU : NUMERO DE LA SURFACE DANS LE LEXIQUE DES SURFACES
C LADEFI : TABLEAU ENTIER DE DEFINITION DE LA SURFACE
C          CF ~/td/d/a_surface__definition
C RADEFI : TABLEAU REEL   DE DEFINITION DE LA SURFACE
C          CES 2 TABLEAUX LADEFI ET RADEFI ONT MEME ADRESSE A L'APPEL
C
C SORTIES:
C -------
C NTFASU : NUMERO      DU TMS 'NSEF' DES NUMEROS DES FACES
C MNFASU : ADRESSE MCN DU TMS 'NSEF' DES NUMEROS DES FACES
C          CF ~/td/d/a___nsef
C NTSOFA : NUMERO      DU TMS 'XYZSOMMET' DE LA SURFACE
C MNSOFA : ADRESSE MCN DU TMS 'XYZSOMMET' DE LA SURFACE
C          CF ~/td/d/a___xyzsommet
C IERR   : =0 SI PAS D'ERREUR
C          >0 SI ERREUR
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR :HENRI CHAJMOWICZ  ANALYSE NUMERIQUE UPMC PARIS    JANVIER 1991
C MODIFS :FERS ET PERRONNET ANALYSE NUMERIQUE UPMC PARIS  SEPTEMBRE 1996
C2345X7..............................................................012
      IMPLICIT INTEGER (W)
      include"./incl/gsmenu.inc"
      include"./incl/a_ligne__definition.inc"
      include"./incl/a_surface__definition.inc"
      include"./incl/a_ligne__bspline.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/ntmnlt.inc"
C
      COMMON / UNITES / LECTEU , IMPRIM , INTERA , NUNITE(29)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
C
      INTEGER LADEFI(0:*)
      INTEGER INT11,INT12,INT21,INT22
      INTEGER DISCRET1,DISCRET2
      INTEGER NUTG
C
      REAL RADEFI(0:*),XYZ6PT(3,6)
      REAL VD1(3),VD2(3),VD3(3),VD4(3),PVD1(3),PVD2(3)
      REAL CENTRE1(3),CENTRE2(3)
      REAL VAUX11(3),VAUX21(3),VAUX31(3),VAUX41(3)
      REAL VAUX12(3),VAUX22(3),VAUX32(3),VAUX42(3)
      REAL T11(3),T21(3)
      REAL T12(3),T22(3)
      REAL NB11(3),NB21(3),NB12(3),NB22(3)
C
C     vecteurs intermediaires utilises pour le calcul des tangentes
      REAL AC(3),T(3)
C
      REAL N1,N2,N3,N4
      REAL N11,N21,N31,N41,N12,N22,N32,N42
      REAL S,S1,S2,SAC,AL1,AL2
      REAL NX
      REAL P1,P2,P3,P4
      REAL LONG1,LONG2,LONG11,LONG12,LONG21,LONG22
      REAL AUX1,AUX2,ALPHA1,ALPHA2,AUX12,AUX22
C
      FRACP(R)=R-INT(R)
      PI = ATAN( 1.0 ) * 4
C
C     LES RAYONS
C     ----------
      RAY1CS = RADEFI( WAY1CS )
      IF( RAY1CS .LE. 0. ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'RAYON=<0'
         CALL LEREUR
         IERR = 3
         GOTO 9999
      ENDIF
C
      RAY2CS = RADEFI( WAY2CS )
      IF( RAY2CS .LE. 0. ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'RAYON=<0'
         CALL LEREUR
         IERR = 3
         GOTO 9999
      ENDIF
C
C     LE NOMBRE D'ARETES EN X ET Y DU CONGE
C     -------------------------------------
      NBIYCS = LADEFI(WBIYCS)
      IF( NBIYCS .LE. 2 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1)='NOMBRE INCORRECT <3  D ARETES'
         CALL LEREUR
         IERR = 2
         GOTO 9999
      ENDIF
C
      NBIXCS= LADEFI(WBIXCS)
      IF (NBIXCS .LE. 0) THEN
         NBLGRC(NRERR) = 1
         KERR(1)='NOMBRE INCORRECT <=0 DANS L''EPAISSEUR'
         CALL LEREUR
         IERR=2
         GOTO 9999
      ENDIF
C
C     LE NUMERO UTILISATEUR DES 6 POINTS DE DEFINITION DU CONGE
C     ---------------------------------------------------------
      DO 3 I=1,6
         IF( I .EQ. 1 ) THEN
            N = LADEFI( WUP1CS )
         ELSE IF( I .EQ. 2 ) THEN
            N = LADEFI( WUP2CS )
         ELSE IF( I .EQ. 3 ) THEN
            N = LADEFI( WUP3CS )
         ELSE IF( I .EQ. 4) THEN
            N = LADEFI( WUP4CS )
         ELSE IF( I .EQ. 5 ) THEN
            N = LADEFI( WUP5CS )
         ELSE IF( I .EQ. 6 ) THEN
            N = LADEFI( WUP6CS )
         ENDIF
C
C        RECUPERATION DES COORDONNEES DES POINTS
C        ---------------------------------------
         CALL LXNLOU( NTPOIN , N , NTLXPO , MN )
         IF( NTLXPO .LE. 0 ) THEN
             NBLGRC(NRERR) = 1
             WRITE(KERR(MXLGER)(1:4),'(I4)') I
             KERR(1) = 'POINT CONGE ' // KERR(MXLGER)(1:4) //' INCONNU'
             CALL LEREUR
             IERR = 6
             GOTO 9999
          ENDIF
          CALL LXTSOU( NTLXPO , 'XYZSOMMET' , NT , MN )
          MN = MN + WYZSOM
C
       DO 4 J=1,3
          XYZ6PT(J,I) = RMCN( MN +J-1)
4      CONTINUE
3     CONTINUE
C
      DO 5 J=1,3
          VD1(J)=XYZ6PT(J,1)-XYZ6PT(J,2)
          VD2(J)=XYZ6PT(J,3)-XYZ6PT(J,2)
          VD3(J)=XYZ6PT(J,4)-XYZ6PT(J,5)
          VD4(J)=XYZ6PT(J,6)-XYZ6PT(J,5)
5     CONTINUE
C
C     ANGLE, NORME,PRODUIT SCALAIRE DES VECTEURS DONNES
C     -------------------------------------------------
      N1 = SQRT( PROSCR(VD1,VD1,3) )
      N2 = SQRT( PROSCR(VD2,VD2,3) )
      S1 = PROSCR(VD1,VD2,3)
C
      N3 = SQRT( PROSCR(VD3,VD3,3) )
      N4 = SQRT( PROSCR(VD4,VD4,3) )
      S2 = PROSCR(VD3,VD4,3)
C
      AUX1 = S1/(N1*N2)
      AUX1 = SQRT(0.5*(AUX1+1))
      ALPHA1 = 2*ACOS(AUX1)
C
      AUX2 = S2/(N3*N4)
      AUX2 = SQRT(0.5*(AUX2+1))
      ALPHA2 = 2*ACOS(AUX2)
C
      IF (ALPHA1 .EQ. 0.) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'ANGLE DU CONGE NUL'
         CALL LEREUR
         IERR = 3
         GOTO 9999
      ENDIF
C
      IF (ALPHA2 .EQ. 0.) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'ANGLE DU CONGE NUL'
         CALL LEREUR
         IERR = 3
         GOTO 9999
      ENDIF
C
C     VERIFICATION DE L'ORDRE D'ENTREE DES POINTS
C     -------------------------------------------
      CALL PROVER(VD1,VD2,PVD1)
      CALL PROVER(VD3,VD4,PVD2)
C
      S = PROSCR(PVD1,PVD2,3)
C
      IF (S .LE. 0.) THEN
         NBLGRC(NRERR) = 1
         KERR(1)='P4 ET P6 SONT INVERSES'
         CALL LEREUR
         IERR=7
         GOTO 9999
      ENDIF
C
      DO 70 I=1,3
         PVD1(I) = PVD1(I)/(N1*N2*SIN(ALPHA1))
         PVD2(I) = PVD2(I)/(N3*N4*SIN(ALPHA2))
70    CONTINUE
C
C     CALCUL DES POINTS DE TANGENCE DU CONGE AVEC LES SEGMENTS
C     --------------------------------------------------------
      AUX12 = RAY1CS/TAN(0.5*ALPHA1)
      AUX22 = RAY2CS/TAN(0.5*ALPHA2)
C
      DO 22 I=1,3
          T11(I) = XYZ6PT(I,2) + (AUX12*VD1(I))/N1
          T21(I) = XYZ6PT(I,2) + (AUX12*VD2(I))/N2
          T12(I) = XYZ6PT(I,5) + (AUX22*VD3(I))/N3
          T22(I) = XYZ6PT(I,5) + (AUX22*VD4(I))/N4
22    CONTINUE
C
      CALL PROVER(PVD1,VD1,VAUX11)
      CALL PROVER(PVD2,VD3,VAUX12)
      DO 24 I=1,3
         VAUX21(I)=XYZ6PT(I,1)-T11(I)
         VAUX31(I)=XYZ6PT(I,3)-T21(I)
C
         VAUX22(I)=XYZ6PT(I,4)-T12(I)
         VAUX32(I)=XYZ6PT(I,6)-T22(I)
24    CONTINUE
C
      P1 = PROSCR(VAUX21,VD1,3)
      P2 = PROSCR(VAUX31,VD2,3)
      P3 = PROSCR(VAUX22,VD3,3)
      P4 = PROSCR(VAUX32,VD4,3)
C
      IF (P1 .LE. 0) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'RAY1CS TROP GRAND'
         CALL LEREUR
         IERR = 4
         GOTO 9999
      ENDIF
C
      IF (P2 .LE. 0) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'RAY1CS TROP GRAND'
         CALL LEREUR
         IERR = 4
         GOTO 9999
      ENDIF
C
      IF (P3 .LE. 0) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'RAY2CS TROP GRAND'
         CALL LEREUR
         IERR = 4
         GOTO 9999
      ENDIF
C
      IF (P4 .LE. 0) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'RAY2CS TROP GRAND'
         CALL LEREUR
         IERR = 4
         GOTO 9999
      ENDIF
C
C     CALCUL DES CENTRES DES CERCLES
C     ------------------------------
      N11 = SQRT( PROSCR(VAUX11,VAUX11,3) )
      N12 = SQRT( PROSCR(VAUX12,VAUX12,3) )
      DO 25 I=1,3
         CENTRE1(I)=T11(I)+RAY1CS*VAUX11(I)/N1
         CENTRE2(I)=T12(I)+RAY2CS*VAUX12(I)/N3
         VAUX41(I) = T11(I) - CENTRE1(I)
         VAUX42(I) = T12(I) - CENTRE2(I)
25    CONTINUE
C
      ALPHA1 = PI - ALPHA1
      ALPHA2 = PI - ALPHA2
C
      N21 = SQRT( PROSCR(VAUX21,VAUX21,3) )
      N31 = SQRT( PROSCR(VAUX31,VAUX31,3) )
C
      N22 = SQRT( PROSCR(VAUX22,VAUX22,3) )
      N32 = SQRT( PROSCR(VAUX32,VAUX32,3) )
C
      LONG1 = N21 + N31 + (ALPHA1*RAY1CS)
      LONG1 = LONG1/REAL(NBIYCS)
      LONG11 = LONG1
      LONG21 = LONG1
C
      IF (FRACP(N21/LONG1) .EQ. 0.) THEN
          INT11 = INT(N21/LONG1)
      ELSE IF (INT(N21/LONG1) .NE. 0) THEN
          LONG11 = N21/REAL(INT(N21/LONG1))
          INT11 = INT(N21/LONG1)
      ELSE
         NBLGRC(NRERR) = 1
         KERR(1) = 'PAS ASSEZ D ARETES'
         CALL LEREUR
         IERR = 10
         GOTO 9999
      ENDIF
C
      IF (FRACP(N31/LONG1) .EQ. 0.) THEN
          INT21 = INT(N31/LONG1)
      ELSE IF (INT(N31/LONG1) .NE. 0) THEN
          LONG21 = N31/REAL(INT(N31/LONG1))
          INT21 = INT(N31/LONG1)
      ELSE
         NBLGRC(NRERR) = 1
         KERR(1) = 'PAS ASSEZ D ARETES'
         CALL LEREUR
         IERR = 10
         GOTO 9999
      ENDIF
C
      DISCRET1 = NBIYCS - (INT11+INT21)
C
      ALPHA1 = ALPHA1/REAL(DISCRET1)
C
      N41 = SQRT( PROSCR(VAUX41,VAUX41,3) )
C
      DO 26 I=1,3
         NB11(I) = VAUX41(I)/N41
26    CONTINUE
C
      CALL PROVER(PVD1,NB11,NB21)
C
      LONG2 = N22 + N32 + (ALPHA2*RAY2CS)
      LONG2 = LONG2/REAL(NBIYCS)
      LONG12 = LONG2
      LONG22 = LONG2
C
      IF (FRACP(N22/LONG2) .EQ. 0.) THEN
          INT12 = INT(N22/LONG2)
      ELSE
          LONG12 = N22/REAL(INT(N22/LONG2))
          INT12 = INT(N22/LONG2)
      ENDIF
C
      IF (FRACP(N32/LONG2) .EQ. 0.) THEN
          INT22 = INT(N32/LONG2)
      ELSE
          LONG22 = N32/REAL(INT(N32/LONG2))
          INT22 = INT(N32/LONG2)
      ENDIF
C
      IF(INT11.NE.INT12) THEN
         INT12=INT11
         LONG12=N22/REAL(INT11)
      ENDIF
C
      IF (INT21.NE.INT22) THEN
         INT22=INT21
         LONG22=N32/REAL(INT21)
      ENDIF
C
      DISCRET2 = NBIYCS - (INT12+INT22)
C
      ALPHA2 = ALPHA2/REAL(DISCRET2)
C
      N42 = SQRT( PROSCR(VAUX42,VAUX42,3) )
C
      DO 27 I=1,3
         NB12(I) = VAUX42(I)/N42
27    CONTINUE
C
      CALL PROVER(PVD2,NB12,NB22)
C
C     CONSTRUCTION DU TABLEAU 'XYZSOMMET'
C     ----------------------------------
      NBSOM = (NBIXCS+1)*(NBIYCS+1)
      NBTGS = (NBIXCS+1)*(DISCRET1+1)
      CALL LXTNDC( NTLXSU, 'XYZSOMMET', 'MOTS', WYZSOM+3*NBSOM+3*NBTGS)
      CALL LXTSOU( NTLXSU , 'XYZSOMMET' ,  NTSOFA , MNSOFA )
C
C     LE NOMBRE DE SOMMETS
      MCN( MNSOFA + WNBSOM ) = NBSOM
C
C     LE NOMBRE DE TANGENTES
      MCN( MNSOFA + WNBTGS ) = NBTGS
C
C     POINTEUR SUR LE DEBUT DE STOCKAGE DES SOMMETS
      MN = MNSOFA + WYZSOM
C     POINTEUR SUR LE DEBUT DE STOCKAGE DES TANGENTES
      MNT = MNSOFA + WYZSOM + 3*NBSOM
      NX = REAL (NBIXCS)
C
      DO 31 J=0,NBIYCS
C
        IF (J .LT. INT11) THEN
C
C         partie plane superieure
          DO 41 I=0,NBIXCS
            DO 11 K=1,3
              RMCN(MN+K-1)=(XYZ6PT(K,1)-(VD1(K)*REAL(J)*LONG11/N1))
     %                    *(1.-REAL(I)/NX) + (REAL(I)/NX)*
     %                    (XYZ6PT(K,4)-(VD3(K)*REAL(J)*LONG12/N3))
11            CONTINUE
              MN = MN+3
41        CONTINUE
        ELSE IF (J.GT.(NBIYCS -INT21)) THEN
C
C         partie plane inferieure
          JL=NBIYCS-J
          DO 44 I=0,NBIXCS
            DO 14 K=1,3
              RMCN(MN+K-1)=(XYZ6PT(K,3)-(VD2(K)*REAL(JL)*LONG21/N2))
     %                   *(1.-REAL(I)/NX)+(REAL(I)/NX)*
     %                   (XYZ6PT(K,6)-(VD4(K)*REAL(JL)*LONG22/N4))
14           CONTINUE
             MN = MN+3
44         CONTINUE
         ELSEIF (J.GE.INT11.AND.J.LE.(NBIYCS-INT21)) THEN
C
C          partie courbe
           AL1=REAL(J-INT11)*ALPHA1
           AL2=REAL(J-INT11)*ALPHA2
           DO 350 I=0,NBIXCS
C             calcul des 3 coordonnes au sommet i,j
              OMEGA = REAL(I) / NX
              OMEGA1 = 1.0 - OMEGA
              DO 19 K=1,3
                 RMCN(MN+K-1)=(CENTRE1(K)+N41*COS(AL1)*NB11(K)
     %                 - N41*SIN(AL1)*NB21(K))*OMEGA1
     %                 +(CENTRE2(K)+N42*COS(AL2)*NB12(K)
     %                 -N42*SIN(AL2)*NB22(K))*OMEGA
19            CONTINUE
C
C             calcul de la tangente au sommet i,j
              DO 300 K=1,3
C                Les 3 coordonnees du vecteur : sommet i,j -> centre du cercle i
                 AC(K)=CENTRE1(K) + OMEGA * (CENTRE2(K)-CENTRE1(K))
     %                -RMCN(MN+K-1)
C                Les 3 coordonnees du vecteur normal au plan des 3 points
                 VAUX21(K)=PVD1(K) + OMEGA * (PVD2(K)-PVD1(K))
300           CONTINUE
C             le vecteur unitaire Sommet i,j -> centre du cercle
              SAC=SQRT( PROSCR(AC,AC,3) )
              DO 301 K=1,3
                 AC(K)=AC(K)/SAC
301           CONTINUE
C             le vecteur tangent au cercle au sommet i,j
              CALL PROVER(VAUX21,AC,T)
              SAC = SQRT( PROSCR(T,T,3) )
              AL  = ALPHA1 + OMEGA * (ALPHA2-ALPHA1)
              RAY = RAY1CS + OMEGA * (RAY2CS-RAY1CS)
              SAC = RAY * AL / SAC
              DO 302 K=1,3
                 RMCN(MNT+K-1) = T(K) * SAC
302           CONTINUE
              MN =MN +3
              MNT=MNT+3
350        CONTINUE
C
        ENDIF
31    CONTINUE
C
C     AJOUT DE LA DATE
      CALL ECDATE( MCN(MNSOFA) )
C     LE NOMBRE DE COORDONNEES PAR SOMMET
      MCN( MNSOFA + WBCOOR ) = 3
C     AJOUT DU NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNSOFA + MOTVAR(6) ) = NONMTD( '~>>>XYZSOMMET' )
C
C     CONSTRUCTION DU TABLEAU 'NSEF' SURFACE STRUCTUREE
C     -------------------------------------------------------
      NBEF   = NBIXCS*NBIYCS
      NBEFTG = (NBIYCS - INT11 - INT21) * NBIXCS
      CALL LXTNDC( NTLXSU, 'NSEF', 'ENTIER', WBARYQ+1+NBEF+9*NBEFTG )
      CALL LXTSOU( NTLXSU, 'NSEF',  NTFASU , MNFASU )
C     LE TYPE DE L'OBJET : ICI SURFACE
      MCN( MNFASU + WUTYOB ) = 3
C     LE TYPE DU MAILLAGE : ICI QUADRANGLE STRUCTURE
      MCN( MNFASU + WUTYMA ) = 4
C     LE NOMBRE D'ARETES DU QUADRANGLE STRUCTURE
      MCN( MNFASU + WBARYQ ) = NBIYCS
      MCN( MNFASU + WBARXQ ) = NBIXCS
C     LE NOMBRE D'EF DU MAILLAGE
      MCN( MNFASU + WBEFOB ) = NBIXCS * NBIYCS
C     LE NOMBRE DE SOMMETS PAR ELEMENT FINI
      MCN( MNFASU + WBSOEF ) = 4
C     LE NOMBRE DE TANGENTES PAR ELEMENT FINI
      MCN( MNFASU + WBTGEF ) = 8
C     LE POINTEUR SUR LES EF A TG
      MCN( MNFASU + WBEFAP ) = NBEF
C     LE NOMBRE D'EF A TG
      MCN( MNFASU + WBEFTG ) = NBEFTG
C
C     POINTEUR SUR LES EF A TG
      NBT = 0
      MN  = MNFASU + WBARYQ
      DO 60 J=1,NBIYCS
         DO 50 I=1,NBIXCS
            MN = MN + 1
            IF( J.GT.INT11 .AND. J.LE.(NBIYCS-INT21) ) THEN
C              PARTIE COURBE A TG
               NBT = NBT + 1
               MCN(MN) = NBT
            ELSE
C              PARTIE PLANE
               MCN(MN) = 0
            ENDIF
 50      CONTINUE
 60   CONTINUE
C
C     LE NUMERO GEOMETRIQUE ICI TRONC DE CONE => NO 14
      DO 80 J=1,NBIYCS
         DO 72 I=1,NBIXCS
            IF( J.GT.INT11 .AND. J.LE.(NBIYCS-INT21) ) THEN
C              PARTIE COURBE A TG
               MN = MN + 1
               MCN(MN) = 14
            ENDIF
 72      CONTINUE
 80   CONTINUE
C
C     NUMERO DES TG DES EF A TG
      NUTG = 0
      DO 100 J=INT11+1,NBIYCS-INT21
           DO 105 I=1,NBIXCS
C             RANGEMENT DES TANGENTES PAR SOMMETS DE L'EF
              NUTG = NUTG + 1
              MCN(MN+1) = 0
              MCN(MN+2) = NUTG
              MCN(MN+3) = NUTG+1
              MCN(MN+4) = 0
              MCN(MN+5) = 0
              MCN(MN+6) = -(NUTG+2+NBIXCS)
              MCN(MN+7) = -(NUTG+1+NBIXCS)
              MCN(MN+8) = 0
              MN=MN+8
105        CONTINUE
           NUTG = NUTG + 1
100   CONTINUE
C
C     AJOUT DE LA DATE
      CALL ECDATE( MCN(MNFASU) )
C     AJOUT DU NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNFASU + MOTVAR(6) ) = NONMTD( '~>>>NSEF' )
C
      IF (IERR .NE. 0) THEN
          CALL LXTSDS ( NTLSU, 'XYZSOMMET')
          CALL LXTSDS ( NTLSU, 'NSEF')
          GOTO 9999
      ENDIF
C
C     ERREUR
C     ======
9999  RETURN
      END
