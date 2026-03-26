      SUBROUTINE SUEX05( NTLXSU , LADEFI ,
     %                   NTFASU , MNFASU , NTSOSU , MNSOSU , IERR )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    GENERER LES FACES DU TRIANGLE BEZIER
C -----    D'APRES DES POINTS DE CONTROLE
C
C ENTREES:
C --------
C NTLXSU : NUMERO DU TABLEAU TS DU LEXIQUE DE LA SURFACE
C LADEFI : TABLEAU ENTIER DE DEFINITION DE LA SURFACE
C
C SORTIES:
C --------
C NTFASU : NUMERO      DU TMS 'NSEF' DES NUMEROS DES FACES
C MNFASU : ADRESSE MCN DU TMS 'NSEF' DES NUMEROS DES FACES
C NTSOSU : NUMERO      DU TMS 'XYZSOMMET' DE LA SURFACE
C MNSOSU : ADRESSE MCN DU TMS 'XYZSOMMET' DE LA SURFACE
C IERR   : =0 SI PAS D'ERREUR
C          >0 SI ERREUR
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEURS: CLEVEDE ZELMANSE ANALYSE NUMERIQUE UPMC PARIS    FEVRIER 1991
C2345X7..............................................................012
      IMPLICIT INTEGER (W)
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU , IMPRIM , INTERA , NUNITE(29)
      include"./incl/a_ligne__definition.inc"
      include"./incl/a_ligne__bspline.inc"
      include"./incl/a_surface__definition.inc"
      include"./incl/a_surface__bspline.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
      INTEGER           LADEFI(0:*)
C
C     VERIFICATIONS DES DONNEES
C     =========================
      IERR   = 0
C
C     LE NOMBRE D'ARETES
      NBARTB = LADEFI( WBARTB )
      IF( NBARTB .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'NOMBRE INCORRECT D''ARETE'
         CALL LEREUR
         IERR = 2
         GOTO 9999
      ENDIF
C
C     LE DEGRE DES POLYNOMES
      LDEGTB = LADEFI( WDEGTB )
      IF( LDEGTB .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'DEGRE INCORRECT DES POLYNOMES'
         CALL LEREUR
         IERR = 2
         GOTO 9999
      ENDIF
C
C     LE NOMBRE DE POINTS DE CONTROLE :
      NBPTTB = ( (LDEGTB+1)*(LDEGTB+2) )/2
C
C     LECTURE DES 3 COORDONNEES DES POINTS DE CONTROLE
C     TABLEAU POINTB(1:3,1:NBPTTB)
      CALL TNMCDC( 'REEL' , 3*NBPTTB , MNPTSB )
      MNSB = MNPTSB
      DO 10 I=1,NBPTTB
C        LE NUMERO DU POINT
         N = LADEFI( WUPTTB - 1 + I )
C        OUVERTURE DU LEXIQUE DU POINT
         CALL LXNLOU( NTPOIN , N , NTLXPO , MN )
         IF( NTLXPO .LE. 0 ) THEN
            NBLGRC(NRERR) = 1
            WRITE(KERR(MXLGER)(1:4),'(I4)') I
            KERR(1) = 'POINT INCONNU ' // KERR(MXLGER)(1:4)
            CALL LEREUR
            IERR = 7
            GOTO 10
         ENDIF
C        OUVERTURE DU TABLEAU SOMMETS
         CALL LXTSOU( NTLXPO , 'XYZSOMMET' , NT , MNS )
         IF( NT .LE. 0 ) THEN
            NBLGRC(NRERR) = 1
            WRITE(KERR(MXLGER)(1:4),'(I4)') I
            KERR(1) = 'POINT ' // KERR(MXLGER)(1:4) //
     %                ' SANS COORDONNEES'
            CALL LEREUR
            IERR = 7
            GOTO 10
         ENDIF
C        LES 3 COORDONNEES DU POINT
         MNS = MNS + WYZSOM
         RMCN( MNSB     ) = RMCN( MNS )
         RMCN( MNSB + 1 ) = RMCN( MNS+1 )
         RMCN( MNSB + 2 ) = RMCN( MNS+2 )
         MNSB = MNSB + 3
 10   CONTINUE
      IF( IERR .NE. 0 ) GOTO 9999
C
C     LE TABLEAU DES N!
      CALL TNMCDC( 'REEL' , LDEGTB+1 , MNFACM )
C
C     CONSTRUCTION DU TABLEAU 'XYZSOMMET'
C     ---------------------------------
      NBSOSU = ( (NBARTB+2)*(NBARTB+1) )/2
      CALL LXTNDC( NTLXSU , 'XYZSOMMET' , 'ENTIER' , WYZSOM + 3*NBSOSU )
      CALL LXTSOU( NTLXSU , 'XYZSOMMET' ,  NTSOSU  , MNSOSU )
C     LE NOMBRE DE SOMMETS
      MCN( MNSOSU + WNBSOM ) = NBSOSU
C
C     CALCUL DES POINTS DE LA TRIANGULATION
C     =====================================
      CALL BEZIER( NBARTB, LDEGTB, NBSOSU, RMCN(MNSOSU+WYZSOM), NBPTTB,
     %             RMCN(MNPTSB), RMCN(MNFACM) )
C     DESTRUCTION DES TABLEAUX DEVENUS INUTILES
      CALL TNMCDS( 'REEL' , 3*NBPTTB, MNPTSB )
      CALL TNMCDS( 'REEL' , LDEGTB+1, MNFACM )
C     AJOUT DE LA DATE
      CALL ECDATE( MCN(MNSOSU) )
C
C     LE NOMBRE DE COORDONNEES D'UN SOMMET
      MCN( MNSOSU + WBCOOR ) = 3
C     AJOUT DU NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNSOSU + WNBTGS ) = 0
      MCN( MNSOSU + MOTVAR(6) ) = NONMTD( '~>>>XYZSOMMET' )
C
C     CONSTRUCTION DU TABLEAU 'NSEF' SURFACE STRUCTUREE
C     -------------------------------------------------------
      CALL LXTNDC( NTLXSU , 'NSEF' , 'ENTIER' , 1 + WBARTR )
      CALL LXTSOU( NTLXSU , 'NSEF' ,  NTFASU  , MNFASU )
C     LE TYPE DE L'OBJET : ICI SURFACE
      MCN( MNFASU + WUTYOB ) = 3
C     LE TYPE DU MAILLAGE : ICI TRIANGLE STRUCTURE
      MCN( MNFASU + WUTYMA ) = 3
C     LE NOMBRE DE SOMMETS PAR FACE
      MCN ( MNFASU + WBSOEF ) = 4
C     LE NOMBRE DE TANGENTES STOCKEES PAR EF : SURFACE C0
      MCN ( MNFASU + WBTGEF ) = 0
      MCN ( MNFASU + WBEFAP ) = 0
      MCN ( MNFASU + WBEFTG ) = 0
C     LE NOMBRE D'EF DE LA SURFACE
      MCN ( MNFASU + WBEFOB ) = NBARTB * NBARTB
C     LE NOMBRE D'ARETES DU TRIANGLE STRUCTURE
      MCN( MNFASU + WBARXQ ) = NBARTB
C     LE TYPE NON-FERME DE FERMETURE DU MAILLAGE
      MCN( MNFASU + WUTFMA ) = 0
C     AJOUT DE LA DATE
      CALL ECDATE( MCN(MNFASU) )
C     AJOUT DU NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNFASU + MOTVAR(6) ) = NONMTD( '~>>>NSEF' )
 9999 RETURN
      END
      SUBROUTINE BEZIER ( NBARTB, LDEGTB, NBSOSU, SOMTB, NBPTTB,
     %                    PTCTTB, FACT )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    GENERER UNE TRIANGULATION D'UN TRIANGLE DE BEZIER.
C------
C ENTREE :
C --------
C NBARTB : NOMBRE D'ARETES D'UN COTE DE LA TRIANGULATION.
C LDEGTB : DEGRE DES POLYNOMES
C NBSOSU : NOMBRE DE SOMMETS
C NBPTTB : NOMBRE DE POINTS DE CONTROLE
C PTCTTB : TABLEAU DES POINTS DE CONTROLE
C FACT   : TABLEAU DES N!
C
C SORTIE :
C --------
C SOMTB  : TABLEAU DES SOMMETS DE LA TRIANGULATION
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEURS: CLEVEDE ZELMANSE ANALYSE NUMERIQUE UPMC PARIS    FEVRIER 1991
C2345X7..............................................................012
      INTEGER  NBARTB, LDEGTB, NBSOSU, NBPTTB
      REAL     SOMTB(3,NBSOSU), PTCTTB(3,NBPTTB), FACT(0:LDEGTB)
      INTEGER  NUMPCT, I, J, K, COORD, PT, L
      REAL     COEFF, BETA(1:3), R, S, T, GAMMA
C
      FACT(0) = 1.
      DO 90 N=1,LDEGTB
         FACT(N) = N*FACT(N-1)
90    CONTINUE
C
      DO 50 PT=1,NBSOSU
         DO 60 COORD=1,3
            SOMTB(COORD,PT)=0.
60       CONTINUE
50    CONTINUE
C
      DO 10 NUMPCT=1,NBPTTB
C        BOUCLE SUR LES POINTS DE CONTROLE.
C
C        CALCUL DU TRIPLET (I,J,K) CORRESPONDANT AU NUMPCT IEME PT DE CTRL:
         CALL ATOIJK (NUMPCT,LDEGTB,I,J,K)
C
         COEFF=FACT(LDEGTB)/( FACT(I)*FACT(J)*FACT(K) )
C
         DO 20 COORD=1,3
            BETA(COORD)=COEFF*PTCTTB(COORD,NUMPCT)
20       CONTINUE
C
         DO 30 PT=1,NBSOSU
C
C           CALCUL DES COORDONNEES PARAMETRIQUES DU PT IEME SOMMET:
            CALL COPARA (PT,NBARTB,R,S,T)
C
            GAMMA = 1.
            DO 110 L=1,I
               GAMMA = GAMMA*R
110         CONTINUE
            DO 120 L=1,J
               GAMMA = GAMMA*S
120         CONTINUE
            DO 130 L=1,K
               GAMMA = GAMMA*T
130         CONTINUE
C
            DO 40 COORD=1,3
               SOMTB(COORD,PT)=SOMTB(COORD,PT)+BETA(COORD)*GAMMA
40          CONTINUE
30      CONTINUE
10    CONTINUE
      RETURN
      END
      SUBROUTINE  COPARA (NUMPT,NBARTB,R,S,T)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : DETERMINATION DES COORDONNEES BARYCENTRIQUES DU POINT NUMPT
C ----- DANS UN TRIANGLE DE REFERENCE (TRIANGULATION DE DELAUNAY).
C
C           6  ...
C           3   5  ...
C           1   2   4  ...
C
C ENTREE :
C --------
C NUMPT : NUMERO DU POINT (NUMEROTATION PAR DIAGONALES, DE
C         GAUCHE A DROITE ET DE BAS EN HAUT [CF SUEXT3]).
C NBARTB : NOMBRE D'ARETE(S) PAR COTE.
C
C SORTIE :
C --------
C R,S,T  : COORDONNEES BARYCENTRIQUES.
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEURS: CLEVEDE ZELMANSE ANALYSE NUMERIQUE UPMC PARIS    FEVRIER 1991
C2345X7..............................................................012
      INTEGER  NUMPT,NBARTB
      REAL     R,S,T
      INTEGER  NDIAG
C
      NDIAG=0
 10   NDIAG=NDIAG+1
      IF (NUMPT.GT.((NDIAG*(NDIAG+1))/2)) GOTO 10
C
C     NDIAG EST LE NUMERO DE LA DIAGONALE QUI PORTE LE POINT NUMPT
      S=( NUMPT - ((NDIAG-1)*NDIAG)/2 - 1.0 ) / NBARTB
      R=( NDIAG-1.0 ) / NBARTB - S
      T= 1.0 - R - S
C
      RETURN
      END
      SUBROUTINE  ATOIJK (NUMPCT,LDEGTB,I,J,K)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : DETERMINATION DU TRIPLET (I,J,K) ASSOCIE AU POINT DE CONTROLE
C ----- NUMPCT CLASSE PAR I DECROISSANTS PUIS J DECROISSANTS.
C
C      PAR EXEMPLE, POUR LDEGTB=3, ON A :
C
C   1       2       3       4       5       6       7       8       9      10
C(3,0,0) (2,1,0) (2,0,1) (1,2,0) (1,1,1) (1,0,2) (0,3,0) (0,2,1) (0,1,2) (0,0,3)
C
C ENTREE :
C --------
C NUMPCT : NUMERO DU POINT DE CONTROLE.
C LDEGTB : DEGRE DES POLYNOMES (I+J+K=LDEGTB).
C
C SORTIE :
C --------
C I,J,K : TRIPLET INDICE DU POINT NUMPCT.
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEURS: CLEVEDE ZELMANSE ANALYSE NUMERIQUE UPMC PARIS    FEVRIER 1991
C2345X7..............................................................012
      INTEGER  NUMPCT,LDEGTB,I,J,K
      INTEGER  U
C
      U = 0
  10  U = U+1
      IF ( ((U*(U+1))/2).LT.(NUMPCT) ) GOTO 10
C
      I = LDEGTB-U+1
      K = NUMPCT-( (U*(U-1))/2 + 1 )
      J = LDEGTB-I-K
C
      RETURN
      END
