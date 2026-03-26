      SUBROUTINE VOE114( NUTYEX, LETORE, NBTREX, MNTREX,
     %                   NOFONC,
     %                   VECTOR, Z,
     %                   NBSOMS, COOR2, COOR3,
     %                   NBTGSU, XYZTGS, XYZTGV,
     %                   IERR   )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : GENERATION DES COORDONNEES DES DIFFERENTES COUCHES D UN MAILLAGE
C ----- ET DES TANGENTES DANS LE CAS NUTYEX=0 OU 1
C
C ENTREES :
C ---------
C NUTYEX : 0 => 3 POINTS UTILISATEUR DANS CHAQUE SURFACE
C          1 => 3 LIGNES UTILISATEUR => TYPE 0
C          2 => VECTOR(3,NBTREX)=Z MAX DECOUPEE REGULIEREMENT
C          3 => COTE MAX DECOUPEE IRREGULIEREMENT SELON VECTOR(3,I)
C               1 <= I <= NBTREX
C          4 => VECTOR(J,1) 1 =< J =< 3 VECTEUR CONSTANT DE TRANSLATION
C          5 => VECTOR(J,I) 1 =< J =< 3  1 =< I =< NBTREX
C               VECTEUR DE TRANSLATION DE LA I+1 - EME COUCHE
C               LA PREMIERE COUCHE EST DANS LE PLAN XOY ( Z = 0 )
C          6 => 3 NUMEROS DES FONCTIONS DE TRANSFORMATION DE LA SURFACE
C               INITIALE EN CHAQUE COUCHE
C LETORE : 1 SI DERNIERE SURFACE = PREMIERE SURFACE
C          0 SINON
C NBTREX : NOMBRE DE COUCHES OU TRANCHES DU VOLUME A MAILLER
C MNTREX : L'ADRESSE RMCN DE LA 1-ERE COORDONNEE DE CHACUN DES 3 SOMMETS
C          DE DEFINITION DU PLAN DE CHAQUE SURFACE
C VECTOR : DECLARE (3,NBTREX) . IL EST REMPLI SELON NUTYEX
C Z      : DECLARE ( NBTREX + 1 ) TABLEAU AUXILIAIRE
C NBSOMS : NOMBRE DE SOMMETS DE LA SURFACE INITIALE
C COOR2  : COORDONNEES DES SOMMETS DE LA SURFACE INITIALE
C NBTGSU : NOMBRE DE TANGENTES DE LA SURFACE
C XYZTGS : LES 3 COMPOSANTES DES NBTGSU TANGENTES DE LA SURFACE
C
C SORTIES :
C ---------
C COOR3  : COORDONNEES DES NBSOMS * ( NBTREX + 1 ) SOMMETS DU VOLUME
C XYZTGV : LES 3 COMPOSANTES DES TANGENTES DU VOLUME
C IERR   : CODE D'ERREUR EN EXECUTION 0 PAS D'ERREUR, NON NUL SINON
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS UPMC    NOVEMBRE 1996
C23456...............................................................012
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      include"./incl/pp.inc"
      COMMON            RMCN(MOTMCN)
      INTEGER           MNTREX(1:3,0:NBTREX)
      INTEGER           NOFONC(3)
      DOUBLE PRECISION  DP(3,3),DPM1(3,3),P1P2(3),P1P3(3),P1P4(3)
      DOUBLE PRECISION  A,DPARAF(4)
      REAL              VECTOR(3,NBTREX),Z(NBTREX),
     %                  COOR2(3,NBSOMS),COOR3(3,*),
     %                  XYZTGS(3,NBTGSU), XYZTGV(3,*),
     %                  V(3)
C
C     LE TRAITEMENT SELON L'OPTION NUTYEX
C     ===================================
      NBTRE1 = NBTREX + 1
 1    GOTO( 2 , 2 , 110 , 200 , 300 , 300 , 500 ) , NUTYEX+1
C
C     NUTYEX 0 : 3 POINTS SUR CHAQUE SURFACE DEFINISSENT LA TRANSFORMATION
C     --------------------------------------------------------------------
C     LES 3 COORDONNEES DES NBSOMS SOMMETS DE LA SURFACE INITIALE
 2    DO 5 J=1,NBSOMS
         DO 3 I=1,3
            COOR3(I,J) = COOR2(I,J)
 3       CONTINUE
 5    CONTINUE
      DO 8 J=1,NBTGSU
         DO 6 I=1,3
            XYZTGV(I,J) = XYZTGS(I,J)
 6       CONTINUE
 8    CONTINUE
C
C     LE TRIEDRE DE LA SURFACE INITIALE
C     CALCUL DE LA MATRICE  PK-P1 POUR K=2,3,4
      DO 10 I=0,2
C        LA DIFFERENCE PK-P1 POUR K=2,3
         P1P2(I+1) = RMCN(MNTREX(2,0)+I) - RMCN(MNTREX(1,0)+I)
         P1P3(I+1) = RMCN(MNTREX(3,0)+I) - RMCN(MNTREX(1,0)+I)
 10   CONTINUE
C     CALCUL DE P1P4 = PRODUIT VECTORIEL P1P2 P1P3
      CALL PROVEC( P1P2, P1P3, P1P4 )
C     HOMOGENEISATION A UNE LONGUEUR DE VECTEUR
      A = 2D0 / ( SQRT( P1P2(1)**2+P1P2(2)**2+P1P2(3)**2 )
     %          + SQRT( P1P3(1)**2+P1P3(2)**2+P1P3(3)**2 ) )
      P1P4(1) = P1P4(1) * A
      P1P4(2) = P1P4(2) * A
      P1P4(3) = P1P4(3) * A
C
C     LA BOUCLE SUR LES SURFACES 1 A NBTREX-LETORE
      NBSOM0 = 0
      NBSOM1 = NBSOMS
      NBTGS0 = 0
      NBTGS1 = NBTGSU
      DO 100 NBTR = 1 , NBTREX - LETORE
C
C        CALCUL DE LA MATRICE DU DEPLACEMENT AMENANT P10 EN P1
C        P20 EN P2 , P30 EN P3 ET P40 EN P4 AVEC P4 OBTENU
C        PAR PRODUIT VECTORIEL P1P2 AVEC P1P3
C        =====================================================
C        REMPLISSAGE DE LA MATRICE A INVERSER
         DO 20 J=1,3
            DP(1,J) = P1P2(J)
            DP(2,J) = P1P3(J)
            DP(3,J) = P1P4(J)
 20      CONTINUE
C
C        INVERSION DE LA MATRICE DP
         CALL DFM1R3( DP, DPM1, IERR )
         IF( IERR .NE. 0 ) THEN
            NBLGRC(NRERR) = 2
            KERR(1) = 'PASSAGE REPERE A REPERE SUIVANT IMPOSSIBLE'
            KERR(2) = 'REVOYEZ LES REPERES FORMES PAR LES 3 POINTS'
            CALL LEREUR
            RETURN
         ENDIF
C
C        CALCUL DE LA MATRICE  PK-P1 POUR K=2,3,4
         DO 30 I=0,2
C           LA DIFFERENCE PK-P1 POUR K=2,3
            P1P2(I+1) = RMCN(MNTREX(2,NBTR)+I) - RMCN(MNTREX(1,NBTR)+I)
            P1P3(I+1) = RMCN(MNTREX(3,NBTR)+I) - RMCN(MNTREX(1,NBTR)+I)
 30      CONTINUE
C        CALCUL DE P1P4 = PRODUIT VECTORIEL P1P2 P1P3
         CALL PROVEC( P1P2, P1P3, P1P4 )
C        HOMOGENEISATION A UNE LONGUEUR DE VECTEUR
         A = 2D0 / ( SQRT( P1P2(1)**2+P1P2(2)**2+P1P2(3)**2 )
     %             + SQRT( P1P3(1)**2+P1P3(2)**2+P1P3(3)**2 ) )
         P1P4(1) = P1P4(1) * A
         P1P4(2) = P1P4(2) * A
         P1P4(3) = P1P4(3) * A
C
C        LE CALCUL DES 9 COEFFICIENTS DE LA MATRICE DE DEPLACEMENT
         DO 40 I=1,3
          DP(I,1)=DPM1(1,1)*P1P2(I)+DPM1(1,2)*P1P3(I)+DPM1(1,3)*P1P4(I)
          DP(I,2)=DPM1(2,1)*P1P2(I)+DPM1(2,2)*P1P3(I)+DPM1(2,3)*P1P4(I)
          DP(I,3)=DPM1(3,1)*P1P2(I)+DPM1(3,2)*P1P3(I)+DPM1(3,3)*P1P4(I)
 40      CONTINUE
C
C        LES 3 COORDONNEES DES SOMMETS DE LA SURFACE NBTR
         NBT0 = NBTR - 1
         DO 90 J=1,NBSOMS
            DO 80 I=1,3
               DPARAF(I) = COOR3(I,J+NBSOM0) - RMCN(MNTREX(1,NBT0)+I-1)
 80         CONTINUE
            DO 86 I=1,3
C              LA TRANSLATION POUR REVENIR SUR LE POINT P1 ACTUEL
               A = RMCN(MNTREX(1,NBTR)+I-1)
               DO 84 K=1,3
                  A = A + DP(I,K) * DPARAF(K)
 84            CONTINUE
               COOR3(I,J+NBSOM1) = REAL( A )
 86         CONTINUE
 90      CONTINUE
C
C        LES TANGENTES DE LA TRANCHE NBTR
         IF( NBTGSU .GT. 0 ) THEN
            DO 99 J=1,NBTGSU
               DO 95 I=1,3
C                 LE DEPLACEMENT (TRANSFORMATION) DU VECTEUR TANGENTE
                  A = 0
                  DO 92 K=1,3
                     A = A + DP(I,K) * XYZTGV(K,J+NBTGS0)
 92               CONTINUE
                  XYZTGV(I,J+NBTGS1) = REAL( A )
 95            CONTINUE
 99         CONTINUE
C           PASSAGE A LA SURFACE SUIVANTE
            NBTGS0 = NBTGS1
            NBTGS1 = NBTGS1 + NBTGSU
         ENDIF
C
C        PASSAGE A LA SURFACE SUIVANTE
         NBSOM0 = NBSOM1
         NBSOM1 = NBSOM1 + NBSOMS
 100  CONTINUE
      RETURN
C
C     NUTYEX 2 : DECOUPAGE EQUIDISTANT DE LA COTE MAXIMALE
C     ----------------------------------------------------
 110  H  = VECTOR(3,NBTREX) / NBTREX
      Z1 = -H
      N  = 0
      DO 130 J=1,NBTRE1
           Z1 = Z1 + H
           IF( NUTYEX .EQ. 3 ) Z1 = Z( J )
           DO 120 I=1,NBSOMS
                N = N + 1
                COOR3( 1 , N ) = COOR2( 1 , I )
                COOR3( 2 , N ) = COOR2( 2 , I )
                COOR3( 3 , N ) = COOR2( 3 , I ) + Z1
 120       CONTINUE
 130  CONTINUE
      RETURN
C
C     NUTYEX 3 : DECOUPAGE IRREGULIER AVEC VECTOR( 3 , . )
C     ---------------------------------------------------
 200  Z(1) = 0.
      DO 210 I=1,NBTREX
         Z( I + 1 ) = ABS( VECTOR( 3 , I ) )
 210  CONTINUE
      CALL TRIREE( NBTRE1 , Z )
      IF( Z(2) .GT. 0. ) GOTO 230
      NBLGRC(NRERR) = 1
      KERR(1) = 'ERREUR VOE114: Z<=0'
      CALL LEREUR
  220 WRITE(IMPRIM,10220) (Z(I),I=1,NBTRE1)
10220 FORMAT('ERREUR VOE114:1-ERE COTE Z<0 OU Z(I)=Z(I+1).
     %NUTYEX FORCEE A 2'/(5G13.5))
      NUTYEX = 2
      GOTO 1
C
 230  DO 240 I=1,NBTREX
         IF( Z(I) .GE. Z(I+1) ) GOTO 220
 240  CONTINUE
      GOTO 110
C
C     NUTYEX 4 : VECTEUR DE TRANSLATION CONSTANT
C     NUTYEX 5 : VECTEUR DE TRANSLATION NON CONSTANT
C     ----------------------------------------------
 300  IF( ABS(VECTOR(1,1))+ABS(VECTOR(2,1))+ABS(VECTOR(3,1)) .GT. 0. )
     %    GOTO 310
      NBLGRC(NRERR) = 1
      KERR(1) = 'TRANSLATION D''UN VECTEUR NUL'
      CALL LEREUR
      IERR = 5
      RETURN
C
 310  V(1) = 0.
      V(2) = 0.
      V(3) = 0.
      N    = 0
      DO 350 J=1,NBTRE1
         DO 320 I=1,NBSOMS
            N = N + 1
            COOR3( 1 , N ) = COOR2( 1 , I ) + V( 1 )
            COOR3( 2 , N ) = COOR2( 2 , I ) + V( 2 )
            COOR3( 3 , N ) = COOR2( 3 , I ) + V( 3 )
 320     CONTINUE
         IF( J .EQ. NBTRE1 ) GOTO 350
         J1 = 1
         IF( NUTYEX .EQ. 5 ) J1 = J
         DO 330 I=1,3
C           ATTENTION AU RISQUE DE MAUVAISE ORIENTATION DES EF GENERES!
            V( I ) = V( I ) + VECTOR( I , J1 )
 330     CONTINUE
 350  CONTINUE
      RETURN
C
C     NUTYEX 6 : COORDONNEES FOURNIES PAR 3 FONCTIONS
C     -----------------------------------------------
 500  DO 590 J=0,NBTREX
C        TRANSFORMATION DES COORDONNEES DES SOMMETS DE LA SURFACE INITIALE
C        POUR LES AMENER SUR LA COUCHE J
         DPARAF(1) = J
         DO 580 K=1,3
C           LE NUMERO DE LA FONCTION EST NOFONC(K)
            DO 510 I=1,NBSOMS
               DPARAF(2) = COOR2(1,I)
               DPARAF(3) = COOR2(2,I)
               DPARAF(4) = COOR2(3,I)
               CALL FONVAL( NOFONC(K) , 4 , DPARAF , IERR , A )
C              TRANSFORMATION REEL DOUBLE EN REEL SIMPLE
               COOR3(K,I+J*NBSOMS) = REAL( A )
               IF( IERR .EQ. 0 ) THEN
C                 ERREUR A L'EXECUTION DEJA EXPLIQUEE DANS FONVAL
                  IERR = 1
                  RETURN
               ENDIF
 510        CONTINUE
 580     CONTINUE
 590  CONTINUE
      IERR = 0
      END
