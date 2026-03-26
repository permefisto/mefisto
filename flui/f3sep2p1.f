      SUBROUTINE F3SEP2P1( DT,  X, TP2P2,
     %                     NOOBSF, NUMISU, NUMASU, LTDESU,
     %                     NOVOLU, NUMIVO, NUMAVO, LTDEVO,
     %                     BE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DU SECOND MEMBRE DU TETRAEDRE TAYLOR HOOD
C -----    P2 CONTINU POUR LES 3 COMPOSANTES DE LA VITESSE
C          P1 CONTINU POUR LA PRESSION
C          INTERPOLATION DES FORCES AUX 10 POINTS DU TETRAEDRE
C
C ENTREES:
C --------
C DT     : PAS DE TEMPS DU SCHEMA D'EULER D'INTEGRATION EN TEMPS
C X      : 3 COORDONNEES DES 10 POINTS DU TETRAEDRE
C TP2P2  : TP2P2(i,j) = Integrale sur e chapeau P2i P2j dx dy dz
C
C NOOBSU : NUMERO DES OBJETS SURFACES DES FACES DE L'ELEMENT FINI
C FORCE  : FORCE IMPOSEE SUR UNE SURFACE FORCE(2) REMPLI DANS CE SP
C
C NOOBSF : NUMERO DE L'OBJET SURFACE DES 4 FACES DE L'EF
C FORCE  : EFFORTS DE SURFACE (2)  REMPLI DANS CE SOUS PROGRAMME
C NUMISU : NUMERO MINIMAL DES OBJETS SURFACES
C NUMASU : NUMERO MAXIMAL DES OBJETS SURFACES
C LTDESU : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES DU FLUIDE
C          DES OBJETS SURFACES
C
C NOVOLU : NUMERO DE L'OBJET VOLUME DE CET ELEMENT FINI
C FORCE  : EFFORTS DE VOLUME REMPLI DANS CE SOUS PROGRAMME
C NUMIVO : NUMERO MINIMAL DES OBJETS VOLUMES
C NUMAVO : NUMERO MAXIMAL DES OBJETS VOLUMES
C LTDEVO : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES VOLUMIQUES DU FLUIDE
C
C SORTIES:
C --------
C BE     : BE(34) LE SECOND MEMBRE ELEMENTAIRE (10+10+10+4)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY      Mai 2010
C23456---------------------------------------------------------------012
C     NO TYPE TETRAEDRE TAYLOR-HOOD
      INTEGER    NUTYEL
      PARAMETER (NUTYEL=20)
      include"./incl/donflu.inc"
      include"./incl/ctemps.inc"
      include"./incl/ponoel.inc"
      REAL              DT, X(10,3)
      DOUBLE PRECISION  TP2P2(10,10), BE(34)
      DOUBLE PRECISION  FORCE(3,10), VITANG(3,2)
      INTEGER           LTDEVO( 1:MXDOFL, NUMIVO:NUMAVO )
      INTEGER           LTDESU( 1:MXDOFL, NUMISU:NUMASU )
      INTEGER           NOOBSF(4), NUMISU, NUMASU,
     %                  NOVOLU,    NUMIVO, NUMAVO
C
      INTEGER           NOPOF(6), I, J, K, L, M, N, MN
      DOUBLE PRECISION  DETM33, ABS, DELTA, A, B, C, D, E, F, G
      DOUBLE PRECISION  XD, YD, ZD
      DOUBLE PRECISION  GL(3), DGL(2,3), DGLN, VN(3)
      REAL              TEMPS0
C
C     CALCUL DU JACOBIEN DU TETRAEDRE
C     ===============================
      XD = X(1,1)
      YD = X(1,2)
      ZD = X(1,3)
      DELTA = ABS( DETM33( X(2,1)-XD, X(3,1)-XD, X(4,1)-XD,
     %                     X(2,2)-YD, X(3,2)-YD, X(4,2)-YD,
     %                     X(2,3)-ZD, X(3,3)-ZD, X(4,3)-ZD ) )
C
C     1) LES EFFORTS VOLUMIQUES: F Omega est interpolee P2 Lagrange
C     =============================================================
      IF( LTDEVO(LPFORC,NOVOLU) .GT. 0 ) THEN
C
C        LES COEFFICIENTS DE LA MATRICE DE MASSE
C        ---------------------------------------
         A =  DELTA / 420.D0
         G =  DELTA / 315.D0
         B =  G * 4.D0
         C =  DELTA / 2520.D0
         D = -DELTA /  630.D0
         E = -A
         F =  G * 2.D0
C
C        RECHERCHE DE LA VALEUR DES EFFORTS VOLUMIQUES AUX 10 NOEUDS
         DO K=1,10
            XD = X(K,1)
            YD = X(K,2)
            ZD = X(K,3)
            CALL REFORC(4, NOVOLU, 3, XD,YD,ZD, 0D0,0D0,0D0,
     %                  LTDEVO(LPFORC,NOVOLU), FORCE(1,K) )
C           FORCE(1:3,K) EST INITIALISE
         ENDDO
C
         M = 0
         DO K=1,3
C
            BE(M+1) = A * FORCE(K,1)
     %         + C * ( FORCE(K,2) + FORCE(K,3) + FORCE(K,4)  )
     %         + D * ( FORCE(K,5) + FORCE(K,7) + FORCE(K,8)  )
     %         + E * ( FORCE(K,6) + FORCE(K,9) + FORCE(K,10) )
C
            BE(M+2) = A * FORCE(K,2)
     %         + C * ( FORCE(K,1) + FORCE(K,3) + FORCE(K,4)  )
     %         + D * ( FORCE(K,5) + FORCE(K,6) + FORCE(K,9)  )
     %         + E * ( FORCE(K,7) + FORCE(K,8) + FORCE(K,10) )
C
            BE(M+3) = A * FORCE(K,3)
     %         + C * ( FORCE(K,1) + FORCE(K,2) + FORCE(K,4)  )
     %         + D * ( FORCE(K,6) + FORCE(K,7) + FORCE(K,10) )
     %         + E * ( FORCE(K,5) + FORCE(K,8) + FORCE(K,9)  )
C
            BE(M+4) = A * FORCE(K,4)
     %         + C * ( FORCE(K,1) + FORCE(K,2) + FORCE(K,3)  )
     %         + D * ( FORCE(K,8) + FORCE(K,9) + FORCE(K,10) )
     %         + E * ( FORCE(K,5) + FORCE(K,6) + FORCE(K,7)  )
C
            BE(M+5) = B * FORCE(K,5)
     %         + D * ( FORCE(K,1) + FORCE(K,2) )
     %         + E * ( FORCE(K,3) + FORCE(K,4) )
     %         + F * ( FORCE(K,6) + FORCE(K,7)+ FORCE(K,8)+ FORCE(K,9) )
     %         + G *   FORCE(K,10)
C
            BE(M+6) = B * FORCE(K,6)
     %         + D * ( FORCE(K,2) + FORCE(K,3) )
     %         + E * ( FORCE(K,1) + FORCE(K,4) )
     %         + F * ( FORCE(K,5) + FORCE(K,7)+ FORCE(K,9)+ FORCE(K,10))
     %         + G *   FORCE(K,8)
C
            BE(M+7) = B * FORCE(K,7)
     %         + D * ( FORCE(K,1) + FORCE(K,3) )
     %         + E * ( FORCE(K,2) + FORCE(K,4) )
     %         + F * ( FORCE(K,5) + FORCE(K,6)+ FORCE(K,8)+ FORCE(K,10))
     %         + G *   FORCE(K,9)
C
            BE(M+8) = B * FORCE(K,8)
     %         + D * ( FORCE(K,1) + FORCE(K,4) )
     %         + E * ( FORCE(K,2) + FORCE(K,3) )
     %         + F * ( FORCE(K,5) + FORCE(K,7)+ FORCE(K,9)+ FORCE(K,10))
     %         + G *   FORCE(K,6)
C
            BE(M+9) = B * FORCE(K,9)
     %         + D * ( FORCE(K,2) + FORCE(K,4) )
     %         + E * ( FORCE(K,1) + FORCE(K,3) )
     %         + F * ( FORCE(K,5) + FORCE(K,6)+ FORCE(K,8)+ FORCE(K,10))
     %         + G *   FORCE(K,7)
C
            BE(M+10) = B * FORCE(K,10)
     %         + D * ( FORCE(K,3) + FORCE(K,4) )
     %         + E * ( FORCE(K,1) + FORCE(K,2) )
     %         + F * ( FORCE(K,6) + FORCE(K,7)+ FORCE(K,8)+ FORCE(K,9))
     %         + G *   FORCE(K,5)
C
            M = M + 10
         ENDDO
C
         DO K=31,34
            BE(K) = 0D0
         ENDDO
C
      ELSE
C
C        MISE A ZERO DU SECOND MEMBRE
         DO K=1,34
            BE(K) = 0D0
         ENDDO
C
      ENDIF
C
C     VECTEUR(3) VITESSE ANGULAIRE AUX 10 NOEUDS DU TETRAEDRE
C     -------------------------------------------------------
      IF( LTDEVO(LPVIAN,NOVOLU) .GT. 0 ) THEN
C        CALCULS AU TEMPS-DT et TEMPS
         TEMPS0 = TEMPS
         DO J=1,10
C
C           AU NOEUD J DU TETRAEDRE
            TEMPS = TEMPS0 - DT
            XD = X(J,1)
            YD = X(J,2)
            ZD = X(J,3)
            DO K=1,2
               CALL REVIAN( 4, NOVOLU, XD, YD, ZD,
     %                      LTDEVO(LPVIAN,NOVOLU), VITANG(1,K) )
               TEMPS = TEMPS0
            ENDDO
C
C           CALCUL DES FORCES DUES A LA ROTATION
C           AU SOMMET J DU TETRAEDRE   MODIFIE LE 22/2/2011
C           - dOmega/dt x r  -g vecteur z  - Omega x ( Omega x r )
            FORCE(1,J) = ( - (VITANG(2,2)-VITANG(2,1))/DT * ZD
     %                     + (VITANG(3,2)-VITANG(3,1))/DT * YD
     %           + VITANG(3,2) * ( VITANG(3,2) * XD - VITANG(1,2) * ZD )
     %           - VITANG(2,2) * ( VITANG(1,2) * YD - VITANG(2,2) * XD )
     %                   ) * DELTA
C
            FORCE(2,J) = ( - (VITANG(3,2)-VITANG(3,1))/DT * XD
     %                     + (VITANG(1,2)-VITANG(1,1))/DT * ZD
     %           - VITANG(3,2) * ( VITANG(2,2) * ZD - VITANG(3,2) * YD )
     %           + VITANG(1,2) * ( VITANG(1,2) * YD - VITANG(2,2) * XD )
     %                   ) * DELTA
C
            FORCE(3,J) = ( - (VITANG(1,2)-VITANG(1,1))/DT * YD
     %                     + (VITANG(2,2)-VITANG(2,1))/DT * XD
     %           + VITANG(2,2) * ( VITANG(2,2) * ZD - VITANG(3,2) * YD )
     %           - VITANG(1,2) * ( VITANG(3,2) * XD - VITANG(1,2) * ZD )
     %                   ) * DELTA
         ENDDO
C
C        COEFFICIENTS DU SECOND MEMBRE LIES AUX EFFORTS DE ROTATION
         DO I=1,10
            DO J=1,10
               D = TP2P2(I,J)
               BE(I   ) = BE(I   ) + D * FORCE(1,J)
               BE(I+10) = BE(I+10) + D * FORCE(2,J)
               BE(I+20) = BE(I+20) + D * FORCE(3,J)
            ENDDO
         ENDDO
      ENDIF
C
C     2) CONTRIBUTION DES EFFORTS SURFACIQUES SUR LES 4 FACES
C     =======================================================
      DO 100 K=1,4
C
C        NO DE SURFACE DE LA FACE K
         NOOB = NOOBSF(K)
         IF( NOOB .GT. 0 ) THEN
C
C           LA FACE EST SUR UNE SURFACE SUPPORT DE FORCE?
            MN = LTDESU( LPFORC, NOOB )
            IF( MN .GT. 0 ) THEN
C
C              UN TABLEAU FORCE EXISTE POUR CETTE FACE K
C              CALCUL DE LA CONTRIBUTION DE LA FACE K A BE
C              ...........................................
C              RECHERCHE DU NUMERO DES POINTS=NOEUDS DE LA FACE K
C              TETRAEDRE P2 NUTYEL=20
               CALL ELNOFA( NUTYEL, K, N, NOPOF )
C              NOPOF(I) NUMERO DU NOEUD I DE LA FACE TRIANGULAIRE K
C              LES 3 PREMIERS NOEUDS SONT LES SOMMETS DE LA FACE K
C              SUPPOSEE DROITE C'EST A DIRE INTERPOLEE P1
               N = NOPOF(1)
               XD = X(N,1)
               YD = X(N,2)
               ZD = X(N,3)
C
               N = NOPOF(2)
               DGL(1,1) = X(N,1) - XD
               DGL(1,2) = X(N,2) - YD
               DGL(1,3) = X(N,3) - ZD
C
               N = NOPOF(3)
               DGL(2,1) = X(N,1) - XD
               DGL(2,2) = X(N,2) - YD
               DGL(2,3) = X(N,3) - ZD
C
C              CALCUL DU VECTEUR NORMAL A LA FACE K SUPPOSEE PLANE
               CALL VECNOR( DGL, DGLN, VN )
C
               DO 40 L=1,6
C                 CALCUL DES FORCES EN CE NOEUD L DE LA FACE K
                  N = NOPOF(L)
                  GL(1) = X(N,1)
                  GL(2) = X(N,2)
                  GL(3) = X(N,3)
                  CALL REFORC( 3, NOOB, 3,
     &                         GL(1), GL(2), GL(3), VN(1), VN(2), VN(3),
     &                         LTDESU(LPFORC,NOOB), FORCE(1,L) )
 40            CONTINUE
C
C              PRODUIT DE LA MATRICE DE MASSE SURFACIQUE * FORCEm
C              ==================================================
C              NUMERO ELEMENTAIRE DU NOEUD 1 DE LA FACE K
               NU = NOPOF(1)
C              LES 3 COMPOSANTES DE LA FORCE ELEMENTAIRE
               N = 0
               DO 51 M=1,3
C                 LE COEFFICIENT DU SECOND MEMBRE ELEMENTAIRE DE LA COMPOSANTE M
                  BE(NU+N) = BE(NU+N) + DGLN *
     %              ( FORCE(M,1)/60D0 - (FORCE(M,2)+FORCE(M,3))/360D0
     %              - FORCE(M,5)/90D0 )
                  N = N + 10
 51            CONTINUE
C
C              NUMERO ELEMENTAIRE DU NOEUD 2 DE LA FACE K
               NU = NOPOF(2)
C              LES 3 COMPOSANTES DE LA FORCE ELEMENTAIRE
               N = 0
               DO 52 M=1,3
C                 LE COEFFICIENT DU SECOND MEMBRE ELEMENTAIRE DE LA COMPOSANTE M
                  BE(NU+N) = BE(NU+N) + DGLN *
     %              ( FORCE(M,2)/60D0 - (FORCE(M,1)+FORCE(M,3))/360D0
     %              - FORCE(M,6)/90D0 )
                  N = N + 10
 52            CONTINUE
C
C              NUMERO ELEMENTAIRE DU NOEUD 3 DE LA FACE K
               NU = NOPOF(3)
C              LES 3 COMPOSANTES DE LA FORCE ELEMENTAIRE
               N = 0
               DO 53 M=1,3
C                 LE COEFFICIENT DU SECOND MEMBRE ELEMENTAIRE DE LA COMPOSANTE M
                  BE(NU+N) = BE(NU+N) + DGLN *
     %              ( FORCE(M,3)/60D0 - (FORCE(M,1)+FORCE(M,2))/360D0
     %              - FORCE(M,4)/90D0 )
                  N = N + 10
 53            CONTINUE
C
C
C              NUMERO ELEMENTAIRE DU NOEUD 4 DE LA FACE K
               NU = NOPOF(4)
C              LES 3 COMPOSANTES DE LA FORCE ELEMENTAIRE
               N = 0
               DO 54 M=1,3
C                 LE COEFFICIENT DU SECOND MEMBRE ELEMENTAIRE DE LA COMPOSANTE M
                  BE(NU+N) = BE(NU+N) + DGLN *
     %               (-FORCE(M,3)/90D0 + 4D0/45D0 * FORCE(M,4)
     %               + 2D0/45D0 * ( FORCE(M,5) + FORCE(M,6) ) )
                  N = N + 10
 54            CONTINUE
C
C              NUMERO ELEMENTAIRE DU NOEUD 5 DE LA FACE K
               NU = NOPOF(5)
C              LES 3 COMPOSANTES DE LA FORCE ELEMENTAIRE
               N = 0
               DO 55 M=1,3
C                 LE COEFFICIENT DU SECOND MEMBRE ELEMENTAIRE DE LA COMPOSANTE M
                  BE(NU+N) = BE(NU+N) + DGLN *
     %               (-FORCE(M,1)/90D0 + 4D0/45D0 * FORCE(M,5)
     %               + 2D0/45D0 * ( FORCE(M,4) + FORCE(M,6) ) )
                  N = N + 10
 55            CONTINUE
C
C              NUMERO ELEMENTAIRE DU NOEUD 6 DE LA FACE K
               NU = NOPOF(6)
C              LES 3 COMPOSANTES DE LA FORCE ELEMENTAIRE
               N = 0
               DO 56 M=1,3
C                 LE COEFFICIENT DU SECOND MEMBRE ELEMENTAIRE DE LA COMPOSANTE M
                  BE(NU+N) = BE(NU+N) + DGLN *
     %               (-FORCE(M,2)/90D0 + 4D0/45D0 * FORCE(M,6)
     %               + 2D0/45D0 * ( FORCE(M,4) + FORCE(M,5) ) )
                  N = N + 10
 56            CONTINUE
            ENDIF
         ENDIF
C
 100  CONTINUE
C
      RETURN
      END
