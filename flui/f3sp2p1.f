      SUBROUTINE F3SP2P1( X,
     %                    NOOBSF, NUMISU, NUMASU, LTDESU,
     %                    NOVOLU, NUMIVO, NUMAVO, LTDEVO,
     %                    BE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DU SECOND MEMBRE DU TETRAEDRE TAYLOR HOOD
C -----    P2 CONTINU POUR LES 3 COMPOSANTES DE LA VITESSE
C          P1 CONTINU POUR LA PRESSION
C          INTERPOLATION DES FORCES AUX 10 POINTS DU TETRAEDRE
C
C ENTREES:
C --------
C X      : 3 COORDONNEES DES 10 POINTS DU TETRAEDRE
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
C AUTEUR: ALAIN PERRONNET Laboratoire J-L.LIONS UPMC Paris     Juin 2007
C AUTEUR: ALAIN PERRONNET Laboratoire J-L.LIONS UPMC Paris Novembre 2008
C23456---------------------------------------------------------------012
C     NO TYPE TETRAEDRE TAYLOR-HOOD
      PARAMETER (NUTYEL=20)
      include"./incl/donflu.inc"
      include"./incl/ponoel.inc"
      REAL              X(10,3)
      DOUBLE PRECISION  FORC(3,10), BE(34)
      INTEGER           NOOBSF(4), NUMISU, NUMASU,
     %                  NOVOLU,    NUMIVO, NUMAVO
      INTEGER           LTDEVO( 1:MXDOFL, NUMIVO:NUMAVO )
      INTEGER           LTDESU( 1:MXDOFL, NUMISU:NUMASU )
C
      INTEGER           K, L, M, MN, N, NOPOF(6), NOOB, NU
      DOUBLE PRECISION  DETM33, ABS, DELTA, A, B, C, D, E, F, G
      DOUBLE PRECISION  XD, YD, ZD
      DOUBLE PRECISION  FORCE(3,6)
      DOUBLE PRECISION  GL(3), DGL(2,3), DGLN, VN(3)
C
C     LES EFFORTS VOLUMIQUES: F Omega est interpolee P2 Lagrange
C     ==========================================================
      IF( LTDEVO(LPFORC,NOVOLU) .GT. 0 ) THEN
C
C        CALCUL DU JACOBIEN
         XD = X(1,1)
         YD = X(1,2)
         ZD = X(1,3)
         DELTA = ABS( DETM33( X(2,1)-XD, X(3,1)-XD, X(4,1)-XD,
     %                        X(2,2)-YD, X(3,2)-YD, X(4,2)-YD,
     %                        X(2,3)-ZD, X(3,3)-ZD, X(4,3)-ZD ) )
C
C        LES COEFFICIENTS DE LA MATRICE DE MASSE
C        =======================================
         A =  DELTA / 420.D0
         G =  DELTA / 315.D0
         B =  G * 4.D0
         C =  DELTA / 2520.D0
         D = -DELTA /  630.D0
         E = -A
         F =  G * 2.D0
C
C        RECHERCHE DE LA VALEUR DES EFFORTS VOLUMIQUES AUX 10 POINTS
C        ===========================================================
         DO K=1,10
            XD = X(K,1)
            YD = X(K,2)
            ZD = X(K,3)
            CALL REFORC(4, NOVOLU, 3, XD,YD,ZD, 0D0,0D0,0D0,
     %                  LTDEVO(LPFORC,NOVOLU), FORC(1,K) )
C           FORC(1:3,K) EST INITIALISE
         ENDDO
C
         M = 0
         DO K=1,3
C
            BE(M+1) = A * FORC(K,1)
     %         + C * ( FORC(K,2) + FORC(K,3) + FORC(K,4)  )
     %         + D * ( FORC(K,5) + FORC(K,7) + FORC(K,8)  )
     %         + E * ( FORC(K,6) + FORC(K,9) + FORC(K,10) )
C
            BE(M+2) = A * FORC(K,2)
     %         + C * ( FORC(K,1) + FORC(K,3) + FORC(K,4)  )
     %         + D * ( FORC(K,5) + FORC(K,6) + FORC(K,9)  )
     %         + E * ( FORC(K,7) + FORC(K,8) + FORC(K,10) )
C
            BE(M+3) = A * FORC(K,3)
     %         + C * ( FORC(K,1) + FORC(K,2) + FORC(K,4)  )
     %         + D * ( FORC(K,6) + FORC(K,7) + FORC(K,10) )
     %         + E * ( FORC(K,5) + FORC(K,8) + FORC(K,9)  )
C
            BE(M+4) = A * FORC(K,4)
     %         + C * ( FORC(K,1) + FORC(K,2) + FORC(K,3)  )
     %         + D * ( FORC(K,8) + FORC(K,9) + FORC(K,10) )
     %         + E * ( FORC(K,5) + FORC(K,6) + FORC(K,7)  )
C
            BE(M+5) = B * FORC(K,5)
     %         + D * ( FORC(K,1) + FORC(K,2) )
     %         + E * ( FORC(K,3) + FORC(K,4) )
     %         + F * ( FORC(K,6) + FORC(K,7) + FORC(K,8) + FORC(K,9) )
     %         + G *   FORC(K,10)
C
            BE(M+6) = B * FORC(K,6)
     %         + D * ( FORC(K,2) + FORC(K,3) )
     %         + E * ( FORC(K,1) + FORC(K,4) )
     %         + F * ( FORC(K,5) + FORC(K,7) + FORC(K,9) + FORC(K,10) )
     %         + G *   FORC(K,8)
C
            BE(M+7) = B * FORC(K,7)
     %         + D * ( FORC(K,1) + FORC(K,3) )
     %         + E * ( FORC(K,2) + FORC(K,4) )
     %         + F * ( FORC(K,5) + FORC(K,6) + FORC(K,8) + FORC(K,10) )
     %         + G *   FORC(K,9)
C
            BE(M+8) = B * FORC(K,8)
     %         + D * ( FORC(K,1) + FORC(K,4) )
     %         + E * ( FORC(K,2) + FORC(K,3) )
     %         + F * ( FORC(K,5) + FORC(K,7) + FORC(K,9) + FORC(K,10) )
     %         + G *   FORC(K,6)
C
            BE(M+9) = B * FORC(K,9)
     %         + D * ( FORC(K,2) + FORC(K,4) )
     %         + E * ( FORC(K,1) + FORC(K,3) )
     %         + F * ( FORC(K,5) + FORC(K,6) + FORC(K,8) + FORC(K,10) )
     %         + G *   FORC(K,7)
C
            BE(M+10) = B * FORC(K,10)
     %         + D * ( FORC(K,3) + FORC(K,4) )
     %         + E * ( FORC(K,1) + FORC(K,2) )
     %         + F * ( FORC(K,6) + FORC(K,7) + FORC(K,8) + FORC(K,9) )
     %         + G *   FORC(K,5)
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
C     CONTRIBUTION DES EFFORTS SURFACIQUES
C     ====================================
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
