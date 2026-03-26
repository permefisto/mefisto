      SUBROUTINE AFVIPR( NUTYEL, NDIM,   NBNOVI, MNXYZN, NODDL,
     %                   TEMPS,  NTDLVP, NBVPAF, VXYZPN,
     %                   VITMIN, VITMAX, VITMOY, PREMIN, PREMAX, PREMOY)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    AFFICHER LE VECTEUR VITESSE+PRESSION CALCULEES AUX NOEUDS
C -----    et au TEMPS
C         (LES ERREURS AUX NOEUDS SONT AFFICHEES DANS afviprer.f)

C ENTREES:
C --------
C NUTYEL : 13 TRIANGLE  BREZZI FORTIN 2D
C          15 TRIANGLE  TAYLOR-HOOD   2D
C          19 TETRAEDRE BREZZI FORTIN 3D
C          20 TETRAEDRE TAYLOR-HOOD   3D
C NDIM   : DIMENSION 2 OU 3 DE L'ESPACE DE L'OBJET
C NBNOVI : NOMBRE DE NOEUDS SUPPORT DE LA VITESSE
C          SOMMETS + BARYCENTRES        des TETRAEDRES POUR BREZZI-FORTIN
C          SOMMETS + MILIEUX DES ARETES des TETRAEDRES POUR TAYLOR-HOOD
C MNXYZN : ADRESSE MCN DU TMS 'XYZNOEUD' DU MAILLAGE DE L'OBJET
C NODDL  : TABLEAU DU NUMERO DU DERNIER D.L. DE CHAQUE NOEUD

C TEMPS  : TEMPS DU VECTEUR VITESSE+PRESSION
C NTDLVP : NOMBRE TOTAL DE DEGRES DE LIBERTES EN VITESSES+PRESSIONS
C NBVPAF : NOMBRE DE NOEUDS DE VITESSE-PRESSION A AFFICHER
C VXYZPN : LA CARTE DES VITESSES-PRESSIONS AUX NOEUDS DU MAILLAGE AU TEMPS

C SORTIES:
C --------
C VITMIN : MINIMUM DE LA NORME DE LA VITESSE CALCULEE EN UN NOEUD de VXYZPN
C VITMAX : MAXIMUM DE LA NORME DE LA VITESSE CALCULEE EN UN NOEUD
C VITMOY : MOYENNE DE LA NORME DE LA VITESSE CALCULEE EN TOUS LES NOEUDS
C PREMIN : MINIMUM DE LA PRESSION CALCULEE EN UN NOEUD DU MAILLAGE
C PREMAX : MAXIMUM DE LA PRESSION CALCULEE EN UN NOEUD DU MAILLAGE
C PREMOY : MOYENNE DE LA PRESSION CALCULEE EN TOUS LES SOMMETS DU MAILLAGE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 2000
C MODIFS : ALAIN PERRONNET Laboratoire JL LIONS UPMC PARIS FEVRIER  2007
C MODIFS : ALAIN PERRONNET LJLL UPMC & Saint Pierre du Perray  Mars 2010
C MODIFS : ALAIN PERRONNET Veulettes sur mer                   Aout 2020
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/a___xyznoeud.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE      (MCN(1), RMCN(1), DMCN(1))
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)

      INTEGER           NODDL(0:NBNOVI), NVHIST(0:7)
      DOUBLE PRECISION  VXYZPN(NTDLVP)
      DOUBLE PRECISION  VITNOR, VITMIN, VITMAX, VITMOY,
     %                  PREMIN, PREMAX, PREMOY, PRECALC,
     %                  XP, YP, ZP, VNHIST(0:7)
      INTRINSIC         SQRT

C     L'ADRESSE MCN DES COORDONNEES DES NOEUDS VITESSE
C     ICI LE TMS XYZNOEUD A ETE ENRICHI DES BARYCENTRES DES EF
C     DANS LES CAS DES EF DE BREZZI-FORTIN
C     MN = MNXYZN + WYZNOE - 3

C     =========================================================
C     AFFICHAGE DES COMPOSANTES DE LA VITESSE ET DE LA PRESSION
C     AUX NOEUDS DU MAILLAGE
C     =========================================================
ccc      NOEUD1 = NBNOVI/2+1
ccc      NOEUD2 = NOEUD1 + NBVPAF
ccc      NOEUD1 = NBNOVI-NBVPAF
ccc      NOEUD2 = NBNOVI

      NOEUD1 = NBNOVI/2
      NOEUD2 = NOEUD1 + NBVPAF
      NOEUD2 = MIN( NOEUD2, NBNOVI )

C     SUPPRESSION DE L'AFFICHAGE SUIVANT. A RESTAURER EVENTUELLEMENT
C     --------------------------------------------------------------
      IF( NBNOVI .GT. 0 ) GOTO 40

      IF( LANGAG .EQ. 0 ) THEN

C        AFFICHAGE EN FRANCAIS
         PRINT 10019, TEMPS,NOEUD1,NOEUD2,NBNOVI,NTDLVP
10019    FORMAT(' Au Temps ',G14.7,': les VITESSES et PRESSIONS des ',
     %    'NOEUDS',I9,' a ',I9,'/',I9,' NOEUDS (Au total',I9,' DL):')

         DO I = NOEUD1, NOEUD2
            MN  = MNXYZN + WYZNOE - 3 + 3*I
            NDL = NODDL(I-1)
            ND  = NODDL(I) - NDL
            IF( NUTYEL .EQ. 19 .OR.  NUTYEL .EQ. 20 ) THEN

C              3D BF ou TH : VITESSE"PRESSION
               IF( ND .GT. NDIM ) THEN
                  PRINT 10023, I, (RMCN(MN+K),K=0,2),
     %                      (VXYZPN(NDL+K),K=1,ND)
               ELSE
                  PRINT 10033, I, (RMCN(MN+K),K=0,2),
     %                      (VXYZPN(NDL+K),K=1,ND)
               ENDIF

            ELSE
C
C              2D: NOMBRE DE DL EN LE NOEUD I
               IF( ND .GT. NDIM ) THEN
                  PRINT 10022, I, (RMCN(MN+K),K=0,2),
     %                 (VXYZPN(NDL+K),K=1,ND)
               ELSE
                  PRINT 10032, I, (RMCN(MN+K),K=0,2),
     %                 (VXYZPN(NDL+K),K=1,ND)
               ENDIF
            ENDIF
         ENDDO

10022 FORMAT(' Noeud',I9,': X=',G14.7,' Y=',G14.7,' Z=',G14.7,
     %'  VITESSE VX=',G14.7,' VY=',G14.7,' PRESSION=',G14.7)

10032 FORMAT(' Noeud',I9,': X=',G14.7,' Y=',G14.7,' Z=',G14.7,
     %'  VITESSE VX=',G14.7,' VY=',G14.7)

10023 FORMAT(' Noeud',I9,': X=',G14.7,' Y=',G14.7,' Z=',G14.7,
     %'  VITESSE VX=',G14.7,' VY=',G14.7,' VZ=',G14.7,
     %'  PRESSION=',G14.7)

10033 FORMAT(' Noeud',I9,': X=',G14.7,' Y=',G14.7,' Z=',G14.7,
     %'  VITESSE VX=',G14.7,' VY=',G14.7, ' VZ=',G14.7 )

      ELSE

C        ENGLISH PRINTING
         PRINT 20019, TEMPS,NOEUD1,NOEUD2,NBNOVI,NTDLVP
20019 FORMAT(' At Time ',G14.7,': VELOCITIES and PRESSURES of nodes',
     %I9,' to ',I9,'/',I9,' NODES (Total:',I9,' DoF):')

         DO I = NOEUD1, NOEUD2
            MN  = MNXYZN + WYZNOE - 3 + 3*I
            NDL = NODDL(I-1)
            ND  = NODDL(I) - NDL

C           AFFICHAGE DES COMPOSANTES DE LA VITESSE ET
C           EVENTUELLEMENT DE LA PRESSION
            IF( NDIM .EQ. 3 ) THEN

C              ELEMENT FINI 3D
               IF( ND .GT. NDIM ) THEN
C                 NOEUD AVEC PRESSION
                  PRINT 20023, I, (RMCN(MN+K),K=0,2),
     %                      (VXYZPN(NDL+K),K=1,ND)
               ELSE
C                 NOEUD SANS PRESSION
                  PRINT 20033, I, (RMCN(MN+K),K=0,2),
     %                      (VXYZPN(NDL+K),K=1,ND)
               ENDIF

            ELSE

C              ELEMENT FINI 2D
               IF( ND .GT. NDIM ) THEN
C                 NOEUD AVEC PRESSION
                  PRINT 20022, I, (RMCN(MN+K),K=0,2),
     %                 (VXYZPN(NDL+K),K=1,ND)
               ELSE
C                 NOEUD SANS PRESSION
                  PRINT 20032, I, (RMCN(MN+K),K=0,2),
     %                 (VXYZPN(NDL+K),K=1,ND)
               ENDIF
            ENDIF
         ENDDO

      ENDIF

20022 FORMAT(' NODE',I9,': X=',G14.7,' Y=',G14.7,' Z=',G14.7,
     %'  VELOCITY VX=',G14.7,' VY=',G14.7,'  PRESSURE=',G14.7)

20032 FORMAT(' NODE',I9,': X=',G14.7,' Y=',G14.7,' Z=',G14.7,
     %'  VELOCITY VX=',G14.7,' VY=',G14.7)

20023 FORMAT(' NODE',I9,': X=',G14.7,' Y=',G14.7,' Z=',G14.7,
     %'  VELOCITY VX=',G14.7,' VY=',G14.7,' VZ=',G14.7,
     %'  PRESSURE=',G14.7)

20033 FORMAT(' NODE',I9,': X=',G14.7,' Y=',G14.7,' Z=',G14.7,
     %'  VELOCITY VX=',G14.7,' VY=',G14.7,' VZ=',G14.7)

C     ====================================================
C     MIN MAX MOYENNE des |Vitesse| et PRESSIONS CALCULEES
C     ====================================================
 40   VITMOY =  0D0
      VITMIN =  1D100
      VITMAX = -1D100
      PREMIN =  1D100
      PREMAX = -1D100
      PREMOY = 0D0
      NBST   = 0

      DO I=1,NBNOVI
         NDL = NODDL(I-1)

C        MODULE DE LA VITESSE CALCULEE AU NOEUD I
         IF( NDIM .EQ. 2 ) THEN
            VITNOR = SQRT( VXYZPN(NDL+1)**2
     %                   + VXYZPN(NDL+2)**2 )
         ELSE
            VITNOR = SQRT( VXYZPN(NDL+1)**2
     %                   + VXYZPN(NDL+2)**2
     %                   + VXYZPN(NDL+3)**2 )
         ENDIF

C        NORME MOYENNE DU MODULE DE LA VITESSE CALCULEE
         VITMOY = VITMOY + VITNOR
C
C        NORME MIN et MAX DU MODULE DE LA VITESSE CALCULEE
         IF( VITNOR .LT. VITMIN ) VITMIN = VITNOR
         IF( VITNOR .GT. VITMAX ) THEN
            VITMAX = VITNOR
            NOVMAX = I
         ENDIF

C        PRESSION CALCULEE AU NOEUD I
         IF( NODDL(I)-NDL .GT. NDIM ) THEN
            NBST    = NBST + 1
            PRECALC = VXYZPN(NODDL(I))
            PREMOY  = PREMOY + PRECALC
C           PRECALC CALCULEE MIN et MAX EN UN NOEUD
            IF( PRECALC .LT. PREMIN ) PREMIN=PRECALC
            IF( PRECALC .GT. PREMAX ) PREMAX=PRECALC
         ENDIF
      ENDDO

C     MODULE MOYEN AUX NOEUDS DE LA VITESSE
      VITMOY = VITMOY / NBNOVI

C     PRESSION MOYENNE AUX SOMMETS
      PREMOY = PREMOY / NBST

C     COORDONNEES DU NOEUD DE VITESSE MAXIMALE
      MN = MNXYZN + WYZNOE - 3 + 3*NOVMAX
      XP = RMCN(MN)
      YP = RMCN(MN+1)
      ZP = RMCN(MN+2)

      IF( LANGAG .EQ. 0 ) THEN
         PRINT 10050,
     %         TEMPS, VITMOY, VITMAX, NOVMAX, XP, YP, ZP,
     %         TEMPS, PREMOY, PREMAX, PREMIN, PREMAX-PREMIN
      ELSE
         PRINT 20050,
     %         TEMPS, VITMOY, VITMAX, NOVMAX, XP, YP, ZP,
     %         TEMPS, PREMOY, PREMAX, PREMIN, PREMAX-PREMIN
      ENDIF

10050 FORMAT(' Au Temps ', G14.7,
     %' |VITESSE|Moyenne=',G14.7,
     %' |VITESSE|Max=',    G14.7,' au noeud',I9,' XYZ=',3G14.6/
     %       ' Au Temps ', G14.7,
     %' PRESSION Moyenne=',G14.7,' PRESSION Max=',G14.7,
     %' PRESSION Min=',    G14.7,' PRESSION Max-Min=',G14.7 )
C
20050 FORMAT(' At Time ',G14.7,
     %' |VELOCITY|Mean=',G14.7,
     %' |VELOCITY|Max=', G14.7,' at node',I9,' XYZ=',3G14.6/
     %       ' At Time ',G14.7,
     %' PRESSURE  Mean=',G14.7,' PRESSURE Max=',G14.7,
     %' PRESSURE  Min=', G14.7,' PRESSURE Max-Min=',G14.7 )


C     HISTOGRAMME DU NOMBRE DE TETRAEDRES POUR LE MODULE DE LA VITESSE
C     ----------------------------------------------------------------
      DO I=0,7
         NVHIST( I ) = 0
      ENDDO

      VNHIST(0) = 0.
      VNHIST(1) = VITMOY
      VNHIST(2) = MIN( VITMAX, 2*VITMOY )
      VNHIST(3) = MIN( VITMAX, 3*VITMOY )
      VNHIST(4) = VITMAX / 3
      VNHIST(5) = VITMAX / 2
      VNHIST(6) = VITMAX
      VNHIST(7) = VITMAX * 2

      DO I=1,NBNOVI
         NDL = NODDL(I-1)

C        MODULE DE LA VITESSE CALCULEE AU NOEUD I
         IF( NDIM .EQ. 2 ) THEN
            VITNOR = SQRT( VXYZPN(NDL+1)**2
     %                   + VXYZPN(NDL+2)**2 )
         ELSE
            VITNOR = SQRT( VXYZPN(NDL+1)**2
     %                   + VXYZPN(NDL+2)**2
     %                   + VXYZPN(NDL+3)**2 )
         ENDIF

C        PLACE de la NORME de la VITESSE dans l'HISTOGRAMME
         DO NV = 1, 6
            IF( VITNOR .LE. VNHIST(NV) ) GOTO 60
         ENDDO

 60      NVHIST( NV ) = NVHIST( NV ) + 1

      ENDDO

      DO NV = 1, 6
         IF( VNHIST(NV-1) .GT.  VNHIST(NV) ) GOTO 9000
         I = NVHIST(NV)
         P = ( 100.0 * I ) / NBNOVI
         IF( LANGAG .EQ. 0 ) THEN
            PRINT 10060, I, VNHIST(NV-1), VNHIST(NV), P
         ELSE
            PRINT 20060, I, VNHIST(NV-1), VNHIST(NV), P
         ENDIF
      ENDDO

10060 FORMAT(I9,' NOEUDS de VITESSE entre',G14.7,' et',G14.7,
     %          ' ->',F7.2,' %')
20060 FORMAT(I9,' NODES with a VELOCITY between',G14.7,' and',G14.7,
     %          ' ->',F7.2,' %')

 9000 RETURN
      END
