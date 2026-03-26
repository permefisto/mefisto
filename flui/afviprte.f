      SUBROUTINE AFVIPRTE( NUTYEL, NDIM,   NBNOVI, MNXYZN, NDDLNO,
     %                     TEMPS,  NTDLVP, NBVPAF, VXYZPN, TEMPER,
     %                     VITMIN, NOVMIN, VITMAX, NOVMAX, VITMOY,
     %                     PREMIN, NOPMIN, PREMAX, NOPMAX, PREMOY,
     %                     TEMMIN, NOTMIN, TEMMAX, NOTMAX, TEMMOY )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  AFFICHER au TEMPS LES VITESSES+PRESSIONS+TEMPERATURES CALCULEES
C -----
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
C NDDLNO : TABLEAU DU NUMERO DU DERNIER D.L. DE CHAQUE NOEUD

C TEMPS  : TEMPS DU VECTEUR VITESSES+PRESSIONS
C NTDLVP : NOMBRE TOTAL DE DEGRES DE LIBERTES EN VITESSES+PRESSIONS
C NBVPAF : NOMBRE DE NOEUDS DE VITESSE-PRESSION-TEMPERATURE A AFFICHER
C VXYZPN : LA CARTE DES VITESSES-PRESSIONS AUX NOEUDS DU MAILLAGE AU TEMPS
C TEMPER : LA CARTE DES TEMPERATURES       AUX NOEUDS DU MAILLAGE AU TEMPS

C SORTIES:
C --------
C VITMIN : MINIMUM DE LA NORME DE LA VITESSE CALCULEE EN UN NOEUD
C NOVMIN : NUMERO DU NOEUD DE VITESSE MINIMALE
C VITMAX : MAXIMUM DE LA NORME DE LA VITESSE CALCULEE EN UN NOEUD
C NOVMAX : NUMERO DU NOEUD DE VITESSE MAXIMALE
C VITMOY : MOYENNE DE LA NORME DE LA VITESSE CALCULEE DE TOUS LES NOEUDS
C PREMIN : MINIMUM DE LA PRESSION CALCULEE EN UN SOMMET DU MAILLAGE
C NOPMIN : NUMERO DU NOEUD DE PRESSION MINIMALE
C PREMAX : MAXIMUM DE LA PRESSION CALCULEE EN UN SOMMET DU MAILLAGE
C NOPMAX : NUMERO DU NOEUD DE VITESSE MAXIMALE
C PREMOY : MOYENNE DE LA PRESSION CALCULEE EN TOUS LES SOMMETS DU MAILLAGE
C TEMMIN : MINIMUM DE LA TEMPERATURE CALCULEE EN UN NOEUD DU MAILLAGE
C NOTMIN : NUMERO DU NOEUD DE TEMPERATURE MINIMALE
C TEMMAX : MAXIMUM DE LA TEMPERATURE CALCULEE EN UN NOEUD DU MAILLAGE
C NOTMAX : NUMERO DU NOEUD DE TEMPERATURE MAXIMALE
C TEMMOY : MOYENNE DE LA TEMPERATURE CALCULEE EN TOUS LES NOEUDS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Saint Pierre du Perray          Novembre 2022
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

      INTEGER           NDDLNO(0:NBNOVI), NVHIST(0:7)
      DOUBLE PRECISION  VXYZPN(NTDLVP), TEMPER(NBNOVI)
      DOUBLE PRECISION  VITMIN, VITMAX, VITMOY, VITNOR,
     %                  TEMMIN, TEMMAX, TEMMOY, TEM,
     %                  PREMIN, PREMAX, PREMOY, PRECALC,
     %                  VNHIST(0:7)
      INTRINSIC         SQRT

      PRINT *
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
            NDL = NDDLNO(I-1)
            ND  = NDDLNO(I) - NDL
            IF( NUTYEL .EQ. 19 .OR.  NUTYEL .EQ. 20 ) THEN

C              3D BF ou TH : VITESSE"PRESSION
               IF( ND .GT. NDIM ) THEN
                  PRINT 10023, I, (RMCN(MN+K),K=0,2),
     %                      (VXYZPN(NDL+K),K=1,ND), TEMPER(I)
               ELSE
                  PRINT 10033, I, (RMCN(MN+K),K=0,2),
     %                      (VXYZPN(NDL+K),K=1,ND), TEMPER(I)
               ENDIF

            ELSE
C
C              2D: NOMBRE DE DL EN LE NOEUD I
               IF( ND .GT. NDIM ) THEN
                  PRINT 10022, I, (RMCN(MN+K),K=0,2),
     %                 (VXYZPN(NDL+K),K=1,ND), TEMPER(I)
               ELSE
                  PRINT 10032, I, (RMCN(MN+K),K=0,2),
     %                 (VXYZPN(NDL+K),K=1,ND), TEMPER(I)
               ENDIF
            ENDIF
         ENDDO

10022 FORMAT(' Noeud',I9,': X=',G14.7,' Y=',G14.7,' Z=',G14.7,
     %'  VITESSE VX=',G14.7,' VY=',G14.7,' PRESSION=',G14.7,
     %' TEMPER=',G14.7)

10032 FORMAT(' Noeud',I9,': X=',G14.7,' Y=',G14.7,' Z=',G14.7,
     %'  VITESSE VX=',G14.7,' VY=',G14.7,' TEMPER=',G14.7)

10023 FORMAT(' Noeud',I9,': X=',G14.7,' Y=',G14.7,' Z=',G14.7,
     %'  VITESSE VX=',G14.7,' VY=',G14.7,' VZ=',G14.7,
     %'  PRESSION=',G14.7,' TEMPER=',G14.7)

10033 FORMAT(' Noeud',I9,': X=',G14.7,' Y=',G14.7,' Z=',G14.7,
     %'  VITESSE VX=',G14.7,' VY=',G14.7,' VZ=',G14.7,
     %' TEMPER=',G14.7)

      ELSE

C        ENGLISH PRINTING
         PRINT 20019, TEMPS,NOEUD1,NOEUD2,NBNOVI,NTDLVP
20019 FORMAT(' At Time ',G14.7,': VELOCITIES and PRESSURES of nodes',
     %I9,' to ',I9,'/',I9,' NODES (Total:',I9,' DoF):')

         DO I = NOEUD1, NOEUD2
            MN  = MNXYZN + WYZNOE - 3 + 3*I
            NDL = NDDLNO(I-1)
            ND  = NDDLNO(I) - NDL

C           AFFICHAGE DES COMPOSANTES DE LA VITESSE ET
C           EVENTUELLEMENT DE LA PRESSION
            IF( NDIM .EQ. 3 ) THEN

C              ELEMENT FINI 3D
               IF( ND .GT. NDIM ) THEN
C                 NOEUD AVEC PRESSION
                  PRINT 20023, I, (RMCN(MN+K),K=0,2),
     %                      (VXYZPN(NDL+K),K=1,ND), TEMPER(I)
               ELSE
C                 NOEUD SANS PRESSION
                  PRINT 20033, I, (RMCN(MN+K),K=0,2),
     %                      (VXYZPN(NDL+K),K=1,ND), TEMPER(I)
               ENDIF

            ELSE

C              ELEMENT FINI 2D
               IF( ND .GT. NDIM ) THEN
C                 NOEUD AVEC PRESSION
                  PRINT 20022, I, (RMCN(MN+K),K=0,2),
     %                 (VXYZPN(NDL+K),K=1,ND), TEMPER(I)
               ELSE
C                 NOEUD SANS PRESSION
                  PRINT 20032, I, (RMCN(MN+K),K=0,2),
     %                 (VXYZPN(NDL+K),K=1,ND), TEMPER(I)
               ENDIF
            ENDIF
         ENDDO

      ENDIF

20022 FORMAT(' NODE',I9,': X=',G14.7,' Y=',G14.7,' Z=',G14.7,
     %'  VELOCITY VX=',G14.7,' VY=',G14.7,'  PRESSURE=',G14.7,
     %' TEMPER=',G14.7 )

20032 FORMAT(' NODE',I9,': X=',G14.7,' Y=',G14.7,' Z=',G14.7,
     %'  VELOCITY VX=',G14.7,' VY=',G14.7,
     %' TEMPER=',G14.7 )

20023 FORMAT(' NODE',I9,': X=',G14.7,' Y=',G14.7,' Z=',G14.7,
     %'  VELOCITY VX=',G14.7,' VY=',G14.7,' VZ=',G14.7,
     %'  PRESSURE=',G14.7,' TEMPER=',G14.7 )

20033 FORMAT(' NODE',I9,': X=',G14.7,' Y=',G14.7,' Z=',G14.7,
     %'  VELOCITY VX=',G14.7,' VY=',G14.7,' VZ=',G14.7,
     %' TEMPER=',G14.7 )

C     ====================================================================
C     MIN MAX MOYENNE des |Vitesse| et PRESSIONS et TEMPERATURES CALCULEES
C     ====================================================================
 40   VITMOY =  0D0
      VITMIN =  1D100
      VITMAX = -1D100
      PREMIN =  1D100
      PREMAX = -1D100
      PREMOY =  0D0
      TEMMIN =  1D100
      TEMMAX = -1D100
      TEMMOY =  0D0
      NBST   = 0

      DO I=1,NBNOVI
         NDL = NDDLNO(I-1)

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
         IF( VITNOR .LT. VITMIN ) THEN
            VITMIN = VITNOR
            NOVMIN = I
         ENDIF
         IF( VITNOR .GT. VITMAX ) THEN
            VITMAX = VITNOR
            NOVMAX = I
         ENDIF

C        PRESSION CALCULEE AU NOEUD I
         IF( NDDLNO(I)-NDL .GT. NDIM ) THEN
            NBST    = NBST + 1
            PRECALC = VXYZPN( NDDLNO(I) )
            PREMOY  = PREMOY + PRECALC
C           PRECALC CALCULEE MIN et MAX EN UN NOEUD
            IF( PRECALC .LT. PREMIN ) THEN
               PREMIN=PRECALC
               NOPMIN = I
            ENDIF
            IF( PRECALC .GT. PREMAX ) THEN
               PREMAX=PRECALC
               NOPMAX = I
            ENDIF
         ENDIF

C        MIN et MAX de la TEMPERATURE CALCULEE
         TEM = TEMPER(I)
         TEMMOY = TEMMOY + TEM
         IF( TEM .LT. TEMMIN ) THEN
            TEMMIN = TEM
            NOTMIN = I
         ENDIF
         IF( TEM .GT. TEMMAX ) THEN
            TEMMAX = TEM
            NOTMAX = I
         ENDIF

      ENDDO

C     MODULE MOYEN AUX NOEUDS DE LA VITESSE
      VITMOY = VITMOY / NBNOVI

C     PRESSION MOYENNE AUX SOMMETS
      PREMOY = PREMOY / NBST

C     TEMPERATURE MOYENNE AUX NOEUDS
      TEMMOY = TEMMOY / NBNOVI

C     COORDONNEES DU NOEUD DE VITESSE MAXIMALE
      MN = MNXYZN + WYZNOE - 3 + 3*NOVMAX
      XPVMAX = RMCN(MN)
      YPVMAX = RMCN(MN+1)
      ZPVMAX = RMCN(MN+2)

C     COORDONNEES DU NOEUD DE PRESSION MAXIMALE
      MN = MNXYZN + WYZNOE - 3 + 3*NOPMAX
      XPPMAX = RMCN(MN)
      YPPMAX = RMCN(MN+1)
      ZPPMAX = RMCN(MN+2)

C     COORDONNEES DU NOEUD DE TEMPERATURE MAXIMALE
      MN = MNXYZN + WYZNOE - 3 + 3*NOTMAX
      XPTMAX = RMCN(MN)
      YPTMAX = RMCN(MN+1)
      ZPTMAX = RMCN(MN+2)

      IF( LANGAG .EQ. 0 ) THEN
         PRINT 10050,
     %         TEMPS, TEMMIN, TEMMOY, TEMMAX,
     %                NOTMAX, XPTMAX, YPTMAX, ZPTMAX,
     %         TEMPS, PREMIN, PREMOY, PREMAX,
     %                NOPMAX, XPPMAX, YPPMAX, ZPPMAX,
     %         TEMPS, VITMIN, VITMOY, VITMAX,
     %                NOVMAX, XPVMAX, YPVMAX, ZPVMAX
      ELSE
         PRINT 20050,
     %         TEMPS, TEMMIN, TEMMOY, TEMMAX,
     %                NOTMAX, XPTMAX, YPTMAX, ZPTMAX,
     %         TEMPS, PREMIN, PREMOY, PREMAX,
     %                NOPMAX, XPPMAX, YPPMAX, ZPPMAX,
     %         TEMPS, VITMIN, VITMOY, VITMAX,
     %                NOVMAX, XPVMAX, YPVMAX, ZPVMAX
      ENDIF

10050 FORMAT(' Au Temps ', G14.7,
     %' TEMPERATUREMin=',G14.7,' TEMPERATUREMoy=',G14.7,
     %' TEMPERATUREMax=',G14.7,' au noeud',I9,' XYZ=',3G13.6/
     %' Au Temps ', G14.7,
     %' PRESSION   Min=',G14.7,'  PRESSION  Moy=',G14.7,
     %' PRESSION   Max=',G14.7,' au noeud',I9,' XYZ=',3G13.6/
     %' Au Temps ', G14.7,
     %' |VITESSE|  Min=',G14.7,' |VITESSE|  Moy=',G14.7,
     %' |VITESSE|  Max=',G14.7,' au noeud',I9,' XYZ=',3G13.6)

20050 FORMAT(' At Time ',G14.7,
     %' TEMPERATUREMin=',G14.7,' TEMPERATUREMean=',G14.7,
     %' TEMPERATUREMax=',G14.7,' at node',I9,' XYZ=',3G13.6/
     %' At Time ',G14.7,
     %'  PRESSURE  Min=',G14.7,'  PRESSURE  Mean=',G14.7,
     %'  PRESSURE  Max=',G14.7,' at node',I9,' XYZ=',3G13.6/
     %' At Time ',G14.7,
     %' |VELOCITY| Min=',G14.7,' |VELOCITY| Mean=',G14.7,
     %' |VELOCITY| Max=',G14.7,' at node',I9,' XYZ=',3G13.6)


C     HISTOGRAMME DU NOMBRE D'EF POUR LE MODULE DE LA VITESSE
C     -------------------------------------------------------
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
         NDL = NDDLNO(I-1)

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
         IF( NV .GT. 1 .AND. VNHIST(NV-1) .GE.  VNHIST(NV) ) GOTO 9000
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
