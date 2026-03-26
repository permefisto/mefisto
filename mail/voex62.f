      SUBROUTINE VOEX62( NTLXVO, LADEFI, RADEFI,
     %                   NTEF6C, MNEF6C, NTST6C, MNST6C, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    GENERER LE MAILLAGE NON STRUCTURE EN 6-CUBES D'UNE BOITE
C -----    PAR HOMOTHETIE EN PROGRESSION GEOMETRIQUE DES 6-CUBES
C
C ENTREES:
C --------
C NTLXVO : NUMERO DU TABLEAU TS DU LEXIQUE DU VOLUME
C LADEFI : TABLEAU DE DEFINITION DES VOLUMES
C          CF '~td/d/a_volume__definition'
C
C SORTIES:
C --------
C NTEF6C : NUMERO      DU TMS 'NSEF' DES NUMEROS DES 6-CUBES
C MNEF6C : ADRESSE MCN DU TMS 'NSEF' DES NUMEROS DES 6-CUBES
C          CF '~td/d/a___nsef'
C NTST6C : NUMERO      DU TMS 'XYZSOMMET' DU MAILLAGE DES 6-CUBES
C MNST6C : ADRESSE MCN DU TMS 'XYZSOMMET' DU MAILLAGE DES 6-CUBES
C          CF '~td/d/a___xyzsommet'
C IERR   : 0 SI PAS D'ERREUR
C        > 0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: : ALAIN PERRONNET  TEXAS A & M UNIVERSITY         15 JULY 2005
C.......................................................................
      IMPLICIT INTEGER (W)
      COMMON / UNITES / LECTEU , IMPRIM , INTERA , NUNITE(29)
      include"./incl/gsmenu.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/a_volume__definition.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE       (MCN(1),RMCN(1))
      INTEGER           LADEFI(0:*)
      REAL              RADEFI(0:*)
      REAL              COIN(6,2)
      INTEGER           NOSOEF(64)
      include"./incl/nusc5c6.inc"
C
      IERR = 0
C
C     DIMENSION DE L'ESPACE DU 6-CUBES
C     ================================
      NBCOOR = 6
C
C     MIN A L'ORIGINE OU CUBE6 CENTRE?
C     NC6CUM = 0 : origine au minimum des coordonnees du cube6
C            =/0 : origine centree au centre du cube6
C     ========================================================
      NC6CUM = LADEFI(WUBE6M)
C
C     NOMBRE DE COUCHES OU HOMOTHETIES
C     ================================
      NC6CUB = LADEFI(WUBE6C)
      IF( NC6CUB .LT. 2 ) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:4),'(I4)') NC6CUB
         KERR(1)='NOMBRE INCORRECT (<2) DE COUCHES ='//KERR(MXLGER)(1:4)
         CALL LEREUR
         IERR = 3
         RETURN
      ENDIF
C
C     LARGEUR DU NOYAU DU 6-CUBES
C     ===========================
      CUB6LA = RADEFI(WUBE6L)
      IF( CUB6LA .LT. 0.0 ) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:14),'(E14.6)') CUB6LA
         KERR(1) =  'LARGEUR INCORRECTE DU NOYAU DU 6-CUBES='
     %           //  KERR(MXLGER)(1:4)
         CALL LEREUR
         IERR = 3
         RETURN
      ENDIF
C
C     LARGEUR TOTALE DU 6-CUBE DANS UNE DIRECTION
C     ===========================================
      CUB6DI = RADEFI(WUBE6D)
      IF( CUB6DI .LT. 0.0 ) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:14),'(E14.6)') CUB6DI
         KERR(1) =  'LARGEUR TOTALE INCORRECTE DU 6-CUBE='
     %           //  KERR(MXLGER)(1:4)
         CALL LEREUR
         IERR = 3
         RETURN
      ENDIF
C
C     RAISON GEOMETRIQUE D''HOMOTHETIE DANS LES 6 DIMENSIONS
C     ======================================================
      CUB6RG = REAL( ( CUB6DI / CUB6LA ) ** (1d0/NC6CUB) )
      IF( CUB6RG .LE. 1.0 ) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:14),'(E14.6)') CUB6RG
         KERR(1) =  'RAISON GEOMETRIQUE INCORRECTE (<=1.0) ='
     %           //  KERR(MXLGER)(1:14)
         CALL LEREUR
         IERR = 2
         RETURN
      ENDIF
C
C     GENERATION DES SOMMETS DE CES 6-CUBES CENTRES HOMOTHETIQUES
C     ===========================================================
C     CONSTRUCTION DU TABLEAU 'XYZSOMMET' DE CE VOLUME
      NBSOM =( 2 ** NBCOOR ) * NC6CUB
      CALL LXTNDC( NTLXVO, 'XYZSOMMET', 'MOTS' , WYZSOM+NBCOOR*NBSOM )
      CALL LXTSOU( NTLXVO, 'XYZSOMMET',  NTST6C, MNST6C )
C     LA DIMENSION DE L'ESPACE = NOMBRE DES COORDONNEES D'UN SOMMET
      MCN( MNST6C + WBCOOR ) = NBCOOR
C     LE NOMBRE DE SOMMETS
      MCN( MNST6C + WNBSOM ) = NBSOM
C
C     CALCUL DES NBCOOR COORDONNEES DES NBSOM SOMMETS
C     => 64 SOMMETS NOUVEAUX PAR COUCHE
      ARETE0 = CUB6LA / 2
      ARETE  = ARETE0
      MN     = MNST6C + WYZSOM
C
      IF( NC6CUM .NE. 0 ) THEN
C
C        MIN DES COORDONNEES DU CUBE 6 A l'ORIGINE
         X0 = 0.0
C
      ELSE
C
C        CUBE6 CENTRE
         X0 = ARETE0
         DO 5  NBC=1,NC6CUB-1
            X0 = X0 + ARETE0 * CUB6RG
 5       CONTINUE
C
      ENDIF
C
      DO 70 NBC=1,NC6CUB
C
         DO 60 N=-1,1,2
            DO 50 M=-1,1,2
               DO 40 L=-1,1,2
                  DO 30 K=-1,1,2
                     DO 20 J=-1,1,2
                        DO 10 I=-1,1,2
                           RMCN( MN     ) = X0 + I * ARETE
                           RMCN( MN + 1 ) = X0 + J * ARETE
                           RMCN( MN + 2 ) = X0 + K * ARETE
                           RMCN( MN + 3 ) = X0 + L * ARETE
                           RMCN( MN + 4 ) = X0 + M * ARETE
                           RMCN( MN + 5 ) = X0 + N * ARETE
                           MN  = MN + NBCOOR
 10                     CONTINUE
 20                  CONTINUE
 30               CONTINUE
 40            CONTINUE
 50         CONTINUE
 60      CONTINUE
C
C        L'HOMOTHETIE POUR PASSER A LA COUCHE SUIVANTE
CCC         ARETE = ARETE + ARETE0 * CUB6RG
         ARETE = ARETE * CUB6RG
         IF( NBC .EQ. NC6CUB-1 ) ARETE = CUB6DI
 70   CONTINUE
C
C     LA DATE DE CREATION
      CALL ECDATE( MCN(MNST6C) )
C     LE NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNST6C + WNBTGS ) = 0
      MCN( MNST6C + MOTVAR(6) ) = NONMTD( '~>>>XYZSOMMET' )
C
C     MIN ET MAX DES COORDONNEES
      CALL CADEXT( MNST6C, COIN )
C
C     GENERATION DES NUMEROS DES SOMMETS DE CHAQUE EF
C     ===============================================
C     CONSTRUCTION DU TABLEAU 'NSEF' DE CES 6-CUBES NON STRUCTURES
      NBSOEF = 2 ** NBCOOR
C     UNE COUCHE = NBFACES D'UN 6-CUBE => 12 6-CUBES
      NBEFOB = 1 + 2 * NBCOOR * (NC6CUB-1)
      MOTS   = WUSOEF + NBSOEF * NBEFOB
      CALL LXTNDC( NTLXVO, 'NSEF', 'ENTIER', MOTS   )
      CALL LXTSOU( NTLXVO, 'NSEF',  NTEF6C,  MNEF6C )
C     MISE A JOUR DU TABLEAU 'NSEF' DE CE VOLUME
C
C     TYPE DE L'OBJET : VOLUME BIEN QU'EN FAIT CE SOIT UN 6-CUBES!
      MCN( MNEF6C + WUTYOB ) = 4
C     LE TYPE INCONNU DE FERMETURE DU MAILLAGE
      MCN( MNEF6C + WUTFMA ) = -1
C
C     VARIABLE NBSOEF 'NOMBRE DE SOMMETS PAR EF' ENTIER
C     ( 4 : 'FACE' , 8 : 'CUBE', '64 : '6-CUBES' ) ;
      MCN( MNEF6C + WBSOEF ) = 64
C
C     VARIABLE NBTGEF 'NOMBRE DE TANGENTES PAR EF' ENTIER
C     ( 0 : '0 TG PAR EF' , 1 : '1 TG PAR SOMMET' , 2 : '2 TG PAR ARETE' ,
C       8 : '8 TG PAR FACE' ,  24 : '24 TG PAR CUBE' )  ;
      MCN( MNEF6C + WBTGEF ) = 0
      MCN( MNEF6C + WBEFAP ) = 0
      MCN( MNEF6C + WBEFTG ) = 0
C
C     VARIABLE NBEFOB 'NOMBRE DES EF DU PLSV'  ENTIER ;
      MCN( MNEF6C + WBEFOB ) = NBEFOB
C
C     VARIABLE NUTYMA 'NUMERO DE TYPE DU MAILLAGE' ENTIER
C     ( 0 : 'NON STRUCTURE' , ... )
      MCN( MNEF6C + WUTYMA ) = 0
C
C     TABLEAU  NUSOEF(1..NBSOEF,1..NBEFOB)
C    'NUMERO DES NBSOEF SOMMETS DES NBEFOB NSEF DE L''OBJET'
C
C     LA 1-ERE COUCHE = LE NOYAU = EF 1
      MN = MNEF6C + WUSOEF - 1
      DO 100 I=1,64
         MCN( MN + I ) = I
 100  CONTINUE
      MN = MN + 64
C
C     BOUCLE SUR LES 6-CUBES DES NC6CUB-1 COUCHES DU VOLUME
C     NO DES 32 SOMMETS DES 12 FACES D'UN 6-CUBE
      NBSOFA = 32
      NBFACE = 12
CCC      CALL SOFA6C( NBSOFA, NBFACE, NUS5C6C )
C
C     NUDS0 NUMERO DU DERNIER SOMMET DE LA COUCHE 0
C     NUDS1 NUMERO DU DERNIER SOMMET DE LA COUCHE 1
      NUDS0 = 0
      NUDS1 = 64
      DO 160 NBC=1,NC6CUB-1
C
C        LA COUCHE NBC+1 DES EF
         DO 150 NF=1,NBFACE
C           FACE COMPLEMENTAIRE DE LA FACE NF
            IF( MOD(NF,2) .EQ. 0 ) THEN
C              FACE PAIRE, SA COMPLEMENTAIRE EST LA PRECEDENTE
               NF1 = NF - 1
               NF2 = NF
            ELSE
C              FACE IMPAIRE, SA COMPLEMENTAIRE EST LA SUIVANTE
               NF1 = NF
               NF2 = NF + 1
            ENDIF
C
C           LE NO DES SOMMETS DE L'EF
            DO 110 I=1,32
               NS = NUS5C6C(I,NF)
               NOSOEF( NUS5C6C(I,NF1) ) = NUDS0 + NS
               NOSOEF( NUS5C6C(I,NF2) ) = NUDS1 + NS
 110        CONTINUE
C
C           TRANSFERT DANS NSEF
            DO 130 I=1,64
               MCN( MN + I ) = NOSOEF( I )
 130        CONTINUE
            MN = MN + 64
 150     CONTINUE
C
         NUDS0 = NUDS1
         NUDS1 = NUDS0 + 64
C
 160  CONTINUE
C
C     LA DATE DE CREATION
      CALL ECDATE( MCN(MNEF6C) )
C     LE NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNEF6C + MOTVAR(6) ) = NONMTD( '~>>>NSEF' )
C
C     RENUMEROTATION DES SOMMETS POUR REDUIRE LE PROFIL DE LA MATRICE
ccc      CALL RE6CUB( NTLXVO, NTEF6C, MNEF6C, MNST6C )
C     AUGMENTE LEGEREMENT LE PROFIL DE LA MATRICE! => SUPPRESSION
C
      RETURN
      END
