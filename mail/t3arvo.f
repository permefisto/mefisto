      SUBROUTINE T3ARVO( NMVOLU )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACER TOUTES LES ARETES (INTERNES ET FRONTALIERES)
C -----    D'UN VOLUME DE NOM NMVOLU
C
C ENTREES:
C --------
C NMVOLU : NOM DU VOLUME
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC PARIS     OCTOBRE 1997
C MODIFS : ALAIN PERRONNET  TEXAS A M UNIVERSITY & LJLL     OCTOBRE 2005
C ...................................................................012
      IMPLICIT INTEGER (W)
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/gsmenu.inc"
      include"./incl/mecoit.inc"
      include"./incl/sotgar.inc"
      include"./incl/sotgfc.inc"
C
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
      CHARACTER*(*)     NMVOLU
      CHARACTER*24      NMSOMM
      REAL              XYZTG(3,2)
      INTEGER           NOSOEL(64), NBSTEF(5:8), NBAREF(5:8)
      include"./incl/nusc1c6.inc"
      DATA              NBSTEF / 4, 6,  8, 64 /
      DATA              NBAREF / 6, 9, 12, 12 /
C
C     LE VOLUME INITIAL
C     =================
C     LE TABLEAU LEXIQUE DE CE VOLUME
      CALL LXLXOU( NTVOLU, NMVOLU, NTLXVL, MN )
      IF( NTLXVL .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'VOLUME INCONNU :' // NMVOLU
         CALL LEREUR
         RETURN
      ENDIF
C
C     LE TABLEAU 'NSEF' DE CE VOLUME
      CALL LXTSOU( NTLXVL, 'NSEF', NTCUVL, MNCUVL )
      IF( NTCUVL .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'VOLUME SANS NSEF :' // NMVOLU
         CALL LEREUR
         RETURN
      ENDIF
C
C     LE TABLEAU 'XYZSOMMET' DE CE VOLUME
      CALL LXTSOU( NTLXVL, 'XYZSOMMET', NTSOVL, MNSOVL )
      IF( NTSOVL .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'VOLUME SANS XYZSOMMET :' // NMVOLU
         CALL LEREUR
         RETURN
      ENDIF
C     DIMENSION DE L'ESPACE ou NOMBRE DE COORDONNEES 3 en general
      NBCOOR = MCN( MNSOVL + WBCOOR )
C
C     LE NOMBRE DE SOMMETS DU MAILLAGE DU VOLUME
      NBSOVL = MCN( MNSOVL + WNBSOM )
C
C     LE NOMBRE DE TANGENTES DU MAILLAGE DU VOLUME
      NBTGS = MCN( MNSOVL + WNBTGS )
C
      MNST  = MNSOVL + WYZSOM - NBCOOR
      MNTG  = MNST   + NBCOOR * NBSOVL
C
C     LES PARAMETRES DES NO SOMMET DU MAILLAGE DU VOLUME
      CALL NSEFPA( MCN(MNCUVL),
     %             NUTYMA, NBSOEL, NBSOEF, NBTGEF,
     %             LDAPEF, LDNGEF, LDTGEF, NBCUVO,
     %             NX    , NY    , NZ    ,
     %             IERR   )
      IF( IERR .NE. 0 ) RETURN
C
      DO 10 NCOGEL=5,9
         IF( NCOGEL .EQ. 8 ) GOTO 10
C
C        LE NUMERO DES 2 SOMMETS DES ARETES DE L'EF DE TYPE NCOGEL
         CALL SOARFA( NCOGEL, NUSOAR(1,1,NCOGEL), NOSTFA )
C
         IF( NBTGS .GT. 0 ) THEN
C           LE NUMERO DES 2 TGS DES ARETES DE L'EF DE TYPE NCOGEL
            DO 5 K=1,NBAREF(NCOGEL)
               CALL TGAREF( NCOGEL, K, NUTGAR(1,K,NCOGEL) )
 5          CONTINUE
         ENDIF
C
 10   CONTINUE
C
C     LA BOUCLE SUR LES CUBES DU MAILLAGE DU VOLUME
C     ---------------------------------------------
      MNS1 = 0
      MNS2 = 0
      DO 100 N=1,NBCUVO
C
C        LE NUMERO DES SOMMETS DU CUBE N
         CALL NSEFNS( N     , NUTYMA, NBSOEF, NBTGEF,
     %                LDAPEF, LDNGEF, LDTGEF,
     %                MNCUVL, NX, NY, NZ,
     %                NCOGEL, NUGEEF, NUEFTG, NOSOEL, IERR )
C
         IF( NUEFTG .GT. 0 .AND. NBTGS .GT. 0 ) THEN
C
            DO 40 K=1,NBAREF(NCOGEL)
C
C              EF A TG :
C              L'ADRESSE DES COORDONNEES DES 2 SOMMETS DE L'ARETE K
               MNS1 = MNST + NBCOOR * NOSOEL( NUSOAR(1,K,NCOGEL) )
               MNS2 = MNST + NBCOOR * NOSOEL( NUSOAR(2,K,NCOGEL) )
C
C              LE NUMERO DES 2 TGS DE L'ARETE
               MN   = MNCUVL + LDTGEF + NBTGEF * ( NUEFTG - 1 ) - 1
               NTG1 = MCN( MN + NUTGAR(1,K,NCOGEL) )
               NTG2 = MCN( MN + NUTGAR(2,K,NCOGEL) )
               IF( NTG1 .EQ. 0 .AND. NTG2 .EQ. 0 ) GOTO 30
C
C              LES 3 COMPOSANTES DE LA TANGENTE 1
               IF( NTG1 .NE. 0 ) THEN
                  MN = MNTG + NBCOOR * ABS(NTG1)
                  IF( NTG1 .GT. 0 ) THEN
                     XYZTG(1,1) = RMCN( MN     )
                     XYZTG(2,1) = RMCN( MN + 1 )
                     XYZTG(3,1) = RMCN( MN + 2 )
                  ELSE
C                    TANGENTE A INVERSER
                     XYZTG(1,1) =-RMCN( MN     )
                     XYZTG(2,1) =-RMCN( MN + 1 )
                     XYZTG(3,1) =-RMCN( MN + 2 )
                  ENDIF
               ELSE
                  XYZTG(1,1) = RMCN(MNS2  ) - RMCN(MNS1  )
                  XYZTG(2,1) = RMCN(MNS2+1) - RMCN(MNS1+1)
                  XYZTG(3,1) = RMCN(MNS2+2) - RMCN(MNS1+2)
               ENDIF
C
C              LES 3 COMPOSANTES DE LA TANGENTE 2
               IF( NTG2 .NE. 0 ) THEN
                  MN = MNTG + NBCOOR * ABS(NTG2)
                  IF( NTG2 .GT. 0 ) THEN
                     XYZTG(1,2) = RMCN( MN     )
                     XYZTG(2,2) = RMCN( MN + 1 )
                     XYZTG(3,2) = RMCN( MN + 2 )
                  ELSE
C                    TANGENTE A INVERSER
                     XYZTG(1,2) =-RMCN( MN     )
                     XYZTG(2,2) =-RMCN( MN + 1 )
                     XYZTG(3,2) =-RMCN( MN + 2 )
                  ENDIF
               ELSE
                  XYZTG(1,2) = RMCN(MNS1  ) - RMCN(MNS2  )
                  XYZTG(2,2) = RMCN(MNS1+1) - RMCN(MNS2+1)
                  XYZTG(3,2) = RMCN(MNS1+2) - RMCN(MNS2+2)
               ENDIF
C
C              LE TRACE DE L'ARETE P3 ** 3
               CALL TRAR3D( NCOUAL, PREDUA,
     %                      RMCN(MNS1), RMCN(MNS2),
     %                      XYZTG(1,1), XYZTG(1,2) )
               GOTO 40
C
C              EF SANS TG => TRACE DES ARETES P1
 30            CALL TRAIT3D( NCOUAF, RMCN(MNS1), RMCN(MNS2) )
C
 40         CONTINUE
C
         ELSE
C
C           EF SANS TG => TRACE DES ARETES P1
            IF( NBCOOR .NE. 6 ) THEN
               DO 50 K = 1, NBAREF(NCOGEL)
                  MNS1 = MNST + NBCOOR * NOSOEL( NUSOAR(1,K,NCOGEL) )
                  MNS2 = MNST + NBCOOR * NOSOEL( NUSOAR(2,K,NCOGEL) )
                  CALL TRAIT3D( NCOUAF, RMCN(MNS1), RMCN(MNS2) )
 50            CONTINUE
            ELSE
C              TRACE DES 192 ARETES=1-CUBES DU 6-CUBE
               DO 60 K = 1, NB1C6C
                  MNS1 = MNST + NBCOOR * NOSOEL( NUS1C6C(1,K) )
                  MNS2 = MNST + NBCOOR * NOSOEL( NUS1C6C(2,K) )
                  CALL TRAIT3D( NCOUAF, RMCN(MNS1), RMCN(MNS2) )
 60            CONTINUE
            ENDIF
         ENDIF
C
C        TRACE EVENTUEL DU NO DES SOMMETS
         IF( IAVNSO .NE. 0 ) THEN
            DO 70 K=1,NBSTEF(NCOGEL)
               NS = NOSOEL(K)
               WRITE( NMSOMM , '(I8)' ) NS
               CALL SANSBL( NMSOMM, L )
               MNS1 = MNST + NBCOOR * NS
               CALL TEXTE3D( NCONSO, RMCN(MNS1), NMSOMM(1:L) )
 70         CONTINUE
         ENDIF
 100  CONTINUE
C
C     LE TRACE DE LA POIGNEE DU VOLUME AU MILIEU DE LA DERNIERE ARETE TRACEE
      XYZTG(1,1) = ( RMCN(MNS1  ) + RMCN(MNS2  ) ) * 0.5
      XYZTG(2,1) = ( RMCN(MNS1+1) + RMCN(MNS2+1) ) * 0.5
      XYZTG(3,1) = ( RMCN(MNS1+2) + RMCN(MNS2+2) ) * 0.5
      CALL NUOBNM( 'VOLUME', NMVOLU, K )
      CALL ITEMV3(  XYZTG,   NMVOLU, K )
C
      RETURN
      END
