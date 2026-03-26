      SUBROUTINE MOINTG( NTNSEF, MNNSEF, NTXYZS, MNXYZS, IERR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :IDENTIFIER LES TANGENTES EGALES EN CHACUN DES SOMMETS DU MAILLAGE
C -----SUPPRIMER LES TANGENTES DOUBLES EN UN SOMMET
C
C ENTREES:
C --------
C NUTYOB : NUMERO DU TYPE DU PLSV (1:POINT, 2:LIGNE, 3:SURFACE, 4:VOLUME)
C NTNSEF : NO TMS 'NSEF'      DU MAILLAGE DU PLSV
C MNNSEF : ADRESSE MCN DU TMS 'NSEF'
C NTXYZS : NO TMS 'XYZSOMMET' DU MAILLAGE DU PLSV
C MNXYZS : ADRESSE MCN DU TMS 'XYZSOMMET'
C
C SORTIES:
C --------
C IERR   : =0 SI PAS D'ERREUR RENCONTREE
C          >0 SINON
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC PARIS    DECEMBRE 1996
C2345X7..............................................................012
      IMPLICIT INTEGER (W)
      include"./incl/gsmenu.inc"
      include"./incl/a___nsef.inc"
      include"./incl/a___xyzsommet.inc"
C
      COMMON / EPSSSS / EPZERO,EPSXYZ
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
C
      INTEGER           NOSOEL(64)
      REAL              COMP(3)
C
      MNNEWT = 0
      MNTGST = 0
C     LE NOMBRE DE COORDONNEES DES SOMMETS ET TGS DU MAILLAGE DU PLSV
      NBCOOR = 3
C
C     LE TABLEAU 'NSEF' DU PLSV
      IF( NTNSEF .LE. 0 .OR. MNNSEF .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'MOINTG: PLSV SANS NSEF'
         CALL LEREUR
         IERR   = 2
         NTXYZS = 0
         MNXYZS = 0
         RETURN
      ENDIF
C
C     NUTYOB : NUMERO DU TYPE DU PLSV
C    (1:POINT, 2:LIGNE, 3:SURFACE, 4:VOLUME)
      NUTYOB = MCN( MNNSEF + WUTYOB )
      IF( NUTYOB .LE. 1 .OR. NUTYOB .GT. 4 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'MOINTG: TYPE DE PLSV NON EGAL A 2 ou 3 ou 4'
         CALL LEREUR
         IERR = 1
         RETURN
      ELSE
C        LE NOMBRE MAXIMAL DE TANGENTES STOCKABLES PAR SOMMET
         MXTGST = NBCOOR * 2 ** (NUTYOB-1)
      ENDIF
C
C     LES PARAMETRES DES NO SOMMET DU MAILLAGE DU PLSV
      CALL NSEFPA( MCN(MNNSEF),
     %             NUTYMA, NBSOEL, NBSOEF, NBTGEF,
     %             LDAPEF, LDNGEF, LDTGEF, NBEFOB,
     %             NX,     NY,     NZ,
     %             IERR   )
      IF( IERR .NE. 0 ) RETURN
C
C     LE NOMBRE D'EF A TG DU PLSV
      NBEFTG = MCN( MNNSEF + WBEFTG )
      IF( NBEFTG .LE. 0 ) RETURN
C
C     LE TABLEAU 'XYZSOMMET' DE CETTE SURFACE
      IF( NTXYZS .LE. 0 .OR. MNXYZS .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'MOINTG: PLSV SANS XYZSOMMET'
         CALL LEREUR
         IERR = 3
         RETURN
      ENDIF
C
C     LE NOMBRE DE SOMMETS DU MAILLAGE DU PLSV
      NBSOM = MCN( MNXYZS + WNBSOM )
C
C     LE NOMBRE DE TANGENTES DU MAILLAGE DU PLSV
      NBTGS = MCN( MNXYZS + WNBTGS )
      IF( NBTGS .LE. 0 ) RETURN
C
C     ADRESSAGES - 4 DU DEBUT DES XYZ ET DES TANGENTES
      MNXYZ  = MNXYZS + WYZSOM -1 -NBCOOR
      MNXYTG = MNXYZ  + NBCOOR * NBSOM
C
C     DECLARATION DU TABLEAU DES NUMEROS DES TANGENTES AVANT ET APRES
      CALL TNMCDC( 'ENTIER', 1+NBTGS, MNNEWT )
      DO 5 I=0,NBTGS
         MCN(MNNEWT+I) = I
 5    CONTINUE
C
C     DECLARATION DU NUMERO DES TANGENTES EN CHACUN DES SOMMETS
      CALL TNMCDC( 'ENTIER', MXTGST*NBSOM, MNTGST )
      CALL AZEROI( MXTGST*NBSOM, MCN( MNTGST ) )
C
C     BOUCLE SUR LES EF DU MAILLAGE
      DO 100 NEF = 1, NBEFOB
C
C        LE NUMERO DES SOMMETS DE L'EF NEF DU PLSV
         CALL NSEFNS( NEF,    NUTYMA, NBSOEF, NBTGEF,
     %                LDAPEF, LDNGEF, LDTGEF,
     %                MNNSEF, NX, NY, NZ,
     %                NCOGEL, NUGEEF, NUEFTG, NOSOEL, IERR )
         IF( IERR   .NE. 0 ) GOTO 9990
         IF( NUEFTG .LE. 0 ) GOTO 100
C
C        EF A TG . LE NOMBRE DE TG PAR SOMMET
         NBTGST = NBTGEF / NBSOEF
C
C        LE NUMERO QUI PRECEDE LA PREMIERE TG
         NTG = 0
C
         DO 90 NS = 1, NBSOEF
C
C           LE NUMERO DU SOMMET NS DE L'EF DANS LE MAILLAGE
            NUST = NOSOEL(NS)
            IF( NUST .LE. 0 ) GOTO 90
C
C           PARCOURS DES NBTGST TANGENTES DU SOMMET NUST
            DO 80 NT = 1, NBTGST
C
               NTG  = NTG + 1
C              LE NUMERO DE 1 A NBTGS DE LA TANGENTE
               NUTG = ABS( NOSOEL( NBSOEF+NTG ) )
               IF( NUTG .EQ. 0 ) GOTO 80
C
C              LA TANGENTE EST ELLE L'ARETE?
               MNXT = MNXYTG + NBCOOR * NUTG
C              RECHERCHE DU SOMMET FINAL DE L'ARETE POUR CETTE TANGENTE
C              LE NUMERO DEPEND DU TYPE DE L'EF
               NS1   = NSFITG( NS, NT, NCOGEL )
               NUST1 = NOSOEL( NS1 )
               MNX1  = MNXYZ + NBCOOR * NUST1
               MNX   = MNXYZ + NBCOOR * NUST
C              LES NBCOOR COMPOSANTES DU VECTEUR ARETE
               COMP(1) = RMCN(MNX1+1) - RMCN(MNX+1)
               COMP(2) = RMCN(MNX1+2) - RMCN(MNX+2)
               COMP(3) = RMCN(MNX1+3) - RMCN(MNX+3)
               CALL XYZIDE( RMCN(MNXT+1), COMP, IDENTQ )
               IF( IDENTQ .EQ. 1 ) THEN
C                 LA TANGENTE EST L'ARETE
                  MCN(MNNEWT+NUTG) = 0
                  GOTO 80
               ENDIF
C
C              LE NUMERO ACTUALISE DE CETTE TANGENTE
               NEWTG = MCN(MNNEWT+NUTG)
               IF( NEWTG .EQ. 0 ) THEN
C
C                 TANGENTE IDENTIFIEE A UNE ARETE => MISE A ZERO ENSUITE
                  GOTO 80
C
               ELSE IF( ABS(NEWTG) .NE. NUTG ) THEN
C
C                 TANGENTE DEJA IDENTIFIEE A UNE AUTRE NEWTG (>0 ou <0)
C                 LE SIGNE DE LA TANGENTE ACTUELLE DANS L'EF COURANT
                  IF( NOSOEL( NBSOEF+NTG ) .LT. 0 ) THEN
                     LSIGNE = -1
                  ELSE
                     LSIGNE =  1
                  ENDIF
C                 MISE A JOUR DANS NSEF
                  MN = MNNSEF + LDTGEF + NBTGEF * (NUEFTG-1) - 1
                  MCN(MN+NTG) = NEWTG * LSIGNE
C
               ELSE
C
C                 TANGENTE NON IDENTIFIEE A UNE AUTRE
C                 RECHERCHE D'UNE TANGENTE IDENTIQUE P
                  MNT = MNTGST + NUST * MXTGST - MXTGST - 1
                  DO 30 I=1,MXTGST
C                    LE NUMERO DE LA TANGENTE RECENSEE AU SOMMET NUST
                     NUTGR = MCN( MNT + I )
                     IF( NUTGR .EQ. 0 ) GOTO 40
C
C                    TENTATIVE D'IDENTIFICATION DES 2 TANGENTES NUTGR ET NUTG
                     CALL TGSIDE( RMCN(MNXYTG+NBCOOR*NUTGR+1),
     %                            RMCN(MNXT+1),
     %                            LSIGNE )
                     IF( LSIGNE .EQ. 0 ) THEN
C
C                       TANGENTES NON IDENTIFIEES
                        GOTO 30
C
                     ELSE
C
C                       LES 2 TANGENTES SONT IDENTIFIEES AU SIGNE PRES
                        NEWTG = LSIGNE * NUTGR
                        MCN(MNNEWT+NUTG) = NEWTG
C                       LE SIGNE DE LA TANGENTE ACTUELLE DANS L'EF COURANT
                        IF( NOSOEL( NBSOEF+NTG ) .LT. 0 ) THEN
                           LSIGNE = -1
                        ELSE
                           LSIGNE =  1
                        ENDIF
C                       MISE A JOUR DANS NSEF
                        MN = MNNSEF + LDTGEF + NBTGEF * (NUEFTG-1) - 1
                        MCN(MN+NTG) = NEWTG * LSIGNE
                        GOTO 80
                     ENDIF
 30               CONTINUE
C
C                 PAS D'IDENTIFICATION : CETTE TANGENTE EST RECENSEE POUR LE SOM
 40               IF( I .LT. MXTGST ) THEN
C                    IL RESTE DE LA PLACE POUR STOCKER CETTE TANGENTE
                     MCN( MNT + I ) = NUTG
                  ENDIF
               ENDIF
 80         CONTINUE
 90      CONTINUE
 100  CONTINUE
C
C     LE NOUVEAU NUMERO DES TANGENTES RECENSEES
      DO 130 NT=1, NBTGS
C        LE NUMERO DE LA TANGENTE EST RECENSE
         NUTG = MCN(MNNEWT+NT)
         IF( NT .EQ. NUTG ) THEN
C           TANGENTE NON IDENTIFIEE
            MCN(MNNEWT+NT) = NUTG
         ELSE
C           TANGENTE IDENTIFIEE
            MCN(MNNEWT+NT) = 0
         ENDIF
 130  CONTINUE
C
C     LA NOUVELLE NUMEROTATION DES TANGENTES RESTANTES
      NBTGS1 = 0
      DO 140 NT=1,NBTGS
C        LE NUMERO ANCIEN DE LA TANGENTE
         NUTG = MCN(MNNEWT+NT)
         IF( NUTG .NE. 0 ) THEN
C           LE NUMERO NOUVEAU DE LA TANGENTE
            NBTGS1 = NBTGS1 + 1
            MCN(MNNEWT+NT) = NBTGS1
         ENDIF
 140  CONTINUE
C
C     LES NBCOOR COMPOSANTES DES TANGENTES DU MAILLAGE DANS LE NOUVEL
C     ORDRE SONT RANGEES DANS LE TABLEAU COMP
      IF( NBTGS1 .GT. 0 ) THEN
         CALL TNMCDC( 'REEL', NBCOOR*NBTGS1, MNCOMP )
         MN0 = MNXYZS + WYZSOM - NBCOOR + NBCOOR * NBSOM
         MN1 = MNCOMP
         DO 150 NT=1,NBTGS
            NUTG = MCN(MNNEWT+NT)
            IF( NUTG .NE. 0 ) THEN
               MN = MN0 + NBCOOR * NT
               RMCN(MN1  ) = RMCN(MN  )
               RMCN(MN1+1) = RMCN(MN+1)
               RMCN(MN1+2) = RMCN(MN+2)
               MN1 = MN1 + NBCOOR
            ENDIF
 150     CONTINUE
C
C        RECOPIE DANS LE TMS 'XYZSOMMET'
         CALL TRTATA( MCN(MNCOMP), MCN(MN0+NBCOOR), NBCOOR*NBTGS1 )
C        DESTRUCTION DU TABLEAU AUXILIAIRE
         CALL TNMCDS( 'REEL', NBCOOR*NBTGS1, MNCOMP )
      ENDIF
C
C     LE NOMBRE DE TANGENTES DE CE MAILLAGE
      MCN( MNXYZS + WNBTGS ) = NBTGS1
C
C     AJOUT DE LA DATE
      CALL ECDATE( MCN(MNXYZS) )
C     AJOUT DU NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNXYZS + MOTVAR(6) ) = NONMTD( '~>>>XYZSOMMET' )
C
C     REDUCTION DE LA TAILLE DU TMS
      IF( NBTGS1 .LT. NBTGS ) THEN
         CALL TAMSRA( NTXYZS, WYZSOM + NBCOOR * (NBSOM+NBTGS1) )
      ENDIF
C
C     MISE A JOUR DU NUMERO DES TANGENTES DES EF A TG
      MN0 = MNNSEF + LDTGEF - 1
      DO 190 NS = 1, NBEFTG
C        +- LE NUMERO DES NBTGEF TG
         DO 180 NT=1,NBTGEF
C           +- L'ANCIEN NUMERO DE LA TG
            NUTG = MCN(MN0+NT)
            IF( NUTG .LT. 0 ) THEN
               LESIGN = -1
               NUTG   = -NUTG
            ELSE
               LESIGN = 1
            ENDIF
C           +- LE NOUVEAU NUMERO DE LA TG
            MCN(MN0+NT) = LESIGN * MCN(MNNEWT+NUTG)
 180     CONTINUE
         MN0 = MN0 + NBTGEF
 190  CONTINUE
C
C     AJOUT DE LA DATE
      CALL ECDATE( MCN(MNNSEF) )
C     AJOUT DU NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNNSEF + MOTVAR(6) ) = NONMTD( '~>>>NSEF' )
C
C     ERREUR DETECTEE OU NON
 9990 IF( MNNEWT .GT. 0 ) CALL TNMCDS( 'ENTIER', 1+NBTGS, MNNEWT )
      IF( MNTGST .GT. 0 ) CALL TNMCDS( 'ENTIER', MXTGST*NBSOM, MNTGST )
C
C     SUPPRESSION DES FAUX EF A TANGENTES
C     ===================================
      IF( IERR .EQ. 0 ) CALL FAEFTG( NTNSEF, MNNSEF, IERR )
      END
