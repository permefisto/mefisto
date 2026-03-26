      SUBROUTINE SUC0G1( NTLXSF, NBSPSV, NUSPSV,
     %                   NTNSEF, MNNSEF, NTXYZS, MNXYZS, IERR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    PROJECTION DES TANGENTES D'UN MEME SOMMET SUR LE PLAN
C -----    A DISTANCE MINIMALE DES EXTREMITES DES TANGENTES
C          => SURFACE MAILLEE G1-CONTINUE SI TGS EN TOUT SOMMET
C
C REMARQUE: SEULES LES COMPOSANTES DES TANGENTES DU TMS XYZSOMMET
C           SONT EVENTUELLEMENT MODIFIEES
C
C ENTREES:
C --------
C NTLXSF : NUMERO TMS DU LEXIQUE DE LA SURFACE DE MAILLAGE A TRAITER
C NBSPSV : NOMBRE D'ENTIERS PLUS UN DU TABLEAU NUSPSV
C          =0 SI TOUS LES SOMMETS SONT A CONSIDERER POUR LA PROJECTION
C          >0 SI SEULS LES SOMMETS TELS QUE NUSPSV(I)=-1 SONT A CONSIDERER
C NUSPSV : NUSPSV(I)=-1 SI LE SOMMET I DOIT AVOIR SES TANGENTES A PROJETER
C          SUR LE PLAN A DISTANCE MINIMALE DES EXTREMITES DES TANGENTES
C          CE TABLEAU N'EST PAS UTILISE SI NBSPSV EST =0
C
C SORTIES:
C --------
C NTNSEF : NUMERO TMS 'NSEF'      DU MAILLAGE DE LA SURFACE
C MNNSEF : ADRESSE MCN DU TMS 'NSEF'
C NTXYZS : NUMERO TMS 'XYZSOMMET' DU MAILLAGE DE LA SURFACE
C MNXYZS : ADRESSE MCN DU TMS 'XYZSOMMET'
C IERR   : =0 SI PAS D'ERREUR RENCONTREE
C          >0 SINON
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC PARIS        JUIN 1998
C ...................................................................012
      IMPLICIT INTEGER (W)
      include"./incl/gsmenu.inc"
      include"./incl/a___nsef.inc"
      include"./incl/a___xyzsommet.inc"
C
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
C
      INTEGER           NUSPSV(0:NBSPSV)
      INTEGER           NOSOEL(12), NUTG(2)
      DOUBLE PRECISION  PLAN(4)
C
C     LE TABLEAU 'NSEF' DE CETTE SURFACE
      CALL LXTSOU( NTLXSF, 'NSEF', NTNSEF, MNNSEF )
      IF( NTNSEF .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'SUC0G1: SURFACE SANS NSEF'
         CALL LEREUR
         IERR = 2
         RETURN
      ENDIF
C
C     LES PARAMETRES DU NO DES SOMMETS DU MAILLAGE DE LA SURFACE
      CALL NSEFPA( MCN(MNNSEF),
     %             NUTYMA, NBSOEL, NBSOEF, NBTGEF,
     %             LDAPEF, LDNGEF, LDTGEF, NBEFOB,
     %             NX,     NY,     NZ,
     %             IERR   )
      IF( IERR .NE. 0 ) RETURN
C
C     LE NOMBRE D'EF A TG DE LA SURFACE
      IF( MCN( MNNSEF + WBEFTG ) .LE. 0 ) RETURN
C
C     LE TABLEAU 'XYZSOMMET' DE CETTE SURFACE
      CALL LXTSOU( NTLXSF, 'XYZSOMMET', NTXYZS, MNXYZS )
      IF( NTXYZS .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'SUC0G1: SURFACE SANS XYZSOMMET'
         CALL LEREUR
         IERR = 3
         RETURN
      ENDIF
C
C     LE NOMBRE DE TANGENTES DU MAILLAGE DE LA SURFACE
      IF( MCN( MNXYZS + WNBTGS ) .LE. 0 ) RETURN
C
C     LE NOMBRE DE SOMMETS DU MAILLAGE DE LA SURFACE
      NBSOM = MCN( MNXYZS + WNBSOM )
C
C     ADRESSE MCN -3 DE LA 1-ERE COMPOSANTE DE LA PREMIERE TG
      MNXYTG = MNXYZS + WYZSOM + 3 * NBSOM - 3
C
C     CALCUL DU NOMBRE D'EF DE CHAQUE SOMMET DU MAILLAGE
C     ==================================================
      MN3CTG = 0
      MNNUTG = 0
      MNNBES = 0
      MONBES = 1 + NBSOM
      CALL TNMCDC( 'ENTIER', MONBES, MNNBES )
      CALL AZEROI( MONBES, MCN(MNNBES) )
C
      DO 20 NUEF=1,NBEFOB
C
C        LE NUMERO DES SOMMETS DE L'EF NUEF
         CALL NSEFNS( NUEF,   NUTYMA, NBSOEF, NBTGEF,
     %                LDAPEF, LDNGEF, LDTGEF,
     %                MNNSEF, NX, NY, NZ,
     %                NCOGEL, NUGEEF, NUEFTG, NOSOEL, IERR )
         IF( NUEFTG .LE. 0 ) GOTO 20
         DO 10 I=1,NCOGEL
C           UN EF DE PLUS CONTIENT LE SOMMET NS
            NS = NOSOEL(I)
            MCN(MNNBES+NS) = MCN(MNNBES+NS) + 1
 10      CONTINUE
 20   CONTINUE
C
C     POINTEUR SUR LA DERNIERE TANGENTE DE CHAQUE SOMMET
C     3 TROUS PAR SOMMET SONT PERMIS
      DO 30 I=1,NBSOM
         MCN(MNNBES+I) = MCN(MNNBES+I-1) + 2*MCN(MNNBES+I) + 3
 30   CONTINUE
C
C     CALCUL DU NUMERO DES TANGENTES DE CHAQUE SOMMET DU MAILLAGE
C     ===========================================================
      MXTG1S = 0
      MONUTG = 1+MCN(MNNBES+NBSOM)
      CALL TNMCDC( 'ENTIER', MONUTG, MNNUTG )
      CALL AZEROI( MONUTG, MCN(MNNUTG) )
C
      DO 50 NUEF=1,NBEFOB
C
C        LE NUMERO DES SOMMETS DE L'EF NUEF
         CALL NSEFNS( NUEF,   NUTYMA, NBSOEF, NBTGEF,
     %                LDAPEF, LDNGEF, LDTGEF,
     %                MNNSEF, NX, NY, NZ,
     %                NCOGEL, NUGEEF, NUEFTG, NOSOEL, IERR )
         IF( NUEFTG .LE. 0 ) GOTO 50
C        EF A TG
         DO 40 I=1,NCOGEL
C           NUMERO DU SOMMET I
            NS = NOSOEL(I)
C           ADRESSE DES TG DE L'EF A TG
            MN = MNNSEF + LDTGEF + NBTGEF * (NUEFTG-1) - 1
C           LE NUMERO DES 2 TANGENTES DU SOMMET I
            NUTG2 = 2 * I
            NUTG(1) = MCN( MN + NUTG2 - 1)
            NUTG(2) = MCN( MN + NUTG2 )
            J0 = MCN(MNNBES+NS-1)+1
            J1 = MCN(MNNBES+NS)
            DO 36 L=1,2
C              TRAITEMENT DE LA TG L DU SOMMET I DE L'EF NUEF
               IF( NUTG(L) .EQ. 0 ) GOTO 36
               DO 34 J=J0,J1
C                 RECHERCHE DE LA DERNIERE TG RECENSEE
                  IF( MCN(MNNUTG+J) .EQ. 0 ) THEN
C                    C'EST J-1
                     DO 32 K=J0,J-1
                        IF( MCN(MNNUTG+K) .EQ. NUTG(L) ) GOTO 36
 32                  CONTINUE
C                    TG NON RETROUVEE DONC AJOUTEE
                     MCN(MNNUTG+J) = NUTG(L)
C                    NOMBRE MAXIMAL DE TGS-1 EN UN SOMMET
                     MXTG1S = MAX( MXTG1S, J-J0 )
                     GOTO 36
                  ENDIF
 34            CONTINUE
 36         CONTINUE
 40      CONTINUE
 50   CONTINUE
      MXTG1S = MXTG1S + 1
C
C     PROJECTION DES TGS EN CHACUN DES SOMMETS DU MAILLAGE
C     SOIT, TOUS LES SOMMETS
C     ====================================================
      MN3CTG = 0
      CALL TNMCDC( 'REEL', 3*MXTG1S, MN3CTG )
C
      J1 = MCN(MNNBES)
      DO 100 NS=1,NBSOM
C
C        POINTEUR 1ERE TG
         J0 = J1 + 1
C
C        POINTEUR DERNIERE TG
         J1 = MCN(MNNBES+NS)
C
         IF( NBSPSV .GT. 0 ) THEN
C           LE SOMMET EST IL A PROJETER?
            IF( NUSPSV(NS) .NE. -1 ) THEN
C              NON
               GOTO 100
            ENDIF
         ENDIF
C
         DO 60 J=J0,J1
C           RECHERCHE DE LA DERNIERE TG RECENSEE J-1
            IF( MCN(MNNUTG+J) .EQ. 0 ) GOTO 70
 60      CONTINUE
C        PAS DE TG EN CE SOMMET
         GOTO 100
C
C        LES TGS DU SOMMET NS VONT DE MCN(MNNUTG+J0) A MCN(MNNUTG+J-1)
C        CONSTRUCTION DES 3 COMPOSANTES DES TGS DU SOMMET NS
 70      MN3 = MN3CTG
         MNS = MNXYZS + WYZSOM + 3 * NS - 3
         DO 80 K=J0,J-1
C           LE NUMERO DE LA TANGENTE
            L = MCN(MNNUTG+K)
            IF( L .GT. 0 ) THEN
C              SIGNE DE LA TANGENTE
               N = 1
            ELSE
C              SIGNE DE LA TANGENTE
               L = -L
               N = -1
            ENDIF
C           LES 3 COMPOSANTES DE LA TANGENTE
            MN = MNXYTG + 3 * L
            RMCN(MN3  ) = N * RMCN(MN  )
            RMCN(MN3+1) = N * RMCN(MN+1)
            RMCN(MN3+2) = N * RMCN(MN+2)
            MN3 = MN3 + 3
 80      CONTINUE
C
C        CALCUL DU PLAN A DISTANCE MINIMALE DES EXTREMITES DES TGS AU SOMMET NS
         CALL PLDIM2( RMCN(MNS), J-J0, RMCN(MN3CTG), PLAN, IERR )
         IF( IERR .NE. 0 ) GOTO 100
C
C        PROJECTION DES TGS SUR LE PLAN
         MN3 = MN3CTG
         DO 90 K=J0,J-1
C
C           LE NUMERO DE LA TANGENTE ET SON SIGNE
            L  = MCN(MNNUTG+K)
            IF( L .GT. 0 ) THEN
C              SIGNE DE LA TANGENTE
               N = 1
            ELSE
C              SIGNE DE LA TANGENTE
               L = -L
               N = -1
            ENDIF
C
C           LE POINT EXTREMITE DE LA TG EST PROJETE SUR LE PLAN
            RMCN(MN3  ) = RMCN(MNS  ) + RMCN(MN3  )
            RMCN(MN3+1) = RMCN(MNS+1) + RMCN(MN3+1)
            RMCN(MN3+2) = RMCN(MNS+2) + RMCN(MN3+2)
            CALL PRPTP2( RMCN(MN3), PLAN, RMCN(MN3), IERR )
C
C           LES 3 COMPOSANTES DE LA NOUVELLE TANGENTE
            MN = MNXYTG + 3 * L
            RMCN(MN  ) = N * ( RMCN(MN3  ) - RMCN(MNS  ) )
            RMCN(MN+1) = N * ( RMCN(MN3+1) - RMCN(MNS+1) )
            RMCN(MN+2) = N * ( RMCN(MN3+2) - RMCN(MNS+2) )
            MN3 = MN3 + 3
 90      CONTINUE
C
 100  CONTINUE
C
C     DESTRUCTION DES TABLEAUX AUXILIAIRES
      CALL TNMCDS( 'REEL', 3*MXTG1S, MN3CTG )
      CALL TNMCDS( 'ENTIER', MONUTG, MNNUTG )
      CALL TNMCDS( 'ENTIER', MONBES, MNNBES )
      IERR = 0
C
C     REDUCTION DES TGS
      CALL MOINTG( NTNSEF, MNNSEF, NTXYZS, MNXYZS, IERR )
      RETURN
      END
