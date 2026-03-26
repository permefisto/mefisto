      SUBROUTINE TRPLSV( NUTYOB, NUPLSV,
     %                   NTXYZS, MNXYZS, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    FAIRE SUBIR LA TRANSFORMATION DU PLSV NUPLSV DE TYPE NUTYOB
C -----    CETTE TRANSFORMATION EST
C          SOIT UNE TRANSFORMATION DEFINIE PAR 3 FONCTIONS
C          SOIT UNE ISOMETRIE
C          SOIT UNE PROJECTION SELON X OU Y OU Z OU UN POINT OU
C               SELON LA NORMALE AU PLAN DEFINI PAR 3 POINTS
C               (PAS DE TANGENTES EN SORTIE DANS CE CAS)
C          SOIT UNE TRANSFORMATION DEFINIE PAR 12 FONCTIONS F(3) ET DF(3,3)
C
C ENTREES:
C --------
C NUTYOB : NUMERO DU TYPE DU PLSV (1:POINT 2:LIGNE, ... )
C NUPLSV : NUMERO DU PLSV DANS LE LEXIQUE DES POINT OU LIGNE OU SURFACE OU VOLUM
C          LE TABLEAU DE DEFINITION D'UNE TRANSFORMATION EST DECRIT DANS
C          LE FICHIER  ~/td/d/a_transfo__definition
C
C SORTIES:
C --------
C NTXYZS : NUMERO      DU TABLEAU 'XYZSOMMET' DU PLSV  XYZ DES SOMMETS ET TANGEN
C MNXYZS : ADRESSE MCN DU TABLEAU 'XYZSOMMET' DU PLSV  XYZ DES SOMMETS ET TANGEN
C NTNSEF : NUMERO      DU TMS 'NSEF' DES NUMEROS DES SOMMETS ET TGS DES EF
C MNNSEF : ADRESSE MCN DU TMS 'NSEF' DES NUMEROS DES SOMMETS ET TGS DES EF
C IERR   : 0   PAS D'ERREUR
C          1   PLSV INCONNU
C          2   DEFINITION INCONNUE DU PLSV
C          3   TRANSFORMATION INCONNUE
C          3   3 FONCTIONS IDENTITE
C          4   TYPE DE PROJECTION INCONNU
C          5   FONCTION SURFACE INCONNU
C          6   CENTRE DE PROJECTION INCONNU
C          7   FONCTION X A PROBLEME
C          8   FONCTION Y A PROBLEME
C          9   FONCTION Z A PROBLEME
C         10   3 POINTS INCONNUS POUR DEFINIR LE PLAN
C              DE NORMALE DIRECTION DE LA PROJECTION
C         11   BISSECTION NON CONVERGENTE LORS DE LA PROJECTION
C         12   PROJECTION AVEC FONCTION SANS RESULTAT
C         20   TRANSFORMATION RENDANT UN EF DE VOLUME<0
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     JANVIER 1997
C23456---------------------------------------------------------------012
      PARAMETER        (MXITER=24)
      IMPLICIT          INTEGER (W)
      include"./incl/gsmenu.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/ponoel.inc"
      include"./incl/a_point__definition.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/a_transfo__definition.inc"
      COMMON / EPSSSS / EPZERO,EPSXYZ
C
      include"./incl/pp.inc"
      COMMON             MCN(MOTMCN)
      REAL              RMCN(1)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE      (DMCN(1),RMCN(1),MCN(1))
C
      REAL              P1(3),P2(3),P3(3),XYZN(3),XYZ(3),COIN(6,2),
     %                  DIRTG(3)
      EQUIVALENCE      (DIRTG,P2)
      DOUBLE PRECISION  DMATRI(3,4),DP(3),DBLVAL,D
      INTEGER           NOSOEL(64),NUFODF(3,3),
     %                  NOST(2),NOTGAR(2),NOTG(2)
      CHARACTER*1       KXYZ(3)
      DATA              KXYZ/ 'X', 'Y', 'Z' /
C
      IF( NUTYOB .LE. 0 .OR. NUTYOB .GT. 4 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'trplsv: L''OBJET N''EST PAS UN PLSV'
         CALL LEREUR
         IERR = 1
         GOTO 9999
      ENDIF
C
      NEWGEO = 1
      NTNSEF = 0
      MNNSEF = 0
      MNST1  = 0
      MNDF1  = 0
      D      = 0
C
C     NUMERO TS DU LEXIQUE DES PLSVS DE CE TYPE D'PLSV
      NTLX   = NTMN( NUTYOB )
C
C     LE NUMERO DE TMS DU LEXIQUE DU PLSV FINAL
      CALL LXNLOU( NTLX, NUPLSV, NTLXOB, MNLXOB )
      IF( NTLXOB .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'trplsv:PLSV INCONNU'
         CALL LEREUR
         IERR = 1
         GOTO 9999
      ENDIF
C
C     LE NUMERO DE TMS DU TABLEAU DEFINITION DE CET PLSV
      CALL LXTSOU( NTLXOB, 'DEFINITION', NTDFOB, MNDFOB )
      IF( NTDFOB .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'trplsv:DEFINITION INCONNUE DU ''PLSV'
         CALL LEREUR
         IERR = 2
         GOTO 9999
      ENDIF
C
C     LE NUMERO DE LA TRANSFORMATION DU PLSV
      NUTRAN = MCN( MNDFOB + WTYTRP )
      IF( NUTRAN .LE. 1 ) THEN
C        IDENTITE
         IERR = 0
         RETURN
      ENDIF
C
C     LE LEXIQUE DE LA TRANSFORMATION
      CALL LXNLOU( NTTRAN, NUTRAN, NTLXTR, MNLXTR )
      IF( NTLXTR .LE. 0 ) THEN
         NUTYTR = NUTRAN
         GOTO 90
      ENDIF
C
C     LE NUMERO DE TMS DU TABLEAU DEFINITION DE CETTE TRANSFORMATION
      CALL LXTSOU( NTLXTR, 'DEFINITION', NTDFTR, MNDFTR )
      IF( NTDFTR .LE. 0 ) THEN
         NBLGRC(NRERR) = 2
         WRITE(KERR(MXLGER)(1:4),'(I4)') NUTRAN
         KERR(1) = 'trplsv:DEFINITION INCONNUE DE LA'
         KERR(2) = 'TRANSFORMATION '  // KERR(MXLGER)(1:4)
         CALL LEREUR
         IERR = 3
         GOTO 9999
      ENDIF
C
C     LE TABLEAU 'XYZSOMMET' DE CET PLSV
      CALL LXTSOU( NTLXOB, 'XYZSOMMET', NTXYZS, MNXYZS )
      IF( NTXYZS .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'trplsv:TABLEAU ''XYZSOMMET'' INCONNU'
         CALL LEREUR
         IERR = 2
         GOTO 9999
      ENDIF
C
C     LE TABLEAU 'NSEF' DE CE PLSV
      IF( NUTYOB .GT. 1 ) THEN
         CALL LXTSOU( NTLXOB, 'NSEF', NTNSEF, MNNSEF )
         IF( NTNSEF .LE. 0 ) THEN
            NBLGRC(NRERR) = 1
            KERR(1) = 'trplsv:TABLEAU ''NSEF'' INCONNU'
            CALL LEREUR
            IERR = 2
            GOTO 9999
         ENDIF
      ENDIF
C
C     LE NOMBRE D'EF A TG
      NBEFTG = MCN( MNNSEF + WBEFTG )
C
C     LE NUMERO DU TYPE DE LA TRANSFORMATION
      NUTYTR = MCN( MNDFTR + WUTYTR )
C     NUTYTR 'NUMERO DU TYPE D''UNE TRANSFORMATION'
C     ( 1  : 'identite de R3 -> R3'
C     , 3  : 'isometrie de R3 -> R3'
C     , 4  : 'projection dans R3'
C     , 6  : 'F:R3->R3 definie par 3 fonctions (Fi)'
C     , 7  : 'F:R3->R3 12 fonctions (Fi et [dFi/dxj])' )
C
C     TRANSFORMATION DES SOMMETS
C     ==========================
C     ADRESSE DES COORDONNEES DES SOMMETS DU PLSV
      MNSO = MNXYZS + WYZSOM
C
C     LE NOMBRE DE SOMMETS
      NBSOM = MCN( MNXYZS + WNBSOM )
C
C     LE NOMBRE DE TANGENTES
      NBTGS = MCN( MNXYZS + WNBTGS )
C
      GOTO ( 9999, 90, 150, 180, 90, 100, 100, 90 ), NUTYTR
C
C     ERREUR
C     ------
 90   NBLGRC(NRERR) = 1
      WRITE(KERR(MXLGER)(1:4),'(I4)') NUTYTR
      KERR(1) = 'trplsv:TRANSFORMATION INCONNUE '//KERR(MXLGER)(1:4)
      CALL LEREUR
      IERR = 3
      GOTO 9999
C
C     3 FONCTIONS : LE NOM DES 3 FONCTIONS EST IL CORRECT ?
C     -----------------------------------------------------
 100  IF( MCN( MNDFTR+WUFONX ) .LE. 0 .AND.
     %    MCN( MNDFTR+WUFONY ) .LE. 0 .AND.
     %    MCN( MNDFTR+WUFONZ ) .LE. 0 )  THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'trplsv: 3 FONCTIONS F:R3->R3 IDENTITE'
         CALL LEREUR
         IERR = 3
         GOTO 9999
      ENDIF
C
      IF( NUTYTR .EQ. 7 ) THEN
C        DANS CE CAS LES 3 COORDONNEES INITIALES DES SOMMETS
C        DOIVENT ETRE PROTEGEES POUR PERMETTRE LE CALCUL DES TANGENTES
C        LE TABLEAU MNST1 CONTIENT LES COORDONNEES APRES TRANSFORMATION
         CALL TNMCDC( 'REEL', 3*NBSOM, MNST1 )
         CALL TRTATA( MCN(MNXYZS+WYZSOM), MCN(MNST1), 3*NBSOM )
         MNSO = MNST1
      ENDIF
C
C     TRANSFORMATION DES NBSOM SOMMETS DU MAILLAGE DU PLSV
      DO 120 I=1,NBSOM
C        LES COORDONNEES DU SOMMET I
         P1(1) = RMCN( MNSO )
         P1(2) = RMCN( MNSO + 1 )
         P1(3) = RMCN( MNSO + 2 )
C        TRANSFORMATION DE CE SOMMET PAR LES 3 FONCTIONS
         DO 110 J=1,3
C           LE NUMERO DE LA FONCTION
            NUFONC = MCN( MNDFTR + WUFONX - 1 + J )
            IF( NUFONC .LE. 0 ) THEN
C              FONCTION NON DEFINIE => IDENTITE
               P2( J ) = P1( J )
            ELSE
C              FONCTION DEFINIE
               DP(1) = P1(1)
               DP(2) = P1(2)
               DP(3) = P1(3)
               CALL FONVAL( NUFONC, 3, DP, NCODEV, DBLVAL )
               IF( NCODEV .LE. 0 ) THEN
                  IERR = 6+J
                  NBLGRC(NRERR) = 1
                  KERR(1) =  KXYZ(J)//' FONCTION INCALCULABLE'
                  CALL LEREUR
                  GOTO 9999
               ENDIF
               P2( J ) = REAL( DBLVAL )
            ENDIF
 110     CONTINUE
C        LA SAUVEGARDE DES COORDONNEES APRES TRANSFORMATION
C        TABLEAU P2 UTILE POUR LE CAS D'UN PB DE CALCUL D'UNE FONCTION
         RMCN( MNSO     ) = P2( 1 )
         RMCN( MNSO + 1 ) = P2( 2 )
         RMCN( MNSO + 2 ) = P2( 3 )
         MNSO = MNSO + 3
 120  CONTINUE
C     ICI LE TABLEAU A L'ADRESSE MNSO CONTIENT LES XYZ DES F(SOMMET)
C
C     3 FONCTIONS NE PERMETTENT PAS DE DECIDER DU NUMERO GEOMETRIQUE
C     DU PLSV CAR LA TRANSFORMATION CHANGE LA GEOMETRIE DU PLSV
      NEWGEO = 0
C
      IF( NUTYOB .EQ. 1 ) THEN
C        SI UN POINT, IL N'Y A PAS DE DIRECTION DES TGS
C        POUR L'INSTANT PAS DE TANGENTES DANS XYZSOMMET
         MCN( MNXYZS + WNBTGS ) = 0
C        LE NOMBRE DE COORDONNEES PAR SOMMET
         MCN( MNXYZS + WBCOOR ) = 3
         GOTO 8000
      ENDIF
C
C     EXISTE T IL UN GRADIENT DE F : DF A APPLIQUER SUR LES TANGENTES ?
      IF( NUTYTR .EQ. 6 ) THEN
C        NON: TANGENTES NON CALCULABLES DONC SUPPRIMEES DANS L'IMAGE
         CALL SUPTGS( NTXYZS, MNXYZS, NTNSEF, MNNSEF )
         GOTO 8000
      ENDIF
C
C     OUI : LE NUMERO DES FONCTIONS DE [DF] POUR CALCULER LES TANGENTES
      NUFODF(1,1) = MCN( MNDFTR+WUDF1U )
      NUFODF(1,2) = MCN( MNDFTR+WUDF1V )
      NUFODF(1,3) = MCN( MNDFTR+WUDF1W )
C
      NUFODF(2,1) = MCN( MNDFTR+WUDF2U )
      NUFODF(2,2) = MCN( MNDFTR+WUDF2V )
      NUFODF(2,3) = MCN( MNDFTR+WUDF2W )
C
      NUFODF(3,1) = MCN( MNDFTR+WUDF3U )
      NUFODF(3,2) = MCN( MNDFTR+WUDF3V )
      NUFODF(3,3) = MCN( MNDFTR+WUDF3W )
C
C     VERIFICATION : LES 9 FONCTIONS DOIVENT ETRE DEFINIES
      DO 123 J=1,3
         DO 122 K=1,3
            IF( NUFODF(J,K) .LE. 0 ) THEN
               NBLGRC(NRERR) = 2
               WRITE(KERR(4)(1:2),'(I2)') J
               WRITE(KERR(4)(3:4),'(I2)') K
               KERR(1) = 'LA FONCTION DU GRADIENT DE F EST NON DEFINIE'
               KERR(2) = 'POUR LA COMPOSANTE' // KERR(4)(1:2) //
     %                   ' DERIVEE PAR RAPPORT A LA VARIABLE' //
     %                   KERR(4)(3:4)
               CALL LEREUR
               IERR = 3
               GOTO 9999
            ENDIF
 122     CONTINUE
 123  CONTINUE
C
C     DECLARATION DU TABLEAU DU GRADIENT EN CHAQUE SOMMET
      CALL TNMCDC( 'REEL', 9*NBSOM, MNDF1 )
      MNDF = MNDF1
      MNSO = MNXYZS + WYZSOM
C
C     CALCUL DE LA MATRICE GRADIENT EN LES SOMMETS INITIAUX DU MAILLAGE
      DO 130 I=1,NBSOM
C
C        LES COORDONNEES INITIALES DU SOMMET I PASSENT EN DOUBLE PRECISION
         DP(1) = RMCN( MNSO     )
         DP(2) = RMCN( MNSO + 1 )
         DP(3) = RMCN( MNSO + 2 )
         MNSO  = MNSO + 3
C
C        LA MATRICE GRADIENT EN CE SOMMET IMAGE DES 9 FONCTIONS DF(J,K)
         DO 128 J=1,3
            DO 125 K=1,3
C              LA FONCTION DFj / DUk  EST APPLIQUEE
               CALL FONVAL( NUFODF(J,K), 3, DP, NCODEV, DBLVAL )
               IF( NCODEV .LE. 0 ) THEN
C                 PROBLEME LORS DU CALCUL DE LA FONCTION EN CE SOMMET
                  IERR = 6+J
                  NBLGRC(NRERR) = 2
                  WRITE(KERR(4)(1:2),'(I2)') J
                  WRITE(KERR(4)(3:4),'(I2)') K
                  KERR(1) ='LA FONCTION DU GRADIENT DE F NON CALCULABLE'
                  KERR(2) ='POUR LA COMPOSANTE' // KERR(4)(1:2) //
     %                     ' ET DERIVEE PAR RAPPORT A LA VARIABLE' //
     %                     KERR(4)(3:4)
                  CALL LEREUR
                  GOTO 9999
               ENDIF
               RMCN(MNDF) = REAL( DBLVAL )
               MNDF = MNDF + 1
 125        CONTINUE
 128     CONTINUE
 130  CONTINUE
C
C     LES PARAMETRES DU LSV
      CALL NSEFPA( MCN(MNNSEF),
     %             NUTYMA, NBSOEL, NBSOEF, NBTGEF,
     %             LDAPEF, LDNGEF, LDTGEF, NBEFOB,
     %             NX   , NY   , NZ   ,
     %             IERR   )
C
C     TOUTE DIRECTION EST UNE TANGENTE => NOUVELLE VALEUR DE NBTGEF
      IF( NUTYOB .EQ. 2 ) THEN
         NBTGEF = 2
      ELSE IF( NUTYOB .EQ. 3 .OR. NUTYOB .EQ. 4 ) THEN
         NBTGEF = 8
      ELSE
         NBTGEF = 24
      ENDIF
C
C     ALLONGEMENT DU TMS 'XYZSOMMET' POUR LES TANGENTES
      CALL TAMSAU( NTXYZS, WYZSOM + 3 * ( NBSOM + NBTGEF*NBEFOB ) )
      CALL TAMSOU( NTXYZS, MNXYZS )
      MNXYTG = MNXYZS + WYZSOM + 3 * NBSOM
C
C     ALLONGEMENT DU TMS 'NSEF' POUR LES EF AVEC TANGENTES
C     ICI TOUT EF EST DECLARE ETRE AVEC TG
      NBT = LDAPEF + NBEFOB * ( 2 + NBTGEF )
      CALL TAMSAU( NTNSEF, NBT )
      CALL TAMSOU( NTNSEF, MNNSEF )
      MNEFAP = MNNSEF + LDAPEF - 1
      MNCGEF = MNEFAP + NBEFOB
      MNNUTG = MNCGEF + NBEFOB
      MNT    = MNNUTG
C
C     LE NOMBRE ACTUEL DE TANGENTES DU MAILLAGE APRES TRANSFORMATION
      NBTG1 = 0
C     ADRESSE DES 3 COMPOSANTES DE LA NOUVELLE TG A STOCKER
      MNXYT1 = MNXYTG
      DO 140 NBT=1, NBEFOB
C
C        LE NUMERO DES SOMMETS ET TGS DU EF NBT
         CALL NSEFNS( NBT   , NUTYMA, NBSOEF, NBTGEF,
     %                LDAPEF, LDNGEF, LDTGEF,
     %                MNNSEF, NX, NY, NZ,
     %                NCOGEL, NUGEEF, NUEFTG, NOSOEL, IERR )
C
C        LE NOMBRE D'ARETES DU EF NBT
         NARET = NBARET( NCOGEL )
C        LE NUMERO DES 2 SOMMETS DES NARET DU EF DE TYPE NCOGEL
         CALL SOARFA( NCOGEL, NOSOAR, NOSOFA )
C
C        LA BOUCLE SUR LES ARETES DU EF NBT
         DO 135 I=1,NARET
C
C           LE NUMERO DANS XYZSOMMET DES 2 SOMMETS DU ARETE I
            NOST(1) = NOSOEL( NOSOAR(1,I) )
            NOST(2) = NOSOEL( NOSOAR(2,I) )
C
C           LE NUMERO DES 2 TANGENTES DE L'ARETE I D'UN EF DE TYPE NCOGEL
            CALL TGAREF( NCOGEL, I, NOTGAR )
C
            IF( NUEFTG .LE. 0 ) THEN
C              EF SANS TG : NUMERO NUL DES 2 TGS DE L'ARETE INITIALE
               NOTG(1) = 0
               NOTG(2) = 0
            ELSE
C              EF AVEC TGS : NUMERO NUL OU NON DES 2 TGS DE L'ARETE INITIALE
C              LE NUMERO DES 2 TGS DANS XYZSOMMET
               NOTG(1) = NOSOEL( NBSOEF + NOTGAR(1) )
               NOTG(2) = NOSOEL( NBSOEF + NOTGAR(2) )
            ENDIF
C
C           LES COMPOSANTES DE LA DIRECTION DES 2 TANGENTES DANS LE MAILLAGE INI
            DO 132 K=1,2
C
               IF( NOTG(K) .EQ. 0 ) THEN
C                 LA DIRECTION EST CELLE DE +- L'ARETE DROITE
                  MN1 = MNXYZS + WYZSOM + 3*NOST(1) -3
                  MN2 = MNXYZS + WYZSOM + 3*NOST(2) -3
                  IF( K .EQ. 1 ) THEN
C                    TANGENTE DE DIRECTION S1 VERS S2
                     LSIGNE = 1
                  ELSE
C                    TANGENTE DE DIRECTION S2 VERS S1
                     LSIGNE = -1
                  ENDIF
                  DIRTG(1) = LSIGNE * (RMCN(MN2  ) - RMCN(MN1  ))
                  DIRTG(2) = LSIGNE * (RMCN(MN2+1) - RMCN(MN1+1))
                  DIRTG(3) = LSIGNE * (RMCN(MN2+2) - RMCN(MN1+2))
               ELSE
C                 LA DIRECTION EST LA TANGENTE NOTG(K) AU SOMMET NOST(K)
                  MN1 = MNXYZS + WYZSOM + 3*(NBSOM+ABS(NOTG(K))) - 3
                  IF( NOTG(K) .LT. 0 ) THEN
                     LSIGNE = -1
                  ELSE
                     LSIGNE = 1
                  ENDIF
                  DIRTG(1) = LSIGNE * RMCN(MN1  )
                  DIRTG(2) = LSIGNE * RMCN(MN1+1)
                  DIRTG(3) = LSIGNE * RMCN(MN1+2)
               ENDIF
C
C              CALCUL DE LA TANGENTE SELON CETTE DIRECTION
C              EN FAISANT LE PRODUIT DU [GRADIENT] {DIRECTION}
C
C              UNE TANGENTE DE PLUS
               NBTG1 = NBTG1 + 1
C              LE NUMERO DE LA TG K DE CETTE ARETE
               MCN( MNT + NOTGAR(K) ) = NBTG1
C
C              ADRESSE GRADIENT DU SOMMET
               MN2 = MNDF1 + 9*NOST(K) - 10
C
C              LES 3 COORDONNEES DE LA TANGENTE
               RMCN(MNXYT1  ) = RMCN(MN2+1) * DIRTG(1)
     %                        + RMCN(MN2+2) * DIRTG(2)
     %                        + RMCN(MN2+3) * DIRTG(3)
               RMCN(MNXYT1+1) = RMCN(MN2+4) * DIRTG(1)
     %                        + RMCN(MN2+5) * DIRTG(2)
     %                        + RMCN(MN2+6) * DIRTG(3)
               RMCN(MNXYT1+2) = RMCN(MN2+7) * DIRTG(1)
     %                        + RMCN(MN2+8) * DIRTG(2)
     %                        + RMCN(MN2+9) * DIRTG(3)
               MNXYT1 = MNXYT1 + 3
 132        CONTINUE
 135     CONTINUE
C
C        LE NUMERO DE POINTEUR SUR L'EF A TG
         MCN( MNEFAP+NBT ) = NBT
C
C        LE CODE GEOMETRIQUE DU EF A TG EST INCONNU => 0
         MCN( MNCGEF+NBT ) = 0
C
C        L'ADRESSE DES FUTURES TANGENTES DE L'EF SUIVANT
         MNT = MNT + NBTGEF
 140  CONTINUE
C
C     COPIE DANS XYZSOMMET DES COORDONNEES DES SOMMETS APRES TRANSFORMATION
      CALL TRTATA( MCN(MNST1), MCN(MNXYZS+WYZSOM), 3*NBSOM )
C
C     LE NOMBRE DE TANGENTES STOCKEES DANS XYZSOMMET
      MCN( MNXYZS + WNBTGS ) = NBTG1
C     LE NOMBRE DE COORDONNEES PAR SOMMET
      MCN( MNXYZS + WBCOOR ) = 3
C     LA DATE
      CALL ECDATE( MCN( MNXYZS ) )
C
C     LE NOMBRE D'EF A TG
      MCN( MNNSEF + WBTGEF ) = NBTGEF
C     LE NOMBRE D'EF A TG ET POINTEURS
      MCN( MNNSEF + WBEFTG ) = NBEFOB
      MCN( MNNSEF + WBEFAP ) = NBEFOB
C     LE TYPE NUTYMA DU MAILLAGE RESTE INCHANGE
C     LA DATE
      CALL ECDATE( MCN( MNNSEF ) )
C
C     DESTRUCTION DES TABLEAUX DEVENUS INUTILES
      CALL TNMCDS( 'REEL', 3*NBSOM, MNST1 )
      CALL TNMCDS( 'REEL', 9*NBSOM, MNDF1 )
C
C     IDENTIFIER LES TANGENTES EGALES EN CHACUN DES SOMMETS DU MAILLAGE
C     SUPPRIMER LES TANGENTES DOUBLES EN UN SOMMET
      CALL MOINTG( NTNSEF, MNNSEF, NTXYZS, MNXYZS, IERR )
      GOTO 8000
C
C     ISOMETRIE SUR LES SOMMETS DU PLSV INITIAL
C     --------------------------------------------
C     CONSTRUCTION DE LA MATRICE D'ISOMETRIE
  150 CALL ISOMEX( MNDFTR, DMATRI, IERR )
      IF( IERR .NE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) =  'trplsv:ISOMETRIE INCONNUE'
         CALL LEREUR
         IERR = 1
         GOTO 9999
      ENDIF
      DO 160 I=1,NBSOM
C        LES COORDONNEES DU SOMMET I
         P1(1) = RMCN( MNSO )
         P1(2) = RMCN( MNSO + 1 )
         P1(3) = RMCN( MNSO + 2 )
C        TRANSFORMATION DE CE SOMMET PAR L'ISOMETRIE
         CALL ISOMVA( DMATRI, P1, P2 )
C        LA SAUVEGARDE DES COORDONNEES APRES TRANSFORMATION
         RMCN( MNSO     ) = P2( 1 )
         RMCN( MNSO + 1 ) = P2( 2 )
         RMCN( MNSO + 2 ) = P2( 3 )
         MNSO = MNSO + 3
 160  CONTINUE
C
C     LA TRANSFORMATION DES TANGENTES
      IF( NBTGS  .LE. 0 ) GOTO 8000
C     ADRESSE DES COMPOSANTES DES TANGENTES
      MNSO = MNXYZS + WYZSOM + 3 * NBSOM
      DO 170 I=1,NBTGS
C        LES COMPOSANTES DE LA TG I
         P1(1) = RMCN( MNSO )
         P1(2) = RMCN( MNSO + 1 )
         P1(3) = RMCN( MNSO + 2 )
C        TRANSFORMATION DE CETTE TANGENTE PAR L'ISOMETRIE
C        TG2(3) = DMATRI(3,3) * TG1(3)
         DO 168 J=1,3
            D = 0.D0
            DO 164 K=1,3
               D = D + DMATRI(J,K) * P1(K)
 164        CONTINUE
            P2(J) = REAL( D )
 168     CONTINUE
C        LA SAUVEGARDE DES COORDONNEES APRES TRANSFORMATION
         RMCN( MNSO     ) = P2( 1 )
         RMCN( MNSO + 1 ) = P2( 2 )
         RMCN( MNSO + 2 ) = P2( 3 )
         MNSO = MNSO + 3
 170  CONTINUE
C
C     LE NOUVEAU NUMERO GEOMETRIQUE
      D = ABS( DMATRI(1,1) + DMATRI(1,2) + DMATRI(1,3) )
      IF( ABS( ABS( DMATRI(2,1) + DMATRI(2,2) + DMATRI(2,3) ) - D )
     %    .GT. D*1D-7 ) THEN
C        DILATATION OU REDUCTION
         NEWGEO = 2
         GOTO 8000
      ENDIF
      IF( ABS( ABS( DMATRI(3,1) + DMATRI(3,2) + DMATRI(3,3) ) - D )
     %    .GT. D*1D-7 ) THEN
C        DILATATION OU REDUCTION
         NEWGEO = 2
         GOTO 8000
      ENDIF
C     ISOMETRIE REELLE
      NEWGEO = 1
      GOTO 8000
C
C     PROJECTION DES SOMMETS DU PLSV INITIAL
C     -----------------------------------------
C     LE TYPE DE LA PROJECTION EST IL CORRECT ?
 180  NUTYPR = MCN( MNDFTR + WUTYPR )
      NTYPRA = ABS( NUTYPR )
      LESIGN = 1
      NUFONC = MCN( MNDFTR + WUFONC )
      IF( NUTYPR .LT. -3 .OR. NUTYPR .GT. 6 .OR. NUTYPR .EQ. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'trplsv:PROJECTION INCONNUE'
         CALL LEREUR
         IERR = 4
         GOTO 9999
      ELSE IF( MCN( MNDFTR+WUFONC ) .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'trplsv:SURFACE F(X,Y,Z)=0 INCONNUE'
         CALL LEREUR
         IERR = 5
         GOTO 9999
      ENDIF
C
      IF( NUTYPR .EQ. 4 ) THEN
         NUPTCP = MCN( MNDFTR + WUPTCP )
         IF( NUPTCP .LE. 0 ) THEN
            NBLGRC(NRERR) = 1
            KERR(1) = 'trplsv:CENTRE DE PROJECTION POINT INCONNU'
            CALL LEREUR
            IERR = 6
            GOTO 9999
         ENDIF
      ELSE IF( NUTYPR .EQ. 5 ) THEN
C        RECHERCHE DE LA DIRECTION XYZ ORTHOGONALE AU PLAN
C        DEFINI PAR  3 POINTS
         IF( MCN( MNDFTR+WUPLP1) .LE. 0 ) THEN
            NBLGRC(NRERR) = 1
            KERR(1) = 'trplsv: PLAN PROJECTION:3 POINTS INCONNUS'
            CALL LEREUR
            IERR = 10
            GOTO 9999
         ENDIF
      ELSE IF( NUTYPR .EQ. 6 ) THEN
C        RECHERCHE DE LA DIRECTION XYZ ORTHOGONALE A LA DROITE
C        DEFINIE PAR  2 POINTS
         IF( MCN( MNDFTR+WUDRP1) .LE. 0 ) THEN
            NBLGRC(NRERR) = 1
            KERR(1) =
     %     'trplsv:PROJECTION ORTHOGONALE A LA DROITE:2 POINTS INCONNUS'
            CALL LEREUR
            IERR = 11
            GOTO 9999
         ENDIF
      ENDIF
C
C     UNE LONGUEUR DE REFERENCE
C     RECHERCHE DU MIN_MAX DES SOMMETS DES 2 PLSVS
      CALL CADEXT( MNXYZS, COIN )
      DMARET = ( COIN(1,2) - COIN(1,1) +
     %           COIN(2,2) - COIN(2,1) +
     %           COIN(3,2) - COIN(3,1) ) * 0.1
      IF( DMARET .LE. 0. ) DMARET = 1.
C
C     LA DIRECTION DE PROJECTION SUIVANT LES AXES
      IF( NUTYPR .GE. 0 ) THEN
         LADIR = 1
      ELSE
         LADIR = -1
      ENDIF
C
C     LE MODE DE PROJECTION
      GOTO( 220, 220, 220, 240, 250, 255 ), NTYPRA
C
C     PROJECTION SELON OX OY OZ
 220  DO 230 I=1,3
         XYZN(I) = 0.
 230   CONTINUE
      XYZN( NTYPRA ) = DMARET * LADIR
      GOTO 260
C
C     LE POINT CENTRE DEFINISSANT LE PLAN AVEC LES 2 EXTREMITES
 240  CALL XYZPOI( MCN(MNDFTR+WUPTCP), MN, IERR )
      IF( IERR .NE. 0 ) GOTO 9999
      XYZN(1) = RMCN( MN     )
      XYZN(2) = RMCN( MN + 1 )
      XYZN(3) = RMCN( MN + 2 )
      DMARET  = 1.
      GOTO 260
C
C     CALCUL DES 3 COORDONNEES XYZN DE LA NORMALE AU PLAN
C     DEFINI PAR SES 3 POINTS DE NUMERO NUPLP1,2,3
 250  CALL NOPLAN( MCN( MNDFTR+WUPLP1), XYZN, IERR )
      IF( IERR .NE. 0 ) GOTO 9999
      GOTO 260
C
C     CALCUL DES 3 COORDONNEES DES 2 POINTS DEFINISSANT LA DROITE
C     LE PREMIER POINT
 255  CALL XYZPOI( MCN(MNDFTR+WUDRP1), MN, IERR )
      IF( IERR .NE. 0 ) GOTO 9999
      XYZN(1) = RMCN( MN     )
      XYZN(2) = RMCN( MN + 1 )
      XYZN(3) = RMCN( MN + 2 )
C     LE SECOND POINT
      CALL XYZPOI( MCN(MNDFTR+WUDRP2), MN, IERR )
      IF( IERR .NE. 0 ) GOTO 9999
C     CES 2 POINTS SONT ILS IDENTIQUES ?
      CALL XYZIDE( XYZN, RMCN(MN), I )
      IF( I .EQ. 1 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'LA DROITE EST DEFINIE PAR 2 POINTS CONFONDUS'
         CALL LEREUR
         IERR = 13
         GOTO 9999
      ENDIF
C     LES 2 POINTS SONT DIFFERENTS
      XYZ(1) = RMCN( MN     ) - XYZN(1)
      XYZ(2) = RMCN( MN + 1 ) - XYZN(2)
      XYZ(3) = RMCN( MN + 2 ) - XYZN(3)
C     LE CARRE DE LA LONGUEUR DU SEGMENT DEFINISSANT L'AXE
      D      = XYZ(1)**2 + XYZ(2)**2 + XYZ(3)**2
      DMARET = 2
C
C     LA BOUCLE SUR LES NBS SOMMETS DU PLSV
C     DMARET LA DISTANCE MOYENNE ENTRE LES SOMMETS DES ARETES
 260  DARET = DMARET
      DO 395 I=1,NBSOM
         DMARET = DARET
C        ITER NOMBRE D'ITERATIONS POUR ENCADRER LA SURFACE DE PROJECTION
         ITER   = 0
C        RECHERCHE DE 2 POINTS P1 ET P2 ENCADRANT LA SURFACE
C        C-A-D DE SIGNES OPPOSES POUR LA FONCTION NUFONC
C        P1 EST INITIALEMENT LE POINT A PROJETER
C        LES COORDONNEES DU SOMMET I
         P1(1) = RMCN( MNSO )
         P1(2) = RMCN( MNSO + 1 )
         P1(3) = RMCN( MNSO + 2 )
C
 305     ITER = ITER + 1
         GOTO( 310, 310, 310, 340, 350, 357 ), NTYPRA
C
C        PROJECTION OX OY OZ
 310     DO 320 J=1,3
            P2(J) = P1(J)
 320     CONTINUE
         P2( NTYPRA ) = P1( NTYPRA ) + ITER * XYZN( NTYPRA )
         GOTO 360
C
C        PROJECTION SELON LE CENTRE XYZN ET LE POINT INITIAL P1
 340     DO 345 J=1,3
            P2(J) = P1(J) + ( P1(J)-XYZN(J) ) * ITER * DMARET
 345     CONTINUE
         GOTO 360
C
C        PROJECTION ORTHOGONALE SELON LA NORMALE XYZN
 350     DO 355 J=1,3
            P2(J) = P1(J) + XYZN(J) * ITER * DMARET
 355     CONTINUE
         GOTO 360
C
C        PROJECTION ORTHOGONALE A UNE DROITE
C        RECHERCHE DU POINT P3 PROJECTION ORTHOGONALE DE P1 SUR LA DROITE
 357     IF( ITER .EQ. 1 ) THEN
            DD = REAL( ( ( XYZN(1)-P1(1) ) * XYZ(1) +
     %                   ( XYZN(2)-P1(2) ) * XYZ(2) +
     %                   ( XYZN(3)-P1(3) ) * XYZ(3) ) / D )
            P3(1) = XYZN(1) - XYZ(1) * DD
            P3(2) = XYZN(2) - XYZ(2) * DD
            P3(3) = XYZN(3) - XYZ(3) * DD
         ENDIF
         DO 359 J=1,3
            P2(J) = P3(J) + (P1(J)-P3(J)) * ITER * DMARET
 359     CONTINUE
C
C        CALCUL DE LA FONCTION DE CES 2 POINTS
 360     DP(1) = P1(1)
         DP(2) = P1(2)
         DP(3) = P1(3)
         CALL FONVAL( NUFONC, 3, DP, NCODEV, DBLVAL )
         IF( NCODEV .LE. 0 ) THEN
            IERR = 12
            GOTO 9999
         ENDIF
         F1 = REAL( DBLVAL )
C
         DP(1) = P2(1)
         DP(2) = P2(2)
         DP(3) = P2(3)
         CALL FONVAL( NUFONC, 3, DP, NCODEV, DBLVAL )
         IF( NCODEV .LE. 0 ) THEN
            IERR = 12
            GOTO 9999
         ENDIF
         F2 = REAL( DBLVAL )
C
         IF( ABS(F1) .LT. EPZERO ) GOTO 390
         IF( ABS(F2) .LT. EPZERO ) THEN
C           P2 EST SOLUTION
            DO 370 J=1,3
               P1(J) = P2(J)
 370        CONTINUE
            GOTO 390
         ENDIF
C
         IF( F1 * F2 .GT. 0. ) THEN
C           LES 2 POINTS NE CONVIENNENT PAS
            IF( ITER .LE. MXITER ) GOTO 305
            IF( LESIGN .GT. 0 .AND. NTYPRA .LE. 3 ) THEN
                XYZN( NTYPRA ) = - LADIR * DMARET / MXITER
                LESIGN = -LESIGN
                ITER   = 0
                GOTO 305
            ENDIF
            IF( (NTYPRA.EQ.4 .OR. NTYPRA.EQ.6) .AND. DMARET.GT.0. ) THEN
C              ON INVERSE LE SENS DE LA RECHERCHE DU POINT P2 PAR RAPPORT A P1
               DMARET = -DMARET / MXITER
               IF( NTYPRA .EQ. 6 ) DMARET = -1. / ( MXITER + 1 )
               ITER   = 0
               GOTO 305
            ENDIF
C           TOUTES LES TENTATIVES ONT ECHOUE
            NBLGRC(NRERR) = 3
            KERR(1) = 'PROJECTION IMPOSSIBLE'
            KERR(2) = 'REVOIR LA FONCTION SURFACE F(X,Y,Z)=0 '
            KERR(3) = 'ET LES POINTS A PROJETER'
            CALL LEREUR
            IERR = 11
            GOTO 9999
         ENDIF
C
C        LES 2 POINTS ENCADRENT LA SURFACE
C        ---------------------------------
C        BISSECTION SUR LA DROITE P1P2
         CALL BISSEC( NUFONC, P1, F1, P2, F2, IERR )
         IF( IERR .NE. 0 ) GOTO 9999
C        LA SAUVEGARDE DES COORDONNEES APRES TRANSFORMATION
 390     RMCN( MNSO     ) = P1( 1 )
         RMCN( MNSO + 1 ) = P1( 2 )
         RMCN( MNSO + 2 ) = P1( 3 )
         MNSO = MNSO + 3
 395  CONTINUE
C
C     LE NUMERO GEOMETRIQUE EST INCONNU => CAS STANDARD
      NEWGEO = 0
C
C     TANGENTES INCALCULABLES DONC SUPPRIMEES
      CALL SUPTGS( NTXYZS, MNXYZS, NTNSEF, MNNSEF )
C
C     MISE A JOUR DE LA DATE DU TABLEAU XYZSOMMET
C     ===========================================
 8000 CALL ECDATE( MCN( MNXYZS ) )
      IF( NUTYOB .LE. 1 ) GOTO 9999
C
C     LES PARAMETRES DU PLSV
      CALL NSEFPA( MCN(MNNSEF),
     %             NUTYMA, NBSOEL, NBSOEF, NBTGEF,
     %             LDAPEF, LDNGEF, LDTGEF, NBEFOB,
     %             NX   , NY   , NZ   ,
     %             IERR   )
C
      IF( NUTYOB .LE. 3 ) GOTO 8500
C
C     VERIFICATION DU VOLUME POSITIF DES TETRAEDRES, PYRAMIDES,
C     PENTAEDRES ET HEXAEDRES
C     =========================================================
      MNC = MNNSEF + WUSOEF - 1
      MNX = MNXYZS + WYZSOM - 4
      DO 8100 NBT=1, NBEFOB
C
C        LE NUMERO DES SOMMETS DU ELEMENT NBT
         CALL NSEFNS( NBT   , NUTYMA, NBSOEF, NBTGEF,
     %                LDAPEF, LDNGEF, LDTGEF,
     %                MNNSEF, NX, NY, NZ,
     %                NCOGEL, NUGEEF, NUEFTG, NOSOEL, IERR )
         IF( NOSOEL(5) .EQ. 0 ) THEN
C           TETRAEDRE
            NS = 4
         ELSE IF( NOSOEL(6) .EQ. 0 ) THEN
C           PYRAMIDE
            NS = 5
         ELSE IF( NOSOEL(7) .EQ. 0 ) THEN
C           PENTAEDRE
            NS = 4
         ELSE
C           HEXAEDRE
            NS = 5
         ENDIF
C        12 PRODUIT VECTORIEL 13
         DO 8050 I=1,3
            V      = RMCN(MNX+3*NOSOEL( 1)+I)
            P1(I)  = RMCN(MNX+3*NOSOEL( 2)+I) - V
            P2(I)  = RMCN(MNX+3*NOSOEL( 3)+I) - V
            XYZ(I) = RMCN(MNX+3*NOSOEL(NS)+I) - V
 8050    CONTINUE
         CALL PROVER( P1, P2, P3 )
C        ( 12 PRODUIT VECTORIEL 13 ) PRODUIT SCALAIRE 14 ou 15
         V = PROSCR( P3, XYZ, 3 )
         IF( V .LT. 0 ) THEN
C           VOLUME NEGATIF
            IF( NUTYMA .GT. 0 ) THEN
C              MAILLAGE STRUCTURE => DESTRUCTURATION DU MAILLAGE
               CALL STNOST( NTLXOB, NTNSEF, MNNSEF, IERR )
               IF( IERR .NE. 0 ) THEN
                  NBLGRC(NRERR) = 3
                  KERR(1) = 'TRANSFORMATION RENDANT UN EF DE VOLUME<0'
                  WRITE(KERR(3)(1:9),'(I9)') NBT
                  KERR(2) = 'EF NUMERO ' // KERR(3)(1:9)
                  KERR(3) = 'MAILLAGE STRUCTURE NON DESTRUCTURABLE'
                  CALL LEREUR
                  IERR = 20
                  GOTO 9999
               ENDIF
C              ON REDEMARRE AVEC LE MAILLAGE NON STRUCTURE
               GOTO 8000
            ELSE
C              MAILLAGE NON STRUCTURE
C              PERMUTATION DES SOMMETS POUR RENDRE LE VOLUME>0
               IF( MCN(MNC+5) .EQ. 0 ) THEN
C                 TETRAEDRE
                  NS         = MCN(MNC+2)
                  MCN(MNC+2) = MCN(MNC+3)
                  MCN(MNC+3) = NS
               ELSE IF( MCN(MNC+6) .EQ. 0 ) THEN
C                 PYRAMIDE
                  NS         = MCN(MNC+2)
                  MCN(MNC+2) = MCN(MNC+4)
                  MCN(MNC+4) = NS
               ELSE IF( MCN(MNC+7) .EQ. 0 ) THEN
C                 PENTAEDRE
                  NS         = MCN(MNC+2)
                  MCN(MNC+2) = MCN(MNC+3)
                  MCN(MNC+3) = NS
                  NS         = MCN(MNC+5)
                  MCN(MNC+5) = MCN(MNC+6)
                  MCN(MNC+6) = NS
               ELSE
C                 HEXAEDRE
                  NS         = MCN(MNC+2)
                  MCN(MNC+2) = MCN(MNC+4)
                  MCN(MNC+4) = NS
                  NS         = MCN(MNC+6)
                  MCN(MNC+6) = MCN(MNC+8)
                  MCN(MNC+8) = NS
               ENDIF
            ENDIF
         ENDIF
         MNC = MNC + 8
 8100 CONTINUE
C
C     TRAITEMENT DU NUMERO GEOMETRIQUE DES EF A TG
C     ============================================
 8500 NBEFTG = MCN( MNNSEF + WBEFTG )
      IF( NBEFTG .GT. 0 ) THEN
C        IL Y A DES EF A TG
         MNC = MNNSEF + LDNGEF - 1
         DO 8600 NBT=1, NBEFTG
C           LE NUMERO GEOMETRIQUE DE CET EF A TG AVANT TRANSFORMATION
            NUGEEF = MCN( MNC + NBT )
            IF( NEWGEO .EQ. 0 ) THEN
C              APRES TRANSFORMATION: LE NUMERO GEOMETRIQUE EST FORCE A LA VALEUR
C              EF A TG STANDARD ( ARC P3 HERMITE ou TRIANGLE HCT ou QUADRANGLE d
               NUGEEF = 0
            ELSE IF( NEWGEO .EQ. 2 ) THEN
C              EF A TG AYANT SUBI UNE DILATATION OU REDUCTION
               IF( NUGEEF .EQ. 1 ) THEN
C                 CERCLE:1 DEVIENT ELLIPSE:2
                  NUGEEF = 2
               ELSE IF( NUGEEF .EQ. 11 ) THEN
C                 SPHERE:11 DEVIENT ELLIPSOIDE:12
                  NUGEEF = 12
C              ELSE
                  NUGEEF = 0
C                 A COMPLETER AU FUR ET A MESURE ...
               ENDIF
            ENDIF
C           LE NOUVEAU NUMERO GEOMETRIQUE APRES TRANSFORMATION
            MCN( MNC + NBT ) = NUGEEF
 8600    CONTINUE
      ENDIF
      IERR = 0
C
C     SUPPRESSION DES TANGENTES SI INCOHERENCE ENTRE LE NOMBRE D'EF A TG
C     ET LE NOMBRE DE COMPOSANTES DES TANGENTES STOCKEES
C     ==================================================================
      NBTGS  = MCN( MNXYZS + WNBTGS )
      NBEFTG = MCN( MNNSEF + WBEFTG )
      IF( (NBEFTG .LE. 0 .AND. NBTGS .GT. 0)  .OR.
     %    (NBEFTG .GT. 0 .AND. NBTGS .LE. 0) ) THEN
         CALL SUPTGS( NTXYZS, MNXYZS, NTNSEF, MNNSEF )
      ENDIF
C
C     ERREUR
C     ======
 9999 RETURN
      END
