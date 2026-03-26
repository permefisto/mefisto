      SUBROUTINE SUTOTG( NUSUIN, NTLXSU,
     %                   NTNSEF, MNNSEF, NTXYZS, MNXYZS, IERR )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    COPIER LES TMS NSEF ET XYZSOMMET D'UNE SURFACE NUSUIN
C -----    DE TELLE SORTE QUE
C          TOUT EF EST UN EF A TG
C          TOUTE TANGENTE EST UNIQUE
C          TOUTE ARETE SANS TG DEVIENT UNE TANGENTE
C
C ENTREES:
C --------
C NUSUIN : NUMERO DE LA SURFACE INITIALE DANS LE LEXIQUE DES SURFACES
C NTLXSU : NUMERO DU TABLEAU TS DU LEXIQUE DE LA SURFACE FINALE
C
C SORTIES:
C --------
C NTNSEF : NUMERO      DU TMS 'NSEF' DES NUMEROS DES FACES DE LA SURFACE
C MNNSEF : ADRESSE MCN DU TMS 'NSEF' DES NUMEROS DES FACES DE LA SURFACE
C          CF $MEFISTO/td/d/a___nsef
C NTXYZS : NUMERO      DU TMS 'XYZSOMMET' DE LA SURFACE
C MNXYZS : ADRESSE MCN DU TMS 'XYZSOMMET' DE LA SURFACE
C          CF $MEFISTO/td/d/a___xyzsommet
C IERR   : 0 SI PAS D'ERREUR
C          > 0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC PARIS        JUIN 1999
C ...................................................................012
      include"./incl/ntmnlt.inc"
      include"./incl/gsmenu.inc"
      include"./incl/a_surface__definition.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
      INTEGER           NOSOEL(12)
C
C     LA SURFACE A TRAITER
C     ====================
C     LE TABLEAU LEXIQUE DE CETTE SURFACE
      CALL LXNLOU( NTSURF, NUSUIN, NTLXSF, MNLXSF )
      IF( NTLXSF .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'SURFACE INCONNUE'
         CALL LEREUR
         IERR = 1
         RETURN
      ENDIF
C
C     LE TABLEAU 'NSEF' DE CETTE SURFACE
      CALL LXTSOU( NTLXSF, 'NSEF', NTNSE0, MNNSE0 )
      IF( NTNSE0 .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'SURFACE SANS NSEF'
         CALL LEREUR
         IERR = 2
         RETURN
      ENDIF
C
C     LES PARAMETRES DES NO SOMMET DU MAILLAGE DE LA SURFACE
      CALL NSEFPA( MCN(MNNSE0),
     %             NUTYMA, NBSOEL, NBSOEF, NBTGEF,
     %             LDAPEF, LDNGEF, LDTGEF, NBEFOB,
     %             NX, NY, NZ,
     %             IERR   )
      IF( IERR .NE. 0 ) RETURN
C
C     LE NOMBRE D'EF A TG DE LA SURFACE
      NBEFT0 = MCN( MNNSE0 + WBEFTG )
C     LE NOMBRE D'EF
      NBEFOB = MCN( MNNSE0 + WBEFOB )
C
C     LE TABLEAU 'XYZSOMMET' DE CETTE SURFACE
      CALL LXTSOU( NTLXSF, 'XYZSOMMET', NTXYZ0, MNXYZ0 )
      IF( NTXYZ0 .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'SURFACE SANS XYZSOMMET'
         CALL LEREUR
         IERR = 3
         RETURN
      ENDIF
C     LE NOMBRE DE SOMMETS
      NBSOM0 = MCN(MNXYZ0+WNBSOM)
C     LE NOMBRE INITIAL DE TGS
      NBTGS0 = MCN(MNXYZ0+WNBTGS)
C
      IF( NBEFT0 .LE. 0 .OR. NBTGS0 .LE. 0 ) THEN
         NBLGRC(NRERR) = 2
         KERR(1) = 'SURFACE INITIALE SANS TANGENTES'
         KERR(2) = 'IMPOSSIBLE DE RENDRE G1-CONTINUE CETTE SURFACE'
         CALL LEREUR
         IERR = 3
         RETURN
      ENDIF
C
C     TOUT EF DEVIENT UN EF A TG POUR RESOUDRE LA G1-CONTINUITE
C     TOUTE TG MULTIPLE REDEVIENT UNE TANGENTE UNIQUE
C     ---------------------------------------------------------
C     LE TABLEAU DES XYZ DES TANGENTES SIMPLES
      CALL TNMCDC( 'REEL', 3*8*NBEFOB, MNCOTG )
C     LE POINTEUR EF A TG, LE CODE GEOMETRIQUE ET LES NO DES TGS DES EF A TG
      CALL TNMCDC( 'ENTIER', NBEFOB,   MNEFAP )
      CALL TNMCDC( 'ENTIER', NBEFOB,   MNEFCG )
      CALL TNMCDC( 'ENTIER', NBEFOB*8, MNNUTG )
C
      MNUTG  = MNNUTG
      MNXYT0 = MNXYZ0 + WYZSOM + 3*NBSOM0 - 4
      MNXYTG = MNCOTG - 1
      MNXYST = MNXYZ0 + WYZSOM - 4
C
      NBTGS  = 0
      DO 90 NUEF=1,NBEFOB
C
C        LE NUMERO DES SOMMETS DE L'EF NUEF
         CALL NSEFNS( NUEF,   NUTYMA, NBSOEF, NBTGEF,
     %                LDAPEF, LDNGEF, LDTGEF,
     %                MNNSE0, NX, NY, NZ,
     %                NCOGEL, NUGEEF, NUEFTG, NOSOEL, IERR )
         IF( IERR .NE. 0 ) GOTO 9000
C
C        UN EF A TG DE PLUS
         MCN( MNEFAP -1 + NUEF ) = NUEF
C        CODE GEOMETRIQUE DE CET EF
         MCN( MNEFCG -1 + NUEF ) = NUGEEF
C
         IF( NUEFTG .LE. 0 ) THEN
C
C           L'EF DEVIENT UN EF A TG
C           -----------------------
C           CREATION DES TANGENTES
            I0 = NCOGEL
            DO 20 I=1,NCOGEL
               IF( I .LT. NCOGEL ) THEN
                  I1 = I + 1
               ELSE
                  I1 = 1
               ENDIF
C              LES 3 COORDONNEES DES ARETES ISSUES DU SOMMET I
               DO 10 K=1,3
                  RMCN(MNXYTG+K  ) = RMCN(MNXYST+3*NOSOEL(I1)+K)
     %                             - RMCN(MNXYST+3*NOSOEL(I )+K)
                  RMCN(MNXYTG+3+K) = RMCN(MNXYST+3*NOSOEL(I0)+K)
     %                             - RMCN(MNXYST+3*NOSOEL(I )+K)
 10            CONTINUE
               MNXYTG = MNXYTG + 6
C
C              LE NUMERO DES TGS DU SOMMET I
               NBTGS = NBTGS + 1
               MCN(MNUTG) = NBTGS
               NBTGS = NBTGS + 1
               MCN(MNUTG+1) = NBTGS
               MNUTG = MNUTG + 2
C
               I0    = I
 20         CONTINUE
C
         ELSE
C
C           L'EF RESTE UN EF A TG A COMPLETER
C           ---------------------------------
C           CREATION DES TANGENTES MANQUANTES
            I0 = NCOGEL
            L  = 4
            DO 70 I=1,NCOGEL
C
C              L'ARETE I->I+1
               L = L + 1
               NOTG = NOSOEL(L)
               IF( NOTG .NE. 0 ) THEN
C
C                 LE NUMERO DE LA TG EST CHANGEE
                  NBTGS = NBTGS + 1
                  MCN(MNUTG) = NBTGS
                  IF( NOTG .GT. 0 ) THEN
                     LSIGNE = 1
                  ELSE
                     LSIGNE = -1
                  ENDIF
                  DO 30 K=1,3
                     RMCN(MNXYTG+K) = LSIGNE*RMCN(MNXYT0+3*ABS(NOTG)+K)
 30               CONTINUE
                  MNXYTG = MNXYTG + 3
               ELSE
C                 LA TG ARETE EST AJOUTEE
                  IF( I .LT. NCOGEL ) THEN
                     I1 = I + 1
                  ELSE
                     I1 = 1
                  ENDIF
C                 LES 3 COORDONNEES DE L'ARETE I->I1
                  DO 40 K=1,3
                     RMCN(MNXYTG+K) = RMCN(MNXYST+3*NOSOEL(I1)+K)
     %                              - RMCN(MNXYST+3*NOSOEL(I )+K)
 40               CONTINUE
                  MNXYTG = MNXYTG + 3
C                 LE NUMERO DE CETTE TG
                  NBTGS = NBTGS + 1
                  MCN(MNUTG) = NBTGS
               ENDIF
               MNUTG = MNUTG + 1
C
C              L'ARETE I->I-1
               L = L + 1
               NOTG = NOSOEL(L)
               IF( NOTG .NE. 0 ) THEN
C
C                 LE NUMERO DE LA TG EST CHANGEE
                  NBTGS = NBTGS + 1
                  MCN(MNUTG) = NBTGS
                  IF( NOTG .GT. 0 ) THEN
                     LSIGNE = 1
                  ELSE
                     LSIGNE = -1
                  ENDIF
                  DO 45 K=1,3
                     RMCN(MNXYTG+K) = LSIGNE*RMCN(MNXYT0+3*ABS(NOTG)+K)
 45               CONTINUE
                  MNXYTG = MNXYTG + 3
C
               ELSE
C                 LES 3 COORDONNEES DE L'ARETE I->I0
                  DO 50 K=1,3
                     RMCN(MNXYTG+K) = RMCN(MNXYST+3*NOSOEL(I0)+K)
     %                              - RMCN(MNXYST+3*NOSOEL(I )+K)
 50               CONTINUE
                  MNXYTG = MNXYTG + 3
C                 LE NUMERO DE CETTE TG
                  NBTGS = NBTGS + 1
                  MCN(MNUTG) = NBTGS
               ENDIF
               MNUTG = MNUTG + 1
C
               I0    = I
 70         CONTINUE
         ENDIF
C
C        COMPLEMENT A ZERO DES TANGENTES 7 8 DU TRIANGLE
         DO 80 I=NCOGEL+1,4
            MCN(MNUTG) = 0
            MNUTG = MNUTG + 1
            MCN(MNUTG) = 0
            MNUTG = MNUTG + 1
 80      CONTINUE
 90   CONTINUE
C
C     CONSTRUCTION DU TABLEAU 'XYZSOMMET' PAR COPIE
C     ---------------------------------------------
      CALL LXTSOU( NTLXSU, 'XYZSOMMET',  NTXYZS, MNXYZS )
      IF( MNXYZS .GT. 0 ) THEN
         CALL LXTSDS( NTLXSU, 'XYZSOMMET' )
      ENDIF
      L = WYZSOM + 3 * ( NBSOM0 + NBTGS )
      CALL LXTNDC( NTLXSU, 'XYZSOMMET', 'MOTS', L )
      CALL LXTSOU( NTLXSU, 'XYZSOMMET',  NTXYZS, MNXYZS )
C
C     COPIE DU TABLEAU 'XYZSOMMET' INITIAL DANS LE NOUVEAU
      CALL TRTATA( MCN(MNXYZ0+WYZSOM), MCN(MNXYZS+WYZSOM), 3*NBSOM0 )
C     COPIE DU TABLEAU DES COORDONNEES DES TGS
      CALL TRTATA( MCN(MNCOTG), MCN(MNXYZS+WYZSOM+3*NBSOM0), 3*NBTGS )
C
C     LE NOMBRE DE SOMMETS
      MCN(MNXYZS+WNBSOM) = NBSOM0
C     LE NOMBRE EFFECTIF DE TGS
      MCN(MNXYZS+WNBTGS) = NBTGS
C
C     AJOUT DU NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNXYZS + MOTVAR(6) ) = NONMTD( '~>>>XYZSOMMET' )
C     MISE A JOUR DE LA DATE
      CALL ECDATE( MCN(MNXYZS) )
C
C     CONSTRUCTION DU TABLEAU 'NSEF' PAR COPIE
C     ----------------------------------------
      CALL LXTSOU( NTLXSU, 'NSEF', NTNSEF, MNNSEF )
      IF( MNNSEF .GT. 0 ) THEN
          CALL LXTSDS( NTLXSU, 'NSEF' )
      ENDIF
      NBEFTG = NBEFOB
      L = LDAPEF + NBEFOB + NBEFTG*9
      CALL LXTNDC( NTLXSU, 'NSEF', 'MOTS', L )
      CALL LXTSOU( NTLXSU, 'NSEF',  NTNSEF, MNNSEF )
C
C     COPIE DU TABLEAU 'NSEF' INITIAL DANS LE NOUVEAU
      CALL TRTATA( MCN(MNNSE0), MCN(MNNSEF), LDAPEF )
C
C     LES POINTEURS SUR LES EF A TG
      CALL TRTATA( MCN(MNEFAP), MCN(MNNSEF+LDAPEF), NBEFOB )
C     LES CODES GEOMETRIQUES
      CALL TRTATA( MCN(MNEFCG), MCN(MNNSEF+LDAPEF+NBEFOB), NBEFTG )
C     LES NUMEROS DES TGS
      CALL TRTATA( MCN(MNNUTG), MCN(MNNSEF+LDAPEF+NBEFOB+NBEFTG),
     %             8*NBEFTG )
C
C     NOMBRE DE TGS PAR EF
      MCN( MNNSEF+WBTGEF ) = 8
C     NOMBRE D'EF A TGS
      MCN( MNNSEF+WBEFTG ) = NBEFTG
C     NOMBRE D'EF A POINTEUR
      MCN( MNNSEF+WBEFAP ) = NBEFOB
C
C     AJOUT DU NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNNSEF + MOTVAR(6) ) = NONMTD( '~>>>NSEF' )
C     MISE A JOUR DE LA DATE
      CALL ECDATE( MCN(MNNSEF) )
C
C     DESTRUCTION DES TABLEAUX AUXILIAIRES
 9000 CALL TNMCDS( 'REEL', 3*8*NBEFOB, MNCOTG )
      CALL TNMCDS( 'ENTIER', NBEFOB,   MNEFAP )
      CALL TNMCDS( 'ENTIER', NBEFOB,   MNEFCG )
      CALL TNMCDS( 'ENTIER', NBEFOB*8, MNNUTG )
C
      RETURN
      END
