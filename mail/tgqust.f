      SUBROUTINE TGQUST( MNSOLI, MNNTGL,
     %                   NTLXSU, MNSOFA, MNFASU )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  AJOUTER LES TANGENTES AU BORD D'UN QUADRANGLE STRUCTURE
C -----  A PARTIR DES TANGENTES DES LIGNES DES 4 COTES PRESENTEES
C        SELON LE SENS C1:S1S2 C2:S2S3  C3:S4S3 C4:S1S4
C
C                       C3
C                S4---->-------S3
C                |             |
C                |             |
C           C4  / \           / \ C2
C                |             |
C                |             |
C                S1---->-------S2
C                       C1
C
C ENTREES:
C --------
C MNSOLI : ADRESSE MCN DU TABLEAU 'xyzsommet' DES 4 LIGNES COTES
C MNNTGL : ADRESSE MCN DU TABLEAU DES 2 NUMEROS DES TANGENTES DES
C          ARETES DE LA LIGNE  (0 SI LE TABLEAU N'EXISTE PAS)
C          TABLEAU A DETRUIRE EN FIN D'UTILISATION
C NTLXSU : NUMERO DU TABLEAU TS DU LEXIQUE DU QUADRANGLE STRUCTURE
C MNSOFA : ADRESSE MCN DU TABLEAU 'xyzsommet' DE LA SURFACE
C MNFASU : ADRESSE MCN DU TABLEAU 'nsef'      DE LA SURFACE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS     OCTOBRE 1996
C2345X7..............................................................012
      IMPLICIT INTEGER (W)
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
      INTEGER           MNSOLI(4), MNNTGL(4)
      INTEGER           NBTG(4), NOTG(8)
C
C     SI LA SURFACE A DES TANGENTES ALORS RETOUR
      IF( MCN(MNSOFA+WNBTGS) .GT. 0 ) RETURN
C
C     SI LA SURFACE N'EST PAS UN QUADRANGLE STRUCTURE ALORS RETOUR
      IF( MCN(MNFASU+WUTYMA) .NE. 4 ) RETURN
C
C     LES PARAMETRES DES NO SOMMET DU MAILLAGE DE LA SURFACE
      CALL NSEFPA( MCN(MNFASU),
     %             NUTYMA, NBSOEL, NBSOEF, NBTGEF,
     %             LDAPEF, LDNGEF, LDTGEF, NBEFOB,
     %             NX    , NY    , NZ    ,
     %             IERR   )
      IF( IERR .NE. 0 ) RETURN
C
C     CALCUL DU NOMBRE DE TANGENTES A STOCKER POUR LE QUADRANGLE STRUCTURE
      NBTGS   = 0
      DO 5 N=1,4
         NBTGS   = NBTGS + MCN(MNSOLI(N)+WNBTGS)
         NBTG(N) = NBTGS
 5    CONTINUE
      IF( NBTGS .LE. 0 ) RETURN
C
C     TRAITEMENT DU TMS 'XYZSOMMET'
C     =============================
C     SAUVEGARDE DU TABLEAU DES COORDONNEES SANS TG
      NBSOM = MCN(MNSOFA+WNBSOM)
      CALL TNMCDC( 'REEL', 3*NBSOM, MNXYZ0 )
      CALL TRTATA( MCN(MNSOFA+WYZSOM), MCN(MNXYZ0), 3*NBSOM )
C
C     DESTRUCTION DU TMS 'XYZSOMMET' DU QUADRANGLE STRUCTURE
      CALL LXTSDS( NTLXSU, 'XYZSOMMET' )
C
C     CONSTRUCTION DU TABLEAU 'XYZSOMMET' DE CETTE SURFACE
      CALL LXTNDC( NTLXSU, 'XYZSOMMET', 'MOTS',
     %             WYZSOM + 3 * (NBSOM+NBTGS) )
      CALL LXTSOU( NTLXSU, 'XYZSOMMET',  NTSOFA, MNSOFA )
C     LE NOMBRE DE COORDONNEES DES SOMMETS
      MCN( MNSOFA + WBCOOR ) = 3
C     LE NOMBRE DE SOMMETS
      MCN( MNSOFA + WNBSOM ) = NBSOM
C     LE NOMBRE DE TANGENTES
      MCN( MNSOFA + WNBTGS ) = NBTGS
C
C     RECOPIE DES 3 COORDONNEES DES NBSOM SOMMETS
      CALL TRTATA( MCN(MNXYZ0), MCN(MNSOFA+WYZSOM), 3*NBSOM )
      CALL TNMCDS( 'REEL', 3*NBSOM, MNXYZ0 )
C
C     RECOPIE DES 3 COMPOSANTES DES TANGENTES DES 4 LIGNES
      MN = MNSOFA+WYZSOM+3*NBSOM
      DO 8 N=1,4
         NBT = MCN( MNSOLI(N) + WNBTGS )
         IF( NBT .GT. 0 ) THEN
            NBS = MCN( MNSOLI(N) + WNBSOM )
            CALL TRTATA( MCN(MNSOLI(N)+WYZSOM+3*NBS),
     %                   MCN(MN), 3*NBT )
            MN = MN + 3 * NBT
         ENDIF
 8    CONTINUE
C
C     LA DATE DE CREATION
      CALL ECDATE( MCN(MNSOFA) )
C     LE NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNSOFA + MOTVAR(6) ) = NONMTD( '~>>>XYZSOMMET' )
C
C     TRAITEMENT DU TMS 'NSEF' DU QUADRANGLE STRUCTURE
C     ================================================
C     RECHERCHE DU NOMBRE D'EF AVEC DES TANGENTES
      CALL TNMCDC( 'ENTIER', NBEFOB, MNEFTG )
      CALL AZEROI( NBEFOB, MCN(MNEFTG) )
C
C     CALCUL DU NOMBRE D'EF A TG DU COTE 1
      MN = MNNTGL(1)
      IF( MN .GT. 0 ) THEN
         MN = MNNTGL(1)
         DO 10 I=1,NX
            IF( MCN(MN) .NE. 0 .OR. MCN(MN+1) .NE. 0 ) THEN
C              EF A TG
               MCN(MNEFTG-1+I) =  1
            ENDIF
            MN = MN + 2
 10      CONTINUE
      ENDIF
C
C     CALCUL DU NOMBRE D'EF A TG DU COTE 2
      MN = MNNTGL(2)
      IF( MN .GT. 0 ) THEN
         DO 20 J=1,NY
            IF( MCN(MN) .NE. 0 .OR. MCN(MN+1) .NE. 0 ) THEN
C              EF A TG
               I = MNEFTG - 1 + NX*J
               MCN(I) = MCN(I) + 1
            ENDIF
            MN = MN + 2
 20      CONTINUE
      ENDIF
C
C     CALCUL DU NOMBRE D'EF A TG DU COTE 3
      MN = MNNTGL(3)
      IF( MN .GT. 0 ) THEN
         MN = MNNTGL(3)
         DO 30 I=1,NX
            IF( MCN(MN) .NE. 0 .OR. MCN(MN+1) .NE. 0 ) THEN
C              EF A TG
               J = MNEFTG - 1 + NX*NY - NX + I
               MCN(J) = MCN(J) + 1
            ENDIF
            MN = MN + 2
 30      CONTINUE
      ENDIF
C
C     CALCUL DU NOMBRE D'EF A TG DU COTE 4
      MN = MNNTGL(4)
      IF( MN .GT. 0 ) THEN
         I  = MNEFTG
         DO 40 J=1,NY
            IF( MCN(MN) .NE. 0 .OR. MCN(MN+1) .NE. 0 ) THEN
C              EF A TG
               MCN(I) = MCN(I) + 1
            ENDIF
            MN = MN + 2
            I  = I + NX
 40      CONTINUE
      ENDIF
C
C     LE NOMBRE D'EF A TG
      NBEFTG = 0
      DO 50 I=MNEFTG,MNEFTG+NBEFOB-1
         IF( MCN(I) .GT. 0 ) NBEFTG = NBEFTG + 1
 50   CONTINUE
C
C     AUGMENTATION DE LA TAILLE DU TMS 'NSEF'
      NBTGEF = 8
C
C     DESTRUCTION DU TMS 'NSEF' DE LA SURFACE
      CALL LXTSDS( NTLXSU, 'NSEF' )
C
C     CONSTRUCTION DU TABLEAU 'NSEF' DE CETTE SURFACE
      CALL LXTNDC( NTLXSU, 'NSEF', 'ENTIER',
     %             WBARYQ+1+NBEFOB+NBEFTG*(1+NBTGEF) )
      CALL LXTSOU( NTLXSU, 'NSEF', NTFASU, MNFASU )
C
C     TYPE DE L'OBJET : SURFACE
      MCN ( MNFASU + WUTYOB ) = 3
C     SURFACE NON FERMEE
      MCN ( MNFASU + WUTFMA ) = 0
C     LE NOMBRE DE SOMMETS PAR FACE
      MCN ( MNFASU + WBSOEF ) = 4
C     LE NOMBRE DE TANGENTES STOCKEES PAR EF : SURFACE C1
      MCN ( MNFASU + WBTGEF ) = NBTGEF
C     LE NOMBRE D'EF DE LA SURFACE
      MCN ( MNFASU + WBEFOB ) = NBEFOB
C     LE NOMBRE D'EF AVEC TANGENTES DE LA SURFACE
      MCN ( MNFASU + WBEFTG ) = NBEFTG
C     LE NOMBRE D'EF AVEC POINTEUR SUR LES EF A TG DE LA SURFACE
      MCN ( MNFASU + WBEFAP ) = NBEFOB
C     NUMERO DU TYPE DU MAILLAGE : QUADRANGLE STRUCTURE
      MCN ( MNFASU + WUTYMA ) = 4
C     NOMBRE D'ARETES EN X DU QUADRANGLE STRUCTURE
      MCN ( MNFASU + WBARXQ ) = NX
C     NOMBRE D'ARETES EN Y DU QUADRANGLE STRUCTURE
      MCN ( MNFASU + WBARYQ ) = NY
C
C     LE POINTEUR SUR LES EF A TG DANS L'ORDRE DES EF (PAR LIGNES)
      NBT = 0
      MN  = MNFASU+WBARYQ
      MNT = MNEFTG - 1
      DO 62 J=1,NY
         DO 61 I=1,NX
            MN  = MN  + 1
            MNT = MNT + 1
            IF( MCN(MNT) .GT. 0 ) THEN
C              EF A TG DE NUMERO NBT
               NBT = NBT + 1
               MCN(MN) = NBT
            ELSE
C              EF SANS TG
               MCN(MN) = 0
            ENDIF
 61      CONTINUE
 62   CONTINUE
C     DESTRUCTION DU TABLEAU EFTG
      CALL TNMCDS( 'ENTIER', NBEFOB, MNEFTG )
      MNEFTG = MNFASU+WBARYQ+1+NBEFOB
C
C     LE CODE GEOMETRIQUE C1 degre 3 => 0
      MN = MNEFTG + NBEFTG
      DO 70 I=MNEFTG, MN-1
         MCN(I) = 0
 70   CONTINUE
C
C     LE NUMERO DES TANGENTES DES ARETES FRONTALIERES
C     -----------------------------------------------
C     LE QUADRANGLE (1,1)
      CALL AZEROI( 8, NOTG )
      I = 0
      MNT = MNNTGL(1)
      IF( MNT .GT. 0 ) THEN
         IF( MCN(MNT) .NE. 0 .OR. MCN(MNT+1) .NE. 0 ) THEN
C           EF A TG PAR LE COTE 1
            NOTG(1) = MCN(MNT)
            NOTG(4) = MCN(MNT+1)
            I = 1
         ENDIF
      ENDIF
      MNT = MNNTGL(4)
      IF( MNT .GT. 0 ) THEN
         IF( MCN(MNT) .NE. 0 .OR. MCN(MNT+1) .NE. 0 ) THEN
C           EF A TG PAR LE COTE 4
            IF( MCN(MNT) .GE. 0 ) THEN
               NOTG(2) = MCN(MNT ) + NBTG(3)
            ELSE
               NOTG(2) = MCN(MNT ) - NBTG(3)
            ENDIF
            IF( MCN(MNT+1) .GE. 0 ) THEN
               NOTG(7) = MCN(MNT+1) + NBTG(3)
            ELSE
               NOTG(7) = MCN(MNT+1) - NBTG(3)
            ENDIF
            I = 1
         ENDIF
      ENDIF
      IF( I .NE. 0 ) THEN
C        EF A TG
         CALL TRTATA( NOTG, MCN(MN), 8 )
         CALL AZEROI( 8, NOTG )
         MN = MN + 8
      ENDIF
C
C     LES QUADRANGLES (I,1)
      MNT = MNNTGL(1)
      IF( MNT .GT. 0 ) THEN
         DO 110 I=2,NX-1
            MNT = MNT + 2
            IF( MCN(MNT) .NE. 0 .OR. MCN(MNT+1) .NE. 0 ) THEN
C              EF A TG PAR LE COTE 1
               NOTG(1) = MCN(MNT)
               NOTG(4) = MCN(MNT+1)
               CALL TRTATA( NOTG, MCN(MN), 8 )
               CALL AZEROI( 8, NOTG )
               MN = MN + 8
            ENDIF
 110     CONTINUE
      ENDIF
C
C     LE QUADRANGLE (NX,1)
      I   = 0
      MNT = MNNTGL(1)
      IF( MNT .GT. 0 ) THEN
         MNT = MNT + 2 * NX - 2
         IF( MCN(MNT) .NE. 0 .OR. MCN(MNT+1) .NE. 0 ) THEN
C           EF A TG PAR LE COTE 1
            NOTG(1) = MCN(MNT)
            NOTG(4) = MCN(MNT+1)
            I = 1
         ENDIF
      ENDIF
      MNT = MNNTGL(2)
      IF( MNT .GT. 0 ) THEN
         IF( MCN(MNT) .NE. 0 .OR. MCN(MNT+1) .NE. 0 ) THEN
C           EF A TG PAR LE COTE 2
            IF( MCN(MNT) .GE. 0 ) THEN
               NOTG(3) = MCN(MNT) + NBTG(1)
            ELSE
               NOTG(3) = MCN(MNT) - NBTG(1)
            ENDIF
            IF( MCN(MNT+1) .GE. 0 ) THEN
               NOTG(6) = MCN(MNT+1) + NBTG(1)
            ELSE
               NOTG(6) = MCN(MNT+1) - NBTG(1)
            ENDIF
            I = 1
         ENDIF
      ENDIF
      IF( I .NE. 0 ) THEN
C        EF A TG
         CALL TRTATA( NOTG, MCN(MN), 8 )
         CALL AZEROI( 8, NOTG )
         MN = MN + 8
      ENDIF
C
      DO 130 J=2,NY-1
C
C        LES QUADRANGLES (1,J)
         IF( MNNTGL(4) .GT. 0 ) THEN
            MNT = MNNTGL(4) + 2 * J - 2
            IF( MCN(MNT) .NE. 0 .OR. MCN(MNT+1) .NE. 0 ) THEN
C              EF A TG PAR LE COTE 4
               IF( MCN(MNT) .GE. 0 ) THEN
                  NOTG(2) = MCN(MNT) + NBTG(3)
               ELSE
                  NOTG(2) = MCN(MNT) - NBTG(3)
               ENDIF
               IF( MCN(MNT+1) .GE. 0 ) THEN
                  NOTG(7) = MCN(MNT+1) + NBTG(3)
               ELSE
                  NOTG(7) = MCN(MNT+1) - NBTG(3)
               ENDIF
               CALL TRTATA( NOTG, MCN(MN), 8 )
               CALL AZEROI( 8, NOTG )
               MN = MN + 8
            ENDIF
         ENDIF
C
C        LES QUADRANGLES (NX,J)
         IF( MNNTGL(2) .GT. 0 ) THEN
            MNT = MNNTGL(2) + 2 * J - 2
            IF( MCN(MNT) .NE. 0 .OR. MCN(MNT+1) .NE. 0 ) THEN
C              EF A TG PAR LE COTE 2
               IF( MCN(MNT) .GE. 0 ) THEN
                  NOTG(3) = MCN(MNT) + NBTG(1)
               ELSE
                  NOTG(3) = MCN(MNT) - NBTG(1)
               ENDIF
               IF( MCN(MNT+1) .GE. 0 ) THEN
                  NOTG(6) = MCN(MNT+1) + NBTG(1)
               ELSE
                  NOTG(6) = MCN(MNT+1) - NBTG(1)
               ENDIF
               CALL TRTATA( NOTG, MCN(MN), 8 )
               CALL AZEROI( 8, NOTG )
               MN = MN + 8
            ENDIF
         ENDIF
 130  CONTINUE
C
C     LE QUADRANGLE (1,NY)
      I = 0
      MNT = MNNTGL(3)
      IF( MNT .GT. 0 ) THEN
         IF( MCN(MNT) .NE. 0 .OR. MCN(MNT+1) .NE. 0 ) THEN
C           EF A TG PAR LE COTE 3
            IF( MCN(MNT) .GE. 0 ) THEN
               NOTG(8) = MCN(MNT) + NBTG(2)
            ELSE
               NOTG(8) = MCN(MNT) - NBTG(2)
            ENDIF
            IF( MCN(MNT+1) .GE. 0 ) THEN
               NOTG(5) = MCN(MNT+1) + NBTG(2)
            ELSE
               NOTG(5) = MCN(MNT+1) - NBTG(2)
            ENDIF
            I = 1
         ENDIF
      ENDIF
      IF( MNNTGL(4) .GT. 0 ) THEN
         MNT = MNNTGL(4) + 2 * NY - 2
         IF( MCN(MNT) .NE. 0 .OR. MCN(MNT+1) .NE. 0 ) THEN
C           EF A TG PAR LE COTE 4
            IF( MCN(MNT) .GE. 0 ) THEN
               NOTG(2) = MCN(MNT) + NBTG(3)
            ELSE
               NOTG(2) = MCN(MNT) - NBTG(3)
            ENDIF
            IF( MCN(MNT+1) .GE. 0 ) THEN
               NOTG(7) = MCN(MNT+1) + NBTG(3)
            ELSE
               NOTG(7) = MCN(MNT+1) - NBTG(3)
            ENDIF
            I = 1
         ENDIF
      ENDIF
      IF( I .NE. 0 ) THEN
C        EF A TG
         CALL TRTATA( NOTG, MCN(MN), 8 )
         CALL AZEROI( 8, NOTG )
         MN = MN + 8
      ENDIF
C
C     LES QUADRANGLES (I,NY)
      MNT = MNNTGL(3)
      IF( MNT .GT. 0 ) THEN
         DO 150 I=2,NX-1
            MNT = MNT + 2
            IF( MCN(MNT) .NE. 0 .OR. MCN(MNT+1) .NE. 0 ) THEN
C              EF A TG PAR LE COTE 3
               IF( MCN(MNT) .GE. 0 ) THEN
                  NOTG(8) = MCN(MNT) + NBTG(2)
               ELSE
                  NOTG(8) = MCN(MNT) - NBTG(2)
               ENDIF
               IF( MCN(MNT+1) .GE. 0 ) THEN
                  NOTG(5) = MCN(MNT+1) + NBTG(2)
               ELSE
                  NOTG(5) = MCN(MNT+1) - NBTG(2)
               ENDIF
               CALL TRTATA( NOTG, MCN(MN), 8 )
               CALL AZEROI( 8, NOTG )
               MN = MN + 8
            ENDIF
 150     CONTINUE
      ENDIF
C
C     LE QUADRANGLE (NX,NY)
      I = 0
      MNT = MNNTGL(3)
      IF( MNT .GT. 0 ) THEN
         MNT = MNT + 2 * NX - 2
         IF( MCN(MNT) .NE. 0 .OR. MCN(MNT+1) .NE. 0 ) THEN
C           EF A TG PAR LE COTE 1
            IF( MCN(MNT) .GE. 0 ) THEN
               NOTG(8) = MCN(MNT) + NBTG(2)
            ELSE
               NOTG(8) = MCN(MNT) - NBTG(2)
            ENDIF
            IF( MCN(MNT+1) .GE. 0 ) THEN
               NOTG(5) = MCN(MNT+1) + NBTG(2)
            ELSE
               NOTG(5) = MCN(MNT+1) - NBTG(2)
            ENDIF
            I = 1
         ENDIF
      ENDIF
      MNT = MNNTGL(2)
      IF( MNT .GT. 0 ) THEN
         MNT = MNT + 2 * NY - 2
         IF( MCN(MNT) .NE. 0 .OR. MCN(MNT+1) .NE. 0 ) THEN
C           EF A TG PAR LE COTE 2
            IF( MCN(MNT) .GE. 0 ) THEN
               NOTG(3) = MCN(MNT) + NBTG(1)
            ELSE
               NOTG(3) = MCN(MNT) - NBTG(1)
            ENDIF
            IF( MCN(MNT+1) .GE. 0 ) THEN
               NOTG(6) = MCN(MNT+1) + NBTG(1)
            ELSE
               NOTG(6) = MCN(MNT+1) - NBTG(1)
            ENDIF
            I = 1
         ENDIF
      ENDIF
      IF( I .NE. 0 ) THEN
C        EF A TG
         CALL TRTATA( NOTG, MCN(MN), 8 )
         CALL AZEROI( 8, NOTG )
         MN = MN + 8
      ENDIF
      END
