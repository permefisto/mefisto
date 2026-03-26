      SUBROUTINE QDEGET( NUTYOB, NMOBJT, MNTSMA, MNSOMM, NBEFMQ, IERR  )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    ESSAYER DE TRANSFORMER LES QUADRANGLES DEGENERES DU MAILLAGE
C------    EN TRIANGLES LORSQUE 2 SOMMETS CONSECUTIFS PEUVENT ETRE
C          IDENTIFIES
C
C ENTREES:
C --------
C NUTYOB : NUMERO DU TYPE DU MAILLAGE A TRACER ( 1:SOMMETS 2:ARETES ...)
C NMOBJT : NOM DE L'OBJET
C
C MODIFIES :
C ----------
C MNTSMA : ADRESSE MCN DU TABLEAU MS 'NSEF' DU MAILLAGE
C MNSOMM : ADRESSE MCN DU TABLEAU MS 'XYZSOMMET' DES SOMMETS
C NBEFMQ : NOMBRE D'EF DE QUALITE INFERIEURE A LA QUALITE MINIMALE
C IERR   : 1 SI TRIANGLE STRUCTURE AVEC UN SOUS-TRIANGLE=ARETE
C          2 SI QUADRANGLE EF AVEC 3 SOMMETS IDENTIQUES
C          INCHANGE SINON
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS       AVRIL 1993
C....................................................................012
      PARAMETER        (QUEFMI=0.001)
      IMPLICIT          INTEGER(W)
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU , IMPRIM , INTERA , NUNITE(29)
      COMMON / MSSFTA /  MSSF(28),NTADAM
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
      CHARACTER*(*)     NMOBJT
      INTEGER           NOSOEL(64)
C
C     SI PAS D'EF DE MAUVAISE QUALITE RETOUR
      IF( NBEFMQ .LE. 0 ) RETURN
C
C     POUR L'INSTANT SEULES LES SURFACES SONT TRAITEES
C     (EN FAIT SEULEMENT LES QUADRANGLES )
      IF( NUTYOB .NE. 3 ) RETURN
C
C     ADRESSE-3 DE LA 1-ERE COORDONNEE DU 1-ER SOMMET DU TMS 'XYZSOMMET'
      MNS   = MNSOMM + WYZSOM
      NBSOM = MCN( MNSOMM + WNBSOM )
C
C     LES PARAMETRES DES NO SOMMET DU MAILLAGE
      CALL NSEFPA( MCN(MNTSMA) ,
     %             NUTYMA , NBSOEL , NBSOEF , NBTGEF,
     %             LDAPEF , LDNGEF , LDTGEF , NBEFOB ,
     %             NX     , NY     , NZ     ,
     %             IERR   )
      IF( IERR .NE. 0 ) RETURN
C
C     LA BOUCLE SUR LES TRIANGLES ET / OU QUADRANGLES DE LA SURFACE
C     =============================================================
      MNNEW  = 0
      MNPILE = 0
      MXPILE = 0
      LHPILE = 0
      MNST   = MNS - 3
      DO 100 N=1,NBEFOB
C
C        LE NUMERO DES NBSOEF SOMMETS DU SOUSOBJET N
         CALL NSEFNS( N      , NUTYMA , NBSOEF , NBTGEF,
     %                LDAPEF , LDNGEF , LDTGEF ,
     %                MNTSMA , NX , NY , NZ ,
     %                NCOGEL , NUGEEF , NUEFTG, NOSOEL , IERR )
C
C        LA QUALITE DE L'EF
         CALL QUALEF( NCOGEL ,  NOSOEL , NBSOM , RMCN(MNS) ,
     %                SURFVOLU, QUALIT , IERR )
C
         IF( QUALIT .LT. QUEFMI .AND. NCOGEL .EQ. 4 ) THEN
C           QUADRANGLE DEGENERE
C           2 SOMMETS CONSECUTIFS SONT ILS IDENTIQUES?
            NBSI = 0
            DO 10 I=1,4
               IF( I .EQ. 4 ) THEN
                  I1 = 1
               ELSE
                  I1 = I + 1
               ENDIF
               MNS1 = MNST + 3 * NOSOEL(I)
               MNS2 = MNST + 3 * NOSOEL(I1)
               CALL XYZIDE( RMCN(MNS1), RMCN(MNS2), ID )
               IF( ID .NE. 0 ) THEN
C                 LES 2 POINTS SONT IDENTIQUES: ILS SONT EMPILES
                  IF( NBSI .NE. 0 ) THEN
C                    AU MOINS 2 POINTS IDENTIFIES A LEUR SUIVANT
C                    RIEN NE PEUT ETRE FAIT
                     NBLGRC(NRERR) = 2
                     WRITE(KERR(MXLGER)(1:10),'(I10)') N
                     KERR(1) ='QUADRANGLE'//KERR(MXLGER)(1:10)
                     KERR(2) ='AVEC 3 SOMMETS IDENTIFIES'
                     CALL LEREUR
                     IERR = 2
                     GOTO 100
                  ENDIF
                  NBSI = 1
                  IF( MXPILE .LE. LHPILE ) THEN
C                    PILE A AUGMENTER
                     CALL TNMCAU( 'ENTIER', MXPILE*2, (MXPILE+128)*2,
     %                             LHPILE*2, MNPILE )
                     MXPILE = MXPILE + 128
                  ENDIF
                  IF( NOSOEL(I) .GT. NOSOEL(I1) ) THEN
                     NS1 = NOSOEL(I1)
                     NS2 = NOSOEL(I)
                  ELSE
                     NS1 = NOSOEL(I)
                     NS2 = NOSOEL(I1)
                  ENDIF
C                 STOCKAGE LE PLUS PETIT D'ABORD
                  MCN( MNPILE + 2 * LHPILE     ) = NS1
                  MCN( MNPILE + 2 * LHPILE + 1 ) = NS2
                  LHPILE = LHPILE + 1
               ENDIF
 10         CONTINUE
         ENDIF
 100  CONTINUE
      IF( LHPILE .LE. 0 ) RETURN
C
C     SUPPRESSION DES DOUBLETS
      MNP  = MNPILE - 2
      MNP1 = MNP
      DO 150 I=2,LHPILE
         MNP1 = MNP1 + 2
C        LES COORDONNEES DU SOMMET 1
         NS1  = MCN( MNP1 )
         MNS1 = MNST + 3 * NS1
         MNP2 = MNP
         DO 110 J=1,I-1
            MNP2 = MNP2 + 2
            NS2  = MCN( MNP2 )
            MNS2 = MNST + 3 * NS2
            CALL XYZIDE( RMCN(MNS1), RMCN(MNS2), ID )
            IF( ID .NE. 0 ) THEN
C              LES 2 POINTS SONT IDENTIQUES
               IF( NS1 .GT. NS2 ) THEN
                  N = NS2
               ELSE
                  N = NS1
               ENDIF
C              LES 2 SOMMETS ONT LE PLUS PETIT NUMERO
               MCN( MNP1 ) = N
               MCN( MNP2 ) = N
               GOTO 150
            ENDIF
 110     CONTINUE
 150  CONTINUE
C
C     GENERATION DU NOUVEAU NUMERO DE SOMMET
      CALL TNMCDC( 'ENTIER', NBSOM+1, MNNEW )
      DO 170 I=0,NBSOM
         MCN( MNNEW + I ) = I
 170  CONTINUE
C
C     IDENTIFICATION DES SOMMETS
      DO 180 I=1,LHPILE
         NS1 = MCN( MNP + 2 * I )
         NS2 = MCN( MNP + 2 * I + 1 )
         MCN( MNNEW + NS2 ) = NS1
 180  CONTINUE
C
C     LES NOUVEAUX NUMEROS DES SOMMETS
      NBS = 0
      DO 200 I=1,NBSOM
         NS1 = MCN( MNNEW + I )
         IF( NS1 .EQ. I ) THEN
C           SOMMET NON IDENTIFIE
            NBS = NBS + 1
            MCN( MNNEW + I ) = NBS
         ELSE
C           SOMMET IDENTIFIE
            MCN( MNNEW + I ) = MCN( MNNEW + NS1 )
         ENDIF
 200  CONTINUE
C
C     LE NUMERO DE LEXIQUE DE LA QUADRANGULATION
      CALL LXLXOU( NTADAM, 'SURFACE', NTSURF, N )
      CALL LXLXOU( NTSURF, NMOBJT, NTLXSU, MNLXSU )
C
C     SI LE MAILLAGE EST UNE QUADRANGULATION STRUCTUREE
C     ALORS ELLE EST DESTUCTUREE
      IF( NUTYMA .EQ. 4 ) THEN
C        LE NOMBRE D'INTERVALLES EN X ET Y
         NBARXQ = MCN ( MNTSMA + WBARXQ )
         NBARYQ = MCN ( MNTSMA + WBARYQ )
C        DESTRUCTION DE NSEF DE L'ANCIENNE QUADRANGULATION
         CALL LXTSDS( NTLXSU, 'NSEF' )
C        CONSTRUCTION DU TABLEAU 'NSEF'
         CALL LXTNDC( NTLXSU , 'NSEF' , 'ENTIER' ,
     %                WUSOEF + NBEFOB * 4 )
         CALL LXTSOU( NTLXSU, 'NSEF', NTTSMA, MNTSMA )
C        LE TYPE DE L'OBJET : ICI SURFACE
         MCN( MNTSMA + WUTYOB ) = 3
C        LE TYPE DU MAILLAGE NUTYMA : ICI NON STRUCTURE
         MCN( MNTSMA + WUTYMA ) = 0
C        LE NOMBRE DE SOMMETS PAR SOUS-OBJET
         MCN( MNTSMA + WBSOEF ) = 4
C        LE NOMBRE DE NSEF
         MCN( MNTSMA + WBEFOB ) = NBEFOB
         MN  = MNTSMA + WUSOEF - 1
         NBX = NBARXQ + 1
         DO 220 NY=1,NBARYQ
            DO 210 NX=1,NBARXQ
               MCN( MN + 1 ) = (NY-1)*NBX+NX
               MCN( MN + 2 ) = (NY-1)*NBX+NX+1
               MCN( MN + 3 ) = NY*NBX+NX+1
               MCN( MN + 4 ) = NY*NBX+NX
               MN            = MN + 4
 210        CONTINUE
 220     CONTINUE
C        LE TYPE INCONNU DE FERMETURE DU MAILLAGE
         MCN( MNTSMA + WUTFMA ) = -1
C        PAS DE TANGENTES STOCKEES
         MCN( MNTSMA + WBTGEF ) = 0
         MCN( MNTSMA + WBEFAP ) = 0
         MCN( MNTSMA + WBEFTG ) = 0
C        AJOUT DE LA DATE
         CALL ECDATE( MCN(MNTSMA) )
C        AJOUT DU NUMERO DU TABLEAU DESCRIPTEUR
         MCN( MNTSMA + MOTVAR(6) ) = NONMTD( '~>>>NSEF' )
      ENDIF
C
C     MODIFICATION DU NUMERO DES SOMMETS DES EF
      MN0    = MNTSMA + WUSOEF - 1
      MN     = MN0
      NBEFO1 = NBEFOB
      DO 300 N=1,NBEFOB
C
C        MISE A JOUR AVEC LE NOUVEAU NUMERO DE SOMMET
         DO 250 I=1,4
C           L'ANCIEN NUMERO DU SOMMET I DE L'EF N
            NS1 = MCN( MN0 + I )
C           LE NOUVEAU NUMERO
            MCN( MN + I ) = MCN( MNNEW + NS1 )
 250     CONTINUE
C
C        LE QUADRANGLE A T IL 2 SOMMETS DE NUMEROS IDENTIQUES?
         DO 280 I=1,4
C           LE NUMERO DE 2 SOMMETS CONSECUTIFS DU QUADRANGLE
            NS1 = MCN( MN + I )
            IF( I .EQ. 4 ) THEN
               I1 = 1
            ELSE
               I1 = I + 1
            ENDIF
            NS2 = MCN( MN + I1 )
            IF( NS1 .EQ. NS2 ) THEN
C
C              2 SOMMETS AVEC LE MEME NUMERO DANS CET EF
C
C              EST CE BIEN UN QUADRANGLE ?
               NS4 = MCN( MN + 4 )
               IF( NS4 .LE. 0 ) THEN
C                 C'EST UN TRIANGLE => IL EST SUPPRIME
                  IF( NUTYMA .EQ. 3 ) THEN
C                    TRIANGLE STRUCTURE => IMPOSSIBLE A TRAITER
                     NBLGRC(NRERR) = 3
                     WRITE(KERR(MXLGER)(1:10),'(I10)') N
                     KERR(1) ='TRIANGLE'//KERR(MXLGER)(1:10)
                     KERR(2) ='AVEC 2 SOMMETS IDENTIFIES'
                     KERR(3) ='DANS UN TRIANGLE STRUCTURE'
                     CALL LEREUR
                     IERR = 1
                  ENDIF
C                 CAS NON STRUCTURE : CE TRIANGLE EST SUPPRIME
C                                     SANS DECLENCHER UNE ERREUR
                  NBEFO1 = NBEFO1 - 1
                  GOTO 290
               ENDIF
C
C              C'EST UN QUADRANGLE => IL DEVIENT UN TRIANGLE
               IF( I .EQ. 1 ) THEN
                  MCN( MN + 2 ) = NS4
               ELSE IF( I .EQ. 2 ) THEN
                  MCN( MN + 3 ) = NS4
               ENDIF
C              C'EST UN TRIANGLE
               MCN( MN + 4 ) = 0
            ENDIF
 280     CONTINUE
C
C        PASSAGE A L'EF SUIVANT
         MN  = MN  + 4
 290     MN0 = MN0 + 4
 300  CONTINUE
C
C     LE NOMBRE DE EF EST REMIS A JOUR
      MCN( MNTSMA + WBEFOB ) = NBEFO1
C
C     DECALAGE DES COORDONNEES DES SOMMETS
      MNS = MNSOMM + WYZSOM - 3
      DO 400 N=1,NBSOM
C        LE NOUVEAU NUMERO
         NS1 = MCN( MNNEW + N )
         IF( NS1 .NE. N ) THEN
C           LES 3 NOUVELLES COORDONNEES
            MNS1 = MNS + 3 * NS1
C           LES 3 ANCIENNES COORDONNEES
            MNS2 = MNS + 3 * N
            RMCN( MNS1   ) = RMCN( MNS2 )
            RMCN( MNS1+1 ) = RMCN( MNS2+1 )
            RMCN( MNS1+2 ) = RMCN( MNS2+2 )
         ENDIF
 400  CONTINUE
C
C     MODIFICATION DU NOMBRE DE SOMMETS
      MCN( MNSOMM + WNBSOM ) = NBS
C
C     MODIFICATION DE LA DATE
      CALL ECDATE( MCN(MNSOMM) )
C
C     DESTRUCTION DES TABLEAUX INUTILES
      IF( MNNEW  .GT. 0 ) CALL TNMCDS( 'ENTIER', NBSOM+1,  MNNEW )
      IF( MNPILE .GT. 0 ) CALL TNMCDS( 'ENTIER', MXPILE*2, MNPILE )
      END
