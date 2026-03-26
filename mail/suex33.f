        SUBROUTINE SUEX33( NTLXSU, LADEFI,
     %                     NTFASU, MNFASU, NTSOFA, MNSOFA, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    GENERER LA TRIANGULATION-QUADRANGULATION D'UNE SURFACE
C -----    PAR JONCTION DES ARETES DE LIGNES 2 A 2
C
C          ATTENTION: LES SOMMETS DES ARETES DOIVENT APPARAITRE
C                     DANS LE MEME SENS SOUS RISQUE D'OBTENIR
C                     DES EF VRILLES, DEGENERES, ...
C ENTREES:
C --------
C NTLXSU : NUMERO DU TABLEAU TS DU LEXIQUE DE LA SURFACE A CREER
C LADEFI : TABLEAU ENTIER DE DEFINITION DE LA SURFACE
C          CF '~/td/d/a_surface__definition'

C SORTIES:
C --------
C NTFASU : NUMERO      DU TMS 'NSEF' DES NUMEROS DES FACES
C MNFASU : ADRESSE MCN DU TMS 'NSEF' DES NUMEROS DES FACES
C          CF '~/td/d/a___nsef'
C NTSOFA : NUMERO      DU TMS 'XYZSOMMET' DE LA SURFACE
C MNSOFA : ADRESSE MCN DU TMS 'XYZSOMMET' DE LA SURFACE
C          CF '~/td/d/a___xyzsommet'
C IERR   : = 0 SI PAS D'ERREUR
C          <>0 SI ERREUR
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS   NOVEMBRE 1996
C23456...............................................................012
      include"./incl/gsmenu.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/a_surface__definition.inc"
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
      INTEGER           LADEFI(0:*)
      CHARACTER*24      KNOMLG

      IERR   = 0
      MNSOLG = 0

C     RESTAURATION DU MAILLAGE DES LIGNES
C     -----------------------------------
      NBLGJD = LADEFI( WBLGJD )
      IF(NBLGJD .LE. 1 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'NOMBRE INCORRECT <2 DE LIGNES A JOINDRE'
         CALL LEREUR
         IERR = 2
         GOTO 9900
      ENDIF

C     LE TABLEAU DES ADRESSES DES TABLEAUX 'XYZSOMMET' ET 'NSEF'
C     DES LIGNES
      CALL TNMCDC( 'ENTIER', NBLGJD+NBLGJD, MNSOLG )
      MNSSLG = MNSOLG + NBLGJD

C     LE NOMBRE D'ARETES DE LA LIGNE QUI PRECEDE
      NBARAV = 0
C     LE NOMBRE DE SOMMETS ET FACES DU MAILLAGE DE LA SURFACE
      NBSOSU = 0
      NBTGSU = 0
      NBTG1  = 0
      NBFASU = 0
      NBEFTG = 0
C     LE TYPE DU MAILLAGE DE LA SURFACE A PRIORI QUADRANGLE STRUCTURE
      NUTYMS = 4

C     BOUCLE SUR LES LIGNES DE LA FUTURE SURFACE
      NLGPRE = 0
      DO 10 N = 1, NBLGJD

C        OUVERTURE DE LA LIGNE N
         NLG = LADEFI( WULGJD-1+N )
         CALL LXNLOU( NTLIGN, NLG, NT, MN )
         IF( NT .LE. 0 ) THEN
            NBLGRC(NRERR) = 1
            WRITE(KERR(MXLGER)(1:5),'(I5)' ) N
            KERR(1) = 'LIGNE '// KERR(MXLGER)(1:5) //' INCONNUE'
            CALL LEREUR
            IERR = IERR + 1000
            GOTO 10
         ENDIF
         CALL NMOBNU( 'LIGNE', NLG, KNOMLG )
         IF( NLGPRE .EQ. NLG ) THEN
            NBLGRC(NRERR) = 2
            KERR(1) = KNOMLG
            KERR(2) = 'IDENTIQUE A LA PRECEDENTE'
            CALL LEREUR
            IERR = IERR + 100
            GOTO 10
         ENDIF

         NLGPRE = NLG
C        RESTAURATION DES TABLEAUX SOMMETS ET NSEF DE LA LIGNE NLG
         CALL LXTSOU( NT, 'XYZSOMMET', NTSOLI, MNS )
         IF( NTSOLI .LE. 0 ) THEN
            NBLGRC(NRERR) = 2
            KERR(1) = KNOMLG
            KERR(2) = 'LIGNE SANS TMS XYZSOMMET'
            CALL LEREUR
            IERR = IERR + 100
            GOTO 10
         ENDIF

         CALL LXTSOU( NT, 'NSEF', NTARLI, MNA )
         IF( NTARLI .LE. 0 ) THEN
            NBLGRC(NRERR) = 2
            KERR(1) = KNOMLG
            KERR(2) = 'LIGNE SANS TMS NSEF'
            CALL LEREUR
            IERR = IERR + 10
            GOTO 10
         ENDIF

C        TENTATIVE DE RESTRUCTURATION DE LA LIGNE N
         CALL LIGSTR( NT, NTARLI, MNA, NTSOLI, MNS, IER )
         IF( IER .GT. 2 ) THEN
            IERR = IERR + IER
            GOTO 10
         ENDIF

C        STOCKAGE DE L'ADRESSE DU TMS XYZSOMMET DE LA LIGNE N
         MCN(MNSOLG-1+N) = MNS

C        LE NOMBRE DE SES SOMMETS
         NBSOM1 = MCN(MNS+WNBSOM)

C        LA LIGNE EXISTE
         MCN(MNSSLG-1+N) = MNA

C        LA LIGNE SUIVANTE DOIT AVOIR UNE DIFFERENCE D'AU PLUS 1 ARETE
         IF( IERR .GT. 0 ) GOTO 10
         IF( MCN(MNA+WUTYMA) .EQ. 0 ) THEN
C           LIGNE NON STRUCTURE
            NBARE1 = MCN( MNA + WBEFOB )
         ELSE
C           LIGNE STRUCTUREE
            NBARE1 = MCN( MNA + WBARSE )
         ENDIF
         IF( NBSOM1 .EQ. NBARE1 ) THEN
C           LIGNE FERMEE => LE MAILLAGE DE LA SURFACE NE PEUT ETRE STRUCTURE
            NUTYMS = 0
         ENDIF

C        LE NOMBRE DE SOMMETS DE LA SURFACE MAILLEE
         NBSOSU = NBSOSU + MCN(MNS+WNBSOM)
C
C        LE NOMBRE DE TANGENTES DE CETTE LIGNE ET DE SA PRECEDENTE
         NBTG0 = NBTG1
         NBTG1 = MCN( MNS + WNBTGS )
C        LA SOMME DU NOMBRE DES TANGENTES CETTE LIGNE N COMPRISE
C        C'EST LE NOMBRE DE TANGENTES DE LA SURFACE
         NBTGSU = NBTGSU + NBTG1

         IF( NBARAV .GT. 0 ) THEN

C           TRAITEMENT DE LA COUCHE N-1 ENTRE LES LIGNES N-1 ET N POUR N>1
            IF( ABS( NBARAV - NBARE1 ) .GT. 1 ) THEN
               NBLGRC(NRERR) = 3
               KERR(1) = KNOMLG
               KERR(2) = 'TROP GRANDE DIFFERENCE ENTRE NOMBRE D''ARETES'
               KERR(3) = 'DE LA LIGNE ET CELUI DE LA LIGNE PRECEDENTE'
               CALL LEREUR
               IERR = IERR + 10000
            ENDIF

C           LE MAILLAGE RESULTANT EST-IL STRUCTURE ?
C           C-A-D TOUTES LES LIGNES ONT MEME NOMBRE DE LIGNES
C           ET AUCUNE N'EST FERMEE
            IF( NBARAV .NE. NBARE1 ) THEN
C              LE MAILLAGE EST NON STRUCTURE
               NUTYMS = 0
C              LE NOMBRE DE TRIANGLES DE CETTE COUCHE
               NBFASU = NBFASU + NBARAV + NBARE1
C              LE NOMBRE D'EF A TG
               NBEFTG = NBEFTG + NBTG0 + NBTG1
            ELSE
C              LE NOMBRE DE QUADRANGLES DE CETTE COUCHE
               NBFASU = NBFASU + NBARE1
C              MAJORATION DU NOMBRE D'EF A TG
C              IL PEUT ARRIVER LE CAS D'UNE TG SUR 1 LIGNE
C              ET 2 AUTRES SUR DES ARETES DIFFERENTES DE l'AUTRE LIGNE!
               NBEFTG = NBEFTG + MIN( NBARAV, NBTG0+NBTG1 )
            ENDIF

         ENDIF

C        LE NOMBRE DES ARETES DE LA LIGNE D'AVANT
         NBARAV = NBARE1

 10   ENDDO
      IF( IERR .NE. 0 ) GOTO 9900

C     CONSTRUCTION DU TABLEAU 'XYZSOMMET' DE CETTE SURFACE
C     ----------------------------------------------------
      CALL LXTNDC( NTLXSU, 'XYZSOMMET', 'MOTS',
     %             WYZSOM + 3 * ( NBSOSU + NBTGSU ) )
      CALL LXTSOU( NTLXSU, 'XYZSOMMET',  NTSOFA, MNSOFA )

C     CONSTRUCTION DU TABLEAU 'NSEF' DE CETTE SURFACE
C     -----------------------------------------------
C     DECLARATION DES NO SOMMET DE LA SURFACE
      IF( NUTYMS .EQ. 0 ) THEN
         MOT = WUSOEF + 4 * NBFASU
      ELSE
         MOT = WBARYQ + 1
      ENDIF
      IF( NBTGSU .GT. 0 ) THEN
         NBTGEF = 8
         NBEFAP = NBFASU
         NBEFTG = NBFASU
      ELSE
         NBTGEF = 0
         NBEFAP = 0
         NBEFTG = 0
      ENDIF
      CALL LXTNDC( NTLXSU, 'NSEF', 'ENTIER',
     %             MOT + NBEFAP + NBEFTG * (1+NBTGEF) )
      CALL LXTSOU( NTLXSU, 'NSEF',  NTFASU, MNFASU )
      MNEFAP = MNFASU + MOT
      MNNUTG = MNEFAP + NBEFAP + NBEFTG
      NBEFT0 = NBEFTG

C     BOUCLE SUR LES COUCHES A MAILLER ENTRE 2 LIGNES SUCCESSIVES
C     -----------------------------------------------------------
      NUDEST = 0
      NUDEEF = 0
      NUSENS = 1
      DO N=1,NBLGJD-1

C        LES ADRESSES DES TABLEAUX 'XYZSOMMET' ET NSEF'
         MNS1   = MCN( MNSOLG + N - 1 )
         NBSOM1 = MCN(MNS1+WNBSOM)
         MNA1   = MCN( MNSSLG + N - 1 )
         IF( MCN(MNA1+WUTYMA) .EQ. 0 ) THEN
             MN1    = MNA1 + WBEFOB
             NBARE1 = MCN(MN1)
             MN1    = MNA1 + WUSOEF + 2 * NBARE1
         ELSE
             MN1    = MNA1 + WBARSE
             NBARE1 = MCN(MN1)
             MN1    = MN1 + 1
         ENDIF
         NBEFA1 = MCN( MNA1 + WBEFAP )
         NBEFT1 = MCN( MNA1 + WBEFTG )
         IF( NBEFT1 .GT. 0 ) THEN
            MNT1 = MN1 + NBEFA1 + NBEFT1
         ELSE
            MNT1 = MN1
         ENDIF
         MNS2   = MCN( MNSOLG + N )
         NBSOM2 = MCN(MNS2+WNBSOM)
         MNA2   = MCN( MNSSLG + N )
         IF( MCN(MNA2+WUTYMA) .EQ. 0 ) THEN
             MN2    = MNA2 + WBEFOB
             NBARE2 = MCN(MN2)
             MN2    = MNA2 + WUSOEF + 2 * NBARE2
         ELSE
             MN2    = MNA2 + WBARSE
             NBARE2 = MCN(MN2)
             MN2    = MN2 + 1
         ENDIF
         NBEFA2 = MCN( MNA2 + WBEFAP )
         NBEFT2 = MCN( MNA2 + WBEFTG )
         IF( NBEFT2 .GT. 0 ) THEN
            MNT2 = MN2 + NBEFA2 + NBEFT2
         ELSE
            MNT2 = MN2
         ENDIF

C        GENERATION DES ELEMENTS FINIS DE CETTE COUCHE
         CALL JOTQ2L( NBSOM1, RMCN(MNS1+WYZSOM),
     %                MCN(MNS1+WNBTGS), RMCN(MNS1+WYZSOM+3*NBSOM1),
     %                MCN(MNA1+WUTYMA), MCN(MNA1+WUTFMA),
     %                NBARE1, MCN(MNA1+WUSOEF),
     %                MCN(MN1), NBEFT1, MCN(MNT1),
     %                MCN(MNS2+WNBSOM), RMCN(MNS2+WYZSOM),
     %                MCN(MNS2+WNBTGS), RMCN(MNS2+WYZSOM+3*NBSOM2),
     %                MCN(MNA2+WUTYMA),  MCN(MNA2+WUTFMA),
     %                NBARE2, MCN(MNA2+WUSOEF),
     %                MCN(MN2), NBEFT2, MCN(MNT2),
     %                NUDEST, NUSENS, NBSOSU, MCN(MNSOFA+WYZSOM),
     %                NUDETG, NBTGSU, MCN(MNSOFA+WYZSOM+3*NBSOSU),
     %                NUDEEF, NUTYMS, NBFASU, MCN(MNFASU+WUSOEF),
     %                MCN(MNEFAP), NBEFTG, MCN(MNNUTG) )

      ENDDO

C     MISE A JOUR DU TABLEAU 'XYZSOMMET'
C     ----------------------------------
C     LE NOMBRE DE SOMMETS
      MCN( MNSOFA + WNBSOM) = NBSOSU
C     LE NOMBRE DE TANGENTES
      MCN( MNSOFA + WNBTGS ) = NBTGSU
C     LE NOMBRE DE COORDONNEES PAR SOMMET
      MCN( MNSOFA + WBCOOR ) = 3
C     LA DATE DE CREATION
      CALL ECDATE( MCN(MNSOFA) )
C     LE NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNSOFA + MOTVAR(6) ) = NONMTD( '~>>>XYZSOMMET' )

C     MISE A JOUR DU TABLEAU 'NSEF'
C     -----------------------------
      IF( NBEFTG .GT. 0 ) THEN
C        LE CODE GEOMETRIQUE
         DO N=MNEFAP+NBEFAP,MNEFAP+NBEFAP+NBEFTG-1
            MCN(N) = 0
         ENDDO
         IF( NBEFT0 .GT. NBEFTG ) THEN
C           TRANSLATION DU TABLEAU DES NUMEROS DES TANGENTES
            MNEFAP = MNEFAP + NBEFAP + NBEFTG
            MN2    = NBEFT0 - NBEFTG
            DO N=MNEFAP,MNEFAP+NBTGEF*NBEFTG-1
               MCN(N) = MCN(N+MN2)
            ENDDO
            CALL TAMSRA( NTFASU, MOT + NBEFAP + NBEFTG * (1+NBTGEF) )
         ENDIF
      ENDIF

C     TYPE DE L'OBJET : SURFACE TRIANGULEE NON STRUCTUREE
      MCN( MNFASU + WUTYOB ) = 3
      MCN( MNFASU + WUTYMA ) = NUTYMS
      IF( NUTYMS .EQ. 4 ) THEN
         MCN( MNFASU + WBARXQ ) = NBARE1
         MCN( MNFASU + WBARYQ ) = NBLGJD - 1
      ENDIF
C     NOMBRE DE SOMMETS PAR EF
      MCN( MNFASU + WBSOEF ) = 4
C     NOMBRE D'EF DU MAILLAGE
      MCN( MNFASU + WBEFOB ) = NBFASU
C     LE TYPE INCONNU DE FERMETURE DU MAILLAGE
      MCN( MNFASU + WUTFMA ) = -1
C     LES TANGENTES STOCKEES
      MCN( MNFASU + WBTGEF ) = NBTGEF
      MCN( MNFASU + WBEFAP ) = NBEFAP
      MCN( MNFASU + WBEFTG ) = NBEFTG
C     LA DATE DE CREATION
      CALL ECDATE( MCN(MNFASU) )
C     LE NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNFASU + MOTVAR(6) ) = NONMTD( '~>>>NSEF' )

 9900 IF( MNSOLG.GT.0 ) CALL TNMCDS( 'ENTIER', NBLGJD+NBLGJD, MNSOLG )

      RETURN
      END
