      SUBROUTINE VOEX04( NTLXVO, LADEFI,
     %                   NTVOOB, MNVOOB, NTSOVO, MNSOVO, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    GENERER LE MAILLAGE D'UN PENTAEDRE DEFINI PAR LES
C -----    SURFACES STRUCTUREES DE CHACUNE DE SES 5 FACES PAR
C          INTERPOLATION TRANSFINIE DE DEGRE 1 (cf CRAS A. PERRONNET)
C ENTREES:
C --------
C NTLXVO : NUMERO DU TABLEAU TS DU LEXIQUE DU VOLUME
C LADEFI : TABLEAU DE DEFINITION DE LE VOLUME PARTITIONNEE
C          cf '~/td/d/a_volume__definition'
C
C SORTIES:
C --------
C NTVOOB : NUMERO      DU TMS 'NSEF' DES NUMEROS DES VOLUMES
C MNVOOB : ADRESSE MCN DU TMS 'NSEF' DES NUMEROS DES VOLUMES
C          cf '~/td/d/a___nsef'
C NTSOVO : NUMERO      DU TMS 'XYZSOMMET' DE L'OBJET
C MNSOVO : ADRESSE MCN DU TMS 'XYZSOMMET' DE L'OBJET
C          cf '~/td/d/a___xyzsommet'
C IERR   : 0 SI PAS D'ERREUR
C        > 0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY      ANALYSE NUMERIQUE UPMC PARIS   MARS     1989
C MODIFS : ALAIN PERRONNET  ANALYSE NUMERIQUE UPMC PARIS   NOVEMBRE 1997
C23456---------------------------------------------------------------012
      IMPLICIT INTEGER (W)
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/ntmnlt.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/a_surface__definition.inc"
      include"./incl/a_volume__definition.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE       (MCN(1),RMCN(1))
      REAL              XYZ(3,4,5)
      CHARACTER*24      KNOMSU
      INTEGER           LADEFI(0:*)
      INTEGER           NBSOFA(5,4),MNSOFA(5),NUSOFA(4),
     %                  NUFACE(5),NSENS(5),NPSOM(5),NUTYPE(5),
     %                  QSUIV(4)
      DATA              QSUIV/2,3,4,1/
C
      IERR = 0
C
C     VERIFICATION DU TYPE DE CHAQUE FACE
C     ===================================
      DO 100 NF=1,5
C
C        LE NUMERO DE LA FACE
         NOFA   = LADEFI(WU5FAC-1+NF)
         IF( NOFA .LE. 0 ) THEN
            NBLGRC(NRERR) = 1
            WRITE(KERR(MXLGER)(1:10),'(I10)') NOFA
            KERR(1) = 'SURFACE DE NUMERO' // KERR(MXLGER)(1:10)
     %             // ' INCORRECT'
            CALL LEREUR
            RETURN
         ENDIF
C        LE TABLEAU LEXIQUE DE CETTE FACE
         CALL LXNLOU( NTSURF, NOFA, NTLXSU, MN )
         IF( NTLXSU .LE. 0 ) THEN
            NBLGRC(NRERR) = 1
            WRITE(KERR(MXLGER)(1:10),'(I10)') NF
            KERR(1) = 'SURFACE INCONNUE SUR LA FACE '
     %             // KERR(MXLGER)(1:10)
            CALL LEREUR
            RETURN
         ENDIF
C
C        LE NOM DE LA SURFACE
         CALL NMOBNU( 'SURFACE', NOFA, KNOMSU )
C
C        LE TABLEAU 'NSEF' DE CETTE FACE
         CALL LXTSOU( NTLXSU, 'NSEF', NTFASU, MNFASU )
         IF( NTFASU .LE. 0 ) THEN
            NBLGRC(NRERR) = 1
            KERR(1) = 'SURFACE SANS NSEF : ' // KNOMSU
            CALL LEREUR
            IERR = 5
            RETURN
         ENDIF
C
C        LA FACE EST ELLE STRUCTUREE ?
         NUTYMA = MCN( MNFASU + WUTYMA )
         NUTYPE(NF) = NUTYMA
         IF  ( NUTYMA .NE. 3 ) THEN
            IF ( NUTYMA .NE. 4 ) THEN
               NBLGRC(NRERR) = 2
               KERR(1) = 'SURFACE NON STRUCTUREE : ' // KNOMSU
               WRITE(KERR(MXLGER)(1:10),'(I10)') NOFA
               KERR(2) =  'FACE' // KERR(MXLGER)(1:10)
     %                    // ' TYPE DIFFERENT DE 3 OU 4'
               CALL LEREUR
               IERR = 1
            ENDIF
         ENDIF
C        LE NOMBRE DE SOMMETS DANS CHAQUE DIRECTION
         IF (NUTYMA .EQ. 3) THEN
            NBSOFA(NF,1) = MCN ( MNFASU + WBARTR ) + 1
            NBSOFA(NF,2) = NBSOFA(NF,1)
            NBSOFA(NF,3) = NBSOFA(NF,2)
            NBSOFA(NF,4) = 0
         ELSE
            NBSOFA(NF,1) = MCN ( MNFASU + WBARXQ ) + 1
            NBSOFA(NF,2) = MCN ( MNFASU + WBARYQ ) + 1
            NBSOFA(NF,3) = MCN ( MNFASU + WBARXQ ) + 1
            NBSOFA(NF,4) = MCN ( MNFASU + WBARYQ ) + 1
         END IF
C
C        LE TABLEAU 'XYZSOMMET' DE CETTE FACE
         CALL LXTSOU( NTLXSU, 'XYZSOMMET', NTSOMM, MNSOMM )
         IF( NTSOMM .LE. 0 ) THEN
            NBLGRC(NRERR) = 1
            KERR(1) = 'SURFACE SANS XYZSOMMET : ' // KNOMSU
            CALL LEREUR
            IERR = 8
            RETURN
         ENDIF
C        LE NOMBRE DE SOMMETS
         NBSO = MCN( MNSOMM + WNBSOM)
C        L'ADRESSE DU TABLEAU DES COORDONNEES DES SOMMETS DE LA FACE
         MNSOFA(NF) = MNSOMM + WYZSOM
C        TEST SUR LE NOMBRE DE SOMMETS
         IF (NUTYMA .EQ. 3) THEN
            NBSF = NBSOFA(NF,1) * ( NBSOFA(NF,1) + 1 ) / 2
            IF( NBSO .NE. NBSF ) THEN
               NBLGRC(NRERR) = 3
               WRITE(KERR(MXLGER)(1:2),'(I2)') NF
               KERR(1) = ' ERREUR : FACE' // KERR(MXLGER)(1:2)
               WRITE(KERR(MXLGER)(1:6),'(I6)') NBSO
               KERR(2) = ' DE ' // KERR(MXLGER)(1:6) // ' SOMMETS'
               WRITE(KERR(MXLGER)(1:6),'(I6)') NBSF
               KERR(1) = ' AU LIEU DE' // KERR(MXLGER)(1:6)
               CALL LEREUR
               IERR = 4
            ENDIF
         ELSE
            NBSF = NBSOFA(NF,1) * NBSOFA(NF,2)
            IF( NBSO .NE. NBSF ) THEN
               NBLGRC(NRERR) = 4
               KERR(1) = 'SURFACE : ' // KNOMSU
               WRITE(KERR(MXLGER)(1:2),'(I2)') NF
               KERR(2) = ' ERREUR : FACE' // KERR(MXLGER)(1:2)
               WRITE(KERR(MXLGER)(1:6),'(I6)') NBSO
               KERR(3) = ' DE ' // KERR(MXLGER)(1:6) // ' SOMMETS'
               WRITE(KERR(MXLGER)(1:6),'(I6)') NBSF
               KERR(4) = ' AU LIEU DE' // KERR(MXLGER)(1:6)
               CALL LEREUR
               IERR = 4
               RETURN
            ENDIF
         END IF
C
 100  ENDDO
C
      IF( IERR .GT. 0 ) RETURN
C
C    VERIFICATION DE LA GEOMETRIE
C    ============================
C
C     1) RANGEMENT DES FACES SUIVANT L'ORDRE USUEL :
C        -------------------------------------------
      DO 200 NF=1,5
C        LE TYPE DE LA FACE
         NUTYMA = NUTYPE(NF)
C        LES NUMEROS DES SOMMETS
         IF (NUTYMA .EQ. 3) THEN
            NBSOM     = 3
            NUSOFA(1) = 1
            NUSOFA(2) = NBSOFA(NF,1) * ( NBSOFA(NF,1) - 1 ) / 2  + 1
            NUSOFA(3) = NBSOFA(NF,1) * ( NBSOFA(NF,1) + 1 ) / 2
            NUSOFA(4) = 0
         ELSE
            NBSOM     = 4
            NUSOFA(1) = 1
            NUSOFA(2) = NBSOFA(NF,1)
            NUSOFA(3) = NBSOFA(NF,1) *   NBSOFA(NF,2)
            NUSOFA(4) = NBSOFA(NF,1) * ( NBSOFA(NF,2) - 1 ) + 1
         END IF
C        LES COORDONNEES DES SOMMETS
         DO NS=1,NBSOM
            IA = MNSOFA(NF) - 1 + 3 * ( NUSOFA(NS) - 1 )
            DO J=1,3
               XYZ(J,NS,NF)=RMCN(IA+J)
            ENDDO
         ENDDO
         DO NS=NBSOM+1,4
            DO J=1,3
               XYZ(J,NS,NF)=0.
            ENDDO
         ENDDO
 200  ENDDO
C
      CALL VOEXP1( NUFACE, NSENS, NPSOM, NUTYPE, XYZ, IERR )
      IF (IERR.NE.0) RETURN
C
C     2) VERIFICATION DU NOMBRE DE SOMMETS SUR CHAQUE ARETE
C        --------------------------------------------------
C     LA FACE 1 TRIANGULAIRE DE SOMMETS S1 S2 S3
      NUM   = NUFACE(1)
      NBSA1 = NBSOFA(NUM,1)
C
C     LA  FACE 2 QUADRANGULAIRE DE SOMMETS S1 S2 S5 S4
      NUM  = NUFACE(2)
      IF( NSENS(2) .EQ. 1 ) THEN
         NBSH2 = NBSOFA(NUM,NPSOM(2))
         NBSA2 = NBSOFA(NUM,QSUIV(NPSOM(2)))
      ELSE
         NBSH2 = NBSOFA(NUM,QSUIV(NPSOM(2)))
         NBSA2 = NBSOFA(NUM,NPSOM(2))
      END IF
C
C     LA  FACE 3 QUADRANGULAIRE DE SOMMETS S2 S3 S6 S5
      NUM  = NUFACE(3)
      IF( NSENS(3) .EQ. 1 ) THEN
         NBSH3 = NBSOFA(NUM,NPSOM(3))
         NBSA3 = NBSOFA(NUM,QSUIV(NPSOM(3)))
      ELSE
         NBSH3 = NBSOFA(NUM,QSUIV(NPSOM(3)))
         NBSA3 = NBSOFA(NUM,NPSOM(3))
      END IF
C
C     LA  FACE 4 QUADRANGULAIRE DE SOMMETS S3 S1 S4 S6
      NUM  = NUFACE(4)
      IF( NSENS(4) .EQ. 1 ) THEN
         NBSH4 = NBSOFA(NUM,NPSOM(4))
         NBSA4 = NBSOFA(NUM,QSUIV(NPSOM(4)))
      ELSE
         NBSH4 = NBSOFA(NUM,QSUIV(NPSOM(4)))
         NBSA4 = NBSOFA(NUM,NPSOM(4))
      END IF
C
C     LA  FACE 5 TRIANGULAIRE DE SOMMETS S4 S5 S6
      NUM   = NUFACE(5)
      NBSA5 = NBSOFA(NUM,1)
C
      NBSA = NBSA1
      NBSH = NBSH2
      IF ( NBSA2  .NE. NBSA .OR.
     %     NBSA3  .NE. NBSA .OR.
     %     NBSA4  .NE. NBSA .OR.
     %     NBSA5  .NE. NBSA .OR.
     %     NBSH3  .NE. NBSH .OR.
     %     NBSH4  .NE. NBSH ) THEN
         NBLGRC(NRERR) = 2
         KERR(1) = 'ERREUR: LE NOMBRE DE SOMMETS DOIT ETRE'
         KERR(2) = 'IDENTIQUE SUR LES ARETES COMMUNES'
         CALL LEREUR
         WRITE (IMPRIM,10100) (NF,(NBSOFA(NF,J),J=1,2),NF=1,5)
10100 FORMAT(' ERREUR : LE NOMBRE DE SOMMETS DOIT ETRE IDENTIQUE SUR'
     %,' LES ARETES COMMUNES',/,5(' FACE',I6,I4,' X ',I4,' SOMMETS'/))
         IERR = 4
         RETURN
      ENDIF
C
C     GENERATION DES SOMMETS ET DES VOLUMES
C     =====================================
      NBST = NBSA * (NBSA+1) / 2
      NBSQ = NBSA * NBSH
C     LE NOMBRE DE SOMMETS DU PENTAEDRE
      NBS  = NBST * NBSH
C
C     1) CONSTRUCTION DU TABLEAU 'XYZSOMMET' DE CE VOLUME
      CALL LXTNDC( NTLXVO, 'XYZSOMMET', 'MOTS',  WYZSOM + 3 * NBS )
      CALL LXTSOU( NTLXVO, 'XYZSOMMET',  NTSOVO, MNSOVO )
C     LE NOMBRE DE SOMMETS
      MCN( MNSOVO + WNBSOM) = NBS
C     L'ADRESSE DU DEBUT DES COORDONNEES DES SOMMETS
      MNCOSO = MNSOVO + WYZSOM
C
C     2) REMPLISSAGE DU TABLEAU 'XYZSOMMET' : LES SOMMETS DU BORD
C     RESERVATION DE LA PLACE NECESSAIRE AU TABLEAU TEMPORAIRE
C     DES 3 COORDONNEES DES SOMMETS DU PENTAEDRE UNITE
      NBSF   = MAX( NBST, NBSQ )
C     LES 2 COORDONNEES DES SOMMETS DE 2 FACES TRIANGULAIRES UNITE
      MNXYST = 0
      CALL TNMCDC( 'REEL', 2*NBST*2, MNXYST )
C     LES 2 COORDONNEES DES SOMMETS DE 2 FACES QUADRANGULAIRES UNITE
      MNXYSC = 0
      CALL TNMCDC( 'REEL', 2*NBSA*NBSH*3, MNXYSC )
C     LES 3 COORDONNEES DES SOMMETS D'UNE FACE
      MNXYZF = 0
      CALL TNMCDC( 'REEL', 3*NBSF, MNXYZF )
      CALL VOEXP2( NBSA, NBSH,
     %             RMCN(MNSOFA(NUFACE(1))), RMCN(MNSOFA(NUFACE(2))),
     %             RMCN(MNSOFA(NUFACE(3))), RMCN(MNSOFA(NUFACE(4))),
     %             RMCN(MNSOFA(NUFACE(5))),
     %             NSENS, NPSOM,
     %             RMCN(MNXYZF), RMCN(MNXYST), RMCN(MNXYST+2*NBST),
     %             RMCN(MNXYSC),
     %             RMCN(MNCOSO) )
C     DESTRUCTION DES TABLEAUX INUTILES
      CALL TNMCDS( 'REEL', 3*NBSF, MNXYZF )
C
C     3) CONSTRUCTION DU TABLEAU 'NSEF' DE CE VOLUME
      CALL LXTNDC( NTLXVO, 'NSEF', 'ENTIER' ,
     %             1 + WBARZP )
      CALL LXTSOU( NTLXVO, 'NSEF',  NTVOOB, MNVOOB )
C     TYPE DE L'OBJET : VOLUME
      MCN ( MNVOOB + WUTYOB ) = 4
C     NUMERO DU TYPE DU MAILLAGE : PENTAEDRE STRUCTURE
      MCN ( MNVOOB + WUTYMA ) = 6
C     VARIABLE NBARTP : NOMBRE D'ARETES DES COTES DES TRIANGLES
      MCN ( MNVOOB + WBARTP ) = NBSA - 1
C     VARIABLE NBARZP : NOMBRE D'ARETES DES COTES DES QUADRILATERES
      MCN ( MNVOOB + WBARZP ) = NBSH - 1
C     LE NOMBRE DE SOMMETS PAR CUBE
      MCN ( MNVOOB + WBSOEF ) = 8
C     LE NOMBRE DE TANGENTES STOCKEES PAR EF : VOLUME C0
      MCN ( MNVOOB + WBTGEF ) = 0
C     LE NOMBRE D'EF DE LA SURFACE
      MCN ( MNVOOB + WBEFOB ) =  (NBSH-1) * ((NBSA-1) ** 2)
C
C     4) INTERPOLATION TRANSFINIE DE DEGRE 1 SUR LE PENTAEDRE
      CALL VOEXP3( NBSA, NBSH,
     %             RMCN(MNXYST), RMCN(MNXYST+2*NBST),
     %             RMCN(MNXYSC),
     %             RMCN(MNCOSO) )
C
C     DESTRUCTION DES TABLEAUX INUTILES
      CALL TNMCDS( 'REEL', 2*NBST*2, MNXYST )
      CALL TNMCDS( 'REEL', 2*NBSA*NBSH*3, MNXYSC )
C
C     FIN DE LA CONSTRUCTION DU VOLUME
C     ================================
C     LA DATE DE CREATION
      CALL ECDATE( MCN(MNSOVO) )
C     LE NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNSOVO + WNBTGS ) = 0
      MCN( MNSOVO + WBCOOR ) = 3
      MCN( MNSOVO + MOTVAR(6) ) = NONMTD( '~>>>XYZSOMMET' )
C
C     LE TYPE INCONNU DE FERMETURE DU MAILLAGE
      MCN( MNVOOB + WUTFMA ) = -1
C     PAS DE TANGENTES STOCKEES
      MCN( MNVOOB + WBTGEF ) = 0
      MCN( MNVOOB + WBEFAP ) = 0
      MCN( MNVOOB + WBEFTG ) = 0
C     LA DATE DE CREATION
      CALL ECDATE( MCN(MNVOOB) )
C     LE NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNVOOB + MOTVAR(6) ) = NONMTD( '~>>>NSEF' )
C
      RETURN
      END
