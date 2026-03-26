      SUBROUTINE GEFASU( MNTSMA, MXMOFA, L1FACE, L2FACE, MNFACE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    FORMER PAR HACHAGE LE TABLEAU DES FACES D'UNE SURFACE
C -----
C ENTREES:
C --------
C MNTSMA : ADRESSE MCN DU TABLEAU 'NSEF' DU MAILLAGE
C MXMOFA : NOMBRE DE MOTS EN PLUS PAR FACE
C          0 SI PAS D'UN TEL STOCKAGE
C
C SORTIES:
C --------
C L1FACE : NOMBRE DE MOTS PAR FACE DU TABLEAU NFACE
C L2FACE : NOMBRE DE FACES DU TABLEAU NFACE
C MNFACE : ADRESSE DANS M DU TABLEAU NFACE DES FACES DU MAILLAGE
C          EN SORTIE NFACE(1,I)= NO DU 1-ER  SOMMET DE LA FACE
C                    NFACE(2,I)= NO DU 2-EME SOMMET > 1-ER  SOMMET
C                    NFACE(3,I)= NO DU 3-EME SOMMET DE LA FACE
C                    NFACE(4,I)= NO DU 4-EME SOMMET DE LA FACE
C                                QUELCONQUE SI TRIANGLE
C                    NFACE(5,I)= CHAINAGE HACHAGE SUR LA FACE SUIVANTE
C                    NFACE(5:5+MXMOFA,I)= NON INITIALISE DANS
C                                         CE SOUS PROGRAMME
C ATTENTION : CES 3 ENTIERS SONT NULS SI UNE ERREUR EST RENCONTREE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS  NOVEMBRE 1988
C ......................................................................
      IMPLICIT INTEGER (W)
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      include"./incl/a___nsef.inc"
      INTEGER           NOSOFA(1:4)
C
      NBARXQ = 0
      NBARX1 = 0
C
C     LE TYPE DU MAILLAGE
      NUTYMA = MCN( MNTSMA + WUTYMA )
C
      IF( NUTYMA .EQ. 0 ) THEN
C
C        SURFACE NON STRUCTUREE( LA SURFACE PEUT ETRE EN PLUSIEURS MORCEAUX)
C        NOMBRE DE SOMMETS ET TGS PAR EF
         NBSOEF = MCN( MNTSMA + WBSOEF )
C        NOMBRE DE NSEF
         NBFASU = MCN( MNTSMA + WBEFOB )
C        ADRESSE DU NUMERO DU 1-ER SOMMET DE LA 1-ERE FACE
         MNSS   = MNTSMA + WUSOEF
C
      ELSE IF( NUTYMA .EQ. 3 .OR. NUTYMA .EQ. 4 ) THEN
C
C        TRIANGLE OU QUADRANGLE STRUCTURE
C        NOMBRE DE FACES
         IF( NUTYMA .EQ. 3 ) THEN
C           TRIANGLE
            NBARTR = MCN ( MNTSMA + WBARTR )
            NBFASU = NBARTR * NBARTR
            NBSOEF = 4
            MNSS   = 0
            CALL TNMCDC( 'ENTIER', NBSOEF * NBFASU , MNSS )
            MNSS0  = MNSS
            CALL SUEXT3(MNTSMA,MNSS)
         ELSE
C           QUADRANGLE
            NBARXQ = MCN( MNTSMA + WBARXQ )
            NBARX1 = NBARXQ + 1
            NBARYQ = MCN( MNTSMA + WBARYQ )
            NBFASU = NBARXQ * NBARYQ
            NBSOEF = 4
            NX     = 0
            NXY    = 0
         ENDIF
C
      ELSE
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:4),'(I4)') NUTYMA
         KERR(1) = 'GEFASU:TYPE INCORRECT DE NSEF'
     %           // KERR(MXLGER)(1:4)
         CALL LEREUR
         L1FACE = 0
         L2FACE = 0
         MNFACE = 0
         RETURN
      ENDIF
C
C     CALCUL DU NOMBRE DES FACES DU MAILLAGE
      L1FACE = 5 + MXMOFA
      L2FACE = MAX( NBFASU , 10 )
C
C     ADRESSAGE DU TABLEAU FACE
      L      = L1FACE * L2FACE
      MNFACE = 0
      CALL TNMCDC( 'ENTIER' , L , MNFACE )
C
C     LE TABLEAU DES FACES EST INITIALISE A ZERO
      CALL AZEROI( L , MCN(MNFACE) )
C
C     LA 1-ERE FACE LIBRE DERRIERE CELLES ADRESSEES PAR LE MINIMUM
      LIBREF = L2FACE
C
C     FORMATION DU TABLEAU DES NO DES SOMMETS DES FACES
C     ==================================================
      DO 100 N=1,NBFASU
C
C        TRAITEMENT SELON LE TYPE DU MAILLAGE
         IF( NUTYMA .EQ. 0 ) THEN
C
C           MAILLAGE NON STRUCTURE
C           ----------------------
            K = MCN( MNSS + 3 )
C           LE NUMERO DES SOMMETS DE LA FACE NF
            IF( K .GT. 0 ) THEN
               NBSOFA = 4
            ELSE
               NBSOFA = 3
            ENDIF
C
C           BOUCLE SUR LES SOMMETS DE LA FACE
C           PERMUTATION CIRCULAIRE DES SOMMETS POUR AMENER
C           LE PLUS PETIT NO DE SOMMET DE LA FACE EN PREMIER
            K     = 1
            L     = MCN( MNSS )
            MNSS1 = MNSS - 1
            DO 30  J=2,NBSOFA
               IF( L .LE. MCN(MNSS1+J) ) GOTO 30
               L = MCN(MNSS1+J)
               K = J
   30       CONTINUE
            CALL TRTATA( MCN(MNSS1+K) , NOSOFA(1) , NBSOFA-K+1 )
            IF( K .NE. 1 ) THEN
               CALL TRTATA( MCN(MNSS) , NOSOFA(NBSOFA-K+2) , K-1 )
            ENDIF
C           PASSAGE A LA FACE SUIVANTE
            MNSS = MNSS + NBSOEF
C
         ELSE IF( NUTYMA .EQ. 3 ) THEN
C
C           TRIANGLE STRUCTURE
C           ------------------
            NBSOFA = 3
            NOSOFA(1) = MCN( MNSS     )
            NOSOFA(2) = MCN( MNSS + 1 )
            NOSOFA(3) = MCN( MNSS + 2 )
C           PASSAGE A LA FACE SUIVANTE
            MNSS = MNSS + 4
C
         ELSE IF( NUTYMA .EQ. 4 ) THEN
C
C           QUADRANGLE STRUCTURE
C           --------------------
            NBSOFA = 4
            NXY    = NXY + 1
            NX     = NX  + 1
            IF( NX .GT. NBARXQ ) THEN
               NX  = 1
               NXY = NXY + 1
            ENDIF
            NOSOFA(1) = NXY
            NOSOFA(2) = NXY + 1
            NOSOFA(3) = NXY + 1 + NBARX1
            NOSOFA(4) = NXY + NBARX1
         ENDIF
C
C        ADJONCTION DE LA FACE SI ELLE N'EXISTE PAS DEJA
C        -----------------------------------------------
         CALL HACHAG( NBSOFA,NOSOFA,L1FACE,L2FACE,MCN(MNFACE),5,
     &                LIBREF,J )
  100 CONTINUE
C     DESTRUCTION DU TABLEAU DES NSEF
      IF (NUTYMA.EQ.3) THEN
         CALL TNMCDS( 'ENTIER' , NBSOEF * NBFASU , MNSS0 )
      END IF
      END
