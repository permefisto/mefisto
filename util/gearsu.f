      SUBROUTINE GEARSU( NBSOEF, NBFASU, NOSOEF, MXMOAR,
     %                   L1ARET, L2ARET, MNARET, IERR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : FORMER PAR HACHAGE LE TABLEAU DES ARETES DU MAILLAGE D'UNE SURFACE
C -----
C ENTREES:
C --------
C NBSOEF : NOMBRE DE SOMMETS ET TANGENTES PAR FACE
C NBFASU : NOMBRE DE EF DE LA SURFACE
C NOSOEF : LES NBSOEF NUMEROS DES SOMMETS DES NBFASU FACES DE LA SURFACE
C MXMOAR : NOMBRE DE MOTS EN PLUS (ou EF) PAR ARETE
C          2 PAR EXEMPLE POUR LE NUMERO DES 2 EF CONTENANT L'ARETE
C          0 SI PAS D'UN TEL STOCKAGE

C SORTIES:
C --------
C L1ARET : NOMBRE DE MOTS PAR ARETE DU TABLEAU LARETE
C L2ARET : NOMBRE DE ARETES DU TABLEAU LARETE
C MNARET : ADRESSE DANS MCN DU TABLEAU LARETE=MCN(MNARET) DES ARETES
C          LARETE(1,I)= NO DU 1-ER  SOMMET DE L'ARETE
C          LARETE(2,I)= NO DU 2-EME SOMMET > 1-ER  SOMMET
C          LARETE(3,I)= CHAINAGE HACHAGE SUR ARETE SUIVANTE
C          LARETE(4,I)=+-NUMERO DU 1-ER EF CONTENANT CETTE ARETE
C                      + SI ARETE   DIRECTE DANS L'EF,
C                      - SI ARETE INDIRECTE DANS L'EF,
C                      0 SI ARETE APPARTENANT A AUCUN EF
C           LARETE(4:3+MXMOAR,I)= +-NUMERO DE L'EF DE CETTE ARETE
C                                = IDEM
C IERR   : =0 SI AUCUNE ARETE APPARTIENT A PLUS DE MXMOAR ELEMENTS FINIS
C          =1 SI    UNE ARETE APPARTIENT A PLUS DE MXMOAR ELEMENTS FINIS
C             ou L1ARET-3 < 2
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 1988
C2345+7...............................................................72
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/ponoel.inc"
      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      INTEGER           NOSOEF(NBSOEF,NBFASU)
      INTEGER           NGS(2), NBARXF(6)

      IERR = 0

C     CALCUL DU NOMBRE DES ARETES DU MAILLAGE
C     ---------------------------------------
C     MXMOAR : ICI NOMBRE MAXIMAL D ELEMENTS CONTENANT UNE MEME ARETE
      L1ARET = 3 + MXMOAR

C     NOMBRE D' ARETES : FORMULE A MODIFIER EVENTUELLEMENT
C     ---------------------------------------------------
      L2ARET = NBFASU * 4

C     ADRESSAGE DES TABLEAUX KARET ET DU POINTEUR DES ARETES REFEREES
C     ---------------------------------------------------------------
      L1     = L1ARET * L2ARET
      MNARET = 0
      CALL TNMCDC( 'ENTIER', L1, MNARET )

C     LE TABLEAU DES ARETES EST INITIALISE A ZERO
      CALL AZEROI( L1, MCN(MNARET) )

C     LA 1-ERE ARETE LIBRE EST LA DERNIERE DU TABLEAU
      LIBREA = L2ARET

C     FORMATION DU TABLEAU DES NO DES SOMMETS DES ARETES
C     ==================================================
      DO 100 N = 1, NBFASU

C        L'ELEMENT FINI N EST IL ACTIF
         IF( NOSOEF(1,N) .EQ. 0 ) GOTO 100

C        LE NOMBRE REEL NBSTEF DE SOMMETS DE L ELEMENT FINI N
C        (QUADRANGLE ou TRIANGLE de 4-eme SOMMET NUL)
         DO NBSTEF=NBSOEF,1,-1
            IF( NOSOEF(NBSTEF,N) .GT. 0 ) GOTO 20
         ENDDO
C
 20      IF( NBSTEF .EQ. 3 ) THEN
            NUTYEL = 13
         ELSE IF( NBSTEF .EQ. 4 ) THEN
            NUTYEL = 16
         ELSE
            GOTO 100
         ENDIF

C        LES CARACTERISTIQUES DE L' ELEMENT FINI N
C        NOTAMMENT: NARET=NOMBRE D'ARETES DE L'EF
         CALL ELTYCA( NUTYEL )

C        BOUCLE SUR LES ARETES DE L ELEMENT FINI
C        ---------------------------------------
         DO I = 1, NARET
            DO J=1,2
               NGS(J) = NOSOEF( NOSOAR(J,I) , N )
            ENDDO

C           PERMUTATION DES SOMMETS POUR AMENER
C           LE PLUS PETIT NO EN PREMIER
            L     = NGS(1)
            ISENS = 1
            IF( L .GT. NGS(2) ) THEN
C               ARETE INDIRECTE
                NGS(1) = NGS(2)
                NGS(2) = L
                ISENS  = -1
            ENDIF

            CALL HACHAG( 2, NGS, L1ARET, L2ARET, MCN(MNARET), 3,
     %                   LIBREA, NOAR )
C           NOAR =0 SI SATURATION DU TABLEAU MNARET
C                >0 SI LE TABLEAU NGS A ETE RETROUVE
C                <0 SI LE TABLEAU NGS A ETE AJOUTE

            MN = MNARET+ L1ARET * ( ABS(NOAR) - 1 )

C           STOCKAGE ADRESSE DE L ELEMENT DANS LARETE(4,...;NOAR)
C           -----------------------------------------------------
            L1 = 2
   35       L1 = L1 + 1
            IF( L1 .LT. L1ARET ) THEN
               IF( MCN(MN + L1) .NE. 0 ) GOTO 35

C              ADRESSE +1 > 0 SI L'ARETE EST   DIRECTE DANS L ELEMENT FINI
C              ADRESSE +1 < 0 SI L'ARETE EST INDIRECTE DANS L ELEMENT FINI
               MCN( MN + L1 ) = N * ISENS

            ELSE

               NBLGRC(NRERR) = 5
               WRITE(KERR(2),'( 2I9 )') NGS
               WRITE(KERR(4),'( 8I9 )') (ABS(MCN(MN+2+L)),L=1,MXMOAR),N
               IF( LANGAG .EQ. 0 ) THEN
                  KERR(1) = 'gearsu: UNE ARETE DE SOMMETS'
                  KERR(3) = 'APPARTIENT AU MOINS AUX ELEMENTS FINIS'
                  KERR(5) = 'CORRIGER LE MAILLAGE DE LA SURFACE'
               ELSE
                KERR(1)='gearsu: AN EDGE of 2 VERTICES'
                KERR(3)='BELONGS TO AT LEAST THE FINITE ELEMENTS'
                KERR(5)='CORRECT THE SURFACE MESH'
               ENDIF
               CALL LEREUR
               PRINT 10060,N,(NOSOEF(L,N),L=1,NBSOEF),(NGS(L),L=1,2),
     %                     MXMOAR

C              AFFICHAGE DU NO DES SOMMETS DES EF
               DO L=1,MXMOAR
                  NF = ABS( MCN(MN+2+L) )
                  PRINT*,'gearsu: EF',NF,':',(NOSOEF(LL,NF),LL=1,NBSOEF)
               ENDDO
               PRINT*,'gearsu: EF',N,':',(NOSOEF(LL,N),LL=1,NBSOEF)

               IERR = 3
ccc               IERR = 3  ON CONTINUE  2020/04/09 ...
ccc               RETURN
C
            ENDIF
         ENDDO
10060    FORMAT(' gearsu: LA FACE',I9,' DE SOMMETS',4I9,
     %          ' POUR L''ARETE',2I9/
     %          ' AUGMENTER MXMOAR=',I1,' DU SP gearsu'/)

C        PASSAGE A L ELEMENT FINI SUIVANT
  100 ENDDO

C     CALCUL DU NOMBRE D'ARETES APPARTENANT A 1, 2, 3, 4, 5, 6 QTANGLES
C     -----------------------------------------------------------------
      CALL NBARXFA( L1ARET, L2ARET, MCN(MNARET),
ccc     %              NBSOEF, NOSOEF,
     %              NBARXF, NOARPB, IERR )

      RETURN
      END
