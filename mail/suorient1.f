      SUBROUTINE SUORIENT1( NBEF,   NOSOEF,
     %                      NBEFAP, NBEFTG, LPEFTG, NUTGEF,
     %                      NBSOMM, XYZSOM,
     %                      L1ARET, L2ARET, LARETE,
     %                      LPTEEF, LPTEAR, LAPIAR, LAPIEF,
     %                      IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    ASSURER L'ORIENTATION DU MAILLAGE DE LA SURFACE
C -----    PAR PARCOURS DES EF A TRAVERS LES ARETES COMMUNES
C          ET PERMUTATION DES SOMMETS 2-NBSTEF DU SECOND EF
C          D'UNE ORIENTATION DIFFERENTE DU PREMIER EF
C ENTREES:
C --------
C NBEF   : NOMBRE D'EF DE LA SURFACE A TRAITER
C NOSOEF : LISTE DES 4 SOMMETS DES TRIANGLES ou QUADRANGLES DE LA SURFACE
C
C NBEFAP : Nombre des EF avec POINTEUR sur EF a TG
C NBEFTG : Nombre des EF avec TG
C LPEFTG : tableau LPEFTG(1..NBEFAP) Numero>0 de l'EF a TG sinon 0
C NUTGEF : tableau NUTGEF(1..8,1..NBEFTG) +-No des TANGENTES de l'EF a TG
C
C NBSOMM : NOMBRE DE SOMMETS DU MAILLAGE DE LA SURFACE
C XYZSOM : 3 COORDONNEES DES SOMMETS     DE LA SURFACE
C L1ARET : NOMBRE DE MOTS PAR ARETE DU TABLEAU LARETE
C L2ARET : NOMBRE DE ARETES DU TABLEAU LARETE
C LARETE : TABLEAU LARETE DES ARETES DU MAILLAGE
C          LARETE(1,I)= NO DU 1-ER  SOMMET DE L'ARETE
C          LARETE(2,I)= NO DU 2-EME SOMMET > 1-ER  SOMMET
C          LARETE(3,I)= CHAINAGE HACHAGE SUR ARETE SUIVANTE
C          LARETE(4,I)= +-NUMERO DU 1-ER EF CONTENANT CETTE ARETE
C                       + SI ARETE DIRECTE DANS L'ELEMENT, - SINON
C                       0 SI PAS DE 1-ER EF
C          LARETE(5,I)= +-NUMERO DU 2-EME EF CONTENANT CETTE ARETE
C                       + SI ARETE DIRECTE DANS L'ELEMENT, - SINON
C                       0 SI PAS DE 2-EME EF
C
C LPTEEF, LPTEAR, LAPIAR, LAPIEF: 4 TABLEAUX AUXILIAIRES
C
C SORTIES:
C --------
C NOSOEF : LISTE DES 4 SOMMETS DES FACES DE LA SURFACE
C IERR   : 0 SI PAS D'ERREUR, >0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN Embrun & St Pierre & LJLL UPMC   Janvier 2010
C2345+7...............................................................72
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      REAL              XYZSOM(3,NBSOMM)
      INTEGER           LPEFTG(NBEFAP), NUTGEF(8,NBEFTG)
      INTEGER           NOSOEF(4,NBEF), LARETE(L1ARET, L2ARET),
     %                  LPTEEF(NBEF),   LPTEAR(L2ARET),
     %                  LAPIAR(L2ARET), LAPIEF(L2ARET)
      INTEGER           NGS(2)
      DOUBLE PRECISION  D
C
ccc      print *,'EN ENTREE de suorient1.f'
ccc      DO NEF1=1, NBEF
ccc         print 10000,(K,NEF1,NOSOEF(K,NEF1),K=1,4)
ccc      ENDDO
ccc10000 FORMAT(4('NOSOEF(',I1,',',I4,')=',I4,'  '))
C
C     MISE A ZERO DU TABLEAU LPTEEF
      CALL AZEROI( NBEF, LPTEEF )
C
C     MISE A ZERO DU TABLEAU LPTEAR
      CALL AZEROI( L2ARET, LPTEAR )
C
C     LE VECTEUR NORMAL DE L'EF 1 INITIALISE L'ORIENTATION DE LA SURFACE
C     ------------------------------------------------------------------
      NEF1  = 1
      NBEFT = 0
C
C     CALCUL DU NOMBRE DES SOMMETS DE NEF1
 10   IF( NOSOEF(4,NEF1) .EQ. 0 ) THEN
         NBSEF1 = 3
      ELSE
         NBSEF1 = 4
      ENDIF
C
C     L'EF NEF1 EST TRAITE
      NBEFT = NBEFT + 1
C     L'EF NEF1 EST MARQUE COMME L'EF NBEFT TRAITE
      LPTEEF(NEF1) = NBEFT
C
C     INITIALISATION DE LA PILE DES ARETES DES EF
      LHPILA = 0
      NS1 = NOSOEF( NBSEF1, NEF1 )
      DO K=1,NBSEF1
C
C        LE NUMERO DES 2 SOMMETS DE L'ARETE K DE L'EF 1
         NS2 = NOSOEF( K, NEF1 )
         IF( NS1 .LT. NS2 ) THEN
            NGS(1) = NS1
            NGS(2) = NS2
         ELSE
            NGS(1) = NS2
            NGS(2) = NS1
         ENDIF
         CALL HACHAG( 2, NGS, L1ARET, L2ARET, LARETE,
     %                3, L2ARET, NOAR )
         IF( NOAR .LE. 0 ) THEN
            NEF2 = NEF1
            GOTO 9900
         ENDIF
C
C        L'ARETE K ET LE NO D'EF SONT EMPILES
         LHPILA = LHPILA + 1
C        L'ARETE NOAR EST EMPILEE
         LAPIAR( LHPILA ) = NOAR
C        LE NUMERO D'EF TRAITE EST EMPILE
         LAPIEF( LHPILA ) = NEF1
C
         NS1 = NS2
C
      ENDDO
C
C     PARCOURS DES ARETES DES TRIANGLES OU QUADRANGLES
C     ------------------------------------------------
 20   IF( LHPILA .GT. 0 ) THEN
C
C        LE NUMERO DE L'ARETE EMPILEE
         NOAR = LAPIAR( LHPILA )
C        LE NUMERO D'EF DE CETTE ARETE EMPILEE
         NEF1 = LAPIEF( LHPILA )
C
C        DEPILAGE
         LHPILA = LHPILA - 1
C
C        CETTE ARETE EST TRAITEE MAINTENANT
         LPTEAR( NOAR ) = 1
C
C        LE NO DES 2 EVENTUELS EF ADJACENTS A L'ARETE
         NEF2 = LARETE(4,NOAR)
         IF( ABS(NEF2) .EQ. ABS(NEF1) ) THEN
            NAR2 = 5
            NEF2 = LARETE(NAR2,NOAR)
C           RECUPERATION DU SIGNE DU TRIANGLE DE L'ARETE
            NEF1 = LARETE(4,NOAR)
         ELSE
            NAR2 = 4
C           RECUPERATION DU SIGNE DU TRIANGLE DE L'ARETE
            NEF1 = LARETE(5,NOAR)
         ENDIF
         IF( NEF2 .EQ. 0 ) THEN
C           ARETE FRONTALIERE
            GOTO 20
         ENDIF
C
C        ABS(NEF2) EST L'EF DE L'AUTRE COTE DE L'ARETE
         IF( LPTEEF( ABS(NEF2) ) .NE. 0 ) THEN
C           EF DEJA TRAITE
            GOTO 20
         ENDIF
C
C        LES 2 SOMMETS DE L'ARETE TRAITEE DANS NEF1
         D = NEF1
         D = D * NEF2
         IF( D .GT. 0D0 ) THEN
C
C           ICI LES 2 EF SONT ORIENTES DIFFEREMMENT.
C           L'ORDRE DES SOMMETS DE NEF2 EST RETOURNE
            NEF2A = ABS(NEF2)
            IF( NOSOEF(4,NEF2A) .EQ. 0 ) THEN
C              PERMUTATION DES SOMMETS 2-3 POUR UN TRIANGLE
               NBSEF2 = 3
            ELSE
C              PERMUTATION DES SOMMETS 2-4 POUR UN QUADRANGLE
               NBSEF2 = 4
            ENDIF
C           PERMUTATION DU NUMERO DES SOMMETS 2-NBSEF2
            NS2                  = NOSOEF(2,     NEF2A)
            NOSOEF(2,     NEF2A) = NOSOEF(NBSEF2,NEF2A)
            NOSOEF(NBSEF2,NEF2A) = NS2
C
C           MISE A JOUR DU SIGNE DE NEF2 DANS LARETE POUR TOUTES SES ARETES
            NS1 = NOSOEF( NBSEF2, NEF2A )
            DO K=1,NBSEF2
C
C              LE NUMERO DES 2 SOMMETS DE L'ARETE K DE L'EF NEF2
               NS2 = NOSOEF( K, NEF2A )
               IF( NS1 .LT. NS2 ) THEN
                  NGS(1) = NS1
                  NGS(2) = NS2
               ELSE
                  NGS(1) = NS2
                  NGS(2) = NS1
               ENDIF
               CALL HACHAG( 2, NGS, L1ARET, L2ARET, LARETE,
     %                      3, L2ARET, NOAR )
               IF( NOAR .LE. 0 ) THEN
                  GOTO 9900
               ENDIF
               IF( NEF2A .EQ. ABS(LARETE(4,NOAR)) ) THEN
                  NAR2 = 4
               ELSE
                  NAR2 = 5
               ENDIF
C              LE SENS DE L'ARETE EST INVERSE
               LARETE(NAR2,NOAR) = -LARETE(NAR2,NOAR)
C
C              PASSAGE A L'ARETE SUIVANTE
               NS1 = NS2
C
            ENDDO
C
            IF( NBEFTG .GT. 0 ) THEN
C              PERMUTATION DES TANGENTES DE NEF2
C              NUMERO DE L'EF A TANGENTES DANS NUTGEF
               NEFATG = LPEFTG( NEF2A )
               IF( NEFATG .GT. 0 ) THEN
C                 TRIANGLE ou QUADRANGLE SOMMET 1: TG1 <-> TG2
                  N                = NUTGEF(1,NEFATG)
                  NUTGEF(1,NEFATG) = NUTGEF(2,NEFATG)
                  NUTGEF(2,NEFATG) = N
                  IF( NOSOEF( 4, NEFATG ) .GT. 0 ) THEN
C                    QUADRANGLE SOMMET 2: TG3 <-> TG8
                     N                = NUTGEF(3,NEFATG)
                     NUTGEF(3,NEFATG) = NUTGEF(8,NEFATG)
                     NUTGEF(8,NEFATG) = N
C                    QUADRANGLE SOMMET 2: TG4 <-> TG7
                     N                = NUTGEF(4,NEFATG)
                     NUTGEF(4,NEFATG) = NUTGEF(7,NEFATG)
                     NUTGEF(7,NEFATG) = N
C                    QUADRANGLE SOMMET 3: TG5 <-> TG6
                     N                = NUTGEF(5,NEFATG)
                     NUTGEF(5,NEFATG) = NUTGEF(6,NEFATG)
                     NUTGEF(6,NEFATG) = N
                  ELSE
C                    TRIANGLE SOMMET 2: TG3 <-> TG6
                     N                = NUTGEF(3,NEFATG)
                     NUTGEF(3,NEFATG) = NUTGEF(6,NEFATG)
                     NUTGEF(6,NEFATG) = N
C                    TRIANGLE SOMMET 2: TG4 <-> TG5
                     N                = NUTGEF(4,NEFATG)
                     NUTGEF(4,NEFATG) = NUTGEF(5,NEFATG)
                     NUTGEF(5,NEFATG) = N
                  ENDIF
               ENDIF
            ENDIF
ccc         WRITE(IMPRIM,*)'PERMUTATION EF',NEF2,' ST',NOSOEF(2,NEF2),NS2
C
         ENDIF
C
C        EMPILAGE DES ARETES NON TRAITEES DE NEF2
         NEF2A = ABS( NEF2 )
         IF( NOSOEF(4,NEF2A) .EQ. 0 ) THEN
C           NEF2A EST UN TRIANGLE
            NBSEF2 = 3
         ELSE
C           NEF2A EST UN QUADRANGLE
            NBSEF2 = 4
         ENDIF
C
C        EMPILAGE DES NBSEF2 ARETES NON TRAITEES DE NEF2
         NS1 = NOSOEF( NBSEF2, NEF2A )
         DO 50 K=1,NBSEF2
C
C           LE NUMERO DES 2 SOMMETS DE L'ARETE K DE L'EF NEF2A
            NS2 = NOSOEF( K, NEF2A )
            IF( NS1 .LT. NS2 ) THEN
               NGS(1) = NS1
               NGS(2) = NS2
            ELSE
               NGS(1) = NS2
               NGS(2) = NS1
            ENDIF
            CALL HACHAG( 2, NGS, L1ARET, L2ARET, LARETE,
     %                   3, L2ARET, NOAR )
            IF( NOAR .LE. 0 ) THEN
               GOTO 9900
            ENDIF
            IF( LPTEAR(NOAR) .NE. 0 ) GOTO 40
C
C           L'ARETE K ET LE NO D'EF SONT EMPILES
            LHPILA = LHPILA + 1
            LAPIAR( LHPILA ) = NOAR
            LAPIEF( LHPILA ) = NEF2A
C
 40         NS1 = NS2
C
 50      CONTINUE
C
C        L'EF NEF2A A ETE TRAITE
         NBEFT = NBEFT + 1
         LPTEEF( NEF2A ) = NBEFT
         GOTO 20
      ENDIF
C
C     VERIFICATION DU NOMBRE EGAL D'EF AVANT et APRES
      IF( NBEF .NE. NBEFT ) THEN
C
C        LA SURFACE  N'EST PAS SIMPLEMENT CONNEXE
C        RECHERCHE D'UNE AUTRE COMPOSANTE CONNEXE
         DO 60 NEF1=1,NBEF
            IF( LPTEEF(NEF1) .EQ. 0 ) GOTO 10
 60      CONTINUE
C
C        PROBLEME!
         NBLGRC(NRERR) = 6
         WRITE(KERR(MXLGER)(1:24),'(2I12)') NBEF, NBEFT
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'SOUS-PROGRAMME SUORIENT1:'
            KERR(2) = 'NOMBRE D''EF INITIAUX' // KERR(MXLGER)( 1:12)
            KERR(3) = 'NOMBRE D''EF FINAUX  ' // KERR(MXLGER)(13:24)
            KERR(4) = 'CES 2 NOMBRES DOIVENT ETRE EGAUX'
            KERR(5) = 'VERIFIER LA CONFORMITE DU MAILLAGE OU'
            KERR(6) = 'L''ORIENTABILITE de la SURFACE'
         ELSE
            KERR(1) = 'SUBROUTINE SUORIENT1:'
            KERR(2) = 'INITIAL NUMBER of FINITE ELEMENTS'
     %                 // KERR(MXLGER)( 1:12)
            KERR(3) = 'FINAL   NUMBER of FINITE ELEMENTS'
     %                 // KERR(MXLGER)(13:24)
            KERR(4) = 'THESE 2 NUMBERS MUST BE EQUAL'
            KERR(5) = 'VERIFY the MESH CONFORMITY or'
            KERR(6) = 'the ORIENTABILITY of the SURFACE'
         ENDIF
         CALL LEREUR
         IERR = 9
      ENDIF
      GOTO 9990
C
C     SORTIES
C     =======
 9900 NBLGRC(NRERR) = 3
      WRITE(KERR(MXLGER  )(1:24),'(2I12)') NGS(1),NGS(2)
      WRITE(KERR(MXLGER-1)(1:12),'(I12)') NEF2
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = 'SUORIENT1: UNE ARETE HACHEE N''EST PAS RETROUVEE'
         KERR(2) = 'SOMMETS ' // KERR(MXLGER)(1:24)
         KERR(3) = 'ELEMENT FINI NO=' // KERR(MXLGER-1)(1:12)
      ELSE
         KERR(1) = 'SUORIENT1: AN EDGE IS NOT RECOVERED in LARETE'
         KERR(2) = 'VERTICES ' // KERR(MXLGER)(1:24)
         KERR(3) = 'FINITE ELEMENT No=' // KERR(MXLGER-1)(1:12)
      ENDIF
      WRITE(IMPRIM,19900) NGS(1), (XYZSOM(K,NGS(1)),K=1,3)
      WRITE(IMPRIM,19900) NGS(2), (XYZSOM(K,NGS(2)),K=1,3)
19900 FORMAT('XYZ',I10,'=',3G15.7)
      CALL LEREUR
      IERR = 8
      WRITE(IMPRIM,*)
      WRITE(IMPRIM,*) 'EN SORTIE DE suorient1.f avec PROBLEME'
      NEF1 = ABS(NEF1)
      NEF2 = ABS(NEF2)
      WRITE(IMPRIM,20000) (K,NEF1,NOSOEF(K,NEF1),K=1,4)
      WRITE(IMPRIM,20000) (K,NEF2,NOSOEF(K,NEF2),K=1,4)
20000 FORMAT(4('NoStEF(',I1,',',I8,')=',I8,'  '))
C
 9990 RETURN
      END
