      SUBROUTINE OBJFER( NUTYOB, NUOBJE, IMPR, NOFERM )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    RECHERCHER SI L'OBJET DE TYPE NMTOBJ EST FERME OU NON
C -----    POUR LE CAS NUTYOB=2 'LIGNE' OU 3 = 'SURFACE'
C          POINTS OU VOLUMES NE SONT PAS CONSIDERES
C          C'EST A DIRE POUR UNE
C LIGNE  : SI TOUT SOMMET D'UNE ARETE APPARTIENT A DEUX ARETES
C SURFACE: SI TOUTE ARETE D'UNE FACE APPARTIENT A 2 FACES (TRIA-QUADRANGLE)
C
C ENTREES:
C --------
C NUTYOB : 2 POUR UNE LIGNE ET 3 POUR UNE SURFACE ET ERREUR SINON
C NUOBJE : NUMERO DANS LE LEXIQUE DE L'OBJET A TRAITER
C IMPR   : 0 PAS D'IMPRESSION DES ARETES SIMPLES
C          1 AVEC  IMPRESSION DES ARETES SIMPLES
C
C SORTIE :
C --------
C NOFERM : =1 SI L'OBJET EST FERME
C          =0 SI L'OBJET N'EST PAS FERME
C          <0 SI SATURATION DU TABLEAU SERVANT AU HACHAGE
C          SI LIGNE OU SURFACE ALORS NUTFMA DU TMS 'NSEF' EST MIS A JOUR
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS    NOVEMBRE 1988
C MODIFS : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS   SEPTEMBRE 1995
C MODIFS : PERRONNET ALAIN TEXAS A&M UNIVERSITY at QATAR       MARS 2011
C MODIFS : ALAIN PERRONNET LJLL UPMC & St Pierre du Perray Novembre 2015
C2345X7..............................................................012
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/trvari.inc"
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/ntmnlt.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
      CHARACTER*10      NMTOBJ,NMTYOB
      CHARACTER*24      KNOM
      INTEGER           NBARXF(3)
C
C     PROTECTION CONTRE UNE SUPPRESSION D'UN TABLEAU NON DECLARE
      MNARFA = 0
      L1ARFA = 0
      L2ARFA = 0

C     NOMBRE D'ARETES APPARTENANT A X FACES
C     NBARXF(3) NOMBRE D'ARETES APPARTENANT A PLUS DE 2 FACES
      DO K=1,3
         NBARXF(K) = 0
      ENDDO
      NOFERM = -1
C
C     LE NUMERO DE L'OBJET
      NMTOBJ = NMTYOB( NUTYOB )
      IF( NUTYOB .NE. 2 .AND. NUTYOB .NE. 3 ) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:4),'(I4)') NUTYOB
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'OBJFER:TYPE '//KERR(MXLGER)(1:4)//' NON TRAITE'
         ELSE
           KERR(1) ='OBJFER:TYPE '//KERR(MXLGER)(1:4)//' NOT PROGRAMMED'
         ENDIF
         CALL LEREUR
         RETURN
      ENDIF
C
C     NTMN(NUTYOB) LE NUMERO DU TAMS ET ADRESSE MCN DU LEXIQUE DE CES OBJETS
C     LE LEXIQUE DE L'OBJET
      CALL LXNLOU( NTMN(NUTYOB), NUOBJE, NTLXOB, MNLXOB )
      IF( NTLXOB .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:10),'(I10)') NUOBJE
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) =  'OBJFER:'//NMTOBJ//' INCONNU: NUMERO '
     %            // KERR(MXLGER)(1:10)
         ELSE
            KERR(1) =  'OBJFER:'//NMTOBJ//' UNKNOWN: NUMBER '
     %            // KERR(MXLGER)(1:10)
         ENDIF
         CALL LEREUR
         RETURN
      ENDIF
      CALL NMOBNU( NMTOBJ, NUOBJE, KNOM )
C
C     LE TABLEAU 'NSEF' DE L'OBJET
      CALL LXTSOU( NTLXOB, 'NSEF', NTNSEF, MNNSEF )
      IF( NTNSEF .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:4),'(I4)') NUOBJE
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'OBJFER:OBJET NON MAILLE '//NMTOBJ
     %           //KERR(MXLGER)(1:4)
         ELSE
            KERR(1) = 'OBJFER: NOT MESHED OBJECT '//NMTOBJ
     %           //KERR(MXLGER)(1:4)
         ENDIF
         CALL LEREUR
         RETURN
      ENDIF
C
C     LE TABLEAU 'XYZSOMMET' DE L'OBJET
      CALL LXTSOU( NTLXOB, 'XYZSOMMET', NTXYZS, MNXYZS )
      IF( NTXYZS .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:4),'(I4)') NUOBJE
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'OBJFER: OBJET SANS SOMMETS '//NMTOBJ//
     %              KERR(MXLGER)(1:4)
         ELSE
            KERR(1) = 'OBJFER: OBJECT without VERTICES '//NMTOBJ//
     %              KERR(MXLGER)(1:4)
         ENDIF
         CALL LEREUR
         RETURN
      ENDIF
C
C     DIMENSION NDIM DE L'ESPACE DES COORDONNEES
      NBSOM = MCN(MNXYZS+WNBSOM)
      IF( NBSOM .LE. 0 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            PRINT*,'objfer: INCORRECT Nombre Sommets=',NBSOM,' de ',KNOM
         ELSE
            PRINT*,'objfer: INCORRECT Vertice Number=',NBSOM,' of ',KNOM
         ENDIF
         NOFERM = 0
         GOTO 9999
      ENDIF
      CALL DIMCOO( NBSOM, MCN(MNXYZS+WYZSOM), NDIM )
C
C     LE TYPE DE FERMETURE DU MAILLAGE DE LA LIGNE OU SURFACE
      NUTFMA = MCN( MNNSEF + WUTFMA )
C     NUTFMA = 1  LA LIGNE OU SURFACE A DEJA ETE TESTEE COMME ETANT FERMEE
C            = 0  LA LIGNE OU SURFACE A DEJA ETE TESTEE COMME ETANT NON-FERMEE
C            =-1  LA FERMETURE DE LA LIGNE EST A TESTER
      NOFERM = NUTFMA
      IF( NUTFMA .GE. 0 ) GOTO 100
C
C     LE TYPE FERMETURE EST INCONNU => IL EST CALCULE
C     ===============================================
C     LE TYPE DU MAILLAGE
      NUTYMA = MCN( MNNSEF + WUTYMA )
C
      IF( NUTYMA .EQ. 0 ) THEN
C        LIGNE OU SURFACE NON STRUCTUREE
C        NOMBRE DE NSEF
         NBEFOB = MCN( MNNSEF + WBEFOB )
C        NOMBRE DE SOMMETS ET TANGENTES PAR EF
         NBSOEF = MCN( MNNSEF + WBSOEF )
      ELSE IF( NUTYMA .LT. 0 ) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:4),'(I4)') NUTYMA
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'OBJFER: TYPE INCORRECT de MAILLAGE '//
     %              KERR(MXLGER)(1:4)
         ELSE
            KERR(1) = 'OBJFER: INCORRECT TYPE of MESH '//
     %              KERR(MXLGER)(1:4)
         ENDIF
         CALL LEREUR
         NOFERM = -1
         GOTO 9000
      ELSE IF( NUTYOB .EQ. 2 .OR. NUTYOB .EQ. 3 ) THEN
C        LIGNE OU SURFACE STRUCTUREE DONC OUVERTE
         NOFERM = 0
         GOTO 100
      ELSE
C        POINT OU VOLUME => OBJET NON FERME
         NOFERM = 0
         GOTO 300
      ENDIF
C
      IF( NBEFOB.LE. 0 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            PRINT*,'objfer: INCORRECT Nombre EF=',NBEFOB,' de ',KNOM
         ELSE
            PRINT*,'objfer: INCORRECT FE Number=',NBEFOB,' of ',KNOM
         ENDIF
         NOFERM = 0
         GOTO 9999
      ENDIF

      IF( NUTYOB .EQ. 2 ) THEN
C
C        LE TABLEAU DES SOMMETS DES ARETES DE LA LIGNE
C        ---------------------------------------------
         MNARFA = 0
         CALL GESOAR( RMCN(MNXYZS+WYZSOM),
     &                NBSOEF, NBEFOB, MCN(MNNSEF+WUSOEF), 2,
     &                L1ARFA, L2ARFA, MNARFA, NBARXF(3) )
         IF( NBARXF(3) .NE. 0 ) THEN
            NOFERM = -1
            GOTO 100
         ENDIF
C
C        RECHERCHE DES SOMMETS APPARTENANT A UNE SEULE ARETE
C        LES SOMMETS APPARTENANT A PLUS DE 2 ARETES ONT DEJA
C        ETE DIAGNOSTIQUES
         MN = MNARFA
         DO 40 NA = 1, L2ARFA
C
            IF( MCN(MN) .NE. 0 ) THEN
C
C              LE SOMMET EST INITIALISE
               IF( MCN(MN+3) .EQ. 0 ) THEN
C
C                 LE SOMMET APPARTIENT A UNE SEULE ARETE
                  NBARXF(1) = NBARXF(1) + 1
C
C                 AFFICHAGE DE L'ARETE SIMPLE
                  IF( IMPR .NE. 0 ) THEN
                     NOS = MCN(MN)
                     IF( LANGAG .EQ. 0 ) THEN
                        PRINT*,'objfer: Le SOMMET',NOS,
     %                      ' APPARTIENT a la SEULE ARETE',NOS,MCN(MN+2)
                     ELSE
                        PRINT*,'objfer: The VERTEX',NOS,
     %                         ' BELONGS to ONLY ONE EDGE',NOS,MCN(MN+2)
                     ENDIF
                     MNST = MNXYZS + WYZSOM + 3 * NOS - 3
                   PRINT 10040, NOS,RMCN(MNST),RMCN(MNST+1),RMCN(MNST+2)
                  ENDIF
C
               ENDIF
            ENDIF
C
            MN = MN + L1ARFA
C
 40      CONTINUE
10040    FORMAT(' objfer: POINT',I9,
     %          ' : X=',G15.7,' Y=',G15.7,' Z=',G15.7)
C
         IF( NBARXF(1) .EQ. 0 ) THEN
            NOFERM = 1
         ELSE
            NOFERM = 0
         ENDIF
         GOTO 105

      ELSE

C        LE TABLEAU DES ARETES DES FACES DE LA SURFACE
C        ---------------------------------------------
         MXFAAR = 4
         CALL GEARFA( MCN(MNXYZS+WYZSOM), NBSOEF, NBEFOB,
     %                MCN(MNNSEF+WUSOEF), MXFAAR,
     %                L1ARFA, L2ARFA, MNARFA, NBARXF, IERR )
         IF( IERR .LT. 0 ) GOTO 110

C        RECHERCHE DES ARETES APPARTENANT A UNE SEULE FACE
         MN = MNARFA
         DO NA = 1, L2ARFA

            IF( MCN(MN) .NE. 0 ) THEN

C              L'ARETE NA EST INITIALISEE
               IF( MCN(MN+4) .EQ. 1 ) THEN

C                 L'ARETE NA APPARTIENT A UNE SEULE FACE
C                 AFFICHAGE DE L'ARETE SIMPLE
                  IF( IMPR .NE. 0 ) THEN

                     WRITE(KERR(MXLGER)( 1:10), '(I10)') MCN(MN)
                     WRITE(KERR(MXLGER)(11:20), '(I10)') MCN(MN+1)
                     WRITE(KERR(MXLGER-1)(1:10),'(I10)') MCN(MN+3)

                     IF( LANGAG .EQ. 0 ) THEN
                        PRINT*,'objfer: L''ARETE',MCN(MN),MCN(MN+1),
     %                         ' APPARTIENT A UNE SEULE FACE',MCN(MN+3)
                     ELSE
                        PRINT*,'objfer: The EDGE',MCN(MN),MCN(MN+1),
     %                         ' BELONGS to ONLY ONE FACE',MCN(MN+3)
                     ENDIF
                     MNST   = MNXYZS + WYZSOM + 3 * MCN(MN) - 3
                     PRINT 10040, MCN(MN),RMCN(MNST),RMCN(MNST+1),
     &                            RMCN(MNST+2)
                     MNST   = MNXYZS + WYZSOM + 3 * MCN(MN+1) - 3
                     PRINT 10040, MCN(MN+1),
     %                            RMCN(MNST),RMCN(MNST+1),RMCN(MNST+2)
                  ENDIF

               ENDIF
            ENDIF

C           PASSAGE A L'ARETE SUIVANTE
            MN = MN + L1ARFA

         ENDDO

         IF( NBARXF(1) .EQ. 0 .AND. NBARXF(3) .EQ. 0 ) THEN
C           SURFACE MAILLEE FERMEE
            NOFERM = 1
         ELSE
C           SURFACE MAILLEE NON FERMEE
            NOFERM = 0
         ENDIF

         GOTO 110
      ENDIF
C
C     AFFICHAGE DU TYPE DE FERMETURE
C     ------------------------------
 100  IF( NUTYOB .EQ. 2 ) THEN
         GOTO 105
      ELSE
         GOTO 110
      ENDIF

C     AFFICHAGE DU TYPE DE FERMETURE DE LA LIGNE
C     ------------------------------------------
 105  IF( LORBITE .NE. 0 ) GOTO 300
      NBC = NUDCNB(KNOM)
      IF( NOFERM .EQ. 1 ) THEN

C        LIGNE FERMEE
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10105) KNOM(1:NBC),NDIM
10105       FORMAT(' ',A,' est une LIGNE FERMEE de R**',I1,
     %             '  .................................')
         ELSE
            WRITE(IMPRIM,20105) KNOM(1:NBC),NDIM
20105       FORMAT(' ',A,' is a CLOSED LINE of R**',I1,
     %             '  .....................................')
         ENDIF

      ELSE IF( NOFERM .EQ. 0 ) THEN

C        LIGNE NON-FERMEE
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10106) KNOM(1:NBC),NDIM
10106        FORMAT(' ',A,' est une LIGNE NON FERMEE de R**',I1,
     %              '  .............................')
         ELSE
            WRITE(IMPRIM,20106) KNOM(1:NBC),NDIM
20106        FORMAT(' ',A,' is NOT a CLOSED LINE of R**',I1,
     %              '  .................................')
         ENDIF

      ELSE IF( NOFERM .LT. 0 ) THEN

C        LIGNE DE FERMETURE INCONNUE
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10107) KNOM(1:NBC),NDIM
10107       FORMAT(' ',A,' est une LIGNE de R**',I1,
     %       ' de FERMETURE INCONNUE  .............................')
         ELSE
            WRITE(IMPRIM,20107) KNOM(1:NBC),NDIM
20107       FORMAT(' ',A,' is a LINE of R**',I1,
     %       ' with UNKNOMN CLOSING  ..............................')
         ENDIF

      ENDIF

C     MISE A JOUR DU TYPE DE FERMETURE DANS LE TMS NSEF
      GOTO 300

C     AFFICHAGE DU TYPE DE FERMETURE DE LA SURFACE
C     --------------------------------------------
 110  IF( LORBITE .NE. 0 ) GOTO 300
      NBC = NUDCNB(KNOM)
      IF( NOFERM .EQ. 1 ) THEN
C
C        SURFACE FERMEE
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10110) KNOM(1:NBC),NDIM
10110       FORMAT(' ',A,' est une SURFACE FERMEE de R**',I1,
     %              '  ...............................')
         ELSE
            WRITE(IMPRIM,20110) KNOM(1:NBC),NDIM
20110       FORMAT(' ',A,' is a CLOSED SURFACE of R**',I1,
     %             '  ..................................')
         ENDIF
C
      ELSE IF( NOFERM .EQ. 0 ) THEN
C
C        SURFACE NON-FERMEE
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10111) KNOM(1:NBC),NDIM
10111       FORMAT(' ',A,' N''est PAS une SURFACE FERMEE de R**',I1,
     %              '  ........................')
         ELSE
            WRITE(IMPRIM,20111) KNOM(1:NBC),NDIM
20111       FORMAT(' ',A,' is NOT a CLOSED SURFACE of R**',I1,
     %             '  ..............................')
         ENDIF
C
      ELSE IF( NOFERM .LT. 0 ) THEN
C
C        SURFACE DE FERMETURE INCONNUE
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10112) KNOM(1:NBC),NDIM
10112       FORMAT(' ',A,' est une SURFACE de R**',I1,
     %       ' de FERMETURE INCONNUE  .............................')
         ELSE
            WRITE(IMPRIM,20112) KNOM(1:NBC),NDIM
20112       FORMAT(' ',A,' is a SURFACE of R**',I1,
     %       ' with UNKNOMN CLOSING  ..............................')
         ENDIF
C
      ENDIF
C
C     MISE A JOUR DU TYPE DE FERMETURE DANS LE TMS NSEF
 300  IF( MNNSEF .GT. 0 ) MCN( MNNSEF + WUTFMA ) = NOFERM
C
C     TRACE DES ARETES, DES ARETES DE 3 FACES OU DES ARETES SIMPLES
C     -------------------------------------------------------------
ccc      IF( NBARXF(1) .GE. 0 .OR. NBARXF(3) .GE. 0 ) THEN
ccc         NBSOM = MCN(MNXYZS+WNBSOM)
ccc         CALL DIMCOO( NBSOM,  MCN(MNXYZS+WYZSOM), NDIMLI )
ccc         CALL TRAR13F( NDIMLI, NBSOM,  MCN(MNXYZS+WYZSOM),
ccc     %                 L1ARFA, L2ARFA, MCN(MNARFA),
ccc     %                 NBARXF(1), NBARXF(3), MCN(MNAR3F) )
ccc      ENDIF
C
C     DESTRUCTION DES TABLEAUX DEVENUS INUTILES
 9000 WRITE(IMPRIM,*)
      IF( MNARFA .GT. 0 ) CALL TNMCDS('ENTIER', L1ARFA*L2ARFA, MNARFA )
C
 9999 RETURN
      END
