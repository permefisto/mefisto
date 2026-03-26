      SUBROUTINE RECUARPE( NBSOMM, MXSOMM, PTXYZD, NPSOFR,
     %                     N1TETS, NOTETR, N1TEVI, NUDTETR,
     %                     INFACO, MXFACO, LEFACO, N1FASC,
     %                     IVOLTE, NVOLTE, MXTE1S, NOTE1S,
     %                     NBFAPE, NOFAPE, NBARRECU, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    RECUPERATION D'ARETES PERDUES PAR 2T -> 3T MULTIPLES
C -----

C ENTREES:
C --------
C NBSOMM : NUMERO DU DERNIER SOMMET AJOUTE A LA TETRAEDRISATION
C MXSOMM : NOMBRE MAXIMAL DE SOMMETS PERMIS POUR LA TETRAEDRISATION
C PTXYZD : X Y Z DISTANCE SOUHAITEE DES POINTS
C NPSOFR : =  0 SI LE POINT EST AJOUTE NON SUR LA SURFACE DE L'OBJET
C          =  1 SI LE POINT EST AJOUTE SUR UNE ARETE DE LA SURFACE
C          =  2 SI LE POINT EST AJOUTE COMME BARYCENTRE SUR LA SURFACE
C          =  1 000 000  + NUMERO DU  POINT INTERNE UTILISATEUR
C          = (1 000 000) * (NO SURFACE + 1) + NUMERO DU SOMMET DANS
C              LA NUMEROTATION DES SOMMETS DE LA SURFACE
C          = -4 SI LE POINT EST SOMMET D'OT NON TROP PROCHE PT OU FACE
C          = -1 SI LE POINT EST SOMMET D'OT TROP PROCHE PT OU FACE
C          = -3 SI LE POINT EST SOMMET D'OT REFUSE DANS LA TETRAEDRISATION
C          = -NPSOFR(I) SI POINT DEPLACE SUR LA SURFACE
C             DANS LEUR SURFACE FERMEE INITIALE ou NO DE POINT INTERNE
C INFACO : =1 LE TABLEAU LEFACO DOIT ETRE PRESENT
C MXFACO : NOMBRE MAXIMAL DE FACES DECLARABLES DANS LEFACO
C LEFACO : LES 3 SOMMETS, 2 MATERIAUX, 3 FACES VOISINES ET CHAINAGE
C          DES FACES TRIANGULAIRES DU CONTOUR ET INTERFACES
C N1FASC : N1FASC(NS)=NUMERO (DANS LEFACO) D'UNE FACE DE SOMMET NS
C IVOLTE : 0 PAS       DE TABLEAU NVOLTE A L'APPEL
C          1 EXISTENCE DU TABLEAU NVOLTE A L'APPEL
C NVOLTE : NUMERO DU VOLUME (1 A NBVOPA) DE CHAQUE TETRAEDRE
C          -1 SI VOLUME INCONNU
C MXTE1S : MAX DE MOTS DU TABLEAU NOTE1S

C MODIFIES:
C ---------
C N1TETS : N1TETS(I) NUMERO D'UN TETRAEDRE AYANT POUR SOMMET I
C NOTETR : LISTE DES TETRAEDRES
C          SOMMET1,    SOMMET2,    SOMMET3,    SOMMET4,
C          TETRAEDRE1, TETRAEDRE2, TETRAEDRE3, TETRAEDRE4
C          DE L'AUTRE COTE DE LA FACE
C          1: 123      2: 234      3: 341      4: 412
C N1TEVI : NUMERO DU 1 PREMIER TETRAEDRE VIDE DANS LE TABLEAU NOTETR
C          LE CHAINAGE DES TETRAEDRES VIDES SE FAIT SUR NOTETR(5,.)
C NUDTETR: PLUS GRAND NUMERO DANS NOTETR DE TETRAEDRE OCCUPE
C NOTE1S : TABLEAU DES FACES CONTENANT UN SOMMET
C NBFAPE : NOMBRE DE FACES PERDUES DE LEFACO DU CONTOUR AVANT ET APRES
C          NON PRESENTES DANS LA TETRAEDRISATION
C NOFAPE : NUMERO DANS LEFACO DE LA FACE PERDUE

C SORTIES:
C --------
C NBSOMM  : NUMERO DU DERNIER MILIEU AJOUTE AU MILIEU D'UNE ARETE PERDUE
C NBARRECU: NOMBRE D'ARETES RETROUVEES PAR CHANGMENTS 2T=>3T
C IERR    : 0 SI PAS D'ERREUR
C           1 SI RETOUR EN ARRIERE A PROGRAMMER
C           2 SI SATURATION DU TABLEAU NOTETR
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC       JUIN 1993
C MODIFS : ALAIN PERRONNET  Saint Pierre du Perray         Decembre 2019
C....................................................................012
      include"./incl/gsmenu.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      DOUBLE PRECISION  PTXYZD(1:4,1:MXSOMM)
      INTEGER           NPSOFR(1:MXSOMM),
     %                  NOTETR(8,*),
     %                  N1TETS(*),
     %                  LEFACO(11,0:MXFACO),
     %                  N1FASC(1:MXSOMM),
     %                  NOFAPE(1:NBFAPE),
     %                  NOTE1S(MXTE1S),
     %                  NVOLTE(*)
      DOUBLE PRECISION  D, DIS2ST, LONARE(3), DARMAX
      INTEGER           NTNOUV(3),
     %                  NOSOTR(3),
     %                  NOSOFATE(3,4)
C     NO DES SOMMETS DES FACES POUR QUE VU DU SOMMET MANQUANT
C     LES SOMMETS SOIENT VUS DE L'INTERIEUR DANS LE SENS DIRECT
      DATA              NOSOFATE/ 1,2,3, 2,4,3, 3,4,1, 4,2,1 /

C     LA BOUCLE SUR LES FACES PERDUES
C     ===============================
      NBSOMM0  = NBSOMM
      NBARRECU = 0
      DO 1000 NFP=1,NBFAPE

C        LE NUMERO LEFACO DE LA FACE PERDUE
         NFPE = NOFAPE(NFP)
         IF( NFPE .LE. 0 ) GOTO 1000

C        BOUCLE SUR LES 3 ARETES DE LA FACE PERDUE NFPE DE LEFACO
C        --------------------------------------------------------
         NBARRE = 0
         DARMAX = -1D0
         DO 900 I=1,3

C           L'ARETE I DE LA FACE PERDUE NFPE
            IF( I .EQ. 3 ) THEN
               I1 = 1
            ELSE
               I1 = I+1
            ENDIF
C           LES NUMEROS DES 2 SOMMETS DE L'ARETE I DE LA FACE PERDUE NFPE
            NS1 = LEFACO(I ,NFPE)
            NS2 = LEFACO(I1,NFPE)

            D = DIS2ST( PTXYZD( 1, NS1 ), PTXYZD( 1, NS2 ) )
            LONARE(I) = D
            IF( D .GT. DARMAX ) THEN
               DARMAX = D
               NARMAX = I
            ENDIF

C           CETTE ARETE I EST ELLE DANS LA TETRAEDRISATION?
            CALL TETR1A( NS1,    NS2,    N1TETS, NOTETR,
     %                   NBTE1A, MXTE1S, NOTE1S, IER )
            IF( NBTE1A .LE. 0 ) THEN

C              L'ARETE N'EST PAS DANS LA TETRAEDRISATION
C              RECHERCHE DES NBF1 TETRAEDRES DE SOMMET NS1
               CALL TETR1S( NS1,  N1TETS, NOTETR,
     %                      NBF1, MXTE1S, NOTE1S, IERR )
               IF( NBF1 .LE. 0 ) GOTO 900

C              PARCOURS DES TETRAEDRES DE SOMMET NS1
C              -------------------------------------
               DO 90 JJ=1,NBF1

C                 LE TETRAEDRE DE SOMMET NS1
                  NT = NOTE1S( JJ )
                  IF( NOTETR(1,NT) .LE. 0 ) GOTO 90
C
C                 RECHERCHE DU TETRAEDRE NTOP OPPOSE A LA FACE OPPOSEE
C                 AU SOMMET NS1 DE NT: RECHERCHE DU NO DE SOMMET NS1 DANS NT
                  DO K=1,4
                     IF( NS1 .EQ. NOTETR(K,NT) ) GOTO 30
                  ENDDO
                  GOTO 90
C
C                 LA FACE OPPOSEE AU SOMMET K EST K MOD 4 + 1
 30               IF( K .NE. 4 ) THEN
                     NF=K+1
                  ELSE
                     NF=1
                  ENDIF
C
C                 NTOP LE TETRAEDRE OPPOSE A NS1 ET LE NUMERO NFOP DE LA FACE
                  CALL NOFAOP( NF, NT, NOTETR, NFOP, NTOP )
C                 LE SOMMET OPPOSE NSOP EST NFOP - 1
                  IF( NFOP .NE. 1 ) THEN
                     NSOP = NFOP - 1
                  ELSE
                     NSOP = 4
                  ENDIF
                  IF( NOTETR(NSOP,NTOP) .EQ. NS2 ) THEN
C
C                    ICI L'ARETE PERDUE NS1 NS2 EST LA "DIAGONALE"
C                    DES 2 TETRAEDRES NT et NTOP : 2T -> 3T?
C                    ===========================================================
C                    L'ARETE NS1-NS2 PERFORE-T-ELLE LA FACE COMMUNE AUX
C                    2 TETRAEDRES ALORS ECHANGER 2 TETRAEDRES => 3 TETRAEDRES
C                    D'ARETE COMMUNE NS1-NS2 QUI APPARTIENT A LA TETRAEDRISATION
                     CALL CH2T3T( INFACO, MXFACO, LEFACO, 1,
     %                            IVOLTE, NVOLTE, NT,    NF,    NTOP,
     %                            PTXYZD, N1TETS, NOTETR,N1TEVI,NUDTETR,
     %                            NTNOUV, IERR )
C                    IERR=1 SI PAS DE TETRAEDRE DE L'AUTRE COTE DE LA FACE
C                         2 SI SATURATION DU TABLEAU NOTETR
C                         3 SI LE VOLUME de NT+NTOP EST DIFFERENT DU VOLUME
C                              DES 3 AUTRES TETRAEDRES
                     IF( IERR .EQ. 2 ) THEN
C                       SATURATION DES TETRAEDRES
                    print*,'recuarpe: SATURATION DES TETRAEDRES N1TEVI='
     %                    ,N1TEVI
                        GOTO 9000
                     ENDIF
                     IF( IERR .NE. 0 ) THEN
                        IERR = 0
                        GOTO 90
                     ELSE
C                       2T=>3T A ETE CORRECTEMENT EXECUTE
C                       L'ARETE NS1-NS2 A ETE CREEE DANS LA TETRAEDRISATION
                        NBARRE = NBARRE + 1
                        GOTO 900
                     ENDIF
                  ENDIF


C                 NS2 N'EST PAS LE SOMMET OPPOSE A LA FACE COMMUNE NF DE NT
C                 RECHERCHE D'UNE SUITE DE TETRAEDRES DONT 2 FACES SONT
C                 PERFOREES PAR L'ARETE NS1-NS2 PERDUE ET TELS QUE
C                 2T => 3T DONNENT LE MEME VOLUME
C                 =========================================================
                  NT0 = NT
                  NF0 = NF
C                 LA FACE EST ELLE PERFOREE PAR L'ARETE NS1-NS2?
C                 NUMERO DES 3 SOMMETS DE LA FACE COMMUNE
                  NOSOTR(1) = NOTETR( NOSOFATE(1,NF), NT )
                  NOSOTR(2) = NOTETR( NOSOFATE(2,NF), NT )
                  NOSOTR(3) = NOTETR( NOSOFATE(3,NF), NT )
C                 INTERSECTION OU NON DE NS1-NS2 AVEC LE TRIANGLE NOSOTR?
                  CALL VOLU2T3T( NS1, NS2, NOSOTR, PTXYZD, NONOUI )
                  IF( NONOUI .NE. 1 ) GOTO 90
C
C                 OUI: NS1-NS2 INTERSECTE LE TRIANGLE NOSOTR
C                      L'ARETE NS1-NOTETR(NSOP,NTOP) PERFORE-T-ELLE
C                      LE TRIANGLE NOSOTR?
 48               CALL VOLU2T3T( NS1, NOTETR(NSOP,NTOP), NOSOTR, PTXYZD,
     %                           NONOUI )
                  IF( NONOUI .NE. 1 ) GOTO 90
C
C                 OUI: RECHERCHE DE LA FACE OPPOSEE DE NTOP PERFOREE PAR NS1-NS2
                  DO 50 J=1,4

                     IF( J .EQ. NFOP ) GOTO 50

C                    LA FACE J DU TETRAEDRE NTOP. SES SOMMETS
                     NOSOTR(1) = NOTETR( NOSOFATE(1,J), NTOP )
                     NOSOTR(2) = NOTETR( NOSOFATE(2,J), NTOP )
                     NOSOTR(3) = NOTETR( NOSOFATE(3,J), NTOP )
C
C                    INTERSECTION OU NON DE NS1-NS2 AVEC LE TRIANGLE NOSOTR?
                     CALL VOLU2T3T( NS1, NS2, NOSOTR, PTXYZD, NONOUI )
                     IF( NONOUI .NE. 1 ) GOTO 50
C
C                    OUI: INTERSECTION NS1-NS2 AVEC LA FACE J
C                    PASSAGE AU TETRAEDRE SUIVANT
                     NT = NTOP
                     NF = J
C
C                    RECHERCHE DU TETRAEDRE OPPOSE
                     CALL NOFAOP( NF, NT, NOTETR, NFOP, NTOP )
                     IF( NFOP .LE. 0 ) THEN
C                       LE TEMPS DE LA MISE AU POINT...
                        PRINT*,'recuarpe 60: PB PAS DE TETRAEDRE OPPOSE 
     %a la FACE',NF,' du TETRAEDRE',NT
                        GOTO 90
                     ENDIF
C
C                    LE SOMMET OPPOSE DANS NTOP EST IL NS2?
                     DO JK=1,4
                        IF( NOTETR(JK,NTOP) .EQ. NS2 ) THEN
C                          OUI: NS2 EST OPPOSE A LA FACE NFOP de NTOP
C                          LE TEMPS DE LA MISE AU POINT...
                           IF( (MOD(JK,4)+1) .NE. NFOP ) THEN
                              PRINT*,'recuarpe: SOMMET NS2=',NS2,
     %                        ' NON OPPOSE A LA FACE',JK
                              GOTO 90
                           ENDIF
C                          NS2 EST OPPOSE A LA FACE NFOP et
C                          C'EST LE DERNIER TETRAEDRE LE LONG DE NS1-NS2
                           GOTO 100
                        ENDIF
                     ENDDO
C
C                    NON : LE SOMMET DE NTOP OPPOSE A LA FACE NFOP
C                    LE SOMMET OPPOSE NSOP EST NFOP - 1
                     IF( NFOP .NE. 1 ) THEN
                        NSOP = NFOP - 1
                     ELSE
                        NSOP = 4
                     ENDIF
                     GOTO 48

 50               ENDDO

 90            ENDDO
               GOTO 900


C              CHANGEMENT DES TETRAEDRES POUR FAIRE APPARAITRE L'ARETE NS1-NS2
C              REDEPART A PARTIR DE LA PREMIERE FACE INTERSECTEE PAR NS1-NS2
C              ===============================================================
 100           NT = NT0
               NF = NF0

C              L'ARETE NS1-NS2 PERFORE LA FACE COMMUNE AUX 2 TETRAEDRES
C              2 TETRAEDRES => 3 TETRAEDRES D'ARETE COMMUNE NS1-NS2
               CALL CH2T3T( INFACO, MXFACO, LEFACO, 1,
     %                      IVOLTE, NVOLTE, NT,     NF,     NTOP,
     %                      PTXYZD, N1TETS, NOTETR, N1TEVI, NUDTETR,
     %                      NTNOUV, IERR )
C              IERR = 1 SI PAS DE TETRAEDRE DE L'AUTRE COTE DE LA FACE
C                     2 SI SATURATION DU TABLEAU NOTETR
C                     3 SI LE VOLUME de NT+NTOP EST DIFFERENT DU VOLUME
C                          DES 3 AUTRES TETRAEDRES
               IF( IERR .EQ. 2 ) THEN
C                 SATURATION DES TETRAEDRES
                  GOTO 9000
               ENDIF
               IF( IERR .NE. 0 ) THEN
                  IF( NT .EQ. NT0 ) THEN
                     IERR = 0
                     GOTO 900
                  ELSE
C                    PROBLEME: 2T=>3T A DEJA ETE FAIT ET TILT AVANT NS2
                     PRINT*,'recuarpe: RETOUR EN ARRIERE A FAIRE'
                     IERR = 1
                     GOTO 9000
                  ENDIF
               ENDIF

C              ICI: 2T=NT+NTOP => 3T=NTNOUV(1:3) A ETE CORRECTEMENT EXECUTE
C              ------------------------------------------------------------
C              RECHERCHE PARMI LES 3 TETRAEDRES NTNOUV GENERES DU TETRAEDRE
C              DONT UNE FACE EST PERFOREE PAR L'ARETE NS1-NS2
               NF0 = 0
  120          DO 150 J=1,3

C                 LE NOUVEAU TETRAEDRE J
                  NT = NTNOUV(J)

                  DO 140 NF=1,4

C                    NUMERO DES 3 SOMMETS DE LA FACE NF DE NT
                     NOSOTR(1) = NOTETR( NOSOFATE(1,NF), NT )
                     NOSOTR(2) = NOTETR( NOSOFATE(2,NF), NT )
                     NOSOTR(3) = NOTETR( NOSOFATE(3,NF), NT )
                     CALL VOLU2T3T( NS1, NS2, NOSOTR, PTXYZD, NONOUI )
                     IF( NONOUI .NE. 1 ) GOTO 140
C
C                    LE TETRAEDRE OPPOSE A LA FACE NF DU TETRAEDRE NT
                     CALL NOFAOP( NF, NT, NOTETR, NFOP, NTOP )
C                    LE SOMMET OPPOSE NSOP EST NFOP - 1
                     IF( NFOP .NE. 1 ) THEN
                        NSOP = NFOP - 1
                     ELSE
                        NSOP = 4
                     ENDIF
                     IF( NOTETR(NSOP,NTOP) .EQ. NS2 ) NF0 = 1
C
C                    L'ARETE NS1-NS2 INTERSECTE LA FACE NF DE NT
C                    2 TETRAEDRES => 3 TETRAEDRES POSSIBLE?
                     CALL CH2T3T( INFACO, MXFACO, LEFACO, 1,
     %                            IVOLTE, NVOLTE, NT, NF, NTOP,
     %                          PTXYZD, N1TETS, NOTETR, N1TEVI, NUDTETR,
     %                            NTNOUV, IERR )
                     IF( IERR .EQ. 2 ) THEN
C                       SATURATION DES TETRAEDRES
                        GOTO 9000
                     ENDIF
                     IF( IERR .NE. 0 ) THEN
                        IF( NT .EQ. NT0 ) THEN
                           IERR = 0
                           GOTO 900
                        ELSE
C                          PROBLEME: 2T=>3T A DEJA ETE FAIT ET TILT AVANT NS2
                       PRINT*,'recuarpe: RETOUR EN ARRIERE A PROGRAMMER'
                           IERR = 1
                           GOTO 9000
                        ENDIF
                     ENDIF

C                    NTOP EST IL LE DERNIER TETRAEDRE?
                     IF( NF0 .NE. 0 ) THEN
C                       OUI: SORTIE POUR CETTE ARETE ENTIEREMENT RETROUVEE
                        NBARRE = NBARRE + 1
                        GOTO 900
                     ELSE
C                       NON: PASSAGE AU TETRAEDRE SUIVANT
                        GOTO 120
                     ENDIF

 140              ENDDO

 150           ENDDO
               PRINT*,'recuarpe: ANOMALIE PAS DE TETRAEDRE PERCE'

            ENDIF

 900     ENDDO

         IF( NBARRE .GT. 0 ) GOTO 990


C        PAS D'ARETE PERDUE RECUPEREE
C        AJOUT DU POINT MILIEU DE LA PLUS GRANDE DES 3 ARETES DE LA FACE NFPE
C        --------------------------------------------------------------------
ccc         CALL VD1F2F(      1, NARMAX, NFPE,   MXFACO, LEFACO, NBFACO,
         CALL VD1F2F(      0, NARMAX, NFPE,   MXFACO, LEFACO, NBFACO,
     %                N1FASC, NBFAPE, NOFAPE,
     %                NBSOMM, MXSOMM, PTXYZD, NPSOFR,
     %                MXTE1S, MNTE1S, IERR )
         IF( IERR .NE. 0 ) GOTO 1000
  

C        NOMBRE TOTAL D'ARETES RETROUVEES
 990     NBARRECU = NBARRECU + NBARRE

1000  ENDDO

C     TRACE DES FACES PERDUES
 9000 CALL TRFAPE( NBFAPE, NOFAPE, MXFACO, LEFACO, PTXYZD )
      PRINT*,'recuarpe:',NBSOMM-NBSOMM0,
     %       ' MILIEUX AJOUTES SUR LES ARETES PERDUES'

      RETURN
      END
