      SUBROUTINE IDSTST( NSNOFR, NSDLFR, MXSOMM, PTXYZD, NPSOFR,
     %                   MXTETR, N1TEVI, N1TETS, NOTETR, NUDTETR,
     %                   IVOLTE, NVOLTE,
     %                   MXFETO, N1FEOC, NFETOI,
     %                   MXTENS, NBTENS, NOTENS, QUAMIN, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : REMPLACER DANS LA TETRAEDRISATION LE SOMMET NSNOFR NON FRONTALIER
C ----- PAR LE SOMMET FRONTALIER NSDLFR EN SUPPRIMANT LES TETRAEDRES
C       D'ARETE CES 2 SOMMETS

C ENTREES:
C --------
C NSNOFR : NUMERO DU SOMMET PTXYZD A IDENTIFIER AU SOMMET NSDLFR
C NSDLFR : NUMERO DU SOMMET PTXYZD DE LA FRONTIERE
C MXSOMM : NOMBRE DE SOMMETS DECLARABLES DANS PTXYZD
C PTXYZD : TABLEAU DES COORDONNEES DES POINTS
C          PAR POINT : X  Y  Z DISTANCE_SOUHAITEE
C NPSOFR : =  0 SI POINT AJOUTE NON SUR LA SURFACE DE L'OBJET
C          =  1 SI LE POINT EST AJOUTE SUR UNE ARETE DE LA SURFACE
C          =  2 SI LE POINT EST AJOUTE COMME BARYCENTRE SUR LA SURFACE
C          =  1 000 000  + NUMERO DU  POINT INTERNE UTILISATEUR
C          = (1 000 000) * (NO SURFACE + 1) + NUMERO DU SOMMET DANS
C              LA NUMEROTATION DES SOMMETS DE LA SURFACE
C          = -4 SI LE POINT EST SOMMET NON TROP PROCHE PT OU FACE
C          = -1 SI LE POINT EST SOMMET TROP PROCHE PT OU FACE
C          = -3 SI LE POINT EST SOMMET REFUSE DANS LA TETRAEDRISATION
C          = -NPSOFR(I) SI POINT DEPLACE SUR LA SURFACE
C            DANS LEUR SURFACE FERMEE INITIALE ou NO DE POINT INTERNE
C MXTETR : NOMBRE DE TETRAEDRES DECLARABLES DANS NOTETR
C MXTENS : NOMBRE MAXIMAL DE TETRAEDRES DECLARABLES DANS NOTENS
C NBTENS : NOMBRE DE  TETRAEDRES DE NOTETR CONTENANT LE SOMMET NSDLFR
C          ou LE SOMMET NSNOFR ou OPPOSE AUX FACES DES TETRAEDRES
C          DE SOMMET NSNOFR
C NOTENS : NUMERO DES TETRAEDRES DE NOTETR CONTENANT LE SOMMET NSDLFR
C          ou LE SOMMET NSNOFR ou OPPOSE AUX FACES DES TETRAEDRES
C          DE SOMMET NSNOFR
C MXFETO : NOMBRE MAXIMAL DE FACES DECLARABLES DANS LE CHAINAGE NFETOI
C NUDTETR: PLUS GRAND NUMERO DANS NOTETR DE TETRAEDRE OCCUPE
C IVOLTE : =0 TABLEAU NVOLTE NON PRESENT
C          =1 TABLEAU NVOLTE     PRESENT

C MODIFIES:
C ---------
C NVOLTE : NUMERO DE VOLUME (1 a NBVOPA) DE CHAQUE TETRAEDRE
C          -1 SI VOLUME INCONNU (Exemple: les TETRAEDRES VIDES DE NOTETR)
C N1FEOC : POINTEUR SUR LA PREMIERE FACE DE L'ETOILE
C          CHAINAGE SUIVANT DANS NFETOI(5,*)
C N1TEVI : NO DU PREMIER TETRAEDRE VIDE DANS NOTETR
C N1TETS : N1TETS(NS) NUMERO D'UN TETRAEDRE AYANT POUR SOMMET NS
C NOTETR : LISTE DES TETRAEDRES
C          SOMMET1,    SOMMET2,    SOMMET3,    SOMMET4,
C          TETRAEDRE1, TETRAEDRE2, TETRAEDRE3, TETRAEDRE4
C          DE L'AUTRE COTE DE LA FACE
C          1: 123      2: 234      3: 341      4: 412
C MXFETO : NOMBRE MAXIMAL DE FACES DECLARABLES DANS LE CHAINAGE NFETOI
C NFETOI : VERSION2 ETOILE DES TETRAEDRES DES SOMMETS NSNOFR et NSDLFR
C          1: NUMERO DU TETRAEDRE DANS NOTETR OPPOSE A CETTE FACE
C          2: NUMERO PTXYZD DU SOMMET 1 DE LA FACE
C          3: NUMERO PTXYZD DU SOMMET 2 DE LA FACE
C          4: NUMERO PTXYZD DU SOMMET 3 DE LA FACE
C             S1S2xS1S3 EST DIRIGE VERS L'INTERIEUR DE L'ETOILE
C          5: CHAINAGE SUIVANT DES FACES OCCUPEES ET VIDES

C SORTIES :
C ---------
C NBTENS : NOMBRE DE  TETRAEDRES DE NOTETR CONTENANT LE SOMMET NSDLFR
C NOTENS : NUMERO DES TETRAEDRES DE NOTETR CONTENANT LE SOMMET NSDLFR
C QUAMIN : SI NBTENS>0 ALORS QUALITE MIN DES NBTENS TETRAEDRES NOTENS
C IERR   : =0 PAS D'ERREUR RENCONTREE ET IDENTIFICATION FAITE
C          >0 PAS D'IDENTIFICATION FAITE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY Novembre 2016
C2345X7..............................................................012
      COMMON / TRTETR / STOPTE, TRACTE
      LOGICAL           STOPTE, TRACTE
ccc      LOGICAL           TRACTE0

      DOUBLE PRECISION  PTXYZD(4,MXSOMM)
      INTEGER           NPSOFR(MXSOMM), NOTETR(8,MXTETR),N1TETS(MXSOMM),
     %                  NOTENS(MXTENS), NFETOI(5,MXFETO),NVOLTE(MXTETR)

      CHARACTER*128     KTITRE
      DOUBLE PRECISION  VOLUT0, VOLUTE, VOLUE0, VOLUE1,
     %                  ARMIN, ARMAX, SURFTR(4)
      INTEGER           NOSOTR(3), NOSOTE(4)

      INTEGER           NOSOFA(3,4)
      DATA              NOSOFA   / 1,3,2,  2,3,4,  3,1,4,  4,1,2 /
C     => LA NORMALE A LA FACE EST ORIENTEE VERS L'EXTERIEUR DU TETRAEDRE

      INTEGER           NOFA2S(1:2, 1:3, 2:4)
      DATA              NOFA2S / 1,4, 0,0, 0,0,
     %                           1,3, 1,2, 0,0,
     %                           3,4, 2,4, 2,3 /
C     => NUMERO DES 2 FACES DES 2 SOMMETS NS1-NS2

C     REMPLACEMENT DU NO DE SOMMET NSNOFR PAR NSDLFR
ccc      PRINT*
ccc      PRINT*,'idstst: Debut d''IDENTIFICATION du sommet',NSNOFR,
ccc     %       ' par le sommet',NSDLFR

      IF( NPSOFR(NSNOFR) .GT. 0 ) THEN
         IERR = 5
         GOTO 9999
      ENDIF

ccc      TRACTE0 = TRACTE

C     REPERAGE DES TETRAEDRES D'ARETE NSNOFR-NSDLFR
C     ---------------------------------------------
      QUAMIN = 2.
      VOLUE1 = 0D0
      VOLUE0 = 0D0
      DO 10 NT=1,NBTENS

C        TETRAEDRE DE SOMMET NSNOFR ou NSDLFR
         NTE = NOTENS( NT )

C        CE TETRAEDRE EST IL ACTIF?
         IF( NTE .LE. 0 ) GOTO 10
         IF( NOTETR(1,NTE) .EQ. 0 ) GOTO 10

C        OUI: SON VOLUME INITIAL ET SA QUALITE INITIALE
         CALL QUATETD( PTXYZD( 1, NOTETR(1,NTE) ),
     %                 PTXYZD( 1, NOTETR(2,NTE) ),
     %                 PTXYZD( 1, NOTETR(3,NTE) ),
     %                 PTXYZD( 1, NOTETR(4,NTE) ),
     %                 ARMIN, ARMAX, SURFTR, VOLUT0, QUALTE )
         VOLUE0 = VOLUE0 + VOLUT0

ccc         PRINT *,'idstst: ',NT,' initial NOTETR(',NTE,')=',
ccc     %           (NOTETR(M,NTE),M=1,8),' Q=',QUALTE,' V=',VOLUT0

         DO N0=1,4

            IF( NOTETR(N0,NTE) .EQ. NSNOFR ) THEN

C              NTE DE SOMMET NSNOFR A T IL AUSSI LE SOMMET NSDLFR?
               DO N1=1,4

                  IF( NOTETR(N1,NTE) .EQ. NSDLFR ) THEN

C                    NTE TETRAEDRE D'ARETE NSNOFR-NSDLFR A SUPPRIMER
C                    -----------------------------------------------
ccc                     PRINT *,'idstst: NSDLFR=',NSDLFR,' NT=',NT,
ccc     %                       ' SUPPRIME NOTETR(',NTE,')=',
ccc     %              (NOTETR(M,NTE),M=1,8),' Q=',QUALTE,' V=',VOLUT0


C                    BOUCLE CI DESSOUS AJOUTEE. A VERIFIER  13/1/2019
C                    SES 4 TETRAEDRES OPPOSES ONT LE TETRAEDRE NTE
C                    OPPOSE QUI DEVIENT INCONNU
                     DO M1=1,4
                        CALL NOFAOP( M1, NTE, NOTETR, NFOP, NTOP )
                        IF( NFOP .GT. 0 ) THEN
                           NOTETR( 4+NFOP, NTOP ) = -1
                        ENDIF
                     ENDDO

C                    LES 2 FACES DE L'ARETE NSNOFR-NSDLFR ONT
C                    DES TETRAEDRES OPPOSES INCONNUS
                     IF( N0 .GT. N1 ) THEN
                        N3 = N1
                        N4 = N0
                     ELSE
                        N3 = N0
                        N4 = N1
                     ENDIF

C                    LE TETRAEDRE OPPOSE A LA FACE NF1 DE NTE
                     NF1    = NOFA2S(1,N3,N4)
                     NTEOP1 = NOTETR(4+NF1,NTE)

C                    LE TETRAEDRE OPPOSE A LA FACE NF2 DE NTE
                     NF2    = NOFA2S(2,N3,N4)
                     NTEOP2 = NOTETR(4+NF2,NTE)

                     IF( NTEOP1.GT.0 .AND. NOTETR(1,NTEOP1).GT.0 ) THEN
                        NOSOTR(1) = NOTETR( NOSOFA(1,NF1), NTE )
                        NOSOTR(2) = NOTETR( NOSOFA(2,NF1), NTE )
                        NOSOTR(3) = NOTETR( NOSOFA(3,NF1), NTE )
                        CALL NUFATRTE( NOSOTR, NOTETR(1,NTEOP1), M1 )
                        IF( M1 .GT. 0 ) THEN
                           IF(NTEOP2.GT.0.AND.NOTETR(1,NTEOP2).GT.0)THEN
                              NOSOTR(1) = NOTETR( NOSOFA(1,NF2), NTE )
                              NOSOTR(2) = NOTETR( NOSOFA(2,NF2), NTE )
                              NOSOTR(3) = NOTETR( NOSOFA(3,NF2), NTE )
                              CALL NUFATRTE(NOSOTR,NOTETR(1,NTEOP2),M2)
                              IF( M2 .GT. 0 ) THEN
                                 NOTETR(4+M1,NTEOP1) = NTEOP2
                                 NOTETR(4+M2,NTEOP2) = NTEOP1
                              ENDIF
                           ENDIF
                        ENDIF
                     ENDIF

C                    LE TETRAEDRE NTE EST MARQUE A SUPPRIMER
                     NOTENS( NT ) = -NTE
                     GOTO 10

                  ENDIF

               ENDDO

C              NTE TETRAEDRE AVEC NSNOFR MAIS PAS NSDLFR
C              DANS NTE LE SOMMET NSNOFR -> NSDLFR
C              -----------------------------------------
               DO M=1,4
                  NOSOTE( M ) = NOTETR( M, NTE )
               ENDDO
               NOSOTE( N0 ) = NSDLFR

C              SON VOLUME ET SA QUALITE
               CALL QUATETD( PTXYZD( 1, NOSOTE(1) ),
     %                       PTXYZD( 1, NOSOTE(2) ),
     %                       PTXYZD( 1, NOSOTE(3) ),
     %                       PTXYZD( 1, NOSOTE(4) ),
     %                       ARMIN, ARMAX, SURFTR, VOLUTE, QUALTE )
               VOLUE1 = VOLUE1 + VOLUTE

ccc                  PRINT *,'idstst: ',NT,' SUPPV<0 NOTETR(',NTE,')=',
ccc     %                     NOSOTE,(NOTETR(M,NTE),M=5,8),
ccc     %                    ' Q1=',QUALTE,' V1=',VOLUTE

               GOTO 10

            ENDIF

         ENDDO

C        ICI: NSNOFR N'EST PAS UN SOMMET DE NTE
C        --------------------------------------
         VOLUE1 = VOLUE1 + VOLUT0

 10   ENDDO

      KTITRE='idstst1:         TETRAEDRES de SOMMET            '
      WRITE(KTITRE(10:16),'(I7)') NBTENS
      WRITE(KTITRE(40:47),'(I8)') NSDLFR
      CALL TRACEETOILE( KTITRE, PTXYZD, NOTETR, NSDLFR, NBTENS, NOTENS )
ccc      CALL TRFETO13( KTITRE, PTXYZD, NBTENS, NOTENS, NOTETR )

C     LES TETRAEDRES DE SOMMETS NSNOFR et NSDLFR SONT SUPPRIMES
C     ET LE SOMMET NSNOFR DEVIENT NSDLFR
C     COMPRESSION DES TETRAEDRES DE SOMMET NSNOFR ou NSDLFR RESTANTS
C     --------------------------------------------------------------
      VOLUE1 = 0D0
      NBT    = 0
      DO 15 NT = 1, NBTENS

C        TETRAEDRE DE SOMMET NSNOFR=NSDLFR
         NTE = NOTENS( NT )

         IF( NTE .LT. 0 ) THEN

C           NTE A SUPPRIMER DEVIENT LE PREMIER TETRAEDRE VIDE
C           -------------------------------------------------
            NTE = -NTE

ccc            PRINT *,'idstst: ',NT,' SUPPRim NOTETR(',NTE,')=',
ccc     %               (NOTETR(M,NTE),M=1,8)

ccc            DO NF = 1,4

cccC              LE TETRAEDRE OPPOSE A LA FACE NF DE NTE
ccc               NTEOP = NOTETR( 4+NF, NTE )

ccc               IF( NTEOP.GT.0 .AND. NOTETR(1,NTEOP).GT.0 ) THEN

cccC                 LES 3 SOMMETS DE LA FACE NF DE NTE
ccc                  NOSOTR(1) = NOTETR( NOSOFA(1,NF), NTE )
ccc                  NOSOTR(2) = NOTETR( NOSOFA(2,NF), NTE )
ccc                  NOSOTR(3) = NOTETR( NOSOFA(3,NF), NTE )
ccc                  CALL NUFATRTE( NOSOTR, NOTETR(1,NTEOP), NFOP )

ccc                  IF( NFOP .GT. 0 ) THEN
cccC                    LE TETRAEDRE NTE SUPPRIME DEVIENT INCONNU
ccc                     NOTETR( 4+NFOP, NTEOP ) = -1
ccc                  ENDIF

ccc               ENDIF

ccc            ENDDO

C           NTE EST UN TETRAEDRE VIDE
            NOTETR( 1, NTE ) = 0
            NOTETR( 2, NTE ) = 0
            NOTETR( 3, NTE ) = 0
            NOTETR( 4, NTE ) = 0
            NOTETR( 5, NTE ) = N1TEVI
            NOTETR( 6, NTE ) = 0
            NOTETR( 7, NTE ) = 0
            NOTETR( 8, NTE ) = 0

C           NUMERO DE VOLUME INCONNU
            IF( IVOLTE .NE. 0 ) NVOLTE( NTE ) = -1

C           NTE DEVIENT LE PREMIER TETRAEDRE VIDE
            N1TEVI = NTE

            GOTO 15

         ENDIF

C        NTE A CONSERVER A T IL POUR SOMMET NSNOFR? => NSDLFR
C        ----------------------------------------------------
         DO N0=1,4

C           LE SOMMET N0 DU TETRAEDRE NTE
            NST = NOTETR( N0, NTE )
            IF( NST .EQ. NSNOFR ) THEN

C              OUI: DANS NTE NSNOFR -> NSDLFR
               NOTETR( N0, NTE ) = NSDLFR
               N1TETS( NSDLFR ) = NTE

            ENDIF

         ENDDO

C        COMPRESSION DES TETRAEDRES VIDES
         NBT = NBT + 1
         NOTENS( NBT ) = NTE

C        VOLUME DE L'ETOILE DES NOUVEAUX TETRAEDRES
C        SON VOLUME ET SA QUALITE
         CALL QUATETD( PTXYZD( 1, NOTETR(1,NTE) ),
     %                 PTXYZD( 1, NOTETR(2,NTE) ),
     %                 PTXYZD( 1, NOTETR(3,NTE) ),
     %                 PTXYZD( 1, NOTETR(4,NTE) ),
     %                 ARMIN, ARMAX, SURFTR, VOLUTE, QUALTE )
         VOLUE1 = VOLUE1 + VOLUTE

ccc         PRINT *,'idstst: ',NT,' SUPPV<0 NOTETR(',NTE,')=',
ccc     %            (NOTETR(M,NTE),M=1,8),
ccc     %           ' Q1=',QUALTE,' V1=',VOLUTE

 15   ENDDO
      NBTENS = NBT

C     MISE A JOUR DU NUMERO DE TETRAEDRE DU SOMMET NSNOFR
      N1TETS( NSNOFR ) = 0

ccc      PRINT *,'idstst: VolET0=',VOLUE0,' VolET1=',VOLUE1,
ccc     %        ' V1/V0=',VOLUE1/VOLUE0,
ccc     %        ' apres IDENTIFICATION de ',NSNOFR,'->',NSDLFR

C     MISE A JOUR DES TETRAEDRES OPPOSES DES FACES
C     DE LEURS TETRAEDRES OPPOSES
C     --------------------------------------------
      CALL MJOPTE( NBTENS, NOTENS, N1TETS, NOTETR, NUDTETR,
     %             N1TEVI, PTXYZD, NBFANR )

      KTITRE='idstst2:          TETRAEDRES de SOMMET            AVANT VE
     %RIFICATION TETRAEDRES MULTIPLES'
      WRITE(KTITRE(10:16),'(I7)') NBT
      WRITE(KTITRE(40:47),'(I8)') NSDLFR
      CALL TRACEETOILE( KTITRE, PTXYZD, NOTETR, NSDLFR, NBTENS, NOTENS )


C     VERIFICATION PAS DE TETRAEDRES DOUBLES?
C     SI OUI: ILS SONT TOUS LES 2 SUPPRIMES
C     NSNOFR A TRAVERSE LA FACE OPPOSEE JUSQU'A NSDLFR
C     ET LES 2 TETRAEDRES ONT DISPARUS
C     ------------------------------------------------
      NBTEDB = 0
      DO NT = 1, NBTENS

         NTE = NOTENS( NT )

         DO 34 MT=NT+1,NBTENS

            NTE2 = NOTENS( MT )
C           NTE2 EST IL UN TETRAEDRE IDENTIQUE A NTE?
            IF( NTE2 .GT. 0 . AND. NOTETR(1,NTE2).GT.0 .AND.
     %          NOTETR(1,NTE) +NOTETR(2,NTE)
     %         +NOTETR(3,NTE) +NOTETR(4,NTE) .EQ.
     %          NOTETR(1,NTE2)+NOTETR(2,NTE2)
     %         +NOTETR(3,NTE2)+NOTETR(4,NTE2) ) THEN

C              TETRAEDRES NTE et NTE2 SONT ILS VRAIMENT IDENTIQUES?
               DO 30 N=1,4
                  NS1 = NOTETR(N,NTE)
                  DO M=1,4
                     NS2 = NOTETR(M,NTE2)
                     IF( NS1 .EQ. NS2 ) GOTO 30
                  ENDDO
C                 NON: LES TETRAEDRES NTE et NTE2 SONT DIFFERENTS
                  GOTO 34
 30            ENDDO

C              LES NUMEROS DES 4 SOMMETS SONT IDENTIQUES
C              => LES TETRAEDRES NTE et NTE2 IDENTIQUES
C              -----------------------------------------
               NBTEDB = NBTEDB + 1
ccc               PRINT *
ccc               PRINT *,'idstst5: 2 TETRAEDRES IDENTIQUES',NTE,NTE2
C              NTE: SON VOLUME ET SA QUALITE
               CALL QUATETD( PTXYZD( 1, NOTETR(1,NTE) ),
     %                       PTXYZD( 1, NOTETR(2,NTE) ),
     %                       PTXYZD( 1, NOTETR(3,NTE) ),
     %                       PTXYZD( 1, NOTETR(4,NTE) ),
     %                       ARMIN, ARMAX, SURFTR, VOLUTE, QUALTE )
ccc               PRINT *,'idstst5: TetraId1 NOTETR(',NTE,')=',
ccc     %                 (NOTETR(M,NTE),M=1,8),
ccc     %                 ' Q=',QUALTE,' V=',VOLUTE

C              NTE2: SON VOLUME ET SA QUALITE
               CALL QUATETD( PTXYZD( 1, NOTETR(1,NTE2) ),
     %                       PTXYZD( 1, NOTETR(2,NTE2) ),
     %                       PTXYZD( 1, NOTETR(3,NTE2) ),
     %                       PTXYZD( 1, NOTETR(4,NTE2) ),
     %                       ARMIN, ARMAX, SURFTR, VOLUTE, QUALTE )
ccc               PRINT *,'idstst5: TetraId2 NOTETR(',NTE2,')=',
ccc     %                 (NOTETR(M,NTE2),M=1,8),
ccc     %                 ' Q=',QUALTE,' V=',VOLUTE

C              LES TETRAEDRES OPPOSES A NTE
               NBT = 0
               DO 31 K=5,8
                  NTEOP=NOTETR(K,NTE)
C                 SON VOLUME ET SA QUALITE
                  CALL QUATETD( PTXYZD( 1, NOTETR(1,NTEOP) ),
     %                          PTXYZD( 1, NOTETR(2,NTEOP) ),
     %                          PTXYZD( 1, NOTETR(3,NTEOP) ),
     %                          PTXYZD( 1, NOTETR(4,NTEOP) ),
     %                          ARMIN, ARMAX, SURFTR, VOLUTE, QUALTE )
ccc                  PRINT *,'idstst5: Oppose',NTE,' NOTETR(',NTEOP,')=',
ccc     %                    (NOTETR(M,NTEOP),M=1,8),
ccc     %                    ' Q=',QUALTE,' V=',VOLUTE
                  IF( NTEOP .EQ. NTE .OR. NTEOP .EQ. NTE2 ) GOTO 31
C                 NTEOP EST IL DEJA DANS LA LISTE?
                  DO L=NBTENS+1,NBT
                     IF( NTEOP .EQ. NOTENS(L) ) GOTO 31
                  ENDDO
C                 NON: IL EST AJOUTE
                  NBT = NBT + 1
                  NOTENS(NBTENS+NBT) = NTEOP
 31            ENDDO

C              LES TETRAEDRES OPPOSES A NTE2
               DO 32 K=5,8
                  NTEOP=NOTETR(K,NTE2)
C                 SON VOLUME ET SA QUALITE
                  CALL QUATETD( PTXYZD( 1, NOTETR(1,NTEOP) ),
     %                          PTXYZD( 1, NOTETR(2,NTEOP) ),
     %                          PTXYZD( 1, NOTETR(3,NTEOP) ),
     %                          PTXYZD( 1, NOTETR(4,NTEOP) ),
     %                          ARMIN, ARMAX, SURFTR, VOLUTE, QUALTE )
ccc                  PRINT *,'idstst5: Oppose',NTE2,' NOTETR(',NTEOP,')=',
ccc     %                    (NOTETR(M,NTEOP),M=1,8),
ccc     %                    ' Q=',QUALTE,' V=',VOLUTE
                  IF( NTEOP .EQ. NTE .OR. NTEOP .EQ. NTE2 ) GOTO 32
C                 NTEOP EST IL DEJA DANS LA LISTE?
                  DO L=NBTENS+1,NBT
                     IF( NTEOP .EQ. NOTENS(L) ) GOTO 32
                  ENDDO
C                 NON: IL EST AJOUTE
                  NBT = NBT + 1
                  NOTENS(NBTENS+NBT) = NTEOP
 32            ENDDO

C              TRACE DES TETRAEDRES OPPOSES AUX 2 TETRAEDRES DOUBLES
               NBT = NBT + 1
               NOTENS(NBTENS+NBT) = NTE
               NBT = NBT + 1
               NOTENS(NBTENS+NBT) = NTE2

        KTITRE='idstst5:          TETRAEDRES de SOMMET            AVANT'
               WRITE(KTITRE(10:16),'(I7)') NBT
               WRITE(KTITRE(40:47),'(I8)') NSDLFR
               CALL TRACEETOILE( KTITRE, PTXYZD, NOTETR, NSDLFR,
     %                           NBT, NOTENS(NBTENS+1) )
ccc               CALL TRFETO13(KTITRE,PTXYZD,NBT,NOTENS(NBTENS+1),NOTETR)

               NBT = NBT - 2

C              SUPPRESSION DES 2 TETRAEDRES IDENTIQUES NTE et NTE2
               DO L=1,8
                  NOTETR(L,NTE) = 0
               ENDDO
               NOTETR( 5, NTE ) = N1TEVI
C              NTE DEVIENT LE PREMIER TETRAEDRE VIDE
               N1TEVI = NTE
               NOTENS( NT ) = -NTE
               IF( IVOLTE .NE. 0 ) NVOLTE( NTE ) = -1

               DO L=1,8
                  NOTETR(L,NTE2) = 0
               ENDDO
               NOTETR( 5, NTE2 ) = N1TEVI
C              NTE2 DEVIENT LE PREMIER TETRAEDRE VIDE
               N1TEVI = NTE2
               NOTENS( MT ) = -NTE2
               IF( IVOLTE .NE. 0 ) NVOLTE( NTE2 ) = -1

C              MISE A JOUR DES TETRAEDRES OPPOSES DES FACES
C              DE LEURS TETRAEDRES OPPOSES
               CALL MJOPTE( NBTENS+NBT, NOTENS, N1TETS, NOTETR, NUDTETR,
     %                      N1TEVI, PTXYZD, NBFANR )

      KTITRE='idstst6:         TETRAEDRES de SOMMET             APRES'
               WRITE(KTITRE(10:16),'(I7)') NBT
               WRITE(KTITRE(40:47),'(I8)') NSDLFR
               CALL TRACEETOILE( KTITRE, PTXYZD, NOTETR, NSDLFR,
     %                           NBTENS+NBT, NOTENS )
ccc               CALL TRFETO13( KTITRE,PTXYZD,NBTENS+NBT, NOTENS, NOTETR )

            ENDIF

 34      ENDDO

      ENDDO
ccc      PRINT*,'idstst: Nombre de TETRAEDRES DOUBLES=',NBTEDB


C     RECHERCHE DES TETRAEDRES OPPOSES AUX FACES DES TETRAEDRES
C     DE SOMMET ANCIENNEMENT NSNOFR
C     ---------------------------------------------------------
      L = NBTENS
      CALL COMPENTP( L, NOTENS,  NBTENS )

C     MISE A JOUR DES TETRAEDRES OPPOSES DANS UNE LISTE DE TETRAEDRES
C     ---------------------------------------------------------------
      CALL MJOPTE( NBTENS, NOTENS, N1TETS, NOTETR, NUDTETR,
     %             N1TEVI, PTXYZD, NBFANR )
C     NBFANR : NOMBRE DE FACES NON RETROUVEES DES TETRAEDRES
      IF( NBFANR .LE. 0 ) GOTO 70

C     LE NUMERO DE TETRAEDRE OPPOSE EST MIS A ZERO
      NBTEPE = 0
      DO NT=1,NBTENS

C        TETRAEDRE DE SOMMET NSNOFR=NSDLFR
         NTE = NOTENS( NT )

         DO 40 NF=1,4

C           LES 3 SOMMETS DE LA FACE NF DU TETRAEDRE NTE
            DO I=1,3
               NOSOTR(I) = NOTETR( NOSOFA(I,NF), NTE )
            ENDDO

C           LE TETRAEDRE OPPOSE A LA FACE NF EST IL A JOUR?
            NTEOP = NOTETR(4+NF,NTE)
            IF( NTEOP .GT. 0 .AND. NOTETR(1,NTEOP) .GT. 0 ) THEN
               CALL NUFATRTE( NOSOTR, NOTETR(1,NTEOP), NF1 )
               IF( NF1 .GT. 0 ) THEN
C                 LA FACE NF DE NTE EST LA FACE NF1 DE NTEOP
                  NOTETR(4+NF1,NTEOP) = NTE
                  NOTETR(4+NF ,NTE  ) = NTEOP
                  GOTO 40
               ENDIF
            ENDIF

C           RECHERCHE DE CETTE FACE PARMI LES TETRAEDRES DE SOMMET NSDLFR
            DO MT=1,NBTENS
               IF( MT .NE. NT ) THEN
                  NTE2 = NOTENS(MT)
                  IF( NTE2 .GT. 0 .AND. NOTETR(1,NTE2) .GT. 0 ) THEN
                     CALL NUFATRTE( NOSOTR, NOTETR(1,NTE2), NF1 )
                     IF( NF1 .GT. 0 ) THEN
C                       LA FACE NF1 DE NTE2 EST LA FACE NF DE NTE
                        NOTETR(4+NF ,NTE ) = NTE2
                        NOTETR(4+NF1,NTE2) = NTE
                        GOTO 40
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO

C           RECHERCHE DE CETTE FACE PARMI LES TETRAEDRES
C           OPPOSES AUX FACES DE L'ETOILE DES TETRAEDRES DE SOMMET NSDLFR
            NF1 = N1FEOC
 38         IF( NF1 .GT. 0 ) THEN
C              LE TETRAEDRE OPPOSE A LA FACE NF1
               NTEOP = NFETOI(1,NF1)
               IF( NTEOP .GT. 0 .AND. NOTETR(1,NTEOP) .GT. 0 ) THEN
                  CALL NUFATRTE( NOSOTR, NOTETR(1,NTEOP), NF2 )
                  IF( NF2 .GT. 0 ) THEN
C                    LA FACE NF1 DE NTE2 EST LA FACE NF DE NTE
                     NOTETR(4+NF ,NTE  ) = NTEOP
                     NOTETR(4+NF2,NTEOP) = NTE
                     GOTO 40
                  ENDIF
               ENDIF
C              FACE SUIVANTE
               NF1 = NFETOI(5,NF1)
               GOTO 38
            ENDIF

            NBTEPE = NBTEPE + 1
            PRINT*
            PRINT*,'idstst: Pb TETRA',NTE,' FACE',NF,' NON RETROUVEE'
            PRINT*,'idstst: NOTETR(',NTE,')=',(NOTETR(I,NTE),I=1,8)

 40      ENDDO

C        TETRAEDRE FINAL
C        SON VOLUME ET SA QUALITE
         CALL QUATETD( PTXYZD( 1, NOTETR(1,NTE) ),
     %                 PTXYZD( 1, NOTETR(2,NTE) ),
     %                 PTXYZD( 1, NOTETR(3,NTE) ),
     %                 PTXYZD( 1, NOTETR(4,NTE) ),
     %                 ARMIN, ARMAX, SURFTR, VOLUTE, QUALTE )
ccc         PRINT *,'idstst: ',NT,' final  NOTETR(',NTE,')=',
ccc     %           (NOTETR(M,NTE),M=1,8),' Q=',QUALTE,' V=',VOLUTE

      ENDDO

      IF( NBTEPE .GT. 0 ) THEN

C        DES TETRAEDRES OPPOSES N'ONT PAS ETE RETROUVES
C        AJOUT DES TETRAEDRES OPPOSES AUX FACES SIMPLES DE L'ETOILE
         NF1 = N1FEOC
 55      IF( NF1 .GT. 0 ) THEN

C           LE TETRAEDRE OPPOSE A LA FACE NF1
            NTEOP = NFETOI(1,NF1)

C           NTEOP EST IL DANS NOTENS?
            DO NT=1,NBTENS
               IF( NTEOP .EQ. NOTENS(NT) ) GOTO 60
            ENDDO
C           NTEOP EST AJOUTE A NOTENS
            NBTENS = NBTENS + 1
            NOTENS( NBTENS ) = NTEOP

C           FACE SUIVANTE
 60         NF1 = NFETOI(5,NF1)
            GOTO 55

         ENDIF

C        MISE A JOUR DES TETRAEDRES OPPOSES DANS UNE LISTE DE TETRAEDRES
C        ---------------------------------------------------------------
         CALL MJOPTE( NBTENS, NOTENS, N1TETS, NOTETR, NUDTETR,
     %                N1TEVI, PTXYZD, NBFANR )
C        NBFANR : NOMBRE DE FACES NON RETROUVEES DES TETRAEDRES
C                 LE NUMERO DE TETRAEDRE OPPOSE EST MIS A ZERO

      ENDIF

cccC     VERIFICATION DES TETRAEDRES OPPOSES PAR LES FACES
cccC     -------------------------------------------------
ccc      CALL VEOPTE( NBTENS, NOTENS, NOTETR, PTXYZD, NBFANR )

C     VERIFICATION DES TETRAEDRES MODIFIES et OPPOSES DE SOMMET NSDLFR
C     ----------------------------------------------------------------
C     RECUPERATION DES NBTENS TETRAEDRES DE SOMMET NSDLFR
 70   CALL TETR1S( NSDLFR, N1TETS, NOTETR,
     %             NBTENS, MXTENS, NOTENS, IERR )
      IF( IERR .NE. 0 ) THEN
         PRINT*,'idstst: PB2 sur les TETRAEDRES de SOMMET',NSDLFR
         IERR = 2
         GOTO 9999
      ENDIF

C     VERIFICATION DES TETRAEDRES OPPOSES PAR LES FACES 
      CALL VEOPTE( NBTENS, NOTENS, NOTETR, PTXYZD, NBFANR )

C     MISE A JOUR QUALITE ET VOLUME
ccc      PRINT *
      VOLUE1 = 0D0
      QUAMOY = 0.0
      QUAMIN = 2.0
      DO NT=1,NBTENS

C        TETRAEDRE DE SOMMET NSDLFR
         NTE = NOTENS( NT )

         DO I=1,4
C           MISE A JOUR DU NO DE TETRAEDRE DE CHAQUE SOMMET
            N1TETS( NOTETR(I,NTE) ) = NTE
         ENDDO

C        SON VOLUME ET SA QUALITE
         CALL QUATETD( PTXYZD( 1, NOTETR(1,NTE) ),
     %                 PTXYZD( 1, NOTETR(2,NTE) ),
     %                 PTXYZD( 1, NOTETR(3,NTE) ),
     %                 PTXYZD( 1, NOTETR(4,NTE) ),
     %                 ARMIN, ARMAX, SURFTR, VOLUTE, QUALTE )
         VOLUE1 = VOLUE1 + VOLUTE

         IF( QUALTE .LE. 0. ) THEN
ccc            TRACTE = .TRUE.
            PRINT *,'idstst: Probleme TETRAEDRE(',NTE,')=',
     %              (NOTETR(M,NTE),M=1,8),' Q=',QUALTE,' V=',VOLUTE
         ENDIF
         QUAMOY = QUAMOY + QUALTE

         IF( QUALTE .LT. QUAMIN ) THEN
            QUAMIN = QUALTE
            NTEMIN = NTE
         ENDIF

      ENDDO
      QUAMOY = QUAMOY / NBTENS

ccc      PRINT *,'idstst : St',NSNOFR,'->',NSDLFR,' avec',NBTENS,
ccc     %        ' TETRAEDRES de Volume=',VOLUE1,
ccc     %        ' Q MIN=',QUAMIN,' Q MOY=',QUAMOY,
ccc     %        ' NBFANR=',NBFANR,' IERR=',IERR

cccC     TETRAEDRE DE QUALITE MINIMALE
ccc         CALL QUATETD( PTXYZD( 1, NOTETR(1,NTEMIN) ),
ccc     %                 PTXYZD( 1, NOTETR(2,NTEMIN) ),
ccc     %                 PTXYZD( 1, NOTETR(3,NTEMIN) ),
ccc     %                 PTXYZD( 1, NOTETR(4,NTEMIN) ),
ccc     %                 ARMIN, ARMAX, SURFTR, VOLUTE, QUAMIN )
ccc      PRINT *,'idstst : St',NSNOFR,'->',NSDLFR,' Q MIN=',QUAMIN,
ccc     %' V MIN=',VOLUTE,' NOTETR(',NTEMIN,')=',(NOTETR(K,NTEMIN),K=1,8)

cccC     TRACE DES TETRAEDRES FINAUX APRES IDENTIFICATION
ccc      TRACTE  = .TRUE.
ccc      KTITRE='idstst9:           TETRAEDRES de SOMMET               FIN'
ccc      WRITE(KTITRE(10:16),'(I7)') NBTENS
ccc      WRITE(KTITRE(40:47),'(I8)') NSDLFR
ccc      CALL TRACEETOILE( KTITRE, PTXYZD, NOTETR, NSDLFR,
ccc     %                  NBTENS, NOTENS )
cccccc      CALL TRFETO13( KTITRE, PTXYZD, NBTENS, NOTENS, NOTETR )
ccc      TRACTE = TRACTE0

 9999 RETURN
      END
