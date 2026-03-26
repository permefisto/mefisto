      SUBROUTINE QUA1ST( QUTEMN, ANTE2P, GRAND,  NS,
     %                   MXSOMM, XYZSOM, NPSOFR, NBDM,   NUDMEF,
     %                   NBSOTE, MXTETR, N1TEVI, NSTETR, MXTETA,
     %                   VOLUMT, QUALIT,
     %                   NO1TSO, MXTESO, N1TESO, NOTESO,
     %                   MXFAET, N1FEOC, N1FEVI, NFETOI, NOSOET,
     %                   IERR  )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   AMELIORER LA QUALITE DE LA TETRAEDRISATION AUTOUR DU SOMMET NS
C -----
C
C ENTREES:
C --------
C QUTEMN : QUALITE MINIMALE AU DESSOUS DE LAQUELLE UN TETRAEDRE DOIT ETRE
C          AMELIORE
C ANTE2P : ANGLE AU DESSOUS DUQUEL 2 FACES ADJACENTES SONT CONSIDEREES
C          COPLANAIRES
C GRAND  : PLUS GRAND REEL STOCKABLE
C NS     : NUMERO DU SOMMET DE QUALITE A AMELIORER
C MXSOMM : NOMBRE MAXIMAL DE SOMMETS DECLARABLES DANS XYZSOM ET NPSOFR
C XYZSOM : COORDONNEES X Y Z DES NBSOMM SOMMETS DES TETRAEDRES
C NBDM   : 0 SI 1 MATERIAU=VOLUME, SINON NOMBRE DE MATERIAUX DU VOLUME
C NBSOTE : VALEUR MAXIMALE DE DECLARATION DU PREMIER INDICE DE NSTETR(>3)
C MXTETR : NOMBRE MAXIMAL DE TETRAEDRES DECLARABLES DANS NSTETR
C MXTESO : NOMBRE MAXIMAL DE NUMERO DE TETRAEDRES DES SOMMETS
C MXFAET : NOMBRE MAXIMAL DE FACES DECLARABLES DANS NFETOI
C
C TABLEAUX AUXILIAIRES :
C ----------------------
C NFETOI(5,MXFAET)  DES ENTIERS
C NOSOET(MXTETR)    DES ENTIERS
C
C MODIFIES :
C ----------
C NBSOMM : NOMBRE ACTUEL DE SOMMETS DE LA TETRAEDRISATION
C NPSOFR : NUMERO 0 SI SOMMET INTERNE
C                 1 SI SOMMET SUR LA FRONTIERE
C                 2 SI SOMMET SUR L'INTERFACE ENTRE 2 MATERIAUX
C                -1 SI SOMMET SUPPRIME LORS DE L'AMELIORATION
C NUDMEF : NUMERO DE MATERIAU DE CHAQUE TETRAEDRE DU MAILLAGE
C          ATTENTION: CE TABLEAU EXISTE SEULEMENT SI NBDM>0
C N1TEVI : NUMERO DANS NSTETR DE LA PREMIERE PLACE VIDE
C          CHAINAGE SUIVANT DANS NSTETR(2,*) ET DERNIER A ZERO
C N1TESO : NUMERO DANS NOTESO DE LA PREMIERE PLACE VIDE
C          CHAINAGE SUIVANT DANS NOTESO(2,*) ET DERNIER A ZERO
C NSTETR : NUMERO DES 4 SOMMETS DE CHAQUE TETRAEDRE
C MXTETA : NUMERO DU PLUS GRAND TETRAEDRE UTILISE
C N1FEOC : NUMERO NFETOI DE LA PREMIERE FACE OCCUPEE
C N1FEVI : NUMERO NFETOI DE LA PREMIERE FACE VIDE
C
C SORTIES:
C --------
C NO1TSO : NUMERO DU 1-ER TETRAEDRE DANS NOTESO DE CHAQUE SOMMET
C NOTESO : NOTESO(1,*) NUMERO DANS NSTETR DU TETRAEDRE
C          NOTESO(2,*) NUMERO DANS NOTESO DU TETRAEDRE SUIVANT
C                      0 SI C'EST LE DERNIER
C QUALIT : QUALITE DES TETRAEDRES DE LA TETRAEDRISATION
C VOLUMT : VOLUME  DES TETRAEDRES DE LA TETRAEDRISATION
C IERR   : =0 SI PAS D'ERREUR
C          >0 EN CAS DE SATURATION D'UN TABLEAU
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC   DECEMBRE 1991
C MODIFS : ALAIN PERRONNET LJLL UPMC et St Pierre du Perray    Juin 2008
C2345X7..............................................................012
      include"./incl/langue.inc"
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      INTEGER           NPSOFR(MXSOMM),
     %                  NUDMEF(MXTETR),
     %                  NSTETR(NBSOTE,MXTETR),
     %                  NO1TSO(MXSOMM),
     %                  NOTESO(2,MXTESO),
     %                  NFETOI(5,MXFAET),
     %                  NOSOET(MXTETR)
      REAL              XYZSOM(3,MXSOMM),VOLUMT(MXTETR),
     %                  QUALIT(MXTETR)
      REAL              BB(3)
      REAL              ARMIN,ARMAX,SURFTR(4)
C
      IF( NPSOFR(NS) .EQ. 2 ) RETURN
C     NS EST SUR UN INTERFACE ENTRE 2 MATERIAUX => IL N'EST PAS MODIFIE
C     =================================================================
C
C     COSINUS DE L'ANGLE DE COPLANEARITE CONVERTI EN RADIANS
      COSE2P = COS( ANTE2P * ATAN(1.0) / 45.0 )
C
C     CALCUL DU VOLUME ET DE LA QUALITE DES TETRAEDRES DE SOMMET NS
      IF( NO1TSO( NS ) .LE. 0 ) GOTO 9000
      CALL QUALST( NS, XYZSOM, NBSOTE, NSTETR, NO1TSO, NOTESO,
     %             VOLUMT, QUALIT,
     %             VOLUM0, QUALI0, NTQMIN, NBTENS )
      IF( QUALI0 .GT. QUTEMN .OR. NBTENS .LE. 0 ) GOTO 9000
      QUAL00 = QUALI0
      QUALI1 = 0
      NBS    = 0
C
      IF( NPSOFR(NS) .EQ. 0 ) THEN
C
C        ===============================================
C        LE POINT NS N'EST PAS SUR LA SURFACE DE L'OBJET
C        ===============================================
C
C        CALCUL DU POINT DE QUALITE MINIMALE MAXIMISEE
C        =============================================
C        RECHERCHE DU TRIANGLE OPPOSE AU SOMMET NS DANS LE TETRAEDRE
C        DE QUALITE MINIMALE
 5       CALL QUATET( XYZSOM(1,NSTETR(1,NTQMIN)),
     %                XYZSOM(1,NSTETR(2,NTQMIN)),
     %                XYZSOM(1,NSTETR(3,NTQMIN)),
     %                XYZSOM(1,NSTETR(4,NTQMIN)),
     %   ARMIN, ARMAX, SURFTR, VOLUMT(NTQMIN), QUALIT(NTQMIN) )
         DO 10 J=1,4
            IF( NSTETR(J,NTQMIN) .EQ. NS ) GOTO 15
 10      CONTINUE
C
         WRITE(IMPRIM,*)
         WRITE(IMPRIM,*) 'QUA1ST: ANOMALIE A TRAITER NS=',NS,
     %   ' NON SOMMET de NSTETR(',NTQMIN,')=',(NSTETR(I,NTQMIN),I=1,4)
         CALL XVPAUSE
         RETURN
C
 15      IF( J .LT. 4 ) THEN
            I = J + 1
         ELSE
            I = 1
         ENDIF
C        LA FACE I DE NTQMIN EST OPPOSEE AU SOMMET J=NS
C        LA HAUTEUR ISSUE DE NS = 6V / 2S = HAUTEUR REELLE/3
CCC         DMIN = VOLUMT( NTQMIN ) / SURFTR(I) * 3.0
C
C        I I1 I2 LES 3 SOMMETS DE LA FACE OPPOSEE AU SOMMET J=NS
         I1 = MOD( I,4) + 1
         I2 = MOD(I1,4) + 1
C
C        LE BARYCENTRE DE LA FACE OPPOSEE ET LA DIRECTION NS-BARYCENTRE
 18      DO 20 K=1,3
            BB(K) =  XYZSOM(K,NS) -
     %             ( XYZSOM(K,NSTETR(I ,NTQMIN))
     %              +XYZSOM(K,NSTETR(I1,NTQMIN))
     %              +XYZSOM(K,NSTETR(I2,NTQMIN)) ) / 3.0
 20      CONTINUE
         XX   = BB(1)
         YY   = BB(2)
         ZZ   = BB(3)
         SS   = SQRT( BB(1)**2 + BB(2)**2 + BB(3)**2 )
         IPAS = 1
C
C        CALCUL DU POINT DE QUALITE MAXIMISEE SUR LA DROITE
C        PASSANT PAR NS ET DE DIRECTION NS-BARYCENTRE DE LA FACE OPPOSEE
 30      CALL QTEMXD( NS,     BB,
     %                XYZSOM, NBSOTE, NSTETR, NO1TSO, NOTESO,
     %                VOLUMT, QUALIT,
     %                VOLUM1, QUALI1, NTQMI1, NBTENS )
         IF( NBTENS .LE. 0 ) GOTO 9000
         IF( NTQMI1 .NE. NTQMIN ) THEN
C           ESSAI D'UNE NOUVELLE DIRECTION
            NTQMIN = NTQMI1
            GOTO 5
         ENDIF
C
         IF( QUALI1-QUALI0  .GT. 0.001  .AND.
     %      (QUALI1-QUALI0) .GT. 0.01*QUALI0 ) THEN
C
C           LA QUALITE A AUGMENTE DE 1 %
            QUALI0 = QUALI1
C           ESSAI DE 2 DIRECTIONS ORTHOGONALES
C           (X Y Z) => (0 Z -Y) => (Y*Y+Z*Z -XY XZ)
            IF( IPAS .EQ. 1 ) THEN
               BB(1) = 0
               BB(2) = ZZ
               BB(3) =-YY
               S = SQRT( BB(2)**2 + BB(3)**2 )
               IF( S .LE. 0 ) GOTO 90
               BB(2) = BB(2) / S * SS
               BB(3) = BB(3) / S * SS
               IPAS  = 2
               GOTO 30
            ELSE IF( IPAS .EQ. 2 ) THEN
               BB(1) = YY*YY+ZZ*ZZ
               BB(2) =-XX*YY
               BB(3) = XX*ZZ
               S = SQRT( BB(1)**2 + BB(2)**2 + BB(3)**2 )
               IF( S .LE. 0 ) GOTO 90
               BB(1) = BB(1) / S * SS
               BB(2) = BB(2) / S * SS
               BB(3) = BB(3) / S * SS
               IPAS  = 3
               GOTO 30
            ELSE
CCC               IPAS = 4
               GOTO 18
            ENDIF
         ENDIF
C
CCC 90      IF( QUAL00 .LT. QUALI1 ) THEN
CCCC           XYZSOM(1:3,NS) EST AU MAXIMUM DE LA QUALITE
CCC            WRITE(IMPRIM,*) 'OUT QTEDMX: NS=',NS,
CCC     %      ' Q0=',QUAL00,' Q1=',QUALI1
CCC         ENDIF
C
 90      QUALI0 = QUALI1
C
C        ESSAI DE TETRAEDRISER L'ETOILE EN SUPPRIMANT LE POINT NS
C        ========================================================
C        LES FACES CHAINEES UNIQUES DE L'ETOILE DU POINT NS
         CALL FACETO( NS,     NBSOTE, NSTETR, NO1TSO, NOTESO,
     %                N1FEVI, N1FEOC, NFETOI, IERR )
         IF( IERR .NE. 0 ) THEN
            IERR = 0
            GOTO 9000
         ENDIF
C
C        LA LISTE DES SOMMETS DES FACES UNIQUES DE L'ETOILE
         NBS = 0
         NF2 = N1FEOC
C        BOUCLE SUR LES FACES DE L'ETOILE
 150     IF( NF2 .GT. 0 ) THEN
            DO 170 I=1,3
C              LE NUMERO K DU SOMMET I DE LA FACE NF2
               K = NFETOI(I,NF2)
               IF( K .EQ. NS ) GOTO 170
C              SOMMET K DIFFERENT DE NS
               DO 160 J=1,NBS
                  IF( K .EQ. NOSOET(J) ) GOTO 170
 160           CONTINUE
C              SOMMET K NON RETROUVE  DONC A AJOUTER
               NBS = NBS + 1
               NOSOET( NBS ) = K
 170        CONTINUE
C           PASSAGE A LA FACE SUIVANTE
            NF2 = NFETOI(5,NF2)
            GOTO 150
         ENDIF
C
      ELSE IF( NPSOFR(NS) .EQ. 1 ) THEN
C
C        ==========================
C        CAS D'UN SOMMET FRONTALIER
C        ==========================
C        LES FACES CHAINEES UNIQUES DE L'ETOILE DU POINT NS
         CALL FACETO( NS,     NBSOTE, NSTETR, NO1TSO, NOTESO,
     %                N1FEVI, N1FEOC, NFETOI, IERR )
         IF( IERR .NE. 0 ) THEN
            IERR = 0
            GOTO 9000
         ENDIF
C
C        CALCUL DU NOMBRE DE NORMALES AUX FACES SIMPLES DE SOMMET NS
         CALL ANSETO( NS,    N1FEOC, NFETOI, XYZSOM, COSE2P,
     %                NBNORM )
         IF( NBNORM .GE. 3 ) THEN
C
C           AU MOINS UN TRIEDRE DE SOMMET NS . PAS DE TRAITEMENT
C           --------------------------------
            GOTO 9000
C
         ELSE IF( NBNORM .EQ. 2 ) THEN
C
C           NS EST SUR LA DROITE D'INTERSECTION DE 2 PLANS
C           ----------------------------------------------
C           IL FAUT RESPECTER CETTE DROITE AFIN DE CONSERVER LA GEOMETRIE
C           RECHERCHE DES 2 SOMMETS EXTREMITES DE CETTE ARETE
C           LA LISTE DES SOMMETS DES FACES DE SOMMET NS
CCC
CCC         LE PROBLEME: 1 FACE DANS 3 TETRAEDRES EST DU a CETTE SEQUENCE
CCC         DONC ELLE EST SUPPRIMEE PAR LE if suivant qui est toujours VRAI
            if( ns .gt. 0 ) return
C
            NBS = 0
            NF2 = N1FEOC
C           BOUCLE SUR LES FACES DE L'ETOILE
 200        IF( NF2 .GT. 0 ) THEN
               DO 205 I=1,3
                  IF( NFETOI(I,NF2) .EQ. NS ) GOTO 208
 205           CONTINUE
               GOTO 225
C              FACE DE SOMMET NS
 208           DO 220 I=1,3
C                 LE NUMERO DU SOMMET I DE LA FACE NF2
                  K = NFETOI(I,NF2)
                  IF( K .EQ. NS ) GOTO 220
                  DO 210 J=1,NBS
                     IF( K .EQ. NOSOET(J) ) GOTO 220
 210              CONTINUE
C                 SOMMET NON RETROUVE  DONC A AJOUTER
                  NBS = NBS + 1
                  NOSOET( NBS ) = K
 220           CONTINUE
C              PASSAGE A LA FACE SUIVANTE
 225           NF2 = NFETOI(5,NF2)
               GOTO 200
            ENDIF
C
C           RECHERCHE DE L'ARETE DES 2 PLANS
            DMIN = GRAND
            I1   = 0
            I2   = 0
            DO 240 I=1,NBS-1
               DO 230 J=I+1,NBS
C                 NS EST IL SUR LA DROITE I-J ?
                  D = DISTPD( XYZSOM(1,NS),
     %                        XYZSOM(1,NOSOET(I)),
     %                        XYZSOM(1,NOSOET(J)) )
                  IF( D .LT. DMIN ) THEN
C                    NS EST LE PLUS PROCHE DE NOSOET(I1) NOSOET(I2)
                     DMIN = D
                     I1   = I
                     I2   = J
                  ENDIF
 230           CONTINUE
 240        CONTINUE
C
C           DMIN EST LA DISTANCE LA PLUS COURTE AUX ARETES POSSIBLES
            D = DIST2P( XYZSOM(1,NOSOET(I1)),
     %                  XYZSOM(1,NOSOET(I2)) )
            IF( DMIN .GE. D * 0.005 ) THEN
               GOTO 9000
            ELSE
C
C              ESSAI DE DEPLACER LE SOMMET NS DANS LE SEGMENT I1-I2
C              LA DIRECTION
               BB(1) = (XYZSOM(1,NOSOET(I1))-XYZSOM(1,NOSOET(I2)))*0.5
               BB(2) = (XYZSOM(2,NOSOET(I1))-XYZSOM(2,NOSOET(I2)))*0.5
               BB(3) = (XYZSOM(3,NOSOET(I1))-XYZSOM(3,NOSOET(I2)))*0.5
               CALL QTEMXD( NS,     BB,
     %                      XYZSOM, NBSOTE, NSTETR, NO1TSO, NOTESO,
     %                      VOLUMT, QUALIT,
     %                      VOLUM1, QUALI1, NTQMI1, NBTENS )
               IF( NBTENS .LE. 0 ) GOTO 9000
               IF( QUALI1 .GT. QUALI0 ) THEN
C                 XYZSOM(1:3,NS) EST AU MAXIMUM DE LA QUALITE
CCC                  WRITE(IMPRIM,*) 'DEPLACEMENT DE NS=',NS,
CCC     %           ' SUR LA DROITE AVEC 2 NORMALES Q0=',QUALI0,
CCC     %           ' Q1=',QUALI1
                  QUALI0 = QUALI1
               ENDIF
            ENDIF
C
C           LA DISTANCE DE NS A LA DROITE I1-I2 EST TRES FAIBLE
C           LES 2 SOMMETS EXTREMITES DE L'ARETE A ESSAYER
            NBS = 2
            NOSOET(1) = NOSOET(I1)
            NOSOET(2) = NOSOET(I2)
C
         ELSE IF( NBNORM .EQ. 1 ) THEN
C
C           TOUTES LES FACES DE SOMMET NS SONT COPLANAIRES
C           ----------------------------------------------
CCCC           ESSAI DE METTRE NS AU BARYCENTRE DES FACES COPLANAIRES
CCC            P(1) = 0.
CCC            P(2) = 0.
CCC            P(3) = 0.
C           LES SOMMETS DES FACES COPLANAIRES SONT ESSAYES
C           EN SUPPRIMANT LE SOMMET NS DE L'ETOILE
            NBS = 0
            NF2 = N1FEOC
C
C           BOUCLE SUR LES FACES DE L'ETOILE
 250        IF( NF2 .GT. 0 ) THEN
C              CETTE FACE CONTIENT ELLE LE SOMMET NS ?
               DO 255 I=1,3
                  IF( NFETOI(I,NF2) .EQ. NS ) GOTO 258
 255           CONTINUE
               GOTO 280
C              OUI
 258           DO 270 I=1,3
C                 LE NUMERO DU SOMMET I DE LA FACE NF2
                  K = NFETOI(I,NF2)
                  IF( K .EQ. NS ) GOTO 270
                  DO 260 J=1,NBS
                     IF( K .EQ. NOSOET(J) ) GOTO 270
 260              CONTINUE
C
C                 SOMMET NON RETROUVE  DONC A AJOUTER
                  NBS = NBS + 1
                  NOSOET( NBS ) = K
CCC                  DO 265 J=1,3
CCC                     P(J) = P(J) + XYZSOM(J,K)
CCC 265              CONTINUE
C
C                 ESSAI D'ECHANGER LES DIAGONALES EN SUPPRIMANT NS-K
C                 --------------------------------------------------
                  CALL QU2T2T( NS,     K,
     %                         XYZSOM, NBSOTE, N1TEVI, NSTETR,
     %                         NBDM,   NUDMEF,
     %                         NO1TSO, N1TESO, NOTESO,
     %                         NBTEDS, NOSOET(NBS+1),
     %                         MXTETA, VOLUMT, QUALIT, IERR )
                  IF( IERR .EQ. -1 ) THEN
C                    REUSSITE
CCC                     WRITE(IMPRIM,*) ' 2T=>2T ',NS,' ',K
                     IERR = 0
                     GOTO 9000
                  ENDIF
C
C                 ESSAI DE DEPLACER LE SOMMET NS SUR LA DROITE
C                 XYZSOM(.,NS)-XYZSOM(.,K)
C                 --------------------------------------------
C                 LA DIRECTION
                  BB(1)= XYZSOM(1,NS) - XYZSOM(1,K)
                  BB(2)= XYZSOM(2,NS) - XYZSOM(2,K)
                  BB(3)= XYZSOM(3,NS) - XYZSOM(3,K)
                  CALL QTEMXD( NS,     BB,
     %                         XYZSOM, NBSOTE, NSTETR, NO1TSO, NOTESO,
     %                         VOLUMT, QUALIT,
     %                         VOLUM1, QUALI1, NTQMI1, NBTENS )
                  IF( NBTENS .LE. 0 ) GOTO 9000
 270           CONTINUE
C
C              PASSAGE A LA FACE SUIVANTE
 280           NF2 = NFETOI(5,NF2)
               GOTO 250
            ENDIF
            IF( QUALI1 .GT. QUAL00 ) THEN
C              XYZSOM(1:3,NS) EST AU MAXIMUM DE LA QUALITE
CCC               WRITE(IMPRIM,*) 'DEPLACEMENT DE NS=',NS,
CCC     %       ' DANS LE PLAN 1 NORMALE Q0=',QUALI0,' Q1=',QUALI1
               QUALI0 = QUALI1
            ENDIF
CCCC
CCCC           LE BARYCENTRE DES FACES COPLANAIRES DE SOMMET NS
CCC            DO 285 J=1,3
CCC               P(J) = P(J) / NBS
CCC 285        CONTINUE
CCCC
CCCC           QUALITE DE L'ETOILE DU SOMMET NS DEPLACE AU POINT P
CCC            CALL VOQUET( P, NS, XYZSOM, NO1TSO, NOTESO,
CCC     %                   NBSOTE, NSTETR,
CCC     %                   VOLUM0, QUALI1 )
CCCC
CCCC           BILAN NS EST LE POINT DE MEILLEURE QUALITE
CCC           IF( QUALI1 .GT. QUALI0 ) THEN
CCCC              LE POINT P EST MEILLEUR QUE NS
CCC               XYZSOM(1,NS) = P(1)
CCC               XYZSOM(2,NS) = P(2)
CCC               XYZSOM(3,NS) = P(3)
CCC               QUALI0 = QUALI1
CCC               WRITE(IMPRIM,*)'QUA1ST 285:NS=',NS,
CCC     %       ' QUALITE AMELIOREE ',NBS,' SOMMETS ESSAYESS Q1=',QUALI1
CCCC              MISE A JOUR DU VOLUME ET QUALITE DES TETRAEDRES
CCC               CALL QUALST( NS,XYZSOM,NBSOTE,NSTETR,NO1TSO,NOTESO,
CCC     %                      VOLUMT, QUALIT,
CCC     %                      VOLUM0, QUALI0, NTQMIN, NBTENS )
CCC               IF( NBTENS .LE. 0 ) GOTO 9000
CCC            ENDIF
         ENDIF
      ENDIF

C     IL Y A NBS SOMMETS A ESSAYER AFIN DE SUPPRIMER LE SOMMET NS
C     ===========================================================
      IMIEUX = 0
      DO I=1,NBS
C        ESSAI DU SOMMET NOSOET(I)
ccc         IF( NOSOET(I) .EQ. NS ) GOTO 320
         CALL VLQLET( NOSOET(I), XYZSOM, N1FEOC, NFETOI,
     %                VOLUM0, QUALI1 )
         IF( QUALI1 .GT. QUALI0 ) THEN
C           QUALITE MEILLEURE AVEC CE SOMMET
            IMIEUX = I
C           LE MAXIMUM DE LA QUALITE POSSIBLE
            QUALI0 = QUALI1
         ENDIF
      ENDDO

C     SI UN POINT OFFRE UNE MEILLEURE QUALITE
C     ALORS SUPPRESSION DU SOMMET NS
C           GENERATION DE LA TETRAEDRISATION DE L'ETOILE
C           A PARTIR DE CE POINT
      IF( IMIEUX .GT. 0 ) THEN
C
C        NO DU SOMMET DE MEILLEUR QUALITE
         NSMIEU = NOSOET(IMIEUX)
         IF( NSMIEU .EQ. NS ) GOTO 9000
C
C        DESTRUCTION DES TETRAEDRES DE L'ETOILE DU SOMMET NS
         CALL DSTETO( NS,     NBSOTE, N1TEVI, NSTETR, NBDM, NUDMEF,
     %                NO1TSO, N1TESO, NOTESO, QUALIT, NUMATE )
         IF( NUMATE .EQ. -2 ) THEN
C           TETRAEDRES DE SOMMET NS DANS DIFFERENTS MATERIAUX
            IERR = 0
            GOTO 9000
         ENDIF
C
C        GENERATION DES NOUVEAUX TETRAEDRES REMPLISSANT L'ETOILE
         CALL GETETO( NSMIEU, XYZSOM, N1FEOC, N1FEVI, NFETOI,
     %                NUMATE, NUDMEF,
     %                NBSOTE, N1TEVI, NSTETR,
     %                NO1TSO, N1TESO, NOTESO,
     %                MXTETA, VOLUMT, QUALIT, QUALI0, IERR )
         IF( IERR .GT. 0 ) GOTO 9000
C
ccc         IF( LANGAG .EQ. 0 ) THEN
ccc            WRITE(IMPRIM,*)'QUA1ST: ST',NSMIEU,
ccc     %    ' de MEILLEURE QUALITE',QUALI0,' que',NS,' de QUALITE',QUAL00
ccc            WRITE(IMPRIM,*)'QUA1ST: Le SOMMET',NS,' est SUPPRIME',
ccc     %    ' POUR',NBS,' SOMMETS ESSAYES'
ccc         ELSE
ccc            WRITE(IMPRIM,*)'QUA1ST: Vertex',NSMIEU,
ccc     %  ' with a better quality',QUALI0,' than ',NS,' of QUALITY',QUAL00
ccc            WRITE(IMPRIM,*)'QUA1ST: The VERTEX',NS,' is DELETED'
ccc         ENDIF
C
C        TEMOIN DE SUPPRESSION DU SOMMET NS
         NPSOFR( NS ) = -1
C
ccc         print *,'qua1st: suppression sommet',NS,
ccc     %           ' pour le sommet', nsmieu
C
      ENDIF
C
C     LES FACES OCCUPEES DE NFETOI DEVIENNENT VIDES
 9000 CALL ETOIOCVI( N1FEOC, N1FEVI, NFETOI )
C
      RETURN
      END
