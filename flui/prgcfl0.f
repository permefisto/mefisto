      SUBROUTINE PRGCFL0( MNTOPO,  MNNPEF, MNXYZN,
     %                    NBNOEUD, NDDLNO, NCODSA,
     %                    MNLPLI,  MNLPCO, IERR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    GENERER LE STOCKAGE MORSE DE LA MATRICE SYMETRIQUE DU MAILLAGE
C -----    POUR NDDLNO DEGRES DE LIBERTE ( NON CONSTANT ) PAR NOEUD
C
C ENTREES:
C --------
C MNTOPO : ADRESSE  MCN DU TABLEAU   'TOPOLOGIE' DE L'OBJET
C MNNPEF : ADRESSES MCN DES TABLEAUX 'NPEF"TYPE_EF' DE L'OBJET
C MNXYZN : ADRESSE  MCN DU TABLEAU   'XYZNOEUD'  DE L'OBJET
C NBNOEUD: NOMBRE DE NOEUDS DU MAILLAGE
C      (Sommets pour BREZZI-FORTIN, Sommets+Milieu aretes pour TAYLOR-HOOD)
C NDDLNO : POINTEUR SUR LE DERNIER DL DE CHAQUE NOEUD (0:NBNOEUD)
C NTDL   : NOMBRE TOTAL DE DEGRES DE LIBERTE (Pour BF SANS CEUX DES BARYCENTRES)
C          REMARQUE NTDL=NDDLNO(NBNOEUD)
C NCODSA : CODE DE STOCKAGE DE LA MATRICE MORSE
C          -1 : MATRICE NON SYMETRIQUE (A PROGRAMMER)
C           1 : MATRICE SYMETRIQUE
C SORTIES:
C --------
C MNLPLI : POINTEUR SUR LES COEFFICIENTS DIAGONAUX DE LA MATRICE MORSE
C          MCN(MNLPLI+NODL)=ADRESSE DU COEF DIAGONAL NODL DANS LA MATRICE
C MNLPCO : MCN(MNLPCO-1+N) =NO DE LA COLONNE DU COEFFICIENT N
C IERR   : 0 SI PAS D'ERREUR
C          1 SI UN TABLEAU 'NPEF"TYPE_EF' EST INCONNU
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Laboratoire J-L.LIONS UPMC Paris Janvier 2009
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/donela.inc"
      include"./incl/donele.inc"
      include"./incl/a___xyznoeud.inc"
      include"./incl/a___npef.inc"
      include"./incl/a_objet__topologie.inc"
      include"./incl/gsmenu.inc"
C
      include"./incl/pp.inc"
      COMMON          MCN(MOTMCN)
      COMMON /UNITES/ LECTEU, IMPRIM, NUNITE(30)
      INTEGER         MNNPEF(1:*), NDDLNO(0:NBNOEUD)
C
      IERR   = 0
      MNLPCO = 0
      IF( NCODSA .EQ. 0 ) RETURN
C
C     NTDL NOMBRE TOTAL DE DEGRES DE LIBERTE DES NOEUDS DU MAILLAGE
C     (Pour BREZZI-FORTIN SANS LES DL DES BARYCENTRES)
      NTDL = NDDLNO( NBNOEUD )
C
C     ADRESSE MCN DES 2 TABLEAUX LIGNES ET COLONNES DE LA MATRICE MORSE
      MNLIVO = 0
      MNLPVO = 0
C
C     FORMATION DE LA LISTE DES NO DES NOEUDS VOISINS DE CHAQUE NOEUD
C     ===============================================================
      CALL LISVOI( MNTOPO, MNNPEF, MNXYZN,
     %             MNLIVO, MNLPVO, MXVOIS, IERR  )
      IF( IERR .NE. 0 ) GOTO 9000
C
C     FORMATION DES TABLEAUX LPLIGN ET LPCOLO
C     =======================================
C
C     ADRESSAGE DU TABLEAU POINTEUR SUR LES LIGNES DES DL
      IF( MNLPLI .EQ. 0 ) THEN
         CALL TNMCDC( 'ENTIER', NTDL+1, MNLPLI )
      ENDIF
C
C     MNELE : ADRESSE DU TABLEAU NPEF"TYPE EF
      MNELE = MNNPEF(1)
C
C     LE NUMERO DU TYPE DE L'ELEMENT FINI
      NUTYEL = MCN( MNELE + WUTYEL )
C
C     LE NOMBRE DE COORDONNEES PAR NOEUD
      NBCOOR = MCN( MNXYZN + WBCOON )
C
C     LE POINTEUR SUR LES COLONNES DES DL DE LA MATRICE MORSE
C     MAJORATION DU NOMBRE DE COEFFICIENTS NON NULS SUR UNE LIGNE
      IF( NBCOOR .EQ. 6 ) THEN
C        6CUBES
         LOLPC0 = NTDL * 729 + 1
      ELSE
         IF( NUTYEL .LE. 15 ) THEN
            NDIM = 2
         ELSE
            NDIM = 3
         ENDIF
C        LE NOMBRE DE COLONNES * NB COEF D'UNE MATRICE NDIM+1 x NDIM+1
         LOLPC0 = MCN(MNLPVO+NBNOEUD) * (NDIM+1)**2
      ENDIF
      CALL TNMCDC( 'ENTIER', LOLPC0, MNLPCO )
C
      NBDLNO = 0
      MCN( MNLPLI ) = 0
      IA = 0
      IF( NCODSA .GT. 0 ) THEN
C
C        MATRICE SYMETRIQUE
C        ------------------
         MN1 = MNLIVO + MCN(MNLPVO)
         DO I=1,NBNOEUD
            MN2 = MNLIVO + MCN(MNLPVO+I) - 1
C
C           LES DEGRES DE LIBERTE DU NOEUD I DEBUTE EN
            NUDDLI = NDDLNO(I-1)
C
C           NOMBRE DE DL DU NOEUD I
            NBDLI = NDDLNO(I) - NUDDLI
C
C           CALCUL DU NOMBRE NBV DE VOISINS DU NOEUD I ET <I
            NBV = 0
            DO MNV = MN1, MN2
C              NUMERO DU NOEUD VOISIN
               NOVO = MCN(MNV)
               IF( NOVO .LT. I ) THEN
                  NBV = NBV + 1
               ENDIF
            ENDDO
C
            DO IN=1,NBDLI
C
C              LA LIGNE IN ASSOCIEE AU DL IN DU NOEUD I
C              LES NO DE COLONNES DES DL DE LA LIGNE IN DU NOEUD I
C
C              PARCOURS DES NOEUDS VOISINS DU NOEUD I SANS I
               DO MNV = MN1, MN1+NBV-1
C                 NUMERO DU NOEUD VOISIN
                  NOVO = MCN(MNV)
C
C                 NOMBRE DE DL DU NOEUD NOVO
                  NBDLV = NDDLNO(NOVO) - NDDLNO(NOVO-1)
C
C                 LES DEGRES DE LIBERTE DU NOEUD NOVO DEBUTE EN
                  NUDDLV = NDDLNO(NOVO-1)
                  DO  JN=1,NBDLV
C                    UNE COLONNE POUR LE DL JN
                     MCN(MNLPCO+IA) = NUDDLV + JN
                     IA = IA + 1
                  ENDDO
               ENDDO
C
C              LES IN DEGRES DE LIBERTE DU NOEUD I
C              LES DEGRES DE LIBERTE DU NOEUD I DEBUTE EN NUDDLI
               DO JN=1,IN
C                 UNE COLONNE POUR LE DL JN
                  MCN(MNLPCO+IA) = NUDDLI + JN
                  IA = IA + 1
               ENDDO
C
C              POINTEUR SUR LA DIAGONALE DE LA MATRICE
               MCN(MNLPLI+NUDDLI+IN) = IA
C
            ENDDO
C
C           PASSAGE AU NOEUD I+1 SUIVANT
            MN1 = MN2 + 1
         ENDDO
C
      ELSE
C
C        MATRICE NON SYMETRIQUE  (CF prgcmc.f)
C   A METTRE A JOUR LE NB NON CONSTANT DE DL PAR NOEUD
C        ----------------------
         MN1 = MNLIVO + MCN(MNLPVO)
         DO I=1,NBNOEUD
            MN2 = MNLIVO + MCN(MNLPVO+I) - 1
            IANVOI = 0
            NBV = MN2 - MN1 + 1
            DO MNV=MN1,MN2
               NOVO = MCN(MNV)
C              LES NBDLNO DEGRES DE LIBERTE DU NOEUD K
               NUDDLV = ( NOVO - 1 ) * NBDLNO
               DO IN=1,NBDLNO
                  DO JN=1,NBDLNO
                     IANV = IA + JN + IANVOI
     &               + ( NBV + 1 )  * NBDLNO * ( IN - 1 )
                     MCN(MNLPCO+IANV-1) = NUDDLV + JN
                  ENDDO
               ENDDO
               IANVOI = IANVOI + NBDLNO
            ENDDO
C           LES NBDLNO DEGRES DE LIBERTE DU NOEUD I
            NUDDLI = ( I - 1 ) * NBDLNO
            DO IN=1,NBDLNO
               KN = 0
               DO JN=1,NBDLNO
                  IF (JN.NE.IN) THEN
                     KN = KN + 1
                     IAN = IA + KN + MCN(MNLPLI+I) * NBDLNO
     &                + ( NBV + 1 )  * NBDLNO * ( IN - 1 )
                     MCN(MNLPCO+IAN-1) = NUDDLI + JN
                  END IF
               ENDDO
               IAN = IA + ( NBV + 1 ) * NBDLNO *  IN
               MCN(MNLPCO+IAN-1) = NUDDLI + IN
               MCN(MNLPLI+NUDDLI+IN) = IAN
            ENDDO
            IA  = IAN
            MN1 = MN2 + 1
         ENDDO
      ENDIF

C     DESTRUCTION DES TABLEAUX AUXILIAIRES
C     ------------------------------------
 9000 IF( MNLIVO .LE. 0 .OR. MNLPVO .LE. 0 ) RETURN
      LOLIVO = MCN(MNLPVO+NBNOEUD)
      CALL TNMCDS( 'ENTIER', LOLIVO,    MNLIVO )
      CALL TNMCDS( 'ENTIER', NBNOEUD+1, MNLPVO )

      LOLPCO = MCN(MNLPLI+NTDL)
      IF( LOLPC0 .LT. LOLPCO ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'prgcfl0: LOLPC0 A AUGMENTER'
         ELSE
            KERR(1) = 'prgcfl0: INCREASE the value of LOLPC0'
         ENDIF
         IERR = 1
         RETURN
      ENDIF

C     REDUCTION DE LA TAILLE DU TABLEAU DES NUMEROS DE COLONNES
      CALL TNMCRA( 'ENTIER', LOLPC0, LOLPCO, MNLPCO )
      RETURN
      END
