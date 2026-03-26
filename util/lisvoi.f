      SUBROUTINE LISVOI( MNTOPO, MNNPEF, MNXYZN,
     %                   MNLIVO, MNLPVO, MXVOIS, IERR  )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   GENERER LA LISTE DES VOISINS DE CHACUN DES NOEUDS DU MAILLAGE
C ----- ( LE NOEUD N'EST PAS PRIS EN COMPTE DANS SA LISTE )
C
C ENTREES:
C --------
C MNTOPO : ADRESSE  MCN DU  TABLEAU  TOPOLOGIE DE L'OBJET
C MNNPEF : ADRESSES MCN DES TABLEAUX NPEF"     DE L'OBJET
C MNXYZN : ADRESSE  MCN DU  TABLEAU  XYZNOEUD  DE L'OBJET
C
C SORTIES:
C --------
C MNLIVO : ADRESSE MCN DU TABLEAU LIVO LISTE DES VOISINS DES NOEUDS
C          ATTENTION: LE NOEUD NE FAIT PAS PARTIE DES VOISINS
C MNLPVO : ADRESSE MCN DU POINTEUR SUR LE DERNIER NUMERO DES VOISINS
C          DE CHACUN DES NOEUDS
C          MCN(MNLPVO+I)=INDICE DANS LISTE VOISINS DU DERNIER NO
C                        DE VOISIN DU NOEUD I
C          MCN(MNLPVO+0)=0
C          =0 SI PAS ASSEZ DE MEMOIRE MCN
C          CES TABLEAUX SONT DECLARES ET REMPLIS PAR LISVOI
C MXVOIS : NOMBRE MAXIMAL DE VOISINS D'UN NOEUD DU MAILLAGE
C IERR   : 0 SI PAS D'ERREUR
C          1 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY  ANALYSE NUMERIQUE UPMC PARIS       OCTOBRE  1989
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/donele.inc"
      include"./incl/a___xyznoeud.inc"
      include"./incl/a___npef.inc"
      include"./incl/a_objet__topologie.inc"
      include"./incl/pp.inc"
      COMMON   MCN(MOTMCN)
      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)
      INTEGER           MNNPEF(1:*)
      DOUBLE PRECISION  D, DBLE

      IERR   = 0
      MNLIVO = 0
      MNLPVO = 0
      MNNOEF = 0
      MXVOIS = 0

C     LE NOMBRE TOTAL DE NOEUDS ET DE DEGRES DE LIBERTE DE L'OBJET
      NBNOE = MCN( MNXYZN + WNBNOE )

C     NOMBRE DE COORDONNEES D'UN NOEUD
      NBCOOR = MCN( MNXYZN + WBCOON )

C     NOMBRE MAXIMUM DE VOISINS D'UN NOEUD DONNE
CCC   modif pour les 6-CUBES  MAXVOI = MXNOEL * 10 devient
      IF( NBCOOR .EQ. 6 ) THEN
         MAXVOI = MXNOEL * 64 / 2
      ELSE
ccc         MAXVOI = MXNOEL * 10 / 2  2/2/2021
         MAXVOI = MXNOEL * 6 / 2
      ENDIF
      NBAGMX = 0

C     ADRESSAGE DES TABLEAUX LIVO ET LPVO
C     -----------------------------------
 10   NBAGMX = NBAGMX + 1
      IF( NBAGMX .GT. 5 ) THEN
C        TROP D'AUGMENTATION DE LA VALEUR DE MAXVOI
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'TROP d''AUGMENTATIONS de MAXVOI dans lisvoi.f'
         ELSE
            KERR(1) = 'TOO AUGMENTATIONS of MAXVOI in lisvoi.f'
         ENDIF
         CALL LEREUR
         IERR = 1
         GOTO 9000
      ELSE
         MAXVOI = MAXVOI * 2
      ENDIF

C     LE POINTEUR ASSOCIE
      CALL TNMCDC( 'ENTIER', NBNOE+1, MNLPVO )
C     LE TABLEAU DES NOEUDS D'UN ELEMENT FINI
      CALL TNMCDC( 'ENTIER', MXNOEL, MNNOEF )

C     LE TABLEAU DES NOEUDS VOISINS
      CALL TNMCMX( 'ENTIER', MAXVAR )
      MOTSNE = NBNOE*MAXVOI + NBNOE + 1 + MXNOEL
      IF( MAXVAR .LT. MOTSNE ) THEN
         NBLGRC(NRERR) = 6
         WRITE(KERR(7),'(I12)') MOTSNE
         WRITE(KERR(8),'(I12)') MAXVAR
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'ERREUR: SUPER TABLEAU MCN A AUGMENTER'
            KERR(2) = 'POUR L''EXECUTION du sp lisvoi.f'
            KERR(3) = KERR(7)(1:12) // ' MOTS MEMOIRE DEMANDES'
            KERR(4) = KERR(8)(1:12) // ' MOTS DISPONIBLES dans MCN'
            KERR(5) = 'ou REDUIRE LE NOMBRE DE NOEUDS'
            KERR(6) = 'ou LE NOMBRE DE NOEUDS VOISINS D''UN NOEUD'
         ELSE
            KERR(1) = 'ERROR: AUGMENT the SIZE of SUPER-ARRAY MCN'
            KERR(2) = 'for the EXECUTION of the subroutine lisvoi.f'
            KERR(3) = KERR(7)(1:12) // ' WISHED MEMORY WORDS'
            KERR(4) = KERR(8)(1:12) // ' WORDS USABLE in MCN'
            KERR(5) = 'or REDUCE the NUMBER of NODES'
            KERR(6) = 'or the NUMBER of NEIGHBORING NODES of a NODE'
         ENDIF
         CALL LEREUR
         IERR = 1
         CALL TNMCDS( 'ENTIER', MXNOEL,  MNNOEF )
         CALL TNMCDS( 'ENTIER', NBNOE+1, MNLPVO )
         GOTO 9000
      ENDIF

C     DECLARATION DES TABLEAUX DES VOISINS
      MXCOEF = NBNOE * MAXVOI
      IF( MXCOEF .LE. 0 ) THEN
         NBLGRC(NRERR) = 3
         D = DBLE(NBNOE) * DBLE(MAXVOI)
         WRITE(KERR(4),'(D25.17)') D
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'ERREUR: TAILLE DE LA LISTE DES VOISINS=' //
     %                 KERR(4)(1:25)
            KERR(2) = 'DEPASSEMENT DU MAX DES ENTIERS > 2**31'
            KERR(3) = 'REDUIRE LA TAILLE DU MAILLAGE'
         ELSE
            KERR(1) = 'ERROR: NEIGHBORS LIST SIZE=' //
     %                 KERR(4)(1:25)
            KERR(2) = 'OVERFLOW of INTEGER > 2**31'
            KERR(3) = 'REDUCE THE SIZE OF THE MESH'
         ENDIF
         CALL LEREUR
         IERR = 1
         GOTO 9000
      ENDIF
      CALL TNMCDC( 'ENTIER', MXCOEF, MNLIVO )
      IF( MNLIVO .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'LISVOI: PAS ASSEZ DE MEMOIRE MCN'
         ELSE
            KERR(1) = 'LISVOI: NOT ENOUGH MEMORY MCN'
         ENDIF
         CALL LEREUR
         IERR = 1
         GOTO 90
      ENDIF
      CALL AZEROI( MXCOEF, MCN(MNLIVO) )

C     LA BOUCLE SUR LES TYPES D'ELEMENTS FINIS DU MAILLAGE
C     ====================================================
      NBTYEL = MCN( MNTOPO + WBTYEL )
      DO 19 NUTYEL=1,NBTYEL

C        LE TABLEAU NPEF"   A POUR ADRESSE MNELE
         MNELE = MNNPEF( NUTYEL )
         IF( MNELE .LE. 0 ) THEN
            NBLGRC(NRERR) = 1
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'ERREUR: TYPE INCONNU D''ELEMENT FINI'
            ELSE
               KERR(1) = 'ERROR: UNKNOW TYPE OF FINITE ELEMENT'
            ENDIF
            CALL LEREUR
            IERR = 2
            GOTO 19
         ENDIF

C        LA BOUCLE SUR LES ELEMENTS FINIS DE CE TYPE
C        ===========================================
         NBELEM = MCN( MNELE + WBELEM )
         NBNDEL = MCN( MNELE + WBNDEL )
         MNNDEL = MNELE + WUNDEL - 1
         DO NUELEM=1,NBELEM
            MNNDE = MNNDEL + NUELEM

C           LES NUMEROS DES NOEUDS DE L'ELEMENT FINI
C           ----------------------------------------
            MN = MNNDE
            DO I=0,NBNDEL-1
               MCN(MNNOEF+I) = MCN(MN)
               MN = MN + NBELEM
            ENDDO

C           RECHERCHE DANS LA LISTE DES VOISINS
C           -----------------------------------
            DO I=0,NBNDEL-1
C              LE NO DU NOEUD I DE L'EF
               NUMI = MCN(MNNOEF+I)
               IAVI = MAXVOI * ( NUMI - 1 ) - 1
               DO 15 J=0,NBNDEL-1
                  IF ( I .EQ. J ) GOTO 15
                  NUMJ = MCN(MNNOEF+J)
                  DO K=1,MAXVOI
                     IAK = MNLIVO + IAVI + K
                     KV  = MCN(IAK)
                     IF (KV.EQ.0) THEN
C                        ON AJOUTE CE NOEUD A LA LISTE
                         MCN(IAK) = NUMJ
                         GOTO 15
                     ELSE  IF (KV.EQ.NUMJ) THEN
C                        LE NOEUD EST DEJA DANS LA LISTE
                         GOTO 15
                     ENDIF
                  ENDDO

C                 IL N'Y A PLUS DE PLACE, ET NUMJ N'EST PAS DANS LA LISTE !
                  NBLGRC(NRERR) = 2
                  WRITE(KERR(MXLGER)(1:15),'(I15)') MAXVOI
                  IF( LANGAG .EQ. 0 ) THEN
                     KERR(1) = 'MAXVOI TROP PETIT '//KERR(MXLGER)(1:15)
                     KERR(2) = 'MAXVOI est AUGMENTE dans le sp lisvoi.f'
                  ELSE
                    KERR(1)='MAXVOI TOO SMALL '//KERR(MXLGER)(1:15)
                    KERR(2)='MAXVOI is AUGMENTED in subroutine lisvoi.f'
                  ENDIF
                  CALL LERESU
C                 DESTRUCTION DES TABLEAUX TROP PETITS
                  CALL TNMCDS( 'ENTIER', MXNOEL,  MNNOEF )
                  CALL TNMCDS( 'ENTIER', NBNOE+1, MNLPVO )
                  CALL TNMCDS( 'ENTIER', MXCOEF,  MNLIVO )
                  GOTO 10

 15            ENDDO
            ENDDO
         ENDDO
 19   ENDDO
      IF( IERR .NE. 0 ) GOTO 90

C     CLASSEMENT DES NOEUDS VOISINS PAR ORDRE CROISSANT
C     COMPRESSION DU TABLEAU LIVO, CREATION DU POINTEUR
C     =================================================
      MCN(MNLPVO)=0
      IC0 = 0
      IC1 = 0
      DO I=1,NBNOE
         IA = MNLIVO + ( I - 1 ) * MAXVOI - 1
C        CLASSEMENT CROISSANT DES VOISINS
         MNLPVI = MNLPVO + I
         MCN(MNLPVI) = 0
         DO 21 K1=1,MAXVOI
            IA1 = IA + K1
            IF( MCN(IA1) .EQ. 0 ) GOTO 23
            MCN(MNLPVI) = MCN(MNLPVI) + 1
            DO K2=K1+1,MAXVOI
               IA2 = IA + K2
               IF( MCN(IA2) .EQ. 0 ) GOTO 21
               IF( MCN(IA2) .LT. MCN(IA1) ) THEN
                  MCNIA1  =MCN(IA1)
                  MCN(IA1)=MCN(IA2)
                  MCN(IA2)=MCNIA1
               ENDIF
            ENDDO
 21      ENDDO
C        COMPRESSION DU TABLEAU
 23      DO K1=1,MCN(MNLPVI)
            MCN(MNLIVO+IC1) = MCN(IA+K1)
            IC1 = IC1 + 1
         ENDDO
         MCN(MNLPVI) = IC1
C        NOMBRE DE VOISINS DU NOEUD I
         ID = IC1 - IC0
         IF( ID .GT. MXVOIS ) MXVOIS=ID
         IC0 = IC1
      ENDDO

C     GESTION DES TABLEAUX
C     --------------------
 90   IF( MNNOEF .GT. 0 ) CALL TNMCDS( 'ENTIER', MXNOEL, MNNOEF )
      IF( MNLIVO .GT. 0 ) THEN
         LOLIVO = MCN(MNLPVO+NBNOE)
         CALL TNMCRA( 'ENTIER', MXCOEF, LOLIVO, MNLIVO )
      ENDIF

      D = IC1
      D = D / NBNOE
      IF( LANGAG .EQ. 0 ) THEN
        WRITE(IMPRIM,*) 'NOMBRE MAXIMAL DE VOISINS D''UN NOEUD=',MXVOIS
        WRITE(IMPRIM,*) 'NOMBRE MOYEN   DE VOISINS D''UN NOEUD=',NINT(D)
        WRITE(IMPRIM,*) 'NOMBRE TOTAL   DE VOISINS DES NOEUDS=',IC1
      ELSE
        WRITE(IMPRIM,*) 'NEIGHBORS NUMBER MAXIMUM of a NODE=',MXVOIS
        WRITE(IMPRIM,*) 'NEIGHBORS NUMBER AVERAGE of a NODE=',NINT(D)
        WRITE(IMPRIM,*) 'NEIGHBORS TOTAL  NUMBER  of  NODES=',IC1
      ENDIF

 9000 RETURN
      END
