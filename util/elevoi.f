      SUBROUTINE ELEVOI ( MNTOPO , MNELEM , MNNOEU ,
     %                    MNLPVO , MNLIVO , NBEL   , IERR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : GENERER LA LISTE DES VOISINS DE CHACUN DES EF DU MAILLAGE
C ----- ( L'EF N'EST PAS PRIS EN COMPTE DANS SA LISTE )
C
C ENTREES :
C ---------
C MNTOPO : ADRESSE  MCN DU TABLEAU TOPOLOGIE  DE L'OBJET
C MNELEM : ADRESSES MCN DES TABLEAUX ELEMENTS DE L'OBJET
C MNNOEU : ADRESSE  MCN DU TABLEAU NOEUDS     DE L'OBJET
C NBEL   : NOMBRE TOTAL D'ELEMENTS DU MAILLAGE
C
C SORTIES :
C ---------
C MNLPVO : ADRESSE MCN DU POINTEUR ASSOCIE AU TABLEAU PRECEDENT
C MNLIVO : ADRESSE MCN DU TABLEAU LIVO LISTE DES VOISINS DES ELEMENTS
C          CES TABLEAUX SONT CREES PAR ELEVOI
C NBEL   : NOMBRE D'ELEMENTS TRAITES DU MAILLAGE
C IERR   : 0 SI PAS D'ERREUR
C          1 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY  ANALYSE NUMERIQUE UPMC PARIS  JUIN 1990
C23456---------------------------------------------------------------012
      IMPLICIT        INTEGER (W)
      include"./incl/gsmenu.inc"
      include"./incl/donele.inc"
      include"./incl/a___xyznoeud.inc"
      include"./incl/a___npef.inc"
      include"./incl/a_objet__topologie.inc"
      include"./incl/pp.inc"
      COMMON          MCN(MOTMCN)
      COMMON /UNITES/ LECTEU,IMPRIM,NUNITE(30)
      INTEGER     MNELEM(1:*)
C
      IERR = 0
C
C     LE NOMBRE TOTAL DE NOEUDS ET DE DEGRES DE LIBERTE DE L'OBJET
      NBNOE  = MCN( MNNOEU + WNBNOE )
C     NOMBRE MAXIMUM DE VOISINS D'UN NOEUD
      MAXVON = MXNOEL
C     NOMBRE MAXIMUM DE VOISINS D'UN ELEMENT
      MAXVOE = MXNOEL
C
C     ADRESSAGE DES TABLEAUX
C     ----------------------
      NBELEM = NBEL
C     LE TABLEAU DES ELEMENTS VOISINS
      CALL TNMCDC( 'ENTIER' , NBELEM*MAXVOE , MNLIVE )
C     LE POINTEUR ASSOCIE
      CALL TNMCDC( 'ENTIER' , NBELEM+1 , MNLPVO )
C     LE TABLEAU DES ELEMENTS CONTENANT CHAQUE NOEUD
      CALL TNMCDC( 'ENTIER' , NBNOE*MAXVON , MNLINO )
C     MISE A ZERO
      CALL AZEROI( NBELEM+1       , MCN(MNLPVO) )
      CALL AZEROI( NBELEM*MAXVOE  , MCN(MNLIVE) )
      CALL AZEROI( NBNOE*MAXVON   , MCN(MNLINO) )
C
C     LA BOUCLE SUR LES TYPES D'ELEMENTS DU MAILLAGE
C     ==============================================
      NBTYEL = MCN( MNTOPO + WBTYEL )
      NBELE  = 0
      DO 1 NUTYEL=1,NBTYEL
C
C        LE TABLEAU DES ELEMENTS A POUR ADRESSE MNELE
         MNELE = MNELEM( NUTYEL )
         IF( MNELE .LE. 0 ) THEN
            NBLGRC(NRERR) = 1
            KERR(1) =  ' ERREUR : ELEMENTS INCONNUS'
            CALL LEREUR
            IERR = 1
            GOTO 1
         ENDIF
C
C        LA BOUCLE SUR LES ELEMENTS DE CE TYPE
         NBELEM = MCN( MNELE + WBELEM )
         NBNDEL = MCN( MNELE + WBNDEL )
         MNNDEL = MNELE + WUNDEL - 1
         DO 2 NUELEM=1,NBELEM
            NBELE  = NBELE  + 1
            MNNDE  = MNNDEL + NUELEM
            MN     = MNNDE
C           LES NUMEROS DES NOEUDS DE L'ELEMENT NBELE
            DO 3 I=0,NBNDEL-1
               NOE = MCN(MN)
               MN  = MN + NBELEM
C              AJOUT DE NBELE DANS LA LISTE DES ELEMENTS CONTENANT NOE
               IA = MNLINO + NOE - 1
               DO 4 K=0,MAXVON-1
                  KV  = MCN(IA)
                  IF (KV.EQ.0) THEN
C                    ON AJOUTE CET ELEMENT A LA LISTE
                     MCN(IA) = NBELE
                     GO TO 3
                  END IF
                  IA = IA + NBNOE
 4             CONTINUE
C              IL N'Y A PLUS DE PLACE, ET NBELE N'EST PAS DANS LA LISTE !
               NBLGRC(NRERR) = 1
               WRITE(KERR(MXLGER)(1:4),'(I4)') MAXVON
               KERR(1) = 'MAXVON TROP PETIT '//KERR(MXLGER)(1:4)
               CALL LEREUR
               IERR = 1
 3          CONTINUE
 2       CONTINUE
 1    CONTINUE
      IF( IERR .NE. 0 ) GOTO 100
C
C     CREATION DU TABLEAU DES ELEMENTS VOISINS
C     ========================================
C     LE TABLEAU DES ELEMENTS CONTENANT UN NOEUD
      CALL TNMCDC( 'ENTIER' , MAXVON , MNNUEL )
      CALL AZEROI( MAXVON   , MCN(MNNUEL) )
      MN = MNLINO
      DO 5 NOE = 1 , NBNOE
C        LE NOEUD NOE
         IAN = MN
         MN  = MN + 1
C        LES ELEMENTS QUI CONTIENNENT LE NOEUD NOE
         NBELCN = 0
         DO 6 K=0,MAXVON-1
            IF (MCN(IAN).NE.0) THEN
               NBELCN = NBELCN + 1
               MCN(MNNUEL+K) = MCN(IAN)
               IAN  = IAN + NBNOE
            ELSE
               GO TO 7
            END IF
 6       CONTINUE
 7       DO 8 K=0,NBELCN-1
            NUEL = MCN(MNNUEL+K)
C           LES VOISINS DE L'ELEMENT NUEL
            IAE = MNLIVE + NUEL -1
            DO 9 L=0,NBELCN-1
               IF (L.NE.K) THEN
                  NUEV = MCN(MNNUEL+L)
                  IAEL = IAE
                  DO 10 LL=0,MAXVOE-1
                     NUEW = MCN(IAEL)
                     IF (NUEW.EQ.0) THEN
                        MCN(IAEL) = NUEV
                        GO TO 9
                     ELSE IF (NUEV.EQ.NUEW) THEN
                        GO TO 9
                     ELSE
                        IAEL = IAEL + NBEL
                     END IF
 10               CONTINUE
C                 IL N'Y A PLUS DE PLACE, ET NUEV N'EST PAS DANS LA LISTE !
                  NBLGRC(NRERR) = 1
                  WRITE(KERR(MXLGER)(1:4),'(I4)') MAXVOE
                  KERR(1) = 'MAXVOE TROP PETIT '//KERR(MXLGER)(1:4)
                  CALL LEREUR
                  IERR = 1
               END IF
 9          CONTINUE
 8       CONTINUE
 5    CONTINUE
C
C     CLASSEMENT DES ELEMENTS VOISINS PAR ORDRE CROISSANT
C     COMPRESSION DU TABLEAU LIVO, CREATION DU POINTEUR
C     ===================================================
      MCN(MNLPVO)=0
      IC = 0
      IA = MNLIVE
      DO 15 NE=1,NBELE
C        CLASSEMENT DES VOISINS PAR ORDRE CROISSANT
         IA1 = IA
         NBELV = 0
         DO 11 K1=0,MAXVOE-1
            IF (MCN(IA1).EQ.0) GO TO 12
            NBELV = NBELV + 1
            IA1 = IA1 + NBEL
 11      CONTINUE
 12      IA1 = IA
         DO 13 K1=1,NBELV
            IA2 = IA1+NBEL
            DO 14 K2=K1+1,NBELV
               IF (MCN(IA2).LT.MCN(IA1)) THEN
                  MCNIA1=MCN(IA1)
                  MCN(IA1)=MCN(IA2)
                  MCN(IA2)=MCNIA1
               END IF
               IA2 = IA2 + NBEL
 14         CONTINUE
            IA1 = IA1 + NBEL
 13      CONTINUE
         IC = IC + NBELV
         MCN(MNLPVO+NE) = IC
         IA = IA + 1
 15   CONTINUE
      LOLIVO = MCN(MNLPVO+NBELE)
C     LE TABLEAU DES ELEMENTS VOISINS
      CALL TNMCDC( 'ENTIER' , LOLIVO , MNLIVO )
      L1 = MCN(MNLPVO)
      IA = MNLIVE
      DO 16 NE=1,NBELE
         L2 = MCN(MNLPVO+NE)
         IAE = IA
         DO 17 L=L1,L2-1
            MCN(MNLIVO+L) = MCN(IAE)
            IAE = IAE + NBEL
 17      CONTINUE
         IA = IA + 1
         L1 = L2
 16   CONTINUE
C
C     GESTION DES TABLEAUX
C     --------------------
 100  CALL TNMCDS( 'ENTIER' , NBNOE*MAXVON  , MNLINO )
      CALL TNMCDS( 'ENTIER' , NBELEM*MAXVOE , MNLIVE )
      CALL TNMCDS( 'ENTIER' , MAXVON        , MNNUEL )
      NBEL = NBELE
      RETURN
      END
