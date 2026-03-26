      SUBROUTINE JOINT1( NBSODO , NBTLIG , NBTNOE , NUTYOB , NUOBJE ,
     %                   XYZEXT , LIGNES , LINUOB , LITYOB , NBJOIN ,
     %          LINUNO , NBNOEU , COOR   ,  SENS  , MXNOJO , IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  RETROUVER LES LIGNES COMMUNES AUX SOUS-DOMAINES
C -----  POUR CREER LE TABLEAU DES JOINTS ( EN DIMENSION DEUX )
C        RESTRICTIONS D'UTILISATION :
C          1) LES SOUS-DOMAINES SONT CREES A PARTIR DE LIGNES
C             QUI SONT LEUR CONTOUR
C          2) LES LIGNES SONT SUPPOSEES AVOIR LES MEMES EXTREMITES
C          3) LES JOINTS NE SONT PAS REUNION DE PLUSIEURS LIGNES
C          4) UNIQUEMENT VALABLE POUR DES ELEMENTS FINIS DE DEGRE 1
C           ( POUR UN DEGRE > 1 TRAVAILLER SUR LE TABLEAU DES NOEUDS )
C
C ENTREES :
C ---------
C NBSODO : NOMBRE D'OBJETS A EXAMINER
C NBTLIG : MAJORANT DU NOMBRE DE NOEUDS DE L'INTERFACE
C NBTNOE : MAJORANT DU NOMBRE DE NOEUDS DE L'INTERFACE
C NUTYOB : LE NUMERO DU TYPE DE L'OBJET
C NUOBJE : LE NUMERO DE L'OBJET DANS LE LEXIQUE
C
C SORTIES :
C ---------
C XYZEXT : LES COORDONNEES DES EXTREMITES DES LIGNES
C LIGNES : LA LISTE DES NUMEROS DES LIGNES DU CONTOUR
C LINUOB : LES NUMEROS DES OBJETS ASSOCIES AUX LIGNES
C LITYOB : LES TYPES DES OBJETS ASSOCIES AUX LIGNES
C NBJOIN : NOMBRE DE JOINTS
C LINUNO : POUR CHAQUE NOEUD DE CHAQUE JOINT,
C          LE NUMERO DU NOEUD DANS L'OBJET ASSOCIE
C NBNOEU : NOMBRE DE NOEUDS DE CHAQUE JOINT
C COOR   : LES COORDONNEES DE CES NOEUDS
C MXNOJO : NOMBRE MAXIMUM DE NOEUDS DANS UN JOINT
C IERR   : 0 SI PAS D'ERREUR , > 0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY  ANALYSE NUMERIQUE UPMC PARIS  MAI 1993
C23456---------------------------------------------------------------012
      IMPLICIT          INTEGER (W)
      include"./incl/a_objet__definition.inc"
      include"./incl/a_objet__topologie.inc"
      include"./incl/a_ligne__definition.inc"
      include"./incl/gsmenu.inc"
      include"./incl/trvari.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/a___xyznoeud.inc"
      include"./incl/a___xyzsommet.inc"
      COMMON / EPSSSS / EPZERO,EPSXYZ
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (RMCN(1),MCN(1))
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
C
      INTEGER  NUTYOB(NBSODO),NUOBJE(NBSODO)
      INTEGER  LIGNES(NBTLIG,2),NBNOEU(NBTLIG,2)
      INTEGER  LINUOB(NBTLIG,2),LITYOB(NBTLIG,2)
      INTEGER  LINUNO(NBTLIG,2,NBTNOE),SENS(NBTLIG)
      REAL     COOR(NBTLIG,2,NBTNOE,3),XYZEXT(3,NBTLIG*2)
      REAL     XYZO1(3),XYZO2(3),XYZE1(3),XYZE2(3)
      CHARACTER*10  NMTYOB,KNOMTY
      CHARACTER*24  KNOM
C
C     INITIALISATIONS
C     ---------------
      IERR   = 0
      NBLICO = 0
C     LE NUMERO DU LEXIQUE DES LIGNES
      NUTYLI = 2
      NTLXLI = NTMN(NUTYLI)
C
C     BOUCLE SUR LES OBJETS SOUS-DOMAINES
C     -----------------------------------
C
      DO 50 NUSD=1,NBSODO
         NUTY = NUTYOB(NUSD)
         NUOB = NUOBJE(NUSD)
C        LE LEXIQUE DE CE TYPE D'OBJETS
         NTLX   = NTMN( NUTY )
C        LE LEXIQUE DE CET OBJET
         CALL LXNLOU( NTLX , NUOB , NTLXOB , NT )
C        RECHERCHE DU TABLEAU DEFINITION DE L'OBJET
         CALL LXTSOU( NTLXOB , 'DEFINITION' , NTDFOB , MNDFOB )
C        LE NOMBRE D'OBJETS COMPOSANT L'OBJET
         NBCOOB = MCN(MNDFOB+WBDOBJ)
C        LES OBJETS COMPOSANT L'OBJET : RECHERCHE DES LIGNES
         MN = MNDFOB + WTYOBJ - 1
         DO 51 NC = 1 , NBCOOB
            NUTYCO = MCN(MN+NC*2-1)
            NUMLIG = MCN(MN+NC*2)
            IF (NUTYCO.EQ.2) THEN
C              RECHERCHE DE CETTE LIGNE DANS LA LISTE ACTUELLE
               DO 52 NBL = 1 , NBLICO
                  IF (NUMLIG.EQ.LIGNES(NBL,1)) THEN
C                    CETTE LIGNE A DEJA ETE TROUVEE !
                     LINUOB(NBL,2) = NUOB
                     LITYOB(NBL,2) = NUTY
                     GOTO 51
                  ENDIF
52             CONTINUE
C              LA LIGNE N'EST PAS DANS LA LISTE
               NBLICO = NBLICO + 1
               LIGNES(NBLICO,1) = NUMLIG
               LIGNES(NBLICO,2) = 0
               LINUOB(NBLICO,1) = NUOB
               LITYOB(NBLICO,1) = NUTY
               LINUOB(NBLICO,2) = 0
               LITYOB(NBLICO,2) = 0
C              RECHERCHE DES EXTREMITES
C              LE LEXIQUE DE LA LIGNE
               CALL LXNLOU( NTLXLI , NUMLIG , NTLX0 , NT )
C              RECHERCHE DU TABLEAU SOMMETS DE LA LIGNE
               CALL LXTSOU( NTLX0 , 'DEFINITION' , NT , MNDFLI )
C              RECHERCHE DU TABLEAU SOMMETS DE LA LIGNE
               CALL LXTSOU( NTLX0 , 'XYZSOMMET' , NT , MNSOLI )
               IF (NT.NE.0) THEN
                   NBSOLI = MCN(MNSOLI+WNBSOM)
                   MNXYZ  = MNSOLI+WNBSOM
                   DO 53 K=1,3
                      XYZEXT(K,NBLICO*2-1)=RMCN(MNXYZ+K)
                      XYZEXT(K,NBLICO*2)=RMCN(MNXYZ+3*(NBSOLI-1)+K)
53                 CONTINUE
                   GOTO 51
               ELSE
                   NUTYLI = MCN(MNDFLI+WUTYLI)
C                  CE CAS N'EST PAS ENCORE PROGRAMME
                   KNOMTY = NMTYOB( NUTYLI )
                   CALL NMOBNU( KNOMTY , NUMLIG , KNOM )
                   NBLGRC(NRERR) = 2
                   KERR(1) = ' ATTENTION : LIGNE '// KNOM
                   KERR(2) = ' SANS TABLEAU XYZSOMMET'
                   CALL LEREUR
                   IERR = 0
                   WRITE(IMPRIM,5000) KNOM
C                  CE CAS N'EST PAS ENCORE RENCONTRE !
               ENDIF
            ENDIF
51       CONTINUE
C
50    CONTINUE
C
C     RECHERCHE DES LIGNES AYANT MEME EXTREMITES
C     ------------------------------------------
C
      DO 60 NB1 = 1 , NBLICO
         IF (LINUOB(NB1,2).EQ.0) THEN
            DO 61 K=1,3
               XYZO1(K)=XYZEXT(K,NB1*2-1)
               XYZE1(K)=XYZEXT(K,NB1*2)
61          CONTINUE
            DO 62 NB2 = NB1+1 , NBLICO
               IF (LINUOB(NB2,2).EQ.0) THEN
C                 COMPARAISON DES EXTREMITES
                  DO 63 K=1,3
                     XYZO2(K)=XYZEXT(K,NB2*2-1)
                     XYZE2(K)=XYZEXT(K,NB2*2)
63                CONTINUE
                  CALL XYZIDE(XYZO1,XYZO2,IDENT)
                  IF (IDENT.EQ.1) THEN
                     CALL XYZIDE(XYZE1,XYZE2,IDENT)
                     IF (IDENT.NE.1) GOTO 62
C                    LES LIGNES ONT MEMES EXTREMITES !
                     SENS(NB1) = 1
                     LIGNES(NB1,2) = LIGNES(NB2,1)
                     LINUOB(NB1,2) = LINUOB(NB2,1)
                     LITYOB(NB1,2) = LITYOB(NB2,1)
                     GOTO 60
                  ELSE
C                    CHANGEMENT DE SENS
                     CALL XYZIDE(XYZE1,XYZO2,IDENT)
                     IF (IDENT.EQ.1) THEN
                        CALL XYZIDE(XYZO1,XYZE2,IDENT)
                        IF (IDENT.NE.1) GOTO 62
C                       LES LIGNES ONT MEMES EXTREMITES !
                        SENS(NB1) = -1
                        LIGNES(NB1,2) = LIGNES(NB2,1)
                        LINUOB(NB1,2) = LINUOB(NB2,1)
                        LITYOB(NB1,2) = LITYOB(NB2,1)
                        GOTO 60
                     ENDIF
                  ENDIF
               ENDIF
62          CONTINUE
         ENDIF
60    CONTINUE
      NBJOIN = 0
      DO 64 NBL = 1 , NBLICO
         IF (LINUOB(NBL,2).NE.0) THEN
            NBJOIN = NBJOIN + 1
            LIGNES(NBJOIN,1) = LIGNES(NBL,1)
            LIGNES(NBJOIN,2) = LIGNES(NBL,2)
            LINUOB(NBJOIN,1) = LINUOB(NBL,1)
            LINUOB(NBJOIN,2) = LINUOB(NBL,2)
            LITYOB(NBJOIN,1) = LITYOB(NBL,1)
            LITYOB(NBJOIN,2) = LITYOB(NBL,2)
            IF (SENS(NBL) .NE. 0) THEN
               SENS(NBJOIN)  = SENS(NBL)
            ELSE
               SENS(NBJOIN)  = 1
            ENDIF
         ENDIF
64    CONTINUE
C
C     CONSTRUCTION DES JOINTS
C     -----------------------
C
      MXNOJO = 0
      DO 70 NBL = 1 , NBJOIN
C        LE NUMERO DE LA LIGNE COTE 1
         NUMLIG = LIGNES(NBL,1)
C        LE LEXIQUE DE LA LIGNE
         CALL LXNLOU( NTLXLI , NUMLIG , NTLX0 , NT )
C        RECHERCHE DU TABLEAU SOMMETS DE LA LIGNE
         CALL LXTSOU( NTLX0 , 'XYZSOMMET' , NT , MNSOLI )
         IF (NT.NE.0) THEN
             NOEUDS = MCN(MNSOLI+WNBSOM)
             IF (MXNOJO.LT.NOEUDS) MXNOJO = NOEUDS
             NBNOEU(NBL,1) = NOEUDS
             MNXYZ = MNSOLI + WYZSOM - 4
             DO 71 NBS = 1 , NOEUDS
                MNXYZ  = MNXYZ + 3
                DO 72 K=1,3
                   COOR(NBL,1,NBS,K) = RMCN(MNXYZ+K)
72              CONTINUE
71          CONTINUE
         ENDIF
C        LE NUMERO DE LA LIGNE COTE 2
         IF (LIGNES(NBL,2).EQ.0) LIGNES(NBL,2) = LIGNES(NBL,1)
         NUMLIG = LIGNES(NBL,2)
C        LE LEXIQUE DE LA LIGNE
         CALL LXNLOU( NTLXLI , NUMLIG , NTLX0 , NT )
C        RECHERCHE DU TABLEAU SOMMETS DE LA LIGNE
         CALL LXTSOU( NTLX0 , 'XYZSOMMET' , NT , MNSOLI )
         IF (NT.NE.0) THEN
             NOEUDS = MCN(MNSOLI+WNBSOM)
             IF (MXNOJO.LT.NOEUDS) MXNOJO = NOEUDS
             NBNOEU(NBL,2) = NOEUDS
             MNXYZ = MNSOLI + WYZSOM - 4
             DO 73 NBS = 1 , NOEUDS
                MNXYZ  = MNXYZ + 3
                IF (SENS(NBL).EQ.1) THEN
                   NBSO = NBS
                ELSE
                   NBSO = NOEUDS + 1 - NBS
                ENDIF
                DO 74 K=1,3
                   COOR(NBL,2,NBSO,K) = RMCN(MNXYZ+K)
74              CONTINUE
73          CONTINUE
         ENDIF
70    CONTINUE
C
C     RECHERCHE DES NUMEROS DES NOEUDS
C     --------------------------------
      DO 80 NUSD=1,NBSODO
         NUTY = NUTYOB(NUSD)
         NUOB = NUOBJE(NUSD)
C        LE LEXIQUE DE CE TYPE D'OBJETS
         NTLX   = NTMN( NUTY )
C        LE LEXIQUE DE CET OBJET
         CALL LXNLOU( NTLX , NUOB , NTLXOB , NT )
C        LE NOM DE CET OBJET
         KNOMTY = NMTYOB( NUTY )
         CALL NMOBNU( KNOMTY , NUOB , KNOM )
C        RECUPERATION DU TABLEAU NOEUDS
         CALL LXTSOU( NTLXOB , 'XYZNOEUD'  , NTNOEU , MNNOEU )
         IF (NTNOEU.NE.0) THEN
            NBNOOB = MCN(MNNOEU+WNBNOE)
            MNXYZO = MNNOEU+WYZNOE-1
         ELSE
C            WRITE(IMPRIM,4000) KNOM
C           RECUPERATION DU TABLEAU SOMMETS
            CALL LXTSOU( NTLXOB , 'XYZSOMMET' , NTSOMM , MNSOMM )
            NBNOOB = MCN(MNSOMM+WNBSOM)
            MNXYZO = MNSOMM+WYZSOM-1
         ENDIF
C        RECHERCHE DES JOINTS LIES A L'OBJET
         DO 81 NBL = 1 , NBJOIN
         DO 82 JI  = 1 , 2
             NUMOB = LINUOB(NBL,JI)
             IF (NUMOB.NE.NUOB) GOTO 82
C            RECHERCHE DES NUMEROS DES NOEUDS DANS L'OBJET NUOB
             NBNOLI = NBNOEU(NBL,JI)
             DO 83 NBS=1,NBNOLI
                DO 84 K=1,3
                   XYZO1(K) = COOR(NBL,JI,NBS,K)
84              CONTINUE
                MNXYZ=MNXYZO
                DO 85 NBT=1,NBNOOB
                   DO 86 K=1,3
                      XYZO2(K)=RMCN(MNXYZ+K)
86                 CONTINUE
                   CALL XYZIDE(XYZO1,XYZO2,IDENT)
                   MNXYZ=MNXYZ+3
                   IF (IDENT.NE.1) GOTO 85
C                  LE NOEUD EST RETROUVE !
                   LINUNO(NBL,JI,NBS) = NBT
                   GOTO 83
85              CONTINUE
83           CONTINUE
82        CONTINUE
81        CONTINUE
80    CONTINUE
C
C     IMPRESSIONS
C     -----------
C      DO 90 NB = 1 , NBJOIN
C         WRITE (IMPRIM,1000) NB,LIGNES(NB,1),LIGNES(NB,1),
C     S                          LITYOB(NB,1),LITYOB(NB,2),
C     S                          LINUOB(NB,1),LINUOB(NB,2),
C     S                          NBNOEU(NB,1),NBNOEU(NB,2)
C         DO 91 JI = 1 , 2
C            WRITE (IMPRIM,2000) JI
C            NBN=NBNOEU(NB,JI)
C            DO 92 NJI = 1 , NBN
C               WRITE (IMPRIM,3000) LINUNO(NB,JI,NJI),
C     S                      (COOR(NB,JI,NJI,K),K=1,3)
C  92        CONTINUE
C  91     CONTINUE
C  90  CONTINUE
C
      RETURN
ccc1000  FORMAT(//,1X,80(1H-)/,' JOINT',I3,' LIGNES',2I5,
ccc     S ' OBJETS',4I5,' NOEUDS',2I5/,1X,80(1H-))
ccc2000  FORMAT(/,' COTE',I4/)
ccc3000  FORMAT(' NUMERO',I5,' COOR',3D15.6)
ccc4000  FORMAT(/,' SP JOINT1 : OBJET ',A24,/,' SANS TABLEAU XYZNOEUD')
5000  FORMAT(/,'SP JOINT1 : LIGNE ',A24,/,' SANS TABLEAU XYZSOMMET')
      END
