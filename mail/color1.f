      SUBROUTINE COLOR1( NBEL   , NBCOMX , LPVOIS , LISVOI ,
     %                   LPCOLO , ELECOU , COUELE , ELELOG ,
     %                   COLLOG , RENUEL ,  IERR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : COLORIER LES ELEMENTS D'UN MAILLAGE
C -----
C
C ENTREES :
C ---------
C NBEL   : NOMBRE D'ELEMENTS DU MAILLAGE
C NBCOMX : NOMBRE MAXIMUM DE COULEURS DU MAILLAGE
C LPVOIS : POINTEUR SUR LES VOISINS D'UN ELEMENT
C LISVOI : LISTE DES VOISINS D'UN ELEMENT
C
C SORTIES :
C ---------
C LPCOLO : POINTEUR SUR LES COULEURS
C ELECOU : LISTE DES ELEMENTS DE COULEUR DONNEE
C COUELE : COULEUR D'UN ELEMENT
C ELELOG : TABLEAU LOGIQUE UTILITAIRE
C COLLOG : TABLEAU LOGIQUE UTILITAIRE
C RENUEL : RENUMEROTATION LOCALE DES ELEMENTS
C IERR   : CODE D'ERREUR
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY  ANALYSE NUMERIQUE UPMC PARIS  JUIN 1990
C23456---------------------------------------------------------------012
      include"./incl/gsmenu.inc"
      COMMON /UNITES/ LECTEU,IMPRIM,NUNITE(30)
      COMMON / MSIMTA / NOIMPR
      INTEGER     LPVOIS(NBEL+1),LISVOI(*),LPCOLO(NBCOMX+1),
     &            ELECOU(*),COUELE(NBEL),RENUEL(NBEL)
      LOGICAL     ELELOG(NBEL),COLLOG(NBEL),CHGT
C
C     INITIALISATION
C     --------------
      IERR = 0
      NBCO = 1
      DO 1 NE = 1 , NBEL
         COUELE(NE) =  NBCO
         ELELOG(NE) = .TRUE.
         COLLOG(NE) = .FALSE.
 1    CONTINUE
C
C     RENUMEROTATION DES ELEMENTS
C     ---------------------------
C
      NOEL = 1
      RENUEL(NOEL) = 1
      COLLOG(1)    = .TRUE.
      DO 11 NEA = 1 , NBEL
         NE = RENUEL(NEA)
         L1 = LPVOIS(NE)+1
         L2 = LPVOIS(NE+1)
         DO 12 L = L1 , L2
            NUELVO = LISVOI(L)
            IF ( .NOT. COLLOG(NUELVO) ) THEN
                 NOEL = NOEL + 1
                 RENUEL(NOEL)   = NUELVO
                 COLLOG(NUELVO) = .TRUE.
            END IF
 12      CONTINUE
 11   CONTINUE
C
C     ITERATIONS
C     ----------
 10   NBCO = NBCO + 1
      IF (NBCO.GT.NBCOMX) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = ' ERREUR : NOMBRE DE COULEURS INSUFFISANT'
         CALL LEREUR
         IERR = 1
         RETURN
      END IF
      CHGT = .FALSE.
      DO 2 NE = 1 , NBEL
         COLLOG(NE) = .FALSE.
 2    CONTINUE
      DO 3 NEA = 1 , NBEL
         NE = RENUEL(NEA)
         IF (.NOT.COLLOG(NE)) THEN
         IF (ELELOG(NE)) THEN
            ELELOG(NE) = .FALSE.
            NUCOL = COUELE(NE)
            L1 = LPVOIS(NE)+1
            L2 = LPVOIS(NE+1)
            DO 4 L = L1 , L2
               NUELVO = LISVOI(L)
C              LES VOISINS SONT-ILS DE MEME COULEUR ?
               IF (COUELE(NUELVO) .EQ. NUCOL ) THEN
                  COUELE(NUELVO) =  NBCO
                  COLLOG(NUELVO) = .TRUE.
                  CHGT = .TRUE.
               END IF
 4          CONTINUE
         END IF
         END IF
 3    CONTINUE
      IF (CHGT) THEN
        GO TO 10
      ELSE
        NBCO = NBCO - 1
        GO TO 20
      END IF
C
C     RANGEMENT PAR COULEUR
C     ---------------------
20    L = 0
      LPCOLO(1) = L
      DO 5 NC = 1 , NBCO
         DO 6 NE = 1 , NBEL
            IF (COUELE(NE) .EQ. NC) THEN
                L = L + 1
                ELECOU(L) = NE
            END IF
 6       CONTINUE
         LPCOLO(NC+1) = L
 5    CONTINUE
C
C     IMPRESSION DES RESULTATS
C     ------------------------
      WRITE(IMPRIM,1000)
      WRITE(IMPRIM,1100) NBCO
      IF (NOIMPR.GE.2) THEN
         DO 7 NC = 1 , NBCO
            NEL = LPCOLO(NC+1) - LPCOLO(NC)
            WRITE(IMPRIM,1200) NC,NEL
 7       CONTINUE
         DO 8 NC = 1 , NBCO
            L1 = LPCOLO(NC)+1
            L2 = LPCOLO(NC+1)
            WRITE(IMPRIM,1300) NC,(ELECOU(L),L=L1,L2)
 8       CONTINUE
         WRITE(IMPRIM,1400) (NE,COUELE(NE),NE=1,NBEL)
      ENDIF
C
      RETURN
1000  FORMAT(/,' RENUMEROTATION DES ELEMENTS PAR COULEURS',/)
1100  FORMAT(' NOMBRE TOTAL DE COULEURS',I8/)
1200  FORMAT(' COULEUR',I4,I6,' ELEMENTS')
1300  FORMAT(/,' COULEUR ',I6,/,1X,7(1H-),100(/,' ELEMENTS',10I6))
1400  FORMAT(' ELEMENT',I6,' COULEUR',I6)
      END
