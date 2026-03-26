      SUBROUTINE PRPRMC ( MNTOPO, MNNPEF, MNXYZN, NDLNOE, NCODSA,
     %                    LPDIAG, IERR)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    GENERER LE PROFIL DE LA MATRICE DU MAILLAGE POUR NDLNOE
C -----    DEGRES DE LIBERTE ( CONSTANT EN TOUT NOEUD ) PAR NOEUD
C
C ENTREES :
C ---------
C MNTOPO : ADRESSE  MCN DU TABLEAU   'TOPOLOGIE' DE L'OBJET
C MNNPEF : ADRESSES MCN DES TABLEAUX 'NPEF"TYPE_EF' DE L'OBJET
C MNXYZN : ADRESSE  MCN DU TABLEAU   'XYZNOEUD'  DE L'OBJET
C NDLNOE : NOMBRE CONSTANT DE DEGRES DE LIBERTE PAR NOEUD
C NCODSA : CODE DE STOCKAGE DE LA MATRICE PROFIL
C          -1 : MATRICE NON SYMETRIQUE
C           1 : MATRICE SYMETRIQUE
C
C SORTIES :
C ---------
C LPDIAG : POINTEUR SUR LES COEFFICIENTS DIAGONAUX DE LA MATRICE PROFIL
C IERR   : 0 SI PAS D'ERREUR
C          1 SI UN TABLEAU 'NPEF"TYPE_EF' EST INCONNU
C          2 SI LE POINTEUR DEVIENT NEGATIF PAR DEBORDEMENT DE 2**31
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS         MAI 1989
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/a___xyznoeud.inc"
      include"./incl/a___npef.inc"
      include"./incl/a_objet__topologie.inc"
      include"./incl/pp.inc"
      COMMON          MCN(MOTMCN)
      COMMON /UNITES/ LECTEU,IMPRIM,NUNITE(30)
C
      INTEGER         MNNPEF(1:*), LPDIAG(0:*)
C
      IERR = 0
      IF( NCODSA .EQ. 0 ) RETURN
C
C     LE NOMBRE TOTAL DE NOEUDS ET DE DEGRES DE LIBERTE DE L'OBJET
      NBNOE  = MCN( MNXYZN + WNBNOE )
      NBLIPR = NBNOE * NDLNOE
C
C     MISE A NBNOE DU TABLEAU DU PLUS PETIT VOISIN DE CHAQUE NOEUD
C     ------------------------------------------------------------
      CALL TNMCDC( 'ENTIER', NBNOE, MNMIVO )
      MNMIV1 = MNMIVO - 1
      DO 10 I=MNMIVO,MNMIV1+NBNOE
         MCN( I ) = NBNOE
 10   CONTINUE
C
C     LA BOUCLE SUR LES TYPES D'EF DU MAILLAGE
C     ========================================
      NBTYEL = MCN( MNTOPO + WBTYEL )
      DO 40 NUTYEL=1,NBTYEL
C
C        LE TABLEAU DES ELEMENTS A POUR ADRESSE MNELE
         MNELE = MNNPEF( NUTYEL )
         IF( MNELE .LE. 0 ) THEN
            WRITE(KERR(4)(1:3),'(I3)') NUTYEL
            NBLGRC(NRERR) = 1
            IF( LANGAG .EQ. 0 ) THEN
              KERR(1)='PRPRMC: TYPE ELEMENT FINI INCONNU '//KERR(4)(1:3)
            ELSE
            KERR(1)='PRPRMC: UNKNOWN FINITE ELEMENT TYPE '//KERR(4)(1:3)
            ENDIF
            IERR = 1
            GOTO 40
         ENDIF
C
C        LA BOUCLE SUR LES EF DE CE TYPE
         NBELEM = MCN( MNELE + WBELEM )
         NBNDEL = MCN( MNELE + WBNDEL )
         MNNDEL = MNELE + WUNDEL - 1
C
         DO 30 NUELEM=1,NBELEM
            MNNDE = MNNDEL + NUELEM
            NMIN  = NBLIPR
C
C           NMIN = NO MINIMUM DES NO DE NOEUDS DE L ELEMENT FINI
C           ----------------------------------------------------
            MN = MNNDE
            DO 15 I=1,NBNDEL
               NMIN = MIN( NMIN , MCN(MN) )
               MN   = MN + NBELEM
 15         CONTINUE
C
C           MIVO(J) = NO DU PLUS PETIT VOISIN DU NOEUD J
C           --------------------------------------------
            MN = MNNDE
            DO 20 I=1,NBNDEL
C              LE NUMERO DU NOEUD I DE L'ELEMENT NUELEM
               N  = MCN(MN)
               MCN(MNMIV1+N) = MIN( MCN(MNMIV1+N) , NMIN )
               MN = MN + NBELEM
 20         CONTINUE
 30      CONTINUE
 40   CONTINUE
      IF( IERR .NE. 0 ) GOTO 9000
C
C     FORMATION DU TABLEAU LPDIAG DES POINTEURS SUR LA DIAGONALE
C     ==========================================================
      LPDIAG( 0 ) = 0
      IDL = 0
      IF( NCODSA .GT. 0 ) THEN
C
C        MATRICE SYMETRIQUE
C        ------------------
         DO 60 J=1,NBNOE
C           NDLNOE DEGRES DE LIBERTE PAR NOEUD
            NBC = ( J - MCN(MNMIV1+J) ) * NDLNOE
            DO 50 I=1,NDLNOE
               NBC = NBC + 1
               N   = IDL + 1
               LPDIAG( N ) = LPDIAG( IDL ) + NBC
C              DEBORDEMENT DES ENTIERS 2**31?
               IF( LPDIAG( N ) .LE. 0 ) GOTO 9000
               IDL = N
 50         CONTINUE
 60      CONTINUE
C
      ELSE
C
C        MATRICE NON SYMETRIQUE
C        ----------------------
         DO 80 J=1,NBNOE
C           NDLNOE DEGRES DE LIBERTE PAR NOEUD
            NBC = ( J - MCN(MNMIV1+J) ) * NDLNOE
            NBC = NBC + NBC - 1
            DO 70 I=1,NDLNOE
               NBC = NBC + 2
               N   = IDL + 1
               LPDIAG( N ) = LPDIAG( IDL ) + NBC
C              DEBORDEMENT DES ENTIERS 2**31?
               IF( LPDIAG( N ) .LE. LPDIAG( IDL ) ) GOTO 9000
               IDL = N
 70         CONTINUE
 80      CONTINUE
      ENDIF
      GOTO 9900
C
C     ERREUR DEBORDEMENT DU MAXIMUM DU STOCKAGE DES ENTIERS
 9000 NBLGRC(NRERR) = 3
      WRITE(KERR(5)(1:12),'(I12)') J
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1)='PRPRMC: LE POINTEUR DEBORDE 2**31 AU NOEUD'
     %            //KERR(5)(1:12)
         KERR(2)='PRPRMC: REDUIRE LE NOMBRE DE NOEUDS'
      ELSE
         KERR(1)='PRPRMC: The POINTER OVERFLOWS 2**31 at NODE'
     %            //KERR(5)(1:12)
         KERR(2)='PRPRMC: REDUCE the NODE NUMBER'
      ENDIF
      CALL LEREUR
      IERR = 2
C
C     DESTRUCTION DU TABLEAU AUXILIAIRE
 9900 CALL TNMCDS( 'ENTIER', NBNOE, MNMIVO )
      RETURN
      END
