      SUBROUTINE RESOUR( NYOBJT, NUOBJT, NBCOOR, XYZPI,
     %                   MNSOUR, SOURCE )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++012
C BUT :    REMPLIR LE TABLEAU SOUR DE LA SOURCE EXERCEE SUR L'OBJET
C -----    EN UN POINT D'INTEGRATION NUMERIQUE A L'INSTANT TEMPS
C
C ENTREES:
C --------
C NYOBJT : NUMERO DU TYPE D'OBJET (1:POINT, 2:LIGNE, ... )
C NUOBJT : NUMERO DU PLSV DANS SON LEXIQUE
C NBCOOR : NOMBRE DE COORDONNEES DES POINTS D'INTEGRATION 3 ou 6
C XYZPI  : LES 3 ou 6 COORDONNEES DU POINT D'INTEGRATION NUMERIQUE
C MNSOUR : ADRESSE MCN DU TABLEAU 'source'
C
C SORTIE :
C --------
C SOURCE : LA SOURCE EXERCEE
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++012
C AUTEUR : ALAIN PERRONNET   ANALYSE NUMERIQUE UPMC PARIS   JUILLET 1989
C MODIF  : ALAIN PERRONNET   TEXAS A & M UNIVERSITY         JUILLET 2005
C....6...............................................................012
      IMPLICIT INTEGER(W)
      include"./incl/donthe.inc"
      include"./incl/a___source.inc"
      include"./incl/ctemps.inc"
      include"./incl/cthet.inc"
      include"./incl/cnonlin.inc"
C
      include"./incl/pp.inc"
      COMMON           MCN(MOTMCN)
      REAL             RMCN(1)
      EQUIVALENCE     (MCN(1),RMCN(1))
C
      DOUBLE PRECISION XYZPI(NBCOOR), XYZ(10), SOURCE
C
C     LE TYPE DES DONNEES DE LA SOURCE
      LTSOUR = MCN( MNSOUR + WTSOUR )
C
C     REMPLISSAGE SELON LE TYPE
      IF( LTSOUR .EQ. 1 ) THEN
C
C        VALEUR CONSTANTE POUR CE PLSV
C        =============================
         SOURCE = RMCN( MNSOUR + WOURCE )
C
      ELSE
C
C        FONCTION DE L'UTILISATEUR
C        =========================
         XYZ(1) = TEMPS
         N = 1
         DO 5 I=1,NBCOOR
            N = N + 1
            XYZ(N) = XYZPI(I)
 5       CONTINUE
         XYZ(N+1) = NYOBJT
         XYZ(N+2) = NUOBJT
C        SOURCE DEPENDANTE DE LA TEMPERATURE (CAS NON LINEAIRE) OU NON
C        LA VALEUR DE LA TEMPERATURE EN CE POINT
C        TEMPERATURE CALCULEE ET STOCKEE DANS LE COMMON de cthet.inc
         N = N + 3
         XYZ(N) = TEMPEL
         CALL FONVAL( MCN(MNSOUR+WFSOUR), N, XYZ,
     %                NCODEV, SOURCE )
C
      ENDIF
      END
