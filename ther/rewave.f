      SUBROUTINE REWAVE( NYOBJT, NUOBJT, NBCOOR, XYZPI,
     %                   MNSOUR, FOMEGA )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++012
C BUT :    REMPLIR LE TABLEAU FOMEGA DE LA SOURCE EXERCEE SUR L'OBJET
C -----    EN UN POINT D'INTEGRATION NUMERIQUE A L'INSTANT TEMPS
C          DANS LE CAS D'UNE ONDE REGIE PAR LES NLSE
C
C ENTREES:
C --------
C NYOBJT : NUMERO DU TYPE D'OBJET (1:POINT, 2:LIGNE, ... )
C NUOBJT : NUMERO DU PLSV DANS SON LEXIQUE
C NBCOOR : NOMBRE DE COORDONNEES DES POINTS D'INTEGRATION 3 ou 6
C XYZPI  : LES 3 COORDONNEES DU POINT D'INTEGRATION NUMERIQUE
C MNSOUR : ADRESSE MCN DU TABLEAU 'source'
C
C SORTIE :
C --------
C FOMEGA : FOMEGA(2) PARTIE REELLE ET IMAGINAIRE DE LA SOURCE
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++012
C AUTEUR : ALAIN PERRONNET Texas A & M University at QATAR  Fevrier 2011
C23456---------------------------------------------------------------012
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
      DOUBLE PRECISION XYZPI(NBCOOR), PARAM(9), FOMEGA(2)
C
C     LE TYPE DES DONNEES DE LA SOURCE
      LTSOUR = MCN( MNSOUR + WTSOUR )
C
C     REMPLISSAGE SELON LE TYPE
      IF( LTSOUR .EQ. 1 ) THEN
C
C        VALEUR CONSTANTE POUR CE PLSV
C        =============================
         FOMEGA(1) = RMCN( MNSOUR + WOURCE )
         FOMEGA(2) = FOMEGA(1)
C
      ELSE
C
C        FONCTION DE L'UTILISATEUR
C        =========================
         PARAM(1) = TEMPS
         N = 1
         DO 5 I=1,NBCOOR
            N = N + 1
            PARAM(N) = XYZPI(I)
 5       CONTINUE
         PARAM(N+1) = NYOBJT
         PARAM(N+2) = NUOBJT
C
C        SOURCE DEPENDANTE DE LA TEMPERATURE (CAS NON LINEAIRE) OU NON
C        LA VALEUR DE LA TEMPERATURE EN CE POINT
C        TEMPERATURE CALCULEE ET STOCKEE DANS LE COMMON de cthet.inc
         PARAM(N+3) = TEMPEL
         PARAM(N+4) = ONDEPI
C
         N = N + 5
C
C        LES 2 PARTIES REELLE et IMAGINAIRE de FOmega AUX 3 SOMMETS
         DO M=1,2
C
C           NUMERO DE LA PARTIE: 1 POUR REELLE, 2 POUR IMAGINAIRE
            PARAM(N) = M
C           LA VALEUR DES SOURCES SURFACIQUES EN CE SOMMET K
            CALL FONVAL( MCN(MNSOUR+WFSOUR), N, PARAM,
     %                   NCODEV, FOMEGA(M) )
         ENDDO
C
      ENDIF
C
      RETURN
      END
