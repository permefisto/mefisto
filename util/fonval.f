      SUBROUTINE FONVAL( NOFONC, NBPAFO, DPARAF,  NCODEV, DBLVAL )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     CALCULER LA VALEUR DE LA FONCTION DE NUMERO NOFONC DANS LE
C -----     LEXIQUE DES FONCTIONS ET DE PARAMETRES DPARAF
C           CF ~/td/g/grammaire_lu DE DEFINITION DU LANGAGE UTILISATEUR
C
C ENTREES :
C ---------
C NOFONC  : NUMERO>0 DE LA FONCTION DANS LE LEXIQUE DES FONCTIONS
C NBPAFO  : NOMBRE DE PARAMETRES DE LA FONCTION
C DPARAF  : VALEURS DES NBPAFO PARAMETRES POUR CALCULER LA FONCTION
C
C SORTIES :
C ---------
C NCODEV  : 0 DBLVAL N'EST PAS INITIALISEE EN SORTIE et DBLVAL=0D0
C           1 DBLVAL   EST     INITIALISEE EN SORTIE
C DBLVAL  : VALEUR REELLE DOUBLE PRECISION
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS    FEVRIER 1990
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/lu.inc"
      include"./incl/gsmenu.inc"
      include"./incl/a_fonction__arbre.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)
      DOUBLE PRECISION  DBLVAL, DPARAF(1:NBPAFO)
      CHARACTER*24      NMFONC

      IF( NOFONC .LE. 0 ) GOTO 9000

C     LES VALEURS DES PARAMETRES SONT EMPILES DANS DCTE
C     =================================================
C     L'INDICE DU PREMIER PARAMETRE DANS DCTE
      ICPAFO = NBCTE

      IF( NBCTE+NBPAFO .GT. MXDCTE ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'LU: dans FONVAL la PILE DCTE est SATUREE'
         ELSE
            KERR(1) = 'LU: in FONVAL the STACK DCTE is SATURATED'
         ENDIF
         CALL LEREUR
         print *,'Fonval: pile des constantes: nbcte=',nbcte,
     %           ' nbpaf0=',NBPAFO
         print *,dcte
         NCODEV = 0
         RETURN
      ENDIF

      DO 10 I=1,NBPAFO
C        LE PARAMETRE I DE LA FONCTION EST EMPILE DANS DCTE
         NBCTE = NBCTE + 1
         DCTE( NBCTE ) = DPARAF( I )
 10   CONTINUE

C     LE CALCUL DE LA FONCTION DE NUMERO NOFONC
C     =========================================
      CALL FONVA0( NOFONC, NBPAFO, ICPAFO+1,  NCODEV, DBLVAL )
      IF( NCODEV .LE. 0 ) GOTO 9000

C     LE NOMBRE DE CONSTANTES EST REMIS A SA VALEUR INITIALE
      NBCTE = ICPAFO
      GOTO 9999

C     ERREUR SUR LE NO DE LA FONCTION
 9000 IF( NOFONC .LE. 0 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            NMFONC = 'FONCTION USAGER NO<=0'
         ELSE
            NMFONC = 'USER FUNCTION of NO<=0'
         ENDIF
      ELSE
         CALL NMOBNU( 'FONCTION', NOFONC, NMFONC )
      ENDIF
      NBLGRC(NRERR) = 1
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = 'fonval: MODIFIEZ les DONNEES et/ou VOTRE FONCTION de
     % NOM ' // NMFONC
      ELSE
         KERR(1)='fonval: MODIFY the DATA and/or your FUNCTION of NAME '
     %              // NMFONC
      ENDIF
      CALL LEREUR
      DBLVAL = 0D0
      NCODEV = 0

 9999 RETURN
      END
