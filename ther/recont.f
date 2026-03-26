      SUBROUTINE RECONT( NYOBJT, NUOBJT, NBCOOR, XYZP,
     %                   MNCONT, CONTAC )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++012
C BUT :    REMPLIR LE TABLEAU CONTAC DE LA TEMPERATURE IMPOSEE
C -----    EN UN NOEUD DE L'OBJET A L'INSTANT TEMPS

C ENTREES :
C ---------
C NYOBJT : NUMERO DU TYPE OBJET (1:POINT, 2:LIGNE, ... )
C NUOBJT : NUMERO DU PLSV DANS SON LEXIQUE
C NBCOOR : NOMBRE DE COORDONNEES DES POINTS D'INTEGRATION 3 ou 6
C XYZP   : LES 3 ou 6 COORDONNEES DU POINT DE CALCUL DE LA TEMPERATURE
C MNCONT : ADRESSE MCN DU TABLEAU 'CONTACT'

C SORTIE :
C --------
C CONTAC : LA TEMPERATURE IMPOSEE
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++012
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    DECEMBRE 1990
C MODIF  : ALAIN PERRONNET   TEXAS A & M UNIVERSITY         JUILLET 2005
C....6...............................................................012
      IMPLICIT INTEGER(W)
      include"./incl/donthe.inc"
      include"./incl/a___contact.inc"
      include"./incl/ctemps.inc"
      include"./incl/cthet.inc"

      DOUBLE PRECISION XYZP(NBCOOR), XYZ(10), CONTAC

      include"./incl/pp.inc"
      COMMON           MCN(MOTMCN)
      REAL             RMCN(1)
      EQUIVALENCE     (MCN(1),RMCN(1))

C     LE TYPE DES DONNEES DU CONTACT OU TEMPERATURE IMPOSEE
      LTCONT = MCN( MNCONT + WTCONT )

C     REMPLISSAGE SELON LE TYPE
      IF( LTCONT .EQ. 1 ) THEN
C
C        VALEURS CONSTANTES POUR CET OBJET
C        =================================
         CONTAC = RMCN( MNCONT + WONTAC )

      ELSE

C        FONCTION UTILISATEUR(t,x,y,z,ntyplsv,nuplsv,tempel)
C        ====================
         XYZ(1) = TEMPS
         N = 1
         DO I=1,NBCOOR
            N = N + 1
            XYZ(N) = XYZP(I)
         ENDDO
         XYZ(N+1) = NYOBJT
         XYZ(N+2) = NUOBJT
C        PROBLEME DEPENDANT DE LA TEMPERATURE (CAS NON LINEAIRE) OU NON
C        LA VALEUR DEJA CALCULEE DE LA TEMPERATURE EN CE POINT
C        TEMPERATURE CALCULEE ET STOCKEE DANS LE COMMON de cthet.inc
         N = N + 3
         XYZ(N) = TEMPEL
         CALL FONVAL( MCN(MNCONT+WFCONT), N, XYZ,
     %                NCODEV, CONTAC )
      ENDIF

      RETURN
      END
