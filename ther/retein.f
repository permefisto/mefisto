      SUBROUTINE RETEIN( NYOBJT, NUOBJT, NBCOOR, XYZNO, MNTEM0, TEMPE0 )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++012
C BUT :    REMPLIR LA VALEUR DE LA TEMPERATURE A L'INSTANT TEMPS
C -----    (EN FAIT C'EST LE PLUS SOUVENT A L'INSTANT INITIAL TEMPS=0!)
C          DANS LE CAS OU ELLE EST CONSTANTE EN TOUS LES NOEUDS
C          OU DEFINIE PAR UNE FONCTION UTILISATEUR EN UN NOEUD
C          DE COORDONNEES XNO,YNO,ZNO A L'INSTANT TEMPS ET CAS 1

C ENTREES :
C ---------
C NYOBJT : NUMERO DU TYPE D'OBJET (1:POINT, 2:LIGNE, ..., 5:OBJET )
C NUOBJT : NUMERO DE L'OBJET DANS SON LEXIQUE
C XNO,YNO,ZNO : LES 3 COORDONNEES DU NOEUD
C MNTEM0 : ADRESSE MCN DU TABLEAU 'TEMPERINIT'
C
C SORTIE :
C --------
C TEMPE0 : TEMPERATURE DU NOEUD A L'INSTANT TEMPS
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++012
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     JANVIER 1997
C....6...............................................................012
      IMPLICIT INTEGER(W)
      include"./incl/donthe.inc"
      include"./incl/a___temperinit.inc"
      include"./incl/ctemps.inc"
      include"./incl/cthet.inc"
      include"./incl/cnonlin.inc"

      DOUBLE PRECISION  XYZNO(NBCOOR), TEMPE0
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      DOUBLE PRECISION  XYZ(9)
      EQUIVALENCE      (MCN(1),RMCN(1))

C     LE TYPE DES DONNEES DE LA TEMPERATURE AU NOEUD
      LTTEM0 = MCN( MNTEM0 + WTTEM0 )

C     REMPLISSAGE SELON LE TYPE
      IF( LTTEM0 .EQ. 1 ) THEN

C        MATERIAU HOMOGENE ISOTROPE
C        ==========================
         TEMPE0 = RMCN( MNTEM0 + WATEM0 )

      ELSE IF( LTTEM0 .EQ. -1 ) THEN

C        FONCTION UTILISATEUR
C        ====================
         XYZ(1) = TEMPS
         N = 1
         DO 5 I=1,NBCOOR
            N = N + 1
            XYZ( N ) = XYZNO(I)
 5       CONTINUE
         XYZ(N+1) = NYOBJT
         XYZ(N+2) = NUOBJT
         N = N + 2
         CALL FONVAL( MCN(MNTEM0+WFTEM0), N, XYZ,
     %                NCODEV, TEMPE0 )
      ENDIF

      RETURN
      END
