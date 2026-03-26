      SUBROUTINE REDEIN( NYOBJT, NUOBJT, NBDEP0, XNO,YNO,ZNO, MNDEP0,
     %                   DEPLA0 )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++012
C BUT :    REMPLIR LA VALEUR DU DEPLACEMENT A L'INSTANT TEMPS
C -----    (EN FAIT C'EST LE PLUS SOUVENT A L'INSTANT INITIAL TEMPS=0!)
C          DANS LE CAS OU IL EST CONSTANT EN TOUS LES NOEUDS
C          OU DEFINI PAR UNE FONCTION UTILISATEUR EN UN NOEUD
C          DE COORDONNEES XNO,YNO,ZNO A L'INSTANT TEMPS ET POUR LE CAS 1
C
C ENTREES :
C ---------
C NYOBJT : NUMERO DU TYPE D'OBJET (1:POINT, 2:LIGNE, ..., 5:OBJET )
C NUOBJT : NUMERO DE L'OBJET DANS SON LEXIQUE
C NBDEP0 : 1 OU 2 OU 3 NOMBRE DE COMPOSANTES DU DEPLACEMENT
C XNO,YNO,ZNO : LES 3 COORDONNEES DU NOEUD
C MNDEP0 : ADRESSE MCN DU TABLEAU 'DEPLACTINIT'
C
C SORTIE :
C --------
C DEPLA0 : TABLEAU (NBDEP0) DES NBDEP0 COMPOSANTES DU DEPLACEMENT
C          AU NOEUD XNO, YNO, ZNO A L'INSTANT TEMPS
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++012
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     JANVIER 1998
C....6...............................................................012
      IMPLICIT INTEGER(W)
      include"./incl/donela.inc"
      include"./incl/a___deplactinit.inc"
      include"./incl/ctemps.inc"
      include"./incl/gsmenu.inc"
C
      DOUBLE PRECISION  XNO, YNO, ZNO, DEPLA0(1:NBDEP0)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      DOUBLE PRECISION  XYZ(7)
      EQUIVALENCE      (MCN(1),RMCN(1))
C
C     LE TYPE DES DONNEES DU DEPLACEMENT AU NOEUD
      LTDEP0 = MCN( MNDEP0 + WTDEP0 )
C
C     REMPLISSAGE SELON LE TYPE
      IF( LTDEP0 .EQ. 1 ) THEN
C
C        DEPLACEMENT CONSTANT
C        ====================
         IF( MCN( MNDEP0 + WBDEP0 ) .NE. NBDEP0 ) THEN
            WRITE(KERR(4)(1:1),'(I1)') NBDEP0
            NBLGRC(NRERR) = 3
            KERR(1) ='ERREUR: NOMBRE COMPOSANTES DEPLACEMENT INITIAL'
            KERR(2) ='NON EGAL A ' // KERR(4)(1:1)
            KERR(3) ='DEPLACEMENT INITIAL NUL IMPOSE'
            CALL LEREUR
            DO 10 I=1,NBDEP0
               DEPLA0(I) = 0D0
 10         CONTINUE
         ELSE
            DO 20 I=1,NBDEP0
               DEPLA0(I) = RMCN( MNDEP0 + WADEP0 - 1 + I )
 20         CONTINUE
         ENDIF
C
      ELSE IF( LTDEP0 .EQ. -1 ) THEN
C
C        FONCTION UTILISATEUR
C        ====================
         XYZ(1) = TEMPS
         XYZ(2) = XNO
         XYZ(3) = YNO
         XYZ(4) = ZNO
         XYZ(5) = NYOBJT
         XYZ(6) = NUOBJT
         DO 30 I=1,NBDEP0
C           LE NUMERO DE LA COMPOSANTE DU DEPLACEMENT
            XYZ(7) = I
            CALL FONVAL( MCN(MNDEP0+WFDEP0), 7, XYZ,
     %                   NCODEV, DEPLA0(I) )
 30      CONTINUE
      ENDIF
      END
