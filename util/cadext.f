      SUBROUTINE CADEXT( MNSOMM, COIN )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     MINIMUM ET MAXIMUM DES SOMMETS D'UN OBJET DEFINI
C ----      PAR SON TABLEAU D'ADRESSE MCN MNSOMM
C
C ENTREE  :
C ---------
C MNSOMM  : ADRESSE MCN DU TABLEAU 'XYZSOMMET' A TRAITER
C
C ENTREE ET SORTIE :
C ------------------
C COIN    : XYZUVW   MIN ET MAX DES NBCOOR DES SOMMETS
C           COIN(I,1) = MIN DE LA I-EME COORDONNEE
C           COIN(I,2) = MAX DE LA I-EME COORDONNEE
C                       I=1,...,NBCOOR
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS   SEPTEMBRE 1988
C MODIFS : ALAIN PERRONNET TEXAS A & M UNIVERSITY           JUILLET 2005
C23456---------------------------------------------------------------012
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL             RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      include"./incl/a___xyzsommet.inc"
      include"./incl/langue.inc"
      REAL             COIN(6,2)
C
C     LE NOMBRE DE XYZ
      NBSOM  = MCN( MNSOMM + WNBSOM )
C
C     LE NOMBRE DE COORDONNEES D'UN SOMMET
      NBCOOR = MCN( MNSOMM + WBCOOR )
C
C     LE MIN MAX DES COORDONNEES DES SOMMETS
      CALL MIMXPT( NBCOOR, NBSOM, RMCN(MNSOMM + WYZSOM), COIN )
C
cccC     AFFICHAGE DES COORDONNEES MIN ET MAX
ccc      WRITE(IMPRIM,*)
ccc      IF( LANGAG .EQ. 0 ) THEN
ccc         DO 10 I=1,NBCOOR
cccC           LA COORDONNEE I
ccc            WRITE(IMPRIM,*) 'PLSVO de COORDONNEE',I,
ccc     %                      ': MIN=', COIN(I,1),' MAX=',COIN(I,2)
ccc 10      CONTINUE
ccc      ELSE
ccc         DO 20 I=1,NBCOOR
cccC           LA COORDONNEE I
ccc            WRITE(IMPRIM,*) 'PLSVO of COORDINATE',I,
ccc     %                      ': MIN=', COIN(I,1),' MAX=',COIN(I,2)
ccc 20      CONTINUE
ccc      ENDIF
C
      RETURN
      END
