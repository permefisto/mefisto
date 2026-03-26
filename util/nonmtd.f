      INTEGER FUNCTION NONMTD( NMTD )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     RETOUVER LE NUMERO D'UN TD A PARTIR DE SON NOM
C -----
C
C ENTREES :
C ---------
C NMTD    : NOM DU TABLEAU DESCRIPTEUR A CHERCHER DANS DICOTD
C
C SORTIES :
C ---------
C NONMTD  : NUMERO DU TABLEAU TD DANS DICOTD
C           0 SI NON RETROUVE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS SEPTEMBRE 1988
C23456---------------------------------------------------------------012
      include"./incl/td.inc"
      include"./incl/gsmenu.inc"
C.......................................................................
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      CHARACTER*(*)     NMTD
      CHARACTER*(NCKTD) NOM
      CHARACTER*(NCKTD) NOMMIN
C
C     RECHERCHE DU DERNIER CARACTERE NON BLANC
      NOM = NMTD
      DO 10 I=NCKTD,1,-1
         IF( NOM(I:I) .NE. ' ' ) GOTO 20
 10   CONTINUE
      NBLGRC(NRERR) = 1
      KERR(1) = 'NOM BLANC DANS NONMTD'
      CALL LEREUR
      GOTO 9000
C
C     SUPPRESSION DES SUFFIXES DE SEPARATEUR KSUFIX
 20   J = 1
C
 25   L = INDEX( NOM(J:I) , KSUFIX )
      IF( L .GT. 0 ) THEN
C        RECHERCHE DE >
         L  = J + L - 1
         LL = INDEX( NOM(L+1:I) , '>' )
         IF( LL .GT. 0 ) THEN
            NOM = NOM(1:L-1) // NOM(L+LL:I)
C           LE 1-ER CARACTERE A EXAMINER
            J        = L
C           LE DERNIER CARACTERE A EXAMINER
            I        = I - LL
         ENDIF
         GOTO 25
      ENDIF
C
C     IDENTIFICATION CARACTERE PAR CARACTERE EN MINUSCULES
      NOMMIN = NOM
      CALL MINUSC( NOMMIN )
      DO 30 NONMTD=1,NBTD
         IF( DICOTD(NONMTD)(1:I) .EQ. NOMMIN(1:I) ) RETURN
 30   CONTINUE
      NBLGRC(NRERR) = 1
      KERR(1) =  NMTD//' NON RETROUVE DANS DICOTD'
      CALL LEREUR
C
9000  NONMTD = 0
      RETURN
      END
