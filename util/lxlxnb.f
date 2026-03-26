       SUBROUTINE LXLXNB( MNLX , NBLXLX )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    COMPTER LE NOMBRE ACTUEL DE LEXIQUES DANS LE LEXIQUE
C -----    OUVERT DANS LE SUPER-TABLEAU MCN A L'ADRESSE MNLX
C ENTREE :
C --------
C MNLX   : ADRESSE MCN DU TABLEAU TAMS DU LEXIQUE PERE
C
C SORTIE :
C --------
C NBLXLX : NOMBRE ACTUEL DE LEXIQUES DANS LE LEXIQUE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS   MARS 1987
C.......................................................................
      include"./incl/pp.inc"
      COMMON  MCN(MOTMCN)
C
C     LE DEBUT DU CHAINAGE DES LEXIQUES
      MN   = MNLX + 5
C
C     LE NOMBRE DE LEXIQUES
      NBLXLX = 0
      IF( MCN( MN ) .LE. 0 ) RETURN
C
C     LE NOMBRE DE MOTS PAR NOM ET ATTRIBUTS
      M1LX = MCN( MNLX )
C
C     LE NOMBRE D'ENTIERS POUR UN NOM
      NBENNM = MCN( MNLX + 2 )
C
C     'LEXI' TRADUIT EN ENTIER
      ILEXI = ICHARX( 'LEXI' )
C
C     POSITION DANS LE LEXIQUE AU DEBUT DU NOM
 10   MN = MNLX + M1LX * MCN( MN )
      IF( MCN(MN+NBENNM+1) .EQ. ILEXI ) THEN
C        C'EST UN LEXIQUE
         NBLXLX = NBLXLX + 1
      ENDIF
C
C     LE NOM SUIVANT
      MN = MCN( MN + NBENNM )
      IF( MN .GT. 0 ) GOTO 10
      END
