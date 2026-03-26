      INTEGER FUNCTION NUMTYP( KTYPE )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : RETOURNER LE NUMERO DU TYPE EN FONCTION DE SON NOM
C -----
C ENTREE :
C --------
C KTYPE  : NOM DE 9 CARACTERES
C
C SORTIE :
C --------
C NUMTYP : LOGIQUE  => 1 CARACTERE=> 2  ENTIER/2 => 3 ENTIER   => 4
C          REEL     => 5 REEL2    => 6  REEL4    => 7 COMPLEXE => 8
C          COMPLEXE2=> 9 MOTS     => 10 ^LEXIQUE =>11 XYZ      =>12
C          TYPEOBJET=>13
C          SI LE TYPE N'EST PAS RETROUVE ALORS DIAGNOSTIC ET ARRET
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : PERRONNET ALAIN ANALYSE NUMERIQUE PARIS  NOVEMBRE 1983
C.......................................................................
      include"./incl/msvaau.inc"
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      CHARACTER*9       K
      CHARACTER*(*)     KTYPE
C
C     PROTECTION DE KTYPE EN ENTREE ET COMPLETION OU TRONCATURE PAR
C     DES BLANCS
      K = KTYPE
C
C     COMPARAISON DE K AVEC LES TYPES
      DO 20 NUMTYP=1,NBTYPV
         IF( K .EQ. KTYPES(NUMTYP) ) RETURN
 20   CONTINUE
C
C     TYPE NON RETROUVE . ERREUR
      NBLGRC(NRERR) = 1
      KERR(1) ='TYPE INCONNU DE VARIABLE '// KTYPE
      CALL LEREUR
      WRITE(IMPRIM,10020) KTYPE,KTYPES
10020 FORMAT('ERREUR NUMTYP : TYPE ',A9,' INCONNU PARMI'/5(2X,A9))
      NUMTYP = 0
      END
