      INTEGER FUNCTION MOTVAR( NUMTYP )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : RETOURNER LE NOMBRE DE MOTS D UNE VARIABLE DE NO DE TYPE NUMTYP
C -----
C
C ENTREE :
C --------
C NUMTYP : NUMERO COMPRIS ENTRE 1 ET 13    OU  21 ( TMS )
C
C SORTIE :
C --------
C MOTVAR : NOMBRE ENTIER DE MOTS OCCUPES PAR UNE VARIABLE DE TYPE NUMTYP
C          LOGIQUE  => 1 CARACTERE=> 1 ENTIER/2 => 1  ENTIER   => 1
C          REEL     => 1 REEL2    => 2 REEL4    => 4  COMPLEXE => 2
C          COMPLEXE2=> 4 MOTS     => 1 ^LEXIQUE => 1  XYZ      => 3
C          TMS      => 1 XYZ      => 3 TYPEOBJET=> 2
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : PERRONNET ALAIN ANALYSE NUMERIQUE PARIS  NOVEMBRE 1983
C.......................................................................
      include"./incl/msvaau.inc"
C
C     L AFFECTATION DE LA VALEUR
      IF( NUMTYP .NE. 21 ) THEN
C        VARIABLE CLASSIQUE
         MOTVAR = MOTSVA( NUMTYP )
      ELSE
C        TABLEAU TMS
         MOTVAR = 1
      ENDIF
      END
