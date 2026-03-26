      SUBROUTINE MSINIT
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : INITIALISER QUELQUES VARIABLES DEPENDANTES DE L ORDINATEUR
C -----
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS DECEMBRE 1983
C.......................................................................
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      include"./incl/msvaau.inc"
C
C     VARIABLES NON DEPENDANTES DE L'ORDINATEUR.A NE PAS MODIFIER
C     ===========================================================
C
C     LE NOMBRE DE TYPES DE VARIABLES
C          LOGIQUE  => 1 CARACTERE=> 2  ENTIER/2 => 3 ENTIER   => 4
C          REEL     => 5 REEL2    => 6  REEL4    => 7 COMPLEXE => 8
C          COMPLEXE2=> 9 MOTS     => 10 ^LEXIQUE =>11 XYZ      =>12
C          TYPEOBJET=>13
C
C     LES DIFFERENTS TYPES EN CARACTERES
      KTYPES( 1 ) = 'LOGIQUE  '
      KTYPES( 2 ) = 'CARACTERE'
      KTYPES( 3 ) = 'ENTIER/2 '
      KTYPES( 4 ) = 'ENTIER   '
      KTYPES( 5 ) = 'REEL     '
      KTYPES( 6 ) = 'REEL2    '
      KTYPES( 7 ) = 'REEL4    '
      KTYPES( 8 ) = 'COMPLEXE '
      KTYPES( 9 ) = 'COMPLEXE2'
      KTYPES(10 ) = 'MOTS     '
      KTYPES(11 ) = '^LEXIQUE '
      KTYPES(12 ) = 'XYZ      '
      KTYPES(13 ) = 'TYPEOBJET'
C
C     VARIABLES DEPENDANTES DE L'ORDINATEUR.A MODIFIER
C     ================================================
C
C     LES UNITES DE LECTURE ET D'IMPRESSION.
C     LES COMMENTAIRES SONT DANS LE CADRE D'UNE EXECUTION PREALABLE
C     DU SP INITIA QUI INITIALISE LUI AUSSI LECTEU ET IMPRIM
C     CES 2 INITIALISATIONS SONT A RESTAURER SINON
C
C     LECTEU = IINFO( 'LECTEUR INITIAL' )
C     IMPRIM = IINFO( 'IMPRIMANTE INITIALE' )
C
C     LE TABLEAU VAD1MO : NOMBRE DE VARIABLES DANS UN MOT MEMOIRE
C                         POUR CHAQUE TYPE DE VARIABLE
      VAD1MO(1) = 1.
      VAD1MO(2) = 4.
      VAD1MO(3) = 2.
      VAD1MO(4) = 1.
      VAD1MO(5) = 1.
      VAD1MO(6) = 0.5
      VAD1MO(7) = 0.25
      VAD1MO(8) = 0.5
      VAD1MO(9) = 0.25
      VAD1MO(10)= 1.
      VAD1MO(11)= 1.
      VAD1MO(12)= 1./3.
      VAD1MO(13)= 0.5
C
C     LE TABLEAU MOTSVA : NOMBRE DE MOTS D UNE VARIABLE SELON LE TYPE
      DO 10 I=1,5
         MOTSVA( I ) = 1
10    CONTINUE
      MOTSVA( 6 ) = 2
      MOTSVA( 7 ) = 4
      MOTSVA( 8 ) = 2
      MOTSVA( 9 ) = 4
      MOTSVA(10 ) = 1
      MOTSVA(11 ) = 1
      MOTSVA(12 ) = 3
      MOTSVA(13 ) = 2
C
C     LE TABLEAU MOTYPV : LE NOMBRE DE VARIABLES DANS UN MOT POUR I=1,3
C                         LE NOMBRE DE MOTS PAR VARIABLES    POUR I=4,9
      MOTYPV( 1 ) = 1
      MOTYPV( 2 ) = 4
      MOTYPV( 3 ) = 2
C
      MOTYPV( 4 ) = 1
      MOTYPV( 5 ) = 1
      MOTYPV( 6 ) = 2
      MOTYPV( 7 ) = 4
      MOTYPV( 8 ) = 2
      MOTYPV( 9 ) = 4
      MOTYPV( 10) = 1
      MOTYPV( 11) = 1
      MOTYPV( 12) = 3
      MOTYPV( 13) = 2
      END
