      SUBROUTINE NVARID( NOIDEN , NBVARI , NOTYPE )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LE NOMBRE NBVARI DE VARIABLES DEFINIES PAR
C -----    L'IDENTIFICATEUR NOIDEN
C
C ENTREE :
C --------
C NOIDEN : NUMERO DE L'IDENTIFICATEUR A AFFICHER
C
C SORTIES:
C --------
C NBVARI : NOMBRE DE VARIABLES DE CET IDENTIFICATEUR
C NOTYPE : NUMERO DE TYPE DE LA VARIABLE
C          LOGIQUE  => 1 CARACTERE=> 2  ENTIER/2 => 3 ENTIER   => 4
C          REEL     => 5 REEL2    => 6  REEL4    => 7 COMPLEXE => 8
C          COMPLEXE2=> 9 MOTS     => 10 ^LEXIQUE =>11 XYZ      =>12
C          TMS      => 21
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS SEPTEMBRE 1988
C23456---------------------------------------------------------------012
C.......................................................................
      include"./incl/td.inc"
C.......................................................................
C     RECHERCHE DU NOMBRE DE VARIABLES SI C'EST UN TABLEAU
      NBVARI = 1
C
C     POSITION DANS LE TAS
      N = IDENT(2,NOIDEN)
      IF( N .GT. 0 ) THEN
C        C'EST UN TABLEAU
C        LE NOMBRE DE SES INDICES
         M = LETAS( N )
         DO 10 I=1,M
            NBVARI = NBVARI * ( LETAS(N+2) - LETAS(N+1) + 1 )
            N      = N + 2
 10      CONTINUE
      ENDIF
C
C     LE TYPE DE LA VARIABLE A AFFICHER
      NOTYPE = IDENT(1,NOIDEN)
      END
