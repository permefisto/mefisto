      PROGRAM  PPPARA
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CREER LES PARAMETER DES TABLEAUX DESCRIPTEURS DES TMS
C -----
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     OCTOBRE 1988
C.......................................................................
      COMMON /UNITES/ LECTEU,IMPRIM,NUNITE(30)
      include"./incl/td.inc"
      include"./incl/gsmenu.inc"
C
C     LECTEUR
      LECTEU = IINFO('LECTEUR INITIAL')

C     L'IMPRIMANTE
      IMPRIM = IINFO('IMPRIMANTE INITIALE')
      NBLGRC(NRERR) = 1
      KERR(1)='GENERATION DES PARAMETRES DES TABLEAUX DESCRIPTEURS'
      CALL LERESU
C
C     LECTURE DU DICTIONNAIRE DES TABLEAUX DESCRIPTEURS
      LHKTD = 0
      CALL LEDICO
C
C     INITIALISATION DE VARIABLES DE LA MS
      CALL MSINIT
C
C     BOUCLE SUR TOUS LES TABLEAUX DESCRIPTEURS
C     =========================================
      DO I=1,NBTD
         print *, 'pppara: dicotd(',i,')=', DICOTD(I)
         CALL PATSTD( DICOTD(I) )
      ENDDO
C
C     FIN
      NBLGRC(NRERR) = 1
      KERR(1) = 'pppara: FIN GENERATION DES PARAMETRES'
      CALL LERESU
C
      STOP
      END
