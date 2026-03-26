      SUBROUTINE TYARCP( NOTYAR, NBP, COPOAR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    DONNER LES NBP COORDONNEES CURVILIGNES SUR (0.,1.)
C -----    DES NBP POINTS NON SOMMETS DE L ARETE DE TYPE NOTYAR
C
C ENTREES:
C --------
C NOTYAR : NO DU TYPE DE L ARETE DONT LES COORDONNEES DES POINTS
C          NON SOMMETS SONT DEMANDES
C NBP    : NOMBRE DE POINTS NON SOMMETS DE CE TYPE D ARETE
C
C SORTIE :
C --------
C COPOAR : COORDONNEE DES NBP POINTS NON SOMMETS DE L ARETE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE UPMC PARIS   DECEMBRE 1980
C ......................................................................
      include"./incl/gsmenu.inc"
      REAL              COPOAR(NBP)
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
C
C     ******************************************************************
      GOTO ( 1000 , 20 , 1000 , 1000 , 1000 ) , NOTYAR
C     ******************************************************************
C
C     ==================================================================
C     NOTYAR = 1  ARETE DE TYPE P1      S1 *.......................* S2
C     ==================================================================
C
C     ==================================================================
C     NOTYAR = 2  ARETE DE TYPE P2      S1 *...........*...........* S2
C                                                    MILIEU
C     ==================================================================
 20   N = 1
      IF(NBP .NE. N) GOTO 9999
      COPOAR(1) = 0.5
      RETURN
C
C     ERREUR: NON CONCORDANCE ENTRE LE NOMBRE DE POINTS ATTENDU
C             ET CELUI DONNE
 9999 NBLGRC(NRERR) = 2
      WRITE(KERR(MXLGER)(1:10),'(I10)') NBP
      KERR(1) = 'TYARCP:NOMBRE DE POINTS EN ENTREE='//KERR(MXLGER)(1:10)
      WRITE(KERR(MXLGER)(1:10),'(I10)') N
      KERR(2) = ' AU LIEU DE ' // KERR(MXLGER)(1:10)
      CALL LEREUR
      COPOAR(1) = 0.
C
 1000 RETURN
      END
