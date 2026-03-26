      SUBROUTINE AFCARX( M )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : AFFICHER LES 4 CARACTERES CODES DANS L'ENTIER M
C -----
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS   SEPTEMBRE 1988
C23456---------------------------------------------------------------012
      include"./incl/td.inc"
C.......................................................................
      include"./incl/nbcamo.inc"
C     NOMBRE DE CARACTERES POUR REPRESENTER LE PLUS GRAND CARX
C     VOIR LES FORMAT NON PARAMETRES POUR ACCELERER LE CALCUL
      CHARACTER*(NBCAMO) KCARX
      CHARACTER*(NBCAMO) CHARX
      COMMON / UNITES /  LECTEU,IMPRIM,NUNITE(30)
C
C     LA LIGNE EST BUFFERISEE A L'AIDE DU POINTEUR LCLIGN DERNIER
C     CARACTERE ENTRE DANS LA LIGNE
C
      IF( IMPRES*NOMUET .LE. 0 ) RETURN
C
      KCARX = CHARX( M )
C
C     RECHERCHE DU 1-ER NON BLANC
      N1 = 0
 10   N1 = N1 + 1
      IF( KCARX(N1:N1) .EQ. ' ' ) THEN
         IF( N1 .GE. NBCAMO ) THEN
C           TOUS LES CARACTERES SONT BLANCS
            RETURN
         ENDIF
         GOTO 10
      ENDIF
C
C     RECHERCHE DU DERNIER NON BLANC
      DO 20 N2=NBCAMO,N1,-1
         IF( KCARX(N2:N2) .NE. ' ' ) GOTO 30
 20   CONTINUE
C
C     LE NOMBRE DE CARACTERES A AFFICHER
 30   NBC = N2 + 1 - N1
      IF( NCLIGN-LCLIGN .LT. NBC ) THEN
         CALL AFLIGN
      ENDIF
C
      KLIGNE(LCLIGN+1:LCLIGN+NBC) = KCARX(N1:N2)
      LCLIGN = LCLIGN + NBC
      END
