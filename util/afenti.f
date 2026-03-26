      SUBROUTINE AFENTI( ENTIER )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : AFFICHER ENTIER
C -----
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS SEPTEMBRE 1988
C23456---------------------------------------------------------------012
      include"./incl/td.inc"
C.......................................................................
      PARAMETER (LENTIE=12)
C     NOMBRE DE CARACTERES POUR REPRESENTER LE PLUS GRAND ENTIER
C     VOIR LES FORMAT NON PARAMETRES POUR ACCELERER LE CALCUL
      CHARACTER*(LENTIE) KENTIE
      INTEGER ENTIER
C
C     LA LIGNE EST BUFFERISEE A L'AIDE DU POINTEUR LCLIGN DERNIER
C     CARACTERE ENTRE DANS LA LIGNE
C
      IF( IMPRES*NOMUET .LE. 0 ) RETURN
C
      WRITE(KENTIE,'(I12)' ) ENTIER
C
C     RECHERCHE DU 1-ER NON BLANC
      N = 0
 10   N = N + 1
      IF( KENTIE(N:N) .EQ. ' ' ) GOTO 10
C
C     LE NOMBRE DE CARACTERES A AFFICHER
      NBC = LENTIE + 2 - N
      IF( NCLIGN-LCLIGN .LT. NBC ) THEN
         CALL AFLIGN
      ENDIF
C
      KLIGNE(LCLIGN+1:LCLIGN+NBC) = ' ' // KENTIE(N:LENTIE)
      LCLIGN = LCLIGN + NBC
      END
