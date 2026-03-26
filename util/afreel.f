      SUBROUTINE AFREEL( REEL )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : AFFICHER LE REEL
C -----
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS   SEPTEMBRE 1988
C23456---------------------------------------------------------------012
      include"./incl/td.inc"
C.......................................................................
      PARAMETER (LREEL=14)
C     NOMBRE DE CARACTERES POUR REPRESENTER LE PLUS GRAND REEL
C     VOIR LES FORMAT NON PARAMETRES POUR ACCELERER LE CALCUL
      CHARACTER*(LREEL) KREEL
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
C
C     LA LIGNE EST BUFFERISEE A L'AIDE DU POINTEUR LCLIGN DERNIER
C     CARACTERE ENTRE DANS LA LIGNE
C
      IF( IMPRES*NOMUET .LE. 0 ) RETURN
C
      WRITE( KREEL , '(G14.7)' ) REEL
C
C     RECHERCHE DU 1-ER NON BLANC
      N1 = 0
 10   N1 = N1 + 1
      IF( KREEL(N1:N1) .EQ. ' ' ) GOTO 10
C
C     RECHERCHE DU DERNIER NON BLANC
      DO 20 N2=LREEL,N1,-1
         IF( KREEL(N2:N2) .NE. ' ' ) GOTO 30
 20   CONTINUE
C
C     LE NOMBRE DE CARACTERES A AFFICHER
 30   NBC = N2 + 2 - N1
      IF( NCLIGN-LCLIGN .LT. NBC ) THEN
         CALL AFLIGN
      ENDIF
C
      KLIGNE(LCLIGN+1:LCLIGN+NBC) = ' ' // KREEL(N1:N2)
      LCLIGN = LCLIGN + NBC
      END
