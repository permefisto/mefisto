      SUBROUTINE AUGMENTER
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   AUGMENTER LE TRACE EN REDUISANT L'HEXAEDRE ENGLOBANT
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 2002
C ......................................................................
      include"./incl/xyzext.inc"
C
      IF( COOEXT(1,1) .GE. RINFO( 'GRAND' ) ) THEN
C
C        COOEXT INCORRECT POUR LA SUITE => C'EST LE CARRE UNITE
         COOEXT(1,1) = 0.0
         COOEXT(1,2) = 1.0
         COOEXT(2,1) = 0.0
         COOEXT(2,2) = 1.0
         COOEXT(3,1) = 0.0
         COOEXT(3,2) = 0.0
C
      ELSE
C
         DO 10 I=1,3
            R = (COOEXT(I,2) - COOEXT(I,1)) * 0.111
            COOEXT(I,1) =  COOEXT(I,1) + R
            COOEXT(I,2) =  COOEXT(I,2) - R
 10      CONTINUE
C
      ENDIF
      END
