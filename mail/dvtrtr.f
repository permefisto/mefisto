      SUBROUTINE DVTRTR( PXYD, NOTRIA, NUMTRI, NCTRIA, NCARET )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACER LE TRIANGLE NUMTRI
C -----
C ENTREES:
C --------
C PXYD   : TABLEAU DES COORDONNEES 2D DES POINTS
C          PAR POINT : X  Y  DISTANCE_SOUHAITEE
C NOTRIA : LISTE DES TRIANGLES
C                 ------- ------- ------- -------- -------- --------
C  PAR TRIANGLE : SOMMET1 SOMMET2 SOMMET3 TR_VOIS1 TR_VOIS2 TR_VOIS3
C                 ------- ------- ------- -------- -------- --------
C                 SOMMET    EST LE NUMERO DU SOMMET
C                 TR_VOIS i EST LE NUMERO DANS NOTRIA DU TRIANGLE
C                          ADJACENT PAR L'ARETE i
C          SOMMET1 = 0 => TRIANGLE VIDE
C NUMTRI : NUMERO DANS NOTRIA DU TRIANGLE
C NCTRIA : NUMERO DE LA COULEUR DU TRIANGLE A TRACER
C NCARET : NUMERO DE LA COULEUR DES ARETES DU TRIANGLE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC  SEPTEMBRE 1994
C....................................................................012
      include"./incl/trvari.inc"
      LOGICAL           TRATRI
      COMMON / DV2DCO / TRATRI
C     TRACE OU NON DES TRIANGLES GENERES DANS LA TRIANGULATION(CF DVTR2D)
      INTEGER           NOTRIA(6,*)
      DOUBLE PRECISION  PXYD(3,*)
      REAL              X(3),Y(3)
C
C     TEST DE VALIDITE
      IF( TRATRI ) THEN
C        LE TRIANGLE DOIT ETRE TRACE
         IF( NUMTRI .LE. 0 ) RETURN
         IF( NOTRIA(1,NUMTRI) .LE. 0 ) RETURN
C
         DO 10 J=1,3
            X(J) = REAL( PXYD(1,NOTRIA(J,NUMTRI)) )
            Y(J) = REAL( PXYD(2,NOTRIA(J,NUMTRI)) )
 10      CONTINUE
C
C        LE TRACE DE LA FACE ET DES ARETES
         CALL FACE2D( NCTRIA, NCARET, 3 , X , Y )
C
CCCC        TRACE DU NUMERO DES EXTREMITES
CCC         DO 20 J=1,3
CCC            CALL ENTIER2D( NCMAGE, X(J), Y(J), NOTRIA(J,NUMTRI) )
CCC 20      CONTINUE
      ENDIF
      END
