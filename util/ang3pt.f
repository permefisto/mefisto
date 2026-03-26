      SUBROUTINE ANG3PT( PT1 , PT2 , PT3 , ANGLE )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCULER L'ANGLE EN RADIANS DE -PI A PI
C ----- DU SECTEUR CENTRE EN PT2 (PT2-PT1 , PT2-PT3)
C
C              PT3                   PT1
C                 \   ANGLE       /
C                  \  <====    /
C                   \       /
C                    \ PT2
C ENTREES:
C --------
C PT1 PT2 PT3 : LES 2 COORDONNEES DES 3 POINTS
C
C SORTIE :
C --------
C ANGLE  : EN RADIANS DE -PI A PI DE L'ANGLE PT2 PT1,PT2 PT3
c          0 SI 2 DES 3 POINTS SONT CONFONDUS
C            AVEC IMPRESSION D'UN DIAGNOSTIC DANS LE SP ANRAL1
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS SEPTEMBRE 1986
C..............................................................................
      PARAMETER ( PI2 = 3.14159265358979312 * 2. )
      REAL        PT1(2),PT2(2),PT3(2)
C
C     ANGLE ( OX , PT2-PT1 )  DE -PI A PI
      CALL ANRAL1( PT2 , PT1 , ANGLE , S1 )
C     ANGLE ( OX , PT2-PT3 )  DE -PI A PI
      CALL ANRAL1( PT2 , PT3 , ANGLE , S2 )
C     ANGLE ( PT2-PT1 , PT2-PT3 ) EVENTUELLEMENT NEGATIF
      ANGLE = S2 - S1
      END
