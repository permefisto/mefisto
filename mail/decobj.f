      SUBROUTINE DECOBJ( NCOGEL, MODECO, TETASS, NBTETR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : DECOUPER PYRAMIDE PENTAEDRE OU HEXAEDRE SELON LE TYPE DE DECOUPAGE
C ----- DES FACES NE CONTENANT PAS LE PLUS PETIT NUMERO DE SOMMET DE L'EF
C       ALGORITHME:
C       LE SOMMET DE NO MINIMAL DE L'EF EST SOMMET DE DECOUPE DE SES 3 FACES
C       LE SOMMET DE NO MINIMAL DES 1 OU 3 FACES SANS LE SOMMET MINIMAL
C       EST SOMMET DE DECOUPE DE SES FACES QUADRANGULAIRES
C
C ENTREE :
C --------
C NCOGEL : CODE GEOMETRIQUE DE L'EF
C          5:TETRAEDRE, 6:PENTAEDRE, 7:HEXAEDRE, 9:PYRAMIDE A BASE CARREE
C MODECO : MODE DE DECOUPAGE DE LA FACE
C
C SORTIE :
C --------
C TETASS : NUMEROS DES 4 SOMMETS DES NBTETR TETRAEDRES
C NBTETR : NOMBRE  DE TETRAEDRES DE L'EF
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEURS: F.CHOUKROUN O.RICOU UPMC ANALYSE NUMERIQUE PARIS JANVIER 1990
C MODIFS : PERRONNET ALAIN ANALYSE NUMERIQUE     UPMC PARIS JANVIER 1999
C MODIFS : PERRONNET ALAIN Laboratoire J-L. LIONS UPMC PARIS   Mars 2007
C--------------------------------------------------------------------012
      INTEGER MODECO(3),TETASS(24),PY(8,2),PE(12,2),HE(24,2,2,2),NBTETR
C
      DATA (PY(I,1),I=1,8)  /1,2,3,5,    1,3,4,5 /
      DATA (PY(I,2),I=1,8)  /1,2,4,5,    2,3,4,5 /
C
      DATA (PE(I,1),I=1,12) /1,6,2,3,    1,5,2,6,    1,6,4,5/
      DATA (PE(I,2),I=1,12) /1,5,2,3,    1,5,3,6,    1,6,4,5/
C
      DATA (HE(I,1,1,1),I=1,24) /1,7,2,3,    1,6,2,7,    1,7,3,8,
     +                           1,8,3,4,    1,5,7,8,    1,6,7,5/
      DATA (HE(I,2,1,1),I=1,24) /1,7,2,3,    1,6,2,7,    1,7,3,8,
     +                           1,8,3,4,    1,6,7,8,    1,6,8,5/
C
      DATA (HE(I,1,2,1),I=1,24) /1,6,2,3,    1,7,3,8,    1,6,3,7,
     +                           1,8,3,4,    1,5,7,8,    1,6,7,5/
      DATA (HE(I,2,2,1),I=1,24) /1,6,2,3,    1,6,3,8,    1,8,3,4,
     +                           1,6,8,5,    3,6,7,8,    0,0,0,0/
C
      DATA (HE(I,1,1,2),I=1,24) /1,7,2,3,    1,6,2,7,    1,7,3,4,
     +                           1,5,7,8,    1,8,7,4,    1,6,7,5/
      DATA (HE(I,2,1,2),I=1,24) /1,7,2,3,    1,6,2,7,    1,7,3,4,
     +                           1,6,7,8,    1,8,7,4,    1,6,8,5/
C
      DATA (HE(I,1,2,2),I=1,24) /1,6,2,3,    1,6,3,7,    1,7,3,4,
     +                           1,5,7,8,    1,8,7,4,    1,6,7,5/
      DATA (HE(I,2,2,2),I=1,24) /1,6,2,3,    1,6,3,7,    1,7,3,4,
     +                           1,6,7,8,    1,8,7,4,    1,6,8,5/
C
      IF( NCOGEL .EQ. 9 ) THEN
C
C        PYRAMIDE
         NBTETR=2
         DO 5 I=1,8
            TETASS(I)=PY(I,MODECO(1))
 5       CONTINUE
C
      ELSE IF( NCOGEL .EQ. 6 ) THEN
C
C        PENTAEDRE
         NBTETR=3
         DO 10 I=1,12
            TETASS(I)=PE(I,MODECO(1))
 10      CONTINUE
C
      ELSE
C
C        HEXAEDRE
         IF( MODECO(1).EQ.2.AND.MODECO(2).EQ.2.AND.MODECO(3).EQ.1 ) THEN
            NBTETR=5
         ELSE
            NBTETR=6
         ENDIF
         DO 30 I=1,NBTETR*4
            TETASS(I)=HE(I,MODECO(1),MODECO(2),MODECO(3))
30       CONTINUE
      ENDIF
C
      RETURN
      END
