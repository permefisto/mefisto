      SUBROUTINE DVAMQT( MXSOMM, PXYD,
     %                   MXTRIA, NDTRIA, NOTRIA, CETRIA, NOTRSO,
     %                   NBCHGT )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    ECHANGE DES DIAGONALES DE 2 TRIANGLES ADJACENTS POUR
C -----    RENDRE LA TRIANGULATION DELAUNAY
C
C ENTREES:
C --------
C MXSOMM : NOMBRE MAXIMAL DE SOMMETS PERMIS POUR LA TRIANGULATION
C PXYD   : X Y DISTANCE SOUHAITEE DES SOMMETS
C MXTRIA : NOMBRE MAXIMAL DE TRIANGLES DECLARABLES DANS NOTRIA
C NDTRIA : NUMERO DU PLUS GRAND TRIANGLE UTILISE DANS NOTRIA
C NOTRIA : LISTE CHAINEE DES TRIANGLES
C                 ------- ------- ------- -------- -------- --------
C  PAR TRIANGLE : SOMMET1 SOMMET2 SOMMET3 TR_VOIS1 TR_VOIS2 TR_VOIS3
C                 ------- ------- ------- -------- -------- --------
C                 SOMMET    EST LE NUMERO DU SOMMET
C                 TR_VOIS i EST LE NUMERO DANS NOTRIA DU TRIANGLE
C                          ADJACENT PAR L'ARETE i
C CETRIA : COORDONNEES DU CENTRE DU CERCLE CIRCONSCRIT ET
C          CARRE DU RAYON
C NOTRSO : NOTRSO(I) NUMERO D'UN TRIANGLE AYANT POUR SOMMET I
C
C SORTIE :
C --------
C NBCHGT : NOMBRE DE CHANGEMENT DE DIAGONALE POUR 2 TRIANGLES ADJACENTS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC       MARS 1991
C....................................................................012
      INTEGER           NOTRIA(6,MXTRIA),
     %                  NOTRSO(MXSOMM)
      DOUBLE PRECISION  PXYD(3,MXSOMM),
     %                  CETRIA(1:3,1:MXTRIA)
C
C     PARCOURS DES TRIANGLES
C     ======================
      NBCHGT = 0
C
      DO 100 NT=1,NDTRIA
         IF( NOTRIA(1,NT) .LE. 0 ) GOTO 100
C
C        PARCOURS DES 3 ARETES DU TRIANGLE NT
         DO 50 I=1,3
C
C           LE TRIANGLE OPPOSE A L'ARETE I
            NT1 = NOTRIA(I+3,NT)
            IF( NT1 .LE. 0 ) GOTO 50
            IF( NOTRIA(1,NT1) .LE. 0 ) GOTO 50
C
C           L'ARETE I DE NT EST COMMUNE A 2 TRIANGLES
C           TENTATIVE POUR ECHANGER LA DIAGONALE => DELAUNAY
            CALL DV2T2T( I,      NT,
     %                   NOTRIA, CETRIA, NOTRSO, PXYD,
     %                   NT1,    NBCHG )
            NBCHGT = NBCHGT + NBCHG
C
 50      CONTINUE
C
 100  CONTINUE
      END
