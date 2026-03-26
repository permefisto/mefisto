         SUBROUTINE  VOEXH3( NBX,   NBY,   NBZ, NBXYZ,
     S                       HEXYZ, CUXYZ, COSO )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    MAILLER UN HEXAEDRE DEFINI PAR SES FACES
C -----    PAR LA METHODE DE COOK-GORDON-HALL
C         (INTERPOLATION TRANSFINIE DE DEGRE 1)
C
C ENTREES :
C ---------
C NBX,NBY,NBZ : NOMBRE DE POINTS DANS CHAQUE DIRECTION
C NBXYZ       : NOMBRE DE POINTS DANS CHAQUE FACE (MAJORATION)
C HEXYZ       : COORDONNEES DES SOMMETS DES 6 FACES DE L'HEXAEDRE
C CUXYZ       : COORDONNEES DES SOMMETS DES 6 FACES DU CUBE UNITE
C
C ENTREES ET SORTIES :
C --------------------
C COSO        : COORDONNEES DES SOMMETS DU MAILLAGE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY      ANALYSE NUMERIQUE UPMC PARIS    JANVIER 1989
C MODIFS : ALAIN  PERRONNET ANALYSE NUMERIQUE UPMC PARIS   NOVEMBRE 1997
C23456---------------------------------------------------------------012
      REAL         HEXYZ(3,6,NBXYZ), CUXYZ(3,6,NBXYZ)
      REAL         COSO(3,NBX*NBY*NBZ)
      REAL         F(3,6), A(3,12)
      INTEGER      N(3)
      REAL         XYZ(3,0:6)
      EQUIVALENCE (XX,XYZ(1,0)), (YY,XYZ(2,0)), (ZZ,XYZ(3,0))
C
C     LES 8 SOMMETS DE L'HEXAEDRE
      NDECAL=NBX*NBY*(NBZ-1)
      N3 = NBX*NBY
      N4 = NBX*(NBY-1)+1
      N5 = NDECAL+1
      N6 = NDECAL+NBX
      N7 = NDECAL+NBX*NBY
      N8 = NDECAL+NBX*(NBY-1)+1
C
C     BOUCLE SUR LES POINTS INTERNES DE L'HEXAEDRE
      DO 1 NZ = 2 , NBZ - 1
         DO 2 NY = 2 , NBY - 1
            DO 3 NX = 2 , NBX - 1
C
C              LES COORDONNEES DES POINTS M1 ET M4 DES FACES 1 ET 4
               N(1) = NBX*(NY-1) + NX
C              LES COORDONNEES DES POINTS M2 ET M5 DES FACES 2 ET 5
               N(2) = NBX*(NZ-1) + NX
C              LES COORDONNEES DES POINTS M3 ET M6 DES FACES 3 ET 6
               N(3) = NBY*(NZ-1) + NY
               DO 4 NF=1,3
                  DO 6 J=1,3
                     XYZ( J, NF   ) = CUXYZ( J,  NF , N(NF) )
                     XYZ( J, NF+3 ) = CUXYZ( J, NF+3, N(NF) )
 6                CONTINUE
 4             CONTINUE
C
C              NUMERO DU POINT M A CREER DANS LE DOMAINE
               NUSO = NBX * NBY * (NZ-1) + NBX * (NY-1) + NX
C
C              CALCUL DU POINT LE PLUS PROCHE DES DROITES M1M4, M2M5, ET M3M6
               CALL MINDIS( 3, XYZ )
C
C              LE POINT D'INTERPOLATION SUR CHACUNE DES SIX FACES
               CALL VOEXH6( XYZ, CUXYZ, HEXYZ, NBX, NBY, NBZ, NBXYZ, F )
C
C              LE POINT D'INTERPOLATION SUR CHACUNE DES 12 ARETES
               CALL VOEXH8( XYZ, CUXYZ, HEXYZ, NBX, NBY, NBZ, A )
C
C              LE SOMMET INTERNE PAR LA FORMULE DE COOK, GORDON ET HALL
               XX1 = 1 - XX
               YY1 = 1 - YY
               ZZ1 = 1 - ZZ
               DO 5 J=1,3
                  COSO(J,NUSO) =
     %  XX1*( F(J,3)+ YY1*( ZZ1*COSO(J, 1) + ZZ*COSO(J,N5) - A(J,9)  )
     %              + YY *( ZZ1*COSO(J,N4) + ZZ*COSO(J,N8) - A(J,12) ))
     %+ YY1*( F(J,2)- ZZ1*A(J,1) - ZZ * A(J,5) )
     %+ ZZ1*( F(J,1)- XX1*A(J,4) - XX * A(J,2) )
     %+ XX *( F(J,6)+ YY1*( ZZ1*COSO(J,NBX) + ZZ*COSO(J,N6) - A(J,10) )
     %              + YY *( ZZ1*COSO(J,N3 ) + ZZ*COSO(J,N7) - A(J,11) ))
     %+ YY *( F(J,5)- ZZ1 * A(J,3) - ZZ * A(J,7) )
     %+ ZZ *( F(J,4)- XX1 * A(J,8) - XX * A(J,6) )
 5             CONTINUE
 3          CONTINUE
 2       CONTINUE
 1    CONTINUE
C
      RETURN
      END
