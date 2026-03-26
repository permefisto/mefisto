      SUBROUTINE DPTSOT( PT,     PTXYZD, NUOT, LARBRO, LARBRT,
     %                   DISMIN, NUSMIN )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    RETROUVER LE PLUS PROCHE SOMMET DE L'OT NUOT POUR LE POINT PT
C -----
C
C ENTREES:
C --------
C PT     : X Y Z DU POINT
C PTXYZD : PAR POINT : X  Y  Z  DISTANCE_SOUHAITEE
C LARBRO : ARBRE-14 DES OCTAEDRES ( FOND DE LA TETRAEDRISATION )
C NUOT   : NUMERO DE L'OT (>0 DANS LARBRO, <0 DANS LARBRT)
C      LARBRO(-1:20,1) : RACINE DE L'ARBRE (OCTAEDRE SANS PERE)
C      LARBRO(-1,J) : NO DU PERE DE L'OCTAEDRE J DANS UN DES 2 ARBRES
C                     >0 => DANS LARBRO
C                     <0 => DANS LARBRT
C      LARBRO(0,J)  : 1 A 14 NO DE FILS DE L'OCTAEDRE J POUR SON PERE
C                     + 100 * NO TYPE DE L'OT J
C                     NO TYPE DE L'OT : 0 SI OCTAEDRE
C                                       1 SI TETRAEDRE T RONDE (T1)
C                                       2 SI TETRAEDRE T       (T2)
C   SI LARBRO(0,J)>0 ALORS J EST UN OCTAEDRE OCCUPE
C      SI LARBRO(1,.)>0 ALORS
C         LARBRO(1:14,J): NO (>0) LARBRO DES 14 SOUS-OCTA-TETRAEDRES
C      SINON
C         LARBRO(1:14,J):-NO PTXYZD DES 0 A 14 POINTS INTERNES DE L'OCTA J
C                         0  SI PAS DE POINT
C                        ( J EST ALORS UNE FEUILLE DE L'ARBRE )
C
C      LARBRO(15:20,J) : NO PTXYZD DES 6 SOMMETS DE L'OCTAEDRE J
C   SINON
C      LARBRO(0,J): -ADRESSE DANS LARBRO DE L'OCTAEDRE VIDE SUIVANT
C
C LARBRT : ARBRE-5 DES TETRAEDRES ( FOND DE LA TETRAEDRISATION )
C      LARBRT(0,0) : NO DU 1-ER TETRAEDRE VIDE DANS LARBRT
C      LARBRT(1,0) : MAXIMUM DU 1-ER INDICE DE LARBRT (ICI -1:9)
C      LARBRT(2,0) : MAXIMUM DECLARE DU 2-EME INDICE DE LARBRT
C                     (ICI = MXARBT)
C
C      LARBRT(-1,J) : NO DU PERE DU TETRAEDRE J DANS UN DES 2 ARBRES
C                     >0 => DANS LARBRO
C                     <0 => DANS LARBRT
C      LARBRT(0,J) : 0 A 4 NO DE FILS DU TETRAEDRE J POUR SON PERE
C                    + 100 * NO TYPE DE L'OT J
C                     NO TYPE DE L'OT : 0 SI OCTAEDRE
C                                       1 SI TETRAEDRE T RONDE (T1)
C                                       2 SI TETRAEDRE T       (T2)
C
C   SI LARBRT(0,J)>0 ALORS J EST UN TETRAEDRE OCCUPE
C      SI LARBRT(1,J)>0 ALORS
C         LARBRT(1:5,J): NO (>0) LARBRT DES 5 SOUS-OCTA-TETRAEDRES
C      SINON
C         LARBRT(1:5,J):-NO PTXYZD DES 0 A 5 POINTS INTERNES AU TETRA J
C                         0  SI PAS DE POINT
C                        ( J EST ALORS UNE FEUILLE DU ARBRE )
C
C      LARBRT(6:9,J) : NO PTXYZD DES 4 SOMMETS DU TETRAEDRE J
C   SINON
C      LARBRT(0,J): ADRESSE DANS LARBRT DU TETRAEDRE VIDE SUIVANT
C
C
C SORTIES:
C --------
C DISMIN : DISTANCE MINIMALE DES SOMMETS AU POINT PT
C NUSMIN : NUMERO (1 A 4) DU SOMMET DE L'OT DE DISTANCE MIN AU POINT PT
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 1992
C2345X7..............................................................012
      DOUBLE PRECISION  PT(4),PTXYZD(4,*)
      INTEGER           LARBRO(-1:20,0:*),
     &                  LARBRT(-1:9,0:*)
C
C     INITIALISATION
      DISMIN = RINFO('GRAND')
C
      IF( NUOT .GT. 0 ) THEN
C        OCTAEDRE
         DO 10 I=1,6
            NP = ABS( LARBRO(14+I,NUOT) )
            D = REAL( ( PT(1) - PTXYZD(1,NP) ) ** 2
     %              + ( PT(2) - PTXYZD(2,NP) ) ** 2
     %              + ( PT(3) - PTXYZD(3,NP) ) ** 2 )
            IF( D .LT. DISMIN ) THEN
C              NOUVEAU MINIMUM
               DISMIN = D
               NUSMIN = I
            ENDIF
 10      CONTINUE
      ELSE
C        TETRAEDRE
         DO 20 I=1,4
            NP = LARBRT(5+I,-NUOT)
            D = REAL( ( PT(1) - PTXYZD(1,NP) ) ** 2
     %              + ( PT(2) - PTXYZD(2,NP) ) ** 2
     %              + ( PT(3) - PTXYZD(3,NP) ) ** 2 )
            IF( D .LT. DISMIN ) THEN
C              NOUVEAU MINIMUM
               DISMIN = D
               NUSMIN = I
            ENDIF
 20      CONTINUE
      ENDIF
C
C     LA RACINE CARREE
      DISMIN = SQRT( DISMIN )
C
      RETURN
      END
