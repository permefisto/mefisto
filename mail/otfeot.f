      SUBROUTINE OTFEOT( NUOTIN, LARBRO, LARBRT,
     %                   MXFILS, MNFILS,
     %                   MXOTFL, NBOTFL, MNOTFL )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    RETROUVER LES FEUILLES DE L'OT NUOTIN
C -----
C
C ENTREES:
C --------
C NUOTIN : NUMERO DE L'OT DE FEUILLES A TROUVER
C PTXYZD : PAR POINT : X  Y  Z  DISTANCE_SOUHAITEE
C LARBRO : ARBRE-14 DES OCTAEDRES ( FOND DE LA TETRAEDRISATION )
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
C MODIFIES:
C ---------
C MXFILS : NOMBRE D'ENTIERS DECLARES DANS LE TABLEAU FILS
C MNFILS : ADRESSE MCN DES FILS DE L'OT
C MXOTFL : NOMBRE D'ENTIERS DECLARES DANS LE TABLEAU OTFL
C MNOTFL : ADRESSE MCN DES NBOTFL FEUILLES  DE L'OT
C
C SORTIE :
C --------
C NBOTFL : NOMBRE DE FEUILLES DE L'OT
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 1992
C2345X7..............................................................012
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      INTEGER           LARBRO(-1:20,0:*),
     &                  LARBRT(-1:9,0:*)
C
C     INITIALISATION DE LA PILE :  L'OT EST EMPILE
      NBOTFL = 0
      MCN(MNFILS) = NUOTIN
      LHFILS = 0
C
C     TANT QU'IL EXISTE UN OT FILS
 1    IF( LHFILS .GE. 0 ) THEN
C        L'OT A TRAITER
         NUOT   = MCN( MNFILS + LHFILS )
         LHFILS = LHFILS - 1
         IF( NUOT .GT. 0 ) THEN
C           OCTAEDRE
            IF( LARBRO(1,NUOT) .GT. 0 ) THEN
C              LES FILS SONT EMPILES
               IF( LHFILS+14 .GE. MXFILS ) THEN
C                 PILE SATUREE => ELLE EST AGRANDIE
                  CALL TNMCAU( 'ENTIER', MXFILS, MXFILS+128, LHFILS,
     &                          MNFILS )
                  MXFILS = MXFILS + 128
               ENDIF
               DO 10 I=1,14
                  LHFILS = LHFILS + 1
                  MCN( MNFILS + LHFILS ) = LARBRO( I, NUOT )
 10            CONTINUE
               GOTO 1
            ENDIF
            GOTO 50
         ELSE IF( NUOT .LT. 0 ) THEN
C           TETRAEDRE
            NUOTT = -NUOT
            IF( LARBRT(1,NUOTT) .GT. 0 ) THEN
C              LES FILS SONT EMPILES
               IF( LHFILS+5 .GE. MXFILS ) THEN
C                 PILE SATUREE => ELLE EST AGRANDIE
                  CALL TNMCAU( 'ENTIER', MXFILS, MXFILS+128, LHFILS,
     &                          MNFILS )
                  MXFILS = MXFILS + 128
               ENDIF
               DO 20 I=1,5
                  LHFILS = LHFILS + 1
                  MCN( MNFILS + LHFILS ) = LARBRT( I, NUOTT )
 20            CONTINUE
               GOTO 1
            ENDIF
            GOTO 50
         ELSE
            GOTO 1
         ENDIF
C
C        L'OT EST UNE FEUILLE A AJOUTER
 50      IF( NBOTFL .GE. MXOTFL ) THEN
C           PILE SATUREE => ELLE EST AGRANDIE
            CALL TNMCAU( 'ENTIER', MXOTFL, MXOTFL+128, NBOTFL,
     &                    MNOTFL )
            MXOTFL = MXOTFL + 128
         ENDIF
         MCN( MNOTFL + NBOTFL ) = NUOT
         NBOTFL = NBOTFL + 1
         GOTO 1
      ENDIF
      END
