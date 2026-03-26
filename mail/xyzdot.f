      SUBROUTINE XYZDOT( PT,    PTXYZD, NUOT, LARBRO, LARBRT,
     %                   NPMIN, DISMI2 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CHERCHER SI L'OT NUOT A UN SOMMET OU UN POINT S'IL
C -----    S'AGIT D'UN OT FEUILLE PLUS PROCHE QUE CELUI ACTUEL

C ENTREES:
C --------
C PT     : X Y Z DU POINT
C PTXYZD : PAR POINT : X  Y  Z  DISTANCE_SOUHAITEE
C NUOT   : LE NUMERO D'OT (>0 DANS LARBRO, <0 DANS LARBRT)
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

C LARBRT : ARBRE-5 DES TETRAEDRES ( FOND DE LA TETRAEDRISATION )
C      LARBRT(0,0) : NO DU 1-ER TETRAEDRE VIDE DANS LARBRT
C      LARBRT(1,0) : MAXIMUM DU 1-ER INDICE DE LARBRT (ICI -1:9)
C      LARBRT(2,0) : MAXIMUM DECLARE DU 2-EME INDICE DE LARBRT
C                     (ICI = MXARBT)

C      LARBRT(-1,J) : NO DU PERE DU TETRAEDRE J DANS UN DES 2 ARBRES
C                     >0 => DANS LARBRO
C                     <0 => DANS LARBRT
C      LARBRT(0,J) : 0 A 4 NO DE FILS DU TETRAEDRE J POUR SON PERE
C                    + 100 * NO TYPE DE L'OT J
C                     NO TYPE DE L'OT : 0 SI OCTAEDRE
C                                       1 SI TETRAEDRE T RONDE (T1)
C                                       2 SI TETRAEDRE T       (T2)

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

C MODIFIES:
C ---------
C NPMIN  : NUMERO PTXYZD DU POINT (SOMMET OU NON DE L'OT) LE PLUS PROCHE
C          0 S'IL N'EN EXISTE PAS
C DISMI2 : CARRE DE LA DISTANCE DE NPMIN AU POINT PT SI NPMIN>0
C          0 SI PT A ETE IDENTIFIE AU POINT NPMIN>0 DE PTXYZD
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 1992
C2345X7..............................................................012
      DOUBLE PRECISION  PT(4), PTXYZD(4,*), D
      INTEGER           LARBRO(-1:20,0:*),
     &                  LARBRT(-1:9,0:*)

      IF( NUOT .GT. 0 ) THEN

C        OCTAEDRE
C        ========
         IF( LARBRO(1,NUOT) .LT. 0 ) THEN
C           OCTAEDRE FEUILLE PARCOURS DES POINTS ET SOMMETS
            I1 = 1
         ELSE
C           OCTAEDRE NON FEUILLE PARCOURS DES SEULS SOMMETS
            I1 = 15
         ENDIF

         DO 10 I=I1,20
            NP = ABS( LARBRO(I,NUOT) )
            IF( NP .EQ. 0 ) GOTO 10
            D = ( PT(1) - PTXYZD(1,NP) ) ** 2
     %        + ( PT(2) - PTXYZD(2,NP) ) ** 2
     %        + ( PT(3) - PTXYZD(3,NP) ) ** 2
            IF( D .LT. DISMI2 ) THEN
C              NOUVEAU MINIMUM
               DISMI2 = REAL( D )
               NPMIN  = NP
            ENDIF
 10      CONTINUE

      ELSE

C        TETRAEDRE
C        =========
         IF( LARBRT(1,-NUOT) .LT. 0 ) THEN
C           TETRAEDRE FEUILLE PARCOURS DES POINTS ET SOMMETS
            I1 = 1
         ELSE
C           TETRAEDRE NON FEUILLE PARCOURS DES SEULS SOMMETS
            I1 = 6
         ENDIF

         DO 20 I=I1,9
            NP = ABS( LARBRT(I,-NUOT) )
            IF( NP .EQ. 0 ) GOTO 20
            D = ( PT(1) - PTXYZD(1,NP) ) ** 2
     %        + ( PT(2) - PTXYZD(2,NP) ) ** 2
     %        + ( PT(3) - PTXYZD(3,NP) ) ** 2
            IF( D .LT. DISMI2 ) THEN
C              NOUVEAU MINIMUM
               DISMI2 = REAL( D )
               NPMIN  = NP
            ENDIF
 20      CONTINUE
      ENDIF

      IF( NPMIN .GT. 0 ) THEN
C        IDENTIFICATION POSSIBLE DE PT AVEC PTXYZD(NPMIN)?
         CALL XYZIDD( PT, PTXYZD(1,NPMIN), NONOUI )
         IF( NONOUI .GT. 0 ) THEN
C           OUI: LA DISTANCE ENTRE LES 2 POINTS EST DONC NULLE
            DISMI2 = 0
         ENDIF
      ENDIF
C
      RETURN
      END
