      SUBROUTINE CRTETR( NS1, NS2, NS3, NS4,
     &                   NUOTPE, NUFILS, NUTYFI, LARBRT,
     &                   NUTETR, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CREATION D'UN TETRAEDRE DANS L'ARBRE LARBRT
C -----
C
C ENTREES:
C --------
C NS1,...,NS4 : LE NUMERO PTXYZD DES 4 SOMMETS DU TETRAEDRE
C NUOTPE : NUMERO LARBRT SI NUOTPE>0 ou LARBRT SI <0 DU OT PERE
C NUFILS : NUMERO DE FILS DE CE TETRAEDRE
C NUTYFI : NUMERO 1 ou 2 DU TYPE DE CE TETRAEDRE
C
C ENTREES ET SORTIES:
C -------------------
C LARBRT : ARBRE-5 DES TETRAEDRES ( FOND DE LA TETRAEDRISATION )
C      LARBRT(0,0) : NO DU 1-ER TETRAEDRE VIDE DANS LARBRT
C      LARBRT(1,0) : MAXIMUM DU 1-ER INDICE DE LARBRT (ICI -1:9)
C      LARBRT(2,0) : MAXIMUM DECLARE DU 2-EME INDICE DE LARBRT
C                     (ICI = MXARBT)
C
C      LARBRO(-1,J) : NO DU PERE DU TETRAEDRE J DANS UN DES 2 ARBRES
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
C SORTIES:
C --------
C NUTETR : NUMERO (<0) LARBRT DU TETRAEDRE CREE
C IERR   : 0  SI PAS D'ERREUR,  >0 SI ERREUR
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC   NOVEMBRE 1992
C2345X7..............................................................012
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      INTEGER LARBRT(-1:9,0:*)
C
C     RECHERCHE D'UN TETRAEDRE VIDE
C     -----------------------------
      NUTETR = LARBRT(0,0)
      IF( NUTETR .EQ. 0 ) THEN
         WRITE(IMPRIM,*)
         WRITE(IMPRIM,*)
     %  'crtetr: SATURATION DE L''ARBRE DES TETRAEDRES'
         WRITE(IMPRIM,*)
     %  'Nombre maximal de sommets de la tetraedrisation a AUGMENTER'
         WRITE(IMPRIM,*)
         IERR = 5
         RETURN
      ENDIF
C     LE TETRAEDRE VIDE SUIVANT
      LARBRT(0,0) = ABS( LARBRT(0,NUTETR) )
C
C     LE PERE DANS LARBRT
C     -------------------
      LARBRT(-1,NUTETR) = NUOTPE
C
C     LE NUMERO DE FILS POUR SON PERE ET SON TYPE D'OT (TETRA => 0)
C     -------------------------------------------------------------
      LARBRT(0,NUTETR) = NUFILS + 100 * NUTYFI
C
C     LE TETRAEDRE EST SUPPOSE NON DECOMPOSE :  0 POINT INTERNE
C     ---------------------------------------------------------
      DO 10 I=1,5
         LARBRT(I,NUTETR) = 0
 10   CONTINUE
C
C     LE NUMERO PTXYZD DES 4 SOMMETS
C     ------------------------------
      LARBRT(6,NUTETR) = NS1
      LARBRT(7,NUTETR) = NS2
      LARBRT(8,NUTETR) = NS3
      LARBRT(9,NUTETR) = NS4
C
C     TETRAEDRE STOCKE DANS LARBRT
      NUTETR = -NUTETR
      END
