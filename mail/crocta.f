      SUBROUTINE CROCTA( NS1, NS2, NS3, NS4, NS5, NS6,
     &                   NUOTPE, NUFILS, LARBRO,
     &                   NUOCTA, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CREATION D'UN OCTAEDRE DANS L'ARBRE LARBRO
C -----
C
C ENTREES:
C --------
C NS1,...,NS6 : LE NUMERO PTXYZD DES 6 SOMMETS DE L'OCTAEDRE
C NUOTPE : NUMERO LARBRO SI NUOTPE>0 ou LARBRT SI <0 DE L'OT PERE
C NUFILS : NUMERO DE FILS DE CET OCTAEDRE
C
C ENTREES ET SORTIES:
C -------------------
C LARBRO : ARBRE-14 DES OCTAEDRES ( FOND DE LA TETRAEDRISATION )
C      LARBRO(0,0) : NO DU 1-ER OCTAEDRE VIDE DANS LARBRO
C      LARBRO(1,0) : MAXIMUM DU 1-ER INDICE DE LARBRO (ICI -1:20)
C      LARBRO(2,0) : MAXIMUM DECLARE DU 2-EME INDICE DE LARBRO
C                    (ICI = MXARBO)
C
C      LARBRO(-1:20,1) : RACINE DE L'ARBRE (OCTAEDRE SANS PERE)
C
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
C SORTIES:
C --------
C NUOCTA : NUMERO (>0) LARBRO DE L'OCTAEDRE CREE
C IERR   : 0  SI PAS D'ERREUR,  >0 SI ERREUR
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC   NOVEMBRE 1992
C2345X7..............................................................012
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      INTEGER LARBRO(-1:20,0:*)
C
C     RECHERCHE D'UN OCTAEDRE VIDE
C     ----------------------------
      NUOCTA = LARBRO(0,0)
      IF( NUOCTA .EQ. 0 ) THEN
         WRITE(IMPRIM,*)
         WRITE(IMPRIM,*)
     %  'crocta: SATURATION DE L''ARBRE DES OCTAEDRES'
         WRITE(IMPRIM,*)
     %  'Nombre maximal de sommets de la tetraedrisation a AUGMENTER'
         WRITE(IMPRIM,*)
         IERR = 4
         RETURN
      ENDIF
C     L'OCTAEDRE VIDE SUIVANT
      LARBRO(0,0) = ABS( LARBRO(0,NUOCTA) )
C
C     LE PERE DANS LARBRO
C     -------------------
      LARBRO(-1,NUOCTA) = NUOTPE
C
C     LE NUMERO DE FILS POUR SON PERE ET SON TYPE D'OT (OCTA => 0)
C     ------------------------------------------------------------
      LARBRO(0,NUOCTA) = NUFILS
C
C     L'OCTAEDRE EST SUPPOSE NON DECOMPOSE :  0 POINT INTERNE
C     -------------------------------------------------------
      DO 10 I=1,14
         LARBRO(I,NUOCTA) = 0
 10   CONTINUE
C
C     LE NUMERO PTXYZD DES 6 SOMMETS
C     ------------------------------
      LARBRO(15,NUOCTA) = NS1
      LARBRO(16,NUOCTA) = NS2
      LARBRO(17,NUOCTA) = NS3
      LARBRO(18,NUOCTA) = NS4
      LARBRO(19,NUOCTA) = NS5
      LARBRO(20,NUOCTA) = NS6
      END
