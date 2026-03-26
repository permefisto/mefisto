      SUBROUTINE COMMO1( NCODSA, LPLIGN, LPCOLO, AG )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    COMPRIMER UNE STRUCTURE MORSE
C -----
C ENTREES :
C ---------
C NCODSA : CODDE DE STOCKAGE DE LA MATRICE
C NTDL   : NOMBRE DE LIGNES DE LA MATRICE
C LPLIGN : POINTEUR SUR LES LIGNES
C LPCOLO : POINTEUR SUR LES COLONNES
C AG     : LA MATRICE MORSE
C NBNOSD : NOMBRE DE NOEUDS SUR L'INTERFACE
C
C SORTIES :
C ---------
C LPLIGN : POINTEUR SUR LES LIGNES MODIFIE
C LPCOLO : POINTEUR SUR LES COLONNES MODIFIE
C AG     : LA MATRICE MORSE MODIFIEE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY     ANALYSE NUMERIQUE UPMC PARIS  MAI 1990
C23456---------------------------------------------------------------012
      INTEGER           LPLIGN(NTDL+1),LPCOLO(1)
      DOUBLE PRECISION  AG(1)
      COMMON / SOUSDO / NORESO,NTDL,NDSM,NDIM,NPIMAX
C
      IF( NCODSA .EQ. 0 ) THEN
C
C        A  DIAGONALE
C        ------------
C
C      CAS NON TRAITE
C
C
      ELSE IF ( NCODSA .EQ. 1 ) THEN
C
C        A SYMETRIQUE
C        ------------
C
         JAM=0
         JAP=0
         DO 2 J=1,NTDL
            JAF=LPLIGN(J+1)
            DO 3 JA=JAP+1,JAF
C               DO 4 NB=1,NBNOSD
C                  NUNOK=NDIM*(LINOSD(NB)-1)
C                  DO 4 K=1,NDIM
C                     IF (LPCOLO(JA) .EQ. NUNOK+K) THEN
                        JAM=JAM+1
                        LPCOLO(JAM) = LPCOLO(JA)
                        AG(JAM) = AG(JA)
C                        GO TO 3
C                     ENDIF
C  4           CONTINUE
  3         CONTINUE
            JAP=JAF
            LPLIGN(J+1)=JAM
  2      CONTINUE
         LPLIGN(1)=0
C
      ELSE IF (NCODSA.EQ.-1) THEN
C
C      A NON-SYMETRIQUE
C      ----------------
C
C      CAS NON TRAITE
C
      ENDIF
C
      RETURN
      END
