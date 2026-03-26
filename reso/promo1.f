      SUBROUTINE PROMO1( NCODSA , NTDL , MUDL , AGP ,
     &                   LPLIGN , LPCOLO , AGM )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CONDENSER UNE STRUCTURE PROFIL EN STRUCTURE MORSE
C -----
C ENTREES :
C ---------
C NCODSA : CODDE DE STOCKAGE DE LA MATRICE
C NTDL   : NOMBRE DE LIGNES DE LA MATRICE
C MUDL   : POINTEUR SUR LES COEFFICIENTS DIAGONAUX
C AGP    : LA MATRICE PROFIL
C
C SORTIES :
C ---------
C LPLIGN : POINTEUR SUR LES LIGNES
C LPCOLO : POINTEUR SUR LES COLONNES
C AGM    : LA MATRICE MORSE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY     ANALYSE NUMERIQUE UPMC PARIS  MAI 1990
C23456---------------------------------------------------------------012
      DIMENSION         MUDL(NTDL+1), LPLIGN(NTDL+1), LPCOLO(1)
      DOUBLE PRECISION  AGP(1), AGM(1)
      COMMON / EPSSSS / EPZERO,EPSXYZ
      INTRINSIC         ABS
C
C     AIGUILLAGE SELON NCODSA
C     -----------------------
C
      IF(NCODSA.EQ.0) THEN
C
C        A  DIAGONALE
C        ------------
C
         LPLIGN(1)=0
         DO 1 I=1,NTDL
            LPLIGN(I+1) = I
            LPCOLO(I)   = I
            AGM(I)      = AGP(I)
 1       CONTINUE
C
      ELSE IF (NCODSA.EQ.1) THEN
C
C        A SYMETRIQUE
C        ------------
C
         JAP=0
         JAM=0
         LPLIGN(1)=0
         DO 2 J=1,NTDL
           JMI=MUDL(J+1)-MUDL(J)
           IF(JMI.LE.1) GOTO 3
           I=J-JMI
           JMI=JMI-1
           DO 4 JD=1,JMI
              JAP=JAP+1
              I=I+1
              IF( ABS( AGP(JAP) ) .GT. EPZERO ) THEN
                 JAM=JAM+1
                 LPCOLO(JAM) = I
                 AGM(JAM)    = AGP(JAP)
              ENDIF
  4      CONTINUE
  3      JAP=JAP+1
         JAM=JAM+1
         LPLIGN(J+1)=JAM
         LPCOLO(JAM)=J
         AGM(JAM)=AGP(JAP)
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
