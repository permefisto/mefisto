      SUBROUTINE FILESKYLINE( FILENAME, NBROWS, LPDIAG, SKYLINE, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  PUT ON an ASCII FILE of NAME FILENAME A SKYLINE SYMMETRIC MATRIX
C -----
C
C INPUT :
C -------
C FILENAME: NAME OF THE ASCII FILE
C NBROWS  : NUMBER OF ROWS of THE SYMMETRIC MATRIX
C LPDIAG  : LPDIAG(0)=0
C           LPDIAG(i)= i-th diagonal coefficient pointer on vector SKYLINE
C           Only the coefficients between the first non zero and
C           the diagonal of each row are stored
C SKYLINE(1..LPDIAG(NBROWS)):  coefficients of the matrix SKYLINE
C
C OUTPUT:
C -------
C IERR  : 0 IF ZERO ERROR, ELSE  not ZERO => ERROR NUMBER
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR:ALAIN PERRONNET NCTS Tsing Hua University HSINCHU TAIWAN 9/2006
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      CHARACTER*(*)     FILENAME
      INTEGER           LPDIAG(0:NBROWS)
      DOUBLE PRECISION  SKYLINE(1:*)
C
C     OPEN the FILE
      CALL TRUNIT( NF )
      OPEN(FILE=FILENAME,UNIT=NF,ACCESS='SEQUENTIAL',
     %     FORM='FORMATTED',STATUS='UNKNOWN',
     %     IOSTAT= IERR )
      IF( IERR .NE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'NOM de FICHIER POUR LA MATRICE PROFIL'// FILENAME
         ELSE
            KERR(1) = 'UNCORRECT SKYLINE FILE NAME' // FILENAME
         ENDIF
         CALL LEREUR
         RETURN
      ENDIF
C
      IF( NBROWS .LE. 0 ) THEN
         IERR = 1
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'MATRICE PROFIL avec un NOMBRE de LIGNES <=0'
         ELSE
            KERR(1) = 'SKYLINE with NUMBER OF ROWS <=0'
         ENDIF
         CALL LEREUR
         RETURN
      ENDIF
C
      IF( LPDIAG(0) .NE. 0 ) THEN
         IERR = 2
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'MATRICE PROFIL avec LPDIAG(0) NON ZERO'
         ELSE
            KERR(1) = 'SKYLINE with LPDIAG(0) NOT EQUAL ZERO'
         ENDIF
         CALL LEREUR
         RETURN
      ENDIF
C
C     WRITE on the FILE
      WRITE(NF,10000)
10000 FORMAT('{ A SYMMETRIC SKYLINE MATRIX }'/
     %'{ NUMBER of ROWS, NUMBER of VARIABLES in the SKYLINE STORAGE }'/
     %'{ NUMBER of every DIAGONAL COEFFICIENT in SKYLINE }'/
     %'{ COEFFICIENTS of the SKYLINE STORAGE }')
      NBVAR = LPDIAG(NBROWS)
      WRITE(NF,10001) NBROWS, NBVAR
10001 FORMAT(2I12)
      WRITE(NF,10002) (LPDIAG(K),K=0,NBROWS)
10002 FORMAT(10I12)
      WRITE(NF,10003) (SKYLINE(K),K=1,NBVAR)
10003 FORMAT(5D25.17)
C
      CLOSE( NF )
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,*)'FICHIER',FILENAME,' EST CONSTRUIT'
      ELSE
         WRITE(IMPRIM,*)'FILE',FILENAME,' is BUILT'
      ENDIF
C
      RETURN
      END
