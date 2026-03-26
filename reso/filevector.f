      SUBROUTINE FILEVECTOR( FILENAME, NBCOMP, NBVECT, VECTORS, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  PUT ON an ASCII FILE of NAME FILENAME NBVECT GLOBAL VECTORS
C -----
C
C INPUT :
C -------
C FILENAME: NAME OF THE ASCII FILE
C NBCOMP  : NUMBER OF COMPONENTS of the VECTORS
C NBVECT  : NUMBER OF VECTORS (>=1)
C VECTORS(1:NBCOMP,1:NBVECT):  coefficients of VECTORS
C
C OUTPUT:
C -------
C IERR  : 0 IF ZERO ERROR, ELSE  not ZERO => ERROR NUMBER
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET TIMS NATIONAL TAIWAN UNIVERSITY    August 2007
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      CHARACTER*(*)     FILENAME
      DOUBLE PRECISION  VECTORS(1:NBCOMP,1:NBVECT)
C
C     OPEN the FILE
      CALL TRUNIT( NF )
      OPEN(FILE=FILENAME,UNIT=NF,ACCESS='SEQUENTIAL',
     %     FORM='FORMATTED',STATUS='UNKNOWN',
     %     IOSTAT= IERR )
      IF( IERR .NE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'NOM INCORRECT de FICHIER pour les VECTEURS'
     %              // FILENAME
         ELSE
            KERR(1) = 'UNCORRECT VECTORS FILE NAME' // FILENAME
         ENDIF
         CALL LEREUR
         RETURN
      ENDIF
C
      IF( NBCOMP .LE. 0 ) THEN
         IERR = 1
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'NOMBRE de COMPOSANTES d''un VECTEUR <=0'
         ELSE
            KERR(1) = 'NUMBER of COMPONENTS of one VECTOR <=0'
         ENDIF
         CALL LEREUR
         RETURN
      ENDIF
C
      IF( NBVECT .LE. 0 ) THEN
         IERR = 1
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'NOMBRE de VECTEURS <=0'
         ELSE
            KERR(1) = 'NUMBER of VECTORS <=0'
         ENDIF
         CALL LEREUR
         RETURN
      ENDIF
C
C     WRITE on the FILE
      WRITE(NF,10000)
10000 FORMAT('{ GLOBAL VECTORS }'/
     %'{ NUMBER of COMPONENTS, NUMBER of VECTORS in the STORAGE }'/
     %'{ COEFFICIENTS of the VECTORS(1:NBCOMP,1:NBVECT) }')
      WRITE(NF,10001) NBCOMP, NBVECT
10001 FORMAT(2I12)
C
C     COEFFICIENTS of VECTORS(1:NBCOMP,1:NBVECT)
      DO 10 J=1,NBVECT
         WRITE(NF,10010) (VECTORS(I,J),I=1,NBCOMP)
 10   CONTINUE
10010 FORMAT(5D25.17)
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
