      SUBROUTINE EFSURFNE( XYZSOM, NBSEF, NOSOEL, SURFACE )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    AFFICHER L'EF AVEC UNE SURFACE NEGATIVE
C -----
C
C ENTREES:
C --------
C XYZSOM : COORDONNEES DES SOMMETS DU MAILLAGE
C NBSEF  : NOMBRE DE SOMMETS DE L'EF
C NOSOEL : NO DES NBSEF SOMMETS DE L'EF
C SURFACE: SURFACE<0 DE L'EF
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Laboratoire J.L. LIONS UPMC Paris   Mars 2007
C23456+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++012
      include"./incl/langue.inc"
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      INTEGER           NOSOEL(NBSEF)
      REAL              XYZSOM(3,*)

      IF( NBSEF .EQ. 3 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,*) 'TRIANGLE ',(NOSOEL(I),I=1,3),
     %                      ' de SURFACE ',SURFACE,' <=0!'
         ELSE
            WRITE(IMPRIM,*) 'TRIANGLE ',(NOSOEL(I),I=1,3),
     %                      ' of SURFACE ',SURFACE,' <=0!'
         ENDIF
      ELSE
         IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,*) 'QUADRANGLE ',(NOSOEL(I),I=1,NBSEF),
     %               ' de TRIANGLES 123+134 de SURFACE ',SURFACE,' <=0!'
         ELSE
         WRITE(IMPRIM,*) 'QUADRANGLE ',(NOSOEL(I),I=1,NBSEF),
     %              ' has TRIANGLES 123+134 of SURFACE ',SURFACE,' <=0!'
         ENDIF
      ENDIF

      DO I=1,NBSEF
         WRITE(IMPRIM,*) 'Pt',NOSOEL(I),
     %                   ' XYZ=',(XYZSOM(K,NOSOEL(I)),K=1,3)
      ENDDO

      RETURN
      END
