      SUBROUTINE EFVONE( XYZSOM, NBSEF, NOSOEL, VOLUME )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    AFFICHER l'EF NOSOEL de VOLUME NEGATIF ou NUL
C -----
C
C ENTREES:
C --------
C XYZSOM : COORDONNEES DES SOMMETS DU MAILLAGE
C NBSEF  : NOMBRE DE SOMMETS DE L'EF
C NOSOEL : NO DES NBSEF SOMMETS DE L'EF
C VOLUME: VOLUME<0 DE L'EF DE SOMMETS NOSOEL
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Laboratoire J.L. LIONS UPMC Paris   Mars 2007
C23456+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++012
      include"./incl/langue.inc"
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      INTEGER           NOSOEL(NBSEF)
      REAL              XYZSOM(3,*)

      IF( NBSEF .EQ. 4 ) THEN

C        TETRAEDRE
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,*)'efvone: TETRAEDRE de SOMMETS',
     %                     (NOSOEL(I),I=1,NBSEF),
     %                     ' de VOLUME ',VOLUME,' <=0!'
         ELSE
            WRITE(IMPRIM,*)'efvone: TETRAHEDRON VERTICES',
     %                      (NOSOEL(I),I=1,NBSEF),
     %                     ' of VOLUME ',VOLUME,' <=0!'
         ENDIF

      ELSE IF( NBSEF .EQ. 5 ) THEN

C        PYRAMIDE
         IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,*) 'efvone: PYRAMIDE ',(NOSOEL(I),I=1,NBSEF),
     %                 ' de TETRAEDRE 1 2 3 4 de VOLUME ',VOLUME,' <=0!'
         ELSE
         WRITE(IMPRIM,*) 'efvone: PYRAMID ',(NOSOEL(I),I=1,NBSEF),
     %              ' of TETRAHEDRON 1 2 3 4 of VOLUME ',VOLUME,' <=0!'
         ENDIF

      ELSE IF( NBSEF .EQ. 6 ) THEN

C        PENTAEDRE
         IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,*) 'efvone: PENTAEDRE ',(NOSOEL(I),I=1,NBSEF),
     %                 ' de TETRAEDRE 1 2 3 4 de VOLUME ',VOLUME,' <=0!'
         ELSE
         WRITE(IMPRIM,*) 'efvone: PENTAHEDRON ',(NOSOEL(I),I=1,NBSEF),
     %               ' of TETRAHEDRON 1 2 3 4 of VOLUME ',VOLUME,' <=0!'
         ENDIF

      ELSE

C        HEXAEDRE
         IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,*) 'efvone: HEXAEDRE ',(NOSOEL(I),I=1,NBSEF),
     %                 ' de TETRAEDRE 1 2 3 4 de VOLUME ',VOLUME,' <=0!'
         ELSE
         WRITE(IMPRIM,*) 'efvone: HEXAHEDRON ',(NOSOEL(I),I=1,NBSEF),
     %               ' of TETRAHEDRON 1 2 3 4 of VOLUME ',VOLUME,' <=0!'
         ENDIF
      ENDIF

      DO I=1,NBSEF
         WRITE(IMPRIM,*) 'efvone: XYZ ',NOSOEL(I),' =',
     %                   (XYZSOM(K,NOSOEL(I)),K=1,3)
      ENDDO

      RETURN
      END
