      SUBROUTINE NOTRI1ST( NSt,    M1TRIA, NBSTRI, NBTRIA, NOTRIA,
     %                     MXTRST, NBTRST, NOTRST )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    RETROUVER 'EN FORCE' TOUS LES TRIANGLES DE SOMMET NSt A
C -----    PARTIR DU TABLEAU NOTRIA

C ENTREES:
C --------
C M1TRIA : NOMBRE MAXIMAL DU PREMIER INDICE DU TABLEAU NOTRIA (=4 ou =6)
C NBSTRI : NUMBRE DE SOMMETS DE L'EF A RETROUVER (TRIANGLE=3,QUADRANGLE=4)
C NBTRIA : NOMBRE DE TRIANGLES
C MXTRST : NOMBRE MAXIMAL DE TRIANGLES DE SOMMET NSt
C NOTRIA : NUMERO XYZSOM DES 3 SOMMETS et NUMERO NOTRIA DES 3 TRIANGLES
C          ADJACENTS

C SORTIES:
C --------
C NBTRST : NOMBRE DE TRIANGLES de SOMMET NSt
C NOTRST : NUMERO NOTRIA DES NBTRST TRIANGLES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET  VEULETTES SUR MER                Fevrier 2020
C....................................................................012
      INTEGER  NOTRIA(M1TRIA,NBTRIA), NOTRST(MXTRST)

      NBTRST = 0
      DO 10 NTR=1,NBTRIA

         IF( NOTRIA(1,NTR) .GT. 0 ) THEN
C           LE TRIANGLE EST ACTIF
            DO K=1,NBSTRI
               IF( NOTRIA(K,NTR) .EQ. NSt ) THEN

C                 LE K-EME SOMMET DE NTR EST LE SOMMET NSt
                  IF( NBTRST .GE. MXTRST ) THEN
                     PRINT*,'notri1st: AUGMENTER MXTRST=',MXTRST
                     GOTO 9999
                  ENDIF

C                 NS APPARTIENT A NBTRST TRIANGLES
                  NBTRST = NBTRST + 1
                  NOTRST( NBTRST ) = NTR
                  GOTO 10

               ENDIF
            ENDDO
         ENDIF

 10   ENDDO

 9999 RETURN
      END
