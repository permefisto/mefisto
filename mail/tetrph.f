      SUBROUTINE TETRPH( NBEFI,  NUSSOS, NBDM, NUMAVI,
     %                   NUMAVF, NBTETR, NBPYRA, NUSTET )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    METTRE EN NUMERO 1 LE PLUS PETIT NUMERO DE SOMMET DE L'EF
C -----
C ENTREES:
C --------
C NBEFI  : NOMBRE DES EF INITIAUX
C NUSSOS : NO DES SOMMETS DES EF INITIAUX
C NBDM   : NOMBRE DE MATERIAUX
C NUMAVI : NUMERO DE MATERIAU DES EF INITIAUX
C
C SORTIES:
C --------
C NUMAVF : NUMERO DE MATERIAU DES EF FINAUX
C NBTETR : NOMBRE DE TETRAEDRES
C NBPYRA : NOMBRE DE PYRAMIDES
C NUSTET : NO DES SOMMETS DES TETRAEDRES FINAUX
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEURS: F.CHOUKROUN O.RICOU UPMC ANALYSE NUMERIQUE PARIS JANVIER 1990
C MODIFS : ALAIN PERRONNET Laboratoire J-L. LIONS UPMC PARIS   Mars 2007
C MODIFS : ALAIN PERRONNET LJLL UPMC et St Pierre du Perray    Juin 2008
C2345X7..............................................................012
      INTEGER NUSSOS(8,*), NUSTET(5,*), TE1SOB(4,6), NSEFJ(8), PLUTET
      INTEGER NUMAVI(NBEFI), NUMAVF(*)
C
      NOMATE = 0
      NBTETR = 0
      NBPYRA = 0
      DO 10 NEF=1,NBEFI
C
C        NO DU MATERIAU DE L'EF INITIAL
         IF( NBDM .GT. 0 ) NOMATE = NUMAVI(NEF)
C
         IF(NUSSOS(5,NEF).EQ.0) THEN
C
C  ------------------    TETRAEDRE ---------
            NCOGEL = 5
            NBTETR = NBTETR+1
            DO 11 J=1,4
               NUSTET(J,NBTETR)=NUSSOS(J,NEF)
 11         CONTINUE
C           NO DU MATERIAU DU TETRAEDRE FINAL
            IF( NBDM .GT. 0 ) NUMAVF(NBTETR) = NOMATE
C
         ELSE IF(NUSSOS(6,NEF).EQ.0) THEN
C
C  ------------------    PYRAMIDE ---------
            NCOGEL = 9
            NBTETR = NBTETR + 1
            NBPYRA = NBPYRA + 1
            DO 12 J=1,5
               NUSTET(J,NBTETR)=NUSSOS(J,NEF)
 12         CONTINUE
C           NO DU MATERIAU DU TETRAEDRE FINAL
            IF( NBDM .GT. 0 ) NUMAVF(NBTETR) = NOMATE
C
         ELSE IF(NUSSOS(7,NEF).EQ.0) THEN
C
C  ------------------    PENTAEDRE ---------
            NCOGEL = 6
            DO 30 I=1,6
               NSEFJ(I)=NUSSOS(I,NEF)
 30         CONTINUE
            CALL TETSPH( NCOGEL, NSEFJ, TE1SOB, PLUTET )
            DO 21 NUTE=1,3
               DO 22 J=1,4
                  NUSTET(J,NBTETR+NUTE)=TE1SOB(J,NUTE)
 22            CONTINUE
 21         CONTINUE
C           NO DU MATERIAU DES 3 TETRAEDRES FINAUX
            IF( NBDM .GT. 0 ) THEN
               NUMAVF(NBTETR+1) = NOMATE
               NUMAVF(NBTETR+2) = NOMATE
               NUMAVF(NBTETR+3) = NOMATE
            ENDIF
            NBTETR=NBTETR+3
C
         ELSE
C
C  ------------------    HEXAEDRE ---------
            NCOGEL = 7
            DO 40 I=1,8
               NSEFJ(I)=NUSSOS(I,NEF)
 40         CONTINUE
            CALL TETSPH( NCOGEL, NSEFJ, TE1SOB, PLUTET )
            DO 31 NUTE=1,PLUTET
               DO 32 J=1,4
                  NUSTET(J,NBTETR+NUTE)=TE1SOB(J,NUTE)
 32            CONTINUE
C              NO DU MATERIAU DU TETRAEDRE FINAL
               IF( NBDM .GT. 0 ) NUMAVF(NBTETR+NUTE) = NOMATE
 31         CONTINUE
            NBTETR=NBTETR+PLUTET
         ENDIF
 10   CONTINUE

      RETURN
      END
