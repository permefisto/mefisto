      SUBROUTINE COEXPD( LADERO , NOINST , NCINST )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  COPIE DE L'ARBRE D'UNE EXPRESSION ARITHMETIQUE
C -----  NCOPER(1:4,1:LADERO) DANS LE TABLEAU NCINST D'UNE FONCTION
C        A PARTIR DE L'INSTRUCTION NOINST
C
C ENTREES :
C ---------
C LADERO : LA DERNIERE OPERATION DE NCOPER A COPIER
C
C ENTREE ET SORTIE :
C ------------------
C NOINST : LA DERNIERE INSTRUCTION ENREGISTREE ACTUELLEMENT
C NCINST : LE TABLEAU DES INSTRUCTIONS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS    FEVRIER 1990
C23456---------------------------------------------------------------012
      include"./incl/lu.inc"
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      INTEGER           NCINST(4,*)
C
      DO 10 J=1,LADERO
C
C        LE CODE OPERATION
         NCOP = NCOPER(1,J)
C        L'INSTRUCTION DANS NCINST
         K    = J + NOINST
C
C        LE CODE OPERATION
         NCINST(1,K) = NCOPER(1,J)
C
C        L'INSTRUCTION PERE
         IF( NCOPER(2,J) .GT. 0 ) THEN
            NCINST(2,K) = NCOPER(2,J) + NOINST
         ELSE
            NCINST(2,K) = NCOPER(2,J)
         ENDIF
C
C        LE FILS GAUCHE
         IF( ( NCOP .GE. 3  .AND. NCOP .LE. 5 ) .OR.
     %       ( NCOP .GE. 10 .AND. NCOP .LE. 12) .OR.
     %         NCOP .EQ. 14 ) THEN
            NCINST(3,K) = NCOPER(3,J)
         ELSE
            IF( NCOPER(3,J) .GT. 0 ) THEN
               NCINST(3,K) = NCOPER(3,J) + NOINST
            ELSE
               NCINST(3,K) = NCOPER(3,J)
            ENDIF
         ENDIF
C
C        LE FILS DROIT
         IF( NCOPER(4,J) .GT. 0 ) THEN
            NCINST(4,K) = NCOPER(4,J) + NOINST
         ELSE
            NCINST(4,K) = NCOPER(4,J)
         ENDIF
 10   CONTINUE
C
C     MISE A JOUR DE LA DERNIERE INSTRUCTION ENREGISTREE ACTUELLEMENT
      NOINST = NOINST + LADERO
      END
