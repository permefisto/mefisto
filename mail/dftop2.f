      SUBROUTINE DFTOP2( MOELEM, NBELEM, MNELEM,
     %                   MOELTG, NBELTG, MNELTG )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    DECLARER EN MCN LES TABLEAUX DE HACHAGES DE L'OBJET
C -----    ET LES TABLEAUX DES NUMEROS DES TG DES EF A TG

C ENTREES:
C --------
C MOELEM : NOMBRE DE MOTS DE CHAQUE TYPE D'ELEMENT FINI
C NBELEM : NBELEM(1)=NOMBRE DE NOEUDSOMMETS
C          NBELEM(2)=SEGMENTS
C          NBELEM(3)=TRIANGLES
C          NBELEM(4)=QUADRANGLES
C          NBELEM(5)=TETRAEDRES
C          NBELEM(6)=PENTAEDRES
C          NBELEM(7)=HEXAEDRES
C          NBELEM(8)=6-CUBES
C MOELTG : NOMBRE DE     TG PAR TYPE D'EF
C NBELTG : NOMBRE D'EF A TG PAR TYPE D'EF
C
C SORTIES:
C --------
C MNELEM : ADRESSE MCN DES EVENTUELS 8 TABLEAUX DES ELEMENTS FINIS
C          0 SI PAS DE D'ELEMENTS FINIS DE CE TYPE
C MNELTG : ADRESSE MCN DES EVENTUELS 8 TABLEAUX DES ELEMENTS FINIS A TG
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE UPMC PARIS   DECEMBRE 1996
C MODIF  : ALAIN PERRONNET  TEXAS A & M UNIVERSITY            JULY  2005
C2345X7..............................................................012
      include"./incl/langue.inc"
      include"./incl/a___npef.inc"
      include"./incl/pp.inc"
      COMMON       MCN(MOTMCN)
      INTEGER      MOELEM(9), NBELEM(9), MNELEM(9)
      INTEGER      MOELTG(9), NBELTG(9), MNELTG(9)

      DO I=1,9
         IF( NBELEM(I) .LE. 0 ) THEN
            MNELEM(I) = 0
         ELSE
            NBV = MOELEM(I) * NBELEM(I)
            PRINT*,'dftop2: TYPE d''EF',I,
     %             ' reservation de',NBV,' ENTIERS pour le HACHAGE'
            CALL TNMCDC( 'ENTIER', NBV, MNELEM(I) )
            CALL AZEROI(  NBV, MCN(MNELEM(I)) )
         ENDIF

         IF( NBELTG(I) .LE. 0 ) THEN
            MNELTG(I) = 0
         ELSE
            NBV = MOELTG(I) * NBELTG(I)
            CALL TNMCDC( 'ENTIER', NBV, MNELTG(I) )
            CALL AZEROI(  NBV, MCN(MNELTG(I)) )
         ENDIF
      ENDDO

      RETURN
      END
