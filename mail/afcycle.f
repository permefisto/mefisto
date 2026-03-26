      SUBROUTINE AFCYCLE( NBCIAS, N1CIAS, NSASFVP )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    AFFICHER LE TABLEAU NSAFVP DES NBCIAS CYCLES D'ARETE
C -----

C ENTREES:
C --------
C NBCIAS : NOMBRE DE CYCLES D'ARETES
C N1CIAS : NO DANS NSASFVP DE LA PREMIERE ARETE DES NBCIAS CYCLES
C NBASFVP: NOMBRE DE FACES DES ARETES SIMPLES DES NBCIAS CYCLES
C NSASFVP: 1: NUMERO DU SOMMET 1  DE L'ARETE SIMPLE
C          2: NUMERO DU SOMMET 2  DE L'ARETE SIMPLE
C          3: NUMERO DANS NSASFVP DE L'ARETE SIMPLE SUIVANTE
C             SI =0 FIN DU CYCLE NON TRAITEE COMME UNE BOUCLE
C             SI JAMAIS =0 C'EST UNE BOUCLE A ARRETER AVEC N1CIAS(NOC)
C          4: NUMERO DANS NFVPSI DE LA FACE DE CETTE ARETE SIMPLE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  Veulettes sur mer               Janvier 2020
C2345X7..............................................................012
      INTEGER  N1CIAS(*), NSASFVP(4,*)

      DO NOC=1,NBCIAS

C        AFFICHAGE DU CYCLE NOC
         N1AC = N1CIAS( NOC )
         NAC  = N1AC

 10      PRINT*,'ciasfaet: CYCLE',NOC,' ARETE',NAC,
     %          ' NS1=',NSASFVP(1,NAC),' NS2=',NSASFVP(2,NAC),
     %          ' Suivant',NSASFVP(3,NAC),' FACE=',NSASFVP(4,NAC)
C        L'ARETE SUIVANTE
         NAC = NSASFVP( 3, NAC )
         IF( NAC .GT. 0 ) THEN
C           L'ARETE SUIVANTE EXISTE
            IF( NAC .NE. N1AC ) THEN
C              L'ARETE SUIVANTE N'EST PAS LA PREMIERE
               GOTO 10
            ENDIF
         ENDIF

      ENDDO

      RETURN
      END
