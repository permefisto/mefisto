      SUBROUTINE GETRI6( L1ARET, L2ARET,  LARETE,
     %                   NBTRIA, M1TRIA4, NOTRIA4,
     %                           M1TRIA6, NOTRIA6, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    GENERER SUR NOTRIA4(4,NBTRIA) LE TABLEAU NOTRIA6(6,NBTRIA)
C -----
C ATTENTION: NOTRIA4 ET NOTRIA6 PEUVENT AVOIR MEME ADRESSE A L'APPEL

C ENTREES:
C --------
C NBTRIA : NOMBRE INITIAL ET FINAL DE TRIANGLES
C M1TRIA4: MAXIMUM DU 1-ER INDICE DU TABLEAU NOTRIA4
C          ATTENTION: OBLIGATION de M1TRIA4<=M1TRIA6
C NOTRIA4: LISTE DES 3 SOMMETS ET 0 DES NBTRIA TRIANGLES
C                 ------- ------- ------- --------
C  PAR TRIANGLE : SOMMET1 SOMMET2 SOMMET3    0     ou AUTRES CAS
C                 ------- ------- ------- --------
C L1ARET : NOMBRE DE MOTS PAR ARETE DU TABLEAU NARET
C L2ARET : NOMBRE DE ARETES DU TABLEAU NARET
C LARETE : TABLEAU DES ARETES DU MAILLAGE
C          LARETE(1,I)= NO DU 1-ER  SOMMET DE L'ARETE
C          LARETE(2,I)= NO DU 2-EME SOMMET > 1-ER  SOMMET
C          LARETE(3,I)= CHAINAGE HACHAGE SUR ARETE SUIVANTE
C          LARETE(4,I)= NUMERO DU 1-ER TRIANGLE CONTENANT CETTE ARETE
C                       0 SI PAS DE 1-ER  TRIANGLE
C          LARETE(5,I)= NUMERO DU 2-EME TRIANGLE CONTENANT CETTE ARETE
C                       0 SI PAS DE 2-EME TRIANGLE
C          SI IL EXISTE PLUS DE 2 TRIANGLES PAR ARETE,
C             LES 3,... SONT IGNORES

C SORTIES :
C ---------
C M1TRIA6 : MAXIMUM DU 1-ER INDICE DU TABLEAU NOTRIA6
C NOTRIA6 : LISTE CHAINEE DES TRIANGLES  (GENERALEMENT SUR NOTRIA4 !)
C                 ------- ------- ------- -------- -------- --------
C  PAR TRIANGLE : SOMMET1 SOMMET2 SOMMET3 TR_VOIS1 TR_VOIS2 TR_VOIS3
C                 ------- ------- ------- -------- -------- --------
C                 SOMMET    EST LE NUMERO DU SOMMET
C                 TR_VOIS i EST LE NUMERO DANS NOTRIA6 DU TRIANGLE
C                          ADJACENT PAR L'ARETE i
C IERR   : 0 SI PAS D'ERREUR, >0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC    OCTOBRE 1993
C....................................................................012
      INTEGER        LARETE(L1ARET,L2ARET),
     %               NOTRIA4(M1TRIA4,NBTRIA),
     %               NOTRIA6(M1TRIA6,NBTRIA)
      INTEGER        NGS(2)

C     REPARTITION DE NOTRIA4 DANS NOTRIA6
C     TRANSFERT DU NO DES 3 SOMMETS ET MISE A ZERO DU NO DES TRIANGLES VOISINS
      DO J = NBTRIA, 1, -1
         DO I=3,1,-1
C           LE SOMMET I (EN ARRIERE EST IMPORTANT. SINON BOGUE!)
            NOTRIA6( I, J ) = NOTRIA4( I, J )
C           LE TRIANGLE VOISIN PAR L'ARETE I
            NOTRIA6( 3+I,J ) = 0
         ENDDO
      ENDDO

C     LE CHAINAGE SUR LES TRIANGLES VOISINS
      LIBREF = L2ARET
      DO J = 1, NBTRIA

         NS1 = NOTRIA6(1,J)
         DO 40 I=1,3

C           ARETE NS1-NS2 DU TRIANGLE J
            IF( I .NE. 3 ) THEN
               NS2 = I + 1
            ELSE
               NS2 = 1
            ENDIF
            NS2 = NOTRIA6(NS2,J)
            IF( NS1 .LT. NS2 ) THEN
               NGS(1) = NS1
               NGS(2) = NS2
            ELSE
               NGS(1) = NS2
               NGS(2) = NS1
            ENDIF

C           RECHERCHE DE L'ARETE NS1-NS2 DANS LES ARETES
C           DE LA TRIANGULATION
            CALL HACHAG( 2, NGS, L1ARET, L2ARET, LARETE, 3,
     &                   LIBREF, NOAR )

            IF( NOAR .LE. 0 ) THEN
C              PROBLEME: ARETE I DU TRIANGLE J NON RETROUVEE
               PRINT*,'getri6: Probleme ARETE de SOMMETS',NS1,NS2,
     %                ' NON RETROUVEE dans la LISTE des ARETES'
               IERR = 1
               GOTO 40
            ENDIF

C           LES 2 TRIANGLES ADJACENTS PAR L'ARETE NS1-NS2
            NT1 = ABS( LARETE(4,NOAR) )
            NT2 = ABS( LARETE(5,NOAR) )
            IF( NT1 .EQ. J ) THEN
               NOTRIA6(3+I,J) = NT2
            ELSE IF( NT2 .EQ. J ) THEN
               NOTRIA6(3+I,J) = NT1
            ENDIF

C           ARETE SUIVANTE
            NS1 = NS2

 40      ENDDO

      ENDDO

      RETURN
      END
