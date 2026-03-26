      SUBROUTINE AJTRET( NOTRI1, NOTRIA, N1AEVI, N1AEOC, NAETOI )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    AJOUTER LES 3 ARETES DU TRIANGLE NOTRI1 DE NOTRIA
C -----    A L'ETOILE D'ARETES STOCKEES DANS NAETOI
C
C ENTREE :
C --------
C NOTRI1 : NUMERO DANS NOTRIA DU TRIANGLE D'ARETES A AJOUTER
C NOTRIA : LISTE DES TRIANGLES
C                 ------- ------- ------- -------- -------- --------
C  PAR TRIANGLE : SOMMET1 SOMMET2 SOMMET3 TR_VOIS1 TR_VOIS2 TR_VOIS3
C                 ------- ------- ------- -------- -------- --------
C                 SOMMET    EST LE NUMERO DU SOMMET
C                 TR_VOIS I EST LE NUMERO DANS NOTRIA DU TRIANGLE
C                               ADJACENT PAR L'ARETE I
C
C ENTREES ET SORTIES:
C -------------------
C N1AEVI : NUMERO DANS NAETOI DE LA PREMIERE ARETE VIDE
C N1AEOC : NUMERO DANS NAETOI DE LA PREMIERE ARETE OCCUPEE
C NAETOI : CHAINAGE SUR LE SUIVANT DANS NAETOI(4,*)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC    FEVRIER 1992
C2345X7..............................................................012
      INTEGER   NOTRIA(6,*),NAETOI(4,*)
C
C     RECHERCHE DES ARETES DE L'ETOILE PARMI LES 3 DU TRIANGLE NOTRI1
C     ---------------------------------------------------------------
      DO 30 I=1,3
C
C        SI    ( L'ARETE I DU TRIANGLE NOTRI1 N'APPARTIENT PAS
C                AUX ARETES DE L'ETOILE NAETOI )
C        ALORS ELLE EST AJOUTEE A L'ETOILE DANS NAETOI
C        SINON ELLE EST EMPILEE DANS NPILE POUR ETRE DETRUITE ENSUITE
C              ELLE EST SUPPRIMEE DE L'ETOILE NAETOI
C
C        LE NUMERO NTOP DU TRIANGLE DE L'AUTRE COTE DE L'ARETE I
         NTOP = NOTRIA(I+3,NOTRI1)
         IF( NTOP .LE. 0 ) THEN
C           ARETE TRAITEE
            NC = -I
            GOTO 20
         ENDIF
         NC = I
C
c        LE NAOP NUMERO DE L'ARETE CONCERNEE DANS LE TRIANGLE OPPOSE
         IF( NOTRIA(4,NTOP) .EQ. NOTRI1 ) THEN
            NAOP = 1
         ELSE IF( NOTRIA(5,NTOP) .EQ. NOTRI1 ) THEN
            NAOP = 2
         ELSE
            NAOP = 3
         ENDIF
C
C        LE DEBUT DU CHAINAGE DES ARETES DE L'ETOILE
C        ET BOUCLE SUR LES ARETES DE L'ETOILE
         NA1 = 0
         NA2 = N1AEOC
 10      IF( NA2 .GT. 0 ) THEN
            IF( NAETOI(1,NA2) .NE. NTOP .OR.
     %          ABS(NAETOI(2,NA2)) .NE. NAOP ) THEN
C              L'ARETE EST DIFFERENTE.PASSAGE A LA SUIVANTE
               NA1 = NA2
               NA2 = NAETOI(4,NA2)
               GOTO 10
            ELSE
C              L'ARETE NAOP DU TRIANGLE NTOP EST IDENTIQUE A NA2
C              DESTRUCTION DE L'ARETE DE L'ETOILE
               IF( NA1 .NE. 0 ) THEN
C                 L'ARETE NA2,PRECEDEE DE NA1,N'EST PAS LA
C                 PREMIERE DE L'ETOILE
                  NAETOI(4,NA1) = NAETOI(4,NA2)
C                 L'ARETE NA2 DEVIENT LA PREMIERE VIDE DE L'ETOILE
                  NAETOI(4,NA2) = N1AEVI
                  N1AEVI = NA2
               ELSE
C                 L'ARETE NA2=N1AEOC EST LA 1-ERE DE L'ETOILE
                  NA1    = NAETOI(4,N1AEOC)
                  NAETOI(4,N1AEOC) = N1AEVI
                  N1AEVI = N1AEOC
                  N1AEOC = NA1
               ENDIF
               GOTO 30
            ENDIF
         ENDIF
C
C        L'ARETE NAOP N'EST PAS RETROUVEE
C        ELLE EST AJOUTEE A L'ETOILE AU DEBUT DES ARETES OCCUPEES
 20      NA1           = N1AEVI
         N1AEVI        = NAETOI(4,N1AEVI)
         NAETOI(1,NA1) = NOTRI1
         NAETOI(2,NA1) = NC
         NAETOI(4,NA1) = N1AEOC
         N1AEOC        = NA1
 30   CONTINUE
      END
