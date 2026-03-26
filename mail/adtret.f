      SUBROUTINE ADTRET( NT,     MOSOAR, NOSOAR, MOARTR, NOARTR,
     %                   N1AEVI, N1AEOC, NAETOI )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    AJOUTER LES 3 ARETES DU TRIANGLE NT DE NOARTR
C -----    A L'ETOILE D'ARETES STOCKEES DANS LE TABLEAU NAETOI
C
C ENTREE :
C --------
C NT     : NUMERO DANS NOARTR DU TRIANGLE D'ARETES A AJOUTER
C MOSOAR : NOMBRE MAXIMAL D'ENTIERS PAR ARETE ET
C          INDICE DANS NOSOAR DE L'ARETE SUIVANTE DANS LE HACHAGE
C NOSOAR : NUMERO DES 2 SOMMETS , NO LIGNE, 2 TRIANGLES DE L'ARETE,
C          CHAINAGE DES ARETES FRONTALIERES, CHAINAGE DU HACHAGE DES ARETES
C MOARTR : NOMBRE MAXIMAL D'ENTIERS PAR ARETE DU TABLEAU NOARTR
C NOARTR : LES 3 ARETES DES TRIANGLES +-ARETE1, +-ARETE2, +-ARETE3
C
C ENTREES ET SORTIES:
C -------------------
C N1AEVI : NUMERO DANS NAETOI DE LA PREMIERE ARETE VIDE
C N1AEOC : NUMERO DANS NAETOI DE LA PREMIERE ARETE OCCUPEE
C NAETOI : CHAINAGE SUR LE SUIVANT DANS NAETOI(4,*)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC       MARS 1997
C2345X7..............................................................012
      INTEGER           NOARTR(MOARTR,*),
     %                  NOSOAR(MOSOAR,*),
     %                  NAETOI(4,*)
C
C     RECHERCHE DES ARETES DE L'ETOILE PARMI LES 3 DU TRIANGLE NT
C     -----------------------------------------------------------
      DO 30 I=1,3
C
C        SI    ( L'ARETE I DU TRIANGLE NT N'APPARTIENT PAS
C                AUX ARETES DE L'ETOILE NAETOI )
C        ALORS ELLE EST AJOUTEE A L'ETOILE DANS NAETOI
C        SINON ELLE EST EMPILEE DANS NPILE POUR ETRE DETRUITE ENSUITE
C              ELLE EST SUPPRIMEE DE L'ETOILE NAETOI
C
C        LE NUMERO DE L'ARETE NOAR DANS NOSOAR DE L'ARETE I DU TRIANGLE NT
         NOAR = ABS( NOARTR(I,NT) )
C
C        LE NUMERO NTOP DU TRIANGLE DE L'AUTRE COTE DE L'ARETE I
         IF( NOSOAR(4,NOAR) .EQ. NT ) THEN
            NTOP = NOSOAR(5,NOAR)
         ELSE
            NTOP = NOSOAR(4,NOAR)
         ENDIF
         IF( NTOP .LE. 0 ) THEN
C           ARETE TRAITEE
            NC = -I
            GOTO 20
         ENDIF
         NC = I
C
C        LE NUMERO NAOP DE L'ARETE CONCERNEE DANS LE TRIANGLE OPPOSE
         IF( ABS( NOARTR(1,NTOP) ) .EQ. NOAR ) THEN
            NAOP = 1
         ELSE IF( ABS( NOARTR(2,NTOP) ) .EQ. NOAR ) THEN
            NAOP = 2
         ELSE
            NAOP = 3
         ENDIF
C
C        LE DEBUT DU CHAINAGE DES ARETES DE L'ETOILE
C        ET BOUCLE SUR LES ARETES DE L'ETOILE
C        ===========================================
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
C
C              L'ARETE NAOP DU TRIANGLE NTOP EST IDENTIQUE A NA2
C              DESTRUCTION DE L'ARETE NA2 DE L'ETOILE NAETOI
C              -------------------------------------------------
               IF( NA1 .NE. 0 ) THEN
C                 L'ARETE NA2, PRECEDEE DE NA1, N'EST PAS LA
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
C        --------------------------------------------------------
C        LE NUMERO DE LA PREMIERE ARETE VIDE DANS NAETOI
 20      NA1           = N1AEVI
C        LE NUMERO DE LA NOUVELLE PREMIERE ARETE VIDE DANS NAETOI
         N1AEVI        = NAETOI(4,N1AEVI)
C        LE NUMERO NOARTR DU TRIANGLE
         NAETOI(1,NA1) = NT
C        LE NUMERO DE L'ARETE DANS NT
         NAETOI(2,NA1) = NC
C        LE CHAINAGE SUR L'ARETE OCCUPEE SUIVANTE
         NAETOI(4,NA1) = N1AEOC
         N1AEOC        = NA1
C
 30   CONTINUE
      END
