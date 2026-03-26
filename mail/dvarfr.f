      SUBROUTINE DVARFR( NS1,NS2,NLSOFR,NBLFTR,NDARLF,NOSOAR,
     %                   NONOUI )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    ETABLIR SI L'ARETE NS1 NS2 EST UNE ARETE D'UNE LIGNE FERMEE
C -----
C
C ENTREES:
C --------
C NS1 NS2: LES 2 NUMEROS DE SOMMETS DE L'ARETE
C NLSOFR : TABLEAU DU NUMERO DE LIGNE(1 A NBLFTR) DES SOMMETS
C         -NUMERO DE POINT INTERNE UTILISATEUR IMPOSE, 0 SINON
C NBLFTR : NOMBRE DE LIGNES FERMEES LIMITANT LA SURFACE A TRIANGULER
C NDARLF : NUMERO DE LA PREMIERE ARETE DE CHAQUE LIGNE FERMEE DANS
C          LE TABLEAU NOSOAR
C NOSOAR : NUMERO DES 2 SOMMETS DE CHAQUE ARETE ET
C          POINTEUR SUR L'ARETE SUIVANTE
C
C SORTIES:
C --------
C NONOUI  : 1 SI NS1 NS2 A ETE RETROUVEE DANS LA LIGNE FERMEE L
C           0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC    JANVIER 1991
C....................................................................012
      INTEGER           NDARLF(NBLFTR),
     %                  NLSOFR(*),
     %                  NOSOAR(3,*)
C
C     LE NUMERO EVENTUEL DE LA LIGNE FERMEE DU SOMMET 1
      L = NLSOFR( NS1 )
      IF( L .GT. 0 .AND. L .EQ. NLSOFR(NS2) ) THEN
C
C        ---------------------------------------------------
C        LES 2 SOMMETS SONT SUR LA MEME LIGNE FERMEE L
C        CETTE ARETE EST ELLE FRONTALIERE DE LIGNE FERMEE L?
C        ---------------------------------------------------
         NA3 = NDARLF(L)
C
C        BOUCLE SUR LES ARETES DE LA LIGNE FERMEE L
C        ==========================================
 10      IF( NA3 .GT. 0 ) THEN
C           L'ARETE DE LA LIGNE FERMEE L
            NS3 = NOSOAR(1,NA3)
            NS4 = NOSOAR(2,NA3)
            IF(( NS3 .NE. NS1 .AND. NS3 .NE. NS2 ).OR.
     %         ( NS4 .NE. NS1 .AND. NS4 .NE. NS2 )) THEN
C
C              ARETE NON RETROUVEE PASSAGE A LA SUIVANTE
               NA3 = NOSOAR(3,NA3)
               GOTO 10
C
            ELSE
C
C              ARETE RETROUVEE SUR LA FRONTIERE L
C              LE POINT N EST DETRUIT AFIN DE CONSERVER
C              CETTE ARETE DANS LA TRIANGULATION ACTUELLE
               NONOUI = 1
               RETURN
            ENDIF
         ENDIF
      ENDIF
C
C     ARETE NON RETROUVEE
      NONOUI = 0
      END
