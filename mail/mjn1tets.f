      SUBROUTINE MJN1TETS( MXTETR, NOTETR, MXSOMM,
     %                     N1TETS, NUDTETR, NBSOMM )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    MISE A JOUR DU TABLEAU N1TETS
C -----

C ENTREES:
C --------
C MXTETR : NOMBRE DE TETRAEDRES DECLARABLES DANS NOTETR
C NOTETR : LISTE DES TETRAEDRES
C          SOMMET1,    SOMMET2,    SOMMET3,    SOMMET4,
C          TETRAEDRE1, TETRAEDRE2, TETRAEDRE3, TETRAEDRE4
C          DE L'AUTRE COTE DE LA FACE
C          1: 123      2: 234      3: 341      4: 412
C MXSOMM : NOMBRE DE SOMMETS DECLARABLES DANS PTXYZD

C SORTIES:
C --------
C N1TETS : N1TETS(NS) NUMERO D'UN TETRAEDRE AYANT POUR SOMMET NS
C NUDTETR: NUMERO NOTETR DU DERNIER TETRAEDRE ACTIF
C NBSOMM : NUMERO MAXIMAL DES SOMMETS DE NOTETR
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St Pierre du Perray  Octobre 2017
C2345X7..............................................................012
      INTEGER  NOTETR(8,MXTETR), N1TETS(MXSOMM)

ccc      PRINT*,'mjn1tets: NBSOMM INITIAL =',NBSOMM

C     MISE A JOUR DU TABLEAU N1TETS
C     -----------------------------
      CALL AZEROI( MXSOMM, N1TETS )

      NBSOMM = 0
      DO NTE = 1, MXTETR
         IF( NOTETR(1,NTE) .GT. 0 ) THEN
            NUDTETR = NTE
            DO K=1,4
               NS = NOTETR(K,NTE)
               N1TETS( NS ) = NTE
               IF( NS .GT. NBSOMM ) NBSOMM = NS
            ENDDO
         ENDIF
      ENDDO

      PRINT*,'mjn1tets: Fin avec NBSOMM=',NBSOMM,' NUDTETR=',NUDTETR

      RETURN
      END
