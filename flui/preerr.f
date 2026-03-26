      SUBROUTINE PREERR( NBNOPR, NBVECT,  PREEXA, PRECAL,
     %                   NUTYEL, NBELEM,  NONOEF,
     %                   ERRPRE, ERRPMIN, NOEMIN, NCAMIN,
     %                           ERRPMAX, NOEMAX, NCAMAX )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LE TABLEAU DE L'ERREUR SUR LES PRESSIONS
C -----    AUX NOEUDS DU MAILLAGE C-A-D
C          AUX SOMMETS POUR  BREZZI-FORTIN et TAYLOR-HOOD +
C          AUX MILIEUX DES ARETES POUR TAYLOR-HOOD
C             (SOMME ERREURS AUX 2 SOMMETS)/2
C
C ENTREES:
C --------
C NBNOPR : NOMBRE DE NOEUDS SUPPORT DE LA PRESSION
C NBVECT : NOMBRE TOTAL DE VECTEURS PRESSION PRESSION
C PREEXA : LA PRESSION EXACTE   EN CHAQUE NOEUD TEMPS
C PRECAL : LA PRESSION CALCULEE EN CHAQUE NOEUD TEMPS
C NUTYEL : NUMERO DU TYPE D'EF
C NBELEM : NOMBRE D'ELEMENTS FINIS
C NONOEF : NUMERO DES NOEUDS P2 DES NBELEM EF
C
C SORTIES:
C --------
C ERRPRE : ERREUR (PRESSIONExacte -PressionCalculee)(Temps,Noeuds)
C ERRPMIN: ERREUR MINIMALE sur la PRESSION(Temps,Noeuds)
C NOEMIN : NUMERO DU NOEUD   DE L'ERREUR MINIMALE sur la PRESSION(Temps,Noeuds)
C NCAMIN : NUMERO DU VECTEUR DE L'ERREUR MINIMALE sur la PRESSION(Temps,Noeuds)
C
C ERRPMAX: ERREUR MAXIMALE sur la PRESSION(Temps,Noeuds)
C NOEMAX : NUMERO DU NOEUD   DE L'ERREUR MAXIMALE sur la PRESSION(Temps,Noeuds)
C NCAMAX : NUMERO DU VECTEUR DE L'ERREUR MAXIMALE sur la PRESSION(Temps,Noeuds)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET TEXAS A & M UNIVERSITY at QATAR  FEVRIER 2012
C234567...............................................................12
      DOUBLE PRECISION  PREEXA(NBNOPR,NBVECT),
     %                  PRECAL(NBNOPR,NBVECT),
     %                  ERRPRE(NBNOPR,NBVECT)
      DOUBLE PRECISION  ERRPR, ERRPMIN, ERRPMAX
      INTEGER           NONOEF(NBELEM,*)
      INTRINSIC         SQRT
C
C     INITIALISATION DU MIN MAX
      NOEMIN = 1
      NCAMIN = 1
      ERRPMIN= 1D100
C
      NOEMAX = 1
      NCAMAX = 1
      ERRPMAX=-1D100
C
      DO K=1,NBVECT
C
C        LE TEMPS K
C        TEMPS = TIMES(K)
C
         DO I=1,NBNOPR
C
C           L'ERREUR |PRESSION|Exact - |PRESSION|Calcul (Noeud I)
            ERRPR = PREEXA(I,K) - PRECAL(I,K)
            ERRPRE(I,K) = ERRPR
C
C           MIN ET MAX AU NOEUD
            IF( ERRPR .LT. ERRPMIN ) THEN
               ERRPMIN = ERRPR
               NOEMIN = I
               NCAMIN = K
            ELSE IF( ERRPR .GT. ERRPMAX ) THEN
               ERRPMAX = ERRPR
               NOEMAX = I
               NCAMAX = K
            ENDIF
C
         ENDDO
C
      ENDDO
C
C     CORRECTION: ERREUR AU MILIEU D'UNE ARETE=(ERREURS AUX 2 EXTREMITES)/2
C     SINON LA VALEUR AU MILIEU EST GRANDE ET PRODUIT DES PICS
      IF( NUTYEL .EQ. 15 ) THEN
C
C        TRIANGLE P2 P1 de TAYLOR-HOOD
C        -----------------------------
         DO NEF=1,NBELEM
C
C           BOUCLE SUR LES 3 ARETES DU TRIANGLE NEF
            DO I=1,3
C
C              NO DU NOEUD MILIEU DE L'ARETE I
               NM = NONOEF(NEF,I+3)
C
C              NO DES 2 SOMMETS DE L'ARETE I
               N1 = NONOEF(NEF,I)
               IF( I .LT. 3 ) THEN
                  J = I+1
               ELSE
                  J = 1
               ENDIF
               N2 = NONOEF(NEF,J)
C
               DO K=1,NBVECT
                  ERRPRE(NM,K)=(ERRPRE(N1,K) + ERRPRE(N2,K)) * 0.5D0
               ENDDO
C
            ENDDO
         ENDDO
C
      ELSE IF( NUTYEL .EQ. 20 ) THEN
C
C        TETRAEDRE P2 P1 de TAYLOR-HOOD
C        ------------------------------
         DO NEF=1,NBELEM
C
C           NO DU SOMMET 4 DU TETRAEDRE
            N4 = NONOEF(NEF,4)
C
C           BOUCLE SUR LES 6 ARETES DU TETRAEDRE NEF
            DO I=1,3
C
C              NO DU NOEUD MILIEU DE L'ARETE I DU BAS
               NMB = NONOEF(NEF,I+4)
C
C              NO DU NOEUD MILIEU DE L'ARETE I+3 DU HAUT
               NMH = NONOEF(NEF,I+7)
C
C              NO DES 2 SOMMETS DE L'ARETE I DU BAS
               N1 = NONOEF(NEF,I)
               IF( I .LT. 3 ) THEN
                  J = I+1
               ELSE
                  J = 1
               ENDIF
               N2 = NONOEF(NEF,J)
C
               DO K=1,NBVECT
                  ERRPRE(NMB,K)=(ERRPRE(N1,K) + ERRPRE(N2,K)) * 0.5D0
                  ERRPRE(NMH,K)=(ERRPRE(N1,K) + ERRPRE(N4,K)) * 0.5D0
               ENDDO
C
            ENDDO
         ENDDO
C
      ENDIF
C
      RETURN
      END
