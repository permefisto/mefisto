      SUBROUTINE ENTMUL( BASE , A , B , R )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  R EST LE PRODUIT DE 2 ENTIERS MULTI MOTS A ET B TELS QUE
C -----  A>=0 B>=0
C
C        ATTENTION : BASE DOIT VERIFIER
C        =========== BASE**2 < AU MAXIMUM DOUBLE PRECISION
C                   2*BASE-2 < AU MAXIMUM ENTIER SUR UN MOT
C                    MIN(NA,NB) + 1 < BASE
C
C ENTREES :
C ---------
C BASE    : LA BASE DE TRAITEMENT
C A       : ENTIER MULTI-MOTS  CODE SOUS LA FORME
C                I=NA
C           A = SOMME ( A(I) * BASE ** I )   ET  A(-1) = SIGNE( A )
C                I=0                             A(-2) = NA
C
C B       : ENTIER MULTI-MOTS  CODE SOUS LA FORME
C                I=NB
C           B = SOMME ( B(I) * BASE ** I )   ET  B(-1) = SIGNE( B )
C                I=0                             B(-2) = NB
C
C SORTIES :
C ---------
C NR      : NOMBRE DE MOTS - 1 DE L'ENTIER R
C R       : ENTIER MULTI-MOTS  A * B CODE SOUS LA FORME
C                I=NR
C           R = SOMME ( R(I) * BASE ** I )   ET  R(-1) = SIGNE( R )
C                I=0                             R(-2) = NR
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS   FEVRIER 1988
C23456...............................................................012
      INTEGER          A(-2:*),B(-2:*),R(-2:*)
      INTEGER          BASE,S0,S1,S2,R1
      DOUBLE PRECISION D0,DS0,D1
C
C     DETECTION D'UN ENTIER NUL
      NA = A(-2)
      NB = B(-2)
      IF( (NA.EQ.0 .AND. A(0).EQ.0) .OR.
     %    (NB.EQ.0 .AND. B(0).EQ.0) ) THEN
C
C        R = A * B = 0
C        -------------
         R(-2 ) = 0
         R(-1 ) = 1
         R( 0 ) = 0
         RETURN
      ENDIF
C
C     ==============================================================
C               K=NA+NB   I=MIN(NA,K)
C     R = A*B = SOMME ( ( SOMME       A(I) * B(K-I) ) * BASE ** K  )
C               K=0       I=MAX(0,K-NB)
C
C     LA SOMME INTERNE POUR K FIXE S'ECRIT
C     S = S0 + BASE * S1 + BASE**2 * S2
C
C     EN EFFET
C     S =< ( MIN(NA,K) - MAX(0,K-NB) ) * (BASE-1)**2
C       =< ( MIN(NA,NB) + 1 ) * (BASE-1)**2
C              ET POUR MIN(NA,NB)+1 < BASE
C        < BASE * (BASE-1)**2 = BASE**2 * (BASE-2) + BASE
C
C     ENTRAINE L'ABSENCE DE TERME SUR BASE**3
C     ==============================================================
C
      S1 = 0
      S2 = 0
      NR = NA + NB
C
      DO 20 K=0,NR
         S0 = S1
         S1 = S2
         S2 = 0
C
C        LE COEFFICIENT DE BASE**K
C        =========================
         DO 10 I = MAX(0,K-NB) , MIN(NA,K)
C
C           A(I) EST MIS DANS UN DOUBLE PRECISION POUR EFFECTUER
C           LES CALCULS EN REEL DOUBLE PRECISION
            D0  = A(I)
            DS0 = S0 + D0 * B(K-I)
C
C           DS0 =< (BASE-1) + (BASE-1)**2 = BASE * (BASE-1)
C                   S0 + R1 * BASE    AVEC  R1 =< BASE-1
C
            IF( DS0 .GE. BASE ) THEN
C
C              R1 EST LA PARTIE ENTIERE DE LA DIVISION
C              ET NON PAS LE PLUS PROCHE ENTIER
               R1 = INT( DS0 / BASE )
C              R1 EST MIS DANS UN DOUBLE PRECISION POUR EFFECTUER
C              LES CALCULS EN REEL DOUBLE PRECISION
               D1 = R1
C
C              ICI LE CALCUL DOIT ETRE EXACT => CONTRAINTE SUR BASE
C              ET S0 EST LE PLUS PROCHE ENTIER DE LA SOUSTRACTION
               S0 = NINT( DS0 - D1 * BASE )
C
C              LA RETENUE EST AJOUTEE A S1
               S1 = S1 + R1
C
C              S1 =< BASE-1 + BASE-1
C                 =< BASE   + BASE-2
C              ENTRAINE UNE RETENUE D'AU PLUS 1
C
               IF( S1 .GE. BASE ) THEN
                  S1 = S1 - BASE
                  S2 = S2 + 1
               ENDIF
C
            ELSE
C
C              DS0 =< BASE-1
               S0 = NINT( DS0 )
            ENDIF
 10      CONTINUE
C
         R(K) = S0
 20   CONTINUE
C
C     ENCADREMENT DE R = A * B
C     ========================
C 1*BASE**NA * 1*BASE**NB =< A*B =< (BASE**(NA+1)-1) * (BASE**(NB+1)-1)
C
C 1 * BASE ** ( NA + NB ) =< A*B =< (BASE**(NA+NB+2)-1)
C
C                                   I=NA+NB+1
C 1 * BASE ** ( NA + NB ) =< A*B =<  SOMME (BASE-1) * BASE**I
C                                   I=0
C     CET ENCADREMENT ENTRAINE QUE S2=0 EN FIN DE BOUCLE K
C
      IF( S1 .GT. 0 ) THEN
         NR = NR + 1
         R( NR ) = S1
      ENDIF
      R(-2) = NR
C
C     LE SIGNE DU PRODUIT R = A * B
      R(-1) = A(-1) * B(-1)
      IF( NR.EQ.0 .AND. R(0).EQ.0 ) THEN
         R(-1 ) = 1
      ENDIF
      END
