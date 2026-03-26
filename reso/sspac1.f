      SUBROUTINE SSPAC1( NOSTRA, MOMANA, MANAG,
     %                   NC, N,  INITX0, IFACTO, TETATO,
     %                   MOTSRG, RG,
     %                   NCODSA, MUA, A,
     %                   NCODSB, MUB, B, COMXMG,
     %                   NCVALP0, NBVALP0, VECP0,
     %                   R, D, VECJAC, AR, BR, W, V, MUTA, DAUX,
     %                   NBVALC, VALP, VECP, VALPES, IERR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   CALCUL DES VALEURS ET VECTEURS PROPRES VALP ET VECP DU PB :
C -----   ( A - VALP * B ) * VECP = 0
C         PAR LA METHODE DE L ITERATION SIMULTANEE INVERSE (SOUS-ESPACE)
C             LA B-ORTHONORMALISATION
C             L ANALYSE DE RAYLEIGH-RITZ SUR L ESPACE ITERE
C
C         CF RAPPORT SUR LE CALCUL DES VALEURS ET VECTEURS PROPRES
C
C ENTREES:
C --------
C NOSTRA : NOMBRE D ITERATIONS DEMANDEES . 0 => STRATEGIE AUTOMATIQUE
C MOMANA : NOMBRE DE MOTS DU TABLEAU D'ENTIERS MANAG
C MANAG  : TABLEAU DESCRIPTIF DE CHAQUE ITERATION DEMANDEE
C
C NCVALP0: NOMBRE DE VALEURS ET VECTEURS PROPRES DECLARES DEJA CALCULES
C NBVALP0: NOMBRE DE VALEURS ET VECTEURS PROPRES DEJA CALCULES
C NC     : NOMBRE DE VALEURS ET VECTEURS PROPRES A CALCULER CETTE FOIS
C N      : NOMBRE DE LIGNES ET COLONNES DES MATRICES A ET B
C          NOMBRE DES COMPOSANTES DE CHACUN DES VECTEURS PROPRES
C INITX0 : 0 X0 ESPACE DE DEPART DEJA DANS VECP
C          1 X0 METHODE BATHE
C          2,3,4 X0 GENERE ALEATOIREMENT
C          5 DIAGNOSTIC ET ARRET
C IFACTO : 0 A N EST PAS FACTORISEE,1 SINON
C TETATO : FACTEUR DE TRANSLATION ( A => A - TETATO * B )
C
C MOTSRG : LE NOMBRE DE MOTS DU TABLEAU RG
C RG     : TABLEAU DE LA MATRICE RG NON FACTORISEE
C
C NCODSA : CODE DE STOCKAGE DE LA MATRICE PROFIL A
C MUA    : POINTEUR SUR LES COEFFICIENTS DIAGONAUX DE A
C A      : MATRICE PROFIL REELLE DOUBLE PRECISION
C NCODSB : CODE DE STOCKAGE DE LA MATRICE PROFIL B
C MUB    : POINTEUR SUR LES COEFFICIENTS DIAGONAUX DE B
C B      : MATRICE PROFIL REELLE DOUBLE PRECISION
C COMXMG : MAXIMUM ACTUEL DES COEFFICIENTS DE LA MATRICE B
C
C VECP0  : VECTEUR(NCVALP0,N) DES VECTEURS PROPRES DEJA CALCULES EN CAS DE
C          VALEURS PROPRES OUBLIEES . VECP SINON .
C R,D,VECJAC,G,AR,BR,W,V,MUTA,DAUX : TABLEAUX AUXILIAIRES
C
C MODIFIES:
C --------
C NBVALC : NOMBRE DE VALEURS PROPRES CONVERGEES
C VALP   : TABLEAUX DES NC VALEURS  PROPRES CONVERGEES
C VECP   : TABLEAU  DES NC VECTEURS PROPRES CALCULES  TABLEAU (NC,N)!
C VALPES : ESTIMATION DE LA VALEUR PROPRE SUIVANT LA DERNIERE CONVERGEE
C IERR   : 0 SI PAS D'ERREUR RENCONTREE, >0 SINON
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEURS: BENAZETH, GOURDIN-SERVENIERE, PERRONNET  ANUM  SEPTEMBRE 1976
C MODIFS : ALAIN PERRONNET LABO ANALYSE NUMERIQUE UPMC PARIS   AOUT 1998
C MODIFS : ALAIN PERRONNET TEXAS A & M UNIVERSITY           JUILLET 2003
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/epsvvp.inc"
      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)
      DOUBLE PRECISION  U,T,RT,S,TETATO,DD
      DOUBLE PRECISION  COMXMG,A(*),B(*),RG(*),
     +                  VECP(NC,N),VALP(NC),R(NC,N),D(NC),
     +                  VECJAC(NC,NC),AR(NC,NC),BR(NC,NC),W(N),V(N),
     +                  DAUX(NC),VECP0(NCVALP0,N)
      INTEGER           MANAG(MOMANA),MUTA(NC),MUA(*),MUB(*),ISIGNE
      DOUBLE PRECISION  TRANSL, VALPES
C
 100  FORMAT(/'Debut SSPAC1: APPROXIMATION de X0 avec INITX0=',I6,
     % ' ',50('-'))
2100  FORMAT(/'Begin SSPAC1: APPROXIMATION of X0 with INITX0=',I6,
     % ' ',50('-'))
 110  FORMAT('ERREUR SSPAC1: 5 METHODES d''INITIALISATION de X0',
     +       ' INSUFFISANTES')
2110  FORMAT('ERROR SSPAC1 : 5 METHODS of INITIALISATION X0',
     +       ' INSUFFICIENT ')
 120  FORMAT('NOMBRE D''ITERATIONS MAXIMUM (NITEM) : ',I6/
     +     5(' MANAG(',I4,') : ',I1,2X))
2120  FORMAT('ITERATION NUMBER MAXIMUM (NITEM) : ',I6/
     +     5(' MANAG(',I4,') : ',I1,2X))
 130  FORMAT('ITERATION ',I4,':  SIMULTANEE SIMPLE')
2130  FORMAT('ITERATION ',I4,':  SIMULTANEOUS SIMPLE')
 140  FORMAT('ITERATION ',I4,':  M-ORTHONORMALISATION')
2140  FORMAT('ITERATION ',I4,':  M-ORTHONORMALIZATION')
 150  FORMAT('ITERATION ',I4,':  TEST sur les ANGLES entre VECTEURS')
2150  FORMAT('ITERATION ',I4,':  TEST on the ANGLES between VECTORS')
 160  FORMAT('ITERATION ',I4,':  TOUS LES ANGLES SUPERIEURS AU SEUIL')
2160  FORMAT('ITERATION ',I4,':  ALL ANGLES SUPERIORS to the LIMIT')
 170  FORMAT('ITERATION ',I4,':  ANALYSE de RAYLEIGHT-RITZ sur l ESPACE
     %ITERE')
2170  FORMAT('ITERATION ',I4,':  RAYLEIGHT-RITZ ANALYSIS on the ITERATED
     % SPACE')
 180  FORMAT(/80('*')/
     +       'ITERATION',I4,': MAXIMUM des ITERATIONS ATTEINT. ',
     +       'SORTIE avec des VALEURS PROPRES NON CONVERGEES'/
     +       'VOIR LES VALEURS PROPRES CALCULEES.'/
     +'CHOISIR UNE BORNE PLUS PROCHE (mais PLUS PETITE) DU MINIMUM DES V
     +ALEURS PROPRES'/
     +'ET      UNE BORNE PLUS GRANDE DU MAXIMUM DES VALEURS PROPRES'/
     +       80('*')/)
2180  FORMAT(/80('*')/
     +     'ITERATION',I4,':  MAXIMUM ITERATIONS REACHED.',
     +     ' EXIT with NOT CONVERGED EIGENVALUES.'/
     +     'SEE the COMPUTED EIGENVALUES.'/
     +     'CHOOSE A NEARER (but SMALLER) BOUND of EIGENVALUE MINIMUM'/
     +     'CHOOSE A GREATER BOUND of EIGENVALUE MAXIMUM'/
     +      80('*')/)
1190  FORMAT('ITERATION ',I4,':  Les VALEURS PROPRES CALCULEES sont '/
     +        5(I4,' : ',G15.7))
2190  FORMAT('ITERATION ',I4,':  The COMPUTED EIGENVALUES are '/
     +        5(I4,' : ',G15.7))
 200  FORMAT('ERREUR SSPAC1: MATRICE RAIDEUR NON INVERSIBLE apres ',
     +       '10 ESSAIS'/)
2200  FORMAT('ERREUR SSPAC1: STIFFNESS MATRIX NOT INVERSIBLE after ',
     +       '10 ESSAIS'/)
1210  FORMAT('FIN SSPAC1 a l''ITERATION',I4)
2210  FORMAT('END SSPAC1 at ITERATION',I4)
C
C     --------------------------------------------------------
C     INITIALISATIONS DES SEUILS DE PRECISION OU CONVERGENCE
C     RTOL EST LA TOLERANCE DE CONVERGENCE DES VALEURS PROPRES
C     --------------------------------------------------------
      IERR   = 0
      NIVEAU = 3
      RTOL   = 5.E-4
      TRANSL = 1.D-4
      ALPHA  = 0.5
      IPASS  = 0
      MOREE2 = MOTVAR(6)
C
C     ***********
C     * ETAPE 1 *  CALCUL DU SOUS ESPACE X0 DE DEPART (NC VECTEURS)
C     ***********
C
C     SI ON A UNE APPROXIMATION DES NC PREMIERS V.P. DANS VECP
C     ON L'UTILISE DIRECTEMENT (C'EST LE CAS INITX0=0)
      IF( LANGAG .EQ. 0 ) THEN
         WRITE (IMPRIM,100) INITX0
      ELSE
         WRITE (IMPRIM,2100) INITX0
      ENDIF
C
      IF( INITX0 .EQ. 1 ) THEN
C
C        INITX0=1 : METHODE BATHE
C        ------------------------
C        CALCUL DES QUOTIENTS DE RAYLEIGH V(I)=ABS( B(I,I)/A(I,I) )
C        INITIALISATION DE LA PREMIERE COLONNE DE VECP AVEC B(I,I)
C
         IF( NCODSB .NE. 0 ) THEN
C
C           A ET B ONT MEME PROFIL NON DIAGONAL
            DO J=1,N
               MUJ  = MUA(J+1)
CCC               V(J) = ABS( B(MUJ) / A(MUJ) )  remplace le 9/8/2006
               V(J) = B(MUJ) / A(MUJ)
               VECP(1,J) = B(MUJ)
            ENDDO
C
         ELSE
C
C           A EST NON DIAGONALE B EST DIAGONALE
            DO J=1,N
               MUJ  = MUA(J+1)
CCC               V(J) = ABS( B(J) / A(MUJ) ) remplace le 9/8/2006
               V(J) = B(J) / A(MUJ)
               VECP(1,J) = B(J)
            ENDDO
C
         ENDIF
C
C        INITIALISATION A 0 DES AUTRES COLONNES
         DO I=2,NC
            DO J=1,N
               VECP(I,J) = 0D0
            ENDDO
         ENDDO
C
C        ON CHOISIT PARMI LES VECTEURS DE BASE CEUX QUI ONT
C        LES PLUS GRANDS QUOTIENTS DE RAILEIGH
         DO I=2,NC
C           RECHERCHE DU PLUS GRAND QUOTIENT ACTUEL
            IJ = 0
            RT = 0D0
            DO J=1,N
               IJ = I
               IF( V(J) .GE. RT ) THEN
                  RT = V(J)
                  IJ = J
               ENDIF
            ENDDO
            V(IJ) = -1D11
            VECP(I,IJ) = 1D0
         ENDDO

      ELSE IF( INITX0 .EQ. 2 ) THEN

C        INITX0=2 : GENERATION ALEATOIRE DU SOUS ESPACE X0
C        -------------------------------------------------
         CALL ALEAD(3,N*NC,VECP)
C
      ELSE IF( INITX0 .EQ. 3 ) THEN
C
C        INITX0=3 : CONSTANTE SUR TOUTES LES COMPOSANTES SAUF SON NO
C        -----------------------------------------------------------
         DO I=1,NC
            DO J=1,N
               VECP(I,J) = 0.1D0
            ENDDO
            DO J=I,N,NC
               VECP(I,J) = 1D0
            ENDDO
         ENDDO
C
      ELSE IF( INITX0 .EQ. 4 ) THEN
C
C        INITX0=4 : GENERATION ALEATOIRE DU SOUS ESPACE X0
C        -------------------------------------------------
         CALL ALEAD(5,N*NC,VECP)
C
      ELSE IF( INITX0 .EQ. 5 ) THEN
C
C        INITX0=5 : GENERATION ALEATOIRE DU SOUS ESPACE X0
C        -------------------------------------------------
         CALL ALEAD(7,N*NC,VECP)
C
      ELSE IF( INITX0 .EQ. 6 ) THEN
C
C        INITX0=6 : GENERATION ALEATOIRE DU SOUS ESPACE X0
C        -------------------------------------------------
         CALL ALEAD(11,N*NC,VECP)
C
      ELSE IF( INITX0 .EQ. 7 ) THEN
C
C        INITX0=7 : ARRET DU TRAITEMENT
C        ------------------------------
         IF( LANGAG .EQ. 0 ) THEN
            WRITE (IMPRIM,110)
         ELSE
            WRITE (IMPRIM,2110)
         ENDIF
         IERR = 1
         RETURN
C
      ENDIF
C
C     ***********
C     * ETAPE 2 *  ITERATIONS DE SOUS ESPACES
C     ***********
C
C     DEBUT DU TRAITEMENT: VECP CONTIENT LES VECTEURS DE DEPART
C     ON PREPARE LE CYCLE DES ITERATIONS  R = B * VECP
C     -------------------------------------------------------
      CALL APBEBD( NCODSB,MUB,B,VECP,NC,N, R )
C
      IF( IFACTO .EQ. 0 ) THEN
 6       IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10006) -TETATO
         ELSE
            WRITE(IMPRIM,20006) -TETATO
         ENDIF
10006    FORMAT(/'FACTORISATION L*D*tL DE K +',G15.7,' M')
20006    FORMAT(/'L*D*tL FACTORIZATION of K +',G15.7,' M')
C
         CALL TRTATA( RG, A, MOTSRG )
         CALL MUA2PD( N,   1D0, NCODSA, MUA, A,
     +                 -TETATO, NCODSB, MUB, B,
     +                                  MUA, A )
         CALL CRMC1D( MUA,A,N,EPS,0, A,NRETOU )
         IFACTO = 1
         NB     = 0
         DO 57 I=1,N
            IF( A(MUA(I+1)) .LT. 0D0 ) NB = NB + 1
 57      ENDDO
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10057) NB,TETATO
         ELSE
            WRITE(IMPRIM,20057) NB,TETATO
         ENDIF
10057    FORMAT(I4,' VALEURS PROPRES AVANT',G15.7)
20057    FORMAT(I4,' EIGENVALUES before ',G15.7)
         IF( NRETOU .NE. 0 ) THEN
            IPASS = IPASS + 1
            IF( IPASS .EQ. 10 ) THEN
               IF( LANGAG .EQ. 0 ) THEN
                  WRITE(IMPRIM,200)
               ELSE
                  WRITE(IMPRIM,2200)
               ENDIF
               IERR = 2
               RETURN
            ENDIF
            TETATO = TETATO * (1D0+TRANSL)
            IF( TETATO .EQ. 0D0 ) TETATO = 1D0
            GOTO 6
         ENDIF
      ENDIF
C
C     INITIALISATION DES VALEURS PROPRES
C     ----------------------------------
      IF( INITX0 .NE. 0 ) THEN
         DO 7 I=1,NC
            VALP(I) = 0D0
 7       ENDDO
      ENDIF
C
C     LES VALEURS PROPRES A L'ITERATION PRECEDENTE
      DO 77 I=1,NC
         D(I) = VALP(I)
 77   ENDDO
C
C     STRATEGIE DES ITERATIONS GEREE PAR LE TABLEAU MANAG
C     ---------------------------------------------------
      IF( NOSTRA .LE. 0 ) THEN
C
C        STRATEGIE AUTOMATIQUE
C        ---------------------
         K        = 0
         NITEM    = MOMANA
         MANAG(1) = 1
         MANAG(2) = 3
         MANAG(3) = 1
         MANAG(4) = 1
CCC         MANAG(5) = 4   => MOINS RAPIDE QUE 5 CAR MOINS D' ITERATIONS
         MANAG(5) = 5
         MANAG(6) = 1
         MANAG(7) = 2
         MANAG(8) = 3
         MANAG(9) = 5
         MANAG(10) = 1
         MANAG(11) = 2
         MANAG(12) = 3
         MANAG(13) = 5
         MANAG(14) = 1
         MANAG(15) = 2
         MANAG(16) = 3
         DO 8 I=17,NITEM
            MANAG(I) = 5
 8       ENDDO
         IF( INITX0 .EQ. 0 ) THEN
C           REPRISE D'UN PRECEDENT CALCUL
            K = 10
         ENDIF
C
      ELSE
C
C        STRATEGIE DE L UTILISATEUR
C        --------------------------
         K     = 0
         NITEM = NOSTRA
         DO 155 I=1,NITEM
            IF( MANAG(I).LE.0 .OR. MANAG(I).GT.5 ) MANAG(I)=5
 155     ENDDO
         MANAG(NITEM) = 5
      ENDIF
C
      IF( LANGAG .EQ. 0 ) THEN
         WRITE (IMPRIM,120) NITEM,(I,MANAG(I),I=1,NITEM)
      ELSE
         WRITE (IMPRIM,2120) NITEM,(I,MANAG(I),I=1,NITEM)
      ENDIF
      WRITE (IMPRIM,*)
C
C     -----------------------------
C     ITERATIONS DIRIGEES PAR MANAG
C     -----------------------------
 1020 K = K + 1
      IF ( MANAG(K) .LT. 5 ) THEN
C
C        ++++++++++++++++++++++++++++++
C        1: ITERATION SIMULTANEE SIMPLE
C        ++++++++++++++++++++++++++++++
         IF( LANGAG .EQ. 0 ) THEN
            WRITE (IMPRIM,130) K
         ELSE
            WRITE (IMPRIM,2130) K
         ENDIF
C
C        ON CALCULE XK+1=VECP tel que L*D*tL*XK+1 = A*XK+1 = B*XK = R
         CALL DRCR1D( N,NCODSA,MUA,A,NC,R,NIVEAU, VECP,IERR )
C
 1030    IF( MANAG(K) .NE. 1 .AND. MANAG(K) .NE. 4 ) THEN
C
C           ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C           2:ITERATION SIMULTANEE AVEC B ORTHONORMALISATION+VALEURS PROPRES
C           3:IDEM MAIS SANS CALCUL DES VALEURS PROPRES
C           ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            IF( LANGAG .EQ. 0 ) THEN
               WRITE (IMPRIM,140) K
            ELSE
               WRITE (IMPRIM,2140) K
           ENDIF
C           ON REMPLACE R PAR R = VECP B-ORTHONORME
C           R = XK+1 / NORME(XK+1)
            CALL ORTHOD( N,MUB,B,NCODSB,COMXMG, VECP0,NCVALP0,NBVALP0,
     +                   VECP,NC, V,W, VALP,R )
C
            IF( MANAG(K) .EQ. 3 ) THEN
C              SI MANAG(K)=3 SAUT DU CALCUL DES VALEURS PROPRES
               CALL TRTATA( R, VECP, MOREE2*NC*N )
               GOTO 1040
            ENDIF
C
C           CALCUL DES VALEURS PROPRES A PARTIR DU QUOTIENT DE RAYLEIGHT
C           VALP  = ( tR * A * R ) / (tR * B * R )
C                 = ( tR * L * D * tL * R ) = ( D * tL * R , tL * R)
            DO 91 J=1,N
               DO 90 I=1,NC
                  VECP(I,J) = 0D0
 90            ENDDO
 91         ENDDO
C
C           CALCUL DE tL * R
            IA = 0
            DO 11 J=2,N+1
               I2 = J-1
               I  = J-MUA(J)+MUA(I2)
               IJ = I2-1
               DO 10 I3=I,IJ
                  IA = IA+1
                  S  = A(IA)
                  DO 94 IL=1,NC
                     VECP(IL,I3) = VECP(IL,I3) + S * R(IL,I2)
 94               ENDDO
 10            ENDDO
C
C              LE COEFFICIENT DIAGONAL DE tL=1
               IA = IA+1
               DO 95 IL=1,NC
                  VECP(IL,I2) = VECP(IL,I2) + R(IL,I2)
 95            ENDDO
 11         ENDDO
C
C           CALCUL DE ( D * tL * R , tL * R )
            DO 13 J=1,NC
               S = 0D0
               DO 12 I=1,N
                  S = S + A(MUA(I+1)) * VECP(J,I) ** 2
 12            ENDDO
               VALP(J) = S
 13         ENDDO
C           IMPRESSION DES VALEURS PROPRES
            IF( LANGAG .EQ. 0 ) THEN
               WRITE (IMPRIM,1190) K,(I,VALP(I)+TETATO,I=1,NC)
            ELSE
               WRITE (IMPRIM,2190) K,(I,VALP(I)+TETATO,I=1,NC)
            ENDIF
C
C           D CONTIENT LES VALEURS PROPRES DE L ITERATION PRECEDENTE
C           ---------------------------------------------------------
            DO 15 I=1,NC
               D(I) = VALP(I)
 15         ENDDO
         ENDIF
C
C        ON PREPARE L'ETAPE SUIVANTE : R = B * VECP
C        ------------------------------------------
 1040    CALL APBEBD( NCODSB,MUB,B,VECP,NC,N, R )
C
         IF ( MANAG(K) .NE. 4 ) GOTO 1020
C
C        +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C        4: TEST SUR LES ANGLES (VECP(1),VECP(I))>ALPHA PAS DE B-ORTHONORMAL
C        +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C        CALCUL DE  tV1*B*V1=T
C
         IF( LANGAG .EQ. 0 ) THEN
            WRITE (IMPRIM,150) K
         ELSE
            WRITE (IMPRIM,2150) K
         ENDIF
         T = 0D0
         DO 16 I=1,N
            T = T + VECP(1,I) * R(1,I)
 16      ENDDO
         DO 18 I=2,NC
C           CALCUL DE ANGLE = 1-ABS(tV1*B*VI) / sqrt((tVI*B*VI)*T))
            S = 0D0
            U = 0D0
            DO 17 J=1,N
               S = S + VECP(1,J) * R(I,J)
               U = U + VECP(I,J) * R(I,J)
 17         ENDDO
            U = SQRT( U * T )
C           ANGLE < ALPHA?
            IF( 1D0 - ABS( S/U ) .LT. ALPHA ) THEN
               MANAG(K) = 3
               GOTO 1030
            ENDIF
 18      ENDDO
         IF( LANGAG .EQ. 0 ) THEN
            WRITE (IMPRIM,160) K
         ELSE
            WRITE (IMPRIM,2160) K
         ENDIF
         GOTO 1020
      ENDIF
C
C     +++++++++++++++++++++++++++++++++++++++++++++
C     5:ANALYSE DE RAYLEIGH-RITZ SUR L ESPACE ITERE
C     +++++++++++++++++++++++++++++++++++++++++++++
      IF( LANGAG .EQ. 0 ) THEN
         WRITE (IMPRIM,170) K
      ELSE
         WRITE (IMPRIM,2170) K
      ENDIF
C
C     CALCUL DES OPERATEURS PROJETES :AR=tVECP*A*VECP et BR=tVECP*B*VECP
C     D'ABORD  A PARTIR DE R ON TROUVE VECP TEL QUE A*VECP=R
C     C.A.D: A*XK+1BARRE = YK = B*XK = R ET ON AURA:VECP = XK+1BARRE
C     ------------------------------------------------------------------
      CALL DRCR1D( N, NCODSA, MUA, A, NC, R, NIVEAU,  VECP, IERR )
C
C     MAINTENANT AR = tVECP * R = t(XK+1BARRE) * A * XK+1BARRE
      DO 20 I=1,NC
         DO 87 J=1,I
            S = 0D0
            DO 19 L=1,N
               S = S + R(J,L) * VECP(I,L)
 19         ENDDO
            AR(J,I) = S
            AR(I,J) = S
 87      ENDDO
 20   ENDDO
C
C     ON REMPLACE R PAR R = B * VECP = B * XK+1
      CALL APBEBD( NCODSB, MUB, B, VECP, NC, N, R )
C
C     MAINTENANT BR = tVECP * R = t(XK+1BARRE) * B * XK+1BARRE
      DO 22 I=1,NC
         DO 85 J=1,I
            S = 0D0
            DO 21 L=1,N
               S = S + R(J,L) * VECP(I,L)
 21         ENDDO
            BR(I,J) = S
            BR(J,I) = S
 85      ENDDO
 22   ENDDO
C
C     RESOLUTION DU PROBLEME PROJETE AR(X)=L*BR(X)
C     VECJAC=QK+1 EST FORME DES NC VECTEURS PROPRES
C     ---------------------------------------------
      CALL JACOBD( NC,   RTOL, AR, BR, DAUX,
     %             VALP, VECJAC, IERR )
C
C     REORDONNANCEMENT DES VALEURS ET VECTEURS PROPRES SUIVANT
C     LEUR MODULE CROISSANT
      CALL ORDRED( NC, VALP, MUTA )
      IF( LANGAG .EQ. 0 ) THEN
         WRITE (IMPRIM,1190) K,(I,VALP(I)+TETATO,I=1,NC)
      ELSE
         WRITE (IMPRIM,2190) K,(I,VALP(I)+TETATO,I=1,NC)
      ENDIF
C     POUR ECONOMISER UN PEU ON ECRIT VECJAC REORDONNEE DANS AR
      DO 23 I=1,NC
         IMA = MUTA(I)
         DO 55 J=1,NC
            AR(J,I) = VECJAC(J,IMA)
 55      ENDDO
 23   ENDDO
C
C     CALCULONS LES VECTEURS PROPRES: XK=(XK+1BARRE)*QK+1
C     C.A.D: R = VECP * AR
      DO 25 I=1,N
         DO 58 J=1,NC
            S = 0D0
            DO 24 L=1,NC
               S = S + VECP(L,I) * AR(L,J)
 24         ENDDO
            R(J,I) = S
 58      ENDDO
 25   ENDDO
C
C     B-ORTHOGONALISATION DES VECP PAR RAPPORT A VECP0 ET EUX-MEMES
      CALL ORTHOD( N,MUB,B,NCODSB,COMXMG, VECP0,NCVALP0,NBVALP0,
     %             R,NC, V,W, VALP,VECP )
C
C     TEST DE CONVERGENCE: CALCUL DU NOMBRE DE VALEURS PROPRES CONVERGEES
C     -------------------
      NBVALC = 0
      DO 26 I=1,NC-5
         IF( ABS( VALP(I)-D(I) ) .GT. RTOL*ABS( VALP(I) ) ) GOTO 1050
         NBVALC = I
 26   ENDDO
C     CONVERGENCE ASSUREE
      GOTO 2000
C
C     LIMITATION DES ITERATIONS
C     -------------------------
 1050 IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,11050) K,NBVALC
      ELSE
         WRITE(IMPRIM,21050) K,NBVALC
      ENDIF
11050 FORMAT('ITERATION',I4,':  NOMBRE de VALEURS PROPRES CONVERGEES='
     %,I4)
21050 FORMAT('ITERATION',I4,':  NUMBER of CONVERGED EIGENVALUES=',I4)
      IF( K .GE. NITEM ) THEN
C        ICI: NOMBRE MAXIMAL D'ITERATIONS ATTEINT => NON CONVERGENCE
         IF( LANGAG .EQ. 0 ) THEN
            WRITE (IMPRIM,180) NITEM
         ELSE
            WRITE (IMPRIM,2180) NITEM
         ENDIF
         IERR = 3
C        SAUVEGARDE DANS VVPR DES VALEURS ET VECTEURS PROPRES CONVERGES
         GOTO 2000
      ENDIF
C
C     ON PREPARE L ETAPE SUIVANTE
C     ---------------------------
      WRITE(IMPRIM,*)
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10027) K+1
      ELSE
         WRITE(IMPRIM,20027) K+1
      ENDIF
10027 FORMAT('DEBUT DE L''ITERATION',I4)
20027 FORMAT('START of the ITERATION',I4)
C
      DO 27 I=1,NC
         D(I) = VALP(I)
 27   ENDDO
C
C     PERTURBATION DES VECTEURS PROPRES NC-3,...,NC
C     POUR INJECTER DE NOUVELLES DIRECTIONS ET CELA
C     POUR RETROUVER DES VALEURS PROPRES OUBLIEES
C     ---------------------------------------------
      ISIGNE = 1
      DO 50 I=NC-3,NC
C        CALCUL DE LA COMPOSANTE MAX
         S = ABS( VECP(I,1) )
         DO 40 J=2,N
            T = ABS( VECP(I,J) )
            IF( T .GT. S ) S = T
 40      ENDDO
C        CES VALEURS POUR NE PAS AVOIR DE VALEURS REPETITIVES
         S  = S / ( 2 + MOD(K,5) )
         DD = (NC-I+100) * S / N
         DO 45 J=1,N
            IF( ABS(VECP(I,J)) .LT. S ) THEN
               ISIGNE = - ISIGNE
               VECP(I,J)= ISIGNE * J * DD
            ENDIF
 45      ENDDO
 50   ENDDO
C
C     ON CALCULE  YK+1 = B * XK+1 C.A.D  R = B * VECP
      CALL APBEBD( NCODSB,MUB,B,VECP,NC,N, R )
      GOTO 1020
C
C     ====================
C     CONVERGENCE REALISEE
C     ====================
 2000 IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,12000) K,NBVALC
      ELSE
         WRITE(IMPRIM,22000) K,NBVALC
      ENDIF
12000 FORMAT('ITERATION',I4,': CONVERGENCE des',I4,' VALEURS PROPRES')
22000 FORMAT('ITERATION',I4,': CONVERGENCE of',I4,' EIGENVALUES')
C
C     REORDONNANCEMENT DES VALEURS ET VECTEURS PROPRES
C     ------------------------------------------------
      CALL REORDD( 1, NC, N, VALP, VECP )
C
C     TRANSLATION INVERSE
C     -------------------
      DO I=1,NC
         VALP(I) = VALP(I) + TETATO
      ENDDO
C     VALEUR PROPRE ESTIMEE SUIVANT LA DERNIERE CONVERGEE
      VALPES = VALP(NBVALC+1)
C
C     IMPRESSIONS DIVERSES
C     --------------------
      IF( LANGAG .EQ. 0 ) THEN
         WRITE (IMPRIM,1190) K, (I,VALP(I),I=1,NBVALC)
         WRITE (IMPRIM,1210) K
      ELSE
         WRITE (IMPRIM,2190) K, (I,VALP(I),I=1,NBVALC)
         WRITE (IMPRIM,2210) K
      ENDIF
      WRITE (IMPRIM,*)

      RETURN
      END
