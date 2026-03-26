      SUBROUTINE BLDLPC( NTDL,   NDSM, NBDLFX, NODLFX, VADLFX,
     %                   NCODSA, LPDIAG, AG,  BG )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : METTRE DANS LE SECOND MEMBRE LA CONTRIBUTION DES DL FIXES
C ----- DE LA LISTE NODLFX SUR UNE MATRICE PROFIL TOUT EN MC
C       LES COEFFICIENTS DES LIGNES DES DL FIXES SONT ANNULES SAUF
C       LE COEFFICIENT DIAGONAL QUI DEVIENT 1
C
C       | AG 11  AG 12 | | XG 1 |   | BG 1 |
C                        |      | =            DEVIENT
C                        | XD 2 |
C
C       | AG 11    0   | | XG 1 |   | BG 1 - AG 12 x XD 2 |
C       |              | |      | = |                     |
C       |  0  Identite | | XG 2 |   | XD 2                |
C
C DONNEES:
C --------
C NTDL   : NOMBRE TOTAL DE DEGRES DE LIBERTE
C NDSM   : NOMBRE DE CAS DE CHARGE
C NBDLFX : NOMBRE DE DEGRES DE LIBERTE FIXES
C NODLFX : NO DES NBDLFX DEGRES DE LIBERTE FIXES
C VADLFX : VALEUR DES DEGRES DE LIBERTE FIXES
C NCODSA : CODE DE STOCKAGE DE LA MATRICE
C          1 : MATRICE SYMETRIQUE
C          0 : MATRICE DIAGONALE => RIEN A FAIRE CAR AG 12 = 0
C         -1 : MATRICE NON SYMETRIQUE
C LPDIAG : POINTEUR SUR LES COEFFICIENTS DIAGONAUX DE LA MATRICE AG
C
C MODIFIES:
C ---------
C AG     : MATRICE PROFIL EN MEMOIRE CENTRALE
C BG     : SECOND MEMBRE GLOBAL EN MEMOIRE CENTRALE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LABORATOIRE ANALYSE NUMERIQUE PARIS  MAI 1989
C23456---------------------------------------------------------------012
      DOUBLE PRECISION AG(*), BG(NTDL, NDSM), VADLFX(NDSM, NBDLFX), X
      INTEGER          LPDIAG(0:NTDL), NODLFX(NBDLFX)
C
      IF( NBDLFX .LE. 0 .OR. NCODSA .EQ. 0 ) RETURN
      JA1 = 0
C
C     CALCUL DE LBDPDL : DEMI LARGEUR DE BANDE DE LA MATRICE PROFIL
      LBDPDL = 0
      DO 5 I=1, NTDL
         LBDPDL = MAX( LBDPDL,  LPDIAG(I)-LPDIAG(I-1) )
 5    CONTINUE
C
      DO 100 I=1, NBDLFX
C
C        LE NUMERO GLOBAL DU DL FIXE I
         NODLF = NODLFX(I)
C        LE NO DU COEFFICIENT DIAGONAL PRECEDENT
         IA0 = LPDIAG(NODLF-1)
C        LE NO DU COEFFICIENT DIAGONAL
         ID  = LPDIAG( NODLF )
         IH  = ID - IA0 - 1
C
         IF( NCODSA .LT. 0 ) THEN
C           PARTICULARITES POUR AG NON SYMETRIQUE
            JA1 = IA0
            IH  = IH / 2
            IA0 = IA0 + IH
         ENDIF
C
         IA  = IA0
         IMI = NODLF - IH
C
C        COEFFICIENTS AG(L,NODLF) ET AG(NODLF,L)  LPDIAG(NODLF-1)=<L<NODLF
C        -----------------------------------------------------------------
         DO 20 L=IMI, NODLF-1
            IA  = IA + 1
            X   = AG(IA)
            IF( X .NE. 0D0 ) THEN
C
C              LE COEFFICIENT MATRICIEL EST NON NUL
C              SA CONTRIBUTION EST REPORTEE DANS LE SECOND MEMBRE
               DO 10 K=1,NDSM
                  BG( L, K ) = BG( L, K ) - X * VADLFX( K, I )
 10            CONTINUE
C              LE COEFFICIENT DANS LA MATRICE DEVIENT NUL
               AG(IA) = 0D0
            ENDIF
C
            IF( NCODSA .LT. 0 ) THEN
C              CAS D'UNE MATRICE NON SYMETRIQUE
               JA1 = JA1 + 1
               AG(JA1) = 0D0
            ENDIF
 20      CONTINUE
C
C        LE COEFFICIENT DIAGONAL DEVIENT L'IDENTITE
         AG(ID) = 1D0
         DO 30 K=1, NDSM
            BG( NODLF, K ) = VADLFX( K, I )
 30      CONTINUE
C
C        LES COEFFICIENTS AG(NODLF,L) AG(L,NODLF) NODLF+1=<L=<NTDL
C        ---------------------------------------------------------
         IF( NODLF .GE. NTDL ) GOTO 100
         JA0 = LPDIAG( NODLF )
         DO 60 J = NODLF+1, MIN(NTDL, NODLF+LBDPDL)
C
C           DEPART DU COEFFICIENT DIAGONAL J
            JA1 = LPDIAG(J)
C           LA HAUTEUR DES COEFFICIENTS
            JH  = JA1 - JA0 - 1
            IF( NCODSA .LT. 0 ) JH = JH / 2
            JMI = J - JH
            IF( JMI .LE. NODLF ) THEN
C
C              ADRESSAGE DU COEFFICIENT DANS LE PROFIL
               IF( NCODSA .LT. 0 ) THEN
C                 CALCUL DE L ADRESSE DE AG( J, NODLF )
C                 MATRICE NON SYMETRIQUE
                  IA = JA0 - J + NODLF + JH + 1
               ELSE
C                 MATRICE SYMETRIQUE
                  IA = JA1 - J + NODLF
               ENDIF
C
C              LA CONTRIBUTION DU COEFFICIENT MATRICIEL NON NUL
C              EST REPORTEE DANS LE SECOND MEMBRE
               DO 40 K=1, NDSM
                  BG( J, K ) = BG( J, K ) - AG( IA ) * VADLFX( K, I )
 40            CONTINUE
C              LE COEFFICIENT DANS LA MATRICE DEVIENT NUL
               AG(IA) = 0D0
C              CAS D'UNE MATRICE NON SYMETRIQUE
               IF( NCODSA .LT. 0 ) AG( JA1 - J + NODLF ) = 0D0
            ENDIF
            JA0 = JA1
 60      CONTINUE
C
 100  CONTINUE
      RETURN
      END
