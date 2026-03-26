      SUBROUTINE BLDLP0( NTDL, NBDLFX, NODLFX, NCODSA, LPDIAG, AG )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : IMPOSER DES FIXATIONS NULLES AUX DEGRES DE LIBERTE
C ----- DE LA LISTE NODLFX SUR UNE MATRICE PROFIL TOUT EN MC
C       LES COEFFICIENTS DE LA LIGNE ET DE LA COLONNE DE CHACUN DES
C       DL FIXES SONT ANNULES SAUF LE COEFFICIENT DIAGONAL QUI DEVIENT 1
C
C       | AG 11  AG 12 | | XG 1 |   | BG 1 |
C                        |      | =           DEVIENT POUR LA MATRICE AG
C                        |  0   |                    SEULEMENT
C
C       | AG 11    0   | | XG 1 |   | BG 1           |
C       |              | |      | = |                |
C       |  0  Identite | |  0   |   | non defini ici |
C
C ENTREES:
C --------
C NTDL   : NOMBRE TOTAL DE DEGRES DE LIBERTE
C NBDLFX : NOMBRE DE DEGRES DE LIBERTE FIXES
C NODLFX : NO DES NBDLFX DEGRES DE LIBERTE FIXES
C NCODSA : CODE DE STOCKAGE DE LA MATRICE
C          1 : MATRICE SYMETRIQUE
C          0 : MATRICE DIAGONALE => RIEN A FAIRE CAR AG 12 = 0
C         -1 : MATRICE NON SYMETRIQUE
C LPDIAG : POINTEUR SUR LES COEFFICIENTS DIAGONAUX DE LA MATRICE AG
C
C MODIFIES:
C ---------
C AG     : MATRICE PROFIL EN MEMOIRE CENTRALE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LABORATOIRE ANALYSE NUMERIQUE PARIS AOUT 1998
C23456---------------------------------------------------------------012
      DOUBLE PRECISION  AG(*)
      INTEGER           LPDIAG(0:NTDL),NODLFX(NBDLFX)
C
      IF( NBDLFX .LE. 0 .OR. NCODSA .EQ. 0 ) RETURN
      JA1 = 0
C
C     CALCUL DE LBDPDL : DEMI LARGEUR DE BANDE DE LA MATRICE PROFIL
      LBDPDL = 0
      DO 5 I=1,NTDL
         LBDPDL = MAX( LBDPDL , LPDIAG(I)-LPDIAG(I-1) )
 5    CONTINUE
C
      DO 100 I=1,NBDLFX
C        LE NUMERO GLOBAL DU DL FIXE I
         NODLF  = NODLFX(I)
         IBM = MIN(NTDL,NODLF+LBDPDL)
C        LE NO DU COEFFICIENT DIAGONAL PRECEDENT
         IA0 = LPDIAG(NODLF-1)
C        LE NO DU COEFFICIENT DIAGONAL
         ID  = LPDIAG( NODLF )
         IH  = ID - IA0 - 1
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
         DO 20 L=IMI,NODLF-1
            IA  = IA + 1
C           LE COEFFICIENT DANS LA MATRICE DEVIENT NUL
            AG(IA) = 0D0
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
C
C        LES COEFFICIENTS AG(NODLF,L) AG(L,NODLF) NODLF+1=<L=<NTDL
C        ---------------------------------------------------------
         IF( NODLF .GE. NTDL ) GOTO 100
         JA0 = LPDIAG( NODLF )
         DO 60 J=NODLF+1,IBM
            JA1 = LPDIAG(J)
C           LA HAUTEUR DES COEFFICIENTS
            JH  = JA1 - JA0 - 1
            IF( NCODSA .LT. 0 ) JH = JH / 2
            JMI = J - JH
            IF( JMI .LE. NODLF ) THEN
C              ADRESSAGE DU COEFFICIENT DANS LE PROFIL
               IF( NCODSA .LT. 0 ) THEN
C                 CALCUL DE L ADRESSE DE AG(J,NODLF)
C                 MATRICE NON SYMETRIQUE
                  IA  = JA0 - J + NODLF + JH + 1
               ELSE
C                 MATRICE SYMETRIQUE
                  IA  = JA1 - J + NODLF
               ENDIF
C              LE COEFFICIENT DANS LA MATRICE DEVIENT NUL
               AG(IA) = 0D0
C              CAS D'UNE MATRICE NON SYMETRIQUE
               IF( NCODSA .LT. 0 ) AG(JA1 - J + NODLF ) = 0D0
            ENDIF
            JA0 = JA1
 60      CONTINUE
 100  CONTINUE
      END
