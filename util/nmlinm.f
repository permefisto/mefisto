      SUBROUTINE NMLINM( NMTOBJ, NBNO, NUMERO, NOSENS, KNOM )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    FORMER UN NOM DE 24 CARACTERES A PARTIR D'UNE LISTE DE NOMS
C -----    OBTENUS PAR LEUR NUMERO DANS UN LEXIQUE ET SELON LE SENS
C          CROISSANT OU DECROISSANT
C
C ENTREES:
C --------
C NMTOBJ : NOM DU TYPE DE L'OBJET 'POINT', 'LIGNE', ...
C NBNO   : NOMBRE DE NUMEROS DANS LE LEXIQUE
C NUMERO : NBNO NUMEROS DANS LE LEXIQUE
C NOSENS : +1 SENS CROISSANT, -1 SENS DECROISSANT
C
C SORTIE :
C --------
C KNOM   : CONCATENATION DES NOMS AVEC REDUCTION POUR TENIR
C          DANS 24 CARACTERES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     FEVRIER 1997
C2345X7..............................................................012
      CHARACTER*(*)  NMTOBJ
      CHARACTER*24   KNOM, NOM
      INTEGER        NUMERO(1:NBNO)
C
      NBTC = 0
      LMIN = 24
C
      IF( NOSENS .GE. 0 ) THEN
         I1 = 1
         I2 = NBNO
      ELSE
         I1 = NBNO
         I2 = 1
      ENDIF
C
      KNOM = '                        '
      DO 10 I = I1, I2, NOSENS
C
C        LE NOM DU NUMERO(I) DANS LE LEXIQUE NMTOBJ
         CALL NMOBNU( NMTOBJ, NUMERO(I), NOM )
C        LE NOMBRE DE CARACTERES NON BLANC DU NOM I
         L    = NUDCNB( NOM )
C        LE NOMBRE MINIMAL DE CARACTERES PAR NOM
         LMIN = MIN( LMIN, L )
C        LE NOM CONCATENE
         L1   = NBTC+1
         L2   = MIN(NBTC+L,24)
         IF( L1.LE.24 .AND. L2.LE.24 ) KNOM(L1:L2) = NOM(1:L)
C        LE NOMBRE TOTAL DE CARACTERES
         NBTC = NBTC + L
C
 10   CONTINUE
      IF( NBTC .LE. 24 ) RETURN
C
C     TROP DE CARACTERES => IL FAUT EN ELIMINER
      IF( LMIN * NBNO .LE. 24 ) THEN
C
C        CHOIX DES MIN OU MOY PREMIERS CARACTERES DES NOMS
         MOY  = 24 / NBNO
         IF( MOY .LE. 0 ) GOTO 30
         MOY  = MOY + 1
         KNOM = '                        '
         NBTC = 0
         DO 20 I = I1, I2, NOSENS
C
C           LE NOM DU NUMERO(I) DANS LE LEXIQUE NMTOBJ
            CALL NMOBNU( NMTOBJ, NUMERO(I), NOM )
C           LE NOMBRE DE CARACTERES NON BLANC DU NOM I
            L = NUDCNB( NOM )
            IF( L .EQ. LMIN ) GOTO 15
            L = MIN( L, MOY )
C           LE NOM CONCATENE
 15         L1 = NBTC+1
            L2 = MIN(NBTC+L,24)
            IF( L1.LE.24 .AND. L2.LE.24 ) KNOM(L1:L2) = NOM(1:L)
C           LE NOMBRE TOTAL DE CARACTERES DE KNOM
            NBTC = NBTC + L
            IF( NBTC .GE. 24 ) RETURN
C
 20      CONTINUE
      ENDIF
      IF( NBTC .LE. 24 ) RETURN
C
C     ARRET DES QUE 24 EST ATTEINT
 30   NBTC = 0
      DO 40 I = I1, I2, NOSENS
C
C        LE NOM DU NUMERO(I) DANS LE LEXIQUE NMTOBJ
         CALL NMOBNU( NMTOBJ, NUMERO(I), NOM )
C        LE NOMBRE DE CARACTERES NON BLANC DU NOM I
         L = NUDCNB( NOM )
C        LE NOM CONCATENE
         L1 = NBTC+1
         L2 = MIN(NBTC+L,24)
         IF( L1.LE.24 .AND. L2.LE.24 ) KNOM(L1:L2) = NOM(1:L)
C        LE NOMBRE TOTAL DE CARACTERES
         NBTC = NBTC + L
         IF( NBTC .GE. 24 ) RETURN
C
 40   CONTINUE
      END
