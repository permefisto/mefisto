      SUBROUTINE REELCA( KVAL0, KVAL1 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT: ECRIRE SANS '  -0.XXXE+-EXPOSANT10 ' UN NOMBRE REEL SANS EXPOSANT
C ---- QUAND C'EST POSSIBLE '  -1.2E+03 ' => '-1200'
C                           '   1.2E-03 ' => '0.0012'
C ENTREE :
C --------
C KVAL0  : CHAINE DE CARACTERES REPRESENTANT UN NOMBRE REEL

C SORTIE :
C --------
C KVAL1  : CHAINE DE CARACTERES REPRESENTANT UN NOMBRE REEL
C          SANS EXPOSANT DE 10 LORSQUE C'EST POSSIBLE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET Alain LJLL UPMC & St Pierre du Perray  Janvier 2012
C2345X7..............................................................012
      CHARACTER*(*)  KVAL0, KVAL1
      CHARACTER*1    CAR
      INTRINSIC      LEN
C
      L1 = LEN( KVAL1 )
      KVAL1 = '                      '
C
      L0 = LEN( KVAL0 )
      L  = NUDCNB( KVAL0 )

ccc      DO M=1,L
ccc         CAR = KVAL0(M:M)
ccc         print *,'Reelca: car0=',car,' M=',M
ccc      ENDDO

C     EXISTE T IL UN CARACTERE 'E' ?
      CAR = KVAL0(L-3:L-3)
      IF( CAR .NE. 'E' ) GOTO 7000

C     QUEL EST L'EXPOSANT DE 10?
      READ( KVAL0(L-2:L), '(I3)') LEXP
C
      IF( LEXP .GT. 4 ) GOTO 7000

C     DECALAGES SI POSSIBLE
      IF( LEXP .GT. 0 ) THEN
C
C        EXPOSANT DE 10 POSITIF < 4
         IF( KVAL0(2:2) .EQ. '-' ) THEN
            KVAL1(1:1) = '-'
            L0 = 3
         ELSE
            L0 = 4
         ENDIF

         DO M=5,L-4
            KVAL1(M-L0:M-L0) = KVAL0(M:M)
         ENDDO

C        NOMBRE DE CHIFFRES SORTIS DU . DECIMAL
         NBC = L-4-5+1
         DO M=L-3,L-4+LEXP-NBC
            KVAL1(M-L0:M-L0) = '0'
         ENDDO

      ELSE IF( LEXP .LT. 0 ) THEN

C        EXPOSANT NEGATIF
         IF( KVAL0(3:3) .EQ. '-' ) THEN
            L0 = 3
            KVAL1(1:3) = '-0.'
         ELSE
            L0 = 2
            KVAL1(1:2) = '0.'
         ENDIF

         DO M=L0+1,L0-LEXP
            KVAL1(M:M) = '0'
         ENDDO

         DO M=5,L-4
            MM=L0-LEXP+M-4
            KVAL1(MM:MM) = KVAL0(M:M)
         ENDDO

      ELSE

C        PAS DE CHANGEMENT
         GOTO 7000

      ENDIF

C     SUPPRESSION DES .000 FINAUX
      L = NUDCNB( KVAL1 )
      K = INDEX( KVAL1(1:L), 'E' )
      IF( K .EQ. 0 ) THEN
         M = INDEX( KVAL1(1:L), '.' )
         IF( M .GT. 0 ) THEN
            DO K = L, M+1 , -1
               IF( KVAL1(K:K) .EQ. '0' ) THEN
                  KVAL1(K:K)=' '
               ELSE
                  GOTO 9000
               ENDIF
            ENDDO
            KVAL1(M:M)=' '
         ENDIF
      ENDIF
      GOTO 8000

C     PAS DE CHANGEMENT DU NOMBRE AFFICHE
 7000 DO M=1,L
         CAR = KVAL0(M:M)
         KVAL1(M:M) = CAR
      ENDDO

C     SUPPRESSION DU DERNIER ZERO SI PAS DE E
 8000 L = NUDCNB( KVAL1 )
      IF( INDEX( KVAL1(1:L), 'E' ) .EQ. 0 ) THEN
         IF( KVAL1(L:L) .EQ. '0' ) THEN
            KVAL1(L:L) = ' '
            GOTO 8000
         ENDIF
      ENDIF

C     SUPPRESSION DES BLANCS DANS LE RESULTAT KVAL1
      LDECAL = 0
      DO M=1,L
         CAR = KVAL1(M:M)
ccc         print *,'Reelca: car1=',car,' M=',M
         IF( CAR .EQ. ' ' ) THEN
            LDECAL = LDECAL + 1
         ELSE
            KVAL1(M-LDECAL:M-LDECAL) = CAR
         ENDIF
      ENDDO

C     MISE A BLANC DES CARACTERES RESTANTS
      DO M=L-LDECAL+1,L1
         KVAL1(M:M) = ' '
      ENDDO

C     SI LE DERNIER CARACTERE EST '.' IL EST SUPPRIME
      M = NUDCNB( KVAL1 )
      IF( KVAL1(M:M) .EQ. '.' ) THEN
         KVAL1(M:M) = ' '
      ENDIF

ccc      print *,'REELCA:',KVAL0, '=>', KVAL1,' NUDCNB=',NUDCNB(KVAL1)
C
 9000 RETURN
      END
