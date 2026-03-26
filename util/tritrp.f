      SUBROUTINE TRITRP( NB, A, NOANC )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRI CROISSANT DU TABLEAU A DE NB REELS PAR LA METHODE DU TAS
C -----    METHODE DUE A WILLIAMS et FLOYD     O( N LOG N )
C          VERSION AVEC UN POINTEUR SUR UN TABLEAU DONT EST EXTRAIT A

C ENTREES:
C --------
C NB     : NOMBRE DE TERMES DU TABLEAU A
C A      : LES NB REELS A TRIER DU TABLEAU A
C NOANC  : NUMERO ANCIENNE POSITION DE L'INFORMATION (SOUVENT NOANC(I)=I)

C SORTIES:
C --------
C A      : LES NB REELS CROISSANTS TRIES DU TABLEAU A
C NOANC  : NUMERO ANCIENNE POSITION DE L'INFORMATION
C          NOANC(1)=NO ANCIENNE POSITION DE L'INFORMATION SUR A(1), ...
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC PARIS     FEVRIER 1991
C ...................................................................012
      INTEGER    NOANC(1:NB)
      INTEGER    PERE, PER, FIL, FILS1, FILS2, FIN
      REAL       A(1:NB), AUX

C     FORMATION DU TAS SOUS FORME D'UN ARBRE BINAIRE
      FIN = NB + 1

      DO PERE = NB/2,1,-1

C        DESCENDRE PERE JUSQU'A N DANS A  DE FACON  A RESPECTER
C        A(PERE)>A(J) POUR J FILS OU PETIT FILS DE PERE
C        C-A-D POUR TOUT J TEL QUE PERE <= E(J/2)<J<NB+1
C                                          A(J/2) >= A(J)
C                                                 >= A(J+1)
C        PROTECTION DU PERE
         PER = PERE

C        LE FILS 1 DU PERE
 10      FILS1 = 2 * PER
         IF( FILS1 .LT. FIN ) THEN
C           IL EXISTE UN FILS1
            FIL   = FILS1
            FILS2 = FILS1 + 1
            IF( FILS2 .LT. FIN ) THEN
C              IL EXISTE 2 FILS . SELECTION DU PLUS GRAND
               IF( A(FILS2) .GT. A(FILS1) ) FIL = FILS2
            ENDIF

C           ICI FIL EST LE PLUS GRAND DES FILS
            IF( A(PER) .LT. A(FIL) ) THEN
C              PERMUTATION DE PER ET FIL
               AUX    = A(PER)
               A(PER) = A(FIL)
               A(FIL) = AUX
C              LE POINTEUR EST AUSSI PERMUTE
               NAUX       = NOANC(PER)
               NOANC(PER) = NOANC(FIL)
               NOANC(FIL) = NAUX
C              LE NOUVEAU PERE EST LE FILS PERMUTE
               PER = FIL
               GOTO 10
            ENDIF
         ENDIF
      ENDDO

C     A CHAQUE ITERATION LA RACINE (PLUS GRANDE VALEUR ACTUELLE DE A)
C     EST MISE A SA PLACE (FIN ACTUELLE DU TABLEAU) ET PERMUTEE AVEC
C     LA VALEUR QUI OCCUPE CETTE PLACE, PUIS DESCENTE DE CETTE NOUVELLE
C     RACINE POUR RESPECTER LE FAIT QUE TOUT PERE EST PLUS GRAND QUE TOUS
C     SES FILS
C     C-A-D POUR TOUT J TEL QUE PERE <= E(J/2)<J<NB+1
C                                          A(J/2) >= A(J)
C                                                 >= A(J+1)
      DO FIN=NB,2,-1
C        LA PERMUTATION PREMIER DERNIER
         AUX    = A(FIN)
         A(FIN) = A(1)
         A(1)   = AUX
C        LE POINTEUR EST AUSSI PERMUTE
         NAUX       = NOANC(FIN)
         NOANC(FIN) = NOANC(1)
         NOANC(1)   = NAUX

C        DESCENDRE A(1) ENTRE 1 ET FIN
         PER = 1

C        LE FILS 1 DU PERE
 30      FILS1 = 2 * PER
         IF( FILS1 .LT. FIN ) THEN
C           IL EXISTE UN FILS1
            FIL   = FILS1
            FILS2 = FILS1 + 1
            IF( FILS2 .LT. FIN ) THEN
C              IL EXISTE 2 FILS . SELECTION DU PLUS GRAND
               IF( A(FILS2) .GT. A(FILS1) ) FIL = FILS2
            ENDIF

C           ICI FIL EST LE PLUS GRAND DES FILS
            IF( A(PER) .LT. A(FIL) ) THEN
C              PERMUTATION DE PER ET FIL
               AUX    = A(PER)
               A(PER) = A(FIL)
               A(FIL) = AUX
C              LE POINTEUR EST AUSSI PERMUTE
               NAUX       = NOANC(PER)
               NOANC(PER) = NOANC(FIL)
               NOANC(FIL) = NAUX
C              LE NOUVEAU PERE EST LE FILS PERMUTE
               PER = FIL
               GOTO 30
            ENDIF
         ENDIF
      ENDDO

      RETURN
      END
