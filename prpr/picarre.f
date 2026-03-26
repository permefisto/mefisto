      program picarre
      double precision picarre
      real vp(64)
      picarre = (datan(1D0) * 4D0) ** 2
      i = 0
      do 5 n=1,4
         do 4 m=1,4
            do 3 k=1,4
               i = i+1
               vp(i) = picarre * ( k*k + m*m + n*n )
               print *, k,m,n, vp(i)
 3          continue
 4       continue
 5    continue
      call tritas( 64, vp)
      do i=1,64
         print *,i,vp(i)
      enddo
      stop
      end



      SUBROUTINE TRITAS( NB, A )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRI CROISSANT DU TABLEAU A DE NB REELS PAR LA METHODE DU TAS
C -----    METHODE DUE A WILLIAMS et FLOYD     O(N LOG N)
C ENTREES:
C --------
C NB     : NOMBRE DE TERMES DU TABLEAU A
C A      : LES NB REELS A TRIER DANS A
C
C SORTIES:
C --------
C A      : LES NB REELS CROISSANTS DANS A
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC PARIS     FEVRIER 1991
C ...................................................................012
      INTEGER    PERE,PER,FIL,FILS1,FILS2,FIN
      REAL       A(1:NB),AUX
C
C     FORMATION DU TAS SOUS FORME D'UN ARBRE BINAIRE
      FIN = NB + 1
C
      DO 20 PERE = NB/2,1,-1
C
C        DESCENDRE PERE JUSQU'A N DANS A  DE FACON  A RESPECTER
C        A(PERE)>A(J) POUR J FILS OU PETIT FILS DE PERE
C        C-A-D POUR TOUT J TEL QUE PERE <= E(J/2)<J<NB+1
C                                          A(J/2) >= A(J)
C                                                 >= A(J+1)
C
C        PROTECTION DU PERE
         PER = PERE
C
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
C
C           ICI FIL EST LE PLUS GRAND DES FILS
            IF( A(PER) .LT. A(FIL) ) THEN
C              PERMUTATION DE PER ET FIL
               AUX    = A(PER)
               A(PER) = A(FIL)
               A(FIL) = AUX
C              LE NOUVEAU PERE EST LE FILS PERMUTE
               PER = FIL
               GOTO 10
            ENDIF
         ENDIF
 20   CONTINUE
C
C     A CHAQUE ITERATION LA RACINE (PLUS GRANDE VALEUR ACTUELLE DE A)
C     EST MISE A SA PLACE (FIN ACTUELLE DU TABLEAU) ET PERMUTEE AVEC
C     LA VALEUR QUI OCCUPE CETTE PLACE, PUIS DESCENTE DE CETTE NOUVELLE
C     RACINE POUR RESPECTER LE FAIT QUE TOUT PERE EST PLUS GRAND QUE TOUS
C     SES FILS
C     C-A-D POUR TOUT J TEL QUE PERE <= E(J/2)<J<NB+1
C                                          A(J/2) >= A(J)
C                                                 >= A(J+1)
      DO 50 FIN=NB,2,-1
C        LA PERMUTATION PREMIER DERNIER
         AUX    = A(FIN)
         A(FIN) = A(1)
         A(1)   = AUX
C
C        DESCENDRE A(1) ENTRE 1 ET FIN
         PER = 1
C
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
C
C           ICI FIL EST LE PLUS GRAND DES FILS
            IF( A(PER) .LT. A(FIL) ) THEN
C              PERMUTATION DE PER ET FIL
               AUX    = A(PER)
               A(PER) = A(FIL)
               A(FIL) = AUX
C              LE NOUVEAU PERE EST LE FILS PERMUTE
               PER = FIL
               GOTO 30
            ENDIF
         ENDIF
 50   CONTINUE
      END
