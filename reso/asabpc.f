      SUBROUTINE ASABPC( NTDL,NDSM,NBDL,NODL,AES,AE,BE,NCODSA,MUDL,
     %                   AG,BG )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    ASSEMBLER LA MATRICE ELEMENTAIRE EN LA MATRICE GLOBALE
C -----    PROFIL SYMETRIQUE OU NON TOUTE EN MEMOIRE CENTRALE
C          ASSEMBLER LE SECOND MEMBRE ELEMENTAIRE EN LE VECTEUR GLOBAL
C
C ENTREES:
C --------
C NDSM   : NOMBRE DE CAS DE CHARGE
C NBDL   : NOMBRE DE   DEGRES DE LIBERTE DE L ELEMENT FINI
C NODL   : NO DES NBDL DEGRES DE LIBERTE DE L ELEMENT FINI
C AES    : MATRICE ELEMENTAIRE SYMETRIQUE (NBDL * (NBDL+1) / 2 )
C AE     : MATRICE ELEMENTAIRE  (NBDL,NBDL)
C BE     : SECOND MEMBRE ELEMENTAIRE  (NBDL)
C NCODSA : CODE DE STOCKAGE DE LA MATRICE PROFIL
C           0 : MATRICE DIAGONALE
C          -1 : MATRICE NON SYMETRIQUE
C           1 : MATRICE SYMETRIQUE
C MUDL   : TABLEAU POINTEUR SUR LES COEFFICIENTS DIAGONAUX
C
C MODIFIES:
C ---------
C AG     : MATRICE PROFIL TOUTE EN MEMOIRE CENTRALE
C BG     : SECOND MEMBRE ELEMENTAIRE TOUT EN MEMOIRE CENTRALE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : A. PERRONNET LABORATOIRE D'ANALYSE NUMERIQUE PARIS   MAI 1989
C ......................................................................
      INTEGER           NODL(NBDL),
     %                  MUDL(1:*)
      DOUBLE PRECISION  AES(1:*),
     %                  AE(NBDL,NBDL),
     %                  BE(NDSM,NBDL),
     %                  AG(1:*),
     %                  BG(NTDL,NDSM)
C
      IF( NCODSA .EQ. 0 ) THEN
C
C        LA MATRICE EST DIAGONALE
C        ========================
         DO 10 I=1,NBDL
C           LE NO GLOBAL DU I-EME DEGRE DE LIBERTE LOCAL
            IEG = NODL( I )
C           ASSEMBLAGE DE AE( IEG )
            AG( IEG ) = AG( IEG ) + AES( I )
 10      CONTINUE
C
      ELSE IF( NCODSA .GT. 0 ) THEN
C
C        LA MATRICE EST SYMETRIQUE
C        =========================
         L = 0
         DO 40 I=1,NBDL
C           LE NO GLOBAL DU I-EME DEGRE DE LIBERTE LOCAL
            IEG = NODL( I )
C
C           ASSEMBLAGE DE LA I-EME COMPOSANTE DE BE DANS BG
            DO 20 J=1,NDSM
               BG( IEG , J ) = BG( IEG , J ) + BE( J , I )
 20         CONTINUE
C
C           ASSEMBLAGE DE AES DANS AG
C           LA BOUCLE SUR LES COLONNES DE LA MATRICE ELEMENTAIRE
            DO 30 J=1,I
               JEG = NODL( J )
               IF( JEG .LT. IEG ) THEN
C                 JEG <  IEG
                  JEG = MUDL( IEG + 1 ) - IEG + JEG
               ELSE
C                 JEG >= IEG
                  JEG = MUDL( JEG + 1 ) - JEG + IEG
               ENDIF
C              ASSEMBLAGE DE AE( I , J )
               L         = L + 1
               AG( JEG ) = AG( JEG ) + AES( L )
 30         CONTINUE
 40      CONTINUE
         RETURN
C
      ELSE
C
C        LA MATRICE EST NON SYMETRIQUE
C        =============================
         DO 140 I=1,NBDL
C           LE NO GLOBAL DU I-EME DEGRE DE LIBERTE LOCAL
            IEG = NODL( I )
C
C           ASSEMBLAGE DE BE( J , I )
            DO 120 J=1,NDSM
               BG( IEG , J ) = BG( IEG , J ) + BE( J , I )
 120        CONTINUE
C
C           LA BOUCLE SUR LES COLONNES DE LA MATRICE ELEMENTAIRE
            DO 130 J=1,NBDL
               JEG = NODL( J )
               IF( JEG .LT. IEG ) THEN
C                 JEG <  IEG
                  JEG = (MUDL(IEG+1)+MUDL(IEG)+1)/2-IEG+JEG
               ELSE
C                 JEG >= IEG
                  JEG = MUDL(JEG+1)-JEG+IEG
               ENDIF
C              ASSEMBLAGE DE AE(I,J)
               AG( JEG ) = AG( JEG ) + AE( I , J )
 130        CONTINUE
 140     CONTINUE
      ENDIF
      END
