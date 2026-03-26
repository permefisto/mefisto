      SUBROUTINE ASMEGC( NBDLEF, NODLEF, NCODAE, AES, AE,
     %                   NCODAG, LPLIGN, LPCOLO, AG )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   ASSEMBLER LA MATRICE ELEMENTAIRE EN LA MATRICE GLOBALE
C -----   MORSE SYMETRIQUE OU NON TOUTE EN MEMOIRE CENTRALE
C         (DIFFERENCE AVEC LE SP ASABGC PAS DE VECTEUR ELEMENTAIRE)

C ENTREES:
C --------
C NBDLEF : NOMBRE DE     DEGRES DE LIBERTE DE L ELEMENT FINI
C NODLEF : NO DES NBDLEF DEGRES DE LIBERTE DE L ELEMENT FINI
C NCODAE : CODE DE STOCKAGE DE LA MATRICE ELEMENTAIRE
C           0 : MATRICE DIAGONALE
C          -1 : MATRICE NON SYMETRIQUE
C           1 : MATRICE SYMETRIQUE
C AES    : MATRICE ELEMENTAIRE SYMETRIQUE (NBDLEF * (NBDLEF+1) / 2 )
C AE     : MATRICE ELEMENTAIRE  (NBDLEF,NBDLEF)
C NCODAG : CODE DE STOCKAGE DE LA MATRICE GLOBALE MORSE
C           0 : MATRICE DIAGONALE
C          -1 : MATRICE NON SYMETRIQUE
C           1 : MATRICE SYMETRIQUE
C LPLIGN : TABLEAU POINTEUR SUR LES COEFFICIENTS DIAGONAUX
C LPCOLO : TABLEAU DES NUMEROS DE COLONNES

C SORTIES:
C --------
C AG     : MATRICE MORSE TOUTE EN MEMOIRE CENTRALE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : A. PERRONNET LABORATOIRE ANALYSE NUMERIQUE PARIS     Mai 1989
C MODIFS : A. PERRONNET LJLL UPMC & St Pierre du Perray    Decembre 2012
C MODIFS : A. PERRONNET             St Pierre du Perray       Avril 2023
C ......................................................................
      INTEGER           NODLEF(1:NBDLEF), LPLIGN(0:*), LPCOLO(1:*)
      DOUBLE PRECISION  AES(1:*), AE(1:NBDLEF,1:NBDLEF), AG(1:*)

      IF( NCODAE .EQ. 0 .AND. NCODAG .EQ. 0 ) THEN

C        LA MATRICE ELEMENTAIRE ET GLOBALE SONT DIAGONALES
C        ====================================================
         DO I=1,NBDLEF

C           LE NO GLOBAL DU I-EME DEGRE DE LIBERTE LOCAL
            IGL = NODLEF( I )
C           ASSEMBLAGE DE AE( IGL )
            AG( IGL ) = AG( IGL ) + AES( I )

         ENDDO
         GOTO 9999

      ELSE IF( NCODAE .GT. 0 .AND. NCODAG .GT. 0 ) THEN

C        LA MATRICE ELEMENTAIRE ET GLOBALE SONT SYMETRIQUES
C        ASSEMBLAGE DE AES ELEMENTAIRE SYMETRIQUE DANS AG GLOBALE SYMETRIQUE
C        ===================================================================
C        LA BOUCLE SUR LES LIGNES DE LA MATRICE ELEMENTAIRE
         L = 0
         DO 40 I=1,NBDLEF

C           LE NO GLOBAL DU I-EME DEGRE DE LIBERTE LOCAL
            IGL = NODLEF( I )

C           LA BOUCLE SUR LES COLONNES DE LA MATRICE ELEMENTAIRE
            DO 30 J=1,I

C              LE NO GLOBAL DU I-EME DEGRE DE LIBERTE LOCAL
               JGL = NODLEF( J )

C              ADRESSE AES(I,J) DANS LA MATRICE MORSE AG
               IF(JGL.LE.IGL) THEN
                  NC1  = LPLIGN(IGL-1) + 1
                  NC2  = LPLIGN(IGL)
                  NCGL = JGL
               ELSE
                  NC1  = LPLIGN(JGL-1) + 1
                  NC2  = LPLIGN(JGL)
                  NCGL = IGL
               ENDIF

C              RECHERCHE BESTIALE DU NUMERO DE COLONNE NCGL
               DO NCM=NC1,NC2
                  IF( LPCOLO(NCM) .EQ. NCGL ) GOTO 29
               ENDDO

            print*,'asmegc.f: no ligne          =',IGL
            print*,'asmegc.f: no colonne        =',JGL
            print*,'asmegc.f: no colonne INCONNU=',NCGL
            print*,'asmegc.f: lpcolo voisins    =',(LPCOLO(k),k=nc1,nc2)
               GOTO 30

 29            L = L + 1
               AG( NCM ) = AG( NCM ) + AES( L )

 30         ENDDO

 40      ENDDO
         GOTO 9999

      ELSE IF( NCODAE .LT. 0 .AND. NCODAG .LT. 0 ) THEN

C        LA MATRICE ELEMENTAIRE ET GLOBALE SONT NON SYMETRIQUES
C        ======================================================
         DO I=1,NBDLEF

C           LE NO GLOBAL DU I-EME DEGRE DE LIBERTE LOCAL
            IGL = NODLEF( I )

C           ADRESSES DANS LA MATRICE MORSE
            NC1 = LPLIGN(IGL-1) + 1
            NC2 = LPLIGN(IGL)

C           LA BOUCLE SUR LES COLONNES DE LA MATRICE ELEMENTAIRE
            DO 130 J=1,NBDLEF

C              LE NO GLOBAL DU J-EME DEGRE DE LIBERTE LOCAL
               JGL = NODLEF( J )

C              RECHERCHE BESTIALE DU NUMERO DE COLONNE NCGL
               DO NCM=NC1,NC2
                  IF( LPCOLO(NCM) .EQ. JGL ) GOTO 129
               ENDDO

            print*,'asmegc.f: NO ligne          =',IGL
            print*,'asmegc.f: NO colonne INCONNU=',JGL
            print*,'asmegc.f: lpcolo voisins    =',(LPCOLO(k),k=nc1,nc2)

C              ASSEMBLAGE DE AE(I,J)
 129           AG( NCM ) = AG( NCM ) + AE( I , J )

 130        ENDDO
         ENDDO
         GOTO 9999

      ELSE IF( NCODAE .EQ. 0 .AND. NCODAG .NE. 0 ) THEN

C        LA MATRICE ELEMENTAIRE EST DIAGONALE MAIS
C        LA MATRICE GLOBALE EST SYMETRIQUE OU NON SYMETRIQUE
C        ===================================================
         DO I=1,NBDLEF

C           LE NO GLOBAL DU I-EME DEGRE DE LIBERTE LOCAL
            IGL = NODLEF( I )
C           SA POSITION SUR LA DIAGONALE DE LA MATRICE GLOBALE MORSE
            IGL = LPLIGN( IGL )
C           ASSEMBLAGE DE AES( IGL )
            AG( IGL ) = AG( IGL ) + AES( I )

         ENDDO

      ENDIF

 9999 RETURN
      END
