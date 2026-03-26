SUBROUTINE MAGCVX( NTDL, NPDLFX, LPLIGN, LPCOLO, AG, X,  Y )
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! BUT:     Y=AG*X avec Y,X VECTEURS (NTDL)  AG MATRICE MORSE SYMETRIQUE
! ----     VARIANTE DE magcve.f AVEC CHAQUE LIGNE D'UN DEGRE DE LIBERTE
!          CONSIDEREE COMME LA LIGNE DE LA MATRICE UNITE MAIS STOCKEE
!          SANS MODIFICATION POUR PRENDRE EN COMPTE LE STOCKAGE
!          SYMETRIQUE DE LA MATRICE
!
!
! ENTREES:
! --------
! NTDL   : NOMBRE DE LIGNES DE LA MATRICE ET DU VECTEUR
! NPDLFX : NPDLFX(N)=0 SI LE DL N EST LIBRE, >0 SI LE DL N EST FIXE
!          NO TEMOIN D'UN DEGRE DE LIBERTE FIXE'
! LPLIGN : LES POINTEURS SUR LES LIGNES DE LA MATRICE MORSE AG
! LPCOLO : LES NUMEROS DES COLONNES DES COEFFICIENTS DE LA MATRICE MORSE
! AG     : LES COEFFICIENTS DE LA MATRICE MORSE AG
! X      : VECTEUR DE NTDL COMPOSANTES
!
! SORTIE ou MODIFIE:
! ------------------
! Y      : VECTEUR DE NTDL COMPOSANTES
!          ATTENTION: Y DOIT ETRE DIFFERENT DE X A L'APPEL
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! AUTEUR : Alain PERRONNET LJLL UPMC & St Pierre du Perray  Octobre 2012
!23456---------------------------------------------------------------012
!$ USE OMP_LIB
   IMPLICIT  NONE
   include"./incl/threads.inc"
   INTEGER           NTDL
   DOUBLE PRECISION  AG(1:*), X(NTDL), Y(NTDL), AJ, XI, YI, P
   INTEGER           NPDLFX(NTDL), LPLIGN(0:NTDL), LPCOLO(1:*)
   INTEGER           I, J, NODIAG0, NODIAG, NOCOL

!///////////////////////////////////////////////////////////////////////
!$OMP PARALLEL PRIVATE( I, J, NODIAG0, NODIAG, NOCOL, AJ, XI, YI, P ) &
!$OMP          SHARED ( NTDL, NPDLFX, LPLIGN, LPCOLO, AG, X,  Y )
!$OMP DO SCHEDULE( STATIC, NTDL/NBTHREADS )
!     INITIALISATION DE TOUS LES Y(I) OBLIGATOIRE AVANT SOMMATION
      DO I = 1, NTDL
         IF( NPDLFX(I) .EQ. 0 ) THEN
!           LIGNE D'UN DEGRE DE LIBERTE LIBRE
            Y(I) = 0D0
         ELSE
!           LA LIGNE I DU DL FIXE EST EN FAIT LA LIGNE IDENTITE
            Y(I) = X(I)
         ENDIF
      ENDDO
!$OMP END DO

!$OMP DO SCHEDULE( STATIC, NTDL/NBTHREADS )
      DO I = 1, NTDL

!        NUMERO +1 DU COEFFICIENT DIAGONAL I-1
         NODIAG0 = LPLIGN(I-1) + 1

!        NUMERO DU COEFFICIENT DIAGONAL I
         NODIAG = LPLIGN(I)

!        I-EME COEFFICIENT DE X
         XI = X(I)

!        I-EME COEFFICIENT DE Y
         YI = Y(I)

         IF( NPDLFX(I) .EQ. 0 ) THEN

!           LIGNE D'UN DEGRE DE LIBERTE LIBRE
!           PARCOURS DES COLONNES AVANT LA DIAGONALE DE LA LIGNE I DE AG
            DO J = NODIAG0, NODIAG-1

!              NO DE LA COLONNE DU COEFFICIENT J
               NOCOL = LPCOLO(J)

!              CONTRIBUTION DU COEFFICIENT AG(J) AVANT LA DIAGONALE
               AJ = AG(J)
               P  = AJ * X(NOCOL)
               YI = YI + P

!              CONTRIBUTION DU COEFFICIENT APRES LA DIAGONALE PAR SYMETRIE
               IF( NPDLFX(NOCOL) .EQ. 0 ) THEN
!                 COLONNE ET LIGNE D'UN DEGRE DE LIBERTE LIBRE
!$OMP ATOMIC
                  Y(NOCOL) = Y(NOCOL) + AJ * XI
               ENDIF
!
            ENDDO

!           CONTRIBUTION DU COEFFICIENT DIAGONAL SUPPOSE 1 SI DL FIXE
            P  = AG(NODIAG) * XI
            YI = YI + P

         ELSE

!           LIGNE D'UN DEGRE DE LIBERTE FIXE
            DO J = NODIAG0, NODIAG-1

!              NO DE LA COLONNE DU COEFFICIENT J
               NOCOL = LPCOLO(J)

               IF( NPDLFX(NOCOL) .EQ. 0 ) THEN
!                 COLONNE D'UN DEGRE DE LIBERTE LIBRE
!$OMP ATOMIC
                  Y(NOCOL) = Y(NOCOL) + AG(J) * XI
                ENDIF

            ENDDO

!           LA LIGNE I DU DL FIXE EST EN FAIT LA LIGNE IDENTITE
            YI = XI

         ENDIF

!        RENDEZ VOUS APRES SOMME DES COEFFICIENTS SOUS-DIAGONAUX DANS YI
         Y(I) = YI

      ENDDO
!$OMP END DO
!$OMP END PARALLEL
!///////////////////////////////////////////////////////////////////////

      RETURN
END SUBROUTINE MAGCVX
