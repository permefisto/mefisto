SUBROUTINE MAXVIT( NDIM,   NBNOVI, NBVECT, VX, VY, VZ, &
                   NOEMIN, NCAMIN, VITMIN, &
                   NOEMAX, NCAMAX, VITMAX, VITMOY )
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! BUT :   CALCULER LA NORME MIN MAX MOYENNE DU VECTEUR VITESSE
! -----   AUX NBNOVI NOEUDS

! ENTREES :
! ---------
! NDIM   : DIMENSION DE L'ESPACE DE L'OBJET (2 OU 3)
! NBNOVI : NOMBRE DE NOEUDS SUPPORT DE LA VITESSE
! NBVECT : NOMBRE TOTAL DE VECTEURS VITESSE PRESSION
! VX     : LA VITESSE EN X EN CHAQUE NOEUD ET NBVECT FOIS
! VY     : LA VITESSE EN Y EN CHAQUE NOEUD ET NBVECT FOIS
! VZ     : LA VITESSE EN Z EN CHAQUE NOEUD ET NBVECT FOIS

! SORTIES:
! --------
! NOEMIN : NUMERO DU NOEUD   DE LA NORME DE LA VITESSE MINIMALE
! NCAMIN : NUMERO DU VECTEUR DE LA NORME DE LA VITESSE MINIMALE
! VITMIN : NORME             DE LA VITESSE MINIMALE
! NOEMAX : NUMERO DU NOEUD   DE LA NORME DE LA VITESSE MAXIMALE
! NCAMAX : NUMERO DU VECTEUR DE LA NORME DE LA VITESSE MAXIMALE
! VITMAX : NORME             DE LA VITESSE MAXIMALE
! VITMAX : NORME MOYENNE DE LA VITESSE
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE DU PERRAY Novembre 2012
!.......................................................................
!$ USE OMP_LIB
      include"./incl/threads.inc"
      DOUBLE PRECISION  VX(NBNOVI,NBVECT)
      DOUBLE PRECISION  VY(NBNOVI,NBVECT)
      DOUBLE PRECISION  VZ(NBNOVI,NBVECT)
      DOUBLE PRECISION  VITMIN, VITMAX, VITMOY, VITE
      INTRINSIC         SQRT

!     LA VITESSE MINIMALE ET MAXIMALE
!     -------------------------------
      NOEMIN = 1
      NCAMIN = 1
      NOEMAX = 1
      NCAMAX = 1
      VITMIN = VX(1,1)**2 + VY(1,1)**2
      IF( NDIM .EQ. 3 ) VITMIN = VITMIN + VZ(1,1)**2
      VITMIN = SQRT(VITMIN)
      VITMAX = VITMIN
      VITMOY = 0D0

!///////////////////////////////////////////////////////////////////////
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE( K, I, VITE )
!$OMP DO REDUCTION(+:VITMOY)
      DO I=1,NBNOVI
         DO K=1,NBVECT

!           LA NORME DE LA VITESSE AU NOEUD
            VITE = VX(I,K)**2 + VY(I,K)**2
            IF( NDIM .EQ. 3 ) VITE = VITE + VZ(I,K)**2
            VITE = SQRT(VITE)

            VITMOY = VITMOY + VITE
            IF( VITE .LT. VITMIN ) THEN
!$OMP CRITICAL
               VITMIN = VITE
               NOEMIN = I
               NCAMIN = K
!$OMP END CRITICAL
            ELSE IF( VITE .GT. VITMAX ) THEN
!$OMP CRITICAL
               VITMAX = VITE
               NOEMAX = I
               NCAMAX = K
!$OMP END CRITICAL
            ENDIF
!
         ENDDO
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
!///////////////////////////////////////////////////////////////////////

!   LA MOYENNE AUX NOEUDS DE LA NORME DE LA VITESSE
!   -----------------------------------------------
    VITMOY = VITMOY / (NBVECT*NBNOVI)

!! PRINT *,'maxvit.f95: VITMIN=',VITMIN,' VITMAX=',VITMAX,' VITMOY=',VITMOY

    RETURN
END SUBROUTINE MAXVIT
