SUBROUTINE VITECRET( NDIM, NBNOVI, VX, VY, VZ, VitECR )
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! BUT :     ECRETAGE DE LA VITESSE AUX NOEUDS
! -----     DES TETRAEDRES de TAYLOR HOOD ou BREZZI-FORTIN
!
! ENTREES :
! ---------
! NDIM   : DIMENSION DE L'ESPACE DE L'OBJET (2 OU 3)
! NBNOVI : NOMBRE DE NOEUDS SUPPORT DE LA VITESSE
!         (POUR TAYLOR HOOD CE SONT LES SOMMETS ET MILIEUX DES ARETES)
!         (POUR BREZZI FORTIN CE SONT LES SOMMETS ET BARYCENTRES DES EF
!          ET IMPLICITEMENT LES SOMMETS SONT NUMEROTES DE 1 A NBNOPR
!             LES AUTRES NOEUDS BARYCENTRE SONT NUMEROTES NBNOPR+NO EF)
!	   VitECR : VITESSE D'ECRETAGE

! ENTREES et SORTIES:
! -------------------
! VX     : LA VITESSE EN X EN CHAQUE NOEUD VITESSE ET NBVECT FOIS
! VY     : LA VITESSE EN Y EN CHAQUE NOEUD VITESSE ET NBVECT FOIS
! VZ     : LA VITESSE EN Z EN CHAQUE NOEUD VITESSE ET NBVECT FOIS
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! AUTEUR : ALAIN PERRONNET  Veulettes sur mer               Octobre 2020
!23456---------------------------------------------------------------012
!$ USE OMP_LIB
   IMPLICIT  NONE
   include"./incl/langue.inc"
   include"./incl/threads.inc"

   DOUBLE PRECISION  VX(NBNOVI), VY(NBNOVI), VZ(NBNOVI), &
                     VitECR, VNORM, D, SQRT
   INTEGER           NDIM, NBNOVI, I, NBECRET

   NBECRET = 0
   IF( VitECR .LE. 0D0 ) GOTO 9999

!  TAYLOR-HOOD ou BREZZI-FORTIN
!  ===========    =============
   IF( NDIM .EQ. 2 ) THEN

!///////////////////////////////////////////////////////////////////////
!$OMP PARALLEL PRIVATE( I, VNORM, D )
!$OMP DO SCHEDULE( STATIC, NBNOVI/NBTHREADS )

      DO I=1,NBNOVI

!        NORME DE LA VITESSE au TEMPS tn+1
         VNORM = SQRT( VX(I)**2 + VY(I)**2 )

         IF( VNORM .GT. VitECR ) THEN
!           ECRETAGE DE LA VITESSE AU NOEUD I
            NBECRET = NBECRET + 1
            D = VitECR / VNORM
            VX(I) = VX(I) * D
            VY(I) = VY(I) * D
         ENDIF

      ENDDO

!$OMP END DO
!$OMP END PARALLEL
!/////////////////////////////////////////////////////////////////////////

   ELSE

!///////////////////////////////////////////////////////////////////////
!$OMP PARALLEL PRIVATE( I, VNORM, D )
!$OMP DO SCHEDULE( STATIC, NBNOVI/NBTHREADS )

      DO I=1,NBNOVI

!        NORME DE LA VITESSE au TEMPS tn+1
         VNORM = SQRT( VX(I)**2 + VY(I)**2 + VZ(I)**2 )

         IF( VNORM .GT. VitECR ) THEN
!           ECRETAGE DE LA VITESSE AU NOEUD I
            NBECRET = NBECRET + 1
            D     = VitECR / VNORM
            VX(I) = VX(I) * D
            VY(I) = VY(I) * D
            VZ(I) = VZ(I) * D
         ENDIF

      ENDDO

!$OMP END DO
!$OMP END PARALLEL
!/////////////////////////////////////////////////////////////////////////

   ENDIF

 9999 IF( LANGAG .EQ. 0 ) THEN
         PRINT*,'vitecret.f95:',NBECRET,'/',NBNOVI,  &
                ' NOEUDS de VITESSE ECRETEE a',VitECR
      ELSE
         PRINT*,'vitecret.f95:',NBECRET,'/',NBNOVI,  &
                ' NODES with REDUCED VELOCITY to',VitECR
      ENDIF

      RETURN
END SUBROUTINE VITECRET
