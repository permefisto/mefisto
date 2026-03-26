SUBROUTINE FLUXVITLG( NBNOVI, XYZNOE, &
                      NUTYEL, NBNOEF, NBELEM, NUNOEF, &
                      NUDDLN, NBLAEF, NULAEF, &
                      VXYZPN, NUMILG, NUMALG, &
                      FLUNEG, FLUPOS, FLUVLG )
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! BUT :    CALCULER LE FLUX DE LA VITESSE A TRAVERS LES ARETES
! -----    FRONTALIERES DECLAREES DANS LA DEFINITION DE L'OBJET
!
! ENTREES:
! --------
! NBNOVI : NOMBRE TOTAL DE NOEUDS VITESSE
! XYZNOE : XYZNOE(3,NBNOVI)  3 COORDONNEES DES NOEUDS DE LA TRIANGULATION
! NUTYEL : NUMERO DU TYPE DE L'EF
! NBNOEF : NOMBRE DE NOEUDS D'UN EF DE CE TYPE
! NBELEM : NOMBRE D'EF DU MAILLAGE
! NUNOEF : NUNOEF(NBELEM,NBNOEF) NO DES NOEUDS DES NBELEM EF
! NUDDLN : NUMERO DU DERNIER DL DE CHAQUE NOEUD
! NBLAEF : NOMBRE DES ARETES DE CE TYPE D'EF
! NULAEF : NUMERO DE LIGNE DES NBLAEF ARETES DES NBELEM EF
! VXYZPN : VITESSE PRESSION EN TOUS LES NOEUDS DU MAILLAGE
! NUMILG : NUMERO MINIMAL DES LIGNES DANS LA DEFINITION DE L'OBJET
! NUMALG : NUMERO MAXIMAL DES LIGNES DANS LA DEFINITION DE L'OBJET
!
! SORTIES:
! --------
! FLUNEG : FLUX NEGATIF DE LA VITESSE A TRAVERS LES ARETES FRONTALIERES
! FLUPOS : FLUX POSITIF DE LA VITESSE A TRAVERS LES ARETES FRONTALIERES
! FLUVLG : FLUX DE LA VITESSE A TRAVERS LES LIGNES DE L'OBJET
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! AUTEUR : ALAIN PERRONNET LJLL UPMC Paris &St PIERRE du PERRAY Mai 2011
!23456---------------------------------------------------------------012
!$ USE OMP_LIB
      IMPLICIT NONE
      include"./incl/threads.inc"
      include"./incl/langue.inc"
      INTEGER           LECTEU, IMPRIM, NUNITE
      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)
      REAL              XYZNOE(3,NBNOVI)
      INTEGER           NUNOEF(NBELEM,NBNOEF), NULAEF(NBLAEF,NBELEM)
      INTEGER           NUDDLN(0:NBNOVI), &
                        NBNOVI, NUTYEL, NBNOEF, NBELEM, NBLAEF, &
                        NUMILG, NUMALG
      DOUBLE PRECISION  VXYZPN(*)
      DOUBLE PRECISION  DN(2), FLUPOS, FLUNEG, FLUXV, VN
      DOUBLE PRECISION  FLUVLG(NUMILG:NUMALG)
      INTEGER           K,I,NONOAK(3),NODL,NUELEM,NUDCNB
      CHARACTER*24      KNOMLG

!     INITIALISATION DES FLUX A ZERO
      FLUNEG = 0D0
      FLUPOS = 0D0
      DO K = NUMILG, NUMALG
         FLUVLG( K ) = 0D0
      ENDDO

!     LA BOUCLE SUR LES ELEMENTS FINIS POUR CALCULER LE FLUX A LA FRONTIERE
!     =====================================================================
!///////////////////////////////////////////////////////////////////////
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE( NUELEM,K,NONOAK,DN,FLUXV,NODL,VN,I )
!$OMP DO SCHEDULE( STATIC, NBELEM/NBTHREADS ) &
!$OMP REDUCTION(+:FLUNEG) REDUCTION(+:FLUPOS)
      DO 100 NUELEM = 1, NBELEM

!        UNE ARETE DU TRIANGLE EST ELLE FRONTALIERE?
!        CALCUL DU VECTEUR NORMAL A L'ARETE K
         DO K = 1, 3

            IF( NULAEF(K,NUELEM) .GT. 0 ) THEN

!              RECHERCHE DU JACOBIEN DE G POUR L'ARETE K
!              NUMERO LOCAL DES NOEUDS DE L'ARETE
               NONOAK(1) = K
               IF( K .NE. 3 ) THEN
                  NONOAK(2) = K+1
               ELSE
                  NONOAK(2) = 1
               ENDIF
!              NO DU NOEUD MILIEU DE L'ARETE
               NONOAK(3) = NONOAK(1) + 3

!              CALCUL DE LA NORMALE * DELTAK: DG/DY, -DG/DX
               NONOAK(1) = NUNOEF( NUELEM, NONOAK(1) )
               NONOAK(2) = NUNOEF( NUELEM, NONOAK(2) )
               DN(1) = XYZNOE(2,NONOAK(2)) - XYZNOE(2,NONOAK(1))
               DN(2) = XYZNOE(1,NONOAK(1)) - XYZNOE(1,NONOAK(2))

!              FLUX NORMAL DE LA VITESSE A TRAVERS L'ARETE K
               FLUXV = 0D0
               IF( NUTYEL .EQ. 13 ) THEN
!                 TRIANGLE BREZZI-FORTIN INTEGRATION AUX 2 SOMMETS
                  DO I = 1, 2
!                    NUMERO DU DL AVANT CEUX DE LA VITESSE AU NOEUD I
                     NODL = NUDDLN( NONOAK(I) - 1 )
!                    Vitesse Normale  DELTAK  au NOEUD I de la ARETE K
                     VN = VXYZPN(NODL+1) * DN(1) &
                        + VXYZPN(NODL+2) * DN(2)
                     FLUXV = FLUXV + VN
                  ENDDO
                  FLUXV = FLUXV / 2D0
               ELSE
!                 TRIANGLE TAYLOR HOOD INTEGRATION AUX 3 NOEUDS DE L'ARETE
!                 POIDS 1/6 AUX 2 SOMMETS et 2/3 AU MILIEU DE L'ARETE
                  DO I = 1, 2
!                    NUMERO DU DL AVANT CEUX DE LA VITESSE AU NOEUD I
                     NODL = NUDDLN( NONOAK(I) - 1 )
!                    Vitesse Normale  DELTAK  au NOEUD I de la ARETE K
                     VN = VXYZPN(NODL+1) * DN(1) &
                        + VXYZPN(NODL+2) * DN(2)
                     FLUXV = FLUXV + VN/6D0
                  ENDDO
!                 LE MILIEU DE L'ARETE
                  NONOAK(3) = NUNOEF( NUELEM, NONOAK(3) )
                  NODL = NUDDLN( NONOAK(3) - 1 )
                  VN = VXYZPN(NODL+1) * DN(1) &
                     + VXYZPN(NODL+2) * DN(2)
                  FLUXV = FLUXV + VN*2D0/3D0
               ENDIF

!              REPARTITION EN FLUX NEGATIF ET POSITIF ET POIDS
               IF( FLUXV .LE. 0D0 ) THEN
                  FLUNEG = FLUNEG + FLUXV
               ELSE
                  FLUPOS = FLUPOS + FLUXV
               ENDIF

!              FLUX DE LA VITESSE SUR LA LIGNE NULAEF(K,NUELEM)
               I = NULAEF(K,NUELEM)
               FLUVLG( I ) = FLUVLG( I ) + FLUXV

            ENDIF

         ENDDO
 100  CONTINUE
!$OMP END DO
!$OMP END PARALLEL
!///////////////////////////////////////////////////////////////////////

!     AFFICHAGE DU FLUX DE LA VITESSE SUR LES LIGNES DE L'OBJET
!     =========================================================
      WRITE(IMPRIM,*)
      DO K = NUMILG, NUMALG

!        LE FLUX DE LA VITESSE SUR LA LIGNE K A-T-ELLE ETE CALCULEE
         FLUXV = FLUVLG( K )

         IF( FLUXV .NE. 0D0 ) THEN
!           LE NOM DE LA LIGNE
            CALL NMOBNU( 'LIGNE', K, KNOMLG )
            I = NUDCNB( KNOMLG )
            IF( LANGAG .EQ. 0 ) THEN
               WRITE(IMPRIM,*)'FLUX de la VITESSE=',FLUXV, &
                              ' sur la LIGNE ',KNOMLG(1:I)
            ELSE
               WRITE(IMPRIM,*)'VELOCITY FLUX=',FLUXV, &
                              ' through the LIGNE ',KNOMLG(1:I)
            ENDIF
         ENDIF

      ENDDO

!     AFFICHAGE DU TOTAL DES FLUX POSITIFS et NEGATIFS DE LA VITESSE
!     SUR LES LIGNES DE L'OBJET
      FLUXV = FLUNEG + FLUPOS
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10100) FLUNEG, FLUPOS, FLUXV
      ELSE
         WRITE(IMPRIM,20100) FLUNEG, FLUPOS, FLUXV
      ENDIF

10100 FORMAT(' FLUX<0 de la VITESSE=', &
      G14.6,'  FLUX>0=',G14.6,'  FLUX TOTAL=',G14.6)
20100 FORMAT(' VELOCITY NEGATIVE FLUX=',G14.6, &
      '  POSITIVE FLUX=',G14.6,'  TOTAL FLUX=',G14.6)
      WRITE(IMPRIM,*)

      RETURN
END SUBROUTINE FLUXVITLG
