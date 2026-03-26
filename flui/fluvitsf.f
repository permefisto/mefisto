      SUBROUTINE FLUVITSF( MNXYZP, MNELE, NODDLN, VXYZPN,
     %                     NUMISF, NUMASF,
     %                     FLUNEG, FLUPOS, FLUVSF )
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! BUT :    CALCULER LE FLUX DE LA VITESSE A TRAVERS LES FACES
! -----    FRONTALIERES DECLAREES DANS LA DEFINITION DE L'OBJET
!
! ENTREES:
! --------
! MNXYZP : ADRESSE MCN DU TABLEAU XYZPOINT DE L'OBJET KNOMOB
! MNELE  : ADRESSE DU TABLEAU NPEF"TYPE EF DE LA TETREDRISATION
! NODDLN : NUMERO DU DERNIER DL DE CHAQUE NOEUD
! VXYZPN : VITESSE PRESSION EN TOUS LES NOEUDS DU MAILLAGE
! NUMISF : NUMERO MINIMAL DES SURFACES DANS LA DEFINITION DE L'OBJET
! NUMASF : NUMERO MAXIMAL DES SURFACES DANS LA DEFINITION DE L'OBJET
!
! SORTIES:
! --------
! FLUNEG : FLUX NEGATIF DE LA VITESSE A TRAVERS LES FACES FRONTALIERES
! FLUPOS : FLUX POSITIF DE LA VITESSE A TRAVERS LES FACES FRONTALIERES
! FLUVSF : FLUX DE LA VITESSE A TRAVERS CHACUNE DES SURFACES DE L'OBJET
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY     Juin 2010
!23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/a___xyzpoint.inc"
      include"./incl/a___npef.inc"
      include"./incl/ponoel.inc"
      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)
      INTEGER           MCN
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      INTEGER           NODDLN(0:*)
      INTEGER           NONOEF(10), NONOFK(6)
      INTEGER           NOOBVC, NOOBSF(6), NOOBLA(12), NOOBPS(8)
      DOUBLE PRECISION  VXYZPN(*)
      DOUBLE PRECISION  DGL(2,3), DN(3), FLUPOS, FLUNEG, FLUXV, VN
      DOUBLE PRECISION  FLUVSF(NUMISF:NUMASF)
      REAL              COORDP(30)
      CHARACTER*24      KNOMSF

!     INITIALISATION DES FLUX A ZERO
      FLUNEG = 0D0
      FLUPOS = 0D0
      DO N = NUMISF, NUMASF
         FLUVSF( N ) = 0D0
      ENDDO

!     NOMBRE D'ELEMENTS FINIS DE LA TETRAEDRISATION
      NBELEM = MCN( MNELE + WBELEM )

!     LE NUMERO DU TYPE DE L'ELEMENT FINI
      NUTYEL = MCN( MNELE + WUTYEL )

!     MNPGEL ADRESSE MCN DES NUMEROS NOEUDS ET POINTS GEOMETRIQUES DES EF
      MNPGEL = MNELE + WUNDEL

!     COEFFICIENTS DE L'INTEGRALE DES POLYNOMES SUR L'EF DE REFERENCE
      IF( NUTYEL .EQ. 19 ) THEN
!        TETRAEDRE BREZZI-FORTIN
         LEDECA = 0
      ELSE
!        TETRAEDRE TAYLOR-HOOD
         LEDECA = 3
      ENDIF

!     LA BOUCLE SUR LES ELEMENTS FINIS POUR CALCULER LE FLUX A LA FRONTIERE
!     =====================================================================
      DO 100 NUELEM = 1, NBELEM

!        NO DES NOEUDS DE L'ELEMENT FINI NUELEM
         CALL EFNOEU( MNELE, NUELEM, NBNOEF, NONOEF )

!        NO DES POINTS LIGNES SURFACES VOLUMES DES SOMMETS ARETES FACES VOLUME
         CALL EFPLSV( MNELE , NUELEM,
     %                NOVCEL, NOSFEL, NOLAEL, NOPSEL,
     %                NOOBVC, NOOBSF, NOOBLA, NOOBPS )

!        UNE FACE DU TETRAEDRE EST ELLE FRONTALIERE?
         DO KF = 1, 4
            IF( NOOBSF(KF) .GT. 0 ) GOTO 10
         ENDDO
         GOTO 100

!        INTERPOLATION DU FLUIDE avec F: e reference->e  P1 ndim
!        COORDONNEES DES POINTS=NOEUDS DE L'ELEMENT FINI
!        INTEGRATION NUMERIQUE AUX 3 SOMMETS DES FACES DES SURFACES
 10      CALL EFXYZP( 3, MNXYZP, NBELEM, NUELEM, MNPGEL, NBNOEF,
     %                COORDP )
         NBNOE2 = 2 * NBNOEF

!        CALCUL DU VECTEUR NORMAL A LA FACE K
         DO K = KF, 4

            IF( NOOBSF(K) .GT. 0 ) THEN
!              RECHERCHE DU NUMERO DES POINTS=NOEUDS DE LA FACE K
               CALL ELNOFA( NUTYEL, K, NBNOFK, NONOFK )
!              NONOFK(I) LES NUMEROS DES NOEUDS DE LA FACE K
!                        CE SONT AUSSI LES POINTS D'INTEGRATION

!              RECHERCHE DU JACOBIEN DE G
               N1 = NONOFK(1)
               N2 = NONOFK(2)
               N3 = NONOFK(3)
               DGL(1,1) = COORDP( N2 )        - COORDP( N1 )
               DGL(2,1) = COORDP( N3 )        - COORDP( N1 )
               DGL(1,2) = COORDP( N2+NBNOEF ) - COORDP( N1+NBNOEF )
               DGL(2,2) = COORDP( N3+NBNOEF ) - COORDP( N1+NBNOEF )
               DGL(1,3) = COORDP( N2+NBNOE2 ) - COORDP( N1+NBNOE2 )
               DGL(2,3) = COORDP( N3+NBNOE2 ) - COORDP( N1+NBNOE2 )

!              CALCUL DE LA NORMALE * DELTAK: PRODUIT VECTORIEL(DG/DX1,DG/DX2)
               DN(1) = DGL(1,2) * DGL(2,3) - DGL(2,2) * DGL(1,3)
               DN(2) = DGL(1,3) * DGL(2,1) - DGL(2,3) * DGL(1,1)
               DN(3) = DGL(1,1) * DGL(2,2) - DGL(2,1) * DGL(1,2)

!              FLUX NORMAL DE LA VITESSE A TRAVERS LA FACE K
               FLUXV = 0D0
               DO I = 1, 3
!                 NUMERO DU DL AVANT CEUX DE LA VITESSE AU NOEUD I
                  NODL = NODDLN( NONOEF( NONOFK(LEDECA+I) ) - 1 )
!                 Vitesse Normale  DELTAK  au NOEUD I de la FACE K
                  VN = VXYZPN(NODL+1) * DN(1)
     %               + VXYZPN(NODL+2) * DN(2)
     %               + VXYZPN(NODL+3) * DN(3)
                  FLUXV = FLUXV + VN
               ENDDO
               FLUXV = FLUXV / 6D0

!              REPARTITION EN FLUX NEGATIF ET POSITIF ET POIDS
               IF( FLUXV .LE. 0D0 ) THEN
                  FLUNEG = FLUNEG + FLUXV
               ELSE
                  FLUPOS = FLUPOS + FLUXV
               ENDIF

!              FLUX DE LA VITESSE SUR LA SURFACE NOOBSF(K)
               FLUVSF( NOOBSF(K) ) = FLUVSF( NOOBSF(K) ) + FLUXV

            ENDIF

         ENDDO
 100  CONTINUE

!     AFFICHAGE DU FLUX DE LA VITESSE SUR LES SURFACES DE L'OBJET
!     ===========================================================
      WRITE(IMPRIM,*)
      DO K = NUMISF, NUMASF
!        LE FLUX DE LA VITESSE SUR LA SURFACE K A-T-ELLE ETE CALCULEE
         FLUXV = FLUVSF( K )
         IF( FLUXV .NE. 0D0 ) THEN
!            OUI: IMPRESSION DE L'INTEGRALE
!            LE NOM DE LA SURFACE
             CALL NMOBNU( 'SURFACE', K, KNOMSF )
             NN = NUDCNB( KNOMSF )
             IF( LANGAG .EQ. 0 ) THEN
                WRITE(IMPRIM,*)'FLUX NORMAL de la VITESSE=',FLUXV,
     %                ' sur la SURFACE ',KNOMSF(1:NN)
             ELSE
                WRITE(IMPRIM,*)'VELOCITY NORMAL FLUX=',FLUXV,
     %                ' through the SURFACE ',KNOMSF(1:NN)
             ENDIF
         ENDIF
      ENDDO

!     AFFICHAGE DU TOTAL DES FLUX POSITIFS et NEGATIFS DE LA VITESSE
!     SUR LES SURFACES DE L'OBJET
      FLUXV = FLUNEG + FLUPOS
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10100) FLUNEG, FLUPOS, FLUXV
      ELSE
         WRITE(IMPRIM,20100) FLUNEG, FLUPOS, FLUXV
      ENDIF

10100 FORMAT(' Totaux sur les surfaces: FLUX<0 de la VITESSE=',
     %G14.6,'  FLUX>0=',G14.6,'  FLUX TOTAL=',G14.6)
20100 FORMAT(' Total on all Surfaces: VELOCITY NEGATIVE FLUX=',G14.6,
     %'  POSITIVE FLUX=',G14.6,'  TOTAL FLUX=',G14.6)
!!!      WRITE(IMPRIM,*)

      RETURN
      END
