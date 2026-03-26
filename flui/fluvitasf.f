      SUBROUTINE FLUVITASF( MNXYZP, MNELE,  NODDLN,
     %                      NCAS0,  NCAS1,  vitx, vity, vitz,
     %                      NUMISF, NUMASF, TIMES,  FLUVSF )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LE FLUX DE LA VITESSE A TRAVERS LES FACES 3D
C -----    FRONTALIERES DECLAREES DANS LA DEFINITION DE L'OBJET 3D
C
C ENTREES:
C --------
C MNXYZP : ADRESSE MCN DU TABLEAU XYZPOINT DE L'OBJET KNOMOB
C MNELE  : ADRESSE DU TABLEAU NPEF"TYPE EF DE LA TETREDRISATION
C NODDLN : NUMERO DU DERNIER DL DE CHAQUE NOEUD

C NCAS0  : NUMERO DU PREMIER VECTEUR VITESSE A TRAITER
C NCAS1  : NUMERO DU DERNIER VECTEUR VITESSE A TRAITER
C NBNOVI : NOMBRE DE NOEUDS SUPPORT DES DL DE LA VITESSE
C         (POUR TAYLOR HOOD   CE SONT LES SOMMETS ET MILIEUX DES ARETES)
C         (POUR BREZZI FORTIN CE SONT LES SOMMETS ET BARYCENTRES DES EF
C          ET IMPLICITEMENT LES SOMMETS SONT NUMEROTES DE 1 A NBSOM=NBNOPR
C          LES NOEUDS BARYCENTRE SONT NUMEROTES NBNOPR+NO EF)
C vitx   : TABLEAU(ncas0:ncas1) de pointeurs sur le tableau vitesse en X(NBNOVI)
C vity   : TABLEAU(ncas0:ncas1) de pointeurs sur le tableau vitesse en Y(NBNOVI)
C vitz   : TABLEAU(ncas0:ncas1) de pointeurs sur le tableau vitesse en Z(NBNOVI)

C NUMISF : NUMERO MINIMAL DES SURFACES DANS LA DEFINITION DE L'OBJET
C NUMASF : NUMERO MAXIMAL DES SURFACES DANS LA DEFINITION DE L'OBJET
C TIMES  : TEMPS DU CALCUL POUR LES VECTEURS STOCKES DANS DES FICHIERS
C
C SORTIES:
C --------
C FLUVSF : 1,NUSF,NCAS FLUX NEGATIF DE LA VITESSE NCAS A TRAVERS LA SURFACE NUSF
C        : 2,NUSF,NCAS FLUX POSITIF DE LA VITESSE NCAS A TRAVERS LA SURFACE NUSF
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St Pierre du Perray    Avril 2012
C MODIFS : ALAIN PERRONNET Veulettes & St Pierre du Perray     Aout 2020
C MODIFS : ALAIN PERRONNET             St Pierre du Perray  Fevrier 2022
C23456---------------------------------------------------------------012
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

      type typ_dptab
         DOUBLE PRECISION, dimension(:), pointer :: dptab
      end type typ_dptab
      type( typ_dptab ), dimension(:), allocatable :: vitx, vity, vitz

      DOUBLE PRECISION  DGL(2,3), DN(3), FLUPOS, FLUNEG, FLUXV, VN
      DOUBLE PRECISION  FLUVSF(2, NUMISF:NUMASF, NCAS0:NCAS1)
      REAL              COORDP(30), TIMES(NCAS0:NCAS1)
      CHARACTER*24      KNOMSF

C     NOMBRE D'ELEMENTS FINIS DE LA TETRAEDRISATION
      NBELEM = MCN( MNELE + WBELEM )
C
C     LE NUMERO DU TYPE DE L'ELEMENT FINI
      NUTYEL = MCN( MNELE + WUTYEL )
C
C     MNPGEL ADRESSE MCN DES NUMEROS NOEUDS ET POINTS GEOMETRIQUES DES EF
      MNPGEL = MNELE + WUNDEL
C
C     COEFFICIENTS DE L'INTEGRALE DES POLYNOMES SUR L'EF DE REFERENCE
      IF( NUTYEL .EQ. 19 ) THEN
C        TETRAEDRE BREZZI-FORTIN
         LEDECA = 0
      ELSE
C        TETRAEDRE TAYLOR-HOOD
         LEDECA = 3
      ENDIF
C
C     INITIALISATION DES FLUX A ZERO
      FLUNEG = 0D0
      FLUPOS = 0D0
      DO NCAS = NCAS0, NCAS1
         DO NUSF = NUMISF, NUMASF
            FLUVSF( 1, NUSF, NCAS ) = 0D0
            FLUVSF( 2, NUSF, NCAS ) = 0D0
         ENDDO
      ENDDO
C
C     LA BOUCLE SUR LES ELEMENTS FINIS POUR CALCULER LE FLUX A LA FRONTIERE
C     =====================================================================
      DO 100 NUELEM = 1, NBELEM
C
C        NO DES NOEUDS DE L'ELEMENT FINI NUELEM
         CALL EFNOEU( MNELE, NUELEM, NBNOEF, NONOEF )
C
C        NO DES POINTS LIGNES SURFACES VOLUMES DES SOMMETS ARETES FACES VOLUME
         CALL EFPLSV( MNELE , NUELEM,
     %                NOVCEL, NOSFEL, NOLAEL, NOPSEL,
     %                NOOBVC, NOOBSF, NOOBLA, NOOBPS )
C
C        UNE FACE DU TETRAEDRE EST ELLE FRONTALIERE?
         DO KF = 1, 4
            IF( NOOBSF(KF) .GT. 0 ) GOTO 10
         ENDDO
         GOTO 100
C
C        INTERPOLATION DU FLUIDE avec F: e reference->e  P1 ndim
C        COORDONNEES DES POINTS=NOEUDS DE L'ELEMENT FINI
C        INTEGRATION NUMERIQUE AUX 3 SOMMETS DES FACES DES SURFACES 3D
 10      CALL EFXYZP( 3, MNXYZP, NBELEM, NUELEM, MNPGEL, NBNOEF,
     %                COORDP )
         NBNOE2 = 2 * NBNOEF
C
C        CALCUL DU VECTEUR NORMAL A LA FACE K
         DO K = KF, 4
C
            NUSF = NOOBSF(K)
            IF( NUSF .GT. 0 ) THEN
C              RECHERCHE DU NUMERO DES POINTS=NOEUDS DE LA FACE K
               CALL ELNOFA( NUTYEL, K, NBNOFK, NONOFK )
C              NONOFK(I) LES NUMEROS DES NOEUDS DE LA FACE K
C                        CE SONT AUSSI LES POINTS D'INTEGRATION
C
C              RECHERCHE DU JACOBIEN DE G
               N1 = NONOFK(1)
               N2 = NONOFK(2)
               N3 = NONOFK(3)
               DGL(1,1) = COORDP( N2 )        - COORDP( N1 )
               DGL(2,1) = COORDP( N3 )        - COORDP( N1 )
               DGL(1,2) = COORDP( N2+NBNOEF ) - COORDP( N1+NBNOEF )
               DGL(2,2) = COORDP( N3+NBNOEF ) - COORDP( N1+NBNOEF )
               DGL(1,3) = COORDP( N2+NBNOE2 ) - COORDP( N1+NBNOE2 )
               DGL(2,3) = COORDP( N3+NBNOE2 ) - COORDP( N1+NBNOE2 )
C
C              CALCUL DE LA NORMALE * DELTAK: PRODUIT VECTORIEL(DG/DX1,DG/DX2)
               DN(1) = DGL(1,2) * DGL(2,3) - DGL(2,2) * DGL(1,3)
               DN(2) = DGL(1,3) * DGL(2,1) - DGL(2,3) * DGL(1,1)
               DN(3) = DGL(1,1) * DGL(2,2) - DGL(2,1) * DGL(1,2)
C
C              FLUX NORMAL 3D DE LA VITESSE A TRAVERS LA FACE K
               DO NCAS = NCAS0, NCAS1
C
                  FLUXV = 0D0
                  DO I = 1, 3
C                    NUMERO DU DL AVANT CEUX DE LA VITESSE AU NOEUD I
                     NODL = NODDLN( NONOEF( NONOFK(LEDECA+I) ) - 1 ) + 1
C                    Vitesse Normale  DELTAK  au NOEUD I de la FACE K
                     VN = vitx( NCAS )%dptab( NODL ) * DN(1)
     %                  + vity( NCAS )%dptab( NODL ) * DN(2)
     %                  + vitz( NCAS )%dptab( NODL ) * DN(3)
                     FLUXV = FLUXV + VN
                  ENDDO
                  FLUXV = FLUXV / 6D0
C
C                 REPARTITION EN FLUX NEGATIF ET POSITIF ET POIDS
                  IF( FLUXV .LE. 0D0 ) THEN
                     FLUNEG = FLUNEG + FLUXV
                  ELSE
                     FLUPOS = FLUPOS + FLUXV
                  ENDIF
C
C                 FLUX DE LA VITESSE SUR LA SURFACE NOOBSF(K)
                  FLUVSF( 1, NUSF, NCAS )=FLUVSF( 1, NUSF, NCAS )+FLUNEG
                  FLUVSF( 2, NUSF, NCAS )=FLUVSF( 2, NUSF, NCAS )+FLUPOS
C
               ENDDO
C
            ENDIF
C
         ENDDO
 100  CONTINUE

C     AFFICHAGE DU FLUX DE LA VITESSE SUR LES SURFACES DE L'OBJET
C     ===========================================================
      WRITE(IMPRIM,*)
      DO NCAS = NCAS0, NCAS1

         WRITE(IMPRIM,*)
         DO NUSF = NUMISF, NUMASF
C
C           LE FLUX DE LA VITESSE SUR LA SURFACE NUSF A-T-IK ETE CALCULE?
            FLUNEG = FLUVSF(1,NUSF,NCAS)
            FLUPOS = FLUVSF(2,NUSF,NCAS)
            FLUXV  = FLUNEG + FLUPOS
            IF( ABS(FLUNEG)+FLUPOS .GT. 0D0 ) THEN
C
C              OUI: IMPRESSION DE L'INTEGRALE
C              LE NOM DE LA SURFACE
               CALL NMOBNU( 'SURFACE', NUSF, KNOMSF )
               NN = NUDCNB( KNOMSF )
C
C              AFFICHAGE DU FLUX POSITIF, NEGATIF et TOTAL DE LA VITESSE
C              SUR LA SURFACE NUSF DE L'OBJET
               FLUXV  = FLUNEG + FLUPOS
               IF( LANGAG .EQ. 0 ) THEN
                  WRITE(IMPRIM,10100) NCAS, TIMES(NCAS), FLUNEG, FLUPOS,
     %                                FLUXV, KNOMSF(1:NN)
               ELSE
                  WRITE(IMPRIM,20100) NCAS, TIMES(NCAS), FLUNEG, FLUPOS,
     %                                FLUXV, KNOMSF(1:NN)
               ENDIF
C
            ENDIF
C
         ENDDO

      ENDDO
10100 FORMAT('Vecteur',I4,' au temps',G14.6,' FLUX<0 de la VITESSE=',
     %G14.6,' FLUX>0=',G14.6,'  FLUX TOTAL=',G14.6,' sur la SURFACE ',A)
20100 FORMAT('Vector',I4,' at time',G14.6, ' VELOCITY NEGATIVE FLUX=',
     %G14.6,' POSITIVE FLUX=',G14.6,'  TOTAL FLUX=',G14.6,' on SURFACE '
     %,A)

      RETURN
      END
