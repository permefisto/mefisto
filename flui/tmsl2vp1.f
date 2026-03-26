      SUBROUTINE TMSL2VP1( NDIM,   NTDLVP, NBVECT,
     %                     NUTYEL, NBELEM, NBNOEF, NONOEF,
     %                     NBNOPR, NBNOVI, NDDLNV,
     %                     XYZEF,  XYZ1EF,
     %                     VXYZPR, VITEXA,   PREEXA,
     %                     VOLUME, INTL2VIT, INTL2PRE,
     %                     NOFOVI, INTL2VER, NOFOPR, INTL2PER, IERR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LE VOLUME DU MAILLAGE
C -----    CALCULER LA NORME L2 DES NBVECT |VITESSE|
C          CALCULER LA NORME L2 DES NBVECT PRESSION
C
C NDIM   : DIMENSION DE L'ESPACE DES COORDONNEES 2 ou 3
C NTDLVP : NOMBRE TOTAL DE DEGRES DE LIBERTE DU MAILLAGE DE L'OBJET
C NBVECT : NOMBRE DE VECTEURS SOLUTION VITESSEPRESSION
C NUTYEL : NUMERO DU TYPE DE L'EF
C NBELEM : NOMBRE D'EF DU MAILLAGE
C NBNOEF : NOMBRE DE NOEUDS D'UN EF
C NONOEF : NONOEF(NBELEM,NBNOEF) NO DES NOEUDS DES EF DU MAILLAGE
C NBNOVI : NOMBRE DE NOEUDS SUPPORT DES PRESSIONS DU MAILLAGE
C NBNOVI : NOMBRE DE NOEUDS SUPPORT DES VITESSES DU MAILLAGE
C NDDLNV : TABLEAU DES POINTEURS SUR LE DERNIER DL POUR CHAQUE NOEUD VITESSE
C          NDDLNV(I)= NO DERNIER DL DU NOEUD (SOMMET ou MILIEU ou BARYCENTRE)
C          CE TABLEAU EST DIMENSIONNE A 0:NBNOVI
C XYZEF  : TABLEAU XYZEF (3,NBNOVI) DES COORDONNEES DES NBNOVI NOEUDS
C XYZ1EF : TABLEAU XYZ1EF(NBNOEF,NDIM) DES COORDONNEES DES NOEUDS D'UN EF
C
C VXYZPR : VXYZPR(NTDLVP,NBVECT) LES NBVECT VECTEURS SOLUTION
C          AVEC DES DL NOEUD PAR NOEUD  NOEUD:VX,VY,VZ + PRESSION SI SOMMET
C VITEXA : VITEXA(NBNOVI,NBVECT,NDIM) VITESSE EXACTE EN CHAQUE NOEUD TEMPS
C PREEXA : PREEXA(NBNOPR,NBVECT) PRESSION EXACTE     EN CHAQUE NOEUD TEMPS
C
C NOFOVI : NUMERO DE LA FONCTION VITESSE_EXACTE(t,x,y,z)
C          0 SI ELLE N'EST PAS DONNEE PAR L'UTILISATEUR
C NOFOPR : NUMERO DE LA FONCTION PRESSION_EXACTE(t,x,y,z)
C          0 SI ELLE N'EST PAS DONNEE PAR L'UTILISATEUR
C
C SORTIES:
C --------
C VOLUME : VOLUME DU MAILLAGE
C          Volume= Som  Delta(e)
C                  e dans E
C                                                          (Vk1C)
C INTL2VIT:Som Jacobien Som (Vk1C,...,VknC) [int Pi Pj DX] ( ...)
C          e de E       k=1,...,d                          (VknC)
C                                                    (P1C)
C INTL2PRE:Som Jacobien (P1C,...,PmC) [int Pi Pj DX] (...)
C          e de E                                    (PmC)
C
C                                                                    (Vk1E-Vk1C)
C INTL2VER:Som Jacobien Som (Vk1E-Vk1C,...,VknE-VknC) [int Pi Pj DX] (    ...  )
C          e de E       k=1,...,d                                    (VknE-VknC)
C          NON CALCULE SI PRESSION_EXACTE(t,x,y,z) NON DONNEE PAR L'UTILISATEUR
C
C                                                            (P1E-P1C)
C INTL2PER:Som Jacobien (P1E-P1C,...,PmE-PmC) [int Pi Pj DX] (  ...  )
C          e de E                                            (PmE-PmC)
C
C          NON CALCULE SI VITESSE_EXACTE(t,x,y,z,nc) NON DONNEE PAR L'UTILISATEU
C IERR   : 0 SI PAS D'ERREUR, 1 SI EF NON PROGRAMME, 2 SI EF DEGENERE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET TEXAS A & M UNIVERSITY at QATAR   Fevrier 2012
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/a___npef.inc"
      include"./incl/a___xyzpoint.inc"
      include"./incl/a___vecteur.inc"
      include"./incl/ponoel.inc"
      include"./incl/gsmenu.inc"
C
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1),DMCN(1))
      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)
C
      REAL              XYZEF(3,NBNOVI), XYZ1EF(NBNOEF,3)
      DOUBLE PRECISION  VXYZPR(NTDLVP,NBVECT),
     %                  VITEXA(NBNOVI,NBVECT,NDIM),
     %                  PREEXA(NBNOPR,NBVECT),
     %                  INTL2PRE(NBVECT), INTL2VIT(NBVECT),
     %                  INTL2PER(NBVECT), INTL2VER(NBVECT),
     %                  VOLUME, VOLEF, PR, PREMIN
      INTEGER           NONOEF(NBELEM,NBNOEF), NDDLNV(0:NBNOVI),
     %                  NONO1EF(10)
C
      IERR   = 0
      MOREE2 = MOTVAR(6)
C
C     CALCUL DE LA PRESSION MINIMALE DE CHACUN DES NBVECT VECTEURS
      MNPRMI = 0
      CALL TNMCDC( 'REEL2', NBVECT, MNPRMI )
      MNDPRMI = (MNPRMI-1) / MOREE2
C
      DO NV=1, NBVECT
C
C        CALCUL DE LA PRESSION MINIMALE DU VECTEUR NV
         PREMIN = 1D100
         NDLI0  = 0
         DO I=1,NBNOVI
C           NUMERO DU DL DE PRESSION AU NOEUD I
            NDLI = NDDLNV( I )
            IF( NDLI-NDLI0 .GT. NDIM ) THEN
                PR = VXYZPR(NDLI,NV)
                IF( PR .LT. PREMIN ) PREMIN=PR
             ENDIF
             NDLI0 = NDLI
         ENDDO
         DMCN( MNDPRMI + NV ) = PREMIN
C
      ENDDO
C
C     CALCUL DU VOLUME DU MAILLAGE
      VOLUME = 0D0
C
      IF( NUTYEL .EQ. 13 .OR. NUTYEL .EQ. 19 .OR.
     %    NUTYEL .EQ. 15 .OR. NUTYEL .EQ. 20 ) GOTO 5
C
C     ERREUR: TYPE D'EF FLUIDE INCONNU
 1    IERR = 1
      GOTO 9999
C
C     LA BOUCLE SUR LES ELEMENTS FINIS DU MAILLAGE DE CE TYPE NUTYEL
C     ==============================================================
 5    DO 100 NUELEM = 1, NBELEM
C
C        LE NUMERO DES NBNOEF NOEUDS DE L'ELEMENT FINI NUELEM
         DO I=1, NBNOEF
C
C           NUMERO DU NOEUD I DE L'EF NUELEM
            NOEUD = NONOEF( NUELEM, I )
            NONO1EF( I ) = NOEUD
C
C           LES NDIM COORDONNEES DU NOEUD
            DO K=1,NDIM
               XYZ1EF( I, K ) = XYZEF( K, NOEUD )
            ENDDO
C
         ENDDO
C
C        CALCUL SELON LE TYPE DE l'EF
         GOTO(   1,   1,  1,  1,  1,   1,  1,  1,   1,  1,
     %           1,   1, 13,  1, 15,   1,  1,  1,  19, 20 ),NUTYEL
         GOTO    1
C
C        *************************
C        2D TRIANGLE BREZZI-FORTIN
C        *************************
C        CALCUL DU VOLUME ET CARRE DE LA NORME L2 DE LA PRESSION P1
 13      NONO1EF( 4 ) = NBNOPR + NUELEM
         CALL FL2PRE2( NBNOEF, XYZ1EF, NONO1EF,
     %                 NTDLVP, NBVECT, VXYZPR, MCN(MNPRMI), NDDLNV,
     %                 VOLEF,  INTL2PRE,
     %                 NOFOPR, NBNOPR, PREEXA,  INTL2PER )
C
C        CARRE DE LA NORME L2 DE |VITESSE|
         CALL FL2VITBF2( NBNOEF,   XYZ1EF, VOLEF,  NONO1EF,
     %                   NTDLVP,   NBVECT, VXYZPR, NDDLNV,
     %                   INTL2VIT, NOFOVI, NBNOVI, VITEXA, INTL2VER )
         GOTO 90
C
C        ***********************
C        2D TRIANGLE TAYLOR-HOOD
C        ***********************
C        CALCUL DU VOLUME ET CARRE DE LA NORME L2 DE LA PRESSION P1
 15      CALL FL2PRE2( NBNOEF, XYZ1EF, NONO1EF,
     %                 NTDLVP, NBVECT, VXYZPR, MCN(MNPRMI), NDDLNV,
     %                 VOLEF,  INTL2PRE,
     %                 NOFOPR, NBNOPR, PREEXA,  INTL2PER )
C
C        CARRE DE LA NORME L2 DE LA VITESSE INTEGREE EXACTEMENT
         CALL FL2VITTH2( NBNOEF, VOLEF,  NONO1EF,
     %                   NTDLVP, NBVECT, VXYZPR, NDDLNV,
     %                   INTL2VIT, NOFOVI, NBNOVI, VITEXA, INTL2VER )
         GOTO 90
C
C        **************************
C        3D TETRAEDRE BREZZI-FORTIN
C        **************************
C        CALCUL DU VOLUME ET CARRE DE LA NORME L2 DE LA PRESSION P1
 19      NONO1EF( 5 ) = NBNOPR + NUELEM
         CALL FL2PRE3( NBNOEF, XYZ1EF, NONO1EF,
     %                 NTDLVP, NBVECT, VXYZPR, MCN(MNPRMI), NDDLNV,
     %                 VOLEF,  INTL2PRE,
     %                 NOFOPR, NBNOPR, PREEXA,  INTL2PER )
C
C        CARRE DE LA NORME L2 DE LA VITESSE INTEGREE EXACTEMENT
         CALL FL2VITBF3( NBNOEF, XYZ1EF, VOLEF,  NONO1EF,
     %                   NTDLVP, NBVECT, VXYZPR, NDDLNV,
     %                   INTL2VIT, NOFOVI, NBNOVI, VITEXA, INTL2VER )
         GOTO 90
C
C        ************************
C        3D TETRAEDRE TAYLOR-HOOD
C        ************************
C        CALCUL DU VOLUME ET CARRE DE LA NORME L2 DE LA PRESSION P1
 20      CALL FL2PRE3( NBNOEF, XYZ1EF, NONO1EF,
     %                 NTDLVP, NBVECT, VXYZPR, MCN(MNPRMI), NDDLNV,
     %                 VOLEF,  INTL2PRE,
     %                 NOFOPR, NBNOPR, PREEXA,  INTL2PER )
C
C        CARRE DE LA NORME L2 DE LA VITESSE INTEGREE EXACTEMENT
         CALL FL2VITTH3( NBNOEF, VOLEF,  NONO1EF,
     %                   NTDLVP, NBVECT, VXYZPR, NDDLNV,
     %                   INTL2VIT, NOFOVI, NBNOVI, VITEXA, INTL2VER )
C
C        CONTRIBUTION DE L'EF AU VOLUME ou SURFACE DU DOMAINE
C        ----------------------------------------------------
 90      VOLUME = VOLUME + VOLEF
C
 100  CONTINUE
C
C     LE CARRE DE LA NORME L2 DEVIENT LA NORME L2
      DO K=1, NBVECT
         INTL2PRE(K) = SQRT( INTL2PRE(K) )
         INTL2VIT(K) = SQRT( INTL2VIT(K) )
         IF( NOFOPR .GT. 0 ) INTL2PER(K) = SQRT( INTL2PER(K) )
         IF( NOFOVI .GT. 0 ) INTL2VER(K) = SQRT( INTL2VER(K) )
      ENDDO
C
 9999 CALL TNMCDS( 'REEL2', NBVECT, MNPRMI )
      RETURN
      END
