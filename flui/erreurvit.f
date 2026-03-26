      SUBROUTINE ERREURVIT( KNOMOB, KNOMFIC, NDIM, NUTYEL, MNELE,
     %                      MNXYZN, NODDL, NBNOVI, NTDLVP, NTDLTE,
     %                      NCAS0,  NCAS1, TIMES,  VXYZPN, TEMPER,
     %                      NOFOVI, ERRVIT, VitErMOY,
     %                      VitErMIN, NOEMIN, NCASMIN,
     %                      VitErMAX, NOEMAX, NCASMAX )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LE TABLEAU DE L'ERREUR SUR LE MODULE DE LA VITESSE
C -----    AUX NOEUDS DU MAILLAGE C-A-D
C          AUX SOMMETS POUR  BREZZI-FORTIN et TAYLOR-HOOD
C          + AUX BARYCENTRES DES EF POUR BREZZI-FORTIN
C          + AUX MILIEUX DES ARETES POUR TAYLOR-HOOD

C ENTREES:
C --------
C KNOMOB : NOM DE L'OBJET CONTENANT LE FLUIDE
C KNOMFIC: NOM DU FICHIER SUPPORT D'UN VECTEUR VITESSE+PRESSION
C NDIM   : DIMENSION 2 OU 3 DE L'ESPACE DE L'OBJET
C NUTYEL : NUMERO DU TYPE D'EF UTILISE
C          13 ou 19 : BREZZI FORTIN <=> TRIA 2P1D ou TETR 3P1D + BARYCENTRE
C          15 ou 20 : TAYLOR HOOD   <=> TRIA 2P2D ou TRIA 2P2C
C MNELE  : ADRESSE MCN DU TABLEAU NPEF"TYPE EF (P1D ou P2C)

C MNXYZN : ADRESSE MCN DU TMS 'XYZNOEUD' DU MAILLAGE DE L'OBJET
C NODDL  : TABLEAU DU NUMERO DU DERNIER D.L. DE CHAQUE NOEUD
C NBNOVI : NOMBRE DE NOEUDS SUPPORT DE LA VITESSE
C          SOMMETS + BARYCENTRES        des TETRAEDRES POUR BREZZI-FORTIN
C          SOMMETS + MILIEUX DES ARETES des TETRAEDRES POUR TAYLOR-HOOD
C NTDLVP : NOMBRE TOTAL DE DEGRES DE LIBERTES EN VITESSES PRESSIONS
C NCAS0  : PREMIER VECTEUR VITESSE-PRESSION A TRAITER
C NCAS1  : DERNIER VECTEUR VITESSE-PRESSION A TRAITER
C TIMES  : LES NCAS0:NCAS1 TEMPS DES CALCULS DE LA VITESSE-PRESSION
C VXYZPN : UNE CARTE NCAS DES VITESSES-PRESSIONS AUX NOEUDS DU MAILLAGE
C NOFOVI : NUMERO DE LA FONCTION VITESSE_EXACTE(t,x,y,z,nc)
C          =0 SI ELLE N'EST PAS DONNEE PAR L'UTILISATEUR -> RETOUR

C SORTIES:
C --------
C ERRVIT   : ERREUR SUR LE MODULE DE LA VITESSE EN TOUS LES NOEUDS
C            ET EN TOUS LES CAS DE TEMPS NCAS0:NCAS1
C VitErMOY : ERREUR MOYENNE SUR |VITESSE|(Temps,Noeuds) POUR TOUS LES NOEUDS ET TEMPS
C VitErMIN : ERREUR MINIMALE SUR |VITESSE|(Temps,Noeuds)
C NOEMIN   : NUMERO DU NOEUD   DE L'ERREUR MINIMALE SUR |VITESSE|(Temps,Noeuds)
C NCASMIN  : NUMERO DU VECTEUR DE L'ERREUR MINIMALE SUR |VITESSE|(Temps,Noeuds)
C VitErMAX : ERREUR MAXIMALE SUR |VITESSE|(Temps,Noeuds)
C NOEMAX   : NUMERO DU NOEUD   DE L'ERREUR MAXIMALE SUR |VITESSE|(Temps,Noeuds)
C NCASMAX  : NUMERO DU VECTEUR DE L'ERREUR MAXIMALE SUR |VITESSE|(Temps,Noeuds)

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Saint Pierre du Perray               Mai 2021
C234567...............................................................12
      include"./incl/langue.inc"
      include"./incl/a___xyznoeud.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE      (MCN(1), RMCN(1), DMCN(1))
      CHARACTER*(*)     KNOMOB
      CHARACTER*80      KNOMFIC
      INTEGER           NODDL(0:NBNOVI), NONOEF(10)
      REAL              TIMES(NCAS0:NCAS1), VitErMIN, VitErMAX
      DOUBLE PRECISION  VXYZPN(1:NTDLVP), ERRVIT(1:NBNOVI,NCAS0:NCAS1),
     %                  TEMPER(*),
     %                  XP, YP, ZP, DPARAF(5), VitErrMOY,
     %                  VitErr, VitErrMIN, VitErrMAX,
     %                  XVitExact, YVitExact, ZVitExact
 
      INTRINSIC         SQRT

      IF( NOFOVI .LE. 0 ) RETURN

      NDIM1 = NDIM + 1
C     INITIALISATION DU MIN MAX
      NOEMIN  = 1
      NCASMIN = 1
      VitErrMIN = 1D100

      NOEMAX  = 1
      NCASMAX = 1
      VitErrMAX =-1D100

      VitErrMOY = 0D0

C     BOUCLE SUR LES TEMPS DES CALCULS
C     ================================
      MNXYZ = MNXYZN + WYZNOE - 3

      DO NCAS = NCAS0, NCAS1

C        LE TEMPS DU CALCUL
         TEMPS = TIMES(NCAS)

C        TRANSFERT DE LA VITESSE-PRESSION AU TEMPS NCAS PAR NOEUDS
C        FICHIER DU REPERTOIRE PROJET DANS LE TABLEAU MNVXYZP0
         CALL LIFIVIPRTE( KNOMOB, TEMPS,  NCAS,   NAVSTO, NBPASDT,
     %                    NDIM,   NBNOVI, NBNOPR, NBNOTE,
     %                    NTDLVP, VXYZPN, NTDLTE, TEMPER,
     %                    KNOMFIC, IERR )

C        L'ADRESSE MCN DES COORDONNEES DES NOEUDS VITESSE
C        ICI LE TMS XYZNOEUD A ETE ENRICHI DES BARYCENTRES DES EF
C        DANS LES CAS DES EF DE BREZZI-FORTIN
         MN = MNXYZ

         DO I=1,NBNOVI

            IF( (NUTYEL .EQ. 13 .OR. NUTYEL .EQ. 19) .AND.
     %           I .GT. NBNOPR ) THEN

C              NOEUD I=BARYCENTRE DE L'EF BREZZI-FORTIN P1+BARYCENTRE
               NUELEM = I - NBNOPR

C              LES NOEUDS DE L'ELEMENT FINI NUELEM
               CALL EFNOEU( MNELE, NUELEM, NBNOEF, NONOEF )

C              CALCUL DES XYZ DU BARYCENTRE DE L'EF NUELEM
               XP = 0D0
               YP = 0D0
               ZP = 0D0
               DO K=1,NDIM1
                  NS = NONOEF( K )
                  MN = MNXYZ + 3 * NS
                  XP = XP + RMCN( MN )
                  YP = YP + RMCN( MN+1 )
                  ZP = ZP + RMCN( MN+2 )
               ENDDO
               XP = XP / NDIM1
               YP = YP / NDIM1
               ZP = ZP / NDIM1

            ELSE
 
C              EF TAYLOR-HOOD: SOMMETS et MILIEUX des ARETES
C              ADRESSE DES COORDONNEES DU NOEUD I
               MN = MN + 3
               XP = RMCN( MN )
               YP = RMCN( MN+1 )
               ZP = RMCN( MN+2 )

            ENDIF

C           NOMBRE DE DL AU NOEUD I
            NDL = NODDL(I-1)
            ND  = NODDL(I) - NDL

C           PARAMETRES D'APPEL DE LA VITESSE EXACTE
            DPARAF(1) = TEMPS
            DPARAF(2) = XP
            DPARAF(3) = YP
            DPARAF(4) = ZP

C           LA VITESSE EXACTE(temps,xp,yp,zp,NoComposante)
C           ----------------------------------------------
C           NUMERO DE LA COMPOSANTE TRAITEE DE LA VITESSE
            DPARAF(5) = 1
C           COMPOSANTE X DE LA VITESSE EXACTE AU NOEUD I AU TEMPS NCAS
            CALL FONVAL( NOFOVI, 5, DPARAF, NCODEV, XVitExact )

            DPARAF(5) = 2
C           COMPOSANTE Y DE LA VITESSE EXACTE AU NOEUD I AU TEMPS NCAS
            CALL FONVAL( NOFOVI, 5, DPARAF, NCODEV, YVitExact )

            IF( NDIM .GE. 3 ) THEN
               DPARAF(5) = 3
C              COMPOSANTE Z DE LA VITESSE EXACTE AU NOEUD I AU TEMPS NCAS
               CALL FONVAL( NOFOVI, 5, DPARAF, NCODEV, ZVitExact )
            ELSE
               ZVitExact = 0D0
            ENDIF

C           CALCUL DES ERREURS LOCALES SUR LA VITESSE
            IF( NDIM .GE. 3 ) THEN
               VitErr  = SQRT( (XVitExact-VXYZPN(NDL+1))**2
     %                       + (YVitExact-VXYZPN(NDL+2))**2
     %                       + (ZVitExact-VXYZPN(NDL+3))**2 )
            ELSE
               VitErr   = SQRT( (XVitExact-VXYZPN(NDL+1))**2
     %                        + (YVitExact-VXYZPN(NDL+2))**2 )
            ENDIF

C           ERREUR SUR LA VITESSE AU NOEUD I ET AU TEMPS NCAS
            ERRVIT( I, NCAS ) = VitErr

C           MODULE DE ERREUR MIN ET MAX AU NOEUD DE LA VITESSE
            IF( VitErr .LT. VitErrMIN ) THEN
               VitErrMIN = VitErr
               NOEMIN  = I
               NCASMIN = NCAS
            ELSE IF( VitErr .GT. VitErrMAX ) THEN
               VitErrMAX = VitErr
               NOEMAX  = I
               NCASMAX = NCAS
            ENDIF

C           AJOUT DE L'ERREUR AU NOEUD I AU TEMPS NCAS
            VitErrMOY = VitErrMOY + VitErr

         ENDDO

      ENDDO

C     ERREUR MOYENNE POUR TOUS LES NOEUDS ET TEMPS
      VitErrMOY = (VitErrMOY / NBNOVI) / (NCAS1-NCAS0+1)
      IF( LANGAG .EQ. 0 ) THEN
         PRINT*,'erreurvit: ERREUR MOYENNE sur la VITESSE=',VitErrMOY
      ELSE
         PRINT*,'erreurvit: VELOCITY MEAN ERROR =',VitErrMOY
      ENDIF

C     CONVERSION EN REEL SIMPLE
      VitErMOY = REAL( VitErrMOY )
      VitErMIN = REAL( VitErrMIN )
      VitErMAX = REAL( VitErrMAX )

      RETURN
      END
