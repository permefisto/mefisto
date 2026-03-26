      SUBROUTINE NLSEINS( KNOMOB, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCULER u(t,X) solution de l'Equation Non Lineaire de
C ----- SCHRODINGER (NLSE) avec ou non un terme de rotation
C       soit U(t,X) est SOLUTION de
C i Rho dU(t,X)/dt - Alfa LAPLACIEN U(t,X) + N(|U|**2) U(t,X)
C                  - i OmegaZ (x dU/dy - y dU/dx) = F

C Exemple: N(|U|**2) = Beta  |U|**2

C CALCUL de la PARTIE REELLE V et IMAGINAIRE W de U(t,X) = V(t,X) + i W(t,X)
C  selon le systeme avec F = Fr(t,X) + i Fi(t,X)

C -Rho dW(t,X)/dt - Alfa LAPLACIEN V(t,X) + N(V**2+W**2) V(t,X)
C                 + OmegaZ ( x dW/dy - y dW/dx ) = Fr(t,X)

C  Rho dV(t,X)/dt - Alfa LAPLACIEN W(t,X) + N(V**2+W**2) W(t,X)
C                 - OmegaZ ( x dV/dy - y dV/dx ) = Fi(t,X)

C La valeur TESTNL du common"./incl/cnonlin.inc" DEFINIT LE SCHEMA DE RESOLUTION

C TESTNL=6: SCHEMA SEMI-IMPLICITE MODIFIE en TEMPS avec Matrice GLOBALE nxn
C -------------------------------------------------------------------------
C VERSION 1:
C   n=0  V=v0      W=w0
C
C(1)m=0  Vn+1m=Vn  Wn+1m=Wn
C
C(2)[N(V0**2+W0**2) -Alfa LAPLACIEN] Vn+1m+1 = Fr(tn+1)
C   + Rho/dt (Wn+1m-Wn) - OmegaZ ( x dWn+1m/dy - y dWn+1m/dx )
C   + ( N(V0**2+W0**2)-N(Vn+1m**2+Wn+1m**2) ) Vn+1m         
C
C   [N(V0**2+W0**2) -Alfa LAPLACIEN] Wn+1m+1 = Fi(tn+1)
C   - Rho/dt (Vn+1m-Vn) + OmegaZ ( x dVn+1m/dy - y dVn+1m/dx )
C   + ( N(V0**2-W0**2)-N(Vn+1m**2+Wn+1m**2) ) Wn+1m
C                       
C    Si ||Vn+1m+1-Vn+1m||>Eps||Vn+1m+1|| ou ||Wn+1m+1-Wn+1m||>Eps||Wn+1m+1||
C    Alors m=m+1  Aller en (2)
C    Sinon Vn+1=Vn+1m+1  Wn+1=Wn+1m+1
C
C   Si TPSINI+(n+1)dt < TPSFIN Alors n=n+1; Aller en (1)
C   Fin
C
C Ajout de N(V0**2+W0**2) pour avoir une matrice definie
C meme en cas de condition de NEUMANN nulle

C VERSION 2:
C    n=0  V=v0      W=w0
C
C(1) m=0  Vn+1m=Vn  Wn+1m=Wn
C
C(2) [Rho/dt -Alfa LAPLACIEN +N(V0**2+W0**2)] Vn+1m+1 = 
C     Fr(tn+1) +Fi(tn+1) + Rho/dt (Vn +Wn+1m-Wn) +Alfa LAPLACIEN Wn+1m
C     +OmegaZ ( x d/dy - y d/dx )(Vn+1m-Wn+1m)
C     -N(Vn+1m**2+Wn+1m**2)(Vn+1m+Wn+1m) + N(V0**2+W0**2) Vn+1m
C
C    [Rho/dt -Alfa LAPLACIEN +N(V0**2+W0**2)] Wn+1m+1 =
C    -Fr(tn+1) +Fi(tn+1) + Rho/dt (-Vn+1m+Vn +Wn) -Alfa LAPLACIEN Vn+1m
C    +OmegaZ ( x d/dy - y d/dx )(Vn+1m+Wn+1m)
C    +N(Vn+1m**2+Wn+1m**2)(Vn+1m-Wn+1m) + N(V0**2+W0**2) Wn+1m
C                       
C    Si ||Vn+1m+1-Vn+1m||>Eps||Vn+1m+1|| ou ||Wn+1m+1-Wn+1m||>Eps||Wn+1m+1||
C    Alors m=m+1  Aller en (2)
C    Sinon Vn+1=Vn+1m+1  Wn+1=Wn+1m+1
C
C   Si TPSINI+(n+1)dt < TPSFIN Alors n=n+1; Aller en (1)
C   Fin
C
C L'IDEE EST DE CALCULER LA MATRICE GLOBALE PROFIL  n x n
C [Rho/dt -Alfa LAPLACIEN +N(V0**2+W0**2)] et de la FACTORISER L D tL CROUT
C POUR RESOUDRE LE SYSTEME
C

C TESTNL=7: SCHEMA IMPLICITE en TEMPS avec Matrice GLOBALE 2nx2n
C --------------------------------------------------------------
C   Integration en temps par EULER a PAS CONSTANT et IMPLICITE MODIFIE
C   n=0  Vn=v0     Wn=w0
C
C(3)m=0  Vn+1m=Vn  Wn+1m=Wn
C
C A L'ETAPE m+1, LE PROBLEME CONSISTE A TROUVER Un+1 SOLUTION DE
C   [-(-Alfa Laplac+N(V0**2+W0**2)); (Rho/dt-OmegaZ(xd/dy-yd/dx)]{Vn+1m+1}
C(4)[(Rho/dt-OmegaZ(xd/dy-yd/dx);   -Alfa Laplac +N(V0**2+W0**2)]{Wn+1m+1} =
C {-(Fr(tn+1)-Rho/dt Wn + (N(V0**2+W0**2)-N(Vn+1m**2+Wn+1m**2)) Vn+1m+1) }
C {  Fi(tn+1)+Rho/dt Vn + (N(V0**2+W0**2)-N(Vn+1m**2+Wn+1m**2)) Wn+1m+1  }

C  Si ||Vn+1m+1-Vn+1m||>Eps||Vn+1m+1|| ou ||Wn+1m+1-Wn+1m||>Eps||Wn+1m+1||
C  Alors m=m+1  Aller en (4)
C  Sinon Vn+1=Vn+1m+1  Wn+1=Wn+1m+1
C
C  Si TPSINI+(n+1)dt < TPSFIN Alors n=n+1; Aller en (3)
C  Fin
C
C L'IDEE EST DE CALCULER UNE SEULE FOIS LA MATRICE GLOBALE PROFIL 2n x 2n
C  [-Alfa Laplac+N(V0**2+W0**2); -(Rho/dt+OmegaZ(xd/dy-yd/dx)]
C  [(Rho/dt+OmegaZ(xd/dy-yd/dx); -Alfa Laplac +N(V0**2+W0**2)]
C et de la FACTORISER L D tL de CROUT POUR RESOUDRE LE SYSTEME LINEAIRE


CCC TESTNL=8: SCHEMA EXPLICITE en TEMPS SANS Matrice GLOBALE Condition CFL: NON TESTE
CCC -----------------------------------------------------------------------
CCC  Integration en temps par EULER a PAS CONSTANT et EXPLICITE
CCC -Wn+1 - dt Alfa LAPLACIEN Vn+1 + dt Beta (Vn**2+Wn**2) Vn+1 = -Wn + dt Fr(t,X)
CCC  Vn+1 - dt Alfa LAPLACIEN Wn+1 + dt Beta (Vn**2+Wn**2) Wn+1 =  Vn + dt Fi(t,X)
CCC
CCC Traitement de la non linearite par des iterations m de Point Fixe:
CCC    n=0  V=v0      W=w0
CCC(5) m=0  Vn+10=Vn  Wn+10=Wn
CCC(6) -Rho(Wn+1m+1-Wn)/dt -Alfa LAPLAC Vn+1m   +Beta*(Vn+1m**2+Wn+1m**2)Vn+1m=Fr(
CCC     Rho(Vn+1m+1-Vn)/dt -Alfa LAPLAC Wn+1m+1 +Beta*(Vm  n+1**2+Wn+1m+1**2)Wn+1m
CCC     Rho(Vn+1m+2-Vn)/dt -Alfa LAPLAC Wn+1m+1 +Beta*(Vn+1m+1**2+Wn+1m+1**2)Wn+1m
CCC    -Rho(Wn+1m+2-Wn)/dt -Alfa LAPLAC Vn+1m+2 +Beta*(Vn+1m+2**2+Wn+1m+1**2)Vn+1m
CCC
CCC     Si ||Vn+1m+2-Vn+1m||>Eps||Vn+1m|| ou ||Wn+1m+2-Wn+1m||>Eps||Wn+1m||
CCC     Alors m=m+2  Aller en (6)
CCC     Sinon Vn+1=Vn+1m+2  Wn+1=Wn+1m+2
CCC   Si TPSINI+(n+1)dt < TPSFIN Alors n=n+1; Aller en (5)
CCC   Fin


C TESTNL=9: SCHEMA SEMI-IMPLICITE MODIFIE en TEMPS avec Matrice GLOBALE nxn
C           pour l'equation de GROSS-PITAEVSKII
C -------------------------------------------------------------------------
C   n=0  V=v0      W=w0
C(7)m=0  V0n+1=Vn  W0n+1=Wn

C(8)(Rho/dt -LAPLACIEN/2 + ( V +N(V0**2+W0**2) ) Vn+1m+1 =
C        Rho Vn/dt + OmegaZ (x dWn+1m/dy - y dWn+1m/dx) 
C        -( N(Vn+1m**2+Wn+1m**2)-N(V0**2-W0**2) ) Vn+1m + Fr

C   (Rho/dt -LAPLACIEN/2 + ( V +N(V0**2+W0**2) ) Wn+1m+1 =
C        Rho Wn/dt - OmegaZ (x dVn+1m/dy - y dVn+1m/dx) 
C        -( N(Vn+1m**2+Wn+1m**2)-N(V0**2-W0**2) ) Wn+1m + Fi

C    NormeL2 = (Integrale |Vn+1m+1**2+Wn+1m+1**2| dX)**1/2
C    Vn+1m+1 = Vn+1m+1 / NormeL2, Wn+1m+1 = Wn+1m+1 / NormeL2

C     Si ||Vn+1m+1-Vn+1m||>Eps||Vn+1m+1|| ou ||Wn+1m+1-Wn+1m||>Eps||Wn+1m+1||
C     Alors m=m+1  Aller en (8)
C     Sinon Vn+1=Vn+1m+1  Wn+1=Wn+1m+1
C   Si TPSINI+(n+1)dt < TPSFIN Alors n=n+1; Aller en (7)
C   Fin


C TESTNL=10: NLSE-GPE APPROXIMATION par ELEMENTS FINIS P1 PRODUISANT
C            UNE EQUAION DU 3-EME DEGRE EN CHACUN DES SOMMETS
C            SANS Matrice GLOBALE nxn
C ------------------------------------------------------------------
C    n=0  V=v0      W=w0
C
C(9) m=0  Vn+1m=Vn  Wn+1m=Wn
C
C(10)-Beta Vn+1m+1**3 
C -  Beta Wn+1m * Vn+1m+1**2
C + {  Rho/dt -V(X) -Beta Wn+1m**2 +[ Alfa d/dx +OmegaZ y, Alfa d/dy -OmegaZ x] [D] } Vn+1m
C + { -Rho/dt -V(X) -Beta Wn+1m**2 +[ Alfa d/dx -OmegaZ y, Alfa d/dy +OmegaZ x] [D] } Wn+1m
C + Rho/dt (Wn-Vn) - Fr -Fi = 0
C   Resolution de l'equation du 3-eme degre en Vn+1m+1(Sommet) en tous les sommets non Dirichlet

C - Beta Wn+1m+1**3 
C + Beta Vn+1m * Wn+1m+1**2
C + {  Rho/dt -V(X) -Beta Vn+1m**2 +[ Alfa d/dx +OmegaZ y, Alfa d/dy -OmegaZ x] [D] } Wn+1m
C - { -Rho/dt -V(X) -Beta Vn+1m**2 +[ Alfa d/dx -OmegaZ y, Alfa d/dy +OmegaZ x] [D] } Vn+1m
C - Rho/dt (Wn+Vn) + Fr -Fi = 0
C   Resolution de l'equation du 3-eme degre en Wn+1m+1(Sommet) en tous les sommets non Dirichlet

C          [d/dx]
C ou [D] = [    ] l'operateur gradient
C          [d/dy]
                    
C    Si ||Vn+1m+1-Vn+1m||>Eps||Vn+1m+1|| ou ||Wn+1m+1-Wn+1m||>Eps||Wn+1m+1||
C    Alors m=m+1  Aller en (10)
C    Sinon Vn+1=Vn+1m+1  Wn+1=Wn+1m+1
C
C   Si TPSINI+(n+1)dt < TPSFIN Alors n=n+1; Aller en (9)
C   Fin


C REMARQUES sur les DONNEES A FOURNIR
C -----------------------------------
C Rho    : Coefficient de DENSITE de MASSE
C ALFA   : Coefficient des -LAPLACIENS doit etre donne comme une
C          CONDUCTIVITE HOMOGENE ISOTROPE
C N ou BETA: Coefficient des TERMES NON LINEAIRES doit etre donne
C          comme une FONCTION DEVANT la TEMPERATURE
C CONTACT: Condition de Dirichlet est la seule condition aux limites permise
C
C TOUTES les autres donnees fournies NE seront PAS PRISES EN COMPTE
C
C ENTREE :
C --------
C KNOMOB : NOM DE L'OBJET A TRAITER
C
C SORTIE :
C --------
C IERR   : 0 SI PAS D'ERREUR D'EXECUTION, NON NUL SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Texas A & M University at QATAR  Fevrier 2011
C MODIFS : ALAIN PERRONNET LJLL UPMC & Saint Pierre du Perray  Aout 2011
C MODIFS : ALAIN PERRONNET LJLL UPMC Saint Pierre du Perray Novembre2013
C MODIFS : ALAIN PERRONNET LJLL UPMC Saint Pierre du Perray Mai     2014
C23456---------------------------------------------------------------012
      PARAMETER     (MXTYEL=7)
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/ctemps.inc"
      include"./incl/cthet.inc"
      include"./incl/cnonlin.inc"
      include"./incl/donthe.inc"
      include"./incl/a___vecteur.inc"
      include"./incl/a___xyznoeud.inc"
C
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE      (MCN(1), RMCN(1), DMCN(1))
C
      CHARACTER*(*)     KNOMOB
      INTEGER           NUMIOB(4), NUMAOB(4), MNDOEL(4), MXDOEL(4)
      REAL              RDeltat, DTSTOC, TPSINI, TPSFIN
      DOUBLE PRECISION  DINFO, D2PI, DCPU, DVECT
      DOUBLE PRECISION  RELMIN
      DATA              RELMIN/-1D28/
C
C     INITIALISATION DU TEMPS CALCUL INITIAL
      DCPU  = DINFO( 'CPU' )
      D2PI  = ATAN(1D0) * 8D0
C
C     INSTANT ACTUEL DU CALCUL
      TEMPS  = 0.0
      MNTIME = 0
C     LE NOMBRE DE MOTS D'UNE VARIABLE REEL2
      MOREE2 = MOTVAR(6)
C     EN CAS DE SORTIE PREMATUREE
      DVECT  = 0D0
      MNX    = 0
      NBCOOR = 0
C
      MNNPEF = 0
      MNDOEL(1) = 0
      MNDOEL(2) = 0
      MNDOEL(3) = 0
      MNDOEL(4) = 0
      MNTPOB = 0
      MNTAUX = 0
      MNTAEL = 0
      MNTHER = 0
      MNNODL = 0
      MONODL = 0
      MNX    = 0
      MNTHDL = 0
C
C     NON LINEAIRE SCHRODINGER EQUATIONS A RESOUDRE
C     ---------------------------------------------
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10000) KNOMOB
      ELSE
         WRITE(IMPRIM,20000) KNOMOB
      ENDIF
10000 FORMAT(/1X,100('*')/
     %' RESOLUTION  i Rho du(t,X)/dt -Alfa LAPLACIEN u(t,X) +N(|u(t,x)|*
     %*2) u(t,X) = F  sur l''objet ',A/1X,100('*'))
20000 FORMAT(/1X,100('*')/
     %' SOLUTION  i Rho du(t,X)/dt -Alfa LAPLACIEN u(t,X) +N(|u(t,x)|**2
     %) u(t,X) = F on the object ',A/1X,100('*'))
C
C     L'ANCIEN HISTORIQUE EST EFFACE
      CALL RECTEF( NRHIST )
      NBLGRC(NRHIST) = 1
      IF( LANGAG .EQ. 0 ) THEN
         KHIST(1) = 'OBJET: ' // KNOMOB
      ELSE
         KHIST(1) = 'OBJECT: ' // KNOMOB
      ENDIF
      CALL LHISTO
C
C     ENTREEE DES DONNEES THERMIQUES POUR L'INSTANT INITIAL POUR
C     UN PROBLEME NON LINEAIRE DE SCHRODINGER INSTATIONNAIRE
C     AVEC UN SCHEMA D'ORDRE 1 EN TEMPS
C     ==========================================================
      CALL THED1T( KNOMOB, MOREE2, NOAXIS,
     %             NTLXOB, MNTOPO, NDIM,   NDPGST, MNXYZP, MNXYZN,
     %             NBTYEL, MNNPEF, NBTTEF,
     %             MNTPOB, NBDLMX, MOAUX,  MNTAUX,
     %             NUMIOB, NUMAOB, NBOBIN, MNOBIN, NBOBCL, MNOBCL,
     %             MXTYEL, MXDOEL, MNDOEL,
     %             IEMAST, IECHMA, IECOND, IEDILA, IEVIFL, IECOET,
     %             IESOIN, IECONT, IEECHA, IESOCL, IESOPO,
     %             MNTHER, MOTAEL, MNTAEL, MNX,
     %             NORESO, MNLPLI, MNLPCO, NIVEAU, NBRDAG,
     %             NCODSK, NCODSM, MNUG0,
     %             RDeltat,DTSTOC, TPSINI, TPSFIN, NBVECT, MXVECT,
     %             NTDL,   NTVECT, MNVECT,
     %             IERR )
      IF( IERR .NE. 0 ) GOTO 9990
C
C     NOMBRE DE COMPOSANTES DES COORDONNEES DE XYZNOEUD
      NBCOOR = MCN(MNXYZN+WBCOON)
C
C     EXECUTION SELON LE CHOIX DE LA METHODE DE RESOLUTION
      IF( TESTNL .EQ. 6 ) THEN
C
C        NLSE: INTEGRATION SEMI-IMPLICITE MODIFIEE AVEC MATRICE GLOBALE nxn
C        Algo2: i du(t)/dt -Alfa Laplacian u(t) + N(|u|**2) u(t)=0
C        ==================================================================
         CALL NLSESIMPL( KNOMOB, MOREE2, NDIM,   NDPGST, MNXYZN,
     %                   NBTYEL, MNNPEF, MNTPOB, MNTAUX,
     %                   NUMIOB, NUMAOB, MNDOEL, IESOPO,
     %                   RELMIN, MNTHER, MNTAEL, MNX,
     %                   NORESO, MNLPLI, MNLPCO, NBRDAG, NCODSK,
     %                   MNUG0,  RDeltat,DTSTOC, TPSINI, TPSFIN,
     %                   NBVECT, MXVECT, NTDL,
     %                   NTVECT, MNVECT, IERR )
C
      ELSE IF( TESTNL .EQ. 7 ) THEN
C
C        NLSE: INTEGRATION IMPLICITE AVEC MATRICE GLOBALE 2nx2n
C        Algo1: i du(t)/dt -Alfa Laplacian u(t) + N(|u|**2) u(t)=0
C        ===========================================================
         CALL NLSEIMPL(  KNOMOB, MOREE2, NDIM,   NDPGST, MNXYZN,
     %                   NBTYEL, MNNPEF, MNTPOB, MNTAUX,
     %                   NUMIOB, NUMAOB, MNDOEL, IESOPO,
     %                   RELMIN, MNTHER, MNTAEL, MNX,
     %                   NORESO, MNLPLI, MNLPCO, NBRDAG, NCODSK,
     %                   MNUG0,  RDeltat,DTSTOC, TPSINI, TPSFIN,
     %                   NBVECT, MXVECT, NTDL,
     %                   NTVECT, MNVECT, IERR )
C
      ELSE IF( TESTNL .EQ. 8 ) THEN
C
C        NLSE: INTEGRATION EN TEMPS EXPLICITE SANS MATRICE GLOBALE
C        =========================================================
         CALL NLSEEXPL( MOREE2, NDIM,   NDPGST, MNXYZN, NBTYEL, MNNPEF,
     %                  NUMIOB, NUMAOB, MNDOEL,
     %                  RELMIN, MNUG0,  RDeltat,DTSTOC, TPSINI, TPSFIN,
     %                  NBVECT, MXVECT, NTDL,
     %                  NTVECT, MNVECT, IERR )
C
      ELSE IF( TESTNL .EQ. 9 ) THEN
C
C        GROSS-PITAEVSKI equation :
C        idu(t,X)/dt +LAPLACIEN u(t,X)/2 - ( V +Beta |u(t,X)|**2 ) u(t,X)
C                    -i OmegaZ (x du/dy - y du/dx) =F
C        avec t->it la I-TIME METHODE et
C        N( |u|**2 ) = r**2/2 + r**4/4 + Beta ( v**2 + w**2 )  devient
C        Rho dU(t,X)/dt + Alfa LAPLACIEN U(t,X) - N(|U**2|) U(t,X)
C                       + i OmegaZ (x dU/dy - y dU/dx) = -F
C        Methode SEMI-IMPLICITE MODIFIEE AVEC MATRICE GLOBALE nxn
C        ================================================================
         CALL NLSEITIME( KNOMOB, MOREE2, NDIM,   NDPGST, MNXYZN,
     %                   NBTYEL, MNNPEF, MNTPOB, MNTAUX,
     %                   NUMIOB, NUMAOB, MNDOEL, IESOPO,
     %                   RELMIN, MNTHER, MNTAEL, MNX,
     %                   NORESO, MNLPLI, MNLPCO, NBRDKG, NCODSK,
     %                   MNUG0,  RDeltat,DTSTOC, TPSINI, TPSFIN,
     %                   NBVECT, MXVECT, NTDL,
     %                   NTVECT, MNVECT, IERR )
C
      ELSE IF( TESTNL .EQ. 10 ) THEN
C
C        GROSS-PITAEVSKI equation :
C        APPROXIMATION ELEMENTS FINIS P1 PRODUISANT UNE EQUATION
C        DU 3-EME DEGRE EN CHACUN DES SOMMETS et SANS Matrice GLOBALE nxn
C        i Rho du(t,X)/dt +Laplacian/2 u(t,X) -( V +Beta |u(t,X)|**2 ) u(t,X)
C                     -i OmegaZ (x d/dy - y d/dx) u(t,X) = F(t,X)
C        ====================================================================
         CALL GPE3DEGP1( KNOMOB, MOREE2, NDIM,   NDPGST, MNXYZN,
     %                   NBTYEL, MNNPEF, MNTAUX,
     %                   NUMIOB, NUMAOB, MNDOEL, RELMIN, MNX,
     %                   MNUG0,  RDeltat,DTSTOC, TPSINI, TPSFIN,
     %                   NBVECT, MXVECT, NTDL,
     %                   NTVECT, MNVECT, IERR )
C
      ENDIF
C
C     DESTRUCTION DES TABLEAUX DEVENUS INUTILES
C     =========================================
 9990 IF( MNNPEF .GT. 0 ) CALL TNMCDS( 'ENTIER', 2*MXTYEL, MNNPEF )
      DO I=1,4
         IF( MNDOEL(I) .GT. 0 .AND. MXDOEL(I) .GT. 0 ) THEN
            CALL TNMCDS( 'ENTIER', MXDOEL(I), MNDOEL(I) )
         ENDIF
      ENDDO
      IF( MNTHDL .GT. 0 ) CALL NLDATADS
      IF( MNTAUX .GT. 0 ) CALL TNMCDS( 'REEL2',  MOAUX,  MNTAUX )
      IF( MNTAEL .GT. 0 ) CALL TNMCDS( 'REEL2',  MOTAEL, MNTAEL )
      IF( MNTHER .GT. 0 ) CALL TNMCDS( 'REEL2',  128,    MNTHER )
      IF( MNNODL .GT. 0 ) CALL TNMCDS( 'ENTIER', MONODL, MNNODL )
      IF( MNX    .GT. 0 ) CALL TNMCDS( 'REEL' ,  NBDLMX*NBCOOR, MNX )
      IF( MNTPOB .GT. 0 ) CALL TNMCDS( 'ENTIER', NBTYEL*MXPOBA, MNTPOB )
C
C     BILAN SUR LE TEMPS CALCUL UTILISE
C     =================================
C     COUT TOTAL DU CALCUL DE LA SOLUTION
      DVECT = DINFO( 'CPU' )
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,19000) DVECT-DCPU
      ELSE
         WRITE(IMPRIM,29000) DVECT-DCPU
      ENDIF
19000 FORMAT(
     %'TEMPS TOTAL  DE LA RESOLUTION =',F12.2,' SECONDES CPU'/)
29000 FORMAT(
     %'CPU TOTAL SOLUTION COMPUTATION=',F12.2,' CPU SECONDS'/)
C
C     SORTIE
      RETURN
      END
