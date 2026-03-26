      SUBROUTINE NLSEEXPL(MOREE2, NDIM,  NDPGST, MNXYZN, NBTYEL, MNNPEF,
     %                   NUMIOB, NUMAOB, MNDOEL,
     %                   RELMIN, MNUG0,  RDeltaT,DTSTOC, TPSINI, TPSFIN,
     %                   NBWAVE, MXVECT, NTDL,   NTVECT, MNVECT, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCULER u(t,X) solution de l'EQUATION NON LINEAIRE de SCHRODINGER
C ----- SOLUTION: i du(t,X)/dt -alfa LAPLACIEN u(t,X) +beta N(u(t,X)**2)U = 0
C       CALCUL de la PARTIE REELLE V et IMAGINAIRE W de u(t,X)=v(t,X)+i w(t,X)
C       selon le systeme
C       -dw(t,X)/dt + alfa LAPLACIEN v(t,X) + beta (v**2+w**2) v(t,X) = Fr(t,X)
C        dv(t,X)/dt + alfa LAPLACIEN w(t,X) + beta (v**2+w**2) w(t,X) = Fi(t,X)
C
C       Integration en temps par EULER a PAS CONSTANT et IMPLICITE, soit:
C       -wn+1 + dt alfa LAPLACIEN vn+1 + dt beta (vn+1**2+wn+1**2) vn+1 = -wn +
C        vn+1 + dt alfa LAPLACIEN wn+1 + dt beta (vn+1**2+wn+1**2) wn+1 =  vn +
C
C   Traitement de la non linearite par des iterations m de Point Fixe:
C   n=0  V=v0  W=w0
C(1)n=n+1
C     m=0
C     V0n+1 = Vn   W0n+1 = Wn
C(2) -(Wm+1n+1-Wn)/dt +alfa LAPLAC Vm  n+1 +beta*(Vm  n+1**2+Wm  n+1**2)Vm  n+1-
C     (Vm+1n+1-Vn)/dt +alfa LAPLAC Wm+1n+1 +beta*(Vm  n+1**2+Wm+1n+1**2)Wm+1n+1-
C     (Vm+2n+1-Vn)/dt +alfa LAPLAC Wm+1n+1 +beta*(Vm+1n+1**2+Wm+1n+1**2)Wm+1n+1-
C    -(Wm+2n+1-Wn)/dt +alfa LAPLAC Vm+2n+1 +beta*(Vm+2n+1**2+Wm+1n+1**2)Vm+2n+1-
C
C     Si ||Vm+2n+1-Vm n+1||>Eps||Vm n+1|| ou ||Wm+2n+1-Wm n+1||>Eps||Wm n+1||
C     Alors m=m+2  Aller en (2)
C     Sinon Vn+1=Vm+2n+1  Wn+1=Wm+2n+1
C   Si TPSINI+(n+1)dt < TPSFIN Alors Aller en (1)
C   Fin
C
C         ===========================================================
C         ATTENTION: Actuellement seul le TRIANGLE 2P1D est PROGRAMME
C         ===========================================================
C
C ENTREES:
C --------
C KNOMOB : NOM DE L'OBJET A TRAITER
C MOREE2 : NOMBRE DE MOTS   D'UNE VARIABLE REELLE DOUBLE PRECISION
C NDIM   : DIMENSION DES COORDONNEES DES POINTS ( 1 OU 2 OU 3 )
C NDPGST : CODE D'IDENTIFICATION DES SOMMETS POINTS NOEUDS
C MNXYZN : ADRESSE MCN DU TABLEAU XYZNOEUD DE L'OBJET KNOMOB
C
C NBTYEL : NOMBRE DE TYPES D'EF DU MAILLAGE DE CET OBJET
C MNNPEF : ADRESSE MCN DU TABLEAU DES ADRESSES MCN DES TMS NPEF"TYPE EF
C NUMIOB : NUMERO MINIMAL DU PLSV DANS LA DEFINITION DE L'OBJET
C NUMAOB : NUMERO MAXIMAL DU PLSV DANS LA DEFINITION DE L'OBJET
C NBOBIN : NOMBRE DE VOLUMES en 3D, SURFACES en 2D DE L'OBJET
C MNOBIN : ADRESSE MCN DU DEBUT DU TABLEAU NUOBIN (PARTIE DE TOPOLOGIE)
C NBOBCL : NOMBRE DE PLS en 3D, PL en 2D DE L'OBJET
C MNOBCL : ADRESSE MCN DU DEBUT DU TABLEAU NUOBCL (PARTIE DE TOPOLOGIE)
C MXTYEL : NOMBRE MAXIMAL DE TYPES D'EF (7)
C MXDOEL : NOMBRE DE MOTS DECLARES DES TABLEAUX D'ADRESSE DANS MNDOEL
C MNDOEL : LES 4 ADRESSES MCN DES TABLEAUX DES ADRESSES DES
C          TABLEAUX DECRIVANT LES DONNEES THERMIQUES DE L'OBJET COMPLET
C
C RELMIN : PLUS PETIT REEL SERVANT DE MARQUEUR DE NON UTILISATION
C MNUG0  : ADRESSE MCN DE U0 SOLUTION  = V0+iW0  A L'INSTANT t
C          MNVG0 = MNUG0 et MNWG0 = MNUG0 + NBSOM * MOREE2
C
C RDeltaT: PAS CONSTANT DU TEMPS REEL SIMPLE PRECISION
C DTSTOC : PAS CONSTANT DU TEMPS ENTRE 2 STOCKAGES DU VECTEUR"ONDENLSE'
C TPSINI : TEMPS INITIAL DU CALCUL
C TPSFIN : TEMPS FINAL   DU CALCUL
C NBWAVE : NUMERO DU DERNIER VECTEUR SOLUTION STOCKE
C MXVECT : NOMBRE DE VECTEUR"ONDENLSE A STOCKER
C NTDL   : NOMBRE TOTAL DE DEGRES DE LIBERTE DE L'OBJET
C          2 FOIS LE NOMBRE DE NOEUDS DU MAILLAGE
C
C SORTIES:
C --------
C NTVECT : NUMERO      DU TMS VECTEUR"ONDENLSE DE L'OBJET
C MNVECT : ADRESSE MCN DU TMS VECTEUR"ONDENLSE DE L'OBJET
C IERR   : 0 SI PAS D'ERREUR D'EXECUTION, NON NUL SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Texas A & M University at QATAR  Fevrier 2011
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St Pierre du Perray  Juillet 2011
C23456---------------------------------------------------------------012
      PARAMETER         (MXTYEL=7)
      PARAMETER         (ITERMX=16)
      PARAMETER         (ITERNLMX=16)
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/donthe.inc"
      include"./incl/donele.inc"
      include"./incl/a_objet__topologie.inc"
      include"./incl/a_objet__definition.inc"
      include"./incl/a_objet__erreurth.inc"
      include"./incl/a___arete.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___xyznoeud.inc"
      include"./incl/a___xyzpoint.inc"
      include"./incl/a___npef.inc"
      include"./incl/a___conductivite.inc"
      include"./incl/a___dilatation.inc"
      include"./incl/a___source.inc"
      include"./incl/a___contact.inc"
      include"./incl/a___vecteur.inc"
      include"./incl/a___dtemperature.inc"
      include"./incl/a___temperinit.inc"
      include"./incl/msvaau.inc"
      include"./incl/ponoel.inc"
      include"./incl/ctemps.inc"
      include"./incl/cthet.inc"
      include"./incl/cnonlin.inc"
C
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE      (MCN(1), RMCN(1), DMCN(1))
C
      DOUBLE PRECISION  RELMIN
      INTEGER           NUMIOB(4), NUMAOB(4), MNDOEL(4)
C
      REAL              RDeltaT, TSTOC,  DTSTOC, TPSINI, TPSFIN
      DOUBLE PRECISION  DT,      ALFA,   DTALFA, BETA,   DTBETA,
     %                  NORDIFV, NORVECV, NORDIFW,  NORVECW
      DOUBLE PRECISION  XYZD(3), DWMIN,   DWMAX
      DOUBLE PRECISION  WCAMIN(2),  WCAMAX(2),  WEXMIN(2), WEXMAX(2)
      REAL              XYZWMIN(3), XYZWMAX(3), ERRMAX(2)
      INTEGER           NONOEF(4),  NOFOWE(2)
      CHARACTER*4       NOMELE(2)
      INTRINSIC         REAL
C
C     L'ADRESSE DE LA SOLUTION A L'ITERATION 0 DU POINT FIXE
      MNTHET = MNUG0
      MNBG   = 0
      MNNOFX = 0
      MNVAFX = 0
      MNMGD  = 0
C
C     PAS DE TEMPS EN DOUBLE PRECISION
      DT = RDeltaT
C
C     NOMBRE DE SOMMETS DU MAILLAGE
      NBSOM  = MCN( MNXYZN + WNBNOE )
C
C     NOMBRE DE DEGRES DE LIBERTE TOTAL
C    (NBSOM POUR PARTIE REELLE et NBSOM POUR PARTIE IMAGINAIRE)
      NTDL = 2 * NBSOM
C
C     AFFICHAGE DE LA SOLUTION AU TEMPS INITIAL
C     =========================================
      CALL AFNLSE( 20,  1, MNXYZN, NTDL/2, 1, MCN(MNUG0),
     %             WCAMIN, WCAMAX, NOFOWE, WEXMIN, WEXMAX, ERRMAX )
C
C     ICI  MAILLAGE DE SEULEMENT DES TRIANGLES 2P1D
C     L'ADRESSE DU TABLEAU NPEF"TRIA 2P1D
      MNELE = MCN( MNNPEF )
C
C     LE NUMERO DU TYPE DE L'ELEMENT FINI
      NUTYEL = MCN( MNELE + WUTYEL )
C
C     ATTENTION: Actuellement seul le TRIANGLE 2P1D est PROGRAMME
      IF( NUTYEL .EQ. 13 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            NBLGRC(NRERR) = 2
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1)='ELEMENT FINI NON TRIA 2P1D NON PROGRAMME'
               KERR(2)='ELEMENT FINI NON ENCORE PROGRAMME'
            ELSE
               KERR(1)='FINITE ELEMENT DIFFERENT OF TRIANGLE 2P1D'
               KERR(2)='FINITE ELEMENT NOT YET PROGRAMMED'
            ENDIF
            CALL LEREUR
            GOTO 9999
         ENDIF
      ENDIF
C
C     LE NOMBRE D'ELEMENTS FINIS DE CE TYPE
      NBELEM = MCN( MNELE + WBELEM )
C
C     L'ADRESSE MCN DES NOEUDS ET POINTS GEOMETRIQUES DES EF DE CE TYPE
      MNPGEL = MNELE + WUNDEL
C
C     LES CARACTERISTIQUES DE L'ELEMENT FINI
      CALL ELNUNM( NUTYEL, NOMELE )
      CALL ELTYCA( NUTYEL )
C
C     LES DONNEES DU PROBLEME DE L'ONDE SONT RETROUVEES SUR LE PREMIER TRIANGLE
C     LES NOEUDS DE L'ELEMENT FINI NUELEM=1
      NUELEM = 1
      CALL EFNOEU( MNELE, NUELEM, NBNDEL, NONOEF )
C     LE NUMERO DE VOLUME  DE L'EF
C     LE NUMERO DE SURFACE DES FACES   DE L'EF
C     LE NUMERO DE LIGNE   DES ARETES  DE L'EF
C     LE NUMERO DE POINT   DES SOMMETS DE L'EF
      CALL EFPLSV( MNELE , NUELEM,
     %             NOVCEL, NOSFEL, NOLAEL, NOPSEL,
     %             NOOBVC, NOOBSF, NOOBLA, NOOBPS )
C
C     RECHERCHE DE ALFA COEFFICIENT DU LAPLACIEN DES NLSE
C     AU 1-ER SOMMET DU 1-ER TRIANGLE
      NS = MCN( MNPGEL )
      MN = MNXYZN + WYZNOE + 3 * (NS-1)
      XYZD(1) = RMCN(MN)
      XYZD(2) = RMCN(MN+1)
      XYZD(3) = RMCN(MN+2)
      CALL REALFA( XYZD, NOOBSF, NUMIOB(3), NUMAOB(3), MCN(MNDOEL(3)),
     %             ALFA )
C
C     RECHERCHE DE BETA COEFFICIENT DU COEFFICIENT DE LA "TEMPERATURE"
C     AU 1-ER SOMMET DU 1-ER TRIANGLE
      CALL REBETA( XYZD, NOOBSF, NUMIOB(3), NUMAOB(3), MCN(MNDOEL(3)),
     %             BETA )
C
C     MULTIPLICATION PAR LE PAS DE TEMPS dt
      DTALFA = DT * ALFA
      DTBETA = DT * BETA
C
C     LISTE DES NUMEROS DES NOEUDS FIXES SUR LA FRONTIERE
C     ET D'UNE VALEUR CORRESPONDANTE A LA PARTIE REELLE V et
C     PARTIE IMAGINAIRE W  FIXEE  VUE ICI COMME UNE TEMPERATURE
C     POUR UNE CONDITION DE DIRICHLET HOMOGENE CETTE VALEUR EST ZERO
      NCDLF0 = 1
      CALL THDLFX( NCDLF0, NTDL,   NDIM,
     %             NBTYEL, MNNPEF, NDPGST,
     %             MNXYZN, NUMIOB, MNDOEL, RELMIN,
     %             NBNOFX, MONOFX, MNNOFX, MOVAFX, MNVAFX, IERR )
C
C     CONSTRUCTION DE LA MATRICE DE MASSE DIAGONALE POUR LE TRIANGLE 2P1D
C     GRACE A UNE INTEGRATION NUMERIQUE AUX 3 SOMMETS DES TRIANGLES
      CALL TNMCDC( 'REEL2', NBSOM, MNMGD )
      CALL TP1P1DIAG( NBSOM,  MCN(MNXYZN+WYZNOE),
     %                NBELEM, MCN(MNPGEL), MCN(MNMGD) )
C
C     LES 2 VECTEURS GLOBAUX DES SOLUTIONS A L'INSTANT 0 et 1 SONT DECOUPES
C     ---------------------------------------------------------------------
C     EN PARTIE REELLE V ET PARTIE IMAGINAIRE W
C     MNUG0: ADRESSE MCN DE U0 SOLUTION V0+iW0  A L'INSTANT INITIAL
      MNVG0 = MNUG0
      MNWG0 = MNUG0 + NBSOM * MOREE2
C
      CALL TNMCDC( 'REEL2', NTDL, MNUG1 )
C     MNUG1 ADRESSE MCN DE U0 SOLUTION  = V0+iW0  A L'INSTANT t+dt
C     MNVG1 = MNUG1 et MNWG1 = MNUG1 + NBSOM * MOREE2
C
C     MNUG1: ADRESSE MCN DE U1 SOLUTION V1+iW1
      MNVG1 = MNUG1
      MNWG1 = MNUG1 + NBSOM * MOREE2
C
C     VECTEUR GLOBAL SECOND MEMBRE
      CALL TNMCDC( 'REEL2', NBSOM, MNBG )
C
C     LE PROCHAIN TEMPS POUR STOCKER LE VECTEUR SOLUTION
      TSTOC  = TPSINI + DTSTOC
      MNWAVE = MNVECT + WECTEU + MOREE2 * NTDL * (NBWAVE-1)
      MNTIME = MNVECT + WECTEU + MOREE2 * NTDL * MXVECT - 1
C
C     LE NOMBRE DE DL D'UN VECTEUR UG = VG + i WG = 2 NBSOM
      MCN( MNVECT + WBCOVE ) = NTDL
C
C     LE NOMBRE DE VECTEURS STOCKES
      MCN( MNVECT + WBVECT ) = NBWAVE
C
C     AU DEPART PAS DE TEMPS INDIQUES
      MCN( MNVECT + WBCPIN ) = 0
C
C     ##############################################################
C     ##                                                          ##
C     ##  LA BOUCLE EN TEMPS AVEC DES PAS DE TEMPS CONSTANTS = DT ##
C     ##                                                          ##
C     ##############################################################
C
C     LE NOUVEAU TEMPS OU SE FAIT LE CALCUL
 100  TEMPS = REAL( TEMPS + DT )
C     --------------------------
C
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10100) TEMPS
      ELSE
         WRITE(IMPRIM,20100) TEMPS
      ENDIF
10100 FORMAT(//'Au TEMPS',G14.6,' CALCUL de l''ONDE ',100('*'))
20100 FORMAT(//'At TIME',G14.6,' COMPUTATION of the WAVE ',100('*'))
C
C     DES ITERATIONS DE POINT FIXE SONT MISES EN OEUVRE
      MNTHET = MNUG0
C
C     MNUG0: ADRESSE MCN DE U0 SOLUTION V0+iW0  A L'INSTANT n
      MNVG0 = MNUG0
      MNWG0 = MNUG0 + NBSOM * MOREE2
C
C     MNUG1: ADRESSE MCN DE U1 SOLUTION V1+iW1  A L'INSTANT n+1
      MNVG1 = MNUG1
      MNWG1 = MNUG1 + NBSOM * MOREE2
C
C     ITERATIONS DE POINT FIXE POUR LA NON LINEARITE
C     ==============================================
C     COPIE DE Vn dans V0 n+1
      CALL TRTATD( MCN(MNVG0), MCN(MNVG1), NBSOM )
C     COPIE DE Wn dans W0 n+1
      CALL TRTATD( MCN(MNWG0), MCN(MNWG1), NBSOM )
C
      ITERNL = 0
C
 200  ITERNL = ITERNL + 1
C     -------------------
      IF( ITERNL .GT. ITERNLMX ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(KERR(MXLGER)(1:8),'(I8)') ITERNLMX
            NBLGRC(NRERR) = 2
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1)='NON CONVERGENCE APRES ' // KERR(MXLGER)(1:8)
     %               //' ITERATIONS DE POINT FIXE'
               KERR(2)='ABANDON des ITERATIONS DE POINT FIXE'
            ELSE
               KERR(1)='NO CONVERGENCE after ' // KERR(MXLGER)(1:8)
     %               //' FIX POINT ITERATIONS'
               KERR(2)='EXIT of FIX POINT ITERATIONS'
            ENDIF
            CALL LEREUR
            GOTO 9999
         ENDIF
      ENDIF
C
C     -(Wm+1n+1-Wn)/dt +alfa LAPLAC Vm  n+1
C                      +beta*(Vmn+1**2+Wmn+1**2)Vmn+1=FOmegaReel => Wm+1n+1
C     c-a-d
C     tPP Wm+1n+1 = tPP Wn + dt alfa tDPDP Vmn+1 + dt beta tPU2P Vm n+1 - dt FOm
      CALL NLSESCH( MCN(MNMGD), MCN(MNWG1), MCN(MNWG0),
     %              DTALFA, MCN(MNVG1),
     %              DTBETA, MCN(MNVG1), MCN(MNWG1), MCN(MNVG1), -DT,
     %              MCN(MNBG),
     %              NBNOFX, MCN(MNNOFX), MCN(MNVAFX), 2,
     %              NBSOM,  MCN(MNXYZN+WYZNOE), NBELEM, MCN(MNPGEL),
     %              NOOBSF, NUMIOB(3), NUMAOB(3), MCN(MNDOEL(3)), 1 )
C     AFFICHAGE DE LA SOLUTION AU TEMPS ACTUEL
      CALL AFNLSE( 20,  1, MNXYZN, NTDL/2, 1, MCN(MNVG1),
     %             WCAMIN, WCAMAX, NOFOWE, WEXMIN, WEXMAX, ERRMAX )
C     NORME DE LA DIFFERENCE VG1-VG0 et NORME de VG0
      CALL NORMDIF( NBSOM, MCN(MNVG0), MCN(MNVG1),  NORDIFV, NORVECV )
C     NORME DE LA DIFFERENCE WG1-WG0 et NORME de WG0
      CALL NORMDIF( NBSOM, MCN(MNWG0), MCN(MNWG1),  NORDIFW, NORVECW )
      WRITE(IMPRIM,10000) 'EQ1:', TEMPS, ITERNL,
     %                     NORVECV, NORDIFV, NORDIFV/NORVECV,
     %                     NORVECW, NORDIFW, NORDIFW/NORVECW
C
C     (Vm+1n+1-Vn)/dt +alfa LAPLAC Wm+1n+1
C             -beta*(Vmn+1**2+Wm+1n+1**2)Wm+1n+1=FOmegaImag => Vm+1n+1
C     c-a-d
C     tPP Vm+1n+1 = tPP Vn - dt alfa tDPDP Wm+1n+1 - dt beta tPU2P Wm+1n+1 + dt
      CALL NLSESCH( MCN(MNMGD), MCN(MNVG1), MCN(MNVG0),
     %             -DTALFA, MCN(MNWG1),
     %             -DTBETA, MCN(MNVG1), MCN(MNWG1), MCN(MNWG1), DT,
     %              MCN(MNBG),
     %              NBNOFX, MCN(MNNOFX), MCN(MNVAFX), 1,
     %              NBSOM,  MCN(MNXYZN+WYZNOE), NBELEM, MCN(MNPGEL),
     %              NOOBSF, NUMIOB(3), NUMAOB(3), MCN(MNDOEL(3)), 2 )
C     AFFICHAGE DE LA SOLUTION AU TEMPS ACTUEL
      CALL AFNLSE( 20,  1, MNXYZN, NTDL/2, 1, MCN(MNVG1),
     %             WCAMIN, WCAMAX, NOFOWE, WEXMIN, WEXMAX, ERRMAX )
C     NORME DE LA DIFFERENCE VG1-VG0 et NORME de VG0
      CALL NORMDIF( NBSOM, MCN(MNVG0), MCN(MNVG1),  NORDIFV, NORVECV )
C     NORME DE LA DIFFERENCE WG1-WG0 et NORME de WG0
      CALL NORMDIF( NBSOM, MCN(MNWG0), MCN(MNWG1),  NORDIFW, NORVECW )
      WRITE(IMPRIM,10000) 'EQ2:', TEMPS, ITERNL,
     %                     NORVECV, NORDIFV, NORDIFV/NORVECV,
     %                     NORVECW, NORDIFW, NORDIFW/NORVECW
C
C     (Vm+2n+1-Vn)/dt +alfa LAPLAC Wm+1n+1
C             +beta*(Vm+1n+1**2+Wm+1n+1**2)Wm+1n+1=FOmegaImag => Vm+2n+1
C     c-a-d
C     tPP Vm+2n+1 = tPP Vn - dt alfa tDPDP Wm+1n+1 - dt beta tPU2P Wm+1n+1 + dt
      CALL NLSESCH( MCN(MNMGD), MCN(MNVG1), MCN(MNVG0),
     %             -DTALFA, MCN(MNWG1),
     %             -DTBETA, MCN(MNVG1), MCN(MNWG1), MCN(MNWG1), DT,
     %              MCN(MNBG),
     %              NBNOFX, MCN(MNNOFX), MCN(MNVAFX), 1,
     %              NBSOM,  MCN(MNXYZN+WYZNOE), NBELEM, MCN(MNPGEL),
     %              NOOBSF, NUMIOB(3), NUMAOB(3), MCN(MNDOEL(3)), 2 )
C     AFFICHAGE DE LA SOLUTION AU TEMPS ACTUEL
      CALL AFNLSE( 20,  1, MNXYZN, NTDL/2, 1, MCN(MNVG1),
     %             WCAMIN, WCAMAX, NOFOWE, WEXMIN, WEXMAX, ERRMAX )
C     NORME DE LA DIFFERENCE VG1-VG0 et NORME de VG0
      CALL NORMDIF( NBSOM, MCN(MNVG0), MCN(MNVG1),  NORDIFV, NORVECV )
C     NORME DE LA DIFFERENCE WG1-WG0 et NORME de WG0
      CALL NORMDIF( NBSOM, MCN(MNWG0), MCN(MNWG1),  NORDIFW, NORVECW )
      WRITE(IMPRIM,10000) 'EQ3:', TEMPS, ITERNL,
     %                     NORVECV, NORDIFV, NORDIFV/NORVECV,
     %                     NORVECW, NORDIFW, NORDIFW/NORVECW
C
C    -(Wm+2n+1-Wn)/dt +alfa LAPLAC Vm+2n+1
C             +beta*(Vm+2n+1**2+Wm+1n+1**2)Vm+2n+1=FOmegaReel => Wm+2n+1
C     c-a-d
C     tPP Wm+2n+1 = tPP Wn + dt alfa tDPDP Vm+2n+1 + dt beta tPU2P Vm+2n+1 - dt
      CALL NLSESCH( MCN(MNMGD), MCN(MNWG1), MCN(MNWG0),
     %              DTALFA, MCN(MNVG1),
     %              DTBETA, MCN(MNVG1), MCN(MNWG1), MCN(MNVG1), -DT,
     %              MCN(MNBG),
     %              NBNOFX, MCN(MNNOFX), MCN(MNVAFX), 2,
     %              NBSOM,  MCN(MNXYZN+WYZNOE), NBELEM, MCN(MNPGEL),
     %              NOOBSF, NUMIOB(3), NUMAOB(3), MCN(MNDOEL(3)), 1 )
C     AFFICHAGE DE LA SOLUTION AU TEMPS ACTUEL
      CALL AFNLSE( 20,  1, MNXYZN, NTDL/2, 1, MCN(MNVG1),
     %             WCAMIN, WCAMAX, NOFOWE, WEXMIN, WEXMAX, ERRMAX )
C     NORME DE LA DIFFERENCE VG1-VG0 et NORME de VG0
      CALL NORMDIF( NBSOM, MCN(MNVG0), MCN(MNVG1),  NORDIFV, NORVECV )
C     NORME DE LA DIFFERENCE WG1-WG0 et NORME de WG0
      CALL NORMDIF( NBSOM, MCN(MNWG0), MCN(MNWG1),  NORDIFW, NORVECW )
      WRITE(IMPRIM,10000) 'EQ4:', TEMPS, ITERNL,
     %                     NORVECV, NORDIFV, NORDIFV/NORVECV,
     %                     NORVECW, NORDIFW, NORDIFW/NORVECW
C
10000 FORMAT(/A,' TEMPS=',G14.6,' ITERNL=',I6,
     %T36,'NORM RP=',G14.6,' NORM DIF RP=',G14.6,'  DIF/RP=',G14.6/
     %T36,'NORM IP=',G14.6,' NORM DIF IP=',G14.6,'  DIF/IP=',G14.6)
C
      IF( NORDIFV .GT. NORVECV * 1D-3   .OR.
     %    NORDIFW .GT. NORVECW * 1D-3 ) GOTO 200
C
C     FIN DES ITERATIONS DE POINT FIXE
C     --------------------------------
C
C     STOCKAGE DE LA SOLUTION A CET INSTANT TEMPS?
      IF( TEMPS .GE. TSTOC*0.9999 .AND. NBWAVE .LT. MXVECT ) THEN
C
C       OUI  STOCKAGE DE LA SOLUTION A CET INSTANT TEMPS
         MNWAVE = MNWAVE + NTDL * MOREE2
         CALL TRTATD( MCN(MNUG1), MCN(MNWAVE), NTDL )
C
C        LE NOMBRE DE VECTEURS SOLUTION STOCKES
         NBWAVE = NBWAVE + 1
C        LE TEMPS DE STOCKAGE
         RMCN(MNTIME+NBWAVE) = TEMPS
C
C        MIN MAX DES SOLUTIONS A CE TEMPS
         CALL MX1VEC( NTDL, 1, 1, MCN(MNUG1),
     %                DWMIN, NOEMIN, DWMAX, NOEMAX )
C        DWMIN  : SOLUTION MINIMALE DU CAS NCAS REEL SIMPLE PRECISION
C        NOEMIN : NUMERO DU NOEUD OU LA SOLUTION EST MINIMALE
C        DWMAX  : SOLUTION MAXIMALE DU CAS NCAS
C        NOEMAX : NUMERO DU NOEUD OU LA SOLUTION EST MAXIMALE
         IF( NOEMIN .GT. NBSOM ) NOEMIN = NOEMIN - NBSOM
         MN = MNXYZN + WYZNOE + NOEMIN*3 -3
         XYZWMIN(1) = RMCN( MN )
         XYZWMIN(2) = RMCN( MN + 1 )
         XYZWMIN(3) = RMCN( MN + 2 )
C
         IF( NOEMAX .GT. NBSOM ) NOEMAX = NOEMAX - NBSOM
         MN = MNXYZN + WYZNOE + NOEMAX*3 -3
         XYZWMAX(1) = RMCN( MN )
         XYZWMAX(2) = RMCN( MN + 1 )
         XYZWMAX(3) = RMCN( MN + 2 )
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10150) TEMPS,   NTDL,  NBWAVE,
     %                   DWMIN, XYZWMIN, WDMAX, XYZWMAX
         ELSE
            WRITE(IMPRIM,20150) TEMPS,   NTDL,  NBWAVE,
     %                   DWMIN, XYZWMIN, WDMAX, XYZWMAX
         ENDIF
10150 FORMAT(' Au TEMPS',G14.6,' STOCKAGE DE',I9,' SOLUTIONS',
     %' en COLONNE',I5,' du TMS VECTEUR"ONDENLSE'/
     %' SOLUTION MIN=',G14.6,' AU NOEUD X=',G14.6,' Y=',G14.6,
     %' Z=',G14.6/
     %' SOLUTION MAX=',G14.6,' AU NOEUD X=',G14.6,' Y=',G14.6,
     %' Z=',G14.6)
20150 FORMAT(' At TIME',G14.6,' STORAGE of',I9,' SOLUTIONS',
     %' in COLUMN',I5,' of TMS VECTEUR"ONDENLSE'/
     %' SOLUTION MIN=',G14.6,' at NODE X=',G14.6,' Y=',G14.6,
     %' Z=',G14.6/
     %' SOLUTION MAX=',G14.6,' at NODE X=',G14.6,' Y=',G14.6,
     %' Z=',G14.6)
C
C        LE PROCHAIN TEMPS DE STOCKAGE
         TSTOC = TSTOC + DTSTOC
      ENDIF
C
C     MISE A JOUR DE MNUG0 PAR PERMUTATION AVEC MNUG1
C     ===============================================
      MN    = MNUG0
      MNUG0 = MNUG1
      MNUG1 = MN
C
C     PASSAGE AU PAS DE TEMPS SUIVANT?
      IF( TEMPS + DT .LT. TPSFIN*1.00001 ) GOTO 100
C
C    ##############################################################
C    ##                                                          ##
C    ##                FIN DE LA BOUCLE EN TEMPS                 ##
C    ##                                                          ##
C    ##############################################################
C
ccc      IF( RMCN(MNTIME+NBWAVE) .NE. TEMPS ) THEN
C     STOCKAGE DES SOLUTIONS A CET INSTANT
      IF( NBWAVE .LT. MXVECT ) THEN
         MNWAVE = MNWAVE + NTDL * MOREE2
C        LE NOMBRE DE VECTEURS SOLUTION STOCKES
         NBWAVE = NBWAVE + 1
      ELSE
C        LE VECTEUR SOLUTION MXVECT
         NBWAVE = MXVECT
      ENDIF
      CALL TRTATD( MCN(MNUG0), MCN(MNWAVE), NTDL )
C     LE TEMPS DE STOCKAGE
      RMCN(MNTIME+NBWAVE) = TEMPS
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10150) TEMPS, NBWAVE, NTDL
      ELSE
         WRITE(IMPRIM,20150) TEMPS, NBWAVE, NTDL
      ENDIF

ccc      ENDIF
C
C     MISE A JOUR DU TMS 'VECTEUR"ONDENLSE'
C     =====================================
      MCN( MNVECT + WBCOVE ) = NTDL
      MCN( MNVECT + WBVECT ) = NBWAVE
      MCN( MNVECT + WBCPIN ) = NBWAVE
      IF( NBWAVE .LT. MXVECT ) THEN
C        LE TMS EST RACOURCI
         L  = MNVECT + WECTEU + NTDL * MXVECT * MOREE2 - 1
         L1 = MNVECT + WECTEU + NTDL * NBWAVE * MOREE2 - 1
         DO I=1,NBWAVE
            RMCN(L1+I) = RMCN(L+I)
         ENDDO
         CALL TAMSRA( NTVECT, WECTEU+NTDL*NBWAVE*MOREE2+NBWAVE )
      ENDIF
C     LA DATE
      CALL ECDATE( MCN(MNVECT) )
C     LE NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNVECT + MOREE2 ) = NONMTD( '~>>>VECTEUR' )
C
C     LES SOLUTIONS A L'INSTANT INITIAL SONT A L'ADRESSE MNWAVE
      MNWAVE = MNVECT + WECTEU
C     L'ADRESSE -1 DU PREMIER TEMPS STOCKE DERRRIERE LES VECTEURS SOLUTION
      MNTIME = MNWAVE + NTDL * NBWAVE * MOREE2 -1
C
C     AFFICHAGE DES SOLUTIONS POUR TOUS LES TEMPS STOCKES
C     ===================================================
      CALL AFNLSE( 20,     NBWAVE, MNXYZN, NTDL/2, NBWAVE, MCN(MNWAVE),
     %             WCAMIN, WCAMAX, NOFOWE, WEXMIN, WEXMAX, ERRMAX )
C
C     DESTRUCTION DES TMC DEVENUS INUTILES
C     ====================================
 9999 IF( MNMGD  .GT. 0 ) CALL TNMCDS( 'REEL2',  NBSOM,  MNMGD  )
      IF( MNUG1  .GT. 0 ) CALL TNMCDS( 'REEL2',  NTDL,   MNUG1  )
      IF( MNUG0  .GT. 0 ) CALL TNMCDS( 'REEL2',  NTDL,   MNUG0  )
      IF( MNBG   .GT. 0 ) CALL TNMCDS( 'REEL2',  NBSOM,  MNBG   )
      IF( MNNOFX .GT. 0 ) CALL TNMCDS( 'ENTIER', MONOFX, MNNOFX )
      IF( MNVAFX .GT. 0 ) CALL TNMCDS( 'ENTIER', MOVAFX, MNVAFX )
C
      RETURN
      END
