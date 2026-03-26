      SUBROUTINE SUEXW3( NBS1, NBS2, NBS12, COSO,  MAT, A, B,
     S                   FCI,  FCJ,  ITRAV, TRAV,
     %                   COSOBO, ALPHA, DISTID, NITCOR, FADIST, IERR )
C***********************************************************************
C BUT :    MAILLER UN QUADRANGLE DANS UN PLAN DEFINI PAR SES BORDS
C          PAR LA METHODE DE WINSLOW
C***********************************************************************
C
C ENTREE:
C           NBS1   : NOMBRE DE POINTS SUR LIGNES "HORIZONTALES"
C           NBS2   : NOMBRE DE POINTS SUR LIGNES "VERTICALES"
C           NBS12  : NBS1 * NBS2 POUR LA DECLARATION STANDARD DES TABLEAUX
C           COSO   : COORDONNEES DES SOMMETS DU BORD DU MAILLAGE
C
C TRAVAIL:
C           COSO   : COORDONNEES DES NOEUDS DU MAILLAGE
C           MAT    : MATRICE DU SYSTEME
C           A      : MATRICE POUR LA RESOLUTION PAR BLOC
C           B      : SECOND MEMBRE
C           FCI    : VECTEUR DES FONCTIONS DE CONTROLE EN I
C           FCJ    : VECTEUR DES FONCTIONS DE CONTROLE EN J
C           ITRAV  : VECTEUR DE TRAVAIL ENTIERS
C           TRAV   : VECTEUR DE TRAVAIL REEELS
C           ALPHA  : VECTEUR DES ANGLES SUR LE BORD
C
C SORTIES:
C           COSO  : COORDONNEES DE TOUS LES SOMMETS DU MAILLAGE
C***********************************************************************
C AUTEUR: CHRISTOPHE DOURSAT ANALYSE NUMERIQUE UPMC NOVEMBRE 1988
C****+7**************************************************************012
C***********************************************************************
C             DECLARATION DES VARIABLES ET DES TABLEAUX
C                  ET INITIALISATION DES VARIABLES
C***********************************************************************
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU , IMPRIM , INTERA , NUNITE(29)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE       (MCN(1),RMCN(1))
C
CCC      DECLARATIONS  AVEC LES VALEURS EFFECTIVES (PB FORTRAN STANDARD)
CCC      DIMENSION         COSO(3,NBS1*NBS2),
CCC     %                  COSOB(3,NBS1*NBS2)
CCC      REAL              MAT(NBS1*NBS2,3)
CCC      DIMENSION         A(NBS1*NBS2,-1:1),
CCC     %                  B(NBS1*NBS2,2)
CCC      DIMENSION         FCI(NBS1*NBS2),
CCC     %                  FCJ(NBS1*NBS2)
CCC      DIMENSION         ITRAV(NBS1*NBS2),
CCC     %                  TRAV(NBS1*NBS2,2)
CCC      DIMENSION         COSOBO(3,2*(NBS1+NBS2)-3),
CCC     %                  ALPHA(2*(NBS1+NBS2)-3)
CCC      DIMENSION         DISTID(2*(NBS1+NBS2-2))
C
      DIMENSION         COSO(3,NBS12)
CCC      DIMENSION         COSOB(3,NBS12)
      REAL              MAT(NBS12,3)
      DIMENSION         A(NBS12,-1:1),
     %                  B(NBS12,2)
      DIMENSION         FCI(NBS12),
     %                  FCJ(NBS12)
      DIMENSION         ITRAV(NBS12),
     %                  TRAV(NBS12,2)
      DIMENSION         COSOBO(3,*),
     %                  ALPHA(*)
      DIMENSION         DISTID(*)
C
      DIMENSION         FADIST(4)
      INTEGER           SENSFR
      DOUBLE PRECISION  DINFO
C
      DOUBLE PRECISION  D2D3(3,3)
      REAL              POINTO(3),POINTI(3),POINTJ(3),PT3D(3),PT2D(3)
C
      IERR = 0
      PI = 4.E0*ATAN(1.E0)
      TESTL1 = 0.00001
      TESTL2 = 0.00001
      TESTL  = TESTL1
      NITC   = 1
C***********************************************************************
C                 VERIFICATION DU NOMBRE DE POINTS
C***********************************************************************
      IF ((NBS1.LE.2).OR.(NBS2.LE.2)) THEN
        IERR = 1
        NBLGRC(NRERR) = 2
        KERR(1) = 'ERREUR SUEXW3 :'
        KERR(2) = 'PAS DE NOEUD INTERIEUR'
        CALL LEREUR
        RETURN
      ENDIF
C***********************************************************************
C                        PASSAGE DANS LE PLAN
C***********************************************************************
      NEUO = 1
      NEUI = 2
      NEUJ = NBS1+1
      DO 10 K=1,3
        POINTO(K) = COSO(K,NEUO)
        POINTI(K) = COSO(K,NEUI)
        POINTJ(K) = COSO(K,NEUJ)
   10 CONTINUE
      CALL DF3D2D(POINTO,POINTI,POINTJ,D2D3,IERR)
      CALL VERPLA(POINTO,D2D3,NBS1,NBS2,COSO,IERR)
      IF (IERR.EQ.1) THEN
        NBLGRC(NRERR) = 2
        KERR(1) = 'ERREUR SUEXW3 :'
        KERR(2) = 'LA SURFACE EST GAUCHE '
        CALL LEREUR
        RETURN
      ELSE
        DO 16 J=1,NBS2
           DO 15 I=1,NBS1
              NEU = (J-1)*NBS1+I
              PT3D(1) = COSO(1,NEU)
              PT3D(2) = COSO(2,NEU)
              PT3D(3) = COSO(3,NEU)
              CALL CH3D3D(POINTO,D2D3,PT3D,PT2D)
              COSO(1,NEU) = PT2D(1)
              COSO(2,NEU) = PT2D(2)
              COSO(3,NEU) = PT2D(3)
 15        CONTINUE
 16     CONTINUE
      ENDIF
C***********************************************************************
C             REMPLISSAGE DE LA MATRICE COSOBO CONTENANT LES
C                   COORDONNEES DES NOEUDS DU BORD
C***********************************************************************
      DO 20 I=1,NBS1
        COSOBO(1,I) = COSO(1,I)
        COSOBO(2,I) = COSO(2,I)
        COSOBO(3,I) = 0.E0
   20 CONTINUE
      DO 30 J=2,NBS2-1
        NEUB = NBS1+J-1
        NEU  = J*NBS1
        COSOBO(1,NEUB) = COSO(1,NEU)
        COSOBO(2,NEUB) = COSO(2,NEU)
        COSOBO(3,NEUB) = 0.E0
   30 CONTINUE
      DO 40 I=NBS1,1,-1
        NEUB = 2*NBS1+NBS2-I-1
        NEU  = (NBS2-1)*NBS1+I
        COSOBO(1,NEUB) = COSO(1,NEU)
        COSOBO(2,NEUB) = COSO(2,NEU)
        COSOBO(3,NEUB) = 0.E0
   40 CONTINUE
      DO 50 J=NBS2-1,1,-1
        NEUB = 2*NBS1+2*NBS2-J-2
        NEU  = (J-1)*NBS1+1
        COSOBO(1,NEUB) = COSO(1,NEU)
        COSOBO(2,NEUB) = COSO(2,NEU)
        COSOBO(3,NEUB) = 0.E0
   50 CONTINUE
C***********************************************************************
C         CALCUL DES ANGLES ENTRE DEUX ARETES DU CONTOUR POUR :
C   1/ CALCULER LE SENS DE ROTATION DU CONTOUR
C        ET VERIFIER QU IL N EST PAS DEGENERE
C   2/ UTILISER CES ANGLES POUR AMELIORER "L'ORTHOGONALITE" AU BORD
C***********************************************************************
      CALL CALANG(NEUB,COSOBO,1,SENSFR,ALPHA)
      IF (SENSFR.EQ.0) THEN
        IERR = 1
        NBLGRC(NRERR) = 2
        KERR(1) = 'ERREUR SUEXW3 :'
        KERR(2) = 'BORD DEGENRE'
        CALL LEREUR
        RETURN
      ENDIF
      DO 70 J=1,NBS2
        DO 60 I=1,NBS1
          NEU = (J-1)*NBS1+I
          FCI(NEU) = 0.
          FCJ(NEU) = 0.
   60   CONTINUE
   70 CONTINUE
CC
C      WRITE (*,*) " "
C      WRITE (*,*) "INVERSION DES FONCTIONS DE CONTROLE ?"
C      WRITE (*,*) "                 NON ........ 0 "
C      WRITE (*,*) "                 OUI ........ 1 "
C      WRITE (*,*) " "
C      READ (*,*) INV
C      WRITE (*,*) " "
C      IF (INV.EQ.1) THEN
C        CALL CALINV(NBS1,NBS2,NBS12,COSO,FCI,FCJ)
C      ENDIF
CC
C      WRITE (*,*) " "
C      WRITE (*,*) "TEST DE SORTIE : NORME 1 .... 1 "
C      WRITE (*,*) "                 NORME SUP .. 2 "
C      WRITE (*,*) " "
C      READ (*,*) NCHT
C      WRITE (*,*) " "
CC
      NCHT = 1
C
      TCPU = REAL(DINFO('DELTA CPU'))
C***********************************************************************
C            CALCUL DU OMEGA OPTIMAL DE LA SURRELAXATION
C             ON PREND CELUI CORRESPONDANT AU LAPLACIEN
C***********************************************************************
      PSXX2 =  COS(PI/(NBS2-1.))/(2-COS(PI/(NBS1-1.)))
      OMEGA = 2.E0/(1E0+SQRT(1.E0-PSXX2**2))
C***********************************************************************
C*                                                                     *
C*                      BOUCLES DES ITERATIONS                         *
C*                                                                     *
C***********************************************************************
      NBITERM = 0
      DO 80 NBITER=1,800
CCC        WRITE (IMPRIM,1000) NBITER
CCC 1000   FORMAT('PROBLEME LINEARISE ',T45,':',I4)
        CALL RESSUR(NBS1,NBS2,COSO,FCI,FCJ,
     S              MAT,A,B,TRAV,OMEGA,TEST,TESTL,NCHT)
        IF (TEST.LE.TESTL) THEN
          WRITE (IMPRIM,1060) NBITER-NBITERM
 1060     FORMAT(' NOMBRE DE PROBLEMES LINEARISES RESOLUS',T45,':',I5)
          NBITERM = NBITER
C***********************************************************************
C*       LA PRECISION DEMANDEE POUR LE CALCUL EST ATTEINTE             *
C*             ON PASSE A DEUX ETAPES SUCCESSIVES                      *
C* 1 VERIFICATION DE L INJECTIVITE DU MAILL., CORRECTION SI NECESSAIRE *
C* 2 AMELIORATION DE LA GEOMETRIE DES MAILLES FRONTIERES               *
C***********************************************************************
C*                      ETAPE 1: INJECTIVITE                           *
C***********************************************************************
C RECHERCHE DES POINTS EXTERIEURS AU MAILLAGE ET DES MAILLES RETOURNEES
C***********************************************************************
          NBELTR = 0
          CALL VERGEO(NBS1,NBS2,COSO,SENSFR,NBELTR,NBELTD,
     S                ITRAV)
          NBELTR = NBELTR+NBELTD
          IF (NBELTR.NE.0) THEN
            CALL VERINJ(NBS1,NBS2,COSO,COSOBO,
     S                  NBPEXT,NBPFRT,ITRAV(NBELTR+1))
            NBPEXT = NBPEXT+NBPFRT
            CALL RESBIJ(NBS1,NBS2,FCI,FCJ,NBELTR,NBPEXT,ITRAV)
            GOTO 80
          ENDIF
C***********************************************************************
C*                ETAPE 2: GEOMETRIE DES MAILLES FRONTIERES            *
C***********************************************************************
          IF (NITC.LE.NITCOR) THEN
            WRITE (IMPRIM,1061) NITC
 1061       FORMAT(' ITERATION DE CORRECTION DE L''ORTHOGONALITE',T45,
     S             ':',I5)
            CALL AMELMA2D(NBS1,NBS2,NITC,COSO,ALPHA,
     S                    SENSFR,FCI,FCJ,DISTID,FADIST,  TRAV )
            NITC = NITC+1
            GOTO 80
          ENDIF
C***********************************************************************
C                        AIGUILLAGE DE SORTIE
C             SI LE TEST EST A LA VALEUR DE TRAVAIL OU
C     SI IL RESTE DES CORRECTIONS SUR LES MAILLES FRONTIERES
C                        ON RENVOIE AU CALCUL
C              SOIT POUR AVOIR UNE MEILLEURE PRECISION
C   SOIT POUR EFFECTUER LES MODIFICATIONS DUES AUX FONC. DE CONT.
C***********************************************************************
          IF (TESTL.GT.TESTL2) THEN
C-----------------------------------------------------------------------
C NOUS NE SOMMES PAS A LA PRECISION FINALE ON VA DONC RENVOYER AU CALCUL
C-----------------------------------------------------------------------
            TESTL = TESTL2
            GOTO 80
          ELSE
C-----------------------------------------------------------------------
C LA PRECISION FINALE EST ATTEINTE ON A DONC FINI LE CALCUL DU MAILLAGE
C                     IMPRESSION DE RESULTATS
C-----------------------------------------------------------------------
            TCPU = REAL(DINFO('DELTA CPU'))
            IF ((NBELTR+NBELTD).EQ.0) THEN
              WRITE (IMPRIM,1070) NBITER
              WRITE (IMPRIM,1072) NINT(TCPU)
 1070         FORMAT(' CONVERGENCE EN ',I4,' ITERATIONS')
 1072         FORMAT(' TEMPS CPU :',I12,' SECONDES')
            ELSE
              WRITE (IMPRIM,1010)
              WRITE (IMPRIM,1036)
              WRITE (IMPRIM,1037)
              WRITE (IMPRIM,1038)
              WRITE (IMPRIM,1010)
            ENDIF
CCC            WRITE (IMPRIM,1033)
CCC            WRITE (IMPRIM,1030) NBITER
CCC            WRITE (IMPRIM,1033)
CCC            WRITE (IMPRIM,1050)
CCC            WRITE (IMPRIM,1051) TCPU
CCC            WRITE (IMPRIM,1010)
 1010       FORMAT('**************************************')
 1011       FORMAT('*             ATTENTION              *')
CCC 1030       FORMAT('*   CONVERGENCE EN',I4,' ITERATIONS    *')
 1033       FORMAT('*                                    *')
 1034       FORMAT('*  MAILLAGE CORRECT TOPOLOGIQUEMENT  *')
 1035       FORMAT('*        BON POUR LE SERVICE         *')
 1036       FORMAT('*         MAUVAIS MAILLAGE           *')
 1037       FORMAT('*    DES ELEMENTS SONT RETOURNES     *')
 1038       FORMAT('*           OU DEGENERES             *')
CCC 1050       FORMAT('*         TEMPS CPU UTILISE          *')
CCC 1051       FORMAT('*           ',1PE11.3,'              *')
            GOTO 90
          ENDIF
        ENDIF
   80 CONTINUE
C***********************************************************************
C   TOUTES LES BOUCLES ONT ETE EFFECTUEES DONC LA METHODE N A PAS
C      CONVERGE EN 100 ITERATIONS. SORTIE AVEC MESSAGE D ERREUR
C***********************************************************************
      CALL VERINJ( NBS1 , NBS2 , COSO , COSOBO , NBPEXT , NBPFRT, ITRAV)
      NBELTR = -1
      CALL VERGEO( NBS1 , NBS2 , COSO , SENSFR , NBELTR , NBELTD, ITRAV)
      WRITE (IMPRIM,1010)
      WRITE (IMPRIM,1011)
      WRITE (IMPRIM,1040)
      WRITE (IMPRIM,1041)
      WRITE (IMPRIM,1042)
      WRITE (IMPRIM,1033)
      NBLGRC(NRERR) = 3
      KERR(1) = 'ATTENTION SUEXW3 :'
      KERR(2) = 'NON CONVERGENCE DU CALCUL'
      IF ((NBPEXT+NBPFRT+NBELTR+NBELTD).EQ.0) THEN
        WRITE (IMPRIM,1034)
        WRITE (IMPRIM,1035)
        KERR(3) = 'MAILLAGE CORRECT'
      ELSE
        WRITE (IMPRIM,1036)
        KERR(3) = 'MAILLAGE DEGENERE'
      ENDIF
      CALL LEREUR
      WRITE (IMPRIM,1010)
 1040 FORMAT('*  PROBLEME DANS LA ROUTINE SUEXW3   *')
 1041 FORMAT('*   TOUTES LES ITERATIONS ONT ETE    *')
 1042 FORMAT('*    EFFECTUEES SANS CONVERGENCE     *')
C***********************************************************************
C                          RETOUR DANS L'ESPACE
C***********************************************************************
   90 CONTINUE
      DO 130 J=1,NBS2
        DO 120 I=1,NBS1
          NEU = (J-1)*NBS1+I
          DO 100 K=1,2
            PT2D(K) = COSO(K,NEU)
  100     CONTINUE
          CALL CH2D3D(POINTO,D2D3,PT2D,PT3D)
          DO 110 K=1,3
            COSO(K,NEU) = PT3D(K)
  110     CONTINUE
  120   CONTINUE
  130 CONTINUE
      RETURN
      END
