      SUBROUTINE FL2NORMES( KNOMOB, IERR, DCPU )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LE VOLUME DU MAILLAGE
C -----    Volume = Som Jacobien(e) / (d+1)!
C                  e dans E
C          CALCULER LA NORME L2 DES NBVECT vecteurs PRESSION
C          norml2P = SQRT( Som  t{Pe} Int t[P] [P] dX {Pe} )
C                          e dans E
C          CALCULER LA NORME L2 DES NBVECT vecteurs |VITESSE|
C          norml2V = SQRT( Som      Som  t{Vei} Int t[P] [P] dX {Vei} )
C                          e dans E i=1...d
C
C          SI LA FONCTION VITESSE_EXACTE(t,x,y,z,nc) existe
C          CALCULER LA NORME L2 DES NBVECT ERREURS SUR LA PRESSION
C          norml2ErP = SQRT( Som  t{Pex-Pcal} Int t[P] [P] dX {Pex-Pcal} )
C                          e dans E
C
C          SI LA FONCTION PRESSION_EXACTE(t,x,y,z) existe
C          CALCULER LA NORME L2 DES NBVECT ERREURS SUR LE |VITESSE|
C          norml2ErV = SQRT( Som Som t{Viex-Vical} Int t[P] [P] dX {Viex-Vical}
C                          e dans E i=1...d
C
C          POUR DES ELEMENTS FINIS DE TAYLOR-HOOD ou BREZZI-FORTIN
C
C ENTREE :
C --------
C KNOMOB : NOM DE L'OBJET A CALCULER
C
C SORTIE :
C --------
C IERR   : 0 SI PAS D'ERREUR D'EXECUTION, NON NUL SINON
C DCPU   : SECONDES DE CPU DE L'EXECUTION DE CE SOUS PROGRAMME
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET TEXAS A & M UNIVERSITY at QATAR   Janvier 2012
C23456---------------------------------------------------------------012
      PARAMETER         (MXTYEL=7)
      PARAMETER         (MOPAGE=512)
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/donflu.inc"
      include"./incl/donele.inc"
      include"./incl/a_objet__topologie.inc"
      include"./incl/a_objet__definition.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___xyznoeud.inc"
      include"./incl/a___xyzpoint.inc"
      include"./incl/a___npef.inc"
      include"./incl/a___contact.inc"
      include"./incl/a___vecteur.inc"
      include"./incl/a___profil.inc"
      include"./incl/a___fluxpt.inc"
      include"./incl/a___tableau1r.inc"
      include"./incl/msvaau.inc"
      include"./incl/ponoel.inc"
      include"./incl/cthet.inc"
      include"./incl/ctemps.inc"
      include"./incl/cnonlin.inc"
C
      include"./incl/pp.inc"
      COMMON             MCN(MOTMCN)
      REAL              RMCN(1)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1),DMCN(1))
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
C
      EXTERNAL          ETTAEL
      DOUBLE PRECISION  RELMIN, VOLUME
      INTEGER           MNDOEL(4), MXDOEL(4)
      DOUBLE PRECISION  DINFO, DCPU
      LOGICAL           AVANT
C
      CHARACTER*(*)     KNOMOB
      DATA              RELMIN/-1D28/
C
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10000) KNOMOB
      ELSE
         WRITE(IMPRIM,20000) KNOMOB
      ENDIF
10000 FORMAT(/100('=')/
     %'Calcul de la NORME L2 de la VITESSE-PRESSION de l''OBJET',A/
     %100('='))
20000 FORMAT(/100('=')/
     %'Computation of L2-NORM of VELOCITY-PRESSURE VECTORS of the OBJECT
     %: ',A/100('='))
C
C     INITIALISATION DU TEMPS CALCUL INITIAL
      DCPU = DINFO( 'CPU' )
C
C     QUELQUES INITIALISATIONS
C     TEMPS DU CALCUL DU FLUIDE
      TEMPS  = 0.0
      IERR   = 0
C     LE NOMBRE DE MOTS D'UNE VARIABLE REEL2
      MOREE2 = MOTVAR(6)
C
C     PROTECTION DES ADRESSES POUR EVITER DES PROBLEMES LORS
C     DE LA DESTRUCTION DES TABLEAUX
      MNNPEF = 0
      DO I=1,4
         MNDOEL(I) = 0
         MXDOEL(I) = 0
      ENDDO
      NBTYEL = 0
      MOAUX  = 0
      MOTAEL = 0
      NBDLMX = 0
      NBNOVI = 0
      MOFLTO = 0
      MOFLPT = 0
      NBVECT = 1
      MNTAUX = 0
      MNINTVIT = 0
      MNINTPRE = 0
C
C     AFFICHAGE ET VERIFICATION DU NOM_DE_L'OBJET
C     ===========================================
C     RECHERCHE DU NOM DE L'OBJET DANS LE LEXIQUE DES OBJETS
      CALL LXLXOU( NTOBJE, KNOMOB, NTLXOB, MNLXOB )
C     S'IL N'EXISTE PAS RETOUR A LA DEMANDE DU NOM DE L'OBJET
      IF( NTLXOB .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'ERREUR: OBJET INCONNU ' // KNOMOB
         ELSE
            KERR(1) = 'ERROR: UNKNOWN OBJECT ' // KNOMOB
         ENDIF
         CALL LEREUR
         CALL LXIM( NTOBJE )
         IERR = 1
         RETURN
      ENDIF
C
C     RECHERCHE DU TABLEAU DEFINITION
      CALL LXTSOU( NTLXOB,'DEFINITION', NTDFOB, MNDFOB )
C     S'IL N'EXISTE PAS RETOUR A LA DEMANDE DU NOM DE L'OBJET
      IF( NTDFOB .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) ='ERREUR: DEFINITION INCONNUE de l''OBJET ' //KNOMOB
         ELSE
            KERR(1) ='ERROR: UNKNOWN DEFINITION for the OBJECT '//KNOMOB
         ENDIF
         CALL LEREUR
         IERR = 2
         RETURN
      ENDIF
C
C     RECHERCHE DU TMS  VECTEUR"VITESSEPRESSION DE L'OBJET
      CALL LXTSOU( NTLXOB, 'VECTEUR"VITESSEPRESSION',  NTVECT, MNVECT )
C     LES DL SONT PAR NOEUDS: VX,VY,VZ, evt PRESSION (aux sommets seulement)
C     ET TEMPS APRES TEMPS
      IF( NTVECT .LE. 0 ) THEN
         NBLGRC(NRERR) = 2
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'ERREUR: OBJET ' // KNOMOB
            KERR(2) = 'VITESSE PRESSION NON CALCULEES'
         ELSE
            KERR(1) = 'ERROR: OBJECT ' // KNOMOB
            KERR(2) = 'VELOCITY PRESSURE NOT COMPUTED'
         ENDIF
         CALL LEREUR
         IERR = 3
         GOTO 9999
      ENDIF
C
C     RECHERCHE DES TABLEAUX SOMMETS NOEUDS POINTS ASSOCIES A L'OBJET
C     ADRESSAGE DES ADRESSES DES TABLEAUX EF DE CET OBJET
      CALL TNMCDC( 'ENTIER', 2*MXTYEL, MNNPEF )
      MNTELE = MNNPEF + MXTYEL
      CALL NDPGEL( NTLXOB, NTTOPO, MNTOPO,
     %             NTXYZP, MNXYZP, NTXYZN, MNXYZN,
     %             NBTYEL, MCN(MNTELE), MCN(MNNPEF), IERR )
C     NTTOPO : NUMERO      DU TMS 'TOPOLOGIE' DE L'OBJET
C     MNTOPO : ADRESSE MCN DU TMS 'TOPOLOGIE' DE L'OBJET
C     NTXYZP : NUMERO      DU TMS 'XYZPOINT'  DE L'OBJET
C     MNXYZP : ADRESSE MCN DU TMS 'XYZPOINT'  DE L'OBJET
C     NTXYZN : NUMERO      DU TMS 'XYZNOEUD'  DE L'OBJET
C     MNXYZN : ADRESSE MCN DU TMS 'XYZNOEUD'  DE L'OBJET
C     NBTYEL : NOMBRE DE TYPES D'ELEMENTS FINIS DU MAILLAGE
C     NTELEM : NUMERO      DU TMS DES NBTYEL TYPES D'ELEMENTS FINIS
C     MNNPEF : ADRESSE MCN DU TMS DES NBTYEL TYPES D'ELEMENTS FINIS
      IF( IERR .NE. 0 ) GOTO 9999
C     L'ADRESSE DU TABLEAU NPEF"TYPE EF
      MNELE = MCN( MNNPEF )
C
C     LE NUMERO DU TYPE DE L'ELEMENT FINI
      NUTYEL = MCN( MNELE + WUTYEL )
C
C     LE NOMBRE D'ELEMENTS FINIS DE CE TYPE
      NBELEM = MCN( MNELE + WBELEM )
C
C     LE NOMBRE DE NOEUDS DE L'ELEMENT FINI
      NBNOEF = MCN( MNELE + WBNDEL )
C
C     LES NORMES L2 DES NBVECT PRESSION, |VITESSE|, VOLUME SONT ELLES
C     DEJA CALCULEES SUR LES VECTEUR"VITESSEPRESSION?
C     ===============================================================
      CALL LXTSOU( NTLXOB, 'VECTEUR"L2VITESSEPRESS', NTVEL2, MNVEL2 )
      IF( NTVEL2 .GT. 0 ) THEN
         IF( AVANT( MCN(MNVECT), MCN(MNVEL2) ) ) THEN
C           L2VITESSES-PRESSIONS POSTERIEURES AUX VITESSES-PRESSIONS
C           DEJA CALCULEES ET RECUPEREES SANS LES RECALCULER
            GOTO 500
         ENDIF
C        CE TABLEAU EST DETRUIT POUR ETRE RECALCULE
         CALL LXTSDS( NTLXOB, 'VECTEUR"L2VITESSEPRESS' )
      ENDIF
C
C     CONSTRUCTION DU TMS 'VECTEUR"L2VITESSEPRESS'
C     ============================================
      CALL TMSL2VP( KNOMOB, NTVEL2, MNVEL2,
     %              NTVIEX, MNVIEX, NTPREX, MNPREX, IERR )
      IF( IERR .NE. 0 .OR. NTVEL2 .LE. 0 ) GOTO 9999
C
C     AFFICHAGE DU TMS 'VECTEUR"L2VITESSEPRESS'
C     =========================================
C     NOMBRE DE NOEUDS SUPPORT DE LA VITESSE DU MAILLAGE
 500  NBNOEU = MCN(MNXYZN+WNBNOE)
C     NDIM DIMENSION EFFECTIVE DE L'ESPACE DES COORDONNEES
      CALL DIMCOO( NBNOEU, MCN(MNXYZN+WYZNOE), NDIM )
C
      NBVECT   = MCN( MNVEL2 + WBCOVE )
      NBCOL2   = MCN( MNVEL2 + WBVECT )
      MNTIMEL2 = MNVEL2 + WECTEU + MOREE2 * NBCOL2 * NBVECT
      VOLUME   = RMCN(MNTIMEL2+NBVECT)
      MNINTVIT = MNVEL2   + WECTEU
      MNINTPRE = MNINTVIT + MOREE2 * NBVECT
C
C     PRESENCE DE L'ERREUR L2 SUR |VITESSE| ?
      NOFOVI = NOFOVITE()
      IF( NOFOVI .GT. 0 ) THEN
         MNINTVER = MNINTPRE + MOREE2 * NBVECT
      ELSE
         MNINTVER = MNINTVIT
      ENDIF
C
C     PRESENCE DE L'ERREUR L2 SUR LA PRESSION ?
      NOFOPR = NOFOPRES()
      IF( NOFOPR .GT. 0 ) THEN
         MNINTPER = MNINTVER + MOREE2 * NBVECT
      ELSE
         MNINTPER = MNINTPRE
      ENDIF
C
C     AFFICHAGE DU TMS 'VECTEUR"L2VITESSEPRESS' des NORMES L2
      CALL AFL2NOVP( NDIM,   NBVECT,  MCN(MNTIMEL2),
     %               VOLUME, MCN(MNINTVIT), MCN(MNINTPRE),
     %               NOFOVI, MCN(MNINTVER), NOFOPR, MCN(MNINTPER) )
C
C     TRACER LES COURBES de la NORME L2 en FONCTION DU TEMPS
C     ||Vitesse||, PRESSION-PRESSION Min,
C     ||Vitesse Exacte-Calculee||, ||Vitesse Exacte-Calculee||/||Vitesse||
C      PRESSION Exacte-Calculee,   PRESSION Exacte-Calculee/PRESSION
C     ====================================================================
      IF( INTERA .GE. 1 ) THEN
         CALL TNMCDC( 'REEL2', NBVECT, MNTAUX )
         CALL TRL2NORM( NDIM, NBVECT, MCN(MNTIMEL2),  VOLUME,
     %                  MCN(MNINTVIT), MCN(MNINTPRE),
     %                  NOFOVI, MCN(MNINTVER), NOFOPR, MCN(MNINTPER),
     %                  MCN(MNTAUX)  )
         IF( MNTAUX .GT. 0 ) CALL TNMCDS( 'REEL2', NBVECT, MNTAUX )
      ENDIF
C
C     DESTRUCTION DES TABLEAUX TEMPORAIRES
9999  IF( MNNPEF .GT. 0 ) CALL TNMCDS( 'ENTIER',2*MXTYEL, MNNPEF )
      DO 11000 I=1,4
         IF( MNDOEL(I) .GT. 0 .AND. MXDOEL(I) .GT. 0 ) THEN
            CALL TNMCDS( 'ENTIER', MXDOEL(I), MNDOEL(I) )
         ENDIF
11000 CONTINUE
C
C     AFFICHAGE DU TEMPS CALCUL QUI PEUT ETRE IMPORTANT PAR L'EXECUTION
C     DES FONCTIONS EXACTES DE L'UTILISATEUR VITESSE_EXACTE, PRESSIONEXACTE
      DCPU = DINFO( 'DELTA CPU' )
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,*) 'TEMPS CALCUL des NORMES L2=',DCPU, ' secondes'
      ELSE
         WRITE(IMPRIM,*) 'CPU TIME to compute L2 NORMS=',DCPU,' seconds'
      ENDIF
C
      RETURN
      END
