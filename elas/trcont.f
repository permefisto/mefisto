      SUBROUTINE TRCONT( KNOMOB, NTLXOB, MODECO, NDIM,
     %                   NBTYEL, MNTOPO, MNNPEF, MNPOGE,
     %                   NBVECT, MNTIME,
     %                   NCAS  , NOPT  , CMFLEC, CMPCON )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACER LES FLECHES REPRESENTANT LES CONTRAINTES ET
C -----    DIRECTIONS PRINCIPALES D'UN OBJET
C
C ENTREES :
C ---------
C KNOMOB : NOM DE L'OBJET
C NTLXOB : NUMERO DU TMS LEXIQUE DE L'OBJET A TRAITER
C MODECO : MODE DE TRACE DES VECTEURS DU TMS D'ADRESSE MNDEPL
C          =1 CE SONT DES DEPLACEMENTS
C          =2 CE SONT DES MODES PROPRES
C NDIM   : DIMENSION DE L'ESPACE =2 OU 3
C NBTYEL : LE NOMBRE DE TYPES D'ELEMENTS DANS CETTE TOPOLOGIE
C MNTOPO : ADRESSE MCN DU TABLEAU TOPOLOGIE
C MNNPEF : ADRESSE MCN DES TABLEAUX ELEMENTS FINIS DE CETTE TOPOLOGIE
C MNPOGE : ADRESSE MCN DU TABLEAU POINTS GEOMETRIQUES DE L'OBJET
C NBVECT : NOMBRE DE VECTEURS DEPLACEMENTS
C MNTIME : ADRESSE MCN DU TABLEAU DES TEMPS OU LES DEPLACEMENTS
C          ONT ETE CALCULES
C          0 SI PAS DE STOCKAGE
C
C MODIFIES EVENTUELLEMENT :
C -------------------------
C NCAS   : NUMERO DU CAS A TRAITER
C NOPT   : NUMERO DE L'OPTION DE TRAITEMENT
C CMFLEC : NOMBRE DE CM DU DE LA FLECHE MAXIMALE
C CMPCON : NOMBRE DE CM POUR L'UNITE DE CONTRAINTE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 1989
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/gsmenu.inc"
      include"./incl/a_objet__topologie.inc"
      include"./incl/a___npef.inc"
      include"./incl/a___contrainte.inc"
      include"./incl/a___xyzpoint.inc"
      include"./incl/ponoel.inc"
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      include"./incl/ctemps.inc"
      include"./incl/xyzext.inc"
      include"./incl/a___face.inc"
      include"./incl/a___aretefr.inc"
C
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      REAL             RMCN(1)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      EQUIVALENCE      (MCN(1),RMCN(1))
C
      REAL              CONTMX,HEXSEC(6,2)
      CHARACTER*(*)     KNOMOB
C
      IERR = 0
      MOREE2 = MOTVAR(6)
C
C     LIMITATION OU NON PAR LA FONCTION 'REGION(t,x,y,z)'
      CALL LXNMNO( NTFONC, 'REGION', NOFORE, I )
C     NOFORE>0 SI CETTE FONCTION EXISTE
C
C     RECHERCHE DU MIN ET MAX DES COORDONNEES DES POINTS
C     RECHERCHE DU MAXIMUM EN VALEUR ABSOLUE DES CONTRAINTES
 10   CALL MXCONT( NTLXOB, NBTYEL, MNNPEF,
     %             HEXSEC, CONTMX, IERR )
      IF( IERR .NE. 0 ) RETURN
      IF( CONTMX .LE. 0 ) CONTMX = MAX( CONTMX, 1.0 )
      CMPCON = CMFLEC / CONTMX
      IF( CMPCON .EQ. 0 ) CMPCON = 1.0
C
C     LE CADRE ENGLOBANT
      DO 14 I=1,3
         DO 12 J=1,2
            COOEXT(I,J) = HEXSEC(I,J)
 12      CONTINUE
 14   CONTINUE
C
C     LECTURE DES DONNEES DU TRACE DES CONTRAINTES
C     ============================================
 100  CALL LIMTCL( 'traccont', NMTCL )
      IF( NMTCL .LE. 0  ) GOTO 9000
      IF( NMTCL .EQ. 50 ) GOTO 250
C
      GOTO( 110, 120, 130, 140, 150, 160, 170, 180, 190, 300,
     %      100, 100, 100, 100, 300, 270 ), NMTCL
C
C     NUMERO DU CAS A VISUALISER
 110  NCVALS = 4
      CALL INVITE( 84 )
      CALL LIRENT( NCVALS, NCAS )
      IF( NCVALS .EQ. -1 ) GOTO 100
C     PROTECTION DU NUMERO DE CAS A TRACER
      NCAS = MAX( 1, ABS(NCAS) )
      NCAS = MIN( NBVECT, NCAS )
C     LES TEMPS ONT ILS ETE STOCKES?
      IF( MNTIME .GT. 0 ) THEN
C        OUI: LE TEMPS INITIAL EST CELUI DU VECTEUR"DEPLACT
         TEMPS  = RMCN( MNTIME + NCAS )
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,*) 'TRACE au TEMPS ',TEMPS
         ELSE
            WRITE(IMPRIM,*) 'DRAWING at TIME ',TEMPS
         ENDIF
      ELSE
C        TEMPS INITIAL SUPPOSE NUL
         TEMPS = 0
      ENDIF
      GOTO 100
C
C     FLECHE MAXIMALE TRACEE EN CM
 120  NOPT   = 1
      NCVALS = 0
      CALL INVITE( 66 )
      CALL LIRRSP(NCVALS,CMFLEC)
      IF( NCVALS .EQ. -1 ) GOTO 100
      IF( CMFLEC .LT. 0. ) CMFLEC = -CMFLEC
      GOTO 100
C
C     1CM  VAUT EN CONTRAINTE
 130  NOPT   = 2
      NCVALS = 0
      CALL INVITE( 97 )
      CALL LIRRSP(NCVALS,CONTMX)
      IF( NCVALS .EQ. -1 ) GOTO 100
      IF( CONTMX .LT. 0. ) CONTMX = -CONTMX
      IF( CONTMX .EQ. 0. ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'DONNEE de CONTRAINTE MAX NULLE. A CORRIGER'
         ELSE
            KERR(1) = 'INPUT of NULL MAXIMUM STRESS. MODIFY IT'
         ENDIF
         CALL LEREUR
         GOTO 10
      ENDIF
      GOTO 100
C
C     ECHELLE PRECEDENTE
 140  IF( CMPCON .LE. 0. ) THEN
         WRITE(KERR(MXLGER)(1:13),'(G13.5)') CMPCON
         NBLGRC(NRERR) = 2
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'ERREUR: 1CM/CONTRAINTE='// KERR(MXLGER)(1:13)
            KERR(2) = 'ECHELLE INTERDITE :'
     %              // ' LA FLECHE MAXIMALE EST FORCEE A 2.5CM'
         ELSE
            KERR(1) = 'ERROR: 1CM/STRESS='// KERR(MXLGER)(1:13)
            KERR(2) = 'INCORRECT VALUE :'
     %              // ' The MAXIMUM ARROW is FORCED to 2.5CM'
         ENDIF
         CALL LEREUR
         NOPT   = 1
         CMFLEC = 2.5
      ELSE
         NOPT   = 3
      ENDIF
      GOTO 100
C
C     COULEUR des ARETES du MAILLAGE
 150  CALL LIMTCL( 'couleur0' , I )
      IF( I .EQ. -1 ) THEN
         GOTO 9000
      ELSE IF( I .EQ. -2 ) THEN
         NCOUAF = -2
      ELSE IF( I .EQ. 0 ) THEN
C        COULEUR NOIRE
         NCOUAF = 0
      ELSE
C        COULEUR RESERVEE
         NCOUAF = N1COEL + I
      ENDIF
      GOTO 100
C
C     TYPE du TRAIT des ARETES du MAILLAGE
 160  CALL LIMTCL( 'typtrait' , I )
      IF( I .EQ. -1 ) GOTO 100
      NTLAFR = I
      GOTO 100
C
C     COULEUR des FLECHES
 170  CALL LIMTCL( 'couleurs' , I )
      IF( I .LT. 0 ) THEN
         GOTO 9000
      ELSE IF( I .EQ. 0 ) THEN
C        COULEUR NOIRE
         NCOUFL = 0
      ELSE
         NCOUFL = N1COEL + I
      ENDIF
      GOTO 100
C
C     NOMBRE D'EPAISSEURS DE TRAIT D'UNE FLECHE
 180  NCVALS = 0
      CALL INVITE( 78 )
      CALL LIRENT( NCVALS, NEPFLE )
      IF( NCVALS .EQ. -1 ) GOTO 100
      IF( NEPFLE .LT. 0 ) NEPFLE=0
      IF( NEPFLE .GT. 5 ) NEPFLE=5
      GOTO 100
C
C     AFFICHAGE DES CONTRAINTES EVENTUELLEMENT DANS UNE REGION
 190  GOTO 300
C
C     EFFACER LE TRACE ACTUEL
 250  CALL EFFACE
C     PLUS D'ITEMS VISIBLES
      CALL ITEMS0
      CALL TRAXES
      GOTO 100
C
C     TRACE des CONTRAINTES dans 1 SECTION PLANE des EF 3D
C     ----------------------------------------------------
 270  IF( NDIM .EQ. 2 ) GOTO 300
      CALL TRCOEFSE( KNOMOB, NTLXOB,
     %               MCN(MNPOGE+WBCOOP), MCN(MNPOGE+WNBPOI),
     %               MCN(MNPOGE+WYZPOI), NOFORE,
     %               MODECO, NCAS, NOPT, CONTMX, CMFLEC, CMPCON )
      GOTO 100
C
C     TRACE des CONTRAINTES dans tous les EF
C     --------------------------------------
 300  CALL TRCOEF( NMTCL,  NOFORE, KNOMOB, NTLXOB, MODECO, NDIM,
     %             NBTYEL, MNTOPO, MNNPEF, MNPOGE,
     %             NCAS,   NOPT,   CONTMX, CMFLEC, CMPCON )
      GOTO 100
C
9000  RETURN
      END
