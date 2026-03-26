      SUBROUTINE DVAJPC( ZEMEPS, UNPEPS,
     %                   NBSOMT, NOPXYD, MXTRIA, PXYD,
     %                   N1TRVI, NOTRIA, CETRIA, NOTRSO,
     %                   NBLFTR, NDARLF, MXSOAR, NOSOAR, NLSOFR,
     %                   MXETOI, NAETOI, NTETOI,
     %                   IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    GENERER LES TRIANGLES DELAUNAY DES SOMMETS 1 A NBSOMT
C -----    DANS L'ORDRE DU TABLEAU NOPXYD QUI POINTE SUR LE TABLEAU PXYD
C          AVEC LA CONTRAINTE DE NE PAS PERDRE D'ARETE FRONTALIERE
C          SI UN POINT PERD UNE ARETE CE POINT EST SUPPRIME
C
C ENTREES:
C --------
C ZEMEPS :  -EPSILON POUR LA COORDONNEE BARYCENTRIQUE
C UNPEPS : 1+EPSILON POUR LA COORDONNEE BARYCENTRIQUE
C NOSTBA : NUMERO PXYD DU BARYCENTRE DE L'ENVELOPPE CONVEXE
C NBSOMT : NOMBRE DE SOMMETS A TRIANGULER
C NOPXYD : NUMERO DANS PXYD DES SOMMETS A TRIANGULER
C MXTRIA : NOMBRE MAXIMAL DE TRIANGLES DECLARABLES
C PXYD   : TABLEAU DES COORDONNEES 2D DES POINTS
C NBLFTR : NOMBRE DE LIGNES FERMEES LIMITANT LA SURFACE A TRIANGULER
C NDARLF : NUMERO DE LA PREMIERE ARETE DE CHAQUE LIGNE FERMEE DANS
C          LE TABLEAU NOSOAR
C NBARFR : NOMBRE D'ARETES  DE LA FRONTIERE OU CONTOUR DU DOMAINE =
C          NOMBRE DE POINTS DE LA FRONTIERE OU CONTOUR DU DOMAINE
C MXSOAR : NOMBRE MAXIMAL D'ARETES FRONTALIERES DECLARABLES
C NOSOAR : NUMERO DES 2 SOMMETS DE CHAQUE ARETE ET
C          POINTEUR SUR L'ARETE SUIVANTE/N1AEVI/
C NLSOFR : TABLEAU DU NUMERO DE LIGNE(1 A NBLFTR) DES SOMMETS
C         -NUMERO DE POINT INTERNE UTILISATEUR IMPOSE, 0 SINON
C
C ENTREES ET SORTIES :
C --------------------
C N1TRVI : NUMERO DU 1 PREMIER TRIANGLE VIDE DANS LE TABLEAU NOTRIA
C          LE CHAINAGE DES TRIANGLES VIDES SE FAIT SUR NOTRIA(4,.)
C NOTRIA : LISTE DES TRIANGLES
C                 ------- ------- ------- -------- -------- --------
C  PAR TRIANGLE : SOMMET1 SOMMET2 SOMMET3 TR_VOIS1 TR_VOIS2 TR_VOIS3
C                 ------- ------- ------- -------- -------- --------
C                 SOMMET    EST LE NUMERO DU SOMMET
C                 TR_VOIS i EST LE NUMERO DANS NOTRIA DU TRIANGLE
C                               ADJACENT PAR L'ARETE i
C
C CETRIA : COORDONNEES DU CENTRE DU CERCLE CIRCONSCRIT ET
C          CARRE DU RAYON
C                         ------- ------- --------
C          PAR TRIANGLE : XCENTRE YCENTRE RAYON**2
C                         ------- ------- --------
C          TABLEAU REEL(3,MXTRIA)
C
C NOTRSO : NOTRSO(I) NUMERO D'UN TRIANGLE AYANT POUR SOMMET I
C
C SORTIES:
C --------
C IERR   : 0 SI PAS D'ERREUR
C          1 SATURATION DES SOMMETS
C          3 SATURATION DES TRIANGLES
C          8 AUCUN CERCLE CIRCONSCRIT AUX TRIANGLES ACTUELS NE CONTIENT
C            LE POINT COURANT
C          9 SATURATION DE NTETOI(MXETOI) => AUGMENTER MXETOI
C         10 SATURATION DE NAETOI(1:4,1:MXETOI) => AUGMENTER MXETOI
C
C TABLEAUX AUXILIAIRES :
C ----------------------
C MXETOI : NOMBRE MAXIMAL D'ARETES OU TRIANGLES DECLARABLES
C          DANS UNE ETOILE
C NAETOI : ENTIER (1:4,1:MXETOI)
C NTETOI : ENTIER (1:MXETOI)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC    JUILLET 1994
C....................................................................012
      include"./incl/trvari.inc"
      include"./incl/gsmenu.inc"
      PARAMETER        (QUADEG=0.01)
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      DOUBLE PRECISION  ZEMEPS, UNPEPS
      INTEGER           NOPXYD(0:*),
     %                  NOTRIA(6,MXTRIA),
     %                  NOTRSO(*),
     %                  NAETOI(4,MXETOI),
     %                  NLSOFR(*),
     %                  NDARLF(NBLFTR),
     %                  NOSOAR(3,MXSOAR),
     %                  NTETOI(MXETOI)
      DOUBLE PRECISION  PXYD(3,*),
     %                  CETRIA(3,MXTRIA)
      DOUBLE PRECISION  CB(3), CB1, CB2, CB3, CMIN
      EQUIVALENCE      (CB(1),CB1),(CB(2),CB2),(CB(3),CB3)
C
      NCOUL  = NDCORE
C
C     REINITIALISATION A VIDE DES ARETES DE L'ETOILE
      N1AEVI = 1
      N1AEOC = 0
      DO 10 I=1,MXETOI
C        NUMERO DU TRIANGLE DE L'AUTRE COTE DE L'ARETE
C         NAETOI(1,I) = 0
C        NUMERO LOCAL AU TRIANGLE DE L'ARETE
C         NAETOI(2,I) = 0
C        NUMERO DANS NAETOI DE L'ARETE SUIVANTE
         NAETOI(4,I) = I+1
  10  CONTINUE
      NAETOI(4,MXETOI) = 0
C
      NP = 0
C
C     TRIANGULATION DU SOMMET NOPXYD(NP+1)
C     ====================================
 100  NP = NP + 1
      IF( NP .LE. NBSOMT ) THEN
C
C        LE NUMERO DANS PXYD DU POINT A TRIANGULER
         N = NOPXYD( NP )
CCC         WRITE(IMPRIM,*) 'SOMMET NP=',NP,' N=',N
         IF( N .LE. 0 ) THEN
            WRITE(IMPRIM,*) 'DVAJPC: ANOMALIE N<=0 N=',N
            CALL XVPAUSE
            GOTO 100
         ENDIF
C
C        REINITIALISATION A VIDE DES ARETES ETOILANT LE POINT N-1
C        --------------------------------------------------------
         NA1 = 0
         NA2 = N1AEOC
 140     IF( NA2 .GT. 0 ) THEN
            NA1 = NA2
            NA2 = NAETOI(4,NA2)
            GOTO 140
         ENDIF
C        FIN DU CHAINAGE DES ARETES
         IF( NA1 .GT. 0 ) THEN
C           LA DERNIERE ARETE OCCUPEE EST CHAINEE SUR LA 1-ERE VIDE
            NAETOI(4,NA1) = N1AEVI
C           LA 1-ERE OCCUPEE DEVIENT LA 1-ERE VIDE
            N1AEVI = N1AEOC
         ENDIF
C
C        IL N'Y A PAS D'ARETES DANS L'ETOILE
         N1AEOC = 0
C
C        LE TABLEAU DES TRIANGLES A DETRUIRE EST VIDE
         NBTRET = 0
C
C        RECHERCHE DU 1-ER TRIANGLE NOTRI1 CONTENANT LE POINT N
C        ======================================================
C        LE 1-ER TRIANGLE A CONSULTER  ( DERNIER CONSTRUIT )
         I = NP
C
 160     I = I - 1
         IF( I .GT. 0 ) THEN
            J      = NOPXYD( I )
            NOTRI1 = NOTRSO( J )
            IF( NOTRI1 .LE. 0 ) GOTO 160
            IF( NOTRIA(1,NOTRI1) .LE. 0 ) GOTO 160
         ELSE
            DO 162 I=1,MXTRIA
               IF( NOTRIA(1,I) .GT. 0 ) THEN
C                 UN TRIANGLE ACTIF DANS L'ENVELOPPE CONVEXE
                  NOTRI1 = I
                  GOTO 165
               ENDIF
 162        CONTINUE
            WRITE(IMPRIM,*) 'DVAJPC: ANOMALIE PT=',N,' DANS 0 TRIANGLE'
            CALL XVPAUSE
         ENDIF
C
C        ***********************************************
C        LE TRIANGLE NOTRI1 CONTIENT IL LE POINT N ?
C        ***********************************************
 165     CALL TRPTDT( N, PXYD, NOTRIA, ZEMEPS, UNPEPS,
     %                NOTRI1, CB1, CB2, CB3, IMIN )
 167     IF( NOTRI1 .LE. 0 ) THEN
C
C           AUCUN TRIANGLE NE CONTIENT CE POINT
            IF( NLSOFR(N) .EQ. 0 ) THEN
C              CE POINT N EST SUPPRIMABLE
C              AUCUN TRIANGLE NE CONTIENT CE POINT : CE POINT EST SUPPRIME
CCC               WRITE(IMPRIM,*)'DVAJPC 1: NP=',NP,' N=',N,' SUPPRIME X=',
CCC     %                         PXYD(1,N),' Y=',PXYD(2,N)
               GOTO 100
            ELSE
C
C              CE POINT N NE PEUT ETRE SUPPRIME CAR IMPOSE
               IF( NOTRI1 .EQ. 0 ) THEN
                  WRITE(IMPRIM,*)'DVAJPC 2: PB PAS DE TRIANGLE OPPOSE ',
     %                            NOTRI1
                  CALL XVPAUSE
                  GOTO 100
               ENDIF
C
C              JONCTION DU POINT AVEC L'ARETE FRONTIERE DU
C              TRIANGLE -(NOTRI1 <0 EN RETOUR)
               NOTRI1 = -NOTRI1
CCC               WRITE(IMPRIM,*) 'DVAJPC 3: TRIANGLE SORTIE ', NOTRI1,
CCC     %         ' PAR LE COTE ',IMIN
C
C              CREATION DU TRIANGLE N, ARETE DE SORTIE DU TRIANGLE NOTRI1
               CALL DV1P1T( N, PXYD, NOTRI1, IMIN,
     %                      N1TRVI, NOTRIA, CETRIA, NOTRSO,
     %                      IERR )
               IF( IERR .EQ. 0 ) THEN
C                 TRAITEMENT TERMINE
                  GOTO 100
               ELSE
C                 ERREUR => SORTIE
                  NBLGRC(NRERR) = 2
                  WRITE(KERR(MXLGER)(1:9),'(I9)') N
                  KERR(1) = 'IMPOSSIBLE AJOUTER LE POINT' //
     %                       KERR(MXLGER)(1:9)
                  KERR(2) = 'MODIFIER LA LIGNE AUTOUR DU POINT'
                  CALL LEREUR
                  GOTO 9999
               ENDIF
            ENDIF
         ENDIF
C
C        ICI CE TRIANGLE NOTRI1 CONTIENT LE POINT N
         NOTRI0 = NOTRI1
C
C        RECHERCHE DE L'ARETE DE NOTRI0 LA PLUS PROCHE DU POINT N
         CMIN = CB(1)
         IMIN = 1
         DO 170 I=2,3
            IF( CB(I) .LT. CMIN ) THEN
               IMIN = I
               CMIN = CB(I)
            ENDIF
 170     CONTINUE
C        CALCUL DE L'ARETE OPPOSEE AU SOMMET IMIN
         IF( IMIN .EQ. 3 ) THEN
            IMIN = 1
         ELSE
            IMIN = IMIN + 1
         ENDIF
C        SI ICI CMIN .LT. 0.01 LE POINT N EST PROCHE DU COTE IMIN DU TRIANGLE NO
         IF( NOTRIA(3+IMIN,NOTRI0) .EQ. 0 .AND. ABS(CMIN) .LT. 0.001D0
     %       .AND. NLSOFR(N) .EQ. 0 ) THEN
C           POINT N NON IMPOSE ET TROP PRES D'UNE ARETE FRONTALIERE
CCC            WRITE(IMPRIM,*)'DVAJPC  : CAS FRONTIERE'
            NOTRI1 = 0
            GOTO 167
         ENDIF
C        ICI LE TRIANGLE NOTRI0 CONTIENT LE POINT N
C
C        ---------------------------------------------------------------
C        CONSTRUCTION DE L'ETOILE DU POINT N A PARTIR DES ARETES
C        DU TRIANGLE NOTRI0 ET RECURSIVITE
C        ---------------------------------------------------------------
C
C        RECHERCHE DES ARETES FRONTALIERES PARMI CELLES DE L'ETOILE
C        ET SUPPRESSION DE L'ETOILE DES TRIANGLES AU DELA DE LA FRONTIERE
C
 172     DO 180 I=1,3
C
C           RECHERCHE DES ARETES FRONTALIERES PARMI CELLES DE L'ETOILE
C           LE NUMERO NTOP DU TRIANGLE DE L'AUTRE COTE DE L'ARETE I
            NTOP = NOTRIA(I+3,NOTRI1)
            IF( NTOP .LE. 0 ) THEN
C              ARETE FRONTIERE SIMPLE
               GOTO 180
            ENDIF
C
C           LE NAOP NUMERO DE L'ARETE CONCERNEE DANS LE TRIANGLE OPPOSE
            IF( NOTRIA(4,NTOP) .EQ. NOTRI1 ) THEN
               NAOP = 1
            ELSE IF( NOTRIA(5,NTOP) .EQ. NOTRI1 ) THEN
               NAOP = 2
            ELSE
               NAOP = 3
            ENDIF
C
            NA1 = 0
            NA2 = N1AEOC
C
C           BOUCLE SUR LES ARETES DE L'ETOILE
C           =================================
 175        IF( NA2 .GT. 0 ) THEN
               IF( NAETOI(1,NA2) .NE. NTOP .OR.
     %             ABS(NAETOI(2,NA2)) .NE. NAOP ) THEN
C
C                 L'ARETE EST DIFFERENTE.PASSAGE A LA SUIVANTE
C                 --------------------------------------------
                  NA1 = NA2
                  NA2 = NAETOI(4,NA2)
                  GOTO 175
C
               ELSE
C
C                 L'ARETE NAOP DU TRIANGLE NTOP EST IDENTIQUE A NA2
C                 CETTE ARETE A DETRUIRE DE L'ETOILE EST ELLE FRONTALIERE ?
C                 ---------------------------------------------------------
C                 LES 2 SOMMETS DE L'ARETE
                  NS1 = NOTRIA(I,NOTRI1)
                  NS2 = I + 1
                  IF( NS2 .EQ. 4 ) NS2 = 1
                  NS2 = NOTRIA(NS2,NOTRI1)
                  CALL DVARFR( NS1, NS2 ,
     %                         NLSOFR, NBLFTR, NDARLF, NOSOAR, NONOUI )
                  IF( NONOUI .NE. 0 ) THEN
C                    ARETE RETROUVEE SUR LA FRONTIERE
C                    SUPPRESSION DU TRIANGLE NOTRI1 DE LA PILE DE L'ETOILE
                     GOTO 290
                  ENDIF
                  GOTO 180
               ENDIF
            ENDIF
C           L'ARETE NAOP N'EST PAS RETROUVEE DANS L'ETOILE
 180     CONTINUE
C
C        ICI AUCUNE DES 3 ARETES DU TRIANGLE NOTRI1 N'EST FRONTALIERE
C        FORMATION EFFECTIVE DE L'ETOILE DU POINT N AVEC LES TRIANGLES RESTANTS
C        PAR SUPPRESSION DES ARETES DOUBLES DE L'ETOILE ET AJOUT DES SIMPLES
C        ----------------------------------------------------------------------
         DO 280 I=1,3
C
C           SI    ( L'ARETE I DU TRIANGLE NOTRI1 N'APPARTIENT PAS
C                   AUX ARETES DE L'ETOILE NAETOI )
C           ALORS ELLE EST AJOUTEE A L'ETOILE DANS NAETOI
C           SINON SI ELLE N'EST PAS FRONTALIERE
C                 ALORS ELLE EST EMPILEE DANS NPILE POUR ETRE DETRUITE
C                       ELLE EST SUPPRIMEE DE L'ETOILE NAETOI
C
C           LE NUMERO NTOP DU TRIANGLE DE L'AUTRE COTE DE L'ARETE I
            NTOP = NOTRIA(I+3,NOTRI1)
            IF( NTOP .LE. 0 ) THEN
C              ARETE TRAITEE
               NC = -I
               GOTO 278
            ENDIF
            NC = I
C
C           LE NAOP NUMERO DE L'ARETE CONCERNEE DANS LE TRIANGLE OPPOSE
            IF( NOTRIA(4,NTOP) .EQ. NOTRI1 ) THEN
               NAOP = 1
            ELSE IF( NOTRIA(5,NTOP) .EQ. NOTRI1 ) THEN
               NAOP = 2
            ELSE
               NAOP = 3
            ENDIF
C
            NA1 = 0
            NA2 = N1AEOC
C
C           BOUCLE SUR LES ARETES DE L'ETOILE
C           =================================
 275        IF( NA2 .GT. 0 ) THEN
               IF( NAETOI(1,NA2) .NE. NTOP .OR.
     %             ABS(NAETOI(2,NA2)) .NE. NAOP ) THEN
C
C                 L'ARETE EST DIFFERENTE.PASSAGE A LA SUIVANTE
C                 --------------------------------------------
                  NA1 = NA2
                  NA2 = NAETOI(4,NA2)
                  GOTO 275
C
               ELSE
C
C                 L'ARETE NAOP DU TRIANGLE NTOP EST IDENTIQUE A NA2
C                 L'ARETE EST DOUBLE DANS L'ETOILE ET NON FRONTALIERE
C                 DESTRUCTION DE CETTE ARETE DE L'ETOILE
C                 ---------------------------------------------------
                  IF( NA1 .NE. 0 ) THEN
C                    L'ARETE NA2,PRECEDEE DE NA1,N'EST PAS LA
C                    PREMIERE DE L'ETOILE
                     NAETOI(4,NA1) = NAETOI(4,NA2)
C                    L'ARETE NA2 DEVIENT LA PREMIERE VIDE DE L'ETOILE
                     NAETOI(4,NA2) = N1AEVI
                     N1AEVI = NA2
                  ELSE
C                    L'ARETE NA2=N1AEOC EST LA 1-ERE DE L'ETOILE
                     NA1    = NAETOI(4,N1AEOC)
                     NAETOI(4,N1AEOC) = N1AEVI
                     N1AEVI = N1AEOC
                     N1AEOC = NA1
                  ENDIF
                  GOTO 280
               ENDIF
            ENDIF
C
C           L'ARETE NAOP N'EST PAS RETROUVEE DANS L'ETOILE
C           ELLE EST AJOUTEE A L'ETOILE AU DEBUT DES ARETES OCCUPEES
 278        IF( N1AEVI .LE. 0 ) THEN
C              SATURATION DES ARETES DE L'ETOILE
               NBLGRC(NRERR) = 2
               KERR(1) = 'SATURATION DES ARETES NAETOI'
               WRITE(KERR(MXLGER)(1:9),'(I9)') MXETOI
               KERR(2) = 'AUGMENTER MXETOI=' // KERR(MXLGER)(1:9)
               CALL LEREUR
               IERR = 10
               GOTO 9999
            ENDIF
            NA1           = N1AEVI
            N1AEVI        = NAETOI(4,N1AEVI)
            NAETOI(1,NA1) = NOTRI1
            NAETOI(2,NA1) = NC
            NAETOI(4,NA1) = N1AEOC
            N1AEOC        = NA1
 280     CONTINUE
C
C        LE TRIANGLE NOTRI1 EST MIS DANS LE TABLEAU NTETOI POUR ETRE DETRUIT
         IF( NBTRET .GE. MXETOI ) THEN
C           SATURATION DES TRIANGLES
            NBLGRC(NRERR) = 2
            KERR(1) = 'SATURATION DES TRIANGLES EMPILES'
            WRITE(KERR(MXLGER)(1:9),'(I9)') MXETOI
            KERR(2) = 'AUGMENTER MXETOI=' // KERR(MXLGER)(1:9)
            CALL LEREUR
            IERR = 9
            GOTO 9999
         ENDIF
         NBTRET = NBTRET + 1
         NTETOI( NBTRET ) = NOTRI1
C
C        PASSAGE AU TRIANGLE SUIVANT PAR LES ARETES SIMPLES DE L'ETOILE
c        --------------------------------------------------------------
 290     NA1 = N1AEOC
C
 300     IF( NA1 .GT. 0 ) THEN
C           LE NO LOCAL DE L'ARETE DE L'ETOILE
            NAOP = NAETOI(2,NA1)
            IF( NAOP .LT. 0 ) THEN
C              ARETE DEJA TRAITEE . PASSAGE A LA SUIVANTE
               NA1 = NAETOI(4,NA1)
               GOTO 300
            ELSE
C              L'ARETE ET SON TRIANGLE SONT ANALYSES ENSUITE
C              LE TRIANGLE DE L'AUTRE COTE
               NOTRI1 = NOTRIA( NAOP+3, NAETOI(1,NA1) )
               IF( NOTRI1 .LE. 0 ) THEN
C                 L'ARETE EST DONC TRAITEE
                  NAETOI(2,NA1) = -NAOP
C                 L'ARETE SUIVANTE
                  NA1 = NAETOI(4,NA1)
                  GOTO 300
               ENDIF
C              LE TEMOIN DE RECHERCHE EFFECTUEE
               NAETOI(2,NA1) = -NAOP
C
C              *****************************************************************
C              LE CERCLE CIRCONSCRIT AU TRIANGLE NOTRI1 CONTIENT IL LE POINT N
C              *****************************************************************
               CALL TRPTDC( UNPEPS, N, PXYD, NOTRI1, CETRIA, NONOUI )
               IF( NONOUI .EQ. 0 ) THEN
C                 NON : PASSAGE A L'ARETE SUIVANTE DE L'ETOILE
                  GOTO 300
               ELSE
C                 OUI : AJOUT DU TRIANGLE
                  GOTO 172
               ENDIF
            ENDIF
         ENDIF
C
C        ICI L'ETOILE DU POINT N EST COMPLETE
C        ==============================================================
C        LE NO DE TRIANGLE            => NO TRIANGLE AU DELA DE L'ARETE
C        LE NO LOCAL DANS LE TRIANGLE => NO 1-ER SOMMET DE L'ARETE
C        LE NO INUTILISE              => NO 2-ME SOMMET DE L'ARETE
C        ==============================================================
         NA1 = N1AEOC
C        BOUCLE SUR LES ARETES
 630     IF( NA1 .GT. 0 ) THEN
C           LE NO DU TRIANGLE ET LOCAL DE L'ARETE
            NOTRI1 = NAETOI(1,NA1)
            I      = ABS( NAETOI(2,NA1) )
C           LE NUMERO DU TRIANGLE AU DELA DE L'ARETE
            NT     = NOTRIA(I+3,NOTRI1)
C
C           LE POINT N EST IL SUR UNE ARETE DE L'ENVELOPPE CONVEXE?
            IF( NOTRI1 .EQ. NOTRI0     .AND. NT .EQ. 0  .AND.
     %          ABS(CMIN) .LT. 0.001D0 .AND. I  .EQ. IMIN ) THEN
C              OUI : TEMOIN D'ARETE SPECIALE SUR L'ENVELOPPE CONVEXE
               NAETOI(1,NA1) = -I
            ELSE
               NAETOI(1,NA1) = NT
            ENDIF
C
C           LE NUMERO DES 2 SOMMETS DANS LE SENS DIRECT
            NS2    = I + 1
            IF( NS2 .EQ. 4 ) NS2 = 1
C           NUMERO DU SOMMET 1 DE L'ARETE
            NAETOI(2,NA1) = NOTRIA(I  ,NOTRI1)
C           NUMERO DU SOMMET 2 DE L'ARETE
            NAETOI(3,NA1) = NOTRIA(NS2,NOTRI1)
CCCC
CCCC           QUALITE DU TRIANGLE A CREER
CCC            CALL QUTR2D( PXYD(1,NAETOI(2,NA1)),
CCC     %                   PXYD(1,NAETOI(3,NA1)),
CCC     %                   PXYD(1,N), Q )
CCC            IF( Q .LT. QUADEG ) THEN
CCCC              TRIANGLE DEGENERE : LE POINT N EST SUPPRIME
CCC               WRITE(IMPRIM,*)'TRIANGLE ', NAETOI(2,NA1), NAETOI(3,NA1),
CCC     %                         N,' DE QUALITE ',Q
CCCCCC               GOTO 100
CCC            ENDIF
C
C           PASSAGE A L'ARETE SUIVANTE
            NA1 = NAETOI(4,NA1)
            GOTO 630
         ENDIF
C
C        ==============================================================
C        DESTRUCTION DES TRIANGLES TABLEAU NTETOI CONTENANT LE SOMMET N
C                         ET ETOILANT CE POINT
C        ==============================================================
         DO 650 I=1,NBTRET
            NT = NTETOI( I )
C           TRACE DU TRIANGLE A DETRUIRE
            CALL DVTRTR( PXYD, NOTRIA, NT, NCBLAN, NCNOIR )
C           NT DEVIENT LE PREMIER TRIANGLE VIDE
            NOTRIA(1,NT) = 0
            NOTRIA(4,NT) = N1TRVI
            N1TRVI       = NT
 650     CONTINUE
C
C        ====================================================================
C        DECLARATION DES TRIANGLES ETOILANT LE POINT N
C        ====================================================================
         NBTRET = 0
         NA2    = N1AEOC
C
C        BOUCLE SUR LES ARETES PERIPHERIQUES DE L'ETOILE
 750     IF( NA2 .GT. 0 ) THEN
            IF( NAETOI(1,NA2) .LT. 0 ) THEN
C              OUI : TEMOIN D'ARETE SPECIALE SUR L'ENVELOPPE CONVEXE
C              LE POINT N EST SUR UNE ARETE DE L'ENVELOPPE CONVEXE
C              LE TRIANGLE (POINT N -> ARETE) NE DOIT PAS ETRE FORME
               GOTO 780
            ENDIF
C
C           L'ARETE NA2 ET LE POINT N FORMENT UN NOUVEAU TRIANGLE A DECLARER
            IF( N1TRVI .LE. 0 ) THEN
C              SATURATION DES TRIANGLES
               NBLGRC(NRERR) = 1
               KERR(1) = 'SATURATION DES TRIANGLES'
               CALL LEREUR
               IERR = 3
               GOTO 9999
            ENDIF
C
C           REMPLISSAGE DU TRIANGLE NT DE SOMMETS CEUX DE L'ARETE NA2 ET N
            NS1 = NAETOI(2,NA2)
            NS2 = NAETOI(3,NA2)
C           LE TRIANGLE EXTERIEUR A L'ARETE AU DELA DE L'ARETE
            NT1 = NAETOI(1,NA2)
C
C           MISE A JOUR DU 1-ER TRIANGLE VIDE
            NT     = N1TRVI
            N1TRVI = NOTRIA(4,N1TRVI)
C
C           LES COORDONNEES DU CENTRE DE LA BOULE CIRCONSCRITE
C           AU TRIANGLE NT ET CARRE DE SON RAYON
            IER = -1
            CALL CENCED( PXYD(1,NS1), PXYD(1,NS2) ,
     %                   PXYD(1,N),   CETRIA(1,NT), IER )
            IF( IER .NE. 0 ) THEN
C              LE TRIANGLE EST DEGENERE => IL EST SUPPRIME
               WRITE(IMPRIM,*) 'TRIANGLE ',NT,' : ST ',
     %                          NS1,' ',NS2,' ',N, ' PLAT'
               GOTO 780
            ENDIF
C
C           NS1 NS2 SOMMETS DE L'ARETE DANS LE SENS DIRECT DU TRIANGLE NT
            NOTRIA(1,NT) = NS1
            NOTRIA(2,NT) = NS2
            NOTRIA(3,NT) = N
C
C           LE CHAINAGE DES TRIANGLES PAR LES ARETES
            NOTRIA(4,NT) = NT1
            NOTRIA(5,NT) = 0
            NOTRIA(6,NT) = 0
C
C           RECHERCHE DU NO LOCAL A NT1 DE L'ARETE NS1 NS2
            IF( NT1 .GT. 0 ) THEN
               DO 760 J=1,3
                  IF( NOTRIA(J,NT1) .EQ. NS2 ) GOTO 770
 760           CONTINUE
 770           NOTRIA(J+3,NT1) = NT
            ENDIF
C
C           TRACE DU TRIANGLE CONSTRUIT AVEC LE POINT N
            NCOUL = NCOUL + 1
            IF( NCOUL .GT. NDCOUL ) NCOUL = NDCORE + 1
            CALL DVTRTR( PXYD, NOTRIA, NT, NCOUL, NCBLAN )
C
C           LA PILE DES TRIANGLES CREES DANS L'ETOILE
            IF( NBTRET .GE. MXETOI ) THEN
C              SATURATION DES TRIANGLES DE L'ETOILE
               NBLGRC(NRERR) = 2
               KERR(1) = 'SATURATION DES TRIANGLES EMPILES NTETOI'
               WRITE(KERR(MXLGER)(1:9),'(I9)') MXETOI
               KERR(2) = 'AUGMENTER MXETOI=' // KERR(MXLGER)(1:9)
               CALL LEREUR
               IERR = 9
               GOTO 9999
            ENDIF
            NBTRET = NBTRET + 1
            NTETOI(NBTRET) = NT
C
C           MISE A JOUR D'UN TRIANGLE CONTENANT CHAQUE SOMMET
C           TOUS LES TRIANGLES ONT PU DISPARAITRE
C           CAS DE L'ETOILE REDUITE AUX 2 TRIANGLES INITIAUX
            NOTRSO(NS1) = NT
            NOTRSO(NS2) = NT
C
C           PASSAGE A L'ARETE SUIVANTE DE L'ETOILE
 780        NA2 = NAETOI(4,NA2)
            GOTO 750
         ENDIF
C
C        LE POINTEUR SOMMET => TRIANGLE
         NOTRSO(N) = NT
C
C        COMPLETION DES CHAINAGES DES TRIANGLES CREES DANS L'ETOILE
C        ----------------------------------------------------------
 800     IF( NBTRET .GT. 0 ) THEN
C            LE HAUT DE LA PILE
             NT = NTETOI(NBTRET)
C            LE TRIANGLE EST DEPILE
             NBTRET = NBTRET - 1
             IF( NT .LE. 0 ) GOTO 800
C
C            QUEL EST LE TRIANGLE OPPOSE A L'ARETE 2 DE SOMMETS NS2 N
             NS2  = NOTRIA(2,NT)
             DO 810 J=1,NBTRET
                NT1 = NTETOI(J)
                IF( NT1 .LE. 0 ) GOTO 810
                IF( NOTRIA(1,NT1) .NE. NS2 ) GOTO 810
C               L'ARETE EST RETROUVEE
                NOTRIA(5,NT ) = NT1
                NOTRIA(6,NT1) = NT
C               SI LE TRIANGLE NT1 EST TOTALEMENT CHAINE
C               IL EST RETIRE DE L'ETOILE
                IF( NOTRIA(5,NT1) .GT. 0 ) NTETOI(J)=0
                GOTO 820
 810         CONTINUE
C
C            SI LE TRIANGLE NT EST TOTALEMENT CHAINE SAUT DE L'ARETE 3
 820         IF( NOTRIA(6,NT) .GT. 0 ) GOTO 800
C
C            QUEL EST LE TRIANGLE OPPOSE A L'ARETE 3 DE SOMMETS N NS1
             NS1  = NOTRIA(1,NT)
             DO 830 J=1,NBTRET
                NT1 = NTETOI(J)
                IF( NT1 .EQ. 0 ) GOTO 830
                IF( NOTRIA(2,NT1) .NE. NS1 ) GOTO 830
C               L'ARETE EST RETROUVEE
                NOTRIA(6,NT ) = NT1
                NOTRIA(5,NT1) = NT
C               SI LE TRIANGLE NT1 EST TOTALEMENT CHAINE
C               IL EST RETIRE DE L'ETOILE
                IF( NOTRIA(6,NT1) .GT. 0 ) NTETOI(J)=0
                GOTO 800
 830         CONTINUE
CCC             WRITE(IMPRIM,*) 'DVAJPC: ANOMALIE EN 830'
CCC             WRITE(IMPRIM,*) 'POINT NP=',NP,' N=',N
CCCCCC             CALL XVPAUSE
C            RETOUR EN HAUT DE PILE
             GOTO 800
         ENDIF
C
C        L'ETOILE EST TRAITEE . PASSAGE AU POINT N SUIVANT
         GOTO 100
      ENDIF
C===============================================================================
C     FIN DU JEU DE POINTS NOSOM1 A NBSOMT
C===============================================================================
 9999 RETURN
      END
