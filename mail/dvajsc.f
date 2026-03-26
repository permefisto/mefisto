      SUBROUTINE DVAJSC( ZEMEPS, UNPEPS,
     %                   MXSOMM, NOSOM1, NOSOM2, MXTRIA, PXYD,
     %                   N1TRVI, NOTRIA, CETRIA, NOTRSO,
     %                   NBLFTR, NDARLF, MXSOAR, NOSOAR,
     %                   NLSOFR,
     %                   MXETRI, NAETOI, NTETOI,
     %                   NBPTSU, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    GENERER LES TRIANGLES DELAUNAY A PARTIR DES POINTS
C -----    NOSOM1 A NOSOM2 DU TABLEAU PXYD
C          AVEC LA CONTRAINTE DE NE PAS PERDRE D'ARETE FRONTALIERE
C          SI UN POINT PERD UNE ARETE CE POINT EST SUPPRIME
C
C ENTREES:
C --------
C ZEMEPS :  -EPSILON POUR LA COORDONNEE BARYCENTRIQUE
C UNPEPS : 1+EPSILON POUR LA COORDONNEE BARYCENTRIQUE
C MXSOMM : NOMBRE MAXIMAL DE SOMMETS PERMIS POUR LA TRIANGULATION
C NOSOM1 : NUMERO DU PREMIER SOMMET A AJOUTER A LA TRIANGULATION
C MXTRIA : NOMBRE MAXIMAL DE TRIANGLES DECLARABLES
C PXYD   : TABLEAU DES COORDONNEES 2D DES POINTS
C NBLFTR : NOMBRE DE LIGNES FERMEES LIMITANT LA SURFACE A TRIANGULER
C NDARLF : NUMERO DE LA PREMIERE ARETE DE CHAQUE LIGNE FERMEE DANS
C          LE TABLEAU NOSOAR
C NBARFR : NOMBRE D'ARETES  DE LA FRONTIERE OU CONTOUR DU DOMAINE =
C          NOMBRE DE POINTS DE LA FRONTIERE OU CONTOUR DU DOMAINE
C MXSOAR : NOMBRE MAXIMAL D'ARETES FRONTALIERES DECLARABLES
C NOSOAR : NUMERO DES 2 SOMMETS DE CHAQUE ARETE ET
C          POINTEUR SUR L'ARETE SUIVANTE
C NLSOFR : TABLEAU DU NUMERO DE LIGNE(1 A NBLFTR) DES SOMMETS
C         -NUMERO DE POINT INTERNE UTILISATEUR IMPOSE, 0 SINON
C
C ENTREES ET SORTIES :
C --------------------
C NOSOM2 : EN ENTREE NUMERO DU DERNIER SOMMET A AJOUTER
C                    A LA TRIANGULATION
C          EN SORTIE NUMERO DU DERNIER SOMMET EFFECTIVEMENT AJOUTE
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
C NBPTSU : NOMBRE DE POINTS SUPPRIMES, NON TRIANGULES
C IERR   : 0 SI PAS D'ERREUR
C          1 SATURATION DES SOMMETS
C          3 SATURATION DES TRIANGLES
C          8 AUCUN CERCLE CIRCONSCRIT AUX TRIANGLES ACTUELS NE CONTIENT
C            LE POINT COURANT
C
C TABLEAUX AUXILIAIRES :
C ----------------------
C MXETRI : NOMBRE MAXIMAL D'ARETES OU TRIANGLES DECLARABLES
C          DANS UNE ETOILE
C NAETOI : ENTIER (1:4,1:MXETRI)
C NTETOI : ENTIER (1:MXETRI)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC    JUILLET 1994
C....................................................................012
      PARAMETER        (QUADEG=0.0001)
      include"./incl/gsmenu.inc"
      include"./incl/trvari.inc"
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      DOUBLE PRECISION  ZEMEPS, UNPEPS
      INTEGER           NOTRIA(6,MXTRIA),
     %                  NOTRSO(MXSOMM),
     %                  NAETOI(4,MXETRI),
     %                  NLSOFR(MXSOMM),
     %                  NDARLF(NBLFTR),
     %                  NOSOAR(3,MXSOAR),
     %                  NTETOI(MXETRI)
      DOUBLE PRECISION  PXYD(3,MXSOMM),
     %                  CETRIA(3,MXTRIA)
      DOUBLE PRECISION  XN, YN, D1, D2, D3, CB(3), CB1, CB2, CB3
      EQUIVALENCE      (CB(1),CB1),(CB(2),CB2),(CB(3),CB3)
C
      NBPTSU = 0
C
C     REINITIALISATION A VIDE DES ARETES DE L'ETOILE
      N1AEVI = 1
      N1AEOC = 0
      DO 10 I=1,MXETRI
C        NUMERO DU TRIANGLE DE L'AUTRE COTE DE L'ARETE
C         NAETOI(1,I) = 0
C        NUMERO LOCAL AU TRIANGLE DE L'ARETE
C         NAETOI(2,I) = 0
C        NUMERO DANS NAETOI DE L'ARETE SUIVANTE
         NAETOI(4,I) = I+1
  10  CONTINUE
      NAETOI(4,MXETRI) = 0
C
      N      = NOSOM1 - 1
      NBSOMM = N
C
C     AJOUT EVENTUEL DES SOMMETS NOSOM1 A NOSOM2
C     ==========================================
 100  N = N + 1
      IF( N .LE. NOSOM2 ) THEN
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
         NBARFT = 0
C
C        RECHERCHE DU 1-ER TRIANGLE NOTRI1
C        CONTENANT LE POINT N DE COORDONNEES ( XN, YN )
C        =======================================================
         XN = PXYD(1,N)
         YN = PXYD(2,N)
C
C        NABOUL = 0 SI PAS DE TRIANGLE RETROUVE CONTENANT LE POINT N
         NABOUL = 0
C
C        LE 1-ER TRIANGLE A CONSULTER  ( DERNIER CONSTRUIT )
         I = NBSOMM
C
 160     I = I - 1
         IF( NOTRSO(I) .LE. 0 ) GOTO 160
         NOTRI1 = NOTRSO(I)
C
 165     IF( NABOUL .EQ. 0 ) THEN
C           ***********************************************
C           LE TRIANGLE NOTRI1 CONTIENT IL LE POINT XN YN ?
C           ***********************************************
            CALL TRPTDT( N, PXYD, NOTRIA, ZEMEPS, UNPEPS,
     %                   NOTRI1, CB1, CB2, CB3, I )
            IF( NOTRI1 .LE. 0 ) THEN
C              AUCUN TRIANGLE NE CONTIENT CE POINT
               IF( NLSOFR(N) .EQ. 0 ) THEN
C                 CE POINT N EST SUPPRIME CAR SUPPRIMABLE
CCC               WRITE(IMPRIM,*)'DVAJSC 1:',N,' SUPPRIME'
                  NBPTSU = NBPTSU + 1
                  GOTO 100
               ELSE
C                 CE POINT N NE PEUT ETRE SUPPRIME CAR IMPOSE
                  IF( NOTRI1 .GE. 0 ) THEN
                     WRITE(IMPRIM,*) 'DVAJSC: PAS DE TRIANGLE OPPOSE',
     %                                NOTRI1
                     CALL XVPAUSE
                  ENDIF
C
C                 JONCTION DU POINT AVEC L'ARETE FRONTIERE DU
C                 TRIANGLE -(NOTRI1 <0 EN RETOUR)
                  NOTRI1 = -NOTRI1
                  WRITE(IMPRIM,*) 'DVAJSC: TRIANGLE DE SORTIE',
     %                             NOTRI1
                  CALL XVPAUSE
               ENDIF
            ENDIF
C
C           LE TRIANGLE NOTRI1 CONTIENT XN YN
            DO 170 I=1,3
               IF( CB(I) .GT. 0.8D0 ) THEN
C                 POINT SUSCEPTIBLE D'ETRE TROP PROCHE DE L'UN DES SOMMETS
C                 RECHERCHE DES LONGUEURS DES COTES POUR RELATIVISER
                  NS1 = NOTRIA(I,NOTRI1)
C                 DISTANCE DU POINT AU SOMMET I
                  D1  = SQRT( (XN-PXYD(1,NS1)) ** 2 +
     %                        (YN-PXYD(2,NS1)) ** 2 )
                  IF( D1 .LT. 0.4D0 * PXYD(3,NS1) ) THEN
C                    POINT SUPPRIME  CAR TROP PRES DU SOMMET NS1
CCC                     WRITE(IMPRIM,*)'DVAJSC 2:',N,' SUPPRIME'
                     NBPTSU = NBPTSU + 1
                     GOTO 100
                  ENDIF
               ENDIF
               IF( CB(I) .LT. 0.2D0 ) THEN
C                 LE POINT EST PROCHE DU COTE I+1
C                 CE COTE EST-IL SUR LA FRONTIERE ?
                  IF( I .NE. 3 ) THEN
                     I2 = I + 1
                  ELSE
                     I2 = 1
                  ENDIF
                  IF( I .NE. 1 ) THEN
                     I1 = I - 1
                  ELSE
                     I1 = 3
                  ENDIF
                  NS1 = NOTRIA(I2,NOTRI1)
                  NS2 = NOTRIA(I1,NOTRI1)
                  CALL DVARFR( NS1,NS2,NLSOFR,NBLFTR,NDARLF,NOSOAR,
     %                         NONOUI )
                  IF( NONOUI .NE. 0 ) THEN
C                    DISTANCE DU POINT AU SOMMET I2
                     D1 = (XN-PXYD(1,NS1)) ** 2 + (YN-PXYD(2,NS1)) ** 2
C                    DISTANCE DU POINT AU SOMMET I1
                     D2 = (XN-PXYD(1,NS2)) ** 2 + (YN-PXYD(2,NS2)) ** 2
C                    DISTANCE ENTRE LES 2 POINTS DE L'ARETE FRONTIERE
                     D3 = (PXYD(1,NS1)-PXYD(1,NS2)) ** 2 +
     %                    (PXYD(2,NS1)-PXYD(2,NS2)) ** 2
                     IF( D3 + D3 .GT. 3 * (D1 + D2) ) THEN
C                       POINT SUPPRIME  CAR TROP PRES DE L'ARETE
C                       L'ANGLE EN XN,YN SERAIT SUPERIEUR A 120 DEGRES
CCC                        WRITE(IMPRIM,*)'DVAJSC 3:',N,' SUPPRIME'
                        NBPTSU = NBPTSU + 1
                        GOTO 100
                     ENDIF
                  ENDIF
               ENDIF
 170        CONTINUE
C           CE TRIANGLE CONTIENT LE POINT N
            NABOUL = 1
         ELSE
C           ********************************************************************
C           LE CERCLE CIRCONSCRIT AU TRIANGLE NOTRI1 CONTIENT IL LE POINT XN YN
C           ********************************************************************
            XG = REAL( CETRIA(1,NOTRI1) - XN )
            YG = REAL( CETRIA(2,NOTRI1) - YN )
            IF( XG*XG + YG*YG .GT. CETRIA(3,NOTRI1)*UNPEPS ) GOTO 290
         ENDIF
C
C        ICI LE CERCLE CIRCONSCRIT DU TRIANGLE NOTRI1 CONTIENT XN YN
C        RECHERCHE DES ARETES DE L'ETOILE PARMI LES 3 DU TRIANGLE NOTRI1
C        ET NON EMPLOI DU TRIANGLE SI UNE ARETE FRONTIERE EST DOUBLE
C        ---------------------------------------------------------------
         DO 180 I=1,3
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
                  CALL DVARFR( NS1, NS2,
     %                         NLSOFR,NBLFTR,NDARLF,NOSOAR,NONOUI )
                  IF( NONOUI .NE. 0 ) THEN
C                    ARETE RETROUVEE SUR LA FRONTIERE
C                    SUPPRESSION DU TRIANGLE NOTRI1 DE LA PILE DE L'ETOILE
                     NBARFT = NBARFT + 1
                     GOTO 290
                  ENDIF
                  GOTO 180
               ENDIF
            ENDIF
C           L'ARETE NAOP N'EST PAS RETROUVEE DANS L'ETOILE
 180     CONTINUE
C
C        ICI LE CERCLE CIRCONSCRIT DU TRIANGLE NOTRI1 CONTIENT XN YN
C        ET AUCUNE DE SES 3 ARETES N'EST FRONTALIERE
C        SUPPRESSION DES ARETES DOUBLES DE L'ETOILE ET AJOUT DES SIMPLES
C        ---------------------------------------------------------------
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
 278        NA1           = N1AEVI
            N1AEVI        = NAETOI(4,N1AEVI)
            NAETOI(1,NA1) = NOTRI1
            NAETOI(2,NA1) = NC
            NAETOI(4,NA1) = N1AEOC
            N1AEOC        = NA1
 280     CONTINUE
C
C        LE TRIANGLE NOTRI1 EST MIS DANS LE TABLEAU NTETOI POUR ETRE DETRUIT
         NBTRET = NBTRET + 1
         NTETOI( NBTRET ) = NOTRI1
C
C        PASSAGE AU TRIANGLE SUIVANT PAR LES ARETES SIMPLES DE L'ETOILE
c        --------------------------------------------------------------
 290     NA1 = N1AEOC
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
               GOTO 165
            ENDIF
         ENDIF
C
C        L'ETOILE DU POINT N EST COMPLETE
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
            NAETOI(1,NA1) = NT
C           LE NUMERO DES 2 SOMMETS DANS LE SENS DIRECT
            NS2    = I + 1
            IF( NS2 .EQ. 4 ) NS2 = 1
C           NUMERO DU SOMMET 1 DE L'ARETE
            NAETOI(2,NA1) = NOTRIA(I  ,NOTRI1)
C           NUMERO DU SOMMET 2 DE L'ARETE
            NAETOI(3,NA1) = NOTRIA(NS2,NOTRI1)
C
C           QUALITE DU TRIANGLE A CREER
            CALL QUTR2D( PXYD(1,NAETOI(2,NA1)),
     %                   PXYD(1,NAETOI(3,NA1)),
     %                   PXYD(1,N), Q )
            IF( Q .LT. QUADEG ) THEN
C              TRIANGLE DEGENERE : LE POINT N EST SUPPRIME
CCC               WRITE(IMPRIM,*)'DVAJSC 4:',N,' SUPPRIME'
               NBPTSU = NBPTSU + 1
               GOTO 100
            ENDIF
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
C           NT DEVIENT LE PREMIER TRIANGLE VIDE
            NOTRIA(1,NT) = 0
            NOTRIA(4,NT) = N1TRVI
            N1TRVI       = NT
 650     CONTINUE
C
C        LE POINT N EST ACCEPTE SOUS LE NUMERO NBSOMM(<=NOSOM2) DE PXYD
C        ==============================================================
C        PAS DE RISQUE DE DEBORDEMENT DE PXYD
         NBSOMM = NBSOMM + 1
         PXYD(1,NBSOMM) = XN
         PXYD(2,NBSOMM) = YN
         PXYD(3,NBSOMM) = PXYD(3,N)
C
C        ====================================================================
C        DECLARATION DES TRIANGLES ETOILANT LE POINT NBSOMM
C        ====================================================================
         NT     = 0
         NBTRET = 0
         NA2    = N1AEOC
C
C        BOUCLE SUR LES ARETES PERIPHERIQUES DE L'ETOILE
 750     IF( NA2 .GT. 0 ) THEN
C           L'ARETE NA2 ET LE POINT N FORMENT UN NOUVEAU TRIANGLE A DECLARER
            IF( N1TRVI .LE. 0 ) THEN
C              SATURATION DES TRIANGLES
               NBLGRC(NRERR) = 1
               KERR(1) = 'SATURATION DES TRIANGLES'
               CALL LEREUR
               IERR = 3
               GOTO 9999
            ENDIF
C           MISE A JOUR DU 1-ER TRIANGLE VIDE
            NT     = N1TRVI
            N1TRVI = NOTRIA(4,N1TRVI)
C
C           REMPLISSAGE DU TRIANGLE NT DE SOMMETS CEUX DE L'ARETE NA2 ET NBSOMM
            NS1 = NAETOI(2,NA2)
            NS2 = NAETOI(3,NA2)
C           NS1 NS2 SOMMETS DE L'ARETE DANS LE SENS DIRECT DU TRIANGLE NT
C           NBSOMM EST TOUJOURS SUPERIEUR A NS1 NS2
            NOTRIA(1,NT) = NS1
            NOTRIA(2,NT) = NS2
            NOTRIA(3,NT) = NBSOMM
C           LE CHAINAGE DES TRIANGLES PAR LES ARETES
C           LE TRIANGLE EXTERIEUR A L'ARETE AU DELA DE L'ARETE
            NT1 = NAETOI(1,NA2)
            NOTRIA(4,NT) = NT1
            NOTRIA(5,NT) = 0
            NOTRIA(6,NT) = 0
C           RECHERCHE DU NO LOCAL A NT1 DE L'ARETE NS1 NS2
            IF( NT1 .GT. 0 ) THEN
               DO 760 J=1,3
                  IF( NOTRIA(J,NT1) .EQ. NS2 ) GOTO 770
 760           CONTINUE
 770           NOTRIA(J+3,NT1) = NT
            ENDIF
C
C           MISE A JOUR D'UN TRIANGLE CONTENANT CHAQUE SOMMET
C           TOUS LES TRIANGLES ONT PU DISPARAITRE
C           CAS DE L'ETOILE REDUITE AUX 2 TRIANGLES INITIAUX
            NOTRSO(NS1) = NT
            NOTRSO(NS2) = NT
C
C           LES COORDONNEES DU CENTRE DE LA BOULE CIRCONSCRITE
C           AU TRIANGLE NT ET CARRE DE SON RAYON
            IER = 1
            CALL CENCED( PXYD(1,NS1), PXYD(1,NS2), PXYD(1,NBSOMM),
     %                   CETRIA(1,NT), IER )
CCCC
CCCC           TRACE DU TRIANGLE
CCC            IF( NBARFT .GT. 0 ) THEN
CCC               CALL DVTRTR( PXYD, NOTRIA, NT, NCBLAN, NCNOIR )
CCC            ENDIF
C
C           LA PILE DES TRIANGLES CREES DANS L'ETOILE
            NBTRET = NBTRET + 1
            IF( NBTRET .GT. MXETRI ) THEN
C              SATURATION DES TRIANGLES DE L'ETOILE
               NBLGRC(NRERR) = 1
               KERR(1) = 'SATURATION DES TRIANGLES DE L''ETOILE'
               CALL LEREUR
               IERR = 3
               GOTO 9999
            ENDIF
            NTETOI(NBTRET) = NT
C
C           PASSAGE A L'ARETE SUIVANTE DE L'ETOILE
            NA2 = NAETOI(4,NA2)
            GOTO 750
         ENDIF
C
C        LE POINTEUR SOMMET => TRIANGLE
         NOTRSO(NBSOMM) = NT
C
C        COMPLETION DES CHAINAGES DES TRIANGLES CREES DANS L'ETOILE
 800     IF( NBTRET .GT. 0 ) THEN
C            LE HAUT DE LA PILE
             NT = NTETOI(NBTRET)
C            LE TRIANGLE EST DEPILE
             NBTRET = NBTRET - 1
             IF( NT .LE. 0 ) GOTO 800
C
C            QUEL EST LE TRIANGLE OPPOSE A L'ARETE 2 DE SOMMETS NS1 NS2
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
C            QUEL EST LE TRIANGLE OPPOSE A L'ARETE 3 DE SOMMETS NS1 NS2
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
             WRITE(IMPRIM,*) 'ANOMALIE EN 830'
C            RETOUR EN HAUT DE PILE
             GOTO 800
         ENDIF
C
C        L'ETOILE EST TRAITEE . PASSAGE AU POINT N SUIVANT
         GOTO 100
      ENDIF
C
C     MISE A JOUR DE NOSOM2
      NOSOM2 = NBSOMM
C===============================================================================
C     FIN DU JEU DE POINTS NOSOM1 A NOSOM2
C===============================================================================
 9999 RETURN
      END
