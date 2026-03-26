      SUBROUTINE ADAP21( NBISO,  VALISO, DFMXIS,
     %                   NCAS  , NDSM,   NTDL,   TEMPER,
     %                   NUTYEL, NBELEM, MNNPEF,
     %                   NBNOEL, NUNDEL, NBPOEL, NUPTEL,
     %                   NBPOIT, XYZPOI,
     %                   MXYZST, NBXYZ,  XYZST,  NOXYZ,
     %                   MXNVTR, NBNVTR, NSNVTR, NUPLIS, NUISOP,
     %                   NUMILF, NUMXLF, L1LGFR,
     %                   MXARCH, NBARCH, N1ARCH, NSARCH,
     %                   NBAEFU, MXSOMM, SOLEL,  COPOE ,
     %                   MXPILE, LTPILE, LAPILE,
     %                   XYVAL,  FBASE,  IERR  )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    ADAPTER LE MAILLAGE 2D D'UN TYPE D'EF EN FONCTION DU TRACE
C -----    DES ISOVALEURS
C
C ENTREES :
C ---------
C NUITER : NUMERO DE L'ITERATION D'ADAPTATION
C KNOMOB : NOM DE L'OBJET (SURFACE) DE MAILLAGE A ADAPTER
C
C NBISO  : NOMBRE D'ISOVALEURS A TRACER
C VALISO : TABLEAU DES VALEURS DES ISO A TRACER
C DFMXIS : DIFFERENCE MAXIMALE ENTRE 2 ISOVALEURS
C
C NCAS   : NUMERO DU CAS A TRAITER
C NDSM   : NOMBRE TOTAL DE VECTEURS TEMPERATURES
C NTDL   : NOMBRE TOTAL DE DEGRES DE LIBERTE EN TEMPERATURE
C TEMPER : LES NTDL * NDSM TEMPERATURES
C
C NUTYEL : NUMERO DU TYPE D'EF A TRAITER ICI
C NBELEM : NOMBRE D'ELEMENTS FINIS DE CE TYPE
C MNNPEF : ADRESSE MCN DU TABLEAU NPEF"... DE CE TYPE D'EF
C NBNOEL : NOMBRE DE  NOEUDS DE L'ELEMENT  FINI  DE CE TYPE
C NUNDEL : NUMERO DES NOEUDS DES  ELEMENTS FINIS DE CE TYPE
C NBPOEL : NOMBRE DE  POINTS DE L'ELEMENT  FINI  DE CE TYPE
C NUPTEL : NUMERO DES POINTS DES  ELEMENTS FINIS DE CE TYPE
C
C NBPOIT : NOMBRE TOTAL DE POINTS DU MAILLAGE
C XYZPOI : COORDONNEES DES POINTS DU MAILLAGE
C
C MXYZST : NOMBRE MAXIMAL DE SOMMETS DECLARABLES DANS XYZST
C MXNVTR : NOMBRE MAXIMAL DE TRIANGLES DECLARABLES DANS NSNVTR
C MXARCH : NOMBRE MAXIMAL D'ARETES DECLARABLES DANS NSARCH
C NBAEFU : NOMBRE DE SUBDIVISIONS EN SOUS-ARETES D'UN COTE DE L'EF UNITE
C MXSOMM : NOMBRE MAXIMAL DE SOMMETS DES SOUS-TRIANGLES DECLARABLES
C MXPILE : NOMBRE MAXIMUM DE SOUS-TRIANGLES DANS LA PILE LTPILE ET LAPILE
C
C MODIFIES :
C ----------
C NBXYZ  : NOMBRE DE SOMMETS ACTUELS DE LA TRIANGULATION
C XYZST  : 3 COORDONNEES DES NBXYZ SOMMETS DE LA TRIANGULATION
C NOXYZ  : NUMERO DANS XYZST DES SOMMETS DES SOUS-TRIANGLES
C          DE L'EF ACTUEL (ILS SONT IDENTIFIES)
C
C NBNVTR : NOMBRE ACTUEL DE TRIANGLES DANS NSNVTR
C NSNVTR : NUMERO DES 3 SOMMETS NO SURFACE DES NBNVTR TRIANGLES ACTUELS
C
C NUMILF : NUMERO MINIMAL DES LIGNES FRONTIERES DE L'OBJET
C NUMXLF : NUMERO MAXIMAL DES LIGNES FRONTIERES DE L'OBJET
C L1LGFR : POUR CHAQUE LIGNE DE NUMILF A NUMXLF POINTE SUR LA PREMIERE
C          ARETE CHAINEE DE LA LIGNE
C
C NBARCH : NOMBRE D'ARETES CHAINEES DANS NSARCH
C N1ARCH : NUMERO DANS NSARCH DE LA PREMIERE ARETE VIDE
C NSARCH : NUMERO DES 2 SOMMETS DE L'ARETE, ARETE PRECEDENTE ET SUIVANTE
C
C FBASE  : VALEUR DES POLYNOMES DE BASE EN UN POINT
C NBXYZ  : NOMBRE DE SOMMETS ACTUELS DE LA TRIANGULATION
C XYZST  : 3 COORDONNEES DES NBXYZ SOMMETS DE LA TRIANGULATION
C NBNVTR : NOMBRE ACTUEL DE TRIANGLES DANS NSNVTR
C NSNVTR : NUMERO DES 3 SOMMETS NO SURFACE DES NBNVTR TRIANGLES ACTUELS
C NUPLIS : NUMERO DE POINT OU LIGNE DE CHACUN DES SOMMETS
C          -2 000 000 - NUPLIS() ANCIEN SI LE SOMMET A DEJA ETE SUPPRIME
C          -NP SI NP EST LE NUMERO DU POINT UTILISATEUR DE CE SOMMET
C          -1 000 000 - NL1 - 1000 * NL2 SI LE SOMMET APPARTIENT
C                     AUX 2 LIGNES NL1 > NL2
C          NU LIGNE SI LE SOMMET EST SUR UNE LIGNE UTILISATEUR
C          0        SINON
C NUISOP : NUMERO DE L'ISO DE CHAQUE POINT ET 0 SINON
C
C L1LGFR : POUR CHAQUE LIGNE DE NUMILF A NUMXLF POINTE SUR LA PREMIERE
C          ARETE CHAINEE DE LA LIGNE
C
C IERR   : 0 SI PAS D'ERREUR DETECTEE
C          1 SATURATION DU TABLEAU NSARCH
C          2 SATURATION DU TABLEAU XYZST
C          3 SATURATION DU TABLEAU NSNVTR
C
C TABLEAUX AUXILIAIRES :
C ----------------------
C SOLEL  : TABLEAU AUXILIAIRE
C COPOE  : TABLEAU AUXILIAIRE
C LTPILE : PILE DES 3 SOMMETS DES SOUS-TRIANGLES
C LAPILE : PILE DES 3 NO DE LIGNES DES ARETES DES SOUS-TRIANGLES
C
C XYVAL  : XYVAL(1 ET 2,.) 2 COORDONNEES DES SOMMETS DES SOUS TRIANGLES
C          XYVAL(3,.) VALEUR DE LA SOLUTION EN CE SOMMET
C FBASE  : TABLEAU AUXILIAIRE (CF SP INTERP)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     JANVIER 1995
C23456---------------------------------------------------------------012
      PARAMETER   ( LIGCON=0, LIGTIR=1 )
      include"./incl/gsmenu.inc"
      include"./incl/donele.inc"
      include"./incl/inteel.inc"
      include"./incl/ponoel.inc"
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
C
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
C
      INTEGER           NUNDEL(NBELEM,NBNOEL),
     %                  NUPTEL(NBELEM,NBPOEL)
      REAL              XYZPOI(3,NBPOIT)
      DOUBLE PRECISION  TEMPER(NTDL,NDSM),
     %                  SOLEL(*),
     %                  COPOE(NBPOEL,2)
C
      DOUBLE PRECISION  XD, YD, ZD, XREF(3), PROSCD, FBASE(30)
C
      REAL              VALISO(1:NBISO),
     %                  XYVAL(3,MXSOMM),
     %                  XYZST(3,MXYZST)
C
      INTEGER           LTPILE(3,MXPILE),
     %                  LAPILE(3,MXPILE),
     %                  NOXYZ(MXSOMM),
     %                  NSNVTR(4,MXNVTR),
     %                  L1LGFR(NUMILF:NUMXLF),
     %                  NUPLIS(MXYZST),
     %                  NUISOP(MXYZST),
     %                  NSARCH(4,MXARCH)
C
      INTEGER           NOSOEL(64)
      INTEGER           NOOBVC, NOOBSF(6), NOOBLA(12), NOOBPS(8)
C
      REAL              XYZ(3),
     %                  COSOTR(2,4)
C
      EQUIVALENCE      (XD,XREF(1)),(YD,XREF(2)),(ZD,XREF(3))
      DATA              ZD /0.D0/
C
C     RETROUVER LE NUMERO DES INTERPOLATIONS DE LA TEMPERATURE
C     ET DES COMPOSANTES DE LA TRANSFORMATION: EF UNITE -> EF COURANT
      CALL ELINTE( 'THERMIQUE', NUTYEL, NDIMF, NOINTF,
     &              NBINVA, NUINVA, NUINTI, NBNDIN )
C     LA TEMPERATURE EST ICI LA SEULE INCONNUE VARIATIONNELLE
      NOINTE = NUINTI(1)
C
C     LA DIMENSION DE L'ESPACE EST 2
      IF( NDIMF .NE. 2 ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'INCOMPATIBILITE DIMENSION ESPACES POUR F ET X'
         CALL LEREUR
         IERR = 10
         RETURN
      ENDIF
C
      IF( NARET .LT. 3  .OR.  NARET .GT. 4 ) THEN
C        ELEMENT NON TRIANGULAIRE OU QUADRANGULAIRE
         NBLGRC(NRERR) = 2
         KERR(1) = 'ERREUR : ELEMENT FINI 2D'
         KERR(2) = 'NON TRIANGLE NON QUADRANGLE'
         CALL LEREUR
         IERR = 11
         RETURN
      ENDIF
C
C     LE TRACE DES ARETES DU MAILLAGE DES EF DE CE TYPE D'EF
C     ------------------------------------------------------
      CALL XVTYPETRAIT( LIGTIR )
      CALL XVEPAISSEUR( 1 )
      CALL TRAREF( 2, NBPOIT, XYZPOI, NBELEM, NUPTEL )
      CALL XVTYPETRAIT( LIGCON )
C
C     OUVERTURE DE LA PILE DU NO DES 3 SOMMETS DES SOUS-TRIANGLES
C     -----------------------------------------------------------
C     LHPIL0 POINTE SUR LE SOMMET DE LA PILE LTPILE
      LHPIL0 = 0
C     NOSOM0 NO DU DERNIER SOMMET GENERE DANS LES SOUS-TRIANGLES
      NOSOM0 = 0
C     NBSOMM NOMBRE DE SOMMETS D'UN TRIANGLE
      NBSOMM = ( NBAEFU + 1 ) * ( NBAEFU + 2 ) / 2
C
C     ===========================================================
C     GENERATION DES SOUS-TRIANGLES DU TRIANGLE (0,0),(1,0),(0,1)
C     ===========================================================
C     LES COORDONNEES DU TRIANGLE DE REFERENCE
C     ----------------------------------------
      I = 1
      COSOTR(1,1) = 0.
      COSOTR(2,1) = 0.
      COSOTR(1,2) = 1.
      COSOTR(2,2) = 0.
      COSOTR(1,3) = 0.
      COSOTR(2,3) = 1.
C
C     LE NO DES 3 SOMMETS DES NBAEFU**2 SOUS-TRIANGLES DANS
C     LE TRIANGLE REFERENCE
C     -----------------------------------------------------
 20   CALL EMPITR( NBAEFU, NOSOM0, MXSOMM, MXPILE, LHPIL0, LTPILE )
C
C     LES COORDONNEES DES (NBAEFU+1)*(NBAEFU+2)/2 SOMMETS DES SOUS-TRIANGLES
C     ----------------------------------------------------------------------
      CALL COPOTR( COSOTR, NBAEFU+1, 3, XYVAL(1,NOSOM0+1) )
      NOSOM0 = NOSOM0 + NBSOMM
C
C     ============================================================
C     SI L ELEMENT EST UN QUADRANGLE GENERATION DES SOUS-TRIANGLES
C     DU TRIANGLE COMPLEMENTAIRE DE SOMMETS (1,1), (0,1), (1,0)
C     ============================================================
      IF( NARET .EQ. 4  .AND.  I .EQ. 1 ) THEN
         I = 2
         COSOTR(1,1) = 1.
         COSOTR(2,1) = 1.
         COSOTR(1,2) = 0.
         COSOTR(2,2) = 1.
         COSOTR(1,3) = 1.
         COSOTR(2,3) = 0.
         GOTO 20
      ENDIF
C
      IF( NOSOM0 .GT. MXSOMM ) THEN
         NBLGRC(NRERR) = 2
         KERR(1) = 'ADAP21: SATURATION DES PILES'
         KERR(2) = 'AUGMENTER LA VALEUR DE MXSOMM'
         CALL LEREUR
         IERR = 4
         RETURN
      ENDIF
C
C     =====================================
C     BOUCLE SUR LES EF COURANTS DE CE TYPE
C     =====================================
      DO 1000 NUELEM=1,NBELEM
C
C        LES NOEUDS DE L'ELEMENT FINI NUELEM
C        -----------------------------------
         CALL EFNOEU( MNNPEF, NUELEM, NBNDEL, NOSOEL )
C
C        LE NUMERO DE SURFACE DE L'EF
C        LE NUMERO DE SURFACE DES FACES   DE L'EF
C        LE NUMERO DE LIGNE   DES ARETES  DE L'EF
C        LE NUMERO DE POINT   DES SOMMETS DE L'EF
C        ----------------------------------------
         CALL EFPLSV( MNNPEF, NUELEM ,
     %                NOVCEL, NOSFEL, NOLAEL, NOPSEL ,
     %                NOOBVC, NOOBSF, NOOBLA, NOOBPS )
C
C        LE NUMERO DE LIGNE DES ARETES DU OU DES 2 TRIANGLES UNITES
C        EST EMPILE POUR LES SOUS-TRIANGLES
         LHPILA = 0
         IF( NARET .EQ. 3 ) THEN
C           UN SEUL TRIANGLE
            CALL EMPITA( NBAEFU, NOOBLA, MXPILE, LHPILA, LAPILE )
         ELSE
C           QUADRANGLE : LE TRIANGLE UNITE
            NOOBLA(5) = NOOBLA(1)
            NOOBLA(6) = 0
            NOOBLA(7) = NOOBLA(4)
            CALL EMPITA( NBAEFU, NOOBLA(5), MXPILE, LHPILA, LAPILE )
C           LE TRIANGLE SYMETRIQUE
            NOOBLA(5) = NOOBLA(3)
            NOOBLA(6) = 0
            NOOBLA(7) = NOOBLA(2)
            CALL EMPITA( NBAEFU, NOOBLA(5), MXPILE, LHPILA, LAPILE )
         ENDIF
C
C        FORMATION DES TABLEAUX COPOE(NBPOEL,2) ET SOLEL(NBNOELI)
C        COORDONNEES DES POINTS ET VALEUR DE LA SOLUTION NOSOL AUX NOEUDS
C        ================================================================
C        EXTRACTION DE LA TEMPERATURE NCAS DE TEMPER AUX NOEUDS DE L'EF
         DO 25 I=1,NBNOEL
            SOLEL( I ) = TEMPER( NUNDEL(NUELEM,I), NCAS )
 25      CONTINUE
C
C        EXTRACTION DES COORDONNEES DES POINTS DE L'EF
         DO 35 I=1,NBPOEL
            DO 30 J=1,2
               COPOE(I,J) = XYZPOI( J, NUPTEL(NUELEM,I) )
 30         CONTINUE
 35      CONTINUE
C
C        ICI LE TABLEAU SOLEL DOUBLE PRECISION CONTIENT LA VALEUR
C        DE LA SOLUTION AUX NOEUDS DE SON INTERPOLATION DANS L EF
C        DE REFERENCE ET LE TABLEAU COPOE(NBPOEL,2) LES COORDONNEES
C        DES POINTS DANS L'EF DE REFERENCE OU L INTERPOLATION EST EFFECTUEE
C
C        BOUCLE SUR LES SOMMETS DES SOUS-TRIANGLES DANS L'EF COURANT
C        VALEUR AUX SOMMETS DES SOUS-TRIANGLES DANS L'EF
C        ===========================================================
C        NOSOM0 DESIGNE LE NOMBRE DE SOMMETS DES SOUS-TRIANGLES
C        ATTENTION: PLUSIEURS SOMMETS ONT MEMES COORDONNEES
C                   CETTE IDENTIFICATION EST FAITE DANS XYZST
         DO 60 I=1,NOSOM0
C
C           LES 2 COORDONNEES DU SOMMET SUR L'EF UNITE
            XD = XYVAL(1,I)
            YD = XYVAL(2,I)
C
C           LA VALEUR DES NBNOEL FONCTIONS DE BASE EN (XD,YD,0D0)
            CALL INTERP( NOINTE, XD, YD, ZD, NBNOEL, FBASE )
C
C           LA VALEUR DE LA SOLUTION EN CE POINT
            XYVAL(3,I) = REAL( PROSCD( FBASE, SOLEL, NBNOEL ) )
C
C           PASSAGE A L'EF COURANT
C           (ATTENTION: ICI NBNOEL DOIT ETRE EGAL A NBPOEL
C                       ET POINTS=NOEUDS PAR EXEMPLE CAS LAGRANGE )
C           CALCUL DES COORDONNEES DES SOMMETS DES ARETES ISOVALEURS
C           A PARTIR DE LA TRANSFORMATION F: EF REFERENCE -> EF COURAN
            DO 40 L=1,2
               XYZ(L) = REAL( PROSCD( FBASE, COPOE(1,L) , NBPOEL ) )
 40         CONTINUE
            XYZ(3) = 0.0
C
C           IDENTIFICATION DE CE POINT PARMI XYZST
            DO 50 L=NBXYZ,1,-1
                CALL XYZIDE( XYZ, XYZST(1,L), IDENTQ )
                IF( IDENTQ .NE. 0 ) THEN
C                   POINT RETROUVE
                    NOXYZ( I ) = L
                    GOTO 60
                ENDIF
 50         CONTINUE
C
C           POINT NON RETROUVE => IL EST AJOUTE
            IF( NBXYZ .GE. MXYZST ) THEN
C              SATURATION DU TABLEAU XYZST
               NBLGRC(NRERR) = 1
               KERR(1) = 'TRSSTR:SATURATION DU TABLEAU XYZST'
               CALL LEREUR
               IERR = 2
               RETURN
            ENDIF
            NBXYZ = NBXYZ + 1
            NOXYZ( I ) = NBXYZ
            XYZST(1,NBXYZ) = XYZ(1)
            XYZST(2,NBXYZ) = XYZ(2)
            XYZST(3,NBXYZ) = XYZ(3)
C           LE SOMMET EST A PRIORI NI UN SOMMET INITIAL NI SUR UNE ARETE ISOVALE
CCC         NUPLIS( NBXYZ ) = 0   (FAIT DANS ADAP2D)
C
 60      CONTINUE
C
C        IDENTIFICATION DES SOMMETS DE L'EF QUI SONT DES POINTS
C        UTILISATEUR ET DONC A IMPOSER SANS DEPLACEMENT OU SUPPRESSION
C        -------------------------------------------------------------
         IF( NARET .EQ. 3 ) THEN
C
C           TRIANGLE
C           LE SOMMET 3 EST LE SOUS-SOMMET 1
            IF( NOOBPS(3) .GT. 0 ) THEN
C              LE NUMERO DE POINT EST IMPOSE
               NUPLIS( NOXYZ(1) ) = -NOOBPS(3)
            ENDIF
C           LE SOMMET 2 EST LE SOUS-SOMMET NBSOMM
            IF( NOOBPS(2) .GT. 0 ) THEN
C              LE NUMERO DE POINT EST IMPOSE
               NUPLIS( NOXYZ(NBSOMM) ) = -NOOBPS(2)
            ENDIF
C           LE SOMMET 1 EST LE SOUS-SOMMET NBSOMM - NBAEFU
            IF( NOOBPS(1) .GT. 0 ) THEN
C              LE NUMERO DE POINT EST IMPOSE
               NUPLIS( NOXYZ(NBSOMM-NBAEFU) ) = -NOOBPS(1)
            ENDIF
C
         ELSE
C
C           QUADRANGLE
C           LE SOMMET 4 EST LE SOUS-SOMMET 1
            IF( NOOBPS(4) .GT. 0 ) THEN
C              LE NUMERO DE POINT EST IMPOSE
               NUPLIS( NOXYZ(1) ) = -NOOBPS(4)
            ENDIF
C           LE SOMMET 3 EST LE SOUS-SOMMET NBSOMM + NBSOMM - NBAEFU
            IF( NOOBPS(3) .GT. 0 ) THEN
C              LE NUMERO DE POINT EST IMPOSE
               NUPLIS( NOXYZ(NBSOMM+NBSOMM-NBAEFU) ) = -NOOBPS(3)
            ENDIF
C           LE SOMMET 2 EST LE SOUS-SOMMET NBSOMM
            IF( NOOBPS(2) .GT. 0 ) THEN
C              LE NUMERO DE POINT EST IMPOSE
               NUPLIS( NOXYZ(NBSOMM) ) = -NOOBPS(2)
            ENDIF
C           LE SOMMET 1 EST LE SOUS-SOMMET NBSOMM-NBAEFU
            IF( NOOBPS(1) .GT. 0 ) THEN
C              LE NUMERO DE POINT EST IMPOSE
               NUPLIS( NOXYZ(NBSOMM-NBAEFU) ) = -NOOBPS(1)
            ENDIF
         ENDIF
C
C        IDENTIFICATION DES SOMMETS DE L'EF QUI APPARTIENNENT A
C        2 LIGNES DE L'UTILISATEUR (FRONTIERE OU INTERFACE) (SOMMET NON SUPPRIMA
C        -----------------------------------------------------------------------
         IF( NARET .EQ. 3 ) THEN
C
C           TRIANGLE
C           LE SOMMET 3 EST LE SOUS-SOMMET 1
            IF( NOOBPS(3) .EQ. 0 .AND.
     %          NOOBLA(2) .GT. 0 .AND. NOOBLA(3) .GT. 0 ) THEN
C              LE NUMERO DE POINT FICTIF -(1000000+NL1+1000*NL2)
               NL1 = NOOBLA(2)
               NL2 = NOOBLA(3)
               IF( NL2 .GT. NL1 ) THEN
                  NL1 = NOOBLA(3)
                  NL2 = NOOBLA(2)
               ENDIF
               NUPLIS( NOXYZ(1) )= -1 000 000 - NL1 - 1000 * NL2
            ENDIF
C           LE SOMMET 2 EST LE SOUS-SOMMET NBSOMM
            IF( NOOBPS(2) .EQ. 0 .AND.
     %          NOOBLA(1) .GT. 0 .AND. NOOBLA(2) .GT. 0 ) THEN
C              LE NUMERO DE POINT EST IMPOSE
               NL1 = NOOBLA(1)
               NL2 = NOOBLA(2)
               IF( NL2 .GT. NL1 ) THEN
                  NL1 = NOOBLA(2)
                  NL2 = NOOBLA(1)
               ENDIF
               NUPLIS( NOXYZ(NBSOMM) )= -1 000 000 - NL1 - 1000 * NL2
            ENDIF
C           LE SOMMET 1 EST LE SOUS-SOMMET NBSOMM - NBAEFU
            IF( NOOBPS(1) .EQ. 0 .AND.
     %          NOOBLA(3) .GT. 0 .AND. NOOBLA(1) .GT. 0 ) THEN
C              LE NUMERO DE POINT EST IMPOSE
               NL1 = NOOBLA(1)
               NL2 = NOOBLA(3)
               IF( NL2 .GT. NL1 ) THEN
                  NL1 = NOOBLA(3)
                  NL2 = NOOBLA(1)
               ENDIF
               NUPLIS( NOXYZ(NBSOMM-NBAEFU) )= -1 000 000 -NL1 -1000*NL2
            ENDIF
C
         ELSE
C
C           QUADRANGLE
C           LE SOMMET 4 EST LE SOUS-SOMMET 1
            IF( NOOBPS(4) .EQ. 0 .AND.
     %          NOOBLA(3) .GT. 0 .AND. NOOBLA(4) .GT. 0 ) THEN
C              LE NUMERO DE POINT EST IMPOSE
               NL1 = NOOBLA(3)
               NL2 = NOOBLA(4)
               IF( NL2 .GT. NL1 ) THEN
                  NL1 = NOOBLA(4)
                  NL2 = NOOBLA(3)
               ENDIF
               NUPLIS( NOXYZ(1) )= -1 000 000 - NL1 - 1000 * NL2
            ENDIF
C           LE SOMMET 3 EST LE SOUS-SOMMET NBSOMM + NBSOMM - NBAEFU
            IF( NOOBPS(3) .EQ. 0 .AND.
     %          NOOBLA(2) .GT. 0 .AND. NOOBLA(3) .GT. 0 ) THEN
C              LE NUMERO DE POINT EST IMPOSE
               NL1 = NOOBLA(2)
               NL2 = NOOBLA(3)
               IF( NL2 .GT. NL1 ) THEN
                  NL1 = NOOBLA(3)
                  NL2 = NOOBLA(2)
               ENDIF
               NUPLIS(NOXYZ(NBSOMM+NBSOMM-NBAEFU))=-1000000-NL1-1000*NL2
            ENDIF
C           LE SOMMET 2 EST LE SOUS-SOMMET NBSOMM
            IF( NOOBPS(2) .EQ. 0 .AND.
     %          NOOBLA(1) .GT. 0 .AND. NOOBLA(2) .GT. 0 ) THEN
C              LE NUMERO DE POINT EST IMPOSE
               NL1 = NOOBLA(1)
               NL2 = NOOBLA(2)
               IF( NL2 .GT. NL1 ) THEN
                  NL1 = NOOBLA(2)
                  NL2 = NOOBLA(1)
               ENDIF
               NUPLIS( NOXYZ(NBSOMM) ) = -1000000-NL1-1000*NL2
            ENDIF
C           LE SOMMET 1 EST LE SOUS-SOMMET NBSOMM-NBAEFU
            IF( NOOBPS(1) .EQ. 0 .AND.
     %          NOOBLA(4) .GT. 0 .AND. NOOBLA(1) .GT. 0 ) THEN
C              LE NUMERO DE POINT EST IMPOSE
               NL1 = NOOBLA(1)
               NL2 = NOOBLA(4)
               IF( NL2 .GT. NL1 ) THEN
                  NL1 = NOOBLA(4)
                  NL2 = NOOBLA(1)
               ENDIF
               NUPLIS(NOXYZ(NBSOMM-NBAEFU))=-1000000-NL1-1000*NL2
            ENDIF
C
         ENDIF
C
C        IDENTIFICATION DES SOMMETS DE L'EF QUI APPARTIENNENT A
C        1 SEULE LIGNE UTILISATEUR (ET/OU FRONTIERE OU INTERFACE)
C        --------------------------------------------------------
         IF( NARET .EQ. 3 ) THEN
C
C           TRIANGLE ARETE 1
            IF( NOOBLA(1) .GT. 0 ) THEN
               N = NBSOMM + 1
               DO 11 I=1,NBAEFU+1
C                 LE NUMERO DU SOMMET DANS L'EF
                  K = N - I
C                 LE NUMERO DANS PXYD
                  K = NOXYZ( K )
                  IF( NUPLIS(K) .EQ. 0  ) NUPLIS(K) = NOOBLA(1)
 11            CONTINUE
            ENDIF
C
C           TRIANGLE ARETE 2
            IF( NOOBLA(2) .GT. 0 ) THEN
               DO 12 I=1,NBAEFU+1
C                 LE NUMERO DU SOMMET DANS L'EF
                  K = I * (I+1) / 2
C                 LE NUMERO DANS PXYD
                  K = NOXYZ( K )
                  IF( NUPLIS(K) .EQ. 0  ) NUPLIS(K) = NOOBLA(2)
 12            CONTINUE
            ENDIF
C
C           TRIANGLE ARETE 3
            IF( NOOBLA(3) .GT. 0 ) THEN
               DO 13 I=1,NBAEFU+1
C                 LE NUMERO DU SOMMET DANS L'EF
                  K = I * (I-1) / 2 + 1
C                 LE NUMERO DANS PXYD
                  K = NOXYZ( K )
                  IF( NUPLIS(K) .EQ. 0  ) NUPLIS(K) = NOOBLA(3)
 13            CONTINUE
            ENDIF
C
         ELSE
C
C           QUADRANGLE ARETE 1
            IF( NOOBLA(1) .GT. 0 ) THEN
               N = NBSOMM + 1
               DO 21 I=1,NBAEFU+1
C                 LE NUMERO DU SOMMET DANS L'EF
                  K = N - I
C                 LE NUMERO DANS PXYD
                  K = NOXYZ( K )
                  IF( NUPLIS(K) .EQ. 0  ) NUPLIS(K) = NOOBLA(1)
 21            CONTINUE
            ENDIF
C
C           QUADRANGLE ARETE 2
            IF( NOOBLA(2) .GT. 0 ) THEN
               DO 22 I=1,NBAEFU+1
C                 LE NUMERO DU SOMMET DANS L'EF
                  K = NBSOMM + I * (I-1) / 2 + 1
C                 LE NUMERO DANS PXYD
                  K = NOXYZ( K )
                  IF( NUPLIS(K) .EQ. 0  ) NUPLIS(K) = NOOBLA(2)
 22            CONTINUE
            ENDIF
C
C           QUADRANGLE ARETE 3
            IF( NOOBLA(3) .GT. 0 ) THEN
               N = NBSOMM + NBSOMM + 1
               DO 23 I=1,NBAEFU+1
C                 LE NUMERO DU SOMMET DANS L'EF
                  K = N - I
C                 LE NUMERO DANS PXYD
                  K = NOXYZ( K )
                  IF( NUPLIS(K) .EQ. 0  ) NUPLIS(K) = NOOBLA(3)
 23            CONTINUE
            ENDIF
C
C           QUADRANGLE ARETE 4
            IF( NOOBLA(4) .GT. 0 ) THEN
               DO 24 I=1,NBAEFU+1
C                 LE NUMERO DU SOMMET DANS L'EF
                  K = I * (I-1) / 2 + 1
C                 LE NUMERO DANS PXYD
                  K = NOXYZ( K )
                  IF( NUPLIS(K) .EQ. 0  ) NUPLIS(K) = NOOBLA(4)
 24            CONTINUE
            ENDIF
C
         ENDIF
C
C        LA HAUTEUR INITIALE DE LA PILE
         LHPILE = LHPIL0
C
C        ********************************************************
C        TANT QUE LA PILE DES SOUS-TRIANGLES N EST PAS VIDE FAIRE
C        ********************************************************
 100     IF( LHPILE .GT. 0 ) THEN
C
C           TRIANGULATION DU SOUS-TRIANGLE LHPILE
            CALL TRSSTR( NBISO,  VALISO, DFMXIS,
     %                   LTPILE(1,LHPILE), NOOBSF(1),
     %                   LAPILE(1,LHPILE), XYVAL,
     %                   NOINTF, NBPOEL, COPOE,  FBASE,
     %                   MXYZST, NBXYZ,  XYZST,  NOXYZ,
     %                   MXNVTR, NBNVTR, NSNVTR, NUPLIS, NUISOP,
     %                   NUMILF, NUMXLF, L1LGFR,
     %                   MXARCH, NBARCH, N1ARCH, NSARCH, IERR )
            IF( IERR .NE. 0 ) RETURN
C
C           LE SOUS-TRIANGLE EST DEPILE
            LHPILE = LHPILE - 1
C           REMONTEE POUR TRAITER LE SOUS-TRIANGLE SUIVANT DE L'EF UNITE
            GOTO 100
C
         ENDIF
C
 1000 CONTINUE
C
C     RETOUR A LA VALEUR PAR DEFAUT
      CALL XVEPAISSEUR( 1 )
      CALL XVTYPETRAIT( LIGCON )
      END
