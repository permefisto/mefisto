      SUBROUTINE SDTRIS( NBISO , VALISO, NTDISO ,
     %                   NDIM  , NCAS  , NDSM  , NTDL  , TEMPER ,
     %                   NUTYEL, NBELEM ,
     %                   NBNOEL, NUNDEL, NBPOEL, NUPTEL ,
     %                   NBPOIT, XYZPOI ,
     %                   MXSOMM, MXPOIS, SOLEL , COPOE  ,
     %                   MXPIL3, NPILE3, NOPOTR ,
     %                   NUINVA, NUINTI, NBNDIN,
     %                   COPOT , VALST , FBASE  ,
     %                   TMI, XMI, YMI, TMA, XMA, YMA )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACER LES ISOTHERMES DES ELEMENTS D'UN TYPE
C -----    TRIANGLE OU QUADRANGLE  ( VERSION SOUS-DOMAINES )
C
C ENTREES :
C ---------
C NBISO  : NOMBRE D'ISOVALEURS A TRACER
C VALISO : TABLEAU DES VALEURS DES ISO A TRACER
C AMPLID : FACTEUR D'AMPLIFICATION DES TEMPERATURES
C NDIM   : ESPACE DE TRAVAIL 2 OU 3
C NCAS   : NUMERO DU CAS A TRAITER
C NDSM   : NOMBRE TOTAL DE VECTEURS TEMPERATURES
C NTDL   : NOMBRE TOTAL DE DEGRES DE LIBERTE EN TEMPERATURE
C TEMPER : LES NTDL * NDSM TEMPERATURES
C NBELEM : NOMBRE D'ELEMENTS DE CE TYPE
C NBNOEL : NOMBRE DE  NOEUDS DE L'ELEMENT
C NUNDEL : NUMERO DES NOEUDS DES  ELEMENTS
C NBPOEL : NOMBRE DE  POINTS DE L'ELEMENT
C NUPTEL : NUMERO DES POINTS DES  ELEMENTS
C NBNOEU : NOMBRE TOTAL DE NOEUDS DU MAILLAGE
C XYZNOE : COORDONNEES DES NOEUDS DU MAILLAGE
C NBPOIT : NOMBRE TOTAL DE POINTS DU MAILLAGE
C XYZPOI : COORDONNEES DES POINTS DU MAILLAGE
C
C TABLEAUX AUXILIAIRES :
C ----------------------
C MXSOMM : NBRE MAXIMUM DE SOMMETS DES SOUS-TRIANGLES GENERES
C MXPOIS : NBRE MAXIMUM DE POINTS DES ISOVALEURS SUR LE TRIANGLE
C SOLEL  : TABLEAU AUXILIAIRE
C COPOE  : TABLEAU AUXILIAIRE
C MXPIL3 : NBRE MAXIMUM DE SOUS-TRIANGLES DANS LA PILE NPILE3
C NPILE3 : PILE DES 3 SOMMETS DES SOUS-TRIANGLES
C NOPOTR : NUMERO DES POINTS SUR LES ISO D'UN TRIANGLE
C NUINVA,NUINTI,NBNDIN,NUNOIN : TABLEAUX INTERMEDIAIRES (CF SP ELINTE)
C COPOT  : COORDONNEES DES POINTS ISO SUR LE  TRIANGLE
C VALST  : VALST(1 ET 2,.) 2 COORDONNEES DES SOMMETS DES SOUS TRIANGLES
C          VALST(3,.) VALEUR DE LA SOLUTION EN CE SOMMET
C FBASE  : TABLEAU AUXILIAIRE (CF SP INTERP)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS   NOVEMBRE 1994
C23456---------------------------------------------------------------012
      PARAMETER      ( LIGCON=0, LIGTIR=1 )
      include"./incl/gsmenu.inc"
      include"./incl/donele.inc"
      include"./incl/inteel.inc"
      include"./incl/ponoel.inc"
      include"./incl/trvari.inc"
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      INTEGER           NUNDEL(NBELEM,NBNOEL),NUPTEL(NBELEM,NBPOEL),
     %                  NPILE3(3,MXPIL3),NOPOTR(3,NBISO),NTDISO(1:NBISO)
      INTEGER           NOS(6)
      REAL              XYZPOI(3,NBPOIT)
      DOUBLE PRECISION  TEMPER(NTDL,NDSM),SOLEL(*),COPOE(NBPOEL,NDIM)
      DOUBLE PRECISION  XD,YD,ZD,PROSCD,FBASE(30),XREF(3)
      REAL              VALISO(1:NBISO),VALST(3,MXSOMM),COPOT(2,MXPOIS)
      REAL              COSOTR(2,3),C(3),C1(3)
      CHARACTER*8       KNOM
      EQUIVALENCE      (XD,XREF(1)),(YD,XREF(2)),(ZD,XREF(3))
      DATA              ZD  /0.D0/
C
C     CALCUL DU NOMBRE DE SUBDIVISIONS DE (0,1)
C     EN FONCTION DU NOMBRE DE NOEUDS DE L ELEMENT
      NBSOUI = 1
      IF( NBNOE .GT. 4 ) NBSOUI = 2
C
C     RETROUVER LE NUMERO DES INTERPOLATIONS DE LA TEMPERATURE
C     ET DES COMPOSANTES DE LA TRANSFORMATION:ELT REFERENCE->ELEMENT
      CALL ELINTE('THERMIQUE',NUTYEL,NDIMF,NOINTF,
     &             NBINVA,NUINVA,NUINTI,NBNDIN)
C     NUINVA(5),NUINTI(5),NBNDIN(5),NUNOIN(30,5) (CF SP ELINTE)
C     LA TEMPERATURE EST ICI LA SEULE INCONNUE VARIATIONNELLE
      NOINTE = NUINTI(1)
C
CCC      IDEM POUR LE GRADIENT DE TEMPERATURE
CCC      CALL ELINT1('THERMIQUE',NUTYEL,NDIMF,NOINFD,NBPOF,
CCC     &            NBINVA,NUINVA,NUINTI,NBNDIN,NUNOIN)
C
C     LA DIMENSION DE TRACE EST NDIM
      IF( NDIM .NE. NDIMF ) THEN
         NBLGRC(NRERR) = 1
         KERR(1) = 'INCOMPATIBILITE DIMENSION ESPACES POUR F ET X'
         CALL LEREUR
         RETURN
      ENDIF
C
      IF( NARET .LT. 3  .OR.  NARET .GT. 4 ) THEN
C        ELEMENT NON TRIANGULAIRE OU QUADRANGULAIRE
         NBLGRC(NRERR) = 2
         KERR(1) = 'ERREUR : ELEMENT FINI 2D'
         KERR(2) = 'NON TRIANGLE NON QUADRANGLE'
         CALL LEREUR
         RETURN
      ENDIF
C
C     LE TRACE DES ARETES
C     -------------------
      CALL XVTYPETRAIT( LIGTIR )
      CALL XVEPAISSEUR( 1 )
      CALL TRAREF( NDIM, NBPOIT, XYZPOI, NBELEM, NUPTEL )
C
C
C     TRACE EFFECTIF DES FLECHES DU GRADIENT DE TEMPERATURE
C     -----------------------------------------------------
C     LIGNES NON EPAISSIES
      CALL XVEPAISSEUR( 2 )
      CALL XVTYPETRAIT( LIGCON )
C
      DO 1000 K=1,NBELEM
C
C        OUVERTURE DE LA PILE DU NO DES 3 SOMMETS DES SOUS-TRIANGLES
C        -----------------------------------------------------------
C        IASPI3 : POINTE SUR LE SOMMET DE LA PILE NPILE3
C        NOSOMM  : NO DU DERNIER SOMMET GENERE DANS LES SOUS-TRIANGLES
C
         IASPI3 = 0
         NOSOMM = 0
C
C        GENERATION DES SOUS-TRIANGLES DU TRIANGLE (0,0),(1,0),(0,1)
C        =============================================================
C        LES COORDONNEES DE L ELEMENT DE REFERENCE
C        -----------------------------------------
         I           = 1
         COSOTR(1,1) = 0.
         COSOTR(2,1) = 0.
         COSOTR(1,2) = 1.
         COSOTR(2,2) = 0.
         COSOTR(1,3) = 0.
         COSOTR(2,3) = 1.
C
C        LE NO DES 3 SOMMETS DE CHAQUE TRIANGLE
C        --------------------------------------
   20    CALL EMPITR( NBSOUI,NOSOMM,MXSOMM,MXPIL3,IASPI3,NPILE3 )
C
C        LES COORDONNEES DES (NBSOUI+1) ** 2 SOMMETS DES SOUS-TRIANGLES
C        --------------------------------------------------------------
         CALL COPOTR( COSOTR, NBSOUI+1, 3, VALST(1,NOSOMM+1) )
         NOSOMM = NOSOMM + (NBSOUI+2) * (NBSOUI+1) / 2
C
C        SI L ELEMENT EST UN QUADRANGLE GENERATION DES SOUS-TRIANGLES
C        DU TRIANGLE (1,1), (0,1), (1,0)
C        ============================================================
         IF( NARET .EQ. 4  .AND.  I .EQ. 1 ) THEN
            I           = 2
            COSOTR(1,1) = 1.
            COSOTR(2,1) = 1.
            COSOTR(1,2) = 0.
            COSOTR(2,2) = 1.
            COSOTR(1,3) = 1.
            COSOTR(2,3) = 0.
            GOTO 20
         ENDIF
C
C        FORMATION DES TABLEAUX COPOE(NBPOEL,NDIM) ET SOLEL(NBNOELI)
C        COORDONNEES DES POINTS ET VALEUR DE LA SOLUTION NOSOL AUX NOEUDS
C        ================================================================
C        EXTRACTION DE LA TEMPERATURE NCAS DE TEMPER
         DO 35 I=1,NBNOEL
            SOLEL( I ) = TEMPER( NUNDEL(K,I), NCAS )
   35    CONTINUE
C
C        EXTRACTION DES COORDONNEES DES POINTS DE L'ELEMENT
         DO 38 I=1,NBPOEL
            DO 37 J=1,NDIM
               COPOE(I,J) = XYZPOI( J, NUPTEL(K,I) )
   37       CONTINUE
   38    CONTINUE
C
C        ICI LE TABLEAU SOLEL DOUBLE PRECISION CONTIENT LA VALEUR
C        DE LA SOLUTION AUX NOEUDS DE SON INTERPOLATION DANS L ELEMENT
C        DE REFERENCE ET LE TABLEAU COPOE(NBPOEL,NDIM) LES COORDONNEES
C        DES POINTS OU L INTERPOLATION EST EFFECTUEE
C
C        BOUCLE SUR LES SOMMETS DES SOUS-TRIANGLES
C        =========================================
         DO 78 I=1,NOSOMM
            XD = VALST(1,I)
            YD = VALST(2,I)
C
C           LA VALEUR DES NBN FONCTIONS DE BASE EN (XD,YD,0D0)
            CALL INTERP(NOINTE,XD,YD,ZD,NBN,FBASE)
C
C           LA VALEUR DE LA SOLUTION EN CE POINT
            VALST(3,I) = REAL( PROSCD( FBASE, SOLEL, NBN ) )
C           RECHERCHE DES TEMPERATURES EXTREMES
            IF (VALST(3,I).LT.TMI) THEN
               TMI = VALST(3,I)
               XREF(1) = XD
               XREF(2) = YD
               CALL FPOXYZ( NDIM,NBPOEL,NOINTF,COPOE,XREF,FBASE,C )
               XMI = C(1)
               YMI = C(2)
            ELSE IF (VALST(3,I).GT.TMA) THEN
               TMA = VALST(3,I)
               XREF(1) = XD
               XREF(2) = YD
               CALL FPOXYZ( NDIM,NBPOEL,NOINTF,COPOE,XREF,FBASE,C )
               XMA = C(1)
               YMA = C(2)
            ENDIF
   78    CONTINUE
C
C        **************************************
C        TANT QUE LA PILE3 N EST PAS VIDE FAIRE
C        **************************************
C
   80    IF( IASPI3 .GT. 0 ) THEN
C
C           LES 3 SOMMETS SONT DEPILES
C           --------------------------
            FMIN =  1.E30
            FMAX = -FMIN
            DO 90 I=1,3
               NOSI = NPILE3(I,IASPI3)
               NOS(I) = NOSI
               F1     = VALST(3,NOSI)
               IF( F1 .GT. FMAX ) FMAX = F1
               IF( F1 .LT. FMIN ) FMIN = F1
   90       CONTINUE
            IASPI3 = IASPI3 - 1
C
C           SI LES 3 VALEURS AUX SOMMETS DU TRIANGLE SONT EXTERIEURES
C           A L INTERVALLE DES ISOS LE TRIANGLE EST ABANDONNE
C           =========================================================
            IF( FMIN .GT. VALISO(NBISO) .OR. FMAX .LT. VALISO(1))GOTO 80
C
C           NOP NOMBRE DE POINTS DES ISOVALEURS SUR CE TRIANGLE
C           ---------------------------------------------------
            NOP = 0
C
C           MISE A ZERO DU TABLEAU NOPOTR
C           -----------------------------
            CALL AZEROI( 3*NBISO, NOPOTR )
C
C           KMIN,KMAX NO DE L ISOVALEUR MIN ET MAX SUR LE TRIANGLE
C           ------------------------------------------------------
            KMIN = NBISO + 1
            KMAX = 0
C
C           BOUCLE SUR LES 3 ARETES DU TRIANGLE
C           ===================================
            DO 200 I=1,3
C
C              NO DES EXTREMITES DE L ARETE
               N1 = NOS(I)
               IF( I .NE. 3 ) THEN
                  N2 = I + 1
               ELSE
                  N2 = 1
               ENDIF
               N2 = NOS(N2)
C
C               LA VALEUR DE LA SOLUTION AUX EXTREMITES
C
                F1 = VALST(3,N1)
                F2 = VALST(3,N2)
cc                IF( F1 - F2 ) 120, 100, 110
                IF( F1 .GT. F2 ) GOTO 110
                IF( F1 .LT. F2 ) GOTO 120
C
C 100           F1 = F2 LA SOLUTION EST ELLE STATIONNAIRE SUR L ELEMENT?
C               --------------------------------------------------------
                N3 = I + 2
                IF(N3 .GT. 3) N3 = N3 - 3
                N3 = NOS(N3)
                F3 = VALST(3,N3)
                IF( F1 .EQ. F3 ) GOTO 80
C
C               SI F1=F2=F3 PAS DE TRACE D ISOVALEUR.PASSAGE AU TRIANGLE
C                           SUIVANT  (=> 80)
C
C               SI F1=F2=/F3 LES 2 EXTREMITES FORMENT EVENTUELLEMENT UN
C                            SEGMENT D ISOVALEUR
C
                DO 105 KK=1,NBISO
                   IF( F1 .EQ. VALISO(KK) ) THEN
C                     F1=F2 EST LA KK-EME ISOVALEUR
C                     ADJONCTION DU SEGMENT N1-N2
C                     -----------------------------
C                     LES EXTREMITES DE L ARETE SONT AJOUTEES
                      CALL POKISO( KK,NOP,VALST(1,N1),COPOT,NOPOTR,IERR)
                      CALL POKISO( KK,NOP,VALST(1,N2),COPOT,NOPOTR,IERR)
                      IF( IERR .NE. 0 ) GOTO 300
                      KMIN = MIN0( KMIN, KK )
                      K2   = KK
                      GOTO 180
                   ENDIF
  105           CONTINUE
                GOTO 200
C
C               F1>F2 PERMUTATION DE (N1,N2) ET (F1,F2)
C               ---------------------------------------
  110           KK = N1
                N1 = N2
                N2 = KK
                F3 = F1
                F1 = F2
                F2 = F3
C
C               F1<F2 . RECHERCHE DE LA PLUS PETITE ISOVALEUR SUR
C                       LE SOUS-TRIANGLE
C               -----   -----------------------------------------
  120           DO 130 KK=1,NBISO
                   K1 = KK
                   IF( VALISO(KK) .GE. F1 ) GOTO 140
  130           CONTINUE
                GOTO 200
C
C               LE TRIANGLE INTERSECTE AU MOINS UNE ISOVALEUR
C               ---------------------------------------------
  140           KMIN = MIN0(KMIN,K1)
                K2   = K1
C
C               BOUCLE SUR LES ISOVALEURS K1 A NBISO
C               ------------------------------------
                DO 170 KK=K1,NBISO
                   IF( VALISO(KK) .GT. F2 ) GOTO 180
C
C                  L ISOVALEUR KK EST COMPRISE ENTRE F1 ET F2
                   K2 = KK
C
C                  ADJONCTION DU POINT D INTERSECTION ISO KK ET ARETE I
C                  ---------------------------------------------------
                   Z    = (VALISO(KK) - F1) / (F2 - F1)
                   C(1) = VALST(1,N1) + Z  * (VALST(1,N2)-VALST(1,N1))
                   C(2) = VALST(2,N1) + Z  * (VALST(2,N2)-VALST(2,N1))
                   CALL POKISO(KK,NOP,C,COPOT,NOPOTR,IERR)
                   IF( IERR .NE. 0 ) GOTO 300
  170           CONTINUE
C
C               NO DE LA PLUS GRANDE ISOVALEUR TRAITEE SUR LE TRIANGLE
  180           KMAX = MAX0(KMAX,K2)
C
C               FIN DE LA BOUCLE SUR LES ARETES DU TRIANGLE
  200        CONTINUE
C
C            INSERTION DES SEGMENTS DE CHAQUE ISOVALEUR DU TRIANGLE
C            DANS CEUX DE L ELEMENT DE REFERENCE COMPLET
C            ======================================================
             IF ( KMIN .LE. KMAX ) THEN
                NCO = NCBLAN
                DO 210 KK=KMIN,KMAX
                   IF(NOPOTR(1,KK) .LE. 1) GOTO 210
C                  LE TRACE DE L'ISOTHERME AVEC LA BONNE COULEUR
                   IF( NOIBLA .EQ. 0 ) THEN
C                     LA COULEUR DE L'ISOTHERME
                      NCO = MOD(KK,NDCOUL) + N1COUL + 1
                   ENDIF
C
C                  IL EXISTE UN SEGMENT A TRACER
C                  -----------------------------
                   N1 = NOPOTR(2,KK)
                   XREF(1) = COPOT(1,N1)
                   XREF(2) = COPOT(2,N1)
                   CALL FPOXYZ( NDIM,NBPOEL,NOINTF,COPOE,XREF,FBASE,C )
C
                   N2 = NOPOTR(3,KK)
                   XREF(1) = COPOT(1,N2)
                   XREF(2) = COPOT(2,N2)
                   CALL FPOXYZ( NDIM,NBPOEL,NOINTF,COPOE,XREF,FBASE,C1 )
                   CALL TRAIT2D( NCO, C(1), C(2), C1(1), C1(2) )
C
C                  TRACE DU NO DE L ISOVALEUR
C                  --------------------------
                   NTDISO( KK ) = NTDISO( KK ) + 1
                   IF( MOD( NTDISO(KK), 16 ) .EQ. 3 ) THEN
                      WRITE( KNOM, '(I8)' ) KK
                      DO 208 J=7,1,-1
                         IF( KNOM(J:J) .EQ. ' ' ) THEN
                            I1 = J + 1
                            GOTO 209
                         ENDIF
  208                 CONTINUE
                      I1 = 1
  209                 CALL TEXTE2D( NCO, (C(1)+C1(1))*0.5,
     %                             (C(2)+C1(2))*0.5, KNOM(I1:8) )
                   ENDIF
  210           CONTINUE
             ENDIF
C
C            REMONTEE POUR TRAITER LE SOUS-TRIANGLE SUIVANT
C            ==============================================
             GOTO 80
C
C
C            IL Y A AU MOINS 3 POINTS D INTERSECTION D UNE ISOVALEUR
C            AVEC LE TRIANGLE . IL EST SUBDIVISE A SON TOUR EN
C            4 SOUS-TRIANGLES
C            =============================================================
C            ACTUELLEMENT CE CAS NE PEUT SE PRODUIRE POUR P1.
C            SI LA RECHERCHE DU ZERO DE L INTERPOLATION EST FAITE
C            IL PEUT Y AVOIR 2 OU PLUSIEURS POINTS D INTERSECTION
C            CETTE RECHERCHE DES ZEROS RESTE A PROGRAMMER
  300        WRITE(IMPRIM,10300) IASPI3+1,(J,NOPOTR(1,J),J=1,NBISO)
10300 FORMAT('ERREUR ELEMENT',I7,' 3 POINTS D INTERSECTION ISO-ARETE'/
     &       5(' NOPOTR(1,',I4,')=',I5))
C
C            GENERATION DU MILIEU DES 3 ARETES DU SOUS-TRIANGLE
C            ==================================================
             DO 320 I=1,3
                NOSOMM = NOSOMM + 1
C
C               CONTROLE DE SATURATION
C               ----------------------
                IF( NOSOMM .GT. MXSOMM ) THEN
C                  TABLEAU VALST TROP PETIT.DIAGNOSTIC.ARRET
                   NBLGRC(NRERR) = 2
                   KERR(1) = 'TAILLE VALST INSUFFISANTE'
                   KERR(2) = 'AUGMENTER MXSOMM=               '
                   WRITE( KERR(2)(18:27),'(I10)') MXSOMM
                   CALL LEREUR
                   RETURN
                ENDIF
C
C               LE NO DU MILIEU DE L ARETEET SES COORDONNEES
C               --------------------------------------------
                NOS(I+3) = NOSOMM
                J = I + 1
                IF( J .GT. 3 ) J = 1
                DO 315 KK=1,2
                   VALST(KK,NOSOMM) = ( VALST(KK,NOS(I))
     &                              +   VALST(KK,NOS(J)) ) * 0.5
  315           CONTINUE
C
C               LA VALEUR DE LA SOLUTION EN CE MILIEU D ARETE
C               ---------------------------------------------
                XD = VALST(1,NOSOMM)
                YD = VALST(2,NOSOMM)
C               LA VALEUR DES FONCTIONS DE BASE EN CE POINT
                CALL INTERP( NOINTE,XD,YD,ZD,NBN,FBASE )
C               LA VALEUR DE LA SOLUTION EN CE POINT
                VALST(3,NOSOMM) = REAL( PROSCD( FBASE, SOLEL, NBN ) )
  320        CONTINUE
C
C            LES 4 SOUS-TRIANGLES SONT EMPILES DANS NPILE3
C            ============================================
             CALL EMPIL3( IASPI3,MXPIL3,NPILE3,NOS(1),NOS(4),NOS(6) )
             CALL EMPIL3( IASPI3,MXPIL3,NPILE3,NOS(4),NOS(2),NOS(5) )
             CALL EMPIL3( IASPI3,MXPIL3,NPILE3,NOS(5),NOS(6),NOS(4) )
             CALL EMPIL3( IASPI3,MXPIL3,NPILE3,NOS(6),NOS(5),NOS(3) )
             GOTO 80
C
C            ****************************************************
C            FIN DU TRAITEMENT DES SOUS-TRIANGLES DE L ELEMENT
C            ****************************************************
         ENDIF
 1000 CONTINUE
C
C     LES LIGNES SONT REMISES A LEUR EPAISSEUR NORMALE
c     ------------------------------------------------
      CALL XVEPAISSEUR( 1 )

      RETURN
      END
