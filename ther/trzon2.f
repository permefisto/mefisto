      SUBROUTINE TRZON2( TMIN,   TMAX,   NDIM,
     %                   NCAS,   NCAS0,  NCAS1, NTDL,
     %                   NTYP,   TEMPER, dptemp,
     %                   NUTYEL, NBELEM,
     %                   NBNOEL, NUNDEL, NBPOEL, NUPTEL,
     %                   NBCOOR, NBNOEU, XYZNOE, NBPOIT, XYZPOI,
     %                   MXSOMM, SOLEL , COPOE ,
     %                   MXPIL3, NPILE3,
     %                   VALST , FBASE , XYZSOM )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACER LA TEMPERATURE PAR ZONES DE COULEURS SUR LES EF 2D
C -----    TRIANGLE OU QUADRANGLE
C
C ENTREES :
C ---------
C TMIN   : TEMPERATURE MINIMALE DES NCAS0:NCAS1 CAS
C TMAX   : TEMPERATURE MAXIMALE DES NCAS0:NCAS1 CAS
C NDIM   : ESPACE DE TRAVAIL 2 OU 3
C NCAS   : NUMERO DU CAS A TRAITER
C NCAS0  : NUMERO DU PREMIER VECTEUR TEMPERATURE
C NCAS1  : NUMERO DU PREMIER VECTEUR TEMPERATURE
C NTDL   : NOMBRE TOTAL DE DEGRES DE LIBERTE EN TEMPERATURE
C NTYP   : =0 EMPLOI de TEMPER
C          =1 EMPLOI de dptemp
C TEMPER : TABLEAU DES NTDL * NCAS0:NCAS1 TEMPERATURES
C dptemp : TABLEAU DES NCAS0:NCAS1 TABLEAU(NTDL) TEMPERATURES
C NUTYEL : NUMERO DU TYPE D'EF A TRAITER ICI
C NBELEM : NOMBRE D'ELEMENTS DE CE TYPE
C NBNOEL : NOMBRE DE  NOEUDS DE L'ELEMENT
C NUNDEL : NUMERO DES NOEUDS DES  ELEMENTS
C NBPOEL : NOMBRE DE  POINTS DE L'ELEMENT
C NUPTEL : NUMERO DES POINTS DES  ELEMENTS
C NBCOOR : NOMBRE DE COORDONNEES D'UN POINT OU NOEUD (3 ou 6)
C NBNOEU : NOMBRE TOTAL DE NOEUDS DU MAILLAGE
C XYZNOE : COORDONNEES DES NOEUDS DU MAILLAGE
C NBPOIT : NOMBRE TOTAL DE POINTS DU MAILLAGE
C XYZPOI : COORDONNEES DES POINTS DU MAILLAGE
C
C TABLEAUX AUXILIAIRES :
C ----------------------
C MXSOMM : NBRE MAXIMUM DE SOMMETS DES SOUS-TRIANGLES GENERES
C SOLEL  : TABLEAU AUXILIAIRE
C COPOE  : TABLEAU AUXILIAIRE
C MXPIL3 : NBRE MAXIMUM DE SOUS-TRIANGLES DANS LA PILE NPILE3
C NPILE3 : PILE DES 3 SOMMETS DES SOUS-TRIANGLES
C VALST  : VALST(1 ET 2,.) 2 COORDONNEES DES SOMMETS DES SOUS TRIANGLES
C          VALST(3,.) VALEUR DE LA SOLUTION EN CE SOMMET
C FBASE  : TABLEAU AUXILIAIRE (CF SP INTERP)
C XYZSOM : TABLEAU AUXILIAIRE (COORDONNEES DES SOMMETS DES SS-TRIANGLES)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 1990
C          PASCAL JOLY     ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 1990
C          ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 1994
C23456---------------------------------------------------------------012
      PARAMETER     ( LIGCON=0 )
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/donele.inc"
      include"./incl/inteel.inc"
      include"./incl/ponoel.inc"
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      COMMON / EPSSSS / EPZERO, EPSXYZ
      INTEGER           NUNDEL(NBELEM, NBNOEL), NUPTEL(NBELEM, NBPOEL),
     %                  NPILE3(3,MXPIL3)
      REAL              XYZNOE(NBCOOR,NBNOEU), XYZPOI(NBCOOR,NBPOIT)
      REAL              XYZSOM(3,MXSOMM)
      type typ_dptab
         DOUBLE PRECISION, dimension(:), pointer :: dptab
      end type typ_dptab
      type( typ_dptab ),dimension(NCAS0:NCAS1) :: dptemp
      DOUBLE PRECISION  TEMPER(NTDL,NCAS0:NCAS1), SOLEL(*),
     %                  COPOE(NBPOEL,NDIM)
      DOUBLE PRECISION  XD, YD, ZD, PROSCD, FBASE(30), XREF(3)
      REAL              VALST(3,MXSOMM)
      REAL              COSOTR(2,3),  C(3),  CS(2,3),  CV(3)
      EQUIVALENCE      (XD,XREF(1)),(YD,XREF(2)),(ZD,XREF(3))
      DATA              ZD  /0.D0/
C
C     CALCUL DU NOMBRE DE SUBDIVISIONS DE (0,1)
C     EN FONCTION DU NOMBRE DE NOEUDS DE L ELEMENT FINI
      IF( NBNOE .GT. 4 ) THEN
         NBSOUI = 2
      ELSE
         NBSOUI = 1
      ENDIF
C
C     RETROUVER LE NUMERO DES INTERPOLATIONS DE LA TEMPERATURE
C     ET DES COMPOSANTES DE LA TRANSFORMATION:ELT REFERENCE->ELEMENT
      CALL ELINTE('THERMIQUE',NUTYEL,NDIMF,NOINTF,
     &             NBINVA,NUINVA,NUINTI,NBNDIN)
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
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'INCOMPATIBILITE DIMENSION ESPACES pour F et X'
         ELSE
            KERR(1) = 'INCOMPATIBILITY DIMENSION SPACES for F and X'
         ENDIF
         CALL LEREUR
         RETURN
      ENDIF
C
      IF( NARET .LT. 3  .OR.  NARET .GT. 4 ) THEN
C        ELEMENT NON TRIANGULAIRE ET NON QUADRANGULAIRE
         NBLGRC(NRERR) = 2
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'ERREUR: ELEMENT FINI 2D'
            KERR(2) = 'NON TRIANGLE NON QUADRANGLE'
         ELSE
            KERR(1) = 'ERROR: 2D FINITE ELEMENT'
            KERR(2) = 'NOT a TRIANGLE NOT a QUADRANGLE'
         ENDIF
         CALL LEREUR
         RETURN
      ENDIF
C
C     PRECISION DES TESTS: TOLERANCE
      C(3) = 0.D0

      TMINMAX = TMAX - TMIN
C
C     LE TRACE DES ISOVALEURS
C     -----------------------
C     LIGNES EPAISSIES ET CONTINUES
      CALL XVEPAISSEUR( 1 )
      CALL XVTYPETRAIT( LIGCON )
C     LE NOMBRE DE COULEURS DISPONIBLES DANS LA PALETTE
      NBCOUL = NDCOUL - N1COUL
C
      DO 1000 K=1,NBELEM
C
C        OUVERTURE DE LA PILE DU NO DES 3 SOMMETS DES SOUS-TRIANGLES
C        -----------------------------------------------------------
C        LHSPI3 : POINTE SUR LE SOMMET DE LA PILE NPILE3
C        NOSOMM : NO DU DERNIER SOMMET GENERE DANS LES SOUS-TRIANGLES
C
         LHSPI3 = 0
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
   20    CALL EMPITR( NBSOUI,NOSOMM,MXSOMM,MXPIL3,LHSPI3,NPILE3 )
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
C        EXTRACTION DE LA TEMPERATURE NCAS DE TEMPER AUX NOEUDS DE L'EF
         DO L=1,NBNOEL
            NDL = NUNDEL(K,L)
            IF( NTYP .EQ. 0 ) THEN
               SOLEL( L ) = TEMPER( NDL, NCAS )
            ELSE
               SOLEL( L ) = dptemp( NCAS )%dptab( NDL )
            ENDIF
         ENDDO
C        EXTRACTION DES COORDONNEES DES POINTS DE L'EF
         DO 28 I=1,NBPOEL
            DO 27 J=1,NDIM
               COPOE(I,J) = XYZPOI( J, NUPTEL(K,I) )
 27         ENDDO
 28      ENDDO
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
            CALL INTERP( NOINTE, XD,YD,ZD, NBN,FBASE )
C
C           LA VALEUR DE LA SOLUTION EN CE POINT
            VALST(3,I) = REAL( PROSCD( FBASE, SOLEL, NBN ) )
C
C           LES COORDONNEES DES SOMMETS DES SOUS-TRIANGLES
            XREF(1) = VALST(1,I)
            XREF(2) = VALST(2,I)
            CALL FPOXYZ( NDIM,NBPOEL,NOINTF,COPOE,XREF,FBASE,C )
            XYZSOM(1,I) = C(1)
            XYZSOM(2,I) = C(2)
   78    ENDDO
C
C        **************************************
C        TANT QUE LA PILE3 N EST PAS VIDE FAIRE
C        **************************************
   80    IF( LHSPI3 .GT. 0 ) THEN
C
C           LES 3 SOMMETS ET LE TRIANGLE SONT DEPILES
C           -----------------------------------------
            DO 90 I=1,3
               NOSI    = NPILE3(I,LHSPI3)
               CS(1,I) = XYZSOM(1,NOSI)
               CS(2,I) = XYZSOM(2,NOSI)
               IF( TMINMAX .NE. 0 ) THEN
                  COULEUR = (VALST(3,NOSI)-TMIN) / TMINMAX
               ELSE
                  COULEUR = 1.0
               ENDIF
C              PROJECTION SI DEPASSEMENT AU DELA DU MIN OU MAX
               IF( COULEUR .GT. 1.0 ) COULEUR = 1.0
               IF( COULEUR .LT. 0.0 ) COULEUR = 0.0
               CV(I) = N1COUL + NBCOUL * COULEUR
   90       ENDDO
            LHSPI3 = LHSPI3 - 1
C
C           TRACE DU TRIANGLE SELON LA COULEUR AUX 3 SOMMETS
            CALL TRIACOUL2D( CS, CV )
C
C           REMONTEE POUR TRAITER LE SOUS-TRIANGLE SUIVANT
            GOTO 80
         ENDIF
 1000 ENDDO
C
      IF( NCOUAF .GE. 0 ) THEN
C
C        LE TRACE DES ARETES
C        -------------------
C        LIGNES NON EPAISSIES ET POINTILLEES
         CALL XVEPAISSEUR( 1 )
         CALL XVTYPETRAIT( NTLAFR )
C
         DO 8 K=1,NBELEM
            DO 5 I=1,NARET
C              LE NUMERO DU 1-ER SOMMET DE L'ARETE
               NS0 = NUNDEL( K, NOSOAR(1,I) )
C              LE NUMERO DU 2-EME SOMMET DE L'ARETE
               NS  = NUNDEL( K, NOSOAR(2,I) )
               IF( NBNOAR(I) .LE. 0 ) THEN
C                 TRACE DE L'ARETE
                  CALL TRAIT2D( NCOUAF, XYZNOE(1,NS0), XYZNOE(2,NS0),
     %                                  XYZNOE(1,NS ), XYZNOE(2,NS ) )
               ELSE
C                 TRACE DU SOMMET1 AU MILIEU
                  NM = NUNDEL( K, NONOAR(1,I) )
                  CALL TRAIT2D( NCOUAF, XYZNOE(1,NS0), XYZNOE(2,NS0),
     %                                  XYZNOE(1,NM ), XYZNOE(2,NM ) )
C                 TRACE DU MILIEU AU SOMMET 2
                  CALL TRAIT2D( NCOUAF, XYZNOE(1,NM), XYZNOE(2,NM),
     %                                  XYZNOE(1,NS), XYZNOE(2,NS) )
               ENDIF
 5          ENDDO
 8       ENDDO
C
      ENDIF
C
C     LES LIGNES SONT REMISES A LEUR EPAISSEUR NORMALE
C     ------------------------------------------------
      CALL XVEPAISSEUR( 1 )
      CALL XVTYPETRAIT( LIGCON )
      RETURN
      END
