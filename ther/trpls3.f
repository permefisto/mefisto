      SUBROUTINE TRPLS3( MODSEC, NOAXE,  POINTS, NBPLAN, COPLAN,
     %                   NBCOOR, NCAS0,  NCAS1,  NTDL,
     %                   NTYP,   TEMPER, dptemp,
     %                   NUTYEL, NBELEM,
     %                   NBNOEL, NUNDEL, NBPOEL, NUPTEL,
     %                   NBPOIT, XYZPOI,
     %                   MXSOMM, SOLEL , COPOE ,
     %                   MXPILE, LAPILE,
     %                   XYZRTC, FBASE,
     %                   NBCOPS, TRIANG, MXFPLA, NBFPLA,
     %                   XYZSFP, TEMSFP )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER EN 3D LES FACES INTERSECTION AVEC UN PLAN
C -----    A X ou Y ou Z=CONSTANTE ou QUELCONQUE
C          AVEC DES ELEMENTS FINIS D'UN TYPE DONNE
C          SOIT TETRAEDRE SOIT PENTAEDRE SOIT HEXAEDRE
C
C ENTREES :
C ---------
C NOAXE  : NUMERO (1 A 3) DE L'AXE OU LA COORDONNEE EST CONSTANTE
C          4 DEFINITION DU PLAN PAR UN POINT ET UN VECTEUR OU AUTRE
C POINTS : 1. UN POINT DU PLAN + 2. LE VECTEUR NORMAL
C TMIN   : TEMPERATURE MINIMALE DES CAS NCAS0 a NCAS1
C TMAX   : TEMPERATURE MAXIMALE DES CAS NCAS0 a NCAS1
C NBPLAN : NOMBRE DE PLANS DE SECTION
C COPLAN : TABLEAU DES VALEURS DE LA COORDONNEE DES NBPLAN PLANS
C NBCOOR : NOMBRE DE COORDONNEES DES NOEUDS et POINTS  ici 3 OU 6
C NCAS0  : NUMERO DU PREMIER CAS A TRAITER
C NCAS1  : NUMERO DU DERNIER CAS A TRAITER
C NTDL   : NOMBRE TOTAL DE DEGRES DE LIBERTE EN TEMPERATURE
C NTYP   : =0 EMPLOI de TEMPER
C         =1  EMPLOI de dptemp
C TEMPER : TABLEAU DES NTDL * NCAS0:NCAS1 TEMPERATURES
C dptemp : TABLEAU DES NCAS0:NCAS1 TABLEAU(NTDL) TEMPERATURES
C
C NUTYEL : NUMERO DU TYPE D'EF A TRAITER ICI
C NBELEM : NOMBRE D'ELEMENTS FINIS DE CE TYPE
C NBNOEL : NOMBRE DE  NOEUDS DE L'ELEMENT  FINI  DE CE TYPE
C NUNDEL : NUMERO DES NOEUDS DES  ELEMENTS FINIS DE CE TYPE
C NBPOEL : NOMBRE DE  POINTS DE L'ELEMENT  FINI  DE CE TYPE
C NUPTEL : NUMERO DES POINTS DES  ELEMENTS FINIS DE CE TYPE
C
C NBNOEU : NOMBRE TOTAL DE NOEUDS DU MAILLAGE
C XYZNOE : COORDONNEES DES NOEUDS DU MAILLAGE
C NBPOIT : NOMBRE TOTAL DE POINTS DU MAILLAGE
C XYZPOI : COORDONNEES DES POINTS DU MAILLAGE
C NBCOPS : 3 COORDONNEES PAR SOMMET
C TRIANG : REEL TEMOIN DE TRIANGLE DANS XYZSFP(3,4,.)
C NBFPLA : NOMBRE DE FACES INITIALES DES PLANS DE SECTION
C MXFPLA : NOMBRE MAXIMAL DE FACES DECLARABLES DANS XYZSFP
C POINTS : POINTS(1:3,1) COORD DU POINT DU PLAN POUR OPTION NOAXE = 4
C          POINTS(1:3,2) VECTEUR NORMAL AU PLAN POUR OPTION NOAXE = 4
C
C TABLEAUX AUXILIAIRES :
C ----------------------
C MXSOMM : NOMBRE MAXIMUM DE SOMMETS DES SOUS-TETRAEDRES DE L'EF REFERENCE
C SOLEL  : TABLEAU AUXILIAIRE
C COPOE  : TABLEAU AUXILIAIRE
C MXPILE : NBRE MAXIMUM DE SOUS-TETRAEDRES DANS LA PILE LAPILE
C LAPILE : PILE DES 4 SOMMETS DES SOUS-TETRAEDRES
C XYZRTC : XYZRTC(1:3,*) 3 COORDONNEES DES SOMMETS DES SOUS-TETRAEDRES
C                        DANS L'EF DE REFERENCE
C          XYZRTC(4,*) VALEUR DE LA COORDONNEE NOAXE EN CE MEME SOMMET
C                      SUR L'EF COURANT
C FBASE  : TABLEAU AUXILIAIRE (CF SP INTERP)
C
C SORTIES:
C --------
C NBFPLA : NOMBRE FINAL DE FACES DES PLANS DE SECTION
C XYZSFP : XYZ DES 4 SOMMETS DES FACES DES PLANS DE SECTION
C TEMSFP : TEMPERATURE DES 4 SOMMETS AU PLUS DES FACES DES PLANS
C          POUR LES CAS NCAS0 A NCAS1
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 1994
C MODIFS : ALAIN PERRONNET LJLL UPMC & ST PIERRE DU PERRAY  OCTOBRE 2010
C23456---------------------------------------------------------------012
C     SEUIL DE PRECISION
      PARAMETER    ( EPS=0.00001,UNPEPS=1.0+EPS)
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/donele.inc"
      include"./incl/inteel.inc"
      include"./incl/ponoel.inc"
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)

      DOUBLE PRECISION  TEMPER(NTDL,NCAS0:NCAS1)
      type typ_dptab
         DOUBLE PRECISION, dimension(:), pointer :: dptab
      end type typ_dptab
      type( typ_dptab ),dimension(NCAS0:NCAS1) :: dptemp

      INTEGER           NUNDEL(NBELEM,NBNOEL),
     %                  NUPTEL(NBELEM,NBPOEL),
     %                  LAPILE(4,MXPILE)

      REAL              XYZPOI(NBCOOR,NBPOIT),
     %                  XYZSFP(NBCOPS,4,MXFPLA),
     %                  TEMSFP(4,NCAS0:NCAS1,MXFPLA),
     %                  COPLAN(1:NBPLAN),
     %                  XYZRTC(4,MXSOMM),
     %                  COSOTE(3,4),
     %                  CEF(4)
C     POUR LE CAS SPECIAL DE PLAN DEFINI PAR UTILISATEUR
      REAL              POINTS(3,2)

      DOUBLE PRECISION  SOLEL(NBNOEL),
     %                  COPOE(NBPOEL,NBCOOR),
     %                  XD, YD, ZD, PROSCD,
     %                  FBASE(30),
     %                  XREF(3,12)

C     DONNEES DES SOUS-TETRAEDRES D'UN TETRAEDRE ou PYRAMIDE ou PENTAEDRE
C     ou HEXAEDRE SUBDIVISE REGULIEREMENT
      include "./incl/nostst.inc"
C
C     CALCUL DU NOMBRE DE SUBDIVISIONS DE (0,1)
C     EN FONCTION DU NOMBRE DE NOEUDS DE L ELEMENT FINI
      IF( NBNOE .LE. 8 ) THEN
         NBSOUI = 1
      ELSE
         NBSOUI = 2
      ENDIF
cccC     TETRAEDRE DE BREZZI-FORTIN
ccc      IF( NUTYEL .EQ. 19 ) NBSOUI = 3 NON PROGRAMME
cccc     TOUS LES SOMMETS SONT SUR LES FACES => INTERPOLATION P1
      NBSTTE = 2 * ( NBSOUI * NBSOUI + 1 )
C
C     RETROUVER LE NUMERO DES INTERPOLATIONS DE LA TEMPERATURE
C     ET DES COMPOSANTES DE LA TRANSFORMATION:ELT REFERENCE->ELEMENT
      CALL ELINTE( 'THERMIQUE', NUTYEL, NDIMF, NOINTF,
     &              NBINVA, NUINVA, NUINTI, NBNDIN )
C     LA TEMPERATURE EST ICI LA SEULE INCONNUE VARIATIONNELLE
      NOINTE = NUINTI(1)
cccC     TETRAEDRE DE BREZZI-FORTIN
ccc      IF( NUTYEL .EQ. 19 ) NOINTE = 3003
cccc     TOUS LES SOMMETS SONT SUR LES FACES => INTERPOLATION P1
C
CCC      IDEM POUR LE GRADIENT DE TEMPERATURE
CCC      CALL ELINT1( 'THERMIQUE', NUTYEL, NDIMF, NOINFD, NBPOF,
CCC     &              NBINVA, NUINVA, NUINTI, NBNDIN, NUNOIN )
C
      IF( NFACE .LT. 4  .OR.  NFACE .GT. 6 ) THEN
C        ELEMENT NON TETRAEDRE OU PENTAEDRE OU HEXAEDRE
         NBLGRC(NRERR) = 3
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'ERREUR: ELEMENT FINI DIFFERENT DE'
            KERR(2) = 'TETRAEDRE ou PENTAEDRE ou HEXAEDRE'
            KERR(3) = 'PLAN-VALEURS NON CALCULABLES'
         ELSE
            KERR(1) = 'ERROR: FINITE ELEMENT DIFFERENT of'
            KERR(2) = 'TETRAHEDRON or PENTAHEDRON or HEXAHEDRON'
            KERR(3) = 'PLANE-VALUES NOT COMPUTABLE'
         ENDIF
         CALL LEREUR
         RETURN
      ENDIF
C
C     FORMATION DES TABLEAUX SUR L'ELEMENT FINI DE REFERENCE
C     ======================================================
C     OUVERTURE DE LA PILE DU NO DES 4 SOMMETS DES SOUS-TETRAEDRES
C     LHPIL0 POINTE SUR LE SOMMET DE LA PILE LAPILE
      LHPIL0 = 0
C     NOSOMM NO DU DERNIER SOMMET GENERE DANS LES SOUS-TETRAEDRE
      NOSOMM = 0
C
C     GENERATION DES SOUS-TETRAEDRES A L'INTERIEUR DE L'EF DE REFERENCE
C     =================================================================
C     LE NO DES 4 SOMMETS DE CHAQUE TETRAEDRE
C     LES COORDONNEES DES SOMMETS DES SOUS-TETRAEDRES DANS L'EF REFERENCE
C     SONT EMPILES DANS LAPILE(4,MXPILE) ET XYZRTC(1:4,MXSOMM)
C
      IF( NFACE .EQ. 4 ) THEN
C
C        *************
C        * TETRAEDRE *
C        *************
C
C        LE NO DES 4 SOMMETS DE CHAQUE SOUS-TETRAEDRE
C        --------------------------------------------
         CALL EMPITE( NBSOUI,NOSOMM,MXSOMM,MXPILE,LHPIL0,LAPILE )
         IF( LHPIL0 .LT. 0 ) GOTO 9999
C
C        LES 3 COORDONNEES DES 2*(NBSOUI**2+1) SOMMETS DES SOUS-TETRAEDRES
C        -----------------------------------------------------------------
         CALL COPTTE( TETR, NBSOUI+1, 4, XYZRTC(1,NOSOMM+1) )
C
         NOSOMM = NOSOMM + NBSTTE
C
C
      ELSE IF( NFACE .EQ. 5 ) THEN
C
         IF( NBNSOM .LE. 5 ) THEN
C
C           *************
C           * PYRAMIDE  *
C           *************
C           BOUCLE SUR 2 SOUS-TETRAEDRES PRIMAIRES
C           ======================================
            DO 14 I = 1, 2
C
C              LE NO DES 4 SOMMETS DE CHAQUE SOUS-TETRAEDRE
C              --------------------------------------------
               CALL EMPITE( NBSOUI,NOSOMM,MXSOMM,MXPILE,LHPIL0,LAPILE )
               IF( LHPIL0 .LT. 0 ) GOTO 9999
C
               DO 12 J = 1, 4
                  NO = NOSY(J, I)
                  COSOTE(1, J) = PYRA(1, NO)
                  COSOTE(2, J) = PYRA(2, NO)
                  COSOTE(3, J) = PYRA(3, NO)
 12            CONTINUE
C
C              LES COORDONNEES DES 2*(NBSOUI**2+1) SOMMETS DES SOUS-TETRAEDRES
C              ---------------------------------------------------------------
               CALL COPTTE( COSOTE, NBSOUI+1, 4, XYZRTC(1,NOSOMM+1) )
C
               NOSOMM = NOSOMM + NBSTTE
 14         CONTINUE
C
         ELSE
C
C           *************
C           * PENTAEDRE *
C           *************
C           BOUCLE SUR 3 SOUS-TETRAEDRES PRIMAIRES
C           ======================================
            DO 18 I = 1, 3
C
C              LE NO DES 4 SOMMETS DE CHAQUE SOUS-TETRAEDRE
C              --------------------------------------------
               CALL EMPITE( NBSOUI,NOSOMM,MXSOMM,MXPILE,LHPIL0,LAPILE )
               IF( LHPIL0 .LT. 0 ) GOTO 9999
C
               DO 16 J = 1, 4
                  NO = NOSP(J, I)
                  COSOTE(1, J) = PENT(1, NO)
                  COSOTE(2, J) = PENT(2, NO)
                  COSOTE(3, J) = PENT(3, NO)
 16            CONTINUE
C
C              LES COORDONNEES DES 2*(NBSOUI**2+1) SOMMETS DES SOUS-TETRAEDRES
C              ---------------------------------------------------------------
               CALL COPTTE( COSOTE, NBSOUI+1, 4, XYZRTC(1,NOSOMM+1) )
C
               NOSOMM = NOSOMM + NBSTTE
 18         CONTINUE
         ENDIF
C
      ELSE
C
C        ************
C        * HEXAEDRE *
C        ************
C
C        BOUCLE SUR 5 SOUS-TETRAEDRES PRIMAIRES
C        ======================================
         DO 26 I = 1, 5
C
C           LE NO DES 4 SOMMETS DE CHAQUE SOUS-TETRAEDRE
C           --------------------------------------------
            CALL EMPITE( NBSOUI,NOSOMM,MXSOMM,MXPILE,LHPIL0,LAPILE )
            IF( LHPIL0 .LT. 0 ) GOTO 9999
C
            DO 24 J = 1, 4
               NO = NOSH(J, I)
               COSOTE(1, J) = HEXA(1, NO)
               COSOTE(2, J) = HEXA(2, NO)
               COSOTE(3, J) = HEXA(3, NO)
 24         CONTINUE
C
C           LES 3 COORDONNEES DES 2*(NBSOUI**2+1) SOMMETS DES SOUS-TETRAEDRES
C           -----------------------------------------------------------------
            CALL COPTTE( COSOTE, NBSOUI+1, 4, XYZRTC(1,NOSOMM+1) )
C
            NOSOMM = NOSOMM + NBSTTE
 26      CONTINUE
      ENDIF
C
C
C     LA BOUCLE SUR LES EF DE CE TYPE
C     ===============================
      DO 1000 NEF=1,NBELEM
C
C        COPOE(NBPOEL,NBCOOR) COORDONNEES DES POINTS DE L'EF COURANT
         DO 38 I=1,NBPOEL
            DO 37 J=1,NBCOOR
               COPOE(I,J) = XYZPOI( J, NUPTEL(NEF,I) )
 37         CONTINUE
 38      CONTINUE
C
C        BOUCLE SUR LES SOMMETS DES SOUS-TETRAEDRES
         DO 78 I=1,NOSOMM
C
C           LES 3 COORDONNEES DU POINT I DANS L'EF DE REFERENCE
            XD = XYZRTC(1,I)
            YD = XYZRTC(2,I)
            ZD = XYZRTC(3,I)
C
C           LA VALEUR DES NBN FONCTIONS DE BASE EN (XD,YD,ZD)
            CALL INTERP( NOINTE, XD, YD, ZD,  NBN, FBASE )
C
C           LA COORDONNEE NOAXE EN LES POINTS DE L'EF COURANT
            IF (NOAXE .EQ. 4) THEN
C
C              PLAN QUELCONQUE
C              CALCUL DE LA COORDONNEE CURVILIGNE A PARTIR DU POINT DU PLAN
C              DANS LA DIRECTION NORMALE
C              PRODUIT SCALAIRE( (SOMMET-PTPLAN) . VECTEUR NORMAL PLAN )
C              PASCAL HAVE DEA ANALYSE NUMERIQUE UPMC PARIS JANVIER 2000
               XYZRTC(4,I) = REAL( 
     %               (PROSCD( FBASE, COPOE(1,1), NBN ) - POINTS(1,1))
     %              * POINTS(1,2)
     %              +(PROSCD( FBASE, COPOE(1,2), NBN ) - POINTS(2,1))
     %              * POINTS(2,2)
     %              +(PROSCD( FBASE, COPOE(1,3), NBN ) - POINTS(3,1))
     %              * POINTS(3,2) )
C
            ELSE
C
C              PLAN A NOAXE CONSTANT
               XYZRTC(4,I) = REAL( PROSCD( FBASE, COPOE(1,NOAXE), NBN ))
C
            ENDIF
 78      CONTINUE
C
C        LE SOMMET DE LA PILE AU DEBUT DE CET EF
         LHPILE = LHPIL0
C
C        ***********************************************************
C        TANT QUE  LA PILE DES SOUS-TETRAEDRES N'ET PAS VIDE   FAIRE
C        ***********************************************************
 80      IF( LHPILE .LE. 0 ) GOTO 1000
C
C        LES 4 SOMMETS DU TETRAEDRE EN HAUT DE LA PILE SONT DEPILES
         DO 200 I=1,4
            K = LAPILE(I,LHPILE)
C           3 COORDONNEES DU SOMMET I DU SOUS TETRAEDRE DANS
C           L'EF DE REFERENCE
            COSOTE(1,I) = XYZRTC(1,K)
            COSOTE(2,I) = XYZRTC(2,K)
            COSOTE(3,I) = XYZRTC(3,K)
C           LA COORDONNEE NOAXE AU SOMMET I DU SOUS-TETRAEDRE COURANT
            CEF( I )    = XYZRTC(4,K)
 200     CONTINUE
         LHPILE = LHPILE - 1
C
C        CALCUL DE LA COORDONNEE NOAXE MIN ET MAX
C        AUX 4 SOMMETS DU SOUS TETRAEDRE
         CMIN = CEF(1)
         CMAX = CEF(1)
         DO 205 K=2,4
            CK = CEF(K)
            IF( CK .GT. CMAX ) THEN
               CMAX = CK
            ELSE IF( CK .LT. CMIN ) THEN
               CMIN = CK
            ENDIF
 205     CONTINUE
C
C        TRAITEMENT DES CAS EXTREMES
         IF( COPLAN(1)      .GT. CMAX   .OR.
     %       COPLAN(NBPLAN) .LT. CMIN ) GOTO 80
C
C        BOUCLE SUR LES PLANS
         DO 300 K=1,NBPLAN
C
C           VALEUR DE LA COORDONNEE CONSTANTE DU PLAN K
            COOPLK = COPLAN( K )
C
C           LE PLAN K EST IL COMPRIS ENTRE CMIN ET CMAX ?
            IF( COOPLK .GE. CMIN .AND. COOPLK .LE. CMAX ) THEN
C
C              RECHERCHE DES POINTS D INTERSECTION DU PLAN K
C              AVEC LES 6 ARETES DU SOUS TETRAEDRE
C              ---------------------------------------------
C              NBINT NOMBRE D INTERSECTION DU PLAN AVEC LES 6 ARETES
               NBINT = 0
C
               DO 250 I=1,6
C
C                 L'ARETE I DU TETRAEDRE DE SOMMETS N1 N2
                  N1 = NOSTAR(1,I)
                  C1 = CEF( N1 )
C
                  N2 = NOSTAR(2,I)
                  C2 = CEF( N2 )
C
                  C21 = C2 - C1
                  IF( ABS( C21 ) .GT. EPS*ABS( C1 ) ) THEN
C
C                    C1 =/ C2
C                    --------
                     CK = ( COOPLK - C1 ) / C21
                     IF( CK .LE. -EPS .OR. CK .GE. UNPEPS ) GOTO 250
C
C                    UN POINT D'INTERSECTION
                     NBINT = NBINT + 1
                     CK0 = 0
                     CK1 = 1
C                    LES COORDONNEES DU POINT D'INTERSECTION
C                    DANS L'EF REFERENCE
                     NBDICHO = 0
 209                 DO 210 LL = 1,3
                        XREF(LL,NBINT) = COSOTE(LL,N1) * (1.0-CK) +
     %                                   COSOTE(LL,N2) * CK
 210                 CONTINUE
C
C                    SI FE N'EST PAS AFFINE IL FAUT TROUVER LE PT
C                    D'INTERSECTION VERIFIANT LA COORDONNEE
                     IF( NBPOEL .GT. 4 .AND. NOAXE .NE. 4 ) THEN
C
C                       LA VALEUR DES FONCTIONS DE BASE AU POINT XREF
                        CALL INTERP( NOINTF, XREF(1,NBINT),
     %                               XREF(2,NBINT), XREF(3,NBINT),
     %                               NBPOEL, FBASE )
C                       LA COORDONNEE NOAXE DU SOMMET DANS L'EF COURANT
                        XC = REAL(PROSCD(FBASE, COPOE(1,NOAXE), NBPOEL))
                        XC = ( XC - COOPLK ) / C21
                        IF( ABS(XC) .GT. 0.001 ) THEN
C                          ITERATION DE DICHOTOMIE
                           NBDICHO = NBDICHO + 1
                           IF( NBDICHO .LT. 64 ) THEN
                              IF( XC .GE. 0 ) THEN
                                 CK1 = CK
                              ELSE
                                 CK0 = CK
                              ENDIF
                              CK = ( CK0 + CK1 ) * 0.5
                              GOTO 209
                           ENDIF
                        ENDIF
                     ENDIF
C
                  ELSE
C
C                    C1 = C2 SONT ILS EGAUX A COOPLK ?
C                    ---------------------------------
                     IF( ABS(COOPLK-C1) .LE. ABS(CMAX-CMIN)*EPS ) THEN
C                       L'ARETE EST DANS LE PLAN K
                        NBINT = NBINT + 2
                        DO 220 LL = 1,3
                           XREF(LL,NBINT-1) = COSOTE(LL,N1)
                           XREF(LL,NBINT  ) = COSOTE(LL,N2)
 220                    CONTINUE
                     ENDIF
C
                  ENDIF
C                 FIN DE TRAITEMENT DE L'ARETE I
 250           CONTINUE
C
C              TRANSFORMATION DES 3 OU 4 SEGMENTS D'INTERSECTION
C              EN UN TRIANGLE OU QUADRANGLE = FACE DU PLAN
C              DANS L'EF DE REFERENCE
C              -------------------------------------------------
               CALL SETRQU( NBINT, XREF )
C
               IF( NBINT .NE. 3 .AND. NBINT .NE. 4 ) GOTO 300
C
C              UNE FACE DE PLUS DANS LE PLAN K
               IF( NBFPLA .GE. MXFPLA ) THEN
                  WRITE(IMPRIM,*) 'MXFPLA=',MXFPLA,
     %                            ' INSUFFISANT DANS trpls3.f'
                  GOTO 9999
               ENDIF
               NBFPLA = NBFPLA + 1
C
C              VALEUR TEMOIN DE FACE=TRIANGLE
C              CETTE VALEUR SERA ECRASEE SI CETTE FACE EST UN QUADRANGLE
               XYZSFP(3,4,NBFPLA) = TRIANG
C
C              COORDONNEES DES NBINT SOMMETS DE LA FACE PLAN DANS L'OBJET
               DO 295 I = 1,NBINT
C
C                 LA VALEUR DES FONCTIONS DE BASE AU POINT I
C                 ..........................................
                  CALL INTERP( NOINTF, XREF(1,I), XREF(2,I), XREF(3,I),
     %                         NBPOEL, FBASE )
C
C                 LES 3 COORDONNEES DU SOMMET I DE LA FACE NBFPLA
C                 ...............................................
                  DO L=1,3
                     XYZSFP(L,I,NBFPLA)=
     %                      REAL( PROSCD(FBASE,COPOE(1,L),NBPOEL) )
                  ENDDO
C                 PLUS LA COORDONNEE CTE DU PLAN K POUR LES PROFILS
                  IF( MODSEC .EQ. 1 ) XYZSFP(4,I,NBFPLA) = COOPLK
c
ccc               print 10292,I,NBFPLA,(XYZSFP(L,I,NBFPLA),L=1,3)
ccc10292          FORMAT('trpls3: XYZSFP(',I1,',',I4,')=',3G14.6)
c
C                 LA SOLUTION AU SOMMET I DE LA FACE POUR LES NCAS0 A NCAS1
C                 .........................................................
C                 ICI POINTS=NOEUDS EST SUPPOSE
                  DO NCAS = NCAS0, NCAS1
C                    EXTRACTION DE LA TEMPERATURE NCAS DE TEMPER
                     DO L=1,NBNOEL
                        NDL = NUNDEL(NEF,L)
                        IF( NTYP .EQ. 0 ) THEN
                           SOLEL( L ) = TEMPER( NDL, NCAS )
                        ELSE
                          SOLEL( L ) = dptemp( NCAS )%dptab( NDL )
                       ENDIF
                     ENDDO

ccc                  IF( NUTYEL .EQ. 19 ) THEN
cccC                    TETRAEDRE DE BREZZI-FORTIN VALEUR AU BARYCENTRE
cccc                    TOUS LES SOMMETS SONT SUR LES FACES => INTERPOLATION P1
ccc                     SOLEL(5) = TEMPER( NBELEM+NEF, NCAS ) INUTILE
ccc                  ENDIF

                     TEMSFP(I,NCAS,NBFPLA)=
     %                      REAL( PROSCD(FBASE, SOLEL, NBPOEL) )
                  ENDDO
 295           CONTINUE
C
            ENDIF
C
 300     CONTINUE
C
C        REMONTEE POUR TRAITER LE SOUS-TETRAEDRE SUIVANT
C        ===============================================
         GOTO 80
C
C        *************************************************
C        FIN DU TRAITEMENT DES SOUS-TETRAEDRES DE L'EF NEF
C        *************************************************
 1000 CONTINUE
C
 9999 RETURN
      END
