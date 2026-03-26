      SUBROUTINE TRISO3( NBISO , VALISO,
     %                   NDIM  , NCAS,   NCAS0,  NCAS1, NTDL,
     %                   NTYP  , TEMPER, dptemp,
     %                   NUTYEL, NBELEM,
     %                   NBNOEL, NUNDEL, NBPOEL, NUPTEL,
     %                   NBCOOR, NBPOIT, XYZPOI,
     %                   MXSOMM, SOLEL , COPOE,
     %                   MXPILE, LAPILE,
     %                   VALST , FBASE,
     %                   TRIANG, MXFISO, NBFISO, XYZISO, NUMISO )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER EN 3D LES FACES ISOTHERMES DES ELEMENTS FINIS D'UN
C -----    TYPE DONNE SOIT TETRAEDRE SOIT PENTAEDRE SOIT HEXAEDRE

C ENTREES :
C ---------
C NBISO  : NOMBRE D'ISOVALEURS A CALCULER
C VALISO : TABLEAU DES VALEURS DES ISOVALEURS A CALCULER
C NDIM   : ESPACE DE TRAVAIL 3
C NCAS   : NUMERO DU CAS A TRAITER
C NCAS0  : NUMERO DU PREMIER VECTEUR TEMPERATURE
C NCAS1  : NUMERO DU PREMIER VECTEUR TEMPERATURE
C NTDL   : NOMBRE TOTAL DE DEGRES DE LIBERTE EN TEMPERATURE
C NTYP   : =0 EMPLOI de TEMPER
C          =1 EMPLOI de dptemp
C TEMPER : TABLEAU DES NTDL * NCAS0:NCAS1 TEMPERATURES
C dptemp : TABLEAU DES NCAS0:NCAS1 TABLEAU(NTDL) TEMPERATURES

C NUTYEL : NUMERO DU TYPE D'EF A TRAITER ICI
C NBELEM : NOMBRE D'ELEMENTS FINIS DE CE TYPE
C NBNOEL : NOMBRE DE  NOEUDS DE L'ELEMENT  FINI  DE CE TYPE
C NUNDEL : NUMERO DES NOEUDS DES  ELEMENTS FINIS DE CE TYPE
C NBPOEL : NOMBRE DE  POINTS DE L'ELEMENT  FINI  DE CE TYPE
C NUPTEL : NUMERO DES POINTS DES  ELEMENTS FINIS DE CE TYPE

C NBCOOR : NOMBRE DE COORDONNEES D'UN POINT (3 ou 6)
C NBNOEU : NOMBRE TOTAL DE NOEUDS DU MAILLAGE
C XYZNOE : COORDONNEES DES NOEUDS DU MAILLAGE
C NBPOIT : NOMBRE TOTAL DE POINTS DU MAILLAGE
C XYZPOI : COORDONNEES DES POINTS DU MAILLAGE
C TRIANG : REEL TEMOIN DE TRIANGLE DANS XYZISO(3,4,.)
C MXFISO : NOMBRE MAXIMAL DE FACES DECLARABLES DANS XYZISO

C TABLEAUX AUXILIAIRES :
C ----------------------
C MXSOMM : NOMBRE MAXIMUM DE SOMMETS DES SOUS-TETRAEDRES DE L'EF REFERENCE
C SOLEL  : TABLEAU AUXILIAIRE
C COPOE  : TABLEAU AUXILIAIRE
C MXPILE : NBRE MAXIMUM DE SOUS-TETRAEDRES DANS LA PILE LAPILE
C LAPILE : PILE DES 4 SOMMETS DES SOUS-TETRAEDRES
C VALST  : VALST(1:3,*) 3 COORDONNEES DES SOMMETS DES SOUS-TETRAEDRES
C          VALST(4,*) VALEUR DE LA SOLUTION EN CE MEME SOMMET
C FBASE  : TABLEAU AUXILIAIRE (CF SP INTERP)

C SORTIES:
C --------
C NBFISO : NOMBRE DE FACES DES ISOVALEURS
C XYZISO : XYZ DES 4 SOMMETS DES FACES DES ISOVALEURS
C NUMISO : NUMERO DE L'ISOVALEUR DE CHAQUE FACE ISOVALEUR
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     OCTOBRE 1994
C23456---------------------------------------------------------------012
C     SEUIL DE PRECISION
      PARAMETER    ( EPS=0.0001,UNPEPS=1.0+EPS)
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/donele.inc"
      include"./incl/inteel.inc"
      include"./incl/ponoel.inc"
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      INTEGER           NUNDEL(NBELEM,NBNOEL),
     %                  NUPTEL(NBELEM,NBPOEL),
     %                  LAPILE(4,MXPILE),
     %                  NUMISO(1:NBISO)
      REAL              XYZPOI(NBCOOR,NBPOIT),
     %                  XYZISO(3,4,MXFISO)
      type typ_dptab
         DOUBLE PRECISION, dimension(:), pointer :: dptab
      end type typ_dptab
      type( typ_dptab ),dimension(NCAS0:NCAS1) :: dptemp

      DOUBLE PRECISION  TEMPER(NTDL,NCAS0:NCAS1),
     %                  SOLEL(*),
     %                  COPOE(NBPOEL,NDIM)
      DOUBLE PRECISION  XD,YD,ZD,PROSCD,
     %                  FBASE(30),
     %                  XYZEFR(3,12)
      REAL              VALISO(1:NBISO),
     %                  VALST(4,MXSOMM),
     %                  COSOTE(3,4),
     %                  V(4)
C
C     DONNEES DES SOUS-TETRAEDRES D'UN TETRAEDRE ou PYRAMIDE ou PENTAEDRE
C     ou HEXAEDRE SUBDIVISE REGULIEREMENT
      include "./incl/nostst.inc"
C
C     NOMBRE DE COULEURS DISPONIBLES
      NBCOUL = NDCOUL - N1COUL
C
C     CALCUL DU NOMBRE DE SUBDIVISIONS DE (0,1)
C     EN FONCTION DU NOMBRE DE NOEUDS DE L ELEMENT FINI
      IF( NBNOE .LE. 8 ) THEN
         NBSOUI = 1
      ELSE
         NBSOUI = 2
      ENDIF
      NBSTTE = 2 * ( NBSOUI * NBSOUI + 1 )
C
C     RETROUVER LE NUMERO DES INTERPOLATIONS DE LA TEMPERATURE
C     ET DES COMPOSANTES DE LA TRANSFORMATION:ELT REFERENCE->ELEMENT
      CALL ELINTE( 'THERMIQUE', NUTYEL, NDIMF, NOINTF,
     &              NBINVA, NUINVA, NUINTI, NBNDIN )
C     LA TEMPERATURE EST ICI LA SEULE INCONNUE VARIATIONNELLE
      NOINTE = NUINTI(1)
C
CCC      IDEM POUR LE GRADIENT DE TEMPERATURE
CCC      CALL ELINT1( 'THER', NUTYEL, NDIMF, NOINFD, NBPOF,
CCC     &              NBINVA, NUINVA, NUINTI, NBNDIN, NUNOIN )
C
C     LA DIMENSION DE L'ESPACE DES COORDONNEES  EST NDIM
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
      IF( NFACE .LT. 4  .OR.  NFACE .GT. 6 ) THEN
C        ELEMENT NON TETRAEDRE OU PENTAEDRE OU HEXAEDRE OU PYRAMIDE
         NBLGRC(NRERR) = 3
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'ERREUR: ELEMENT FINI DIFFERENT DE'
            KERR(2) = 'TETRAEDRE ou PENTAEDRE ou HEXAEDRE ou PYRAMIDE'
            KERR(3) = 'ISOVALEURS NON CALCULABLES'
         ELSE
            KERR(1) = 'ERROR: FINITE ELEMENT DIFFERENT of'
           KERR(2)='TETRAHEDRON or PENTAHEDRON or HEXAHEDRON or PYRAMID'
            KERR(3) = 'ISOVALUES NOT COMPUTABLE'
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
C     SONT EMPILES DANS LAPILE(4,MXPILE) ET VALST(4,MXSOMM)
C
      IF( NFACE .EQ. 4 ) THEN
C
C        *************
C        * TETRAEDRE *
C        *************
C        LE NO DES 4 SOMMETS DE CHAQUE SOUS-TETRAEDRE
C        --------------------------------------------
         CALL EMPITE( NBSOUI,NOSOMM,MXSOMM,MXPILE,LHPIL0,LAPILE )
         IF( LHPIL0 .LT. 0 ) GOTO 9999
C
C        LES COORDONNEES DES 2*(NBSOUI**2+1) SOMMETS DES SOUS-TETRAEDRES
C        ---------------------------------------------------------------
         CALL COPTTE( TETR, NBSOUI+1, 4, VALST(1,NOSOMM+1) )
C
         NOSOMM = NOSOMM + NBSTTE
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
 12            ENDDO
C
C              LES COORDONNEES DES 2*(NBSOUI**2+1) SOMMETS DES SOUS-TETRAEDRES
C              ---------------------------------------------------------------
               CALL COPTTE( COSOTE, NBSOUI+1, 4, VALST(1,NOSOMM+1) )
C
               NOSOMM = NOSOMM + NBSTTE
 14         ENDDO
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
 16            ENDDO
C
C              LES COORDONNEES DES 2*(NBSOUI**2+1) SOMMETS DES SOUS-TETRAEDRES
C              ---------------------------------------------------------------
               CALL COPTTE( COSOTE, NBSOUI+1, 4, VALST(1,NOSOMM+1) )
C
               NOSOMM = NOSOMM + NBSTTE
 18         ENDDO
         ENDIF
C
      ELSE
C
C        ************
C        * HEXAEDRE *
C        ************
C        BOUCLE SUR 5 SOUS-TETRAEDRES PRIMAIRES
C        ======================================
         DO 25 I = 1, 5
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
 24         ENDDO
C
C           LES COORDONNEES DES 2*(NBSOUI**2+1) SOMMETS DES SOUS-TETRAEDRES
C           ---------------------------------------------------------------
            CALL COPTTE( COSOTE, NBSOUI+1, 4, VALST(1,NOSOMM+1) )
C
            NOSOMM = NOSOMM + NBSTTE
 25      ENDDO
      ENDIF
C
C     LA BOUCLE SUR LES EF DE CE TYPE
C     ===============================
      DO 1000 NEF=1,NBELEM
C
C        FORMATION DES TABLEAUX SUR L'EF COURANT
C        SOLEL(NBNOEL)      TEMPERATURE AUX NBNOEL DE L'EF COURANT
C        COPOE(NBPOEL,NDIM) COORDONNEES DES POINTS DE L'EF COURANT
C        =========================================================
C        EXTRACTION DE LA TEMPERATURE NCAS DE TEMPER
         DO L=1,NBNOEL
            NDL = NUNDEL(NEF,L)
            IF( NTYP .EQ. 0 ) THEN
               SOLEL( L ) = TEMPER( NDL, NCAS )
            ELSE
               SOLEL( L ) = dptemp( NCAS )%dptab( NDL )
            ENDIF
         ENDDO

C        EXTRACTION DES COORDONNEES DES POINTS DE L'ELEMENT FINI
         DO 38 I=1,NBPOEL
            DO 37 J=1,NDIM
               COPOE(I,J) = XYZPOI( J, NUPTEL(NEF,I) )
 37         ENDDO
 38      ENDDO
C
C        CALCUL DE LA SOLUTION AUX SOMMETS DES SOUS-TETRAEDRES
C        =====================================================
C        BOUCLE SUR LES SOMMETS DES SOUS-TETRAEDRES
C        ------------------------------------------
         DO 78 I=1,NOSOMM
C           LES 3 COORDONNEES DU POINT I DANS L'EF DE REFERENCE
            XD = VALST(1,I)
            YD = VALST(2,I)
            ZD = VALST(3,I)
C
C           LA VALEUR DES NBN FONCTIONS DE BASE EN (XD,YD,ZD)
            CALL INTERP( NOINTE, XD, YD, ZD,  NBN, FBASE )
C
C           LA VALEUR DE LA SOLUTION EN CE POINT DE L'EF
C          (MEME VALEUR DE LA SOLUTION SUR L'EF REFERENCE OU COURANT)
            VALST(4,I) = REAL( PROSCD( FBASE, SOLEL, NBN ) )
 78      ENDDO
C
C        LE SOMMET DE LA PILE AU DEBUT DE CET EF
         LHPILE = LHPIL0
C
C        ***********************************************************
C        TANT QUE  LA PILE DES SOUS-TETRAEDRES N'ET PAS VIDE   FAIRE
C        ***********************************************************
 80      IF( LHPILE .LE. 0 ) GOTO 1000
C
C        LES 4 SOMMETS DU TETRAEDRE EN SOMMET DE PILE SONT DEPILES
C        ---------------------------------------------------------
         DO 200 I=1,4
            J = LAPILE(I,LHPILE)
            COSOTE(1,I) = VALST(1,J)
            COSOTE(2,I) = VALST(2,J)
            COSOTE(3,I) = VALST(3,J)
            V( I )      = VALST(4,J)
 200     ENDDO
         LHPILE = LHPILE - 1
C
C        CALCUL DU MIN ET MAX AUX 4 SOMMETS DU TETRAEDRE
C        ===============================================
         VMIN = V(1)
         VMAX = V(1)
         DO 205 I=2,4
            VI = V(I)
            IF( VI .GT. VMAX ) THEN
               VMAX = VI
            ELSE IF( VI .LT. VMIN ) THEN
               VMIN = VI
            ENDIF
 205     ENDDO
C
C        TRAITEMENT DES CAS EXTREMES
         IF( VALISO(1)     .GT. VMAX   .OR.
     %       VALISO(NBISO) .LT. VMIN ) GOTO 80
C
C        RECHERCHE DE L INTERSECTION AVEC LES 6 ARETES DU TETRAEDRE
C        ==========================================================
         DO 300 K=1,NBISO
C
C           VALEUR DE L'ISOVALEUR K
            VALISK = VALISO( K )
C
C           L'ISOVALEUR K EST ELLE COMPRISE ENTRE VMIN ET VMAX ?
            IF( VALISK .GE. VMIN .AND. VALISK .LE. VMAX ) THEN
C
C              RECHERCHE DES POINTS D INTERSECTION DE L'ISOVALEUR AVEC LES ARETE
C              -----------------------------------------------------------------
C              NBINT NOMBRE D INTERSECTION DE L'ISOVALEUR AVEC LES 6 ARETES
               NBINT = 0
               DO 250 I=1,6
C
C                 L'ARETE I DU TETRAEDRE DE SOMMETS N1 N2
                  N1 = NOSTAR(1,I)
                  V1 = V( N1 )
C
                  N2 = NOSTAR(2,I)
                  V2 = V( N2 )
C
                  V21 = V2 - V1
                  IF( ABS( V21 ) .GT. EPS*ABS( V1 ) ) THEN
C
C                    V1 =/ V2
C                    --------
                     VI = ( VALISK - V1 ) / V21
                     IF( VI .GE. -EPS .AND. VI .LE. UNPEPS ) THEN
C
C                       UN POINT D'INTERSECTION
                        NBINT = NBINT + 1
                        VI0 = 0
                        VI1 = 1
C                       LES COORDONNEES DU POINT D'INTERSECTION
C                       DANS L'EF REFERENCE
                        NBDICHO = 0
 209                    DO 210 LL = 1,3
                           XYZEFR(LL,NBINT) = COSOTE(LL,N1) * (1.0-VI) +
     %                                        COSOTE(LL,N2) * VI
 210                    ENDDO
C
C                       SI FE N'EST PAS AFFINE IL FAUT TROUVER LE PT
C                       D'INTERSECTION VERIFIANT LA COORDONNEE
                        IF( NBPOEL .GT. 4 ) THEN
C
C                          LA VALEUR DES FONCTIONS DE BASE AU POINT XYZEFR
                           CALL INTERP( NOINTE, XYZEFR(1,NBINT),
     %                                  XYZEFR(2,NBINT),XYZEFR(3,NBINT),
     %                                  NBPOEL, FBASE )
C                          LA TEMPERATURE DU SOMMET DANS L'EF COURANT
                           XC = REAL( PROSCD( FBASE, SOLEL, NBPOEL ) )
                           XC = ( XC - VALISK ) / V21
                           IF( ABS(XC) .GT. 0.001 ) THEN
C                             ITERATION DE DICHOTOMIE
                              NBDICHO = NBDICHO + 1
                              IF( NBDICHO .LT. 64 ) THEN
                                 IF( XC .GE. 0 ) THEN
                                    VI1 = VI
                                 ELSE
                                    VI0 = VI
                                 ENDIF
                                 VI = ( VI0 + VI1 ) * 0.5
                                 GOTO 209
                              ENDIF
                           ENDIF
                        ENDIF
                     ENDIF
C
                  ELSE
C
C                    V1 = V2 SONT ILS EGAUX A VALISK ?
C                    ---------------------------------
                     IF( ABS( VALISK-V1 ) .LE. EPS ) THEN
C                       L'ARETE EST SUR LA SURFACE ISOVALEUR
                        NBINT = NBINT + 2
                        DO 220 LL = 1,3
                           XYZEFR(LL,NBINT-1) = COSOTE(LL,N1)
                           XYZEFR(LL,NBINT  ) = COSOTE(LL,N2)
 220                    ENDDO
                     ENDIF

                  ENDIF

C                 FIN DE TRAITEMENT DE L'ARETE I
 250           ENDDO

C              TRANSFORMATION DES 3 OU 4 SEGMENTS D'INTERSECTION
C              EN UN TRIANGLE OU QUADRANGLE = FACE DE L'ISOVALEUR
C              DANS L'EF DE REFERENCE
               CALL SETRQU( NBINT, XYZEFR )

               IF( NBINT .NE. 3 .AND. NBINT .NE. 4 ) GOTO 300

C              UNE FACE DE PLUS SUR LA SURFACE ISOVALEUR K
               IF( NBFISO .GE. MXFISO ) THEN
                  WRITE(IMPRIM,*) 'TRISO3: AUGMENTER MXFISO=',MXFISO,
     %                            ' NEF=',NEF,' NBELEM=',NBELEM
                  RETURN
               ENDIF
               NBFISO = NBFISO + 1

C              CETTE FACE NBFISO EST SUR LA SURFACE ISOVALEUR K
               NUMISO( NBFISO ) = K

C              VALEUR TEMOIN DE FACE=TRIANGLE
C              CETTE VALEUR SERA ECRASEE SI CETTE FACE EST UN QUADRANGLE
               XYZISO(3,4,NBFISO) = TRIANG

C              LES COORDONNEES DES NBINT SOMMETS DE LA FACE ISO DANS L'OBJET
               DO 295 I = 1,NBINT
                  CALL FPOXYZ( NDIM,NBPOEL,NOINTF,COPOE,XYZEFR(1,I),
     %                         FBASE,XYZISO(1,I,NBFISO) )
 295           ENDDO
            ENDIF
 300     ENDDO

C        REMONTEE POUR TRAITER LE SOUS-TETRAEDRE SUIVANT
C        ===============================================
         GOTO 80

C        *************************************************
C        FIN DU TRAITEMENT DES SOUS-TETRAEDRES DE L'EF NEF
C        *************************************************
 1000 ENDDO
C
 9999 RETURN
      END
